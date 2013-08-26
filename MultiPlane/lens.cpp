/*
 * halo_model.cpp
 *
 *  Created on: Jan 14, 2012
 *      Author: mpetkova, bmetcalf
 */

#include "slsimlib.h"

#include <algorithm>

using namespace std;

#define MIN_PLANE_DIST 1E-8

namespace
{
	class lens_halo_less
	{
	public:
		lens_halo_less(COSMOLOGY* c) : cosmo(c) {}
		
		bool operator()(const LensHalo* a, const LensHalo* b)
		{
			// compare sizes and check that b is not in eps around a
			return (a->getZlens() < b->getZlens()) && std::abs(cosmo->coorDist(a->getZlens(), b->getZlens())) > MIN_PLANE_DIST;
		}
		
	private:
		COSMOLOGY* cosmo;
	};
}

/**
 * \brief Creates an empty lens. Main halos and field halos need to be inserted by hand from the user.
 */
Lens::Lens(long *my_seed)
: seed(my_seed), halo_pos(0)
{
	cosmo = new COSMOLOGY();
	
	if( (cosmo->getOmega_matter() + cosmo->getOmega_lambda()) != 1.0 ){
		printf("ERROR: MultiLens can only handle flat universes at present.  Must change cosmology.\n");
		exit(1);
	}
	
	read_sim_file = false;
	
	charge = 4*pi*Grav;
	std::cout << "charge: " << charge << std::endl;
	
	// initially let source be the one inputed from parameter file
	index_of_new_sourceplane = -1;
	toggle_source_plane = false;
	
	std:: cout << " done " << std:: endl;
}

/**
 * \ingroup Constructor
 * \brief allocates space for the halo trees and the inout lens, if there is any
 */
Lens::Lens(InputParams& params, Source* source, long* my_seed,CosmoParamSet cosmoset)
: seed(my_seed), halo_pos(0)
{
	cosmo = new COSMOLOGY(cosmoset);
	readCosmology(params);

	if( (cosmo->getOmega_matter() + cosmo->getOmega_lambda()) != 1.0 ){
		printf("ERROR: MultiLens can only handle flat universes at present.  Must change cosmology.\n");
		exit(1);
	}

	assignParams(params);

	read_sim_file = false;

	charge = 4*pi*Grav;
	std::cout << "charge: " << charge << std::endl;

	// initially let source be the one inputed from parameter file
	index_of_new_sourceplane = -1;
	toggle_source_plane = false;

	// set up the lens contents
	buildPlanes(params, source);

	std:: cout << " done " << std:: endl;
}

Lens::~Lens()
{
	Utilities::delete_container(lensing_planes);

	Utilities::free_PosTypeMatrix(halo_pos, field_halos.size(), 3);

	Utilities::delete_container(main_halos_created);
	Utilities::delete_container(field_halos);
}

/// read in Cosmological Parameters
void Lens::readCosmology(InputParams& params){
	double tmp;

	if( params.get("Omega_matter",tmp) ) cosmo->setOmega_matter(tmp,true);
	if( params.get("Omega_lambda",tmp) ) cosmo->setOmega_lambda(tmp);
	if( params.get("Omega_baryon",tmp) ) cosmo->setOmega_baryon(tmp);
	if( params.get("Omega_neutrino",tmp) ) cosmo->setOmega_neutrino(tmp);
	if( params.get("hubble",tmp) ) cosmo->sethubble(tmp);
	if( params.get("sigma_8",tmp) ) cosmo->power_normalize(tmp);
}

/// Retrieve input parameters for construction
void Lens::assignParams(InputParams& params)
{
	if(!params.get("main_halo_on",flag_switch_main_halo_on))
	{
		ERROR_MESSAGE();
		cout << "parameter main_halo_on needs to be set in the parameter file " << params.filename() << endl;
		exit(0);
	}
	
	if(flag_switch_main_halo_on)
	{
		if(!params.get("main_DM_halo_type",main_halo_type))
		{
			ERROR_MESSAGE();
			cout << "parameter main_DM_halo_type needs to be set in the parameter file " << params.filename() << endl;
			exit(0);
		}
		if(!params.get("main_galaxy_halo_type",main_galaxy_halo_type))
		{
			main_galaxy_halo_type = null_gal;
		}
	}
	
	if(!params.get("redshift_planes_file",redshift_planes_file))
		read_redshift_planes = false;
	else
		read_redshift_planes = true;
	
	if(!params.get("field_off",flag_switch_field_off))
	{
		flag_switch_field_off = false;
		std::cout << "parameter field_off needs to be set in the parameter file " << params.filename() << std::endl;
		throw runtime_error("need input parameter");
	}
	else
	{
		if(!flag_switch_field_off)
		{
			if(!params.get("field_Nplanes",field_Nplanes))
			{
				ERROR_MESSAGE();
				cout << "parameter field_Nplanes needs to be set in the parameter file " << params.filename() << endl;
				exit(0);
			}
			
			if(!params.get("field_fov",fieldofview))
			{
				ERROR_MESSAGE();
				cout << "parameter field_fov needs to be set in the parameter file " << params.filename() << endl;
				exit(0);
			}
			
			if(!params.get("field_internal_profile",field_int_prof_type))
			{
				ERROR_MESSAGE();
				cout << "parameter field_internal_profile needs to be set in the parameter file " << params.filename() << endl;
				exit(0);
			}
			
			if(!params.get("field_internal_profile_galaxy",field_int_prof_gal_type))
			{
				flag_field_gal_on = false;
				field_int_prof_gal_type = null_gal;
			}
			else
			{
				flag_field_gal_on = true;
				
			}
			
			if(!params.get("field_mass_func_alpha",mass_func_PL_slope))
				mass_func_PL_slope = 1./6.;
			if(!params.get("field_prof_internal_slope_pl",field_prof_internal_slope) && field_int_prof_type == pl_lens)
				field_prof_internal_slope = -1.0;
			if(!params.get("field_prof_internal_slope_pnfw",field_prof_internal_slope) && field_int_prof_type == pnfw_lens)
				field_prof_internal_slope = 2.0;
			
			if(!params.get("field_input_simulation_file",field_input_sim_file))
			{
				// No simulation input file provided
				sim_input_flag = false;
				
				if(!params.get("field_mass_func_type",field_mass_func_type))
				{
					ERROR_MESSAGE();
					cout << "parameter field_mass_func_type needs to be set in the parameter file " << params.filename() << endl;
					exit(0);
				}
				
				if(!params.get("field_min_mass",field_min_mass))
				{
					ERROR_MESSAGE();
					cout << "parameter field_min_mass needs to be set in the parameter file " << params.filename() << endl;
					exit(0);
				}
				
				if(!params.get("field_buffer",field_buffer))
				{
					field_buffer = 0.0;
					cout << "default field buffer of 0 Mpc is being used." << endl;
				}
			}
			else
			{
				field_min_mass = 0.0;
				sim_input_flag = true;
			}
		}
	}
	
	if(!params.get("z_source",zsource))
	{
		ERROR_MESSAGE();
		cout << "parameter z_source needs to be set in the parameter file " << params.filename() << endl;
		exit(0);
	}
	
	if(!params.get("deflection_off",flag_switch_deflection_off))
		flag_switch_deflection_off = false;
	
	// Some checks for valid parameters
	if(flag_switch_field_off == false && field_Nplanes == 0)
	{
		ERROR_MESSAGE();
		cout << "Do you want to run _with_ field halos, but with _without_ field planes? Change field_Nplanes to a bigger number!" << endl;
		exit(1);
	}
	
	if(flag_switch_main_halo_on == false && flag_switch_field_off == true)
	{
		ERROR_MESSAGE();
		cout << "Do you want an empty simulation? Set main_halo_on to true for a main lens, or field_off to false for field lenses." << endl;
		exit(1);
	}
	
	if(!flag_switch_field_off)
	{
		if(field_int_prof_type == pl_lens && field_prof_internal_slope >= 0)
		{
			ERROR_MESSAGE();
			cout << "Power Law internal slope >=0 not possible." << endl;
			exit(1);
		}
		
		if(field_int_prof_type == pnfw_lens && field_prof_internal_slope <= 0)
		{
			ERROR_MESSAGE();
			cout << "Pseudo NFW internal slope <=0 not possible." << endl;
			exit(1);
		}
		
		if(field_int_prof_type == pnfw_lens && (field_prof_internal_slope / floor(field_prof_internal_slope) > 1.0))
		{
			ERROR_MESSAGE();
			cout << "Pseudo NFW internal slope needs to be a whole number." << endl;
			exit(1);
		}
		
		if(field_input_sim_file.size() < 1 && field_int_prof_type == nsie_lens)
		{
			ERROR_MESSAGE();
			cout << "The NSIE internal profile works only for Millennium DM simulations for now." << endl;
			cout << "Set field_input_simulation_file in sample_paramfile." << endl;
			exit(1);
		}
	}
	
	// convert to square degrees
	fieldofview /= 3600. * 3600.;
	
	printMultiLens();
}

void Lens::resetFieldNplanes(std::size_t Np)
{
	Utilities::delete_container(field_planes);
	
	field_Nplanes = Np;
	
	field_plane_redshifts.clear();
	field_Dl.clear();
	
	setupFieldPlanes();
	createFieldPlanes();
	
	combinePlanes();
}

void Lens::resetFieldHalos()
{
	Utilities::delete_container(field_halos);
	Utilities::delete_container(field_planes);
	
	Utilities::free_PosTypeMatrix(halo_pos, field_halos.size(), 3);
	
	if(sim_input_flag){
		if(read_sim_file == false) readInputSimFile();
	}
	else{
		createFieldHalos();
	}
	
	createFieldPlanes();
	
	combinePlanes();
}

void Lens::printMultiLens(){
	cout << endl << "MAIN HALOS" << endl;
	cout << "Main lens profile type:" << endl;
	switch(main_halo_type){
	case null_lens:
		cout << "no lens" << endl;
		break;
	case nfw_lens:
		cout << "NFW lens" << endl;
		break;
	case pnfw_lens:
		cout << "PseudoNFW lens" << endl;
		cout << "slope: " << field_prof_internal_slope << endl;
		break;
	case pl_lens:
		cout << "PowerLaw lens" << endl;
		cout << "slope: " << field_prof_internal_slope << endl;
		break;
	case nsie_lens:
		cout << "NSIE lens" << endl;
		break;
	case ana_lens:
		cout << "AnaNSIE lens" << endl;
		break;
	case uni_lens:
		cout << "UniNSIE lens" << endl;
		break;
	case moka_lens:
		cout << "MOKA lens" << endl;
		break;
	case dummy_lens:
		cout << "Dummy lens" << endl;
		break;
	case hern_lens:
		cout << "Hernquist lens" << endl;
		break;
	case jaffe_lens:
		cout << "Jaffe lens" << endl;
		break;
	}

	cout << endl << "Main galaxies profile type:" << endl;
	switch(main_galaxy_halo_type){
	case null_gal:
		cout << "no galaxy" << endl;
		break;
	case nsie_gal:
		cout << "NSIE galaxy" << endl;
		break;
	}

	if(flag_switch_field_off == false){

		cout << "field of view " << fieldofview << endl;

		cout << endl << "FIELD HALOS" << endl;

		cout << "field Nplanes " << field_Nplanes << endl;

		cout << "min mass " << field_min_mass << endl;
		cout << "Mass function type: "<< endl;

		switch(field_mass_func_type){
		case PS:
			cout << "  Press-Schechter mass function " << endl;
			break;
		case ST:
			cout << "  Sheth-Tormen mass function " << endl;
			break;
		case PL:
			cout << "  Power law mass function " << endl;
			cout << "  slope: " << mass_func_PL_slope << endl;
			break;
		}

		cout << endl << "Field halos profile type:" << endl;
		switch(field_int_prof_type){
		case null_lens:
			cout << "no field type" << endl;
			break;
		case nfw_lens:
			cout << "NFW field type" << endl;
			break;
		case pnfw_lens:
			cout << "PseudoNFW field type" << endl;
			cout << "slope: " << field_prof_internal_slope << endl;
			break;
		case pl_lens:
			cout << "PowerLaw field type" << endl;
			cout << "slope: " << field_prof_internal_slope << endl;
			break;
		case nsie_lens:
			cout << "NSIE field type" << endl;
			break;
		case ana_lens:
			cout << "AnaNSIE field type" << endl;
			break;
		case uni_lens:
			cout << "UniNSIE field type" << endl;
			break;
		case moka_lens:
			cout << "MOKA field type" << endl;
			break;
		case dummy_lens:
			cout << "Dummy field type" << endl;
			break;
		case hern_lens:
			cout << "Hernquist field type" << endl;
		break;
		case jaffe_lens:
			cout << "Jaffe field type" << endl;
		break;
		}

		cout << endl << "Field galaxies profile type:" << endl;
		switch(field_int_prof_gal_type){
		case null_gal:
			cout << "no field galaxy type" << endl;
			break;
		case nsie_gal:
			cout << "NSIE field galaxy type" << endl;
			break;
		}
	}

	cout << endl;
}

void Lens::setupFieldPlanes()
{
	// create spacing of lens planes or read from file
	if(read_redshift_planes)
		setFieldDistFromFile();
	else
		setFieldDist();	
}

/**
 * \brief Populates the planes with field_halos by dividing the space around the planes into
 * equal redshift distances, where the plane with the input lens is excluded
 * since it will not contain any field_halos
 *
 * Then the halo trees are built, depending on the internal profile model that
 * has been chosen in the parameter file
 */
void Lens::createFieldPlanes()
{
	std::cout << "Lens::createFieldPlanes zsource = " << zsource << std::endl;
	
	assert(field_plane_redshifts.size() == field_Nplanes);
	
	// the bounds for sorting field halos onto redshifts
	double z1 = 0, z2 = 0;
	std::size_t k1 = 0, k2 = 0;
	
	// go through planes
	for(std::size_t i = 0; i < field_Nplanes; ++i)
	{
		assert(field_plane_redshifts[i] > 0);
		assert(field_Dl[i] > 0);
		
		// previous upper bound is now lower bound
		z1 = z2;
		k1 = k2;
		
		// find upper bound
		if(i == field_Nplanes-1)
		{
			z2 = zsource;
			k2 = field_halos.size();
		}
		else
		{
			z2 = cosmo->invCoorDist(0.5*(field_Dl[i] + field_Dl[i+1]));
			k2 = Utilities::lower_bound<LensHalo>(field_halos, z2);
		}
		
		/*
		 * finding the average mass surface density in field_halos
		 */
		
		// TODO: Ben: test this
		double sigma_back = cosmo->haloMassInBufferedCone(field_min_mass,z1,z2,fieldofview*pow(pi/180,2),field_buffer,field_mass_func_type,mass_func_PL_slope)
		/(pi*pow(sqrt(fieldofview/pi)*pi*field_Dl[i]/180/(1+field_plane_redshifts[i]) + field_buffer,2));
		
		double sb=0.0;
		//double max_r = 0,tmp;
		for(std::size_t j = k1; j < k2; ++j)
		{
			sb += field_halos[j]->get_mass();
			
			//** test lines********
			//if(max_r < (tmp = halo_pos[j][0]*halo_pos[j][0] + halo_pos[j][1]*halo_pos[j][1])) max_r = tmp;
			
			// convert to proper distance on the lens plane
			halo_pos[j][0] *= field_Dl[i]/(1+field_plane_redshifts[i]);
			halo_pos[j][1] *= field_Dl[i]/(1+field_plane_redshifts[i]);
		}
		
		//max_r=sqrt(max_r);
		
		sb /= (pi*pow(sqrt(fieldofview/pi)*pi*field_Dl[i]/180/(1+field_plane_redshifts[i]) + field_buffer,2));
		
		std::cout << "sigma_back from mass function " << sigma_back << " from sum of halos " << sb << " " << sb/sigma_back - 1 << std::endl;
		if(sim_input_flag) sigma_back = sb;
		
		/*
		 * create the lensing plane
		 */
		
		std::cout << "  Building lensing plane " << i << " number of halos: " << k2-k1 << std::endl;
		
		field_planes.push_back(new LensPlaneTree(&halo_pos[k1], &field_halos[k1], k2-k1, sigma_back));
	}
	
	assert(field_planes.size() == field_Nplanes);
}

void Lens::addMainHaloToPlane(LensHalo* halo)
{
	// the redshift and distance of the halo
	double halo_z = halo->getZlens();
	double halo_Dl = cosmo->coorDist(0, halo_z);
	
	// find the position of the new lens plane
	std::size_t i = std::distance(main_Dl.begin(), std::upper_bound(main_Dl.begin(), main_Dl.end(), halo_Dl));
	
	// go though all options for adding
	if(i > 0 && (halo_Dl - main_Dl[i-1]) < MIN_PLANE_DIST)
	{
		// add to plane at (i-1)
		main_planes[i-1]->addHalo(halo);
	}
	else if(i == main_Dl.size())
	{
		// add new plane at the end
		main_planes.push_back(new LensPlaneSingular(&halo, 1));
		main_plane_redshifts.push_back(halo_z);
		main_Dl.push_back(halo_Dl);
	}
	else if((main_Dl[i] - halo_Dl) < MIN_PLANE_DIST)
	{
		// add to existing plane at position i
		main_planes[i]->addHalo(halo);
	}
	else
	{
		// create new plane at position i
		main_planes.insert(main_planes.begin() + i, new LensPlaneSingular(&halo, 1));
		main_plane_redshifts.insert(main_plane_redshifts.begin() + i, halo_z);
		main_Dl.insert(main_Dl.begin() + i, halo_Dl);
	}
}

void Lens::createMainPlanes()
{
	// clear arrays
	main_planes.clear();
	main_plane_redshifts.clear();
	main_Dl.clear();
	
	// sort halos by redshift
	std::sort(main_halos.begin(), main_halos.end(), lens_halo_less(cosmo));
	
	// put everything with same redshift (within epsilon) onto a plane

	std::vector<LensHalo*>::iterator it = main_halos.begin(); 
  //Utilities::MixedVector<LensHalo*>::iterator<> it = main_halos.begin(); temp_tag
	while(it != main_halos.end())
	{
		// find halos with higher redshift
		std::vector<LensHalo*>::iterator jt = std::upper_bound(it, main_halos.end(), *it, lens_halo_less(cosmo));
		//Utilities::MixedVector<LensHalo*>::iterator<> jt = std::upper_bound(it, main_halos.end(), *it, lens_halo_less(cosmo)); temp_tag
		
		// add halos until higher redshift to plane
		main_planes.push_back(new LensPlaneSingular(&(*it), std::distance(it, jt)));
		
		// add plane to arrays
		main_plane_redshifts.push_back((*it)->getZlens());
		main_Dl.push_back(cosmo->coorDist(0, main_plane_redshifts.back()));
		
		// advance iterator
		it = jt;
	}
}

/**
 * /brief Calculate the coordinate distances of the field planes.
 * 
 * Set the redshifts and distances of the field planes by dividing the
 * coordinate distance space into equal intervals.
 */
void Lens::setFieldDist()
{
	double Dmax = cosmo->coorDist(0, zsource);
	
	std::vector<double> lD;
	std::size_t Np = field_Nplanes + 1;
	
	assert(Np > 1);
	
	// spaces interval equally up to the source, including 0 and Dmax
	// therefore we need field_Nplanes+1 values
	fill_linear(lD, Np, 0.0, Dmax);
	
	// spacing of the distances
	double dlD = lD[1]-lD[0];
	
	// assigns the redshifts and plugs in the input plane
	for(std::size_t i = 1; i < Np; ++i)
	{
		// ensures that the first plane and the last before the source plane
		// have the same volume as all the other planes
		lD[i] -= 0.5*dlD;
		
		field_Dl.push_back(lD[i]);
	}
	
	assert(field_Dl.size() == field_Nplanes);
	
	// assigns the redshifts and plugs in the input plane
	for(std::size_t i = 0; i < field_Nplanes; ++i)
	{
		// get redshift for calculated distance
		double z = cosmo->invCoorDist(field_Dl[i]);
		field_plane_redshifts.push_back(z);
		
		// refit the distances to match the redshift
		field_Dl[i] = cosmo->coorDist(0, z);
	}
	
	assert(field_plane_redshifts.size() == field_Nplanes);
}

void Lens::setFieldDistFromFile()
{
	double value;
	
	std::ifstream file_in(redshift_planes_file.c_str());
	if(!file_in)
		throw std::runtime_error("Can't open file " + redshift_planes_file);
	
	while(file_in >> value)
	{
		if(!value)
			throw std::runtime_error("can't read double from " + redshift_planes_file);
		else
			field_plane_redshifts.push_back(value);
	}
	
	file_in.close();
	
	assert(field_plane_redshifts.size() == field_Nplanes);
	
	for(std::size_t i = 0; i < field_plane_redshifts.size(); ++i)
		field_Dl.push_back(cosmo->coorDist(0, field_plane_redshifts[i]));
}

/**
 * \brief Creates main lens halo as set up in the parmeter file.
 *
 */
void Lens::createMainHalos(InputParams& params, Source* source)
{
	switch(main_halo_type)
	{
	case null_lens:
		break;
	case nfw_lens:
		main_halos.push_back(new LensHaloNFW(params));
		break;
	case pnfw_lens:
		main_halos.push_back(new LensHaloPseudoNFW(params));
		break;
	case pl_lens:
		main_halos.push_back(new LensHaloPowerLaw(params));
		break;
	case nsie_lens:
		main_halos.push_back(new LensHaloSimpleNSIE(params));
		break;
	case ana_lens:
		main_halos.push_back(new LensHaloAnaNSIE(params));
		break;
	case uni_lens:
		main_halos.push_back(new LensHaloUniform(params));
		break;
	case moka_lens:
		{
			LensHaloMOKA* moka = new LensHaloMOKA(params);
			fieldofview = pow(1.5*moka->map->boxlrad*180/pi,2.0);
			main_halos.push_back(moka);
		}
		break;
	case dummy_lens:
		main_halos.push_back(new LensHaloDummy(params));
		break;
	case hern_lens:
		main_halos.push_back(new LensHaloHernquist(params));
		break;
	case jaffe_lens:
		main_halos.push_back(new LensHaloJaffe(params));
		break;
	}

	if(main_galaxy_halo_type!=0){
		switch(main_galaxy_halo_type){
		case null_gal:
			break;
		case nsie_gal:
			main_halos.push_back(new LensHaloSimpleNSIE(params));
			break;
		}
	}

	for(std::size_t i = 0; i < main_halos.size(); ++i)
		main_halos[i]->setInternalParams(cosmo);
}


void Lens::clearMainHalos()

{
	Utilities::delete_container(main_halos_created);
	main_halos.clear();
	
	flag_switch_main_halo_on = false;

	
	Utilities::delete_container(main_planes);
	main_plane_redshifts.clear();
	main_Dl.clear();
	
	combinePlanes();
}

/**
 * \brief Inserts a single main lens halo.
 * Then all lensing planes are updated accordingly.
 */
void Lens::insertMainHalo(Source* source, LensHalo* halo)
{
	halo->setInternalParams(cosmo);
	main_halos.push_back(halo);
	
	flag_switch_main_halo_on = true;
	
	addMainHaloToPlane(halo);
	
	combinePlanes();
}

/**
 * \brief Inserts a sequense of main lens halos and ads them to the existing ones.
 * Then all lensing planes are updated accordingly.
 */
void Lens::insertMainHalos(Source* source, LensHalo** halos, std::size_t Nhalos)
{
	for(std::size_t i = 0; i < Nhalos; ++i)
	{
		halos[i]->setInternalParams(cosmo);
		main_halos.push_back(halos[i]);
		addMainHaloToPlane(halos[i]);
	}
	
	flag_switch_main_halo_on = true;
	
	combinePlanes();
}

/**
 * \brief Inserts a single main lens halo and deletes all previous ones.
 * Then all lensing planes are updated accordingly.
 */
void Lens::replaceMainHalos(Source* source, LensHalo* halo)
{
	Utilities::delete_container(main_halos_created);
	main_halos.clear();
	
	halo->setInternalParams(cosmo);
	main_halos.push_back(halo);
	
	flag_switch_main_halo_on = true;
	
	Utilities::delete_container(main_planes);
	createMainPlanes();
	combinePlanes();
}

/**
 * \brief Inserts a sequense of main lens halos and deletes all previous ones.
 * Then all lensing planes are updated accordingly.
 */
void Lens::replaceMainHalos(Source* source, LensHalo** halos, std::size_t Nhalos)
{
	Utilities::delete_container(main_halos_created);
	main_halos.clear();
	
	for(std::size_t i = 0; i < Nhalos; ++i)
	{
		halos[i]->setInternalParams(cosmo);
		main_halos.push_back(halos[i]);
	}
	
	flag_switch_main_halo_on = true;
	
	Utilities::delete_container(main_planes);
	createMainPlanes();
	combinePlanes();
}

void Lens::createFieldHalos()
{
	const int Nzbins=64;
	const int Nmassbin=64;
	int NZSamples = 50;
	std::vector<double> zbins,Nhalosbin(Nzbins);
	unsigned long i,k,j_max,k1,k2;
	std::vector<double> Logm;
	//double pos_max[2];
    double z_max;
	const double MaxLogm=16.;
	double z1, z2, mass_max,Nhaloestot;
	int np;
	double rr,theta,maxr;
	HALO *halo_calc = new HALO(cosmo,field_min_mass,0.0);
    double mo=7.3113e10,M1=2.8575e10,gam1=7.17,gam2=0.201,be=0.557;
    double field_galaxy_mass_fraction = 0;


    if (field_min_mass < 1.0e5) {
       std::cout << "Are you sure you want the minimum field halo mass to be " << field_min_mass << " Msun?" << std::endl;
       throw;
    }
    
	double aveNhalos = cosmo->haloNumberInBufferedCone(field_min_mass,0,zsource,fieldofview*pow(pi/180,2),field_buffer,field_mass_func_type,mass_func_PL_slope);

	fill_linear(zbins,Nzbins,0.0,zsource);
	// construct redshift distribution table
	Nhalosbin[0] = 1;
	zbins[0] = 0;

	for(k=1;k<Nzbins-1;++k){
		Nhalosbin[k] = cosmo->haloNumberInBufferedCone(field_min_mass,zbins[k],zsource,fieldofview*pow(pi/180,2),field_buffer,field_mass_func_type,mass_func_PL_slope)/aveNhalos;
	}
	zbins[Nzbins-1] = zsource;
	Nhalosbin[Nzbins-1] = 0.0;

	std::size_t Nhalos = static_cast<std::size_t>(poidev(float(aveNhalos), seed));

	std::vector<double> halo_zs_vec;
	std::vector<double *> halo_pos_vec;

	// assign redsshifts to field_halos according to the redshift distribution

	for(i=0;i < Nhalos;++i){
		halo_zs_vec.push_back(InterpolateYvec(Nhalosbin,zbins,ran2(seed)));
	}

    // sort redshifts
	std::sort(halo_zs_vec.begin(),halo_zs_vec.end());

	assert(halo_zs_vec[0] < halo_zs_vec[1]);
	assert(halo_zs_vec[0] < halo_zs_vec[Nhalos-1]);

	// fill the log(mass) vector
	Logm.resize(Nmassbin);
	Nhalosbin.resize(Nmassbin);
	fill_linear(Logm,Nmassbin,log10(field_min_mass),MaxLogm);

	double *theta_pos,*theta2;
	int j = 0;
	k2 = 0;
	std::vector<double>::iterator it1,it2;
	for(np=0,mass_max=0;np<NZSamples;np++){

		z1 = np*zsource/(NZSamples);
		z2 = (np+1)*zsource/(NZSamples);

		it1 = std::lower_bound(halo_zs_vec.begin(),halo_zs_vec.end(),z1);
		it2 = std::lower_bound(halo_zs_vec.begin(),halo_zs_vec.end(),z2);

		k1 = it1 - halo_zs_vec.begin();
		k2 = it2 - halo_zs_vec.begin();

		Nhaloestot = cosmo->haloNumberInBufferedCone(pow(10,Logm[0]),z1,z2,fieldofview*pow(pi/180,2),field_buffer,field_mass_func_type,mass_func_PL_slope);

		Nhalosbin[0] = 1;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(k)
#endif
		for(k=1;k<Nmassbin-1;k++){
			// cumulative number density in one square degree
			Nhalosbin[k] = cosmo->haloNumberInBufferedCone(pow(10,Logm[k]),z1,z2,fieldofview*pow(pi/180,2),field_buffer,field_mass_func_type,mass_func_PL_slope)
					/Nhaloestot;
		}
		Nhalosbin[Nmassbin-1] = 0;

		for(i = k1; i < k2; i++){
			double Ds = cosmo->angDist(0,halo_zs_vec[i]);

			maxr = pi*sqrt(fieldofview/pi)/180. + field_buffer/Ds; // fov is a circle
			rr = maxr*sqrt(ran2(seed));

			assert(rr == rr);

			theta_pos = new double[3];

			theta = 2*pi*ran2(seed);

			theta_pos[0] = rr*cos(theta);//*Ds;
			theta_pos[1] = rr*sin(theta);//*Ds;
			theta_pos[2] = 0.0;

			switch(field_int_prof_type){
			case null_lens:
				ERROR_MESSAGE();
				std::cout << "field_int_prof_type is null!!!!" << std::endl;
				break;
			case nfw_lens:
				field_halos.push_back(new LensHaloNFW);
				break;
			case pnfw_lens:
				field_halos.push_back(new LensHaloPseudoNFW);
				break;
			case pl_lens:
				field_halos.push_back(new LensHaloPowerLaw);
				break;
			case nsie_lens:
				field_halos.push_back(new LensHaloSimpleNSIE);
				break;
			case ana_lens:
				ERROR_MESSAGE();
				std::cout << "AnaNSIE not supported." << std::endl;
				break;
			case uni_lens:
				ERROR_MESSAGE();
				std::cout << "UniNSIE not supported." << std::endl;
				break;
			case moka_lens:
				ERROR_MESSAGE();
				std::cout << "MOKA not supported." << std::endl;
				break;
			case dummy_lens:
				field_halos.push_back(new LensHaloDummy);
				break;
			case hern_lens:
				field_halos.push_back(new LensHaloHernquist);
				break;
			case jaffe_lens:
				field_halos.push_back(new LensHaloJaffe);
				break;
			}

			float mass = pow(10,InterpolateYvec(Nhalosbin,Logm,ran2 (seed)));

			halo_calc->reset(mass,halo_zs_vec[i]);

			float Rmax = halo_calc->getRvir();
			float rscale = Rmax/halo_calc->getConcentration(0);
            assert(rscale < Rmax);
      

			field_halos[j]->setZlens(halo_zs_vec[i]);
			if(flag_field_gal_on){
                field_galaxy_mass_fraction = 2*mo*pow(mass/M1,gam1)
                /pow(1+pow(mass/M1,be),(gam1-gam2)/be)/mass;
                if(field_galaxy_mass_fraction > 1.0) field_galaxy_mass_fraction = 1;

				field_halos[j]->initFromMassFunc(mass*(1-field_galaxy_mass_fraction),Rmax,rscale,field_prof_internal_slope,seed);
			}
			else{
				field_halos[j]->initFromMassFunc(mass,Rmax,rscale,field_prof_internal_slope,seed);
			}

			if(mass > mass_max) {
				mass_max = mass;
				j_max = i;
				//pos_max[0] = theta_pos[0];
				//pos_max[1] = theta_pos[1];
				z_max = halo_zs_vec[i];
			}

			halo_pos_vec.push_back(theta_pos);

			++j;

			if(flag_field_gal_on){
				switch(field_int_prof_gal_type){
				case null_gal:
					ERROR_MESSAGE();
					std::cout << "flag_field_gal_on is true, but field_int_prof_gal_type is null!!!!" << std::endl;
					break;
				case nsie_gal:
					field_halos.push_back(new LensHaloSimpleNSIE);
					break;
				}

				field_halos[j]->setZlens(halo_zs_vec[i]);
				field_halos[j]->initFromMassFunc(mass*field_galaxy_mass_fraction,Rmax,rscale,field_prof_internal_slope,seed);

                // Another copy of this position must be made to avoid rescaling it twice when it is converted into
                // distance on the lens plane in Lens::buildLensPlanes()
                theta2 = new double[3];
                theta2[0]=theta_pos[0]; theta2[1]=theta_pos[1]; theta2[2]=theta_pos[2];

				halo_pos_vec.push_back(theta2);

				++j;
			}

		}

		Nhalosbin.empty();
	}

	assert(k2 == Nhalos);
	delete halo_calc;

	std::cout << Nhalos << " halos created." << std::endl;

	Nhalos = field_halos.size();
	halo_pos = Utilities::PosTypeMatrix(Nhalos,3);

	for(i=0;i<Nhalos;++i){
		halo_pos[i] = halo_pos_vec[i];
	}

	std::cout << "leaving Lens::createFieldHalos()" << std::endl;
}

/**
 * \brief Read in information from a Virgo Millennium Data Base http://gavo.mpa-garching.mpg.de/MyMillennium/
 *
 * query select * from MockLensing.dbo.m21_20_39_021_bc03_Ben_halos
 *
 * This is information on the dark matter field_halos only.  There are 13 entries in each line separated by commas.
 * The comments must be removed from the beginning of the data file and the total number of field_halos must be added
 * as the first line.
 */
void Lens::readInputSimFile()
{
	double ra,dec,z,vmax,vdisp,r_halfmass;
	unsigned long i,j;
	unsigned long haloid,idd,np;
	double mo=7.3113e10,M1=2.8575e10,gam1=7.17,gam2=0.201,be=0.557;
    double field_galaxy_mass_fraction = 0;

	double rmax2=0,rtmp=0;

	std::ifstream file_in(field_input_sim_file.c_str());
	if(!file_in){
		std::cout << "Can't open file " << field_input_sim_file << std::endl;
    ERROR_MESSAGE();
    throw std::runtime_error(" Cannot open file.");
		exit(1);
	}

	// skip through header information in data file
	i=0;
	while(file_in.peek() == '#'){
		file_in.ignore(10000,'\n');
		++i;
	}
	std::cout << "skipped "<< i << " comment lines in file " << field_input_sim_file << std::endl;

	std::vector<double *> halo_pos_vec;

	// read in data
	int j_max;
	double mass_max=0,R_max=0,V_max=0,minmass=1e30;
	double *theta,*theta2;
	int ncolumns = 9;
	//int ncolumns = 13;

	void *addr[ncolumns];
  
	addr[0] = &haloid;
	addr[1] = &idd;
	addr[2] = &ra;
	addr[3] = &dec;
	addr[4] = &z;
	addr[5] = &np;
	addr[6] = &vdisp;
	addr[7] = &vmax;
	addr[8] = &r_halfmass;
/*
  double z_app,m_crit200,m_mean200;
  unsigned long galid;
  
  addr[0] = &galid;
  addr[1] = &haloid;
	addr[2] = &idd;
	addr[3] = &ra;
	addr[4] = &dec;
	addr[5] = &z;
  addr[6] = &z_app;
	addr[7] = &np;
  addr[8] = &m_crit200;
  addr[9] = &m_mean200;
	addr[10] = &vmax;
  addr[11] = &vdisp;
	addr[12] = &r_halfmass;
  */
  
	unsigned long myint;
	double mydouble;
	std::string myline;
	std::string strg;
	std::string f=",";
	std::stringstream buffer;

	for(i=0,j=0 ; ; ++i){
		// read a line of data
		myline.clear();
		std::getline(file_in,myline);

		if(myline[0] == '#')
			break;
		for(int l=0;l<ncolumns; l++){
			int pos = myline.find(f);
			strg.assign(myline,0,pos);
			buffer << strg;
			if(l <= 1 || l == 5){
      //  if(l <= 2 || l == 7){
				buffer >> myint;
				*((unsigned long *)addr[l]) = myint;
			}
			else{
				buffer >> mydouble;
				*((double *)addr[l]) = mydouble;
			}
			myline.erase(0,pos+1);
			strg.clear();
			buffer.clear();
			buffer.str(std::string());
		}

		// position on lens plane
		theta = new double[2];
		//double Ds = cosmo->angDist(0,z);
		theta[0] = -ra*pi/180.;
		theta[1] = dec*pi/180.;

    if(rmax2 < (rtmp = (theta[0]*theta[0]+theta[1]*theta[1]))) rmax2 = rtmp;
    
		if(np > 0.0 && vdisp > 0.0 && z <= zsource){

			switch(field_int_prof_type){
			case null_lens:
				ERROR_MESSAGE();
				std::cout << "field_int_prof_type is null!!!!" << std::endl;
				break;
			case nfw_lens:
				field_halos.push_back(new LensHaloNFW);
				break;
			case pnfw_lens:
				field_halos.push_back(new LensHaloPseudoNFW);
				break;
				ERROR_MESSAGE();
				std::cout << "PseudoNFW not supported." << std::endl;
				break;
			case pl_lens:
				ERROR_MESSAGE();
				std::cout << "PowerLaw not supported." << std::endl;
				break;
			case nsie_lens:
				field_halos.push_back(new LensHaloSimpleNSIE);
				break;
			case ana_lens:
				ERROR_MESSAGE();
				std::cout << "AnaNSIE not supported." << std::endl;
				break;
			case uni_lens:
				ERROR_MESSAGE();
				std::cout << "UniNSIE not supported." << std::endl;
				break;
			case moka_lens:
				ERROR_MESSAGE();
				std::cout << "MOKA not supported." << std::endl;
				break;
			case dummy_lens:
				field_halos.push_back(new LensHaloDummy);
				ERROR_MESSAGE();
				std::cout << "Why would you want dummy file halos???" << std::endl;
				break;
			case hern_lens:
				ERROR_MESSAGE();
				std::cout << "Hernquist not supported." << std::endl;
				break;
			case jaffe_lens:
				ERROR_MESSAGE();
				std::cout << "Jaffe not supported." << std::endl;
				break;
			}

			float mass = np*8.6e8/cosmo->gethubble();

			field_halos[j]->setZlens(z);
			if(flag_field_gal_on){
				field_galaxy_mass_fraction = 2*mo*pow(mass/M1,gam1)
				  /pow(1+pow(mass/M1,be),(gam1-gam2)/be)/mass;
				if(field_galaxy_mass_fraction > 1.0) field_galaxy_mass_fraction = 1;

				field_halos[j]->initFromFile(mass*(1-field_galaxy_mass_fraction),seed,vmax,r_halfmass*cosmo->gethubble());
			}
			else{
        // TODO: test line ****
        //std::cout << "Warning:  This needs to be changed" << std::endl;
        //mass /= 10;
        /***********************/
				field_halos[j]->initFromFile(mass,seed,vmax,r_halfmass*cosmo->gethubble());
			}


			if(field_halos[j]->get_Rmax() > R_max) R_max = field_halos[j]->get_Rmax();
			if(vdisp > V_max) V_max = vdisp;

			halo_pos_vec.push_back(theta);

			if(mass > mass_max) {
				mass_max = mass;
				j_max = j;
			}
			if(mass < minmass) {
				minmass = mass;
			}

			++j;

			if(flag_field_gal_on){
				switch(field_int_prof_gal_type){
				case null_gal:
					ERROR_MESSAGE();
					std::cout << "flag_field_gal_on is true, but field_int_prof_gal_type is null!!!!" << std::endl;
					break;
				case nsie_gal:
					field_halos.push_back(new LensHaloSimpleNSIE);
					break;
				}

				field_halos[j]->setZlens(z);
				field_halos[j]->initFromFile(mass*field_galaxy_mass_fraction,seed,vmax,r_halfmass*cosmo->gethubble());

        // Another copy of this position must be made to avoid rescaling it twice when it is converted into
        // distance on the lens plane in Lens::buildLensPlanes()
        theta2 = new double[2];
        theta2[0]=theta[0]; theta2[1]=theta[1];
				halo_pos_vec.push_back(theta2);

				++j;
			}

		}
	}
	file_in.close();
	std::cout << field_halos.size() << " halos read in."<< std::endl
			<< "Max input mass = " << mass_max << "  R max = " << R_max << "  V max = " << V_max
      << "Min imput mass = " << minmass << std::endl;

	/// setting the minimum halo mass in the simulation
	field_min_mass = minmass;
	if(field_buffer > 0.0){
		std::cout << "Overiding field_buffer to make it 0 because halos are read in." << std::endl;
		field_buffer = 0.0;
	}

	halo_pos = Utilities::PosTypeMatrix(field_halos.size(), 3);

	for(i = 0; i < field_halos.size(); ++i)
	{
		halo_pos[i] = halo_pos_vec[i];
	}

	std::cout << "Overiding input file field of view to make it fit the simulation light cone." << std::endl;
	fieldofview = pi*rmax2*pow(180/pi,2);  // Resets field of view to range of input galaxies

	std::cout << "Setting mass function to Sheth-Tormen." << std::endl;
	field_mass_func_type = ST; // set mass function

	std::cout << "sorting in Lens::readInputSimFile()" << std::endl;
	// sort the field_halos by readshift
	Lens::quicksort(field_halos.data(),halo_pos,field_halos.size());

	std::cout << "leaving Lens::readInputSimFile()" << std::endl;

  field_buffer = 0.0;
	read_sim_file = true;
}

void Lens::combinePlanes()
{
	// clear old plane configuration
	lensing_planes.clear();
	plane_redshifts.clear();
	Dl.clear();
	dDl.clear();
	
	// index of current main/field plane
	std::size_t i_field = 0, i_main = 0;
	
	// always get plane with least redshift from either field or main
	while(i_field < field_planes.size() && i_main < main_planes.size())
	{
		// decide if main or field plane is next
		if(main_plane_redshifts[i_main] <= field_plane_redshifts[i_field])
		{
			// next plane is main
			lensing_planes.push_back(main_planes[i_main]);
			plane_redshifts.push_back(main_plane_redshifts[i_main]);
			Dl.push_back(main_Dl[i_main]);
			
			// advance main index
			++i_main;
		}
		else
		{
			// next plane is field
			lensing_planes.push_back(field_planes[i_field]);
			plane_redshifts.push_back(field_plane_redshifts[i_field]);
			Dl.push_back(field_Dl[i_field]);
			
			// check if planes are too close together
			if(std::abs(field_Dl[i_field] - main_Dl[i_main]) < MIN_PLANE_DIST)
			{
				// move back the inserted field plane
				Dl.back() = main_Dl[i_main] + MIN_PLANE_DIST;
				plane_redshifts.back() = cosmo->invCoorDist(Dl.back());
				
				// TODO: make this more intelligent or make it possible to have all halos on the same planes
			}
			
			// advance field index
			++i_field;
		}
	}
	
	// add rest of planes, one array will already be at end
	for(; i_field < field_planes.size(); ++i_field)
	{
		lensing_planes.push_back(field_planes[i_field]);
		plane_redshifts.push_back(field_plane_redshifts[i_field]);
		Dl.push_back(field_Dl[i_field]);
	}
	for(; i_main < main_planes.size(); ++i_main)
	{
		lensing_planes.push_back(main_planes[i_main]);
		plane_redshifts.push_back(main_plane_redshifts[i_main]);
		Dl.push_back(main_Dl[i_main]);
	}
	
	assert(lensing_planes.size() == field_planes.size() + main_planes.size());
	
	// add the pseudo-plane for rayshooting at the end of the arrays
	plane_redshifts.push_back(zsource);
	Dl.push_back(cosmo->coorDist(0, zsource));
	
	// calculate deltas
	dDl.push_back(Dl[0]);
	for(std::size_t i = 1; i < Dl.size(); ++i)
		dDl.push_back(Dl[i] - Dl[i-1]); // distance from plane i-1 to plane i
	
	// output resulting setup
	std::cout
	<< "\ncombinePlanes()"
	<< "\n---------------"
	<< std::endl;
	std::cout << "\nz:";
	for(std::size_t i = 0, n = plane_redshifts.size(); i < n; ++i)
		std::cout << " " << plane_redshifts[i];
	std::cout << "\nDl:";
	for(std::size_t i = 0, n = Dl.size(); i < n; ++i)
		std::cout << " " << Dl[i];
	std::cout << "\ndDl:";
	for(std::size_t i = 0, n = dDl.size(); i < n; ++i)
		std::cout << " " << dDl[i];
	std::cout << "\n" << std::endl;
}

void Lens::buildPlanes(InputParams& params, Source* source)
{
	// build field
	if(!flag_switch_field_off)
	{
		// set the distances of the field planes
		setupFieldPlanes();
		
		// create or read the field halos
		if(sim_input_flag)
			readInputSimFile();
		else
			createFieldHalos();
		
		// create field planes and sort halos onto them
		createFieldPlanes();
	}
	
	// build main
	if(flag_switch_main_halo_on)
	{
		// create the main halos
		createMainHalos(params, source);
		
		// create the main planes for the halos
		createMainPlanes();
	}
	
	// combine the different planes
	combinePlanes();
}

/**
 * \brief Changes the maximum redshift that the rays are shot to.
 *
 * The multilens must have been initially constructed with a source redshift that is higher
 * than this redshift.  This is used to rayshoot to a source whose line of sight passes through the
 * simulation volume.  The source can be at higher redshift than the simulation volume.
 *
 * To revert the source redshift to its original value use Lens::RevertSourcePlane().
 *
 */
short Lens::ResetSourcePlane(
		double z                 /// redshift of implanted source
		,bool nearest           /** If true, set the source plane to the nearest (in coordinate distance)
			                      * lensing plane that was created already.  This can be used to avoid self-lensing
			                      * by the halo of the source.  If the source is at higher redshift than the simulation
			                      * volume the source will be at its real redshift.
			                      */
		,unsigned long GalID
		,double *xx
		){
	unsigned long j;
	short out=0;

	toggle_source_plane = true;

	if(z<=0.0){
		cout << "Warning: Source redshift can't be set to " << z << " in MultiLens::ResetSourcePlane." << endl;
		return out;
	}


	// distance to new source plane
	double Ds = cosmo->coorDist(0,z);
	// find bounding index
	locateD(Dl.data()-1,lensing_planes.size(),Ds,&j);
	// j is the index of the next plane at higher redshift, This plane will be temporarily replaced and used as a source plane
	assert(j <= lensing_planes.size() && j >=0);

	if(j > 0)
	{
		// check if source plane coincides with previous lens plane
		if(Dl[j-1] == Ds)
			--j;
		// or check if previous plane is nearer when asked to
		else if(nearest)
		{
			double z1 = cosmo->invCoorDist(Dl[j]-0.5*dDl[j]);
			if(z < z1) 
				--j;
		}
	}

	if(nearest && (j < lensing_planes.size()) ){
		zs_implant = plane_redshifts[j];
		Ds_implant = Dl[j];
		if(j > 0) dDs_implant = dDl[j];
		else  dDs_implant = Ds_implant;
	}else{
		// if nearest==false or the source is at higher redshift than the last plane use the real redshift
		Ds_implant = Ds;
		zs_implant = z;
		if(j > 0) dDs_implant = cosmo->coorDist(plane_redshifts[j-1],z); //Ds - Dl[j-1];
		else  dDs_implant = Ds;
	}

	std::cout << "Source on plane " << j << " zs " << zs_implant << " Ds " << Ds << " dDs " << dDs_implant << std::endl;

	index_of_new_sourceplane = j;

	out = j;
	return out;
}

/// Sort field_halos[], brr[][], and id[] by content off arr[]
void Lens::quicksort(LensHaloHndl *halos,double **pos,unsigned long N){
	double pivotvalue;
	unsigned long pivotindex,newpivotindex,i;

	if(N <= 1) return ;

	// pick pivot as the median of the first, last and middle values
	if ((halos[0]->getZlens() >= halos[N/2]->getZlens() && halos[0]->getZlens() <= halos[N-1]->getZlens())
			|| (halos[0]->getZlens() >= halos[N-1]->getZlens() && halos[0]->getZlens() <= halos[N/2]->getZlens())) pivotindex = 0;
	else if ((halos[N/2]->getZlens() >= halos[0]->getZlens() && halos[N/2]->getZlens() <= halos[N-1]->getZlens())
			|| (halos[N/2]->getZlens() >= halos[N-1]->getZlens() && halos[N/2]->getZlens() <= halos[0]->getZlens())) pivotindex = N/2;
	else pivotindex = N-1;
	pivotvalue=halos[pivotindex]->getZlens();

	// move pivet to end of halosay
	std::swap(halos[pivotindex],halos[N-1]);
	std::swap(pos[pivotindex],pos[N-1]);
	newpivotindex=0;

	// partition list and halosay
	for(i=0;i<N;++i){
		if(halos[i]->getZlens() <= pivotvalue){
			std::swap(halos[newpivotindex],halos[i]);
			std::swap(pos[newpivotindex],pos[i]);
			++newpivotindex;
		}
	}
	--newpivotindex;

	quicksort(&halos[0],pos,newpivotindex);
	quicksort(&halos[newpivotindex+1],&pos[newpivotindex+1],N-newpivotindex-1);

	return ;
}
