/*
 * halo_model.cpp
 *
 *  Created on: Jan 14, 2012
 *      Author: mpetkova, bmetcalf
 */

#include "slsimlib.h"

using namespace std;


/*
 * \ingroup Constructor
 * allocates space for the halo trees and the inout lens, if there is any
 */
Lens::Lens(InputParams& params, CosmoHndl cosmo, SourceHndl source, long *my_seed) : seed(my_seed){

	if( (cosmo->getOmega_matter() + cosmo->getOmega_lambda()) != 1.0 ){
		printf("ERROR: MultiLens can only handle flat universes at present.  Must change cosmology.\n");
		exit(1);
	}

	assignParams(params);

	/// makes the oordinate distance table for the calculation of the redshifts of the different planes
	table_set = false;
	make_table(cosmo);

	read_sim_file = false;

	charge = 4*pi*Grav*mass_scale;
	std::cout << "charge: " << charge << std::endl;

	if(flag_input_lens)
		createMainHalos(params,cosmo,source);

	// initially let source be the one inputed from parameter file
	index_of_new_sourceplane = -1;
	toggle_source_plane = false;

	seed = my_seed;

	if(read_redshift_planes){
		setCoorDistFromFile(cosmo);
	}else{
		setCoorDist(cosmo);
	}

	if(flag_switch_field_off == false){
		if(sim_input_flag){
			if(read_sim_file == false) readInputSimFile(cosmo);
		}
		else{
			createFieldHalos(cosmo,seed);
		}
	}

	buildLensPlanes(cosmo);

	std:: cout << " done " << std:: endl;
}

Lens::~Lens(){
	Utilities::delete_container(lensing_planes);

	Dl.clear();
	plane_redshifts.clear();
	dDl.clear();

	Utilities::free_PosTypeMatrix(halo_pos,Nhalos,3);

	main_halos.clear();
	Utilities::delete_container(field_halos);

	coorDist_table.clear();
}

void Lens::make_table(CosmoHndl cosmo){
	int i;
	double x, dx = (zsource+1.0)/(double)NTABLE, Dl;

	for(i = 0 ; i< NTABLE; i++){
		x = i*dx;
		Dl = cosmo->coorDist(0,x);
		coorDist_table.insert ( std::pair<double,double>(Dl,x));
	}
	//std::map<double,double>::iterator it;
	//for (it=coorDist_table.begin(); it!=coorDist_table.end(); ++it)
	//   std::cout << it->first << " => " << it->second << '\n';
	table_set=true;
}

/// Retrieve input parameters for construction
void Lens::assignParams(InputParams& params){

  if(!params.get("outputfile",outputfile)){
		  ERROR_MESSAGE();
		  cout << "parameter outputfile needs to be set in the parameter file " << params.filename() << endl;
		  exit(0);
	}
	if(!params.get("Nplanes",Nplanes)){
		  ERROR_MESSAGE();
		  cout << "parameter Nplanes needs to be set in the parameter file " << params.filename() << endl;
		  exit(0);
	}
	if(!params.get("flag_input_lens",flag_input_lens)){
		  ERROR_MESSAGE();
		  cout << "parameter flag_input_lens needs to be set in the parameter file " << params.filename() << endl;
		  exit(0);
	}
	if(flag_input_lens){
		if(!params.get("DM_halo_type",main_halo_type)){
			ERROR_MESSAGE();
			cout << "parameter flag_input_lens needs to be set in the parameter file " << params.filename() << endl;
			exit(0);
		}
		if(!params.get("galaxy_halo_type",galaxy_halo_type)){
			galaxy_halo_type = null_gal;
		}
	}

	if(!params.get("field_off",flag_switch_field_off)) flag_switch_field_off = false;

	if(!flag_switch_field_off){
		if(!params.get("internal_profile",int_prof_type)){
			ERROR_MESSAGE();
			cout << "parameter internal_profile needs to be set in the parameter file " << params.filename() << endl;
			exit(0);
		}
		if(!params.get("internal_profile_galaxy",int_prof_gal_type)){
			int_prof_gal_type = null_gal;
			flag_galaxy_subhalo = false;
		}
		else{
			flag_galaxy_subhalo = true;
			if(!params.get("galaxy_mass_fraction",galaxy_mass_fraction)){
				ERROR_MESSAGE();
				cout << "to construct a DM + galaxy model the parameter galaxy_mass_fraction needs to be set in the parameter file " << params.filename() << endl;
				exit(0);
			}
		}
		if(!params.get("alpha",mass_func_PL_slope))                mass_func_PL_slope = 1./6.;
		if(!params.get("internal_slope_pl",halo_slope) && int_prof_type == pl_lens)     halo_slope = -1.0;
		if(!params.get("internal_slope_pnfw",halo_slope) && int_prof_type == pnfw_lens) halo_slope = 2.0;

		if(!params.get("redshift_planes_file",redshift_planes_file)) read_redshift_planes = false;
		else read_redshift_planes = true;

		if(!params.get("input_simulation_file",input_sim_file)){
			sim_input_flag = false;
			// No simulation input file provided
			if(!params.get("mass_func_type",mass_func_type)){
				  ERROR_MESSAGE();
				  cout << "parameter mass_func_type needs to be set in the parameter file " << params.filename() << endl;
				  exit(0);
			}
			if(!params.get("min_mass",min_mass)){
				  ERROR_MESSAGE();
				  cout << "parameter min_mass needs to be set in the parameter file " << params.filename() << endl;
				  exit(0);
			}
			if(!params.get("fov",fieldofview)){
				  ERROR_MESSAGE();
				  cout << "parameter fov needs to be set in the parameter file " << params.filename() << endl;
				  exit(0);
			}
			if(!params.get("field_buffer",field_buffer)){
				field_buffer = 0.0;
				cout << "default field buffer of 0 Mpc is being used." << endl;
			}
		}else{
			min_mass = 0.0;
			sim_input_flag = true;
		}
	}

	if(!params.get("z_source",zsource)){
		  ERROR_MESSAGE();
		  cout << "parameter z_source needs to be set in the parameter file " << params.filename() << endl;
		  exit(0);
	}

	if(!params.get("mass_scale",mass_scale)){
		  ERROR_MESSAGE();
		  cout << "parameter mass_scale needs to be set in the parameter file " << params.filename() << endl;
		  exit(0);
	}

	if(!params.get("deflection_off",flag_switch_deflection_off)) flag_switch_deflection_off = false;

	// Some checks for valid parameters
	if(flag_switch_field_off == true && Nplanes != 1){
		ERROR_MESSAGE();
		cout << "Do you want to run _without_ field halos, but with more than one lens planes? Change Nplanes to 1!" << endl;
		exit(1);
	}

	if(flag_switch_field_off == false && Nplanes == 1){
		ERROR_MESSAGE();
		cout << "Do you want to run _with_ field halos, but with _only_ one lens planes? Change Nplanes to a bigger number!" << endl;
		exit(1);
	}

	if(flag_input_lens == 0 && Nplanes == 1){
		ERROR_MESSAGE();
		cout << "Do you want an empty simulation? Set flag_input_lens to > 0 for a main lens." << endl;
		exit(1);
	}

	if(flag_input_lens > 0){
		flag_input_lens = 1;
	}
	if(!flag_switch_field_off){
		if(int_prof_type == pl_lens && halo_slope >= 0){
			ERROR_MESSAGE();
			cout << "Power Law internal slope >=0 not possible." << endl;
			exit(1);
		}

		if(int_prof_type == pnfw_lens && halo_slope <= 0){
			ERROR_MESSAGE();
			cout << "Pseudo NFW internal slope <=0 not possible." << endl;
			exit(1);
		}

		if(int_prof_type == pnfw_lens && (halo_slope / floor(halo_slope) > 1.0)){
			ERROR_MESSAGE();
			cout << "Pseudo NFW internal slope needs to be a whole number." << endl;
			exit(1);
		}

		if(input_sim_file.size() < 1 && int_prof_type == nsie_lens){
			ERROR_MESSAGE();
			cout << "The NSIE internal profile works only for Millenium DM simulations for now." << endl;
			cout << "Set input_simulation_file in sample_paramfile." << endl;
			exit(1);
		}
	}

	  // to compensate for the last plane, which is the source plane
	  Nplanes++;

	  // to compensate for additional lens planes

	  // convert to square degrees
	  fieldofview /= 3600. * 3600.;

	  printMultiLens();
}

void Lens::resetNplanes(CosmoHndl cosmo, int Np){
	Utilities::delete_container(lensing_planes);

	Nplanes = Np;

	Dl.clear();
	plane_redshifts.clear();
	dDl.clear();

	setCoorDist(cosmo);
	buildLensPlanes(cosmo);
}

void Lens::resetFieldHalos(CosmoHndl cosmo){
	Utilities::delete_container(field_halos);
	Utilities::delete_container(lensing_planes);

	Utilities::free_PosTypeMatrix(halo_pos,Nhalos,3);

	if(sim_input_flag){
		if(read_sim_file == false) readInputSimFile(cosmo);
	}
	else{
		createFieldHalos(cosmo,seed);
	}

	buildLensPlanes(cosmo);
}

void Lens::printMultiLens(){
	cout << "Nplanes " << Nplanes << endl;

	cout << "mass scale " << mass_scale << endl;

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
		cout << "slope: " << halo_slope << endl;
		break;
	case pl_lens:
		cout << "PowerLaw lens" << endl;
		cout << "slope: " << halo_slope << endl;
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
	}

	cout << endl << "Main galaxies profile type:" << endl;
	switch(galaxy_halo_type){
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

		cout << "min mass " << min_mass << endl;
		cout << "Mass function type: "<< endl;

		switch(mass_func_type){
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
		switch(int_prof_type){
		case null_lens:
			cout << "no field type" << endl;
			break;
		case nfw_lens:
			cout << "NFW field type" << endl;
			break;
		case pnfw_lens:
			cout << "PseudoNFW field type" << endl;
			cout << "slope: " << halo_slope << endl;
			break;
		case pl_lens:
			cout << "PowerLaw field type" << endl;
			cout << "slope: " << halo_slope << endl;
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
		}

		cout << endl << "Field galaxies profile type:" << endl;
		switch(int_prof_gal_type){
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

/**
 * Populates the planes with field_halos by dividing the space around the planes into
 * equal redshift distances, where the plane with the input lens is excluded
 * since it will not contain any field_halos
 *
 * Then the halo trees are built, depending on the internal profile model that
 * has been chosen in the parameter file
 */
void Lens::buildLensPlanes(
		CosmoHndl cosmo /// the cosmology
		){
	int j, Ntot;
	double z1, z2;
	unsigned long j1,j2;
	std::map<double,double>::iterator ind;

	std::cout << "Lens::buildLensPlanes zsource = " << zsource << std::endl;

	assert(plane_redshifts[Nplanes-1] == zsource);

	for(j=0,Ntot=0;j<Nplanes-1;j++){
		if(flag_input_lens && j == (flag_input_lens % Nplanes)){
			std::cout << "  Building lensing plane " << j << " number of halos: " << main_halos.size() << std::endl;

			lensing_planes.push_back(new SingularLensPlane(main_halos.data(),main_halos.size()));
		}
		else if(flag_switch_field_off == false){


			/*
			 * Setting the redshift range
			 * If there is a plane with an input lens on it, it is skipped over
			 * since it will not contain any field_halos
			 */
			if(j == 0) z1 = 0.0;
			else{
				ind = coorDist_table.lower_bound((Dl[j]-0.5*dDl[j]));
				z1 = ind->second;
			}

			if(flag_input_lens && j-1 == (flag_input_lens % Nplanes)){
				ind = coorDist_table.lower_bound(Dl[j] - 0.5*(Dl[j] - Dl[j-2]));
				z1 = ind->second;
			}

			if(j == Nplanes-2) z2 = zsource;
			else{
				ind = coorDist_table.lower_bound(Dl[j] + 0.5*dDl[j+1]);
				z1 = ind->second;
			}

			if(flag_input_lens && j+1 == (flag_input_lens % Nplanes)){
				ind = coorDist_table.lower_bound(Dl[j] + 0.5*(Dl[j+2] - Dl[j]));
				z1 = ind->second;
			}

			/// Find which field_halos are in redshift range
			j1 = Utilities::lower_bound<LensHalo>(field_halos,z1);
			j2 = Utilities::lower_bound<LensHalo>(field_halos,z2);

			/*
			 * finding the average mass surface density in field_halos
			 */

			// TODO Ben: test this
			double sigma_back = cosmo->haloMassInBufferedCone(min_mass*mass_scale,z1,z2,fieldofview*pow(pi/180,2),field_buffer,mass_func_type,mass_func_PL_slope)
			            		 /(pi*pow(sqrt(fieldofview/pi)*pi*Dl[j]/180/(1+plane_redshifts[j]) + field_buffer,2))/mass_scale;

			double sb=0.0;
			for(int m=0;m<j2-j1;m++){
				sb+=field_halos[j1+m]->get_mass();
			}

			sb /= (pi*pow(sqrt(fieldofview/pi)*pi*Dl[j]/180/(1+plane_redshifts[j]) + field_buffer,2));

			std::cout << sigma_back << " " << sb << " " << sb/sigma_back - 1 << std::endl;
			if(sim_input_flag) sigma_back = sb;

			/// Use other constructor to create halo data
			std::cout << "  Building lensing plane " << j << " number of halos: " << j2-j1 << std::endl;

			lensing_planes.push_back(new TreeLensPlane(&halo_pos[j1],field_halos.data(),j2-j1,sigma_back));
		}
	}
}

/**
 * Set the coordinate distances of the planes by dividing the coordinate distance space into equal intervals
 * and then plugging the analytic input plane in between.
 *
 * After this flag_input_lens will hold the index of the plane it is on
 * In case it is on the first plane, it will hold the index Nplanes, to make
 * sure that it is not zero (i.e. not set)
 */
void Lens::setCoorDist(CosmoHndl cosmo){
	int i;

	double Dlens;
	double Ds = cosmo->coorDist(0,zsource);
	if(flag_input_lens) Dlens = cosmo->coorDist(0,main_halos[0]->getZlens());
	else Dlens = Ds;

	if(flag_input_lens && Nplanes == 2){
		Dl.push_back(Dlens);
		Dl.push_back(Ds);
	}else{

		std:: vector<double> lD;
		int Np;

		if(flag_input_lens == 0 && flag_switch_field_off == false)
			Np = Nplanes;
		if(flag_input_lens && flag_switch_field_off == false)
			Np = Nplanes-1;

		cout << flag_input_lens << " " << flag_switch_field_off << " " << Np << " " << Dlens << " " << Ds << endl;

		/// spaces lD equally up to the source, including 0 and Ds
		/// therefore we need Nplanes+1 values
		/// however, if there is an input plane, we will need Nplanes values, since the input plane will take up a value itself
		fill_linear(lD,Np,0.0,Ds);

		/// ensures that the first plane and the last before the source plane have the same volume
		/// as all ther planes
		double dlD = lD[1]-lD[0];
		for(i=0; i<Np; i++){
			lD[i] -= 0.5*dlD;
		}

		/// puts the input plane first if the case
		int flag=0;
		if(flag_input_lens && Dlens < lD[1]){
			Dl.push_back(Dlens);
			flag_input_lens = Nplanes;
			flag = 1;
		}

		/// assigns the redshifts and plugs in the input plane
		for(i=1; i<Np; i++){
			Dl.push_back(lD[i]);

			if(flag_input_lens && flag == 0)
				if(Dlens > lD[i] && Dlens <= lD[i+1]){
					Dl.push_back(Dlens);
					flag_input_lens = Dl.size()-1;
					flag = 1;
				}
		}

		Dl.push_back(Ds);
	}

	int j;
	dDl.push_back(Dl[0]);  // distance between jth plane and the previous plane
	for(j = 1; j < Nplanes; j++){
		dDl.push_back(Dl[j] - Dl[j-1]); // distance between jth plane and the previous plane
	}

	if(flag_input_lens)
		cout << "zlens " << main_halos[0]->getZlens() << " on plane number " << (flag_input_lens % Nplanes) << endl;

	cout << "Dl: ";
	for(j = 0; j < Nplanes; j++)
		cout << Dl[j] << " ";
	cout << endl;

	cout << "dDl: ";
	for(j = 0; j < Nplanes; j++)
		cout << dDl[j] << " ";
	cout << endl;


	std::map<double,double>::iterator ind;
	// assigns the redshifts and plugs in the input plane
	cout << "z: ";
	for(i=0; i<Nplanes-1; i++){
		ind = coorDist_table.lower_bound(Dl[i]);
		plane_redshifts.push_back(ind->second);
		cout << plane_redshifts[i] << " ";
	}
	plane_redshifts[Nplanes-1] = zsource;
	cout << plane_redshifts[i] << " " << std::endl;

	assert(Dl.size() == Nplanes);
}


void Lens::setCoorDistFromFile(CosmoHndl cosmo){

	double value;

	std::ifstream file_in(redshift_planes_file.c_str());
	if(!file_in){
		std::cout << "Can't open file " << redshift_planes_file << std::endl;
		exit(1);
	}

	while(file_in >> value){
		if(!value){
			ERROR_MESSAGE();
			cout << "can't read double from " << redshift_planes_file << endl;
			exit(1);
		}else{
			plane_redshifts.push_back(value);
		}
	}

	file_in.close();

	plane_redshifts.push_back(zsource);

	cout << "Dl: ";
	int j;
	for(j = 0; j < Nplanes; j++){
		Dl.push_back(cosmo->coorDist(0,plane_redshifts[j]));
		cout << Dl[j] << " ";
	}
	cout << endl;

	cout << "dDl: ";
	dDl.push_back(Dl[0]);  // distance between jth plane and the previous plane
	cout << dDl[0] << " ";
	for(j = 1; j < Nplanes; j++){
		dDl.push_back(Dl[j] - Dl[j-1]); // distance between jth plane and the previous plane
		cout << dDl[j] << " ";
	}
	cout << endl;

	cout << "z: ";
	for(j=0; j<Nplanes; j++){
		cout << plane_redshifts[j] << " ";
	}
	cout << endl;

	assert(plane_redshifts.size() == Nplanes);
}

void Lens::createMainHalos(
		InputParams& params
		,CosmoHndl cosmo     /// cosmology
		,SourceHndl source
){
	switch(main_halo_type){
	case null_lens:
		break;
	case nfw_lens:
		main_halos.push_back(new NFWLensHalo(params));
		break;
	case pnfw_lens:
		main_halos.push_back(new PseudoNFWLensHalo(params));
		break;
	case pl_lens:
		main_halos.push_back(new PowerLawLensHalo(params));
		break;
	case nsie_lens:
		main_halos.push_back(new SimpleNSIELensHalo(params));
		break;
	case ana_lens:
		main_halos.push_back(new AnaNSIELensHalo(params));
		break;
	case uni_lens:
		main_halos.push_back(new UniNSIELensHalo(params));
		break;
	case moka_lens:
		main_halos.push_back(new MOKALensHalo(params));
		//fieldofview = pow(1.5*mokalens->map->boxlrad*180/pi,2.0);
		break;
	case dummy_lens:
		main_halos.push_back(new DummyLensHalo(params));
		break;
	}

	if(galaxy_halo_type!=0){
		switch(galaxy_halo_type){
		case null_gal:
			break;
		case nsie_gal:
			main_halos.push_back(new SimpleNSIELensHalo(params));
			break;
		}
	}

	NmainHalos = main_halos.size();

	for(int i=0; i< NmainHalos; i++)
		main_halos[i]->setInternalParams(cosmo,source);
}

void Lens::createFieldHalos(
		CosmoHndl cosmo     /// cosmology
		,long *seed
	){

	const int Nzbins=64;
	const int Nmassbin=64;
	int NZSamples = 50;
	std::vector<double> zbins,Nhalosbin(Nzbins);
	unsigned long i,k,j_max,k1,k2;
	std::vector<double> Logm;
	double pos_max[2], z_max;
	const double MaxLogm=16.;
	double z1, z2, mass_max,mass_tot,Nhaloestot;
	int np;
	double rr,theta,maxr;
	float dummy;
	HALO *halo_calc = new HALO(cosmo,min_mass*mass_scale,0.0);

	double aveNhalos = cosmo->haloNumberInBufferedCone(min_mass*mass_scale,0,zsource,fieldofview*pow(pi/180,2),field_buffer,mass_func_type,mass_func_PL_slope);

	fill_linear(zbins,Nzbins,0.0,zsource);
	// construct redshift distribution table
	Nhalosbin[0] = 1;
	zbins[0] = 0;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(k)
#endif
	for(k=1;k<Nzbins-1;++k){
		Nhalosbin[k] = cosmo->haloNumberInBufferedCone(min_mass*mass_scale,zbins[k],zsource,fieldofview*pow(pi/180,2),field_buffer,mass_func_type,mass_func_PL_slope)/aveNhalos;
	}
	zbins[Nzbins-1] = zsource;
	Nhalosbin[Nzbins-1] = 0.0;

	Nhalos = (long)(poidev(float(aveNhalos), seed) );

	std::vector<double> halo_zs_vec;
	std::vector<double *> halo_pos_vec;

	// assign redsshifts to field_halos and sort them

	for(i=0;i < Nhalos;++i){
		halo_zs_vec.push_back(InterpolateYvec(Nhalosbin,zbins,ran2(seed)));
	}

	std::sort(halo_zs_vec.begin(),halo_zs_vec.end());

	assert(halo_zs_vec[0] < halo_zs_vec[1]);
	assert(halo_zs_vec[0] < halo_zs_vec[Nhalos-1]);

	// fill the log(mass) vector
	Logm.resize(Nmassbin);
	Nhalosbin.resize(Nmassbin);
	fill_linear(Logm,Nmassbin,log10(min_mass*mass_scale),MaxLogm);

	double *pos;
	int j = 0;
	k2 = 0;
	std::vector<double>::iterator it1,it2;
	for(np=0,mass_max=0;np<NZSamples;np++){

		z1 = np*zsource/(NZSamples);
		z2 = (np+1)*zsource/(NZSamples);

		it1 = std::lower_bound(halo_zs_vec.begin(),halo_zs_vec.end(),z1);
		it2 = std::lower_bound(halo_zs_vec.begin(),halo_zs_vec.end(),z2);

		k1 = it1-halo_zs_vec.begin();
		k2 = it2-halo_zs_vec.begin();

		Nhaloestot = cosmo->haloNumberInBufferedCone(pow(10,Logm[0]),z1,z2,fieldofview*pow(pi/180,2),field_buffer,mass_func_type,mass_func_PL_slope);

		Nhalosbin[0] = 1;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(k)
#endif
		for(k=1;k<Nmassbin-1;k++){
			// cumulative number density in one square degree
			Nhalosbin[k] = cosmo->haloNumberInBufferedCone(pow(10,Logm[k]),z1,z2,fieldofview*pow(pi/180,2),field_buffer,mass_func_type,mass_func_PL_slope)
					/Nhaloestot;
		}
		Nhalosbin[Nmassbin-1] = 0;

		for(i = k1; i < k2; i++){
			double Ds = cosmo->angDist(0,halo_zs_vec[i]);

			maxr = pi*sqrt(fieldofview/pi)/180. + field_buffer/Ds; // fov is a circle
			rr = maxr*sqrt(ran2(seed));

			assert(rr == rr);

			pos = new double[2];

			theta = 2*pi*ran2(seed);

			pos[0] = rr*cos(theta)*Ds;
			pos[1] = rr*sin(theta)*Ds;

			switch(int_prof_type){
			case null_lens:
				ERROR_MESSAGE();
				std::cout << "int_prof_type is null!!!!" << std::endl;
				break;
			case nfw_lens:
				main_halos.push_back(new NFWLensHalo);
				break;
			case pnfw_lens:
				main_halos.push_back(new PseudoNFWLensHalo);
				break;
			case pl_lens:
				main_halos.push_back(new PowerLawLensHalo);
				break;
			case nsie_lens:
				main_halos.push_back(new SimpleNSIELensHalo);
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
				main_halos.push_back(new DummyLensHalo);
				break;
			}

			float mass = pow(10,InterpolateYvec(Nhalosbin,Logm,ran2 (seed)));

			halo_calc->reset(mass,halo_zs_vec[i]);

			float Rmax = halo_calc->getRvir();
			float rscale = Rmax/halo_calc->getConcentration(0);

			field_halos[j]->setZlens(halo_zs_vec[i]);
			if(flag_galaxy_subhalo){
				if(galaxy_mass_fraction > 1.0) galaxy_mass_fraction = 1;
				field_halos[j]->initFromMassFunc(mass*(1-galaxy_mass_fraction),mass_scale,Rmax,rscale,halo_slope,seed);
			}
			else{
				field_halos[j]->initFromMassFunc(mass,mass_scale,Rmax,rscale,halo_slope,seed);
			}

			if(mass > mass_max) {
				mass_max = mass;
				j_max = i;
				pos_max[0] = pos[0];
				pos_max[1] = pos[1];
				z_max = halo_zs_vec[i];
			}

			halo_pos_vec.push_back(pos);

			++j;

			if(flag_galaxy_subhalo){
				switch(int_prof_gal_type){
				case null_gal:
					ERROR_MESSAGE();
					std::cout << "flag_galaxy_subhalo is true, but int_prof_gal_type is null!!!!" << std::endl;
					break;
				case nsie_gal:
					field_halos.push_back(new SimpleNSIELensHalo);
					break;
				}

				field_halos[j]->setZlens(halo_zs_vec[i]);
				field_halos[j]->initFromMassFunc(mass*galaxy_mass_fraction,mass_scale,Rmax,rscale,halo_slope,seed);

				halo_pos_vec.push_back(pos);

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

	std::cout << "leaving MultiLens::createHaloData_buffered()" << std::endl;
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
void Lens::readInputSimFile(CosmoHndl cosmo){

	double ra,dec,z,vmax,vdisp,r_halfmass;
	unsigned long i,j;
	unsigned long haloid,idd,np;
	double mo=7.3113e10,M1=2.8575e10,gam1=7.17,gam2=0.201,be=0.557;

	double rmax=0,rtmp=0;

	std::ifstream file_in(input_sim_file.c_str());
	if(!file_in){
		std::cout << "Can't open file " << input_sim_file << std::endl;
		exit(1);
	}

	// skip through header information in data file
	i=0;
	while(file_in.peek() == '#'){
		file_in.ignore(10000,'\n');
		++i;
	}
	std::cout << "skipped "<< i << " comment lines in file " << input_sim_file << std::endl;

	std::vector<double> halo_zs_vec;
	std::vector<double *> halo_pos_vec;

	// read in data
	int j_max;
	double mass_max=0,R_max=0,V_max=0,minmass=1e30;
	double *theta;
	int ncolumns = 9;

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
			if(l == 0 || l == 1 || l == 5){
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

		/// pos in radians
		theta = new double[2];
		/// pos in physical radians
		double Ds = cosmo->angDist(0,z);
		theta[0] = ra*pi/180.*Ds;
		theta[1] = dec*pi/180.*Ds;

        if(rmax < (rtmp = theta[0]*theta[0]+theta[1]*theta[1])) rmax = rtmp;

		if(np > 0.0 && vdisp > 0.0 && z <= zsource){

			switch(int_prof_type){
			case null_lens:
				ERROR_MESSAGE();
				std::cout << "int_prof_type is null!!!!" << std::endl;
				break;
			case nfw_lens:
				main_halos.push_back(new NFWLensHalo);
				break;
			case pnfw_lens:
				main_halos.push_back(new PseudoNFWLensHalo);
				break;
				ERROR_MESSAGE();
				std::cout << "PseudoNFW not supported." << std::endl;
				break;
			case pl_lens:
				ERROR_MESSAGE();
				std::cout << "PowerLaw not supported." << std::endl;
				break;
			case nsie_lens:
				main_halos.push_back(new SimpleNSIELensHalo);
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
				main_halos.push_back(new DummyLensHalo);
				ERROR_MESSAGE();
				std::cout << "Why would you wand dummy file halos???" << std::endl;
				break;
			}

			float mass = np*8.6e8/cosmo->gethubble();

			field_halos[j]->setZlens(z);
			if(flag_galaxy_subhalo){
				galaxy_mass_fraction = 2*mo*pow(mass/M1,gam1)
				  /pow(1+pow(mass/M1,be),(gam1-gam2)/be)/mass;
				if(galaxy_mass_fraction > 1.0) galaxy_mass_fraction = 1;

				field_halos[j]->initFromFile(mass*(1-galaxy_mass_fraction),mass_scale,seed,vmax,r_halfmass*cosmo->gethubble());
			}
			else{
				field_halos[j]->initFromFile(mass,mass_scale,seed,vmax,r_halfmass*cosmo->gethubble());
			}


			if(field_halos[j]->get_Rmax() > R_max) R_max = field_halos[j]->get_Rmax();
			if(vdisp > V_max) V_max = vdisp;

			halo_zs_vec.push_back(z);
			halo_pos_vec.push_back(theta);

			if(mass > mass_max) {
				mass_max = mass;
				j_max = j;
			}
			if(mass < minmass) {
				minmass = mass;
			}

			++j;

			if(flag_galaxy_subhalo){
				switch(int_prof_gal_type){
				case null_gal:
					ERROR_MESSAGE();
					std::cout << "flag_galaxy_subhalo is true, but int_prof_gal_type is null!!!!" << std::endl;
					break;
				case nsie_gal:
					field_halos.push_back(new SimpleNSIELensHalo);
					break;
				}

				field_halos[j]->setZlens(z);
				field_halos[j]->initFromFile(mass*galaxy_mass_fraction,mass_scale,seed,vmax,r_halfmass*cosmo->gethubble());

				halo_pos_vec.push_back(theta);

				++j;
			}

		}
	}
	file_in.close();
	std::cout << field_halos.size() << " halos read in."<< std::endl
			<< "Max input mass = " << mass_max << "  R max = " << R_max << "  V max = " << V_max << std::endl;

	Nhalos = field_halos.size();

	/// setting the minimum halo mass in the simulation
	min_mass = minmass;
	if(field_buffer > 0.0){
		std::cout << "Overiding field_buffer to make it 0 because halos are read in." << std::endl;
		field_buffer = 0.0;
	}

	halo_pos = Utilities::PosTypeMatrix(Nhalos,3);

	for(i=0;i<Nhalos;++i){
		halo_pos[i] = halo_pos_vec[i];
	}

	std::cout << "sorting in MultiLens::readInputSimFile()" << std::endl;
	// sort the field_halos by readshift
	Lens::quicksort(field_halos.data(),halo_pos,field_halos.size());

	std::cout << "leaving MultiLens::readInputSimFile()" << std::endl;

	read_sim_file = true;
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
		CosmoHndl cosmo           /// cosmology
		,double z                 /// redshift of implanted source
		,bool nearest             /** If true, set the source plane to the nearest (in coordinate distance)
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

	// j is the index of the next plane at higher redshift, This plane will be temporarily replaced and used as a source plane

	double Ds = cosmo->coorDist(0,z);
	locateD(Dl.data(),Nplanes,Ds,&j);
	assert(j <= Nplanes && j >=0);

	if(j >= Nplanes-1){
	  j--;
	}
	else if(j > 0){
		std::map<double,double>::iterator ind;

		ind = coorDist_table.lower_bound((Dl[j]-0.5*dDl[j]));
		double z1 = ind->second;

		if(nearest) j = (z>=z1) ? j : j-1;  // reset j to the nearest plane
	}

	if(nearest && (j < Nplanes-1) ){
		zs_implant = plane_redshifts[j];
		Ds_implant = Dl[j];
		if(j > 0) dDs_implant = dDl[j];
		else  dDs_implant = Ds_implant;
	}else{
		// if nearest==false or the source is at higher redshift than the last plane use the real redshift
		Ds_implant = Ds;
		zs_implant = z;
		if(j > 0) dDs_implant = cosmo->coorDist(plane_redshifts[j-1],z);
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
