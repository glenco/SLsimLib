/*
 * halo_model.cpp
 *
 *  Created on: Jan 14, 2012
 *      Author: mpetkova, bmetcalf
 */

#include "slsimlib.h"
#include <algorithm>
#include "lens_halos.h"
#include <iomanip>      // std::setprecision

using namespace std;

#define MIN_PLANE_DIST 1E-8

namespace
{
	class lens_halo_less
	{
	public:
		lens_halo_less(const COSMOLOGY& c) : cosmo(c) {}
		
		bool operator()(const LensHalo* a, const LensHalo* b)
		{
			// compare sizes and check that b is not in eps around a
			return (a->getZlens() < b->getZlens()) && fabs(cosmo.coorDist(a->getZlens(), b->getZlens())) > MIN_PLANE_DIST;
		}
		
	private:
		const COSMOLOGY& cosmo;
	};
}

/**
 * \brief Creates an empty lens. Main halos and field halos need to be inserted by hand from the user.
 */
Lens::Lens(long* my_seed,PosType z_source, CosmoParamSet cosmoset,bool verbose)
: seed(my_seed), cosmo(cosmoset),zsource(z_source), central_point_sphere(1,0,0), inv_ang_screening_scale(0)
{
  init_seed = 0;
  
	if((cosmo.getOmega_matter() + cosmo.getOmega_lambda()) != 1.0)
	{
		printf("ERROR: Lens can only handle flat universes at present. Must change cosmology.\n");
		exit(1);
	}
	
	read_sim_file = false;
	
  charge = 4*pi*Grav;
	if(verbose) std::cout << "charge: " << charge << std::endl;
	
	// initially let source be the one inputed from parameter file
	index_of_new_sourceplane = -1;
	toggle_source_plane = false;
  flag_switch_deflection_off = false;
  flag_switch_lensing_off = false;
  
  //charge = cosmo.angDist(zsource)/cosmo.angDist(0.3)/cosmo.angDist(0.3,zsource);
  //charge = 4*pi/cosmo.angDist(0.3);
  PosType ztmp = zsource;
  combinePlanes(true);
  if(zsource != ztmp) ResetSourcePlane(ztmp,false);
  std::cout << "number of field halos : " << field_halos.size() << std::endl;
}

Lens::Lens(long* my_seed,PosType z_source, const COSMOLOGY &cosmoset,bool verbose)
: seed(my_seed), cosmo(cosmoset),zsource(z_source), central_point_sphere(1,0,0), inv_ang_screening_scale(0)
{
  init_seed = 0;
  
  if((cosmo.getOmega_matter() + cosmo.getOmega_lambda()) != 1.0)
  {
    printf("ERROR: Lens can only handle flat universes at present. Must change cosmology.\n");
    exit(1);
  }
  
  read_sim_file = false;
  
  charge = 4*pi*Grav;
  if(verbose) std::cout << "charge: " << charge << std::endl;
  
  // initially let source be the one inputed from parameter file
  index_of_new_sourceplane = -1;
  toggle_source_plane = false;
  flag_switch_deflection_off = false;
  flag_switch_lensing_off = false;
  
  //charge = cosmo.angDist(zsource)/cosmo.angDist(0.3)/cosmo.angDist(0.3,zsource);
  //charge = 4*pi/cosmo.angDist(0.3);
  PosType ztmp = zsource;
  combinePlanes(true);
  if(zsource != ztmp) ResetSourcePlane(ztmp,false);
  std::cout << "number of field halos :" << field_halos.size() << std::endl;
}

/**
 * \ingroup Constructor
 * \brief allocates space for the halo trees and the inout lens, if there is any
 */
Lens::Lens(InputParams& params, long* my_seed, CosmoParamSet cosmoset, bool verbose)
: seed(my_seed), cosmo(cosmoset), central_point_sphere(1,0,0), inv_ang_screening_scale(0)
{
  
  //init_params = params;
  init_seed = *my_seed;
  
	readCosmology(params);
	
	if((cosmo.getOmega_matter() + cosmo.getOmega_lambda()) != 1.0)
	{
		printf("ERROR: Lens can only handle flat universes at present. Must change cosmology.\n");
		exit(1);
	}
	
	assignParams(params);
	
	read_sim_file = false;
	
	charge = 4*pi*Grav;
	if(verbose) std::cout << "charge: " << charge << std::endl;
	
	// initially let source be the one inputed from parameter file
	index_of_new_sourceplane = -1;
	toggle_source_plane = false;

  
  if(flag_switch_field_off == false) {
    std::cout << "Nzbins = " << Nzbins << std::endl ;
    
    // Resizing the "number of Halos" binning table :
    zbins.resize(Nzbins) ;
    NhalosbinZ.resize(Nzbins) ;
    Nhaloestot_Tab.resize (NZSamples);
    
    // Initialising the "number of Halos" binning table :
    for(int k=0 ; k<Nzbins ; k++)
    {
      zbins[k] = 0. ;
      NhalosbinZ[k] = 0. ;
    }
    aveNhalosField = 0. ;
    for(int k=0 ; k<NZSamples ; k++) Nhaloestot_Tab[k] = 0. ;
    
    // Computing the number of halos per bins :
    if(sim_input_flag){
      // Do Nothing ! No step is necessary here !
    }
    else {
      if(field_buffer == 0.0){
        field_buffer = pow(3.0e14/800/pi/cosmo.rho_crit(0),1.0/3.);
        std::cout << "    Resetting field buffer to " << field_buffer << " Mpc." << std::endl;
      }
      // Compute the distribution variables :
      ComputeHalosDistributionVariables();
      // for (int i=0; i<Nzbins; i++) std::cout << NhalosbinZ[i] << " " ;
      // std::cout << std::endl ;
    }
    
    std::cout << "Lens : field_Nplanes_original = " << field_Nplanes_original << std::endl ;
  }
  
  // set up the lens contents :
  PosType ztmp = zsource;
	buildPlanes(params, verbose);
  if(zsource != ztmp) ResetSourcePlane(ztmp,false);
  std::cout << "number of field halos : " << field_halos.size() << std::endl;

}
/**
 * \ingroup Constructor
 * \brief allocates space for the halo trees and the inout lens, if there is any
 */
Lens::Lens(InputParams& params, long* my_seed, const COSMOLOGY &cosmoset, bool verbose)
: seed(my_seed), cosmo(cosmoset), central_point_sphere(1,0,0), inv_ang_screening_scale(0)
{
  
  //init_params = params;
  init_seed = *my_seed;
  
  readCosmology(params);
  
  if((cosmo.getOmega_matter() + cosmo.getOmega_lambda()) != 1.0)
  {
    printf("ERROR: Lens can only handle flat universes at present. Must change cosmology.\n");
    exit(1);
  }
  
  assignParams(params);
  
  read_sim_file = false;
  
  charge = 4*pi*Grav;
  if(verbose) std::cout << "charge: " << charge << std::endl;
  
  // initially let source be the one inputed from parameter file
  index_of_new_sourceplane = -1;
  toggle_source_plane = false;
  
  
  if(flag_switch_field_off == false) {
    std::cout << "Nzbins = " << Nzbins << std::endl ;
    
    // Resizing the "number of Halos" binning table :
    zbins.resize(Nzbins) ;
    NhalosbinZ.resize(Nzbins) ;
    Nhaloestot_Tab.resize (NZSamples);
    
    // Initialising the "number of Halos" binning table :
    for(int k=0 ; k<Nzbins ; k++)
    {
      zbins[k] = 0. ;
      NhalosbinZ[k] = 0. ;
    }
    aveNhalosField = 0. ;
    for(int k=0 ; k<NZSamples ; k++) Nhaloestot_Tab[k] = 0. ;
    
    // Computing the number of halos per bins :
    if(sim_input_flag){
      // Do Nothing ! No step is necessary here !
    }
    else {
      if(field_buffer == 0.0){
        field_buffer = pow(3.0e14/800/pi/cosmo.rho_crit(0),1.0/3.);
        std::cout << "    Resetting field buffer to " << field_buffer << " Mpc." << std::endl;
      }
      // Compute the distribution variables :
      ComputeHalosDistributionVariables();
      // for (int i=0; i<Nzbins; i++) std::cout << NhalosbinZ[i] << " " ;
      // std::cout << std::endl ;
    }
    
    std::cout << "Lens : field_Nplanes_original = " << field_Nplanes_original << std::endl ;
  }
  
  // set up the lens contents :
  PosType ztmp = zsource;
  buildPlanes(params, verbose);
  if(zsource != ztmp) ResetSourcePlane(ztmp,false);
  std::cout << "number of field halos :" << field_halos.size() << std::endl;
  
}

Lens::~Lens()
{
	Utilities::delete_container(lensing_planes);
  
	//Utilities::free_PosTypeMatrix(halo_pos, field_halos.size(), 3);
  
	Utilities::delete_container(main_halos_created);
	Utilities::delete_container(field_halos);
  Utilities::delete_container(substructure.halos);
}

/// read in Cosmological Parameters
void Lens::readCosmology(InputParams& params)
{
	PosType tmp;
	if(params.get("Omega_matter", tmp)) cosmo.setOmega_matter(tmp, true);
	if(params.get("Omega_lambda", tmp)) cosmo.setOmega_lambda(tmp);
	if(params.get("Omega_baryon", tmp)) cosmo.setOmega_baryon(tmp);
	if(params.get("Omega_neutrino", tmp)) cosmo.setOmega_neutrino(tmp);
	if(params.get("hubble", tmp)) cosmo.sethubble(tmp);
	if(params.get("sigma_8", tmp)) cosmo.power_normalize(tmp);
}

/// Retrieve input parameters for construction.
void Lens::assignParams(InputParams& params,bool verbose)
{
	if(!params.get("main_halo_on",flag_switch_main_halo_on))
	{
		ERROR_MESSAGE();
		std::cout << "parameter main_halo_on needs to be set in the parameter file " << params.filename() << endl;
		exit(0);
	}
	
	if(flag_switch_main_halo_on)
	{
		if(!params.get("main_DM_halo_type",main_halo_type))
		{
			ERROR_MESSAGE();
			std::cout << "parameter main_DM_halo_type needs to be set in the parameter file " << params.filename() << endl;
			exit(0);
		}
		if(!params.get("main_galaxy_halo_type",main_galaxy_halo_type))
		{
			main_galaxy_halo_type = null_gal;
		}
	}
	else
	{
		main_halo_type = null_lens;
		main_galaxy_halo_type = null_gal;
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
			if(!params.get("field_Nplanes",field_Nplanes_original))
			{
				ERROR_MESSAGE();
				std::cout << "parameter field_Nplanes needs to be set in the parameter file " << params.filename() << endl;
				exit(0);
			}
      field_Nplanes_current = field_Nplanes_original;
			
			if(!params.get("field_internal_profile",field_int_prof_type))
			{
				ERROR_MESSAGE();
				std::cout << "parameter field_internal_profile needs to be set in the parameter file " << params.filename() << endl;
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
        switch ( field_int_prof_gal_type )
        {
          case null_gal:
            break;
          case nsie_gal:
            break;
          case pl_gal:
            break;
          case hern_gal:
            break;
          case jaffe_gal:
            break;
          default :
            ERROR_MESSAGE();
            std::cout<<"field_internal_profile_galaxy must be Null, NSIE,PowerLaw, Hernquist or Jaffe"<<std::endl;
            exit(0);
        }
			}
      
      if (field_int_prof_gal_type == pl_gal){
        if(!params.get("field_internal_profile_galaxy_slope",field_int_prof_gal_slope)){
          ERROR_MESSAGE();
          std::cout << "field_internal_profile_galaxy_slope must be specified for field_internal_profile_galaxy PowerLaw in the parameter file "
          << params.filename() << endl;
          exit(0);
        }
			}
      

			if(!params.get("field_mass_func_alpha",mass_func_PL_slope))
				mass_func_PL_slope = 1./6.;
			if(!params.get("field_prof_internal_slope_pl",field_prof_internal_slope) && field_int_prof_type == pl_lens)
				field_prof_internal_slope = -1.0;
			if(!params.get("field_prof_internal_slope_pnfw",field_prof_internal_slope) && field_int_prof_type == pnfw_lens)
				field_prof_internal_slope = 2.0;
      
      if(params.get("field_input_simulation_file",field_input_sim_file)){
        ERROR_MESSAGE();
        std::cout << "parameter field_input_simulation_file should be changed to field_input_simulation_path "
        << params.filename() << endl;
        exit(0);
      }
      
			if(!params.get("field_input_simulation_path",field_input_sim_file))
			{
				// No simulation input file provided
				sim_input_flag = false;
				
				if(!params.get("field_mass_func_type",field_mass_func_type))
				{
					ERROR_MESSAGE();
					std::cout << "parameter field_mass_func_type needs to be set in the parameter file " << params.filename() << endl;
					exit(0);
				}
				
				if(!params.get("field_min_mass",field_min_mass))
				{
					ERROR_MESSAGE();
					std::cout << "parameter field_min_mass needs to be set in the parameter file " << params.filename() << endl;
					exit(0);
				}
				
				if(!params.get("field_buffer",field_buffer))
				{
					field_buffer = 0.0;
					std::cout << "default field buffer of 0 Mpc is being used for now." << endl;
				}
        if(!params.get("field_fov",fieldofview))
        {
          ERROR_MESSAGE();
          std::cout << "parameter field_fov needs to be set in the parameter file " << params.filename() << endl;
          exit(0);
        }
			}
			else
			{
        field_min_mass = 0.0;
				sim_input_flag = true;
        
        if(!params.get("field_input_simulation_format",field_input_sim_format)){
          ERROR_MESSAGE();
					std::cout << "parameter field_input_simulation_format needs to be set in the parameter file " << params.filename() << endl;
					exit(0);
        }
        
        if(!params.get("field_input_simulation_center_RA",central_point_sphere.phi)){
          central_point_sphere.phi = 0.0;
        }else{
          central_point_sphere.phi *= degreesTOradians;
        }
        if(!params.get("field_input_simulation_center_DEC",central_point_sphere.theta)){
          central_point_sphere.theta = 0.0;
        }else{
          central_point_sphere.theta *= degreesTOradians;
        }
        if(!params.get("field_input_simulation_radius",sim_angular_radius)){
          sim_angular_radius = 0.0;
        }else{
          sim_angular_radius *= degreesTOradians;
        }
			}
      if(!params.get("field_fov",fieldofview));
		}
		else
		{
			// no field
			field_Nplanes_original = field_Nplanes_current = 0;
			fieldofview = 0;
		}
	}
	
	// read Pixeliz parameters if necessary
  if(!params.get("pixelmaps_on",pixel_map_on)) pixel_map_on = 0;
	if(pixel_map_on)
	{
		if(!params.get("pixelmaps_input_file", pixel_map_input_file))
		{
			ERROR_MESSAGE();
			std::cout << "parameter pixelmaps_input_file needs to be set in the parameter file " << params.filename() << endl;
			exit(1);
		}
    if(flag_switch_field_off == false && pixel_map_input_file.find(".fits") == pixel_map_input_file.npos){
      std::cerr << "Warning: Are you sure you want to put both field lenses and pixel maps in the light-cone?" << endl;
    }
    if(!params.get("pixelmaps_padding_factor",pixel_map_zeropad)){
      pixel_map_zeropad = 4;
    }
    if(!params.get("pixelmaps_zeromean",pixel_map_zeromean)){
      pixel_map_zeromean = true;
    }
	}
	
	if(!params.get("z_source",zsource))
	{
		ERROR_MESSAGE();
		std::cout << "parameter z_source needs to be set in the parameter file " << params.filename() << endl;
		exit(0);
	}
	
	if(!params.get("deflection_off",flag_switch_deflection_off))
		flag_switch_deflection_off = false;
  
  if(!params.get("lensing_off",flag_switch_lensing_off))
		flag_switch_lensing_off = false;
  
	// Some checks for valid parameters
	if(flag_switch_field_off == false && field_Nplanes_original == 0 )
	{
		ERROR_MESSAGE();
		std::cout << "Do you want to run _with_ field halos, but with _without_ field planes? Change field_Nplanes to a bigger number!" << endl;
		exit(1);
	}
	
	if(flag_switch_main_halo_on == false && flag_switch_field_off == true && !pixel_map_on)
	{
		ERROR_MESSAGE();
		std::cout << "Do you want an empty simulation? Set main_halo_on to true for a main lens, or field_off to false for field lenses." << endl;
		exit(1);
	}
	
	if(!flag_switch_field_off)
	{
		if(field_int_prof_type == pl_lens && field_prof_internal_slope >= 0)
		{
			ERROR_MESSAGE();
			std::cout << "Power Law internal slope >=0 not possible." << endl;
			exit(1);
		}
		
		if(field_int_prof_type == pnfw_lens && field_prof_internal_slope <= 0)
		{
			ERROR_MESSAGE();
			std::cout << "Pseudo NFW internal slope <=0 not possible." << endl;
			exit(1);
		}
		
		if(field_int_prof_type == pnfw_lens && (field_prof_internal_slope / floor(field_prof_internal_slope) > 1.0))
		{
			ERROR_MESSAGE();
			std::cout << "Pseudo NFW internal slope needs to be a whole number." << endl;
			exit(1);
		}
		
		if(field_input_sim_file.size() < 1 && field_int_prof_type == nsie_lens)
		{
			ERROR_MESSAGE();
			std::cout << "The NSIE internal profile works only for Millennium DM simulations for now." << endl;
			std::cout << "Set field_input_simulation_path in sample_paramfile." << endl;
			exit(1);
		}
	}
	
	// convert to square degrees
	fieldofview /= 3600. * 3600.;
	
	if(verbose) printMultiLens();
}


// get the redshift of a plane
PosType Lens::getFieldPlaneZ(int n)
{
  return field_plane_redshifts[n] ;
}


void Lens::resetFieldNplanes(std::size_t Np, bool verbose)
{
	Utilities::delete_container(field_planes);
	
	field_Nplanes_original = field_Nplanes_current = Np;
	
	field_plane_redshifts.clear();
	field_Dl.clear();
	
	setupFieldPlanes();
	createFieldPlanes(verbose);
	
	combinePlanes(verbose);
  
  field_plane_redshifts_original = field_plane_redshifts;
  field_Dl_original = field_Dl;
}

void Lens::resetFieldHalos(bool verbose)
{
  Utilities::delete_container(field_halos);
	Utilities::delete_container(field_planes);
  field_Nplanes_current = field_Nplanes_original;
  field_plane_redshifts = field_plane_redshifts_original;
  field_Dl = field_Dl_original;
  
	if(sim_input_flag){
		if(read_sim_file == false){
      if(field_input_sim_format == MillenniumObs) readInputSimFileMillennium(verbose);
      if(field_input_sim_format == MultiDarkHalos) readInputSimFileMultiDarkHalos(verbose);
      if(field_input_sim_format == ObservedData) readInputSimFileObservedGalaxies(verbose);
    }
	}
	else{
		createFieldHalos(verbose);
	}
  
  // set up the lens contents :
	createFieldPlanes(verbose);
	combinePlanes(verbose);
  
  if(WasInsertSubStructuresCalled == YES){
    WasInsertSubStructuresCalled = MAYBE ;
  }
  std::cout << "number of field halos : " << field_halos.size() << std::endl;

}

void Lens::printMultiLens(){
	std::cout << endl << "MAIN HALOS" << endl;
	std::cout << "Main lens profile type:" << endl;
	switch(main_halo_type){
    case null_lens:
      std::cout << "no lens" << endl;
      break;
    case nfw_lens:
      std::cout << "NFW lens" << endl;
      break;
    case pnfw_lens:
      std::cout << "PseudoNFW lens" << endl;
      std::cout << "slope: " << field_prof_internal_slope << endl;
      break;
    case pl_lens:
      std::cout << "PowerLaw lens" << endl;
      std::cout << "slope: " << field_prof_internal_slope << endl;
      break;
    case nsie_lens:
      std::cout << "NSIE lens" << endl;
      break;
    case ana_lens:
      std::cout << "AnaNSIE lens" << endl;
      break;
    case uni_lens:
      std::cout << "UniNSIE lens" << endl;
      break;
    case moka_lens:
      std::cout << "MOKA lens" << endl;
      break;
    case dummy_lens:
      std::cout << "Dummy lens" << endl;
      break;
    case hern_lens:
      std::cout << "Hernquist lens" << endl;
      break;
    case jaffe_lens:
      std::cout << "Jaffe lens" << endl;
      break;
	}
  
  if(pixel_map_on) std::cout << "PixelDMap lens" << endl;
    
	std::cout << endl << "Main galaxies profile type:" << endl;
	switch(main_galaxy_halo_type){
    case null_gal:
      std::cout << "no galaxy" << endl;
      break;
    case nsie_gal:
      std::cout << "NSIE galaxy" << endl;
      break;
    case pl_gal:
      std::cout << "PowerLaw galaxy" << endl;
      break;
    case hern_gal:
      std::cout << "Hernquist galaxy" << endl;
      break;
    case jaffe_gal:
      std::cout << "Jaffe galaxy" << endl;
      break;
	}
  
	if(flag_switch_field_off == false){
    
		std::cout << "field of view " << fieldofview << endl;
    
		std::cout << endl << "FIELD HALOS" << endl;
    
    std::cout << "field Nplanes original" << field_Nplanes_original << endl;
    std::cout << "field Nplanes current" << field_Nplanes_current << endl;
    
		std::cout << "min mass " << field_min_mass << endl;
		std::cout << "Mass function type: "<< endl;
    
		switch(field_mass_func_type){
      case PS:
        std::cout << "  Press-Schechter mass function " << endl;
        break;
      case ST:
        std::cout << "  Sheth-Tormen mass function " << endl;
        break;
      case PL:
        std::cout << "  Power law mass function " << endl;
        std::cout << "  slope: " << mass_func_PL_slope << endl;
        break;
		}
    
		std::cout << endl << "Field halos profile type:" << endl;
		switch(field_int_prof_type)
		{
			case null_lens:
				std::cout << "no field type" << endl;
				break;
			case nfw_lens:
				std::cout << "NFW field type" << endl;
				break;
			case pnfw_lens:
				std::cout << "PseudoNFW field type" << endl;
				std::cout << "slope: " << field_prof_internal_slope << endl;
				break;
			case pl_lens:
				std::cout << "PowerLaw field type" << endl;
				std::cout << "slope: " << field_prof_internal_slope << endl;
				break;
			case nsie_lens:
				std::cout << "NSIE field type" << endl;
				break;
			case ana_lens:
				std::cout << "AnaNSIE field type" << endl;
				break;
			case uni_lens:
				std::cout << "UniNSIE field type" << endl;
				break;
      case moka_lens:
        std::cout << "MOKA field type" << endl;
        break;
			case dummy_lens:
				std::cout << "Dummy field type" << endl;
				break;
			case hern_lens:
				std::cout << "Hernquist field type" << endl;
				break;
			case jaffe_lens:
				std::cout << "Jaffe field type" << endl;
				break;
		}
    
    if(pixel_map_on) std::cout << "PixelDMap lens" << endl;

		std::cout << endl << "Field galaxies profile type:" << endl;
		switch(field_int_prof_gal_type){
      case null_gal:
        std::cout << "no field galaxy type" << endl;
        break;
      case nsie_gal:
        std::cout << "NSIE field galaxy type" << endl;
        break;
      case pl_gal:
        std::cout << "PowerLaw field galaxy type" << endl;
        break;
      case hern_gal:
        std::cout << "Hernquist field galaxy type" << endl;
        break;
      case jaffe_gal:
        std::cout << "Jaffe field galaxy type" << endl;
        break;

		}
	}
  
	std::cout << endl;
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
void Lens::createFieldPlanes(bool verbose)
{
	if(verbose) std::cout << "Lens::createFieldPlanes zsource = " << zsource << std::endl;
    
	assert(field_plane_redshifts.size() == field_Nplanes_original);
	
	// the bounds for sorting field halos onto redshifts
	PosType z1 = 0, z2 = 0;
	std::size_t k1 = 0, k2 = 0;
	
	// go through planes
	for(std::size_t i = 0; i < field_Nplanes_original; ++i)
	{
		assert(field_plane_redshifts[i] > 0);
		assert(field_Dl[i] > 0);
		
		// previous upper bound is now lower bound
		z1 = z2;
		k1 = k2;
		
		// find upper bound
		if(i == field_Nplanes_original-1)
		{
			z2 = zsource;
			k2 = field_halos.size();
		}
		else
		{
			z2 = cosmo.invCoorDist(0.5*(field_Dl[i] + field_Dl[i+1]));
			k2 = Utilities::lower_bound<LensHalo>(field_halos, z2);
		}
		
		/*
		 * finding the average mass surface density in field_halos
		 */
    
    PosType sigma_back;
    if(sigma_back_Tab.size() < field_Nplanes_original){
      sigma_back =
        cosmo.haloMassInBufferedCone(field_min_mass,z1,z2,fieldofview*pow(pi/180,2),field_buffer,field_mass_func_type,mass_func_PL_slope)/(pi*pow(sqrt(fieldofview/pi)*pi*field_Dl[i]/180/(1+field_plane_redshifts[i]) + field_buffer,2));
      sigma_back_Tab.push_back(sigma_back);
    }else{
      sigma_back = sigma_back_Tab[i];
    }
		
		PosType sb=0.0;
		//PosType max_r = 0,tmp;
    
		for(std::size_t j = k1; j < k2; ++j)
		{
			sb += field_halos[j]->get_mass();
			
			//** test lines********
			//if(max_r < (tmp = halo_pos[j][0]*halo_pos[j][0] + halo_pos[j][1]*halo_pos[j][1])) max_r = tmp;
			
      PosType tmp[2];
      
			// convert to proper distance on the lens plane
      field_halos[j]->getX(tmp);
      field_halos[j]->setX(tmp[0]*field_Dl[i]/(1+field_plane_redshifts[i])
                           ,tmp[1]*field_Dl[i]/(1+field_plane_redshifts[i]));
			//halo_pos[j][0] *= field_Dl[i]/(1+field_plane_redshifts[i]);
			//halo_pos[j][1] *= field_Dl[i]/(1+field_plane_redshifts[i]);
		}
		
		//max_r=sqrt(max_r);
		
		sb /= (pi*pow(sqrt(fieldofview/pi)*pi*field_Dl[i]/180/(1+field_plane_redshifts[i]) + field_buffer,2));
		
    
    assert(sb == sb);
    
    if(verbose) std::cout << "sigma_back from mass function " << sigma_back
      << " from sum of halos " << sb << " " << sb/sigma_back - 1 << std::endl;
		if(sim_input_flag) sigma_back = sb;
		//sigma_back = sb;
		/*
		 * create the lensing plane
		 */
		
		if(verbose) std::cout << "  Building lensing plane " << i << " number of halos: "
      << k2-k1 << std::endl;
		
    PosType tmp = inv_ang_screening_scale*(1+field_plane_redshifts[i])/field_Dl[i];
    
    if(tmp > 1/2.) tmp = 1/2.;  // TODO: Try to remove this arbitrary minimum distance
		field_planes.push_back(new LensPlaneTree(&field_halos[k1], k2-k1, sigma_back,tmp));
		//field_planes.push_back(new LensPlaneTree(&halo_pos[k1], &field_halos[k1], k2-k1, sigma_back) );
	}
  
	assert(field_planes.size() == field_Nplanes_original);
}



/// * INSERT SUB STRUCTURE * ///


void Lens::insertSubstructures(PosType Rregion,           // in radians
                               PosType center[],
                               PosType NumberDensity,     // in number / unit^2
                               PosType Mass_min,          // in M_sun
                               PosType Mass_max,          // in M_sun
                               PosType redshift,
                               PosType alpha,             // Careful ! alpha is opposite sign wrt Metcalf, Amara 2011.
                               PosType density_contrast,  // dimensionless
                               bool verbose
                               )
{
  
  substructure.alpha = alpha;
  substructure.center.x[0] = center[0];
  substructure.center.x[1] = center[1];
  substructure.Mmax = Mass_max;
  substructure.Mmin = Mass_min;
  substructure.Ndensity = NumberDensity;
  substructure.rho_tidal = density_contrast;
  substructure.Rregion = Rregion;
  substructure.redshift = redshift;
  
  if(WasInsertSubStructuresCalled == YES || WasInsertSubStructuresCalled == MAYBE)
  {
    ERROR_MESSAGE();
    cout << "Lens::insertSubStructures() cannot be called 2 times !" << endl;
    exit(0);
  }
  
  if(alpha == -1) throw std::invalid_argument("alpha must not be -1 in Lens::createOneFieldPlane");

  // non-projected number of halos :
  PosType aveNhalos = NumberDensity * pi*Rregion*Rregion; // in Number/radians^2 * radians^2 = dimensionless
  if(verbose) std::cout << "Lens::insertSubstructures : Average number of Substructures : " << aveNhalos << std::endl;
  // So numberDensity refers to the number density in 3D.
  
  std::cout << "aveNhalos = " << aveNhalos << std::endl;
  
  // projected number of halos (less than the non-projected, ie. 3D, number) :
  substructure.NhalosSub = static_cast<std::size_t>(poidev(float(aveNhalos), seed));
  if(verbose) std::cout << "Lens::insertSubstructures : Actual number of Substructures : " << substructure.NhalosSub << std::endl;
  
  // in case there is none :
  if(substructure.NhalosSub == 0)
  {
    // if(WasInsertSubStructuresCalled==YES || WasInsertSubStructuresCalled==MAYBE) deleteSubstructures();
    WasInsertSubStructuresCalled = MAYBE ;
    return;
  }
  
  // Either we called deleteSub before insertSub, or we do it here :
  if(substructure.halos.size() > 0)
  {
    // Problem: Expanding the vector is a problem if we want to add substructure
    // multiple times because if the vector is copied it will invalidate the pointers
    // on previous planes. To do this we would need to be able to expand the vector
    // without copying it or have multiple substructure_halo vectors.
    
    // throw std::runtime_error("Lens::insertSubstructures : Can only add substructure halos once to a lens.");
    std::cout << "Lens::insertSubstructures : Can only add substructure halos once to a lens. Calling deleteSubstructures." << std::endl ;
    deleteSubstructures();
    insertSubstructures(Rregion,center,NumberDensity,Mass_min,Mass_max,redshift,alpha,density_contrast,verbose);
  }
  
  PosType Dl = cosmo.angDist(redshift),rr,theta; // angular distance of the lens in PhysMpc
  PosType *theta_pos;
  PosType r = Mass_min/Mass_max,f,mass;
  size_t haloid = field_halos.size();
  PosType Rsize;
  PosType mass_min = 0.;
  if(substructure.NhalosSub > 0) mass_min = 1.e50 ;
  PosType mass_max = 0.;
  PosType SumMassSub = 0.;
  PosType rmax_max = 0.;
  
  PosType rho = density_contrast*cosmo.rho_crit(0)*cosmo.getOmega_matter()*(1+redshift)*(1+redshift)*(1+redshift);
  // rho in 1 * (M_sun/Mpc^3) * 1 * (1+z)^3 = M_sun / PhysMpc^3,
  // where Mpc \equiv comoving Mpc.
  
  for(size_t ii=0;ii<substructure.NhalosSub;++ii)
  {
    // random position
    rr = Rregion*sqrt(ran2(seed)); // in radians
    theta_pos = new PosType[3];
    theta = 2*pi*ran2(seed);       // in radians
    
    // position in proper distance
    theta_pos[0] = (rr*cos(theta) + center[0])*Dl; // in radians * angular Distance in PhysMpc = PhysMpc
    theta_pos[1] = (rr*sin(theta) + center[1])*Dl; // same : PhysMpc
    theta_pos[2] = 0.0;

    f = ran2(seed); // dimensionless and between 0 and 1
    
    // mass from power law mass function (inversing the integration of Eq. (9) in Metcalf, Amara 2011 with f \equiv (eta(m)/eta_*)*(sigma/sigma_*) ) :
    mass = Mass_max*pow( f + pow(r,alpha+1)*(1-f) , 1.0/(1+alpha) ); // in Msun, r = Mass_Min / Mass_Max is dimensionless.
    // Variable for the sum of the substructure masses :
    SumMassSub += mass ;
    assert(mass < 1.e50);
    
    // keeping track of the highest/lowest substructure mass :
    mass_max = MAX(mass,mass_max); // in Msun
    mass_min = MIN(mass,mass_min);
    
    // Rsize from tidal truncation
    Rsize = pow(mass/rho/4/pi,1.0/3.); // in [Msun / (Msun / PhysMpc^3)]^(1/3) = PhysMpc
  
    // keeping track of the highest substructure rmax :
    rmax_max = MAX(Rsize,rmax_max); // in PhysMpc
    
    // Adding the randomly-generated halo into the substructure :
    substructure.halos.push_back(new LensHaloPowerLaw(mass,Rsize,redshift,1.0,1.0,0,0));
    substructure.halos.back()->setX(theta_pos);
    ++haloid;
    substructure.halos.back()->setID(haloid);
  }
  substructure.mass_min = mass_min;
  substructure.mass_max = mass_max;
  substructure.SumMassSub = SumMassSub;
  substructure.rmax_max = rmax_max;
  
  // Averaged mass estimated for 1 subhalo from theory :
  PosType AveMassTh = Mass_max * ((1+alpha)/(2+alpha)) * ((1-pow(r,2+alpha))/(1-pow(r,1+alpha))); // Average mass for one sub halo.
  // Estimating the theoretical substructure density :
  PosType SigmaSubTh = AveMassTh * NumberDensity ; // in (Number/radians^2) * Msun = Msun / radians^2.
  // Averaged mass estimated for the whole subhalo population from theory :  //
  AveMassTh *= substructure.NhalosSub ; // Now for the total amount of sub halos.
  substructure.AveMassTh = AveMassTh;
  // Estimating the real substructure density :
  PosType SigmaSubReal = SumMassSub / (pi * Rregion * Rregion); // in (Number/radians^2) * Msun = Msun / radians^2.
  
  // Dimensionless kappa_sub or Sigma_Sub / Sigma_Crit :
  PosType DimLessKappaSub = SumMassSub / (pi * Rregion * Rregion * cosmo.angDist(redshift) * cosmo.angDist(redshift)) / cosmo.SigmaCrit(redshift,zsource);
  
  
  if(verbose)
  {
    std::cout << std::endl ;
    std::cout << "Lens::insertSubstructures : aveNhalos = " << aveNhalos << " , NhalosSub = " << substructure.NhalosSub << " , rho = " << rho << " Msun/PhysMpc^3, Dl = " << Dl << " PhysMpc." << std::endl ;
    std::cout << "Lens::insertSubstructures : Max mass = " << mass_max << " Msun , Max radius = "
    << rmax_max << " PhysMpc, Number of substructure halos = " << substructure.halos.size() << std::endl ;
    std::cout << "Lens::insertSubstructures : SumMassSub = " << SumMassSub << " Msun, Theoretical total averaged mass = " << AveMassTh << " Msun." << std::endl ;
    std::cout << "Lens::insertSubstructures : Sigma_crit = " << cosmo.SigmaCrit(redshift, zsource) << " Msun/PhysMpc^2." << std::endl ;
    std::cout << "Lens::insertSubstructures : 4 pi G = " << 4*pi*Grav << " Mpc/Msun." << std::endl ;
    std::cout << "Lens::insertSubstructures : SigmaSubTh = " << SigmaSubTh << " Msun/radians^2 , SigmaSubReal = " << SigmaSubReal << " Msun/radians^2." << std::endl ;
    std::cout << "Lens::insertSubstructures : SigmaSubTh = " << SigmaSubTh *cosmo.angDist(redshift)*cosmo.angDist(redshift)*(1+redshift)*(1+redshift)*1.e-6 << " Msun/Physkpc^2 , SigmaSubReal = " << SigmaSubReal*cosmo.angDist(redshift)*cosmo.angDist(redshift)*(1+redshift)*(1+redshift)*1.e-6 << " Msun/Physkpc^2." << std::endl ;
    std::cout << "Lens::insertSubstructures : DimLessKappaSub = " << DimLessKappaSub << endl;
  }

  // Test :
  // std::cout << " field_plane_redshifts.begin() : " << field_plane_redshifts.size() << std::endl ;
  // for(int i = 0 ; i < field_plane_redshifts.size() ; i++) std::cout << field_plane_redshifts[i] << " " ;
  // std::cout << std::endl ;
  
  // the new plane must be inserted in order of redshift
  if(field_Nplanes_current != 0)
  {
    std::cout << field_planes.size() << " " << field_plane_redshifts.size() << " " << field_Dl.size() << " " << field_Nplanes_original << std::endl ;
    assert(field_planes.size() == field_Nplanes_original);
    assert(field_plane_redshifts.size() == field_Nplanes_original);
    assert(field_Dl.size() == field_Nplanes_original);

    cout << "field : " << field_planes.size() << " " << field_plane_redshifts.size() << " " << field_Dl.size() << endl ;
    cout << "field redshifts :" << endl;
    for(int i=0;i<field_plane_redshifts.size();i++) cout << field_plane_redshifts[i] << " " ;
    cout << endl;
    
    std::vector<LensPlane*>::iterator it = field_planes.begin();
    std::vector<PosType>::iterator itz = field_plane_redshifts.begin();
    std::vector<PosType>::iterator itd = field_Dl.begin();
    while(*itz < redshift && it != field_planes.end()){
      ++it;
      ++itz;
      ++itd;
    }
    if(verbose) std::cout << "Lens::insertSubstructures : redshift " << redshift << " nearest plane z = " << *itz << std::endl;
    it = field_planes.insert(it, new LensPlaneTree(substructure.halos.data(), substructure.NhalosSub, 0., 0));
    field_plane_redshifts.insert(itz,redshift);
    field_Dl.insert(itd,Dl*(1+redshift));
    substructure.plane = *it;
  }
  else // in the case where no field plane exists
  {
    // Before insertion :
    if(verbose)
    {
      std::cout << "Lens::insertSubstructures : Before insertion of plane :" << std::endl;
      if(field_plane_redshifts.size()==0) std::cout << "X" ;
      for (int i=0;i<field_plane_redshifts.size();i++) std::cout << field_plane_redshifts[i] << " " ;
      std::cout << std::endl ;
      if(field_Dl.size()==0) std::cout << "X" ;
      for (int i=0;i<field_Dl.size();i++) std::cout << field_Dl[i] << " " ;
      std::cout << std::endl ;
    }
    
    // Insertion :
    if(verbose) std::cout << "Lens::insertSubstructures : inserting a new plane at redshift z = " << redshift << std::endl;
    field_planes.push_back(new LensPlaneTree(substructure.halos.data(), substructure.NhalosSub, 0, 0));
    field_plane_redshifts.push_back(redshift);
    field_Dl.push_back(Dl*(1+redshift));
    substructure.plane = field_planes[0];

    // After insertion :
    if(verbose)
    {
      std::cout << "Lens::insertSubstructures : After insertion of plane :" << std::endl;
      for (int i=0;i<field_plane_redshifts.size();i++) std::cout << field_plane_redshifts[i] << " " ;
      std::cout << std::endl ;
      for (int i=0;i<field_Dl.size();i++) std::cout << field_Dl[i] << " " ;
      std::cout << std::endl ;
    }
  }
  ++field_Nplanes_current;
  
  combinePlanes(verbose);
  std::cout << "Lens::insertSubstructures : field_planes.size() = " << field_planes.size() << std::endl;
  // assert(field_planes.size() == field_Nplanes);
  
  WasInsertSubStructuresCalled = YES ;
  std::cout << field_planes.size() << " " << field_Nplanes_current << std::endl ;
  assert(field_planes.size() == field_Nplanes_current);

  // Test :
  // std::cout << "field_plane and field_plane_redshifts : " << std::endl ;
  // for(int k=0 ; k<field_planes.size() ; k++) std::cout << field_planes[k] << " " << field_plane_redshifts[k] << std::endl ;
  // std::cout << std::endl ;

}



/// * RESET SUB STRUCTURE * ///

void Lens::resetSubstructure(bool verbose){

  // new verbose choice :
  substructure.verbose = verbose ;
  
  // Reinitiating varying quantities stored in substructure :
  substructure.NhalosSub = 0.;
  substructure.mass_min = 0.;
  substructure.mass_max = 0.;
  substructure.rmax_max = 0.;
  substructure.SumMassSub = 0.;
  substructure.AveMassTh = 0.;
  
  // test of whether insertSubstructures has been called :
  if(WasInsertSubStructuresCalled == NO)
  {
    ERROR_MESSAGE();
    cout << "Lens::insertSubStructures() has to be called before Lens::resetSubStructure() !" << endl;
    exit(0);
  }
  else if(WasInsertSubStructuresCalled == MAYBE)
  {
    // Clearing the substructure halos :
    Utilities::delete_container(substructure.halos);
    // substructure.halos.clear();
    
    // Some checks of consistency :
    assert(field_Nplanes_current == field_Nplanes_original);
    assert(field_plane_redshifts.size() == field_plane_redshifts_original.size());
    assert(field_Dl.size() == field_Dl_original.size());
    
    // Changing the flag :
    WasInsertSubStructuresCalled = NO ;
    
    // Reconstructing the plane with arguments given to insertSubStructures :
    insertSubstructures(substructure.Rregion,substructure.center.x,substructure.Ndensity,substructure.Mmin,substructure.Mmax,substructure.redshift,substructure.alpha,substructure.rho_tidal,verbose);
    return ;
  }
  // We have WasInsertSubStructuresCalled = YES after this point.
  
  // find which plane has the substructures on it
  int fplane_index = 0,lplane_index = 0;
  while(field_planes[fplane_index] != substructure.plane) ++fplane_index;
  while(lensing_planes[lplane_index] != substructure.plane) ++lplane_index;
  
  // redshift and Dl of substructure plane :
  PosType redshift = field_plane_redshifts[fplane_index];
  assert(redshift == substructure.redshift);
  PosType Dlsub = Dl[fplane_index]/(1+redshift);
  
  // Computing the number of halos in substructure :
  PosType aveNhalos = substructure.Ndensity * pi * substructure.Rregion * substructure.Rregion; // in Number/radians^2 * radians^2 = dimensionless
  substructure.NhalosSub = static_cast<std::size_t>(poidev(float(aveNhalos), seed));
  if(verbose)
  {
    std::cout << "Lens::resetSubstructures : Average number of Substructures : " << aveNhalos << std::endl;
    std::cout << "Lens::resetSubstructures : Actual number of Substructures : " << substructure.NhalosSub << std::endl;
  }
  
  //  Testing if NhalosSub = 0 :
  if (substructure.NhalosSub == 0) std::cout << "Be careful ! NhalosSub = 0 in resetSubstructure !" << std::endl;
  
  PosType rr,theta;
  PosType *theta_pos;
  PosType r = substructure.Mmin/substructure.Mmax,f,mass;
  PosType Rsize;
  size_t haloid = field_halos.size();
  PosType mass_min = 0.;
  if(substructure.NhalosSub > 0) mass_min = 1.e50 ;
  PosType mass_max = 0.;
  PosType SumMassSub = 0.;
  PosType rmax_max = 0.;
  
  PosType rho = substructure.Ndensity *cosmo.rho_crit(0)*cosmo.getOmega_matter()*(1+redshift)*(1+redshift)*(1+redshift);
  if(verbose) std::cout << "Lens::resetSubstructures : rho = " << rho << " , Dlsub = " << Dlsub << std::endl ;
  
  // Clearing the previous halos :
  Utilities::delete_container(substructure.halos);
  assert(substructure.halos.size() == 0);
  
  // Re-creating the new distribution of substructure halos :
  for(size_t ii=0;ii<substructure.NhalosSub;++ii)
  {
    // random position
    rr = substructure.Rregion*sqrt(ran2(seed));
    theta_pos = new PosType[3];
    theta = 2*pi*ran2(seed);
    
    // position in proper distance
    theta_pos[0] = (rr*cos(theta) + substructure.center[0])*Dlsub;
    theta_pos[1] = (rr*sin(theta) + substructure.center[1])*Dlsub;
    theta_pos[2] = 0.0;
    
    f = ran2(seed);
    
    // mass from power law mass function
    mass = substructure.Mmax * pow( f + pow(r,substructure.alpha+1)*(1-f), 1.0/(1+substructure.alpha) );
    // Variable for the sum of the substructure masses :
    SumMassSub += mass ;
    assert(mass < 1.e50);
    
    // keeping track of the highest/lowest substructure mass :
    mass_max = MAX(mass,mass_max); // in Msun
    mass_min = MIN(mass,mass_min);
    
    // Rsize from tidal truncation
    Rsize = pow(mass/rho/4/pi,1.0/3.);
    
    rmax_max = MAX(Rsize,rmax_max);
    
    substructure.halos.push_back(new LensHaloPowerLaw(mass,Rsize,redshift,1.0,1.0,0,0));
    substructure.halos.back()->setX(theta_pos);
    ++haloid;
    substructure.halos.back()->setID(haloid);
  }
  substructure.mass_min = mass_min;
  substructure.mass_max = mass_max;
  substructure.SumMassSub = SumMassSub;
  substructure.rmax_max = rmax_max;
  
  // Averaged mass estimated from theory :
  PosType AveMassTh = substructure.Mmax * ((1+substructure.alpha)/(2+substructure.alpha)) * ((1-pow(r,2+substructure.alpha))/(1-pow(r,1+substructure.alpha))); // Average mass for one sub halo.
  // Estimating the theoretical substructure density :
  // PosType SigmaSubTh = AveMassTh * substructure.Ndensity ; // in (Number/radians^2) * Msun = Msun / radians^2.
  // Averaged mass estimated for the whole subhalo population from theory :  //
  AveMassTh *= substructure.NhalosSub ; // Now for the total amount of sub halos.
  substructure.AveMassTh = AveMassTh;
  // Estimating the real substructure density :
  // PosType SigmaSubReal = SumMassSub / (pi * substructure.Rregion * substructure.Rregion); // in (Number/radians^2) * Msun = Msun / radians^2.
  
  if(verbose){
    std::cout << "Lens::resetSubstructures : Max mass = " << mass_max << " , Max radius = "
    << rmax_max << " , number of substructure halos = " << substructure.halos.size() << std::endl ;
  }
  
  assert(substructure.halos.size() == substructure.NhalosSub);

  // replacing the substructure plane with the new plane :
  delete substructure.plane;
  lensing_planes[lplane_index] = field_planes[fplane_index] = substructure.plane = new LensPlaneTree(substructure.halos.data(), substructure.NhalosSub, 0, 0);
}







/// * DELETE SUB STRUCTURE * ///

void Lens::deleteSubstructures()
{
  if(WasInsertSubStructuresCalled == NO)
  {
    ERROR_MESSAGE();
    cout << "Lens::insertSubStructures() has to be called before Lens::deleteSubstructure() !" << endl;
    exit(0);
  };
  
  /*
   std::cout << " ===== deleteSubstructure ===== " << std::endl ;
   for(int i=0;i<field_planes.size();i++) std::cout << field_planes[i] << " " ;
   std::cout << std::endl;
   for(int i=0;i<field_plane_redshifts.size();i++) std::cout << field_plane_redshifts[i] << " " ;
   std::cout << std::endl;
   for(int i=0;i<field_Dl.size();i++) std::cout << field_Dl[i] << " " ;
   std::cout << std::endl;
   std::cout << " ========================= " << std::endl ;
   */
  
  if(field_planes.size() > 0)
  {
    // Finding the iterators that correspond to the substructure plane :
    std::vector<LensPlane*>::iterator it = field_planes.begin();
    std::vector<PosType>::iterator itz = field_plane_redshifts.begin();
    std::vector<PosType>::iterator itd = field_Dl.begin();
    while(*itz < substructure.redshift && it != field_planes.end()){
      ++it;
      ++itz;
      ++itd;
    }
    
    // Erasing the corresponding quantities :
    field_planes.erase(it);
    field_plane_redshifts.erase(itz);
    field_Dl.erase(itd);
    
    field_Nplanes_current -= 1;
  }
  
  /*
   std::cout << " ===== deleteSubstructure ===== " << std::endl ;
   for(int i=0;i<field_planes.size();i++) std::cout << field_planes[i] << " " ;
   std::cout << std::endl;
   for(int i=0;i<field_plane_redshifts.size();i++) std::cout << field_plane_redshifts[i] << " " ;
   std::cout << std::endl;
   for(int i=0;i<field_Dl.size();i++) std::cout << field_Dl[i] << " " ;
   std::cout << std::endl;
   std::cout << " ========================= " << std::endl ;
   */
  
  // Recombining the planes :
  combinePlanes(false);
  
  // Clearing the substructure from the previous halos :
  substructure.halos.clear();
  // substructure.plane = new LensPlaneTree(substructure.halos.data(), 0, 0, 0);
  delete(substructure.plane);
  
  // Clearing some stored quantities :
  substructure.Rregion = 0;
  substructure.Mmax = 0;
  substructure.Mmin = 0;
  substructure.alpha = 0;
  substructure.Ndensity = 0;
  substructure.center[0] = substructure.center[1] = 0.;
  substructure.rho_tidal = 0;
  substructure.redshift = 0;
  substructure.verbose = false;
  substructure.NhalosSub = 0;
  substructure.mass_min = 0;
  substructure.mass_max = 0;
  substructure.rmax_max = 0;
  substructure.SumMassSub = 0;
  substructure.AveMassTh = 0;
  
  // Resetting the flag :
  WasInsertSubStructuresCalled = NO ;
  
  return ;
}




void Lens::addMainHaloToPlane(LensHalo* halo)
{
	// the redshift and distance of the halo
	PosType halo_z = halo->getZlens();
	PosType halo_Dl = cosmo.coorDist(0, halo_z);
	
	// find the position of the new lens plane
	std::size_t i = std::distance(main_Dl.begin(), std::upper_bound(main_Dl.begin(), main_Dl.end(), halo_Dl));
	
	// go through all options for adding
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
  
	Utilities::MixedVector<LensHalo*>::iterator<> it = main_halos.begin();
	while(it != main_halos.end())
	{
		// find halos with higher redshift
		Utilities::MixedVector<LensHalo*>::iterator<> jt = std::upper_bound(it, main_halos.end(), *it, lens_halo_less(cosmo));
		
		// add halos until higher redshift to plane
		main_planes.push_back(new LensPlaneSingular(&(*it), std::distance(it, jt)));
		
		// add plane to arrays
		main_plane_redshifts.push_back((*it)->getZlens());
		main_Dl.push_back(cosmo.coorDist(0, main_plane_redshifts.back()));
		
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
	PosType Dmax = cosmo.coorDist(0, zsource);
	
	std::vector<PosType> lD;
	std::size_t Np = field_Nplanes_current + 1;
	
	assert(Np > 1);
	
	// spaces interval equally up to the source, including 0 and Dmax
	// therefore we need field_Nplanes+1 values
	Utilities::fill_linear(lD, Np, 0.0, Dmax);
	
	// spacing of the distances
	PosType dlD = lD[1]-lD[0];
	
	// assigns the redshifts and plugs in the input plane
	for(std::size_t i = 1; i < Np; ++i)
	{
		// ensures that the first plane and the last before the source plane
		// have the same volume as all the other planes
		lD[i] -= 0.5*dlD;
		
		field_Dl.push_back(lD[i]);
	}
	
	assert(field_Dl.size() == field_Nplanes_current);
	
	// assigns the redshifts and plugs in the input plane
	for(std::size_t i = 0; i < field_Nplanes_current; ++i)
	{
		// get redshift for calculated distance
		PosType z = cosmo.invCoorDist(field_Dl[i]);
		field_plane_redshifts.push_back(z);
		
		// refit the distances to match the redshift
		field_Dl[i] = cosmo.coorDist(0, z);
	}
	
	assert(field_plane_redshifts.size() == field_Nplanes_current);
  field_plane_redshifts_original = field_plane_redshifts;
  field_Dl_original = field_Dl;
}

void Lens::setFieldDistFromFile()
{
	PosType value;
	
	std::ifstream file_in(redshift_planes_file.c_str());
	if(!file_in)
		throw std::runtime_error("Can't open file " + redshift_planes_file);
	
	while(file_in >> value)
	{
		if(!value)
			throw std::runtime_error("can't read PosType from " + redshift_planes_file);
		else
			field_plane_redshifts.push_back(value);
	}
	
	file_in.close();
	
	assert(field_plane_redshifts.size() == field_Nplanes_current);
	
	for(std::size_t i = 0; i < field_plane_redshifts.size(); ++i)
		field_Dl.push_back(cosmo.coorDist(0, field_plane_redshifts[i]));
  
  field_plane_redshifts_original = field_plane_redshifts;
  field_Dl_original = field_Dl;
}

/**
 * \brief Creates main lens halo as set up in the parmeter file.
 *
 */
void Lens::createMainHalos(InputParams& params)
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
      main_halos.push_back(new LensHaloRealNSIE(params));
      break;
    case ana_lens:
      main_halos.push_back(new LensHaloAnaNSIE(params));
      break;
    case uni_lens:
      main_halos.push_back(new LensHaloUniform(params));
      break;
    case moka_lens:
      main_halos.push_back(new LensHaloMassMap(params, cosmo));
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
  
  if(pixel_map_on) readPixelizedDensity();
  
	if(main_galaxy_halo_type != 0)
	{
		switch(main_galaxy_halo_type)
		{
      case null_gal:
        break;
      case nsie_gal:
        main_halos.push_back(new LensHaloRealNSIE(params));
        break;
      case pl_gal:
        main_halos.push_back(new LensHaloPowerLaw(params));
        break;
      case hern_gal:
        main_halos.push_back(new LensHaloHernquist(params));
        break;
      case jaffe_gal:
        main_halos.push_back(new LensHaloJaffe(params));
        break;
		}
	}
  
	for(std::size_t i = 0; i < main_halos.size(); ++i){
		main_halos[i]->setCosmology(cosmo);
  }
}

void Lens::clearMainHalos(bool verbose)
{
	Utilities::delete_container(main_halos_created);
	main_halos.clear();
	
	flag_switch_main_halo_on = false;
  
	
	Utilities::delete_container(main_planes);
	main_plane_redshifts.clear();
	main_Dl.clear();
	
	combinePlanes(verbose);
}

/**
 * \brief Inserts a single main lens halo.
 * Then all lensing planes are updated accordingly.
 */
void Lens::insertMainHalo(LensHalo* halo,bool verbose)
{
	halo->setCosmology(cosmo);
	main_halos.push_back(halo);
  
	flag_switch_main_halo_on = true;
	
	addMainHaloToPlane(halo);
	
	combinePlanes(verbose);
}

/**
 * \brief Inserts a sequense of main lens halos and ads them to the existing ones.
 * Then all lensing planes are updated accordingly.
 */
void Lens::insertMainHalos(LensHalo** halos, std::size_t Nhalos, bool verbose)
{
	for(std::size_t i = 0; i < Nhalos; ++i)
	{
		halos[i]->setCosmology(cosmo);
		main_halos.push_back(halos[i]);
		addMainHaloToPlane(halos[i]);
	}
	
	flag_switch_main_halo_on = true;
	
	combinePlanes(verbose);
}




/**
 * \brief Inserts a single main lens halo and deletes all previous ones.
 * Then all lensing planes are updated accordingly.
 */
void Lens::replaceMainHalos(LensHalo* halo,bool verbose)
{
	Utilities::delete_container(main_halos_created);
	main_halos.clear();
	
	halo->setCosmology(cosmo);
	main_halos.push_back(halo);
	
	flag_switch_main_halo_on = true;
	
	Utilities::delete_container(main_planes);
	createMainPlanes();
	combinePlanes(verbose);
}

/**
 * \brief Inserts a sequense of main lens halos and deletes all previous ones.
 * Then all lensing planes are updated accordingly.
 */
void Lens::replaceMainHalos(LensHalo** halos, std::size_t Nhalos,bool verbose)
{
	Utilities::delete_container(main_halos_created);
	main_halos.clear();
	
	for(std::size_t i = 0; i < Nhalos; ++i)
	{
		halos[i]->setCosmology(cosmo);
		main_halos.push_back(halos[i]);
	}
	
	flag_switch_main_halo_on = true;
	
	Utilities::delete_container(main_planes);
	createMainPlanes();
	combinePlanes(verbose);
}


/**
 * \brief Computes some quantities necessary for the function createFieldHalos. This material was before computed in createFieldHalos.
 *
 * Especially it computes the quantities related to the mass function. By calling this function in the constructor of the lens, we make that these quantities are stored in the lens and do not have to be recomputed for each new realisation of the field.
 */

void Lens::ComputeHalosDistributionVariables ()
{
  const PosType MaxLogm=16.;
  
  aveNhalosField = cosmo.haloNumberInBufferedCone(field_min_mass,0,zsource,fieldofview*pow(pi/180,2),field_buffer,field_mass_func_type,mass_func_PL_slope);
  
  Utilities::fill_linear(zbins,Nzbins,0.0,zsource);
  // construct redshift distribution table
  NhalosbinZ[0] = 1;
  zbins[0] = 0;
  
  for(int k=1;k<Nzbins-1;++k){
    NhalosbinZ[k] = cosmo.haloNumberInBufferedCone(field_min_mass,zbins[k],zsource,fieldofview*pow(pi/180,2),field_buffer,field_mass_func_type,mass_func_PL_slope)/aveNhalosField;
  }
  // std::cout << std::endl ;
  zbins[Nzbins-1] = zsource;
  NhalosbinZ[Nzbins-1] = 0.0;
  
  // fill the log(mass) vector
  Logm.resize(Nmassbin);
  Utilities::fill_linear(Logm,Nmassbin,log10(field_min_mass),MaxLogm);
  
  // this will be used for the cumulative number density in one square degree
  NhalosbinMass.resize(NZSamples);
  for(int np=0;np<NZSamples;np++) NhalosbinMass[np].resize(Nmassbin);
  
  for(int np=0;np<NZSamples;np++){
    PosType z1, z2;
    z1 = np*zsource/(NZSamples);
    z2 = (np+1)*zsource/(NZSamples);
    
    Nhaloestot_Tab[np] = cosmo.haloNumberInBufferedCone(pow(10,Logm[0]),z1,z2,fieldofview*pow(pi/180,2),field_buffer,field_mass_func_type,mass_func_PL_slope);
    // std::cout << Nhaloestot_Tab[np] << " "  ;
    
    NhalosbinMass[np][0] = 1;
    for(int k=1;k<Nmassbin-1;k++){
      // cumulative number density in one square degree
      NhalosbinMass[np][k] = cosmo.haloNumberInBufferedCone(pow(10,Logm[k]),z1,z2,fieldofview*pow(pi/180,2),field_buffer,field_mass_func_type,mass_func_PL_slope)/Nhaloestot_Tab[np];
    }
    NhalosbinMass[np][Nmassbin-1] = 0;
    
  }
  // std::cout << std::endl ;
  
}

/**
 * \brief Creates the field of halos as specified in the parameter file.
 *
 */
void Lens::createFieldHalos(bool verbose)
{
  std::cout << "Creating Field Halos from Mass Function" << std::endl;
  //verbose = true;
  
	unsigned long i,j_max,k1,k2;
  PosType z_max;
	PosType z1, z2, mass_max,Nhaloestot;
	int np;
	PosType rr,theta,maxr;
	HALOCalculator *halo_calc = new HALOCalculator(&cosmo,field_min_mass,0.0);
  PosType mo=7.3113e10,M1=2.8575e10,gam1=7.17,gam2=0.201,be=0.557;
  PosType field_galaxy_mass_fraction = 0;
  PosType r200=0, r_half_stel_mass=0;
  size_t haloid=0;
  
  if (field_min_mass < 1.0e5) {
    std::cout << "Are you sure you want the minimum field halo mass to be " << field_min_mass << " Msun?" << std::endl;
    throw;
  }
  
  std::vector<PosType> halo_zs_vec;
  //std::vector<PosType *> halo_pos_vec;  /// temporary vector to store angular positions
  
  // assign redshifts to field_halos according to the redshift distribution
  
  int Nhalos = static_cast<std::size_t>(poidev(float(aveNhalosField), seed));

  for(int i=0;i < Nhalos;++i){
    halo_zs_vec.push_back(Utilities::InterpolateYvec(NhalosbinZ,zbins,ran2(seed)));
  }
  
  if(verbose){
    std::cout << "redshift distribution function" << std::endl;
    std::cout << "zbins[]       NhalosbinZ[]" << std::endl;
    
    for(int i=0;i<zbins.size();++i){
      std::cout << zbins[i] << "  " << NhalosbinZ[i] << std::endl;
    }
  }
  
  // sort redshifts
  std::sort(halo_zs_vec.begin(),halo_zs_vec.end());
  
  assert(halo_zs_vec[0] < halo_zs_vec[1]);
  assert(halo_zs_vec[0] < halo_zs_vec[Nhalos-1]);
  
	PosType *theta_pos,*theta2;
	size_t j = 0;
	k2 = 0;
	std::vector<PosType>::iterator it1,it2;

  ///////////////////////////////////////
	for(np=0,mass_max=0;np<NZSamples;np++){

		z1 = np*zsource/(NZSamples);
		z2 = (np+1)*zsource/(NZSamples);

		it1 = std::lower_bound(halo_zs_vec.begin(),halo_zs_vec.end(),z1);
		it2 = std::lower_bound(halo_zs_vec.begin(),halo_zs_vec.end(),z2);
    
		k1 = it1 - halo_zs_vec.begin();
		k2 = it2 - halo_zs_vec.begin();

    Nhaloestot = Nhaloestot_Tab[np] ;
    
    if(verbose){
      std::cout << "   np = " << np << " z: " << z1 << " " << z2 << std::endl;
      std::cout << "          n = " << k2 - k1 << std::endl;
    }
    
		for(i = k1; i < k2; i++){
			PosType Ds = cosmo.angDist(0,halo_zs_vec[i]);
      
			maxr = pi*sqrt(fieldofview/pi)/180. + field_buffer/Ds; // fov is a circle
			rr = maxr*sqrt(ran2(seed));
      
      if(verbose) std::cout << "          maxr = " << maxr << std::endl;

			assert(rr == rr);
      
			theta_pos = new PosType[3];
      
			theta = 2*pi*ran2(seed);
      
			theta_pos[0] = rr*cos(theta);//*Ds;
			theta_pos[1] = rr*sin(theta);//*Ds;
			theta_pos[2] = 0.0;
      
      float mass = pow(10,Utilities::InterpolateYvec(NhalosbinMass[np],Logm,ran2 (seed)));
      
			halo_calc->reset(mass,halo_zs_vec[i]);
      
			float Rsize = halo_calc->getRvir();
			float rscale = Rsize/halo_calc->getConcentration(0);
      assert(rscale < Rsize);
      
      float sigma = 0;
      if(flag_field_gal_on){
        // from Moster et al. 2010ApJ...710..903M
        field_galaxy_mass_fraction = mo*pow(mass/M1,gam1)
        /pow(1+pow(mass/M1,be),(gam1-gam2)/be)/mass;
        if(field_galaxy_mass_fraction > 1.0) field_galaxy_mass_fraction = 1;
        sigma = 126*pow(mass*(1-field_galaxy_mass_fraction)/1.0e10,0.25); // From Tully-Fisher and Bell & de Jong 2001
				//field_halos[j]->initFromMassFunc(mass*(1-field_galaxy_mass_fraction),Rsize,rscale,field_prof_internal_slope,seed);
        r200 = halo_calc->getR200(); // used for Kravtsov 2013 2013ApJ...764L..31K
        r_half_stel_mass = pow(10.,(0.95*log10(0.015*r200*1000.)+0.015))/1000.; // in Mpc
      
      }else{
        field_galaxy_mass_fraction = 0;
      }
      
			switch(field_int_prof_type)
			{
				case null_lens:
					ERROR_MESSAGE();
					std::cout << "field_int_prof_type is null!!!!" << std::endl;
					break;
				case nfw_lens:
					//field_halos.push_back(new LensHaloNFW);
          field_halos.push_back(new LensHaloNFW(mass*(1-field_galaxy_mass_fraction),Rsize,halo_zs_vec[i],Rsize/rscale,1.0,0,0));
					break;
				case pnfw_lens:
					//field_halos.push_back(new LensHaloPseudoNFW);
          field_halos.push_back(new LensHaloPseudoNFW(mass*(1-field_galaxy_mass_fraction),Rsize,halo_zs_vec[i],Rsize/rscale,3,1.0,0,0));
					break;
				case pl_lens:
					//field_halos.push_back(new LensHaloPowerLaw);
          field_halos.push_back(new LensHaloPowerLaw(mass*(1-field_galaxy_mass_fraction),Rsize,halo_zs_vec[i],1.0,1.0,0,0));
          
					break;
				case nsie_lens:
          //std::cout << "Warning: All galaxies are spherical" << std::endl;
					field_halos.push_back(new LensHaloRealNSIE(mass*(1-field_galaxy_mass_fraction),halo_zs_vec[i],sigma,0.0,1.0,0,0));
					//field_halos.push_back(new LensHaloRealNSIE);
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
          field_halos[j]->setZlens(halo_zs_vec[i]);
					break;
				case hern_lens:
					//field_halos.push_back(new LensHaloHernquist);
          field_halos.push_back(new LensHaloHernquist(mass*(1-field_galaxy_mass_fraction),Rsize,halo_zs_vec[i],rscale,1.0,0,0));
					break;
				case jaffe_lens:
					//field_halos.push_back(new LensHaloJaffe);
          field_halos.push_back(new LensHaloJaffe(mass*(1-field_galaxy_mass_fraction),Rsize,halo_zs_vec[i],rscale,1.0,0,0));
					break;
			}
      
			if(mass > mass_max) {
				mass_max = mass;
				j_max = i;
				z_max = halo_zs_vec[i];
			}
      
			//halo_pos_vec.push_back(theta_pos);
      field_halos.back()->setX(theta_pos);
      ++haloid;
      field_halos.back()->setID(haloid);
      
			++j;
      
			if(flag_field_gal_on){
        switch(field_int_prof_gal_type){
          case pl_gal:
            ERROR_MESSAGE();
            std::cout << "field_int_prof_gal_type 2, i.e. PowerLaw not yet implemented!!!" << std::endl;
            break;
          case hern_gal:
            ERROR_MESSAGE();
            std::cout << "field_int_prof_gal_type 3, i.e. Hernquist not yet implemented!!!" << std::endl;
            break;
          case jaffe_gal:
            ERROR_MESSAGE();
            std::cout << "field_int_prof_gal_type 4, i.e. Jaffe not yet implemented!!!" << std::endl;
            break;
          case null_gal:
            ERROR_MESSAGE();
            std::cout << "flag_field_gal_on is true, but field_int_prof_gal_type is null!!!" << std::endl;
            break;
            
          case nsie_gal:
            
            float sigma = 126*pow(mass*field_galaxy_mass_fraction/1.0e10,0.25); // From Tully-Fisher and Bell & de Jong 2001
            //std::cout << "Warning: All galaxies are spherical" << std::endl;
            float fratio = (ran2(seed)+1)*0.5;  //TODO: Ben change this!  This is a kluge.
            float pa = 2*pi*ran2(seed);  //TODO: This is a kluge.
            field_halos.push_back(new LensHaloRealNSIE(mass*field_galaxy_mass_fraction,halo_zs_vec[i],sigma,0.0,fratio,pa,0));
            
            //field_halos[j]->initFromMassFunc(mass*field_galaxy_mass_fraction,Rsize,rscale,field_prof_internal_slope,seed);
            break;
				}
        
        // Another copy of this position must be made to avoid rescaling it twice when it is converted into
        // distance on the lens plane in Lens::buildLensPlanes()
        theta2 = new PosType[3];
        theta2[0]=theta_pos[0]; theta2[1]=theta_pos[1]; theta2[2]=theta_pos[2];
        
        //halo_pos_vec.push_back(theta2);
        field_halos.back()->setX(theta2);
        field_halos.back()->setID(haloid);
        
        ++j;
        
      }
      
		}
	}
  
	assert(k2 == Nhalos);
	delete halo_calc;
  
  if(verbose){
    std::cout << field_halos.size() << " halos created." << std::endl
    << "   largest mass: " << mass_max << "  at redshift redshift: " << z_max << std::endl;
    
  }
  
	if(verbose) std::cout << "leaving Lens::createFieldHalos()" << std::endl;
}


/**
 * \brief Read in information from a Virgo Millennium Data Base http://gavo.mpa-garching.mpg.de/MyMillennium/
 *
 * query select * from MockLensing.dbo.m21_20_39_021_bc03_Ben_halos
 *
 * This is information on the dark matter field_halos only.  There are 13 entries in each line separated by commas.
 *
 *
 * The LensHalo id numbers are set to the id of the parent halo in the simulation.  This
 * is not the FOF parent.  Multiple LensHalos can have the same id because they were derived from the same simulation halo (ex. one for galaxy and one for DM halo).
 *
 */
void Lens::readInputSimFileMillennium(bool verbose)
{
    
  std::cout << "Reading Field Halos from " << field_input_sim_file << std::endl;
	PosType z,vmax,vdisp,r_halfmass;
	unsigned long i,j;
	unsigned long haloid,idd,np;
	const PosType mo=7.3113e10,M1=2.8575e10,gam1=7.17,gam2=0.201,be=0.557;
	HALOCalculator *halo_calc = new HALOCalculator(&cosmo,field_min_mass,0.0);
  PosType field_galaxy_mass_fraction = 0;
	PosType rmax=0,rtmp=0;
  PosType theta[2];
  PosType r_half_stel_mass=0, r200=0;
  
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
	if(verbose) std::cout << "skipped "<< i << " comment lines in file " << field_input_sim_file << std::endl;
  
	//std::vector<PosType *> halo_pos_vec;
  Utilities::Geometry::SphericalPoint tmp_sph_point(1,0,0);
  
	// read in data
	int j_max;
	PosType mass_max=0,R_max=0,V_max=0,minmass=1e30;
	//PosType *theta,*theta2;
	int ncolumns = 9;
	//int ncolumns = 13;
  
	void *addr[ncolumns];
  
	addr[0] = &haloid;
	addr[1] = &idd;
	//addr[2] = &ra;
	//addr[3] = &dec;
	addr[2] = &(tmp_sph_point.phi);
	addr[3] = &(tmp_sph_point.theta);
	addr[4] = &z;
	addr[5] = &np;
	addr[6] = &vdisp;
	addr[7] = &vmax;
	addr[8] = &r_halfmass;
  
	unsigned long myint;
	PosType myPosType;
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
				buffer >> myPosType;
				*((PosType *)addr[l]) = myPosType;
			}
			myline.erase(0,pos+1);
			strg.clear();
			buffer.clear();
			buffer.str(std::string());
		}
    
    if(np == 0) continue;
		//PosType Ds = cosmo->angDist(0,z);
		//theta[0] = -ra*pi/180.;
		//theta[1] = dec*pi/180.;
    tmp_sph_point.theta *= pi/180;
    tmp_sph_point.phi *= -pi/180;
    
    rtmp = Utilities::Geometry::AngleSeporation(central_point_sphere,tmp_sph_point);
    if(sim_angular_radius > 0 && rtmp > sim_angular_radius) continue;
    if(rtmp > rmax) rmax = rtmp;
    
		// position on lens plane
		//theta = new PosType[2];
    tmp_sph_point.OrthographicProjection(central_point_sphere,theta);
    
    float mass = np*8.6e8/cosmo.gethubble(),sigma=0;
    halo_calc->reset(mass,z);
    
    
    if(flag_field_gal_on){
      // from Moster et al. 2010ApJ...710..903M
      field_galaxy_mass_fraction = mo*pow(mass/M1,gam1)
      /pow(1+pow(mass/M1,be),(gam1-gam2)/be)/mass;
      if(field_galaxy_mass_fraction > 1.0) field_galaxy_mass_fraction = 1;
      sigma = 126*pow(mass*(1-field_galaxy_mass_fraction)/1.0e10,0.25); // From Tully-Fisher and Bell & de Jong 2001
      r200 = halo_calc->getR200(); // used for Kravtsov 2013 2013ApJ...764L..31K ;
      r_half_stel_mass = pow(10.,(0.95*log10(0.015*r200*1000.)+0.015))/1000.; // in Mpc
    }else{
      field_galaxy_mass_fraction = 0;
    }
    
    if(np > 0.0 && vdisp > 0.0 && z <= zsource){
      
      switch(field_int_prof_type)
			{
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
          field_halos.push_back(new LensHaloRealNSIE(mass*field_galaxy_mass_fraction,z,sigma,0.0,1.0,0.0,0));
          
					//field_halos.push_back(new LensHaloRealNSIE);
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
      field_halos[j]->setID(haloid);
      
			field_halos[j]->setZlens(z);
      if(field_int_prof_type != nsie_lens){
        
        field_halos[j]->initFromFile(mass*(1-field_galaxy_mass_fraction),seed,vmax,r_halfmass*cosmo.gethubble());
			}
      
      
			if(field_halos[j]->get_Rmax() > R_max) R_max = field_halos[j]->get_Rmax();
			if(vdisp > V_max) V_max = vdisp;
      
			//halo_pos_vec.push_back(theta);
      field_halos.back()->setX(theta);
      
      /****** test **********
      {
        PosType tmpx[2];
        field_halos.back()->getX(tmpx);
        std::cout << "halo " << tmpx[0] << "  " << tmpx[1] << " " << field_halos.back()->get_mass()
        << " " << field_halos.back()->get_Rmax() << " " << field_halos.back()->getID() << std::endl;
      }*/
      
			if(mass > mass_max) {
				mass_max = mass;
				j_max = j;
			}
			if(mass < minmass) {
				minmass = mass;
			}
      
			++j;
      
      if(flag_field_gal_on){
        float sigma = 126*pow(mass*field_galaxy_mass_fraction/1.0e10,0.25); // From Tully-Fisher and Bell & de Jong 2001
        
        //std::cout << "Warning: All galaxies are spherical" << std::endl;
        
        //****** test lines taken out
        //field_galaxy_mass_fraction *= 2;

        float fratio = (ran2(seed)+1)*0.5;  //TODO: Ben change this!  This is a kluge.
        float pa = 2*pi*ran2(seed);  //TODO: This is a kluge.
        
 				switch(field_int_prof_gal_type){
          case null_gal:
            ERROR_MESSAGE();
            std::cout << "flag_field_gal_on is true, but field_int_prof_gal_type is null!!!!" << std::endl;
            break;
          case nsie_gal:
            field_halos.push_back(new LensHaloRealNSIE(mass*field_galaxy_mass_fraction,z,sigma,0.0,fratio,pa,0));
            //std::cout << sigma << std::endl;
            break;
          case pl_gal:
            assert(field_int_prof_gal_slope>0);
          
            field_halos.push_back(new LensHaloPowerLaw(mass*field_galaxy_mass_fraction,rmaxNSIE(sigma, mass*field_galaxy_mass_fraction, fratio, 0.0),z,field_int_prof_gal_slope,fratio,pa+pi/2.,0,Fourier));
            //field_halos.push_back(new LensHaloPowerLaw(mass*field_galaxy_mass_fraction,r_half_stel_mass/1.6*2.5,z,field_int_prof_gal_slope,0.99,pa+pi/2.,0,Fourier)); // explanation for r_half_stel_mass/1.34: relation between r_half_stel_mass and effective radius according to Kravtsev 2013 used!
            //field_halos.push_back(new LensHaloPowerLaw(mass*field_galaxy_mass_fraction,mass*Grav*lightspeed*lightspeed*sqrt(fratio)/pi/sigma/sigma,z,field_int_prof_gal_slope,fratio,pa,0,Fourier));
            
            //std::cout << "PL "<<r_half_stel_mass/1.34 << std::endl;
            break;
          case hern_gal:
            field_halos.push_back(new LensHaloHernquist(mass*field_galaxy_mass_fraction,rmaxNSIE(sigma, mass*field_galaxy_mass_fraction, fratio, 0.0),z,1,fratio,pa,0,Pseudo));
            break;
          case jaffe_gal:
            field_halos.push_back(new LensHaloJaffe(mass*field_galaxy_mass_fraction,rmaxNSIE(sigma, mass*field_galaxy_mass_fraction,fratio,0.0),z,1,fratio,pa,0,Pseudo));
            break;
            
          default:
            throw std::runtime_error("Don't support any but NSIE, PowerLaw, Hernquist and Jaffe galaxies yet!");
            break;
				}
        
        //field_halos[j]->setZlens(z);
				//field_halos[j]->initFromFile(mass*field_galaxy_mass_fraction,seed,vmax,r_halfmass*cosmo.gethubble());
        
        // Another copy of this position must be made to avoid rescaling it twice when it is converted into
        // distance on the lens plane in Lens::buildLensPlanes()
        //theta2 = new PosType[2];
        //theta2[0]=theta[0]; theta2[1]=theta[1];
				//halo_pos_vec.push_back(theta2);
        field_halos.back()->setX(theta);
        field_halos.back()->setID(haloid);
        
        /****** test **********
        {
          PosType tmpx[2];
          field_halos.back()->getX(tmpx);
          std::cout << "gal " << tmpx[0] << "  " << tmpx[1] << " " << field_halos.back()->get_mass()
          << " " << field_halos.back()->get_Rmax() << " " << sigma << " " << fratio << " " << field_halos.back()->getID() << std::endl;
        }*/
				++j;
			}
      
		}
	}
	file_in.close();
	if(verbose) std::cout << field_halos.size() << " halos read in."<< std::endl
    << "Max input mass = " << mass_max << "  R max = " << R_max << "  V max = " << V_max
    << "Min input mass = " << minmass << std::endl;
  
	/// setting the minimum halo mass in the simulation
	field_min_mass = minmass;
	if(field_buffer > 0.0){
		if(verbose) std::cout << "Overiding field_buffer to make it 0 because halos are read in." << std::endl;
		field_buffer = 0.0;
	}
  
	//halo_pos = Utilities::PosTypeMatrix(field_halos.size(), 3);
  
	//for(i = 0; i < field_halos.size(); ++i)
	//{
	//	halo_pos[i] = halo_pos_vec[i];
	//}
  
	std::cout << "Overiding input file field of view to make it fit the simulation light cone." << std::endl;
	fieldofview = pi*rmax*rmax*pow(180/pi,2);  // Resets field of view to range of input galaxies
	std::cout << "    It is now " << fieldofview << " deg^2" << std::endl;
  
	if(verbose) std::cout << "Setting mass function to Sheth-Tormen." << std::endl;
	field_mass_func_type = ST; // set mass function
  
	if(verbose) std::cout << "sorting in Lens::readInputSimFileMillennium()" << std::endl;
	// sort the field_halos by readshift
	//Lens::quicksort(field_halos.data(),halo_pos,field_halos.size());
  std::sort(field_halos.begin(),field_halos.end(),LensHaloZcompare);
  
	if(verbose) std::cout << "leaving Lens::readInputSimFileMillennium()" << std::endl;
  
  field_buffer = 0.0;
	read_sim_file = true;
 	delete halo_calc;
  
}

/**
 * \brief Read in information from a MultiDark Halo Catalog
 
 */
void Lens::readInputSimFileMultiDarkHalos(bool verbose)
{
  std::cout << "Reading Field Halos from " << field_input_sim_file << std::endl;
	PosType z,zob,xpos,ypos,zpos,vx,vy,vz,mass;
	unsigned long i,j;
	const PosType mo=7.3113e10,M1=2.8575e10,gam1=7.17,gam2=0.201,be=0.557;
  PosType field_galaxy_mass_fraction = 0;
//  const PosType masslimit =2.0e12;
  const PosType masslimit = 0.0;
  
  
  Utilities::Geometry::SphericalPoint tmp_sph_point(1,0,0);
  
	PosType rmax=0,rtmp=0,boundary_p1[2],boundary_p2[2],boundary_diagonal[2];
  
  std::vector<std::string> filenames;
  //Utilities::ReadFileNames(field_input_sim_file.c_str(),".dat",filenames);
  filenames.push_back(" ");
  //std::vector<PosType *> halo_pos_vec;
  //std::vector<Utilities::Geometry::SphericalPoint> sph_halo_pos_vec;
  
  int j_max;
	PosType mass_max=0,R_max=0,minmass=1e30;
	PosType theta[2];
  PosType center[2] = {0.0,0.0};
	int ncolumns = 11;
  size_t haloid = 0;
  
  void *addr[ncolumns];
  
	//addr[0] = &ra;
	//addr[1] = &dec;
  
  addr[0] = &(tmp_sph_point.phi);
	addr[1] = &(tmp_sph_point.theta);
  
	addr[2] = &z;
	addr[3] = &zob;
	addr[4] = &xpos;
	addr[5] = &ypos;
	addr[6] = &zpos;
	addr[7] = &vx;
	addr[8] = &vy;
	addr[9] = &vz;
	addr[10] = &mass;
  
  //for(int jj=0;jj<filenames.size();++jj){
  for(int jj=0;jj<1;++jj){
    
    //std::ifstream file_in( field_input_sim_file + filenames[jj].c_str());
    std::ifstream file_in( field_input_sim_file );
    if(!file_in){
      std::cout << "Can't open file " << field_input_sim_file + filenames[jj] << std::endl;
      ERROR_MESSAGE();
      throw std::runtime_error(" Cannot open file.");
      exit(1);
    }
    
    std::cout << "reading halos from file: " << field_input_sim_file + filenames[jj] << std::endl;
    
    // skip through header information in data file
    i=0;
    while(file_in.peek() == '#'){
      file_in.ignore(10000,'\n');
      ++i;
    }
    if(verbose) std::cout << "skipped "<< i << " comment lines in file " <<  filenames[jj]
      << std::endl;
    
    // read in data
    
    PosType myPosType;
    std::string myline;
    std::string strg;
    std::string f=" ";
    std::stringstream buffer;
    
    i=j=0;
    myline.clear();
    while(std::getline(file_in,myline)){
      
      int pos = myline.find_first_not_of(f);
      myline.erase(0,pos);
      
      for(int l=0;l<ncolumns; l++){
        pos = myline.find(f);
        strg.assign(myline,0,pos);
        buffer << strg;
        buffer >> myPosType;
        *((PosType *)addr[l]) = myPosType;
        
        myline.erase(0,pos);
        strg.clear();
        buffer.clear();
        buffer.str(std::string());
        
        pos = myline.find_first_not_of(f);
        myline.erase(0,pos);
      }
      ++haloid;
      
      mass = pow(10,mass)/cosmo.gethubble();
      if(mass > masslimit && z <= zsource){
        
        tmp_sph_point.theta *= pi/180;
        tmp_sph_point.phi *= -pi/180;
        
        rtmp = Utilities::Geometry::AngleSeporation(central_point_sphere,tmp_sph_point);
        
        if(sim_angular_radius > 0 && rtmp > sim_angular_radius) continue;
        
        if(rtmp > rmax) rmax = rtmp;
        
        // position on lens plane
        //theta = new PosType[2];
        tmp_sph_point.OrthographicProjection(central_point_sphere,theta);
        
        if(field_halos.size() > 0 ){
          if(boundary_p1[0] > theta[0]) boundary_p1[0] = theta[0];
          if(boundary_p1[1] > theta[1]) boundary_p1[1] = theta[1];
          
          if(boundary_p2[0] < theta[0]) boundary_p2[0] = theta[0];
          if(boundary_p2[1] < theta[1]) boundary_p2[1] = theta[1];
          
          PosType tmp = theta[0]+theta[1];
          if(boundary_diagonal[0] > tmp) boundary_diagonal[0] = tmp;
          if(boundary_diagonal[1] < tmp) boundary_diagonal[1] = tmp;
        }else{
          boundary_p1[0] = boundary_p2[0] = theta[0];
          boundary_p1[1] = boundary_p2[1] = theta[1];
          boundary_diagonal[0] = boundary_diagonal[1] = theta[0]+theta[1];
        }
        
        //theta[0] = -ra*pi/180.;
        //theta[1] = dec*pi/180.;
        
        center[0] += tmp_sph_point.theta;
        center[1] += tmp_sph_point.phi;
        
        float sigma=0;
        
        if(flag_field_gal_on){
          // from Moster et al. 2010ApJ...710..903M
          field_galaxy_mass_fraction = mo*pow(mass/M1,gam1)
          /pow(1+pow(mass/M1,be),(gam1-gam2)/be)/mass;
          if(field_galaxy_mass_fraction > 1.0) field_galaxy_mass_fraction = 1;
          sigma = 126*pow(mass*(1-field_galaxy_mass_fraction)/1.0e10,0.25); // From Tully-Fisher and Bell & de Jong 2001
        }else{
          field_galaxy_mass_fraction = 0;
        }
        
        switch(field_int_prof_type)
        {
          case null_lens:
            ERROR_MESSAGE();
            std::cout << "field_int_prof_type is null!!!!" << std::endl;
            break;
          case nfw_lens:
            // calculate the average size and concentration of a NFW at the mass and reshift
            if(mass > 0){
              HALOCalculator hcalc(&cosmo,mass*(1-field_galaxy_mass_fraction),z);
              
              field_halos.push_back(new LensHaloNFW(mass*(1-field_galaxy_mass_fraction),hcalc.getRvir(),z,hcalc.getConcentration(),1.0,0.0,0));
            }
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
            field_halos.push_back(new LensHaloRealNSIE(mass*(1-field_galaxy_mass_fraction),z,sigma,0.0,1.0,0.0,0));
            
            //field_halos.push_back(new LensHaloRealNSIE);
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
        
        //if(mass > 0.0) halo_pos_vec.push_back(theta);
        if(mass > 0.0){
          field_halos.back()->setX(theta);
          field_halos.back()->setID(haloid);
        }
        if(mass > mass_max) {
          mass_max = mass;
          j_max = j;
        }
        if(mass < minmass) {
          minmass = mass;
        }
        
        ++j;
        
        if(flag_field_gal_on){
          float sigma = 126*pow(mass*field_galaxy_mass_fraction/1.0e10,0.25); // From Tully-Fisher and Bell & de Jong 2001
          //std::cout << "Warning: All galaxies are spherical" << std::endl;
          float fratio = (ran2(seed)+1)*0.5;  //TODO: Ben change this!  This is a kluge.
          float pa = 2*pi*ran2(seed);  //TODO: This is a kluge.
          
          switch(field_int_prof_gal_type){
            case null_gal:
              ERROR_MESSAGE();
              std::cout << "flag_field_gal_on is true, but field_int_prof_gal_type is null!!!!" << std::endl;
              break;
            case nsie_gal:
              field_halos.push_back(new LensHaloRealNSIE(mass*field_galaxy_mass_fraction,z,sigma,0.0,fratio,pa,0));
              break;
            default:
              throw std::runtime_error("Don't support any but NSIE galaxies yet!");
              break;
          }
          
          //field_halos[j]->setZlens(z);
          //field_halos[j]->initFromFile(mass*field_galaxy_mass_fraction,seed,vmax,r_halfmass*cosmo.gethubble());
          
          // Another copy of this position must be made to avoid rescaling it twice when it is converted into distance on the lens plane in Lens::buildLensPlanes()
          //theta2 = new PosType[2];
          //theta2[0]=theta[0]; theta2[1]=theta[1];
          //halo_pos_vec.push_back(theta2);
          field_halos.back()->setX(theta);
          field_halos.back()->setID(haloid);
          
          ++j;
        }
        
      }
      myline.clear();
      ++i;
    }
    if(verbose){
      std::cout << field_halos.size() << " halos read in."<< std::endl;
      std::cout << "center is : theta:" << center[0]/field_halos.size() << "  phi:" << center[1]/field_halos.size() << std::endl;
    }
    file_in.close();
  }
  
  
	if(verbose) std::cout << field_halos.size() << " halos read in."<< std::endl
    << "Max input mass = " << mass_max << "  R max = " << R_max
    << " Min imput mass = " << minmass << std::endl;
  
	/// setting the minimum halo mass in the simulation
	field_min_mass = minmass;
	if(field_buffer > 0.0){
		if(verbose) std::cout << "Overiding field_buffer to make it 0 because halos are read in." << std::endl;
		field_buffer = 0.0;
	}
  
	//halo_pos = Utilities::PosTypeMatrix(field_halos.size(), 3);
  
	//for(i = 0; i < field_halos.size(); ++i)
	//{
	//	halo_pos[i] = halo_pos_vec[i];
	//}
  
  // determine if the region is a circle or a rectangle
  //PosType diagonal1 = (boundary_diagonal[1] - boundary_diagonal[0])/sqrt(2);
  //PosType diagonal2 = sqrt(pow(boundary_p2[0] - boundary_p1[0],2) + pow(boundary_p2[1] - boundary_p1[1],2));
  
	std::cout << "Overiding input file field of view to make it fit the simulation light cone." << std::endl;
  rmax = (boundary_p2[0] - boundary_p1[0])/2;
  /*if(diagonal1 < diagonal2*0.9){
   // circular region
   rmax = diagonal1/2;
   fieldofview = pi*rmax*rmax*pow(180/pi,2);  // Resets field of view to range of input galaxies
   inv_ang_screening_scale = 0.0;
   }else{
   fieldofview = (boundary_p2[0] - boundary_p1[0])*(boundary_p2[1] - boundary_p1[1])*pow(180/pi,2);
   rmax = diagonal2/2;
   inv_ang_screening_scale = 5.0/(MIN(boundary_p2[0] - boundary_p1[0],boundary_p2[1] - boundary_p1[1]));
   }*/
  
  if(sim_angular_radius == 0.0){
    fieldofview = (boundary_p2[0] - boundary_p1[0])*(boundary_p2[1] - boundary_p1[1])*pow(180/pi,2);
    inv_ang_screening_scale = 5.0/(MIN(boundary_p2[0] - boundary_p1[0],boundary_p2[1] - boundary_p1[1]));
  }else{
    fieldofview = pi*rmax*rmax*pow(180/pi,2);
    inv_ang_screening_scale = 0.0;
  }
  
	if(verbose) std::cout << "Setting mass function to Sheth-Tormen." << std::endl;
	field_mass_func_type = ST; // set mass function
  
	if(verbose) std::cout << "sorting in Lens::readInputSimFileMultiDarkHalos()" << std::endl;
	// sort the field_halos by readshift
  std::sort(field_halos.begin(),field_halos.end(),LensHaloZcompare);
  
  
	if(verbose) std::cout << "leaving Lens::readInputSimFileMultiDarkHalos()" << std::endl;
  
  field_buffer = 0.0;
	read_sim_file = true;
}

/**
 * \brief Read in information from a catalog of observed galaxies
 *
 * The LensHalo id numbers are set to the id of the parent halo in the simulation.  In this case the id 
 * numbers are the order they appear in the halo catalog.
 */
void Lens::readInputSimFileObservedGalaxies(bool verbose)
{
  std::cout << "Reading Field Halos from " << field_input_sim_file << std::endl;
  PosType z,mass,rcut,vdist;
  unsigned long i,j;
  //const PosType mo=7.3113e10,M1=2.8575e10,gam1=7.17,gam2=0.201,be=0.557;
  PosType field_galaxy_mass_fraction = 0;
  //const PosType masslimit =2.0e12;

  Utilities::Geometry::SphericalPoint tmp_sph_point(1,0,0);
  
  PosType rmax=0,rtmp=0,boundary_p1[2],boundary_p2[2],boundary_diagonal[2];
  
  //Utilities::ReadFileNames(field_input_sim_file.c_str(),".dat",filenames);
  //filenames.push_back(" ");
  //std::vector<PosType *> halo_pos_vec;
  //std::vector<Utilities::Geometry::SphericalPoint> sph_halo_pos_vec;
  
  int j_max;
  PosType mass_max=0,R_max=0,minmass=1e30;
  PosType *theta;
  PosType center[2] = {0.0,0.0};
  int ncolumns = 6;
  size_t haloid = 0;
  
  void *addr[ncolumns];
  
  //addr[0] = &ra;
  //addr[1] = &dec;
 
  addr[0] = &haloid;
  addr[1] = &(tmp_sph_point.phi);
  addr[2] = &(tmp_sph_point.theta);
  addr[3] = &z;
  addr[4] = &rcut;
  addr[5] = &vdist;
  
  std::ifstream file_in( field_input_sim_file );
  if(!file_in){
    std::cout << "Can't open file " << field_input_sim_file << std::endl;
    ERROR_MESSAGE();
    throw std::runtime_error(" Cannot open file.");
    exit(1);
  }
  
  std::cout << "reading halos from file: " << field_input_sim_file << std::endl;
  
  // skip through header information in data file
  i=0;
  while(file_in.peek() == '#'){
    file_in.ignore(10000,'\n');
    ++i;
  }
  if(verbose) std::cout << "skipped "<< i << " comment lines in file "
    << std::endl;
  
  // read in data
  
  PosType myPosType;
  std::string myline;
  std::string strg;
  std::string f=" ";
  std::stringstream buffer;
  
  i=j=0;
  myline.clear();
  while(std::getline(file_in,myline)){
    
    int pos = myline.find_first_not_of(f);
    myline.erase(0,pos);
    
    for(int l=0;l<ncolumns; l++){
      //std::cout << myline << std::endl;
      
      pos = myline.find(f);
      strg.assign(myline,0,pos);
      buffer << strg;
      buffer >> myPosType;
      *((PosType *)addr[l]) = myPosType;
      
      myline.erase(0,pos);
      strg.clear();
      buffer.clear();
      buffer.str(std::string());
      
      pos = myline.find_first_not_of(f);
      myline.erase(0,pos);
    }
    
    for(int l=0;l<ncolumns; l++) std::cout << *((PosType *)(addr[l])) << "  ";
    std::cout << std::endl;
      
    tmp_sph_point.theta *= pi/180;
    tmp_sph_point.phi *= pi/180;
    rcut *= 1.0e-3;
    
    rtmp = Utilities::Geometry::AngleSeporation(central_point_sphere,tmp_sph_point);
    
    if(sim_angular_radius > 0 && rtmp > sim_angular_radius) continue;
    
    if(rtmp > rmax) rmax = rtmp;
    
    // position on lens plane
    theta = new PosType[2];
    tmp_sph_point.OrthographicProjection(central_point_sphere,theta);
    
    if(field_halos.size() > 0 ){
      if(boundary_p1[0] > theta[0]) boundary_p1[0] = theta[0];
      if(boundary_p1[1] > theta[1]) boundary_p1[1] = theta[1];
      
      if(boundary_p2[0] < theta[0]) boundary_p2[0] = theta[0];
      if(boundary_p2[1] < theta[1]) boundary_p2[1] = theta[1];
      
      PosType tmp = theta[0]+theta[1];
      if(boundary_diagonal[0] > tmp) boundary_diagonal[0] = tmp;
      if(boundary_diagonal[1] < tmp) boundary_diagonal[1] = tmp;
    }else{
      boundary_p1[0] = boundary_p2[0] = theta[0];
      boundary_p1[1] = boundary_p2[1] = theta[1];
      boundary_diagonal[0] = boundary_diagonal[1] = theta[0]+theta[1];
    }
    
    //theta[0] = -ra*pi/180.;
    //theta[1] = dec*pi/180.;
    
    center[0] += tmp_sph_point.theta;
    center[1] += tmp_sph_point.phi;
    
    /*
     if(flag_field_gal_on){
     // from Moster et al. 2010ApJ...710..903M
     field_galaxy_mass_fraction = mo*pow(mass/M1,gam1)
     /pow(1+pow(mass/M1,be),(gam1-gam2)/be)/mass;
     if(field_galaxy_mass_fraction > 1.0) field_galaxy_mass_fraction = 1;
     sigma = 126*pow(mass*(1-field_galaxy_mass_fraction)/1.0e10,0.25); // From Tully-Fisher and Bell & de Jong 2001
     }else{
     field_galaxy_mass_fraction = 0;
     }
     */
    
    switch(field_int_prof_type)
    {
      case null_lens:
        ERROR_MESSAGE();
        std::cout << "field_int_prof_type is null!!!!" << std::endl;
        break;
      case nfw_lens:
        // calculate the average size and concentration of a NFW at the mass and reshift

        std::cout << "Making an NFW halo, but not of the right mass!!!!" << std::endl;
      {
        mass = pi*vdist*vdist*rcut/Grav/lightspeed/lightspeed/field_galaxy_mass_fraction;
        
        HALOCalculator hcalc(&cosmo,mass,z);
        
        field_halos.push_back(new LensHaloNFW(mass,hcalc.getRvir(),z,hcalc.getConcentration(),1.0,0.0,0));
      }
        break;
      case pnfw_lens:
        ERROR_MESSAGE();
        std::cout << "PseudoNFW not supported." << std::endl;
        break;
      case pl_lens:
        ERROR_MESSAGE();
        std::cout << "PowerLaw not supported." << std::endl;
        break;
      case nsie_lens:
        
        mass = pi*vdist*vdist*rmax/Grav/lightspeed/lightspeed;
        field_halos.push_back(new LensHaloRealNSIE(mass,z,vdist,0.0,1.0,0.0,0));
        
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
    
    //if(mass > 0.0) halo_pos_vec.push_back(theta);
    if(mass > 0.0){
      field_halos.back()->setX(theta);
      field_halos.back()->setID(haloid);
    }
    if(mass > mass_max) {
      mass_max = mass;
      j_max = j;
    }
    if(mass < minmass) {
      minmass = mass;
    }
    
    ++j;
    
    /*
     if(flag_field_gal_on){
     float sigma = 126*pow(mass*field_galaxy_mass_fraction/1.0e10,0.25); // From Tully-Fisher and Bell & de Jong 2001
     //std::cout << "Warning: All galaxies are spherical" << std::endl;
     float fratio = (ran2(seed)+1)*0.5;  //TODO: Ben change this!  This is a kluge.
     float pa = 2*pi*ran2(seed);  //TODO: This is a kluge.
     
     switch(field_int_prof_gal_type){
     case null_gal:
     ERROR_MESSAGE();
     std::cout << "flag_field_gal_on is true, but field_int_prof_gal_type is null!!!!" << std::endl;
     break;
     case nsie_gal:
     field_halos.push_back(new LensHaloRealNSIE(mass*field_galaxy_mass_fraction,z,sigma,0.0,fratio,pa,0));
     break;
     default:
     throw std::runtime_error("Don't support any but NSIE galaxies yet!");
     break;
     }
     
     //field_halos[j]->setZlens(z);
     //field_halos[j]->initFromFile(mass*field_galaxy_mass_fraction,seed,vmax,r_halfmass*cosmo.gethubble());
     
     // Another copy of this position must be made to avoid rescaling it twice when it is converted into
     // distance on the lens plane in Lens::buildLensPlanes()
     theta2 = new PosType[2];
     theta2[0]=theta[0]; theta2[1]=theta[1];
     //halo_pos_vec.push_back(theta2);
     field_halos.back()->setX(theta2);
     field_halos.back()->setID(haloid);
     
     ++j;
     }
     */
    
    myline.clear();
    ++i;
  }
  if(verbose){
    std::cout << field_halos.size() << " halos read in."<< std::endl;
    std::cout << "center is at: theta:" << center[0]*180/pi/field_halos.size() << "  phi:" << center[1]*180/pi/field_halos.size() << " degrees" << std::endl;
  }
  file_in.close();
  
  
  
  if(verbose) std::cout << field_halos.size() << " halos read in."<< std::endl
    << "Max input mass = " << mass_max << "  R max = " << R_max
    << " Min input mass = " << minmass << std::endl;
  
  /// setting the minimum halo mass in the simulation
  field_min_mass = minmass;
  if(field_buffer > 0.0){
    if(verbose) std::cout << "Overiding field_buffer to make it 0 because halos are read in." << std::endl;
    field_buffer = 0.0;
  }
  
  //halo_pos = Utilities::PosTypeMatrix(field_halos.size(), 3);
  
  //for(i = 0; i < field_halos.size(); ++i)
  //{
  //	halo_pos[i] = halo_pos_vec[i];
  //}
  
  // determine if the region is a circle or a rectangle
  //PosType diagonal1 = (boundary_diagonal[1] - boundary_diagonal[0])/sqrt(2);
  //PosType diagonal2 = sqrt(pow(boundary_p2[0] - boundary_p1[0],2) + pow(boundary_p2[1] - boundary_p1[1],2));
  
  std::cout << "Overiding input file field of view to make it fit the simulation light cone." << std::endl;
  rmax = (boundary_p2[0] - boundary_p1[0])/2;
  
  fieldofview = (boundary_p2[0] - boundary_p1[0])*(boundary_p2[1] - boundary_p1[1])*pow(180/pi,2);
  std::cout << "Field of view is " << fieldofview << " sq.deg." << std::endl;
  inv_ang_screening_scale = 0.0;

  
  if(verbose) std::cout << "sorting in Lens::readInputSimFileObservedGalaxies()" << std::endl;
  // sort the field_halos by readshift
  std::sort(field_halos.begin(),field_halos.end(),
            [](LensHalo *lh1,LensHalo *lh2){return (lh1->getZlens() < lh2->getZlens());});
  
  
  if(verbose) std::cout << "leaving Lens::readInputSimFileObservedGalaxies()" << std::endl;
  
  field_buffer = 0.0;
  read_sim_file = true;
}

void Lens::combinePlanes(bool verbose)
{
  if(verbose)
  {
    std::cout << std::endl << "Lens::combinePlanes before clearing." << std::endl ;

    if(plane_redshifts.size() > 0)
    {
      std::cout << "Lens::combinePlanes : plane_redshifts[] = " ;
      for(int i=0;i<plane_redshifts.size();i++) std::cout << plane_redshifts[i] << " " ;
      std::cout << std::endl ;
    }
    else std::cout << "Lens::combinePlanes : plane_redshifts[] size is zero." << std::endl ;

    if(Dl.size() > 0)
    {
      std::cout << "Lens::combinePlanes : Dl[] = " ;
      for(int i=0;i<Dl.size();i++) std::cout << Dl[i] << " " ;
      std::cout << std::endl ;
    }
    else std::cout << "Lens::combinePlanes : Dl[] size is zero." << std::endl ;
    
    if(dDl.size() > 0)
    {
      std::cout << "Lens::combinePlanes : dDl[] = " ;
      for(int i=0;i<dDl.size();i++) std::cout << dDl[i] << " " ;
      std::cout << std::endl << std::endl ;
    }
    else std::cout << "Lens::combinePlanes : dDl[] size is zero." << std::endl ;

    if(field_plane_redshifts.size() > 0)
    {
      std::cout << "Lens::combinePlanes : field_plane_redshifts[] = " ;
      for(int i=0;i<field_plane_redshifts.size();i++) std::cout << field_plane_redshifts[i] << " " ;
      std::cout << std::endl ;
    }
    else std::cout << "Lens::combinePlanes : field_plane_redshifts[] size is zero." << std::endl ;
    
    if(main_plane_redshifts.size() > 0)
    {
      std::cout << "Lens::combinePlanes : main_plane_redshifts[] = " ;
      for(int i=0;i<main_plane_redshifts.size();i++) std::cout << main_plane_redshifts[i] << " " ;
      std::cout << std::endl ;
    }
    else std::cout << "Lens::combinePlanes : main_plane_redshifts[] size is zero." << std::endl ;
  }
  
  std::cout << "Lens::combinePlanes : field_planes.size() = " << field_planes.size() << " , main_planes.size() = " << main_planes.size() << std::endl ;
  
  // clear old plane configuration
  lensing_planes.clear();
  plane_redshifts.clear();
  Dl.clear();
  dDl.clear();

  // copy main and field planes into master Dl
  Dl.resize(field_Dl.size() + main_Dl.size());
  int i;
  for(i =0; i < field_Dl.size(); ++i){
    Dl[i] = field_Dl[i];
  }
  int j = 0;
  while(i < Dl.size()){
    Dl[i++] = main_Dl[j++];
  }
  
  if(verbose && Dl.size() > 0)
  {
    std::cout << "Lens::combinePlanes : before sorting : Dl master = " ;
    for(int i=0;i<Dl.size();i++) std::cout << std::setprecision(13) << Dl[i] << " " ;
    std::cout << std::endl ;
  }
        
  // sort these Dl and keep an index so that the identity of each Dl is remembered
  std::vector<size_t> index(Dl.size());
  Utilities::sort_indexes(Dl,index);
  std::sort(Dl.begin(),Dl.end());
  
  // changing the position of the planes too close to each other
  for(int i=1; i < Dl.size(); ++i){
    if( (Dl[i] - Dl[i-1]) < MIN_PLANE_DIST) Dl[i] = Dl[i-1] + MIN_PLANE_DIST;
  }

  if(verbose && Dl.size() > 0)
  {
    std::cout << "Lens::combinePlanes : after sorting and adjusting distances : Dl master = " ;
    for(int i=0;i<Dl.size();i++) std::cout << std::setprecision(13) << Dl[i] << " " ;
    std::cout << std::endl ;
  }
  
  // filling the tables for redshift and lensing planes
  for(auto i : index){
    if(i<field_Dl.size()){
      plane_redshifts.push_back(field_plane_redshifts[i]);
      lensing_planes.push_back(field_planes[i]);
    }else{
      plane_redshifts.push_back(main_plane_redshifts[i - field_Dl.size()]);
      lensing_planes.push_back(main_planes[i - field_Dl.size()]);
    }
  }
  
  if(verbose && plane_redshifts.size() > 0)
  {
    std::cout << "Lens::combinePlanes : redshifts = " ;
    for(int i=0;i<plane_redshifts.size();i++) std::cout << std::setprecision(13) << plane_redshifts[i] << " " ;
    std::cout << std::endl ;
  }
  
  assert(lensing_planes.size() == field_planes.size() + main_planes.size());
  // std::cout << "assert : " << zsource << " , " << plane_redshifts.back() << std::endl ;
  //assert(zsource > plane_redshifts.back()); // !!!
  
  if(plane_redshifts.size() > 0 && zsource <= plane_redshifts.back()){
    zsource = plane_redshifts.back() + 0.1;
  }

  // add the pseudo-plane for rayshooting at the end of the arrays
  plane_redshifts.push_back(zsource);
  Dl.push_back(cosmo.coorDist(0, zsource));
  
  // calculate deltas
  dDl.push_back(Dl[0]);
  for(std::size_t i = 1; i < Dl.size(); ++i)
  dDl.push_back(Dl[i] - Dl[i-1]); // distance from plane i-1 to plane i
  
  // output resulting setup
  if(verbose)
  {
    std::cout << "\ncombinePlanes()" << "\n---------------" << std::endl;
    std::cout << "\nz:";
    for(std::size_t i = 0, n = plane_redshifts.size(); i < n; ++i) std::cout << " " << plane_redshifts[i];
    std::cout << "\nDl:";
    for(std::size_t i = 0, n = Dl.size(); i < n; ++i) std::cout << std::setprecision(12) << " " << Dl[i];
    std::cout << "\ndDl:";
    for(std::size_t i = 0, n = dDl.size(); i < n; ++i) std::cout << " " << dDl[i];
    std::cout << "\n" << std::endl;
  }
  
}



void Lens::buildPlanes(InputParams& params, bool verbose)
{
	// build field
	if(!flag_switch_field_off)
	{
		// set the distances of the field planes
		setupFieldPlanes();

		// create or read the field halos
		if(sim_input_flag){
      if(field_input_sim_format == MillenniumObs) readInputSimFileMillennium(verbose);
      if(field_input_sim_format == MultiDarkHalos) readInputSimFileMultiDarkHalos(verbose);
      if(field_input_sim_format == ObservedData) readInputSimFileObservedGalaxies(verbose);
    }
    else{
      createFieldHalos(verbose);
		}
    // create field planes and sort halos onto them
		createFieldPlanes(verbose);
	}
  
	// build main
	if(flag_switch_main_halo_on || pixel_map_on)
	{
		// create the main halos
		createMainHalos(params);
		
		// create the main planes for the halos
		createMainPlanes();
	}
	
	// combine the different planes
	combinePlanes(verbose);
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
                             PosType z                 /// redshift of implanted source
                             ,bool nearest           /** If true, set the source plane to the nearest (in coordinate distance)
                                                      * lensing plane that was created already.  This can be used to avoid self-lensing
                                                      * by the halo of the source.  If the source is at higher redshift than the simulation
                                                      * volume the source will be at its real redshift.
                                                      */
                             ,unsigned long GalID
                             ,PosType *xx
                             ,bool verbose
                             ){
	unsigned long j;
	short out=0;
  
	toggle_source_plane = true;
  
	if(z<=0.0){
		std::cout << "Warning: Source redshift can't be set to " << z << " in MultiLens::ResetSourcePlane." << endl;
		return out;
	}
  
  
	// distance to new source plane
	PosType Ds = cosmo.coorDist(0,z);
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
			PosType z1 = cosmo.invCoorDist(Dl[j]-0.5*dDl[j]);
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
		if(j > 0) dDs_implant = cosmo.coorDist(plane_redshifts[j-1],z); //Ds - Dl[j-1];
		else  dDs_implant = Ds;
	}
  
	if(verbose) std::cout << "Source on plane " << j << " zs " << zs_implant << " Ds " << Ds << " dDs " << dDs_implant << std::endl;
  
	index_of_new_sourceplane = j;
  
	out = j;
	return out;
}

/// Sort field_halos[], brr[][], and id[] by content off arr[]
void Lens::quicksort(LensHaloHndl *halos,PosType **pos,unsigned long N){
	PosType pivotvalue;
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
