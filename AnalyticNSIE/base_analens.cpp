/*
 * base_analens.cpp
 *
 *  Created on: Jan 22, 2013
 *      Author: mpetkova
 */

#include "slsimlib.h"

using namespace std;

/**
 * \brief Reads in a parameter file and sets up an analytic lens.
 *
 * Sets many parameters within the lens model, source model and
 * force calculation.
 */
void BaseAnaLens::assignParams(InputParams& params){

	// Host lens parameters
	if(!params.get("z_lens",zlens)) error_message1("z_lens",params.filename());


	// Distortion of host lens parameters
	if(!params.get("NDistortionModes",perturb_Nmodes)) error_message1("NDistortionModes",params.filename());

    // Substructure parameters
    if(!params.get("NdensitySubstruct",sub_Ndensity)) error_message1("NdensitySubstruct",params.filename());
    else if(sub_Ndensity > 0){
    	if(!params.get("beta_sub",sub_beta)) error_message1("beta_sub",params.filename());
    	if(!params.get("alpha_sub",sub_alpha)) error_message1("alpha_sub",params.filename());
    	if(!params.get("R_submax",sub_Rmax)) error_message1("R_submax",params.filename());
    	if(!params.get("mass_max",sub_Mmax)) error_message1("mass_max",params.filename());
    	if(!params.get("mass_min",sub_Mmin)) error_message1("mass_min",params.filename());
    	if(sub_Mmin < 1.0e3){
    		ERROR_MESSAGE();
    		std::cout << "Are you sure the minimum halo mass should be " << sub_Mmin << " Msun?" << std::endl;
    		exit(1);
    	}
    	if(!params.get("sub_type",sub_type)) error_message1("sub_type",params.filename());
    }
	  // Stars parameters
    if(!params.get("Nstars",stars_N)) error_message1("Nstars",params.filename());
    else if(stars_N){
    	if(!params.get("fstars",star_fstars)) error_message1("fstars",params.filename());
    	if(!params.get("stars_mass",star_massscale)) error_message1("stars_mass",params.filename());
    }

}

void BaseAnaLens::error_message1(std::string parameter,std::string file){
		  ERROR_MESSAGE();
		  std::cout << "Parameter " << parameter << " is needed to construct a BaseAnaLens.  It needs to be set in parameter file " << file << "!" << endl;
		  exit(0);
}


double BaseAnaLens::getZlens(){
	return zlens;
}

void BaseAnaLens::setZlens(double z){
	zlens = z;
}

void BaseAnaLens::reNormSubstructure(double kappa_sub){
	/* renomalizes substructure so that
	 * the average surface density it kappa_sub
	 */
	  double avem;
	  avem=sub_Mmax*(sub_alpha+1)
	    /(sub_alpha+2)*(1-pow(sub_Mmin/sub_Mmax,sub_alpha+2))/
	    (1-pow(sub_Mmin/sub_Mmax,sub_alpha+1));

	  sub_Ndensity=kappa_sub*Sigma_crit/avem;

	  return ;
}

void BaseAnaLens::setInternalParams(CosmoHndl cosmo, SourceHndl source){
	double Ds, Dls;

	Dl = cosmo->angDist(0,zlens);
	Ds = cosmo->angDist(0,source->getZ());
	Dls = cosmo->angDist(zlens,source->getZ());

	MpcToAsec = 60*60*180 / pi / Dl;
		// in Mpc
	host_ro=4*pi*pow(host_sigma/2.99792e5,2)*Dl
		*Dls/Ds;
	// find critical density
	Sigma_crit=Ds/Dls/Dl/4/pi/Grav;
	to = (1+zlens)*Ds/Dls/Dl/8.39428142e-10;
}

BaseAnaLens::BaseAnaLens(InputParams& params) : Lens(){

  perturb_rms = new double[6];

  assignParams(params);

  if(perturb_Nmodes){
  	perturb_modes = new double[perturb_Nmodes+1];
  	// zero perturbation modes until use BaseAnaLens::RandomlyDistortLens()
  	for(int i=0;i< perturb_Nmodes+1 ;++i) perturb_modes[i] =  0;
  }

  if(sub_Ndensity > 0){
  	switch(sub_type){
  	case nfw:
		  sub_alpha_func = alphaNFW;
		  sub_kappa_func = kappaNFW;
		  sub_gamma_func = gammaNFW;
		  sub_phi_func = 0;
		  ERROR_MESSAGE();
		  break;
  	case powerlaw:
		  sub_alpha_func = alphaPowLaw;
		  sub_kappa_func = kappaPowLaw;
		  sub_gamma_func = gammaPowLaw;
		  sub_phi_func = phiPowLaw;
		  break;
  	case pointmass:
		  sub_alpha_func = NULL;
		  sub_kappa_func = NULL;
		  sub_gamma_func = NULL;
		  sub_phi_func = NULL;
		  break;
  	default:
		  ERROR_MESSAGE();
		  cout << "ERROR: no submass internal profile chosen" << endl;
		  exit(1);
		  break;
  	}
  }

  // parameters for stars
  stars_implanted = false; // stars are implanted later
  star_theta_force = 0.1;
  sub_theta_force = 0.1;

  sub_sigmaScale = host_sigma = host_pos_angle = host_ro = host_axis_ratio = host_core = 0.0;

  if(sub_Ndensity == 0)
	  sub_N = 0;

  Sigma_crit = 0;

  stars_implanted = false;
  set = true;
}


BaseAnaLens::~BaseAnaLens(){
	cout << "deleting lens" << endl;

	delete[] perturb_rms;

	if(perturb_Nmodes > 0){
		cout << "deleting modes" << endl;
		delete[] perturb_modes;
	}
	if(sub_N > 0 && substruct_implanted){
		cout << "deleting subs" << endl;
		free_PosTypeMatrix(sub_x,sub_N,2);
		delete[] sub_Rcut;
		delete[] sub_mass;
		delete[] sub_substructures;
	}
	if(stars_N > 0 && stars_implanted){
		cout << "deleting stars" << endl;
		delete[] star_masses;
		delete[] stars;
		free_PosTypeMatrix(stars_xp,stars_N,3);
		delete[] star_region;
		delete[] star_kappa;
		free_PosTypeMatrix(star_xdisk,star_Nregions,2);
		delete star_tree;
	}
}

