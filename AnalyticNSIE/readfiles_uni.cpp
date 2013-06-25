/*
 * readfiles_uni.cpp
 *
 *  Created on: Jan 22, 2013
 *      Author: Leier
 */


#include "slsimlib.h"

using namespace std;

/**
 * \brief Reads in a parameter file and sets up an uniform lens.
 *
 * Sets many parameters within the lens model, source model and
 * force calculation.
 */

LensHaloUniform::LensHaloUniform(InputParams& params) : LensHalo(params){

  assignParams(params);

  if(std::numeric_limits<float>::has_infinity){
    Rmax = std::numeric_limits<float>::infinity();
  }else{
    Rmax = std::numeric_limits<float>::max();
  }
  perturb_Nmodes=3;
  perturb_modes = new double[3];

  perturb_modes[0]=kappa_uniform;
  perturb_modes[1]=gamma_uniform[0];
  perturb_modes[2]=gamma_uniform[1];

  PrintLens(false,false);
}

LensHaloUniform::~LensHaloUniform(){

}

void LensHaloUniform::setInternalParams(CosmoHndl cosmo){
		Dl = cosmo->angDist(0,reference_z);
		Ds = cosmo->angDist(0,zlens);
		Dls = cosmo->angDist(zlens,reference_z);
		double norm_factor = 8*pi*Grav*Dls*Dl/Ds;
		kappa_uniform /= norm_factor;
		gamma_uniform[0] /= norm_factor;
		gamma_uniform[1] /= norm_factor;
}

void LensHaloUniform::force_halo(
		double *alpha     /// mass/Mpc
		,KappaType *kappa
		,KappaType *gamma
		,double *xcm
		,bool no_kappa
		,bool subtract_point /// if true contribution from a point mass is subtracted
		)
{
    double alpha_tmp[2];
     KappaType kappa_tmp = 0.0, gamma_tmp[3], dt = 0,tmp = 0;

     gamma_tmp[0] = gamma_tmp[1] = gamma_tmp[2] = 0.0;
     alpha_tmp[0] = alpha_tmp[1] = 0.0;

	 *kappa += lens_expand(perturb_beta,perturb_modes
			 ,perturb_Nmodes,xcm,alpha_tmp,gamma_tmp,&dt);

	 alpha[0] += alpha_tmp[0];
	 alpha[1] += alpha_tmp[1];

	      if(!no_kappa){
	    	  gamma[0] += gamma_tmp[0];
	    	  gamma[1] += gamma_tmp[1];
	      }

	      // add stars for microlensing
	      if(stars_N > 0 && stars_implanted){
	     	 force_stars(alpha,kappa,gamma,xcm,no_kappa);
	      }


}


void LensHaloUniform::assignParams(InputParams& params){

	//if(perturb_Nmodes > 0){
	if(!params.get("kappa_uniform",kappa_uniform)) error_message1("kappa_uniform",params.filename());
	if(!params.get("gamma_uniform_1",gamma_uniform[0])) error_message1("gamma_uniform_1",params.filename());
	if(!params.get("gamma_uniform_2",gamma_uniform[1])) error_message1("gamma_uniform_2",params.filename());

	if(!params.get("stellar_mass_function",imf_type)) error_message1("stellar_mass_function",params.filename());
	if(!params.get("min_mstar",min_mstar)) error_message1("min_mstar",params.filename());
	if(!params.get("max_mstar",max_mstar)) error_message1("max_mstar",params.filename());
	if(!params.get("bending_point",bend_mstar)) error_message1("bending_point",params.filename());
	if(!params.get("slope_1",lo_mass_slope)) error_message1("slope_1",params.filename());
	if(!params.get("slope_2",hi_mass_slope)) error_message1("slope_2",params.filename());

	return;
}

/** \ingroup ImageFinding
 * \brief Prints the parameters of the analytic lens to stdout
 */
void LensHaloUniform::PrintLens(bool show_substruct,bool show_stars){
	int i;


	// uni lens parameters only
	cout << endl << "**Host lens model**" << endl;
	// redshifts
	cout << "zlens " << zlens << endl;
	cout << "kappa " << kappa_uniform << endl;
	cout << "gamma " << gamma_uniform[0] << " " << gamma_uniform[1] << endl;

	if (stars_implanted) PrintStars(show_stars);

}


void LensHalo::implant_stars(double x, double y, unsigned long Nregions,long *seed, IMFtype type){

	if(Nregions <= 0) return;

	double ** centers = new double*[Nregions];
	for (int i = 0; i < Nregions; ++i)
		centers[i] = new double[2];
	centers[0][0] = x;
	centers[0][1] = y;
	implant_stars(centers,Nregions,seed,type);
}


