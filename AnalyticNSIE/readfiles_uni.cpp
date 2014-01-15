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

LensHaloUniform::LensHaloUniform(InputParams& params, const COSMOLOGY& cosmo, bool verbose)
: LensHalo()
{
  
  assignParams(params);
  
  if(std::numeric_limits<float>::has_infinity){
    Rmax = std::numeric_limits<float>::infinity();
  }else{
    Rmax = std::numeric_limits<float>::max();
  }
  perturb_Nmodes=3;
  perturb_modes = new PosType[3];
  
  setCosmology(cosmo);
  if(verbose) PrintLens(false,false);
}

LensHaloUniform::LensHaloUniform(InputParams& params, bool verbose): LensHalo(){
  
  assignParams(params);
  
  if(std::numeric_limits<float>::has_infinity){
    Rmax = std::numeric_limits<float>::infinity();
  }else{
    Rmax = std::numeric_limits<float>::max();
  }
  perturb_Nmodes=3;
  perturb_modes = new PosType[3];
  
  if(verbose) PrintLens(false,false);
}

LensHaloUniform::~LensHaloUniform(){

}

void LensHaloUniform::setCosmology(const COSMOLOGY& cosmo)
{
  PosType zlens = LensHalo::getZlens();
	Dl = cosmo.angDist(0,zlens);
	Ds = cosmo.angDist(0,zsource_reference);
	Dls = cosmo.angDist(zlens,zsource_reference);
  SigmaCrit = Ds/Dl/Dls/(4*pi*Grav);
	
  Sigma_uniform = kappa_uniform*SigmaCrit;
  gammaCrit_uniform[0] = gamma_uniform[0]*SigmaCrit;
  gammaCrit_uniform[1] = gamma_uniform[1]*SigmaCrit;
	
  perturb_modes[0] = Sigma_uniform;
  perturb_modes[1] = gammaCrit_uniform[0];
  perturb_modes[2] = gammaCrit_uniform[1];
}

void LensHaloUniform::force_halo(
		PosType *alpha     /// mass/Mpc
		,KappaType *kappa
		,KappaType *gamma
		,PosType *xcm
		,bool no_kappa
		,bool subtract_point /// if true contribution from a point mass is subtracted
		)
{
    PosType alpha_tmp[2];
     KappaType gamma_tmp[3], dt = 0;

     gamma_tmp[0] = gamma_tmp[1] = gamma_tmp[2] = 0.0;
     alpha_tmp[0] = alpha_tmp[1] = 0.0;

  
    *kappa += lens_expand(perturb_modes,xcm,alpha,gamma,&dt);

    // add stars for microlensing
    if(stars_N > 0 && stars_implanted){
      force_stars(alpha,kappa,gamma,xcm,no_kappa);
    }
}

PosType LensHaloUniform::lens_expand(PosType *mod,PosType *x,PosType *alpha,KappaType *gamma,KappaType *phi){
  PosType theta,r,cosx,sinx,cos2theta,sin2theta;

   // add shear
  alpha[0] +=  x[0]*mod[1] + x[1]*mod[2];
  alpha[1] += -x[1]*mod[1] + x[0]*mod[2];
  
  // add flat kappa
  alpha[0] += -1.0*x[0]*mod[0];
  alpha[1] += -1.0*x[1]*mod[0];
    
  
  gamma[0] += mod[1];
  gamma[1] += mod[2];
 
  r=sqrt(x[0]*x[0] + x[1]*x[1]);
  if(r == 0.0){
    *phi = 0.0;
  }else{
    theta=atan2(x[1],x[0]);
    cosx=x[0]/r;
    sinx=x[1]/r;

    cos2theta=2*cosx*cosx-1;
    sin2theta=2*cosx*sinx;

    // potential
    *phi = r*r*(mod[0] + mod[1]*cos2theta + mod[2]*sin2theta)/2;
  }
  
  //printf("  lens_expand *phi = %e\n",*phi);
  return mod[0];
}

void LensHalo::assignParams_stars(InputParams& params){

   	if(!params.get("main_stars_fraction",star_fstars)) error_message1("main_stars_fraction",params.filename());
   	if(star_fstars < 0 || star_fstars > 1){
   		ERROR_MESSAGE();
    	cout << "main_stars_fraction cannot be less than 0 or larger than 1 in file " << params.filename() <<endl;
    	exit(0);
   	}
   	if(!params.get("main_stars_mass",star_massscale)) error_message1("main_stars_mass",params.filename());

	if(!params.get("main_stars_mass_function",main_stars_imf_type)) error_message1("main_stars_mass_function",params.filename());
  if(main_stars_imf_type != One){
    if(!params.get("main_stars_min_mass",main_stars_min_mass)) error_message1("main_stars_min_mass",params.filename());
    if(main_stars_imf_type != Mono){
      if(!params.get("main_stars_max_mass",main_stars_max_mass)) error_message1("main_stars_max_mass",params.filename());
      if(!params.get("main_stars_bending_point",bend_mstar)) error_message1("main_stars_bending_point",params.filename());
      if(!params.get("main_stars_lo_mass_slope",lo_mass_slope)) error_message1("main_stars_lo_mass_slope",params.filename());
      if(!params.get("main_stars_hi_mass_slope",hi_mass_slope)) error_message1("main_stars_hi_mass_slope",params.filename());
    }
  }else{
    main_stars_min_mass = main_stars_max_mass = 1.0;    
  }

	return;
}

void LensHaloUniform::assignParams(InputParams& params){

	//if(perturb_Nmodes > 0){
	if(!params.get("main_zlens",zlens)) error_message1("main_zlens",params.filename());
	if(!params.get("kappa_uniform",kappa_uniform)) error_message1("kappa_uniform",params.filename());
	if(!params.get("gamma_uniform_1",gamma_uniform[0])) error_message1("gamma_uniform_1",params.filename());
	if(!params.get("gamma_uniform_2",gamma_uniform[1])) error_message1("gamma_uniform_2",params.filename());
	if(!params.get("zsource_reference",zsource_reference)) error_message1("zsource_reference",params.filename());
	if(!params.get("main_stars_N",stars_N)) error_message1("main_stars_N",params.filename());
    else if(stars_N){
    	assignParams_stars(params);
    }

	return;
}

/** \ingroup ImageFinding
 * \brief Prints the parameters of the analytic lens to stdout
 */
void LensHaloUniform::PrintLens(bool show_substruct,bool show_stars) const
{
	// uni lens parameters only
	cout << endl << "**Host lens model**" << endl;
	// redshifts
	cout << "zlens " << getZlens() << endl;
	cout << "kappa " << Sigma_uniform/SigmaCrit << endl;
	cout << "gamma " << gammaCrit_uniform[0]/SigmaCrit << " " << gammaCrit_uniform[1]/SigmaCrit << endl;

	if (stars_implanted) PrintStars(show_stars);
}

/* creates a single star halo in pos (x,y)
void LensHalo::implant_stars(
      PosType x, PosType y,long *seed, IMFtype type){

	if(Nregions <= 0) return;

	PosType ** centers = new PosType*[Nregions];
	for (int i = 0; i < Nregions; ++i){
    centers[i] = new PosType[2];
  
    centers[i][0] = x[i];
    centers[i][1] = y[i];
  }
	implant_stars(centers,Nregions,seed,type);
  
  for (int i = 0; i < Nregions; ++i) delete[] centers[i];
  delete[] centers;
}*/


