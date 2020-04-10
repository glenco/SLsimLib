/*
 * readfiles_uni.cpp
 *
 *  Created on: Jan 22, 2013
 *      Author: Leier
 */


#include "slsimlib.h"

using namespace std;

/*
 * \brief Reads in a parameter file and sets up an uniform lens.
 *
 * Sets many parameters within the lens model, source model and
 * force calculation.
 */
/*
LensHaloUniform::LensHaloUniform(InputParams& params, const COSMOLOGY& cosmo, bool verbose)
: LensHalo(params)
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
*/

//LensHaloUniform::LensHaloUniform(InputParams& params, bool verbose): LensHalo(){
//  
//  assignParams(params);
//  
//  if(std::numeric_limits<float>::has_infinity){
//    Rmax = std::numeric_limits<float>::infinity();
//  }else{
//    Rmax = std::numeric_limits<float>::max();
//  }
//  perturb_Nmodes=3;
//  perturb_modes = new PosType[3];
//  
//  if(verbose) PrintLens(false,false);
//}


LensHaloUniform::LensHaloUniform(double zlens,double Sigma,Point_2d &Shear,COSMOLOGY &cosmo): LensHalo(zlens,cosmo){
  
  if(std::numeric_limits<float>::has_infinity){
    Rmax = std::numeric_limits<float>::infinity();
  }else{
    Rmax = std::numeric_limits<float>::max();
  }
  
  perturb_modes[0] = Sigma;
  perturb_modes[1] = Shear[0];
  perturb_modes[2] = Shear[1];
  
  assert(Shear[0] == Shear[0] && Shear[1] == Shear[1]);
  assert(Sigma == Sigma);
  
}

LensHaloUniform::~LensHaloUniform(){
}

void LensHaloUniform::force_halo(
                                 PosType *alpha     /// mass/Mpc
                                 ,KappaType *kappa
                                 ,KappaType *gamma
                                 ,KappaType *phi
                                 ,PosType const *xcm
                                 ,bool subtract_point /// if true contribution from a point mass is subtracted
                                 ,PosType screening
                                 )
{
  *kappa += lens_expand(perturb_modes,xcm,alpha,gamma,phi);

  alpha[0] *= -1;
  alpha[1] *= -1;
  
  //*phi += phi_tmp ;
  
  assert(alpha[0] == alpha[0] && alpha[1] == alpha[1]);
  assert(gamma[0] == gamma[0] && gamma[1] == gamma[1]);
  assert(kappa == kappa);
}

PosType LensHaloUniform::lens_expand(PosType *mod
                                     ,PosType const *x
                                     ,PosType *alpha
                                     ,KappaType *gamma
                                     ,KappaType *phi
                                     )
{
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

  if(r == 0.0)
  {
    *phi = 0.0;
  }
  else
  {
    theta=atan2(x[1],x[0]);
    cosx=x[0]/r;
    sinx=x[1]/r;

    cos2theta=2*cosx*cosx-1;
    sin2theta=2*cosx*sinx;

    // potential
    *phi = r*r*(mod[0] + mod[1]*cos2theta + mod[2]*sin2theta)/2;
  }
  
  return mod[0];
}

//void LensHalo::assignParams_stars(InputParams& params){
//  
//  if(!params.get("main_stars_fraction",star_fstars)) error_message1("main_stars_fraction",params.filename());
//  if(star_fstars < 0 || star_fstars > 1){
//    ERROR_MESSAGE();
//    cout << "main_stars_fraction cannot be less than 0 or larger than 1 in file " << params.filename() <<endl;
//    exit(0);
//  }
//  if(!params.get("main_stars_mass",star_massscale)) error_message1("main_stars_mass",params.filename());
//  
//  if(!params.get("main_stars_mass_function",main_stars_imf_type)) error_message1("main_stars_mass_function",params.filename());
//  if(main_stars_imf_type != One){
//    if(!params.get("main_stars_min_mass",main_stars_min_mass)) error_message1("main_stars_min_mass",params.filename());
//    if(main_stars_imf_type != Mono){
//      if(!params.get("main_stars_max_mass",main_stars_max_mass)) error_message1("main_stars_max_mass",params.filename());
//      if(!params.get("main_stars_bending_point",bend_mstar)) error_message1("main_stars_bending_point",params.filename());
//      if(!params.get("main_stars_lo_mass_slope",lo_mass_slope)) error_message1("main_stars_lo_mass_slope",params.filename());
//      if(!params.get("main_stars_hi_mass_slope",hi_mass_slope)) error_message1("main_stars_hi_mass_slope",params.filename());
//    }
//  }else{
//    main_stars_min_mass = main_stars_max_mass = 1.0;
//  }
//  
//  return;
//}
