/*
 * uniform_lens.h
 *
 *  Created on: Jan 21, 2013
 *      Author: D. Leier
 */

#ifndef UNIFORM_LENS_H_
#define UNIFORM_LENS_H_

#include "base_analens.h"

/**
 * \brief A uniform surface density and shear lens.
 *
 * Note the units of the input shear and surface density are not Sigma_crit's.
 * This lens will extend to infinity in the angular directions.
 *
 *
 */

class LensHaloUniform: public LensHalo{
public:
 /// Direct constructor without any stars
   LensHaloUniform(double zlens  /// lens redshift
                  ,double Sigma  /// surface density of lens in Msun/Mpc^2
                  ,Point_2d &Shear  /// shear in M_sun/Mpc^2 ie gamma * Sigma_crit
                  ,COSMOLOGY &cosmo
                  );

  ~LensHaloUniform();

  /// surface density in Msun/Mpc^2
  float getSigma() const {return perturb_modes[0];};
  /// shear in Msun/Mpc^2
  Point_2d getShear() const{
    Point_2d p(perturb_modes[1],perturb_modes[2]);
    return p;
  }
  
  /// reset the shear keeping the redshift of the plane fixed
  void resetGamma(const Point_2d &Shear  /// Shear in units of critical density Msun/Mpc^2, ie gamma * Sigma_crit
                  ){
    perturb_modes[1] = Shear[0];
    perturb_modes[2] = Shear[1];
  }
  
  /// reset the convergence keeping the redshift of the plane fixed
  void resetKappa(double Sigma  /// Surface density in Msun/Mpc^2, ie kappa * Sigma_crit
                  ){
    perturb_modes[0] = Sigma;;
  }
  
  void force_halo(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType const *xcm,bool subtract_point=false,PosType screening = 1.0);
  
  LensHaloUniform(const LensHaloUniform &h):LensHalo(h){
    perturb_modes[0] = h.perturb_modes[0];
    perturb_modes[1] = h.perturb_modes[1];
    perturb_modes[2] = h.perturb_modes[2];
  }
  
  LensHaloUniform & operator=(const LensHaloUniform &h){
    if(this == &h) return *this;
    LensHalo::operator=(h);
    perturb_modes[0] = h.perturb_modes[0];
    perturb_modes[1] = h.perturb_modes[1];
    perturb_modes[2] = h.perturb_modes[2];
    return *this;
 }

protected:

  PosType perturb_modes[3];  ///first two are shear
  PosType lens_expand(PosType *mod,PosType const *x,PosType *alpha,KappaType *gamma,KappaType *phi);
};


#endif /* UNIFORM_LENS_H_ */
