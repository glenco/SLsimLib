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
 *  If the magnitude g and position angle, theta,  are specified then
 *
 * Shear[0] = g * cos( 2*theta )
 * Shear[1] = g * sin(  2*theta )
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

class LensHaloKappaDisk: public LensHalo{
public:
 /// Direct constructor without any stars
   LensHaloKappaDisk(double zlens  /// lens redshift
                     ,double Sigma_in  /// surface density of lens in Msun/Mpc^2
                     ,double Rmax   /// radius in Mpc
                     ,const Point_2d ang_center  /// angular center
                     ,COSMOLOGY &cosmo
                     )
  :LensHalo(zlens,cosmo),Sigma(Sigma_in)
  {
  
    if(Rmax <= 0){
      std::cerr << "LensHaloKappaDisk: Rmax cannot be <= 0" << std::endl;
      throw std::invalid_argument("Rmax<=0");
    }
    LensHalo::setMass(PI*Rmax*Rmax*Sigma);
    LensHalo::Rmax = Rmax;
    LensHalo::setRsize(Rmax);
    LensHalo::set_rsize(Rmax);
    LensHalo::setTheta(ang_center);
   }

  ~LensHaloKappaDisk(){};

  /// surface density in Msun/Mpc^2
  float getSigma() const {return Sigma;};
    
  LensHaloKappaDisk(const LensHaloKappaDisk &h):LensHalo(h){
    Sigma = h.Sigma;
  }
  
  LensHaloKappaDisk & operator=(const LensHaloKappaDisk &h){
    if(this == &h) return *this;
    LensHalo::operator=(h);
    Sigma = h.Sigma;
    return *this;
 }

  void force_halo(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType const *xcm,bool subtract_point=false,PosType screening = 1.0){

      float tmp;
      PosType tmp_mass;
      
      PosType r=sqrt(xcm[0]*xcm[0] + xcm[1]*xcm[1]);

      if(r < Rmax){
        alpha[0] -= Sigma*xcm[0];
        alpha[1] -= Sigma*xcm[1];
        *kappa = Sigma;
      }else{
        tmp_mass = mass/PI/r/r; //Sigma * Rmax*Rmax/r/r; 
        alpha[0] -= tmp_mass*xcm[0];
        alpha[1] -= tmp_mass*xcm[1];

        tmp = 2*tmp_mass/r/r;
        gamma[0] -= 0.5*(xcm[0]*xcm[0] - xcm[1]*xcm[1])*tmp;
        gamma[1] -= xcm[0]*xcm[1]*tmp;
      }
    
    return;
  }
protected:
  double Sigma; // surface density
};
#endif /* UNIFORM_LENS_H_ */
