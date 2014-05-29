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
 * \brief A uniform surface density and shear with or without stars.
 *
 *
 *<pre>
 * Input Parameters:
 *
 *  ****  Uniform model parameters
 *	kappa_uniform
 *  gamma_uniform_1
 *  gamma_uniform_2
 *
 *
 *  **** Stars parameters
 *	main_stars_N                 Total number of stars that will be used in the simulation.  If zero the other star parameters are not needed.
 *	main_stars_fraction                 Fraction of surface denity in stars.
 *	main_stars_mass             Mass of stars.
 *
 *	main_stars_mass_function  Choose between: One, Mono, Salpeter, BrokenPowerLaw, SinglePowerLaw, Chabrier.
 *						   One - all stellar masses set to 1 Msol - default value in case one does not give the IMFtype option in function implant_stars()
 *						   Mono - all stellar masses set to the same value - given in main_stars_min_mass and main_stars_max_mass, Note main_stars_min_mass and main_stars_max_mass must be equal.
 *						   Salpeter - Salpeter power law IMF with slope -2.35. Note: main_stars_min_mass and main_stars_max_mass must be defined and different from each other.
 * 						   BrokenPowerLaw - IMF for which two slopes (the low mass slope 'main_stars_lo_mass_slope' and the high mass slope 'main_stars_hi_mass_slope') can be defined as well as a bending point (main_stars_bending_point).
 *                         SinglePowerLaw - Like Salpeter but with main_stars_lo_mass_slope=main_stars_hi_mass_slope as a free parameter.
 * 						   Chabrier - Chabrier IMF with exponential cut-off toward low mass, main_stars_bending_point is fixed to 1.0 Msol, Slope of powerlaw above 1 Msol is fixed to -2.3. Note that the three Chabrier Parameters are defined in implant_stars.cpp and (for now) fixed to [0.086,0.22,0.57].
 *	main_stars_min_mass          	   Lower mass limit of above IMFs (not used for IMFs One and Mono)
 *	main_stars_max_mass			   Upper mass limit of above IMFs (not used for IMFs One and Mono)
 *  main_stars_bending_point      	   Bending mass for BrokenPowerLaw
 *  main_stars_lo_mass_slope            	   index of the powerlaw (this is the index of the low mass end below main_stars_bending_point if BrokenPowerLaw is selected)
 *  main_stars_hi_mass_slope			   	   index of the powerlaw above main_stars_bending_point only used for BrokenPowerLaw. Must be equal to main_stars_lo_mass_slope for SingelPowerLaw.
 *
 *  The stars are not initially present.  They must be implanted later.
 *</pre>
 *
 */

class LensHaloUniform: public LensHalo{
public:
	LensHaloUniform(InputParams& params, const COSMOLOGY& cosmo, bool verbose = false);
	LensHaloUniform(InputParams& params,bool verbose =false);
	~LensHaloUniform();

	void assignParams(InputParams& params);
	void PrintLens(bool show_substruct,bool show_stars) const;
  /// Average surface density in Msun/Mpc^2
	float getSigma_uniform() const {return Sigma_uniform;}
  /// Shear times the critical surface density
	float* getGammaCrit_uniform(){return gammaCrit_uniform;}
  /// Unitless convergence
	float getKappa_uniform() const {return Sigma_uniform/SigmaCrit;}
  /// Shear times the critical surface density
	PosType getAveMag() const { return 1.0/( pow(1-Sigma_uniform/SigmaCrit,2) -
                                  (gammaCrit_uniform[0]*gammaCrit_uniform[0] + gammaCrit_uniform[1]*gammaCrit_uniform[1])/SigmaCrit/SigmaCrit );}
  /// Critical surface density
	PosType getSigmaCrit() const { return SigmaCrit;}

  //PosType getPerturb_beta(){return perturb_beta;}
  //int getPerturb_Nmodes(){return perturb_Nmodes;}    /// this includes two for external shear
  //PosType *perturb_modes;  ///first two are shear

    void force_halo(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType const *xcm,bool no_kappa,bool subtract_point=false);
    
    void setCosmology(const COSMOLOGY& cosmo);

protected:

  PosType SigmaCrit;
  PosType *perturb_modes;  ///first two are shear  
  PosType lens_expand(PosType *mod,PosType const *x,PosType *alpha,KappaType *gamma,KappaType *phi);

  /// Surface density of lens in Msun/ Mpc^2
  float Sigma_uniform;
  /// gamma*Sigma_crit
  float gammaCrit_uniform[3];
  /// input values for kappa and gamma
  float kappa_uniform, gamma_uniform[3];
   // perturbations to host.  These are protected so that in some derived classes they can or cann't be changed.
  int perturb_Nmodes;    /// this includes two for external shear
  PosType perturb_beta;
  PosType *perturb_rms;

  /// redshift for which uniform kappa and gamma are valid
  float zsource_reference;
  PosType Dl, Ds, Dls;
};


#endif /* UNIFORM_LENS_H_ */
