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
 *	Nstars                 Total number of stars that will be used in the simulation.  If zero the other star parameters are not needed.
 *	fstars                 Fraction of surface denity in stars.
 *	stars_mass             Mass of stars.
 *
 *	stellar_mass_function  Choose between: One, Mono, Salpeter, BrokenPowerLaw, SinglePowerLaw, Chabrier.
 *						   One - all stellar masses set to 1 Msol - default value in case one does not give the IMFtype option in function implant_stars()
 *						   Mono - all stellar masses set to the same value - given in min_mstar and max_mstar, Note min_mstar and max_mstar must be equal.
 *						   Salpeter - Salpeter power law IMF with slope -2.35. Note: min_mstar and max_mstar must be defined and different from each other.
 * 						   BrokenPowerLaw - IMF for which two slopes (the low mass slope 'slope_1' and the high mass slope 'slope_2') can be defined as well as a bending point (bending_point).
 *                         SinglePowerLaw - Like Salpeter but with slope_1=slope_2 as a free parameter.
 * 						   Chabrier - Chabrier IMF with exponential cut-off toward low mass, bending_point is fixed to 1.0 Msol, Slope of powerlaw above 1 Msol is fixed to -2.3. Note that the three Chabrier Parameters are defined in implant_stars.cpp and (for now) fixed to [0.086,0.22,0.57].
 *	min_mstar          	   Lower mass limit of above IMFs (not used for IMFs One and Mono)
 *	max_mstar			   Upper mass limit of above IMFs (not used for IMFs One and Mono)
 *  bending_point      	   Bending mass for BrokenPowerLaw
 *  slope_1            	   index of the powerlaw (this is the index of the low mass end below bending_point if BrokenPowerLaw is selected)
 *  slope_2			   	   index of the powerlaw above bending_point only used for BrokenPowerLaw. Must be equal to slope_1 for SingelPowerLaw.
 *
 *  The stars are not initially present.  They must be implanted later.
 *</pre>
 *
 */

class LensHaloUniform: public LensHalo{
public:
	LensHaloUniform(InputParams& params);
	~LensHaloUniform();

	void assignParams(InputParams& params);
	void PrintLens(bool show_substruct,bool show_stars);
	float getKappa_uniform(){return kappa_uniform;}
	float* getGamma_uniform(){return gamma_uniform;}
	double getAveMag(){ return 1.0/( pow(1-kappa_uniform,2) - gamma_uniform[0]*gamma_uniform[0] - gamma_uniform[1]*gamma_uniform[1]);}

	  double getPerturb_beta(){return perturb_beta;}
	  int getPerturb_Nmodes(){return perturb_Nmodes;}    /// this includes two for external shear
	  double *perturb_modes;  ///first two are shear

		/// overridden function to calculate the lensing properties
		void force_halo(double *alpha,KappaType *kappa,KappaType *gamma,double *xcm,bool no_kappa,bool subtract_point=false);
		void setInternalParams(CosmoHndl);

protected:

   float kappa_uniform;
   float gamma_uniform[3];
   // perturbations to host.  These are protected so that in some derived classes they can or cann't be changed.
 int perturb_Nmodes;    /// this includes two for external shear
 double perturb_beta;
 double *perturb_rms;
 float reference_z;
 double Dl, Ds, Dls;


};


#endif /* UNIFORM_LENS_H_ */
