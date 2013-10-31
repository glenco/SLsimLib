/*
 * base_analens.h
 *
 *  Created on: Jan 22, 2013
 *      Author: mpetkova
 */

#ifndef BASE_ANALENS_H_
#define BASE_ANALENS_H_

#include "forceTree.h"
#include "quadTree.h"
#include "source.h"
#include "InputParams.h"

/**
 * \brief An "analytic" model to represent a lens on a single plane.
 *
 * The lens consists of a "host" lens which is a non-singular isothermal ellipsoid (NSIE) plus axial distortion
 * modes, substructures and stars.
 *
 *<pre>
 * Input Parameters:
 *
 *  **** Distortion parameters
 *	main_NDistortionModes       Number of distortion modes to be used.  If zero the other distortion parameters are not needed.
 *	main_perturb_beta
 *	kappa_peturb
 *	gamma_peturb
 *	monopole_peturb
 *	quadrapole_peturb
 *	hexopole_peturb
 *	octopole_peturb
 *
 *  **** Substructure parameters
 *	main_sub_Ndensity      Number density of substructures.  They are distributed uniformly.  If zero the other substructure parameters are not needed.
 *	main_sub_beta               Logarithmic slope of the internal clump profile.  Used if main_sub_type == powerlaw
 *	main_sub_alpha              Logarithmic slope of the mass function.
 *	main_sub_Rmax               Maximum radius of most massive substructure (see Metcalf & Amara 2012)
 *	main_sub_mass_max               Maximum mass
 *	main_sub_mass_min               Minimum mass
 *	main_sub_type               Mass profile of clumps - 0 or nfw,1 or powerlaw, 2 or pointmass
 *
 *  **** Stars parameters
 *	main_stars_N                 Total number of stars that will be used in the simulation.  If zero the other star parameters are not needed.
 *	main_stars_fraction                 Fraction of surface denity in stars.
 *	main_stars_mass             Mass of stars.
 *
 * The stars are not initially present.  They must be implanted later.
 *</pre>
 *
 * TODO: BEN finish this documentation.
 */
class LensHaloBaseNSIE : public LensHalo{
public:
	LensHaloBaseNSIE(InputParams& params);
	virtual ~LensHaloBaseNSIE();

  /// critical surface density
  double getSigma_crit(){return Sigma_crit;}
  /// the time delay scale in days/Mpc^2
  double get_to(){return to;}
   /// Angular size distance to lens plane
  double get_Dl(){return Dl;}
  /// conversion factor from Mpc to Arcsec
  double get_MpcToAsec(){return MpcToAsec;}

  // private derived quantities

	/// get the velocity dispersion
	double get_sigma(){return sigma;};
	/// get the NSIE radius
	double get_Rsize(){return Rsize;};
	/// get the axis ratio
	double get_fratio(){return fratio;};
	/// get the position angle
	double get_pa(){return pa;};
	/// get the core radius
	double get_rcore(){return rcore;};

  /// substructures
  bool AreSubStructImaplated(){return substruct_implanted;}
  double sub_sigmaScale;
  double sub_Ndensity;
  /// actual number of substructures
  int sub_N;
  double **sub_x;
  /// slope of mass profile
  double sub_beta;
  /// slope of mass function
  double sub_alpha;
  /// radius of largest mass substructures
  double sub_Rmax;
  double sub_Mmax;
  double sub_Mmin;
  double sub_theta_force;
  LensHalo *subs;
  TreeQuad *sub_tree;
  IndexType *sub_substructures;
  ClumpInternal main_sub_type;

  void setZlens(double zlens);
//  void setInternalParams(CosmoHndl,SourceHndl);
//  void setInternalParams(CosmoHndl cosmo);
  void assignParams(InputParams& params);
  void PrintLens(bool show_substruct,bool show_stars);
  void error_message1(std::string name,std::string filename);

  virtual void force_halo(double *alpha,KappaType *kappa,KappaType *gamma,double *xcm,bool no_kappa,bool subtract_point=false);

  // in randoimize_lens.c
  double averageSubMass();

  // in readlens_ana.c
  void reNormSubstructure(double kappa_sub);
  //void toggleStars(bool implanted);

  double getHost_Dl(){return Dl;}

  double getEinstein_ro(){return Einstein_ro;}

  double getPerturb_beta(){return perturb_beta;}
  IMFtype getIMF_type(){return main_stars_imf_type;}
  int getPerturb_Nmodes(){return perturb_Nmodes;}    /// this includes two for external shear
  double *perturb_modes;  ///first two are shear
	
	std::size_t Nparams() const;
	double getParam(std::size_t p) const;
	double setParam(std::size_t p, double value);
	
	void printCSV(std::ostream& out, bool header = false) const;

protected:

  /// critical surface density
  double Sigma_crit;
   /// the time delay scale in days/Mpc^2
  double to;
   /// Angular size distance to lens plane
  double Dl;
  /// Einstein radius
  double Einstein_ro;

	/// velocity dispersion of NSIE
	float sigma;
	/// Actual edge of mass distribution in elliptical radius, Rmax is the range beyond which the halo is a point mass
	float Rsize;
	/// axis ratio of surface mass distribution
	float fratio;
	/// position angle on sky, radians
	float pa;
	/// core size of NSIE
	float rcore;


    // perturbations to host.  These are protected so that in some derived classes they can or cann't be changed.
  int perturb_Nmodes;    /// this includes two for external shear
  double perturb_beta;
  double *perturb_rms;

  bool substruct_implanted;

   // private derived quantities
   /// private: conversion factor between Mpc on the lens plane and arcseconds
   double MpcToAsec;

   /// redshift for which the perturbation modes are normalised
   float zsource_reference;
   double Ds, Dls;

};

namespace Utilities{
	double RandomFromTable(double *table,unsigned long Ntable,long *seed);
	void rotation(float *xout,float *xin,double theta);
	void rotation(double *xout,double *xin,double theta);
}

void alphaNSIE(double *alpha,double *xt,double f,double bc,double theta);
KappaType kappaNSIE(double *xt,double f,double bc,double theta);
void gammaNSIE(KappaType *gam,double *xt,double f,double bc,double theta);
KappaType invmagNSIE(double *x,double f,double bc,double theta
                     ,float *gam,float kap);
double rmaxNSIE(double sigma,double mass,double f,double rc );
double ellipticRadiusNSIE(double *x,double f,double pa);
void quadMomNSIE(float mass,float Rmax,float f,float rc,float theta,double *quad);


//  in powerlow.c

void alphaPowLaw(double *alpha,double *x,double R,double mass,double beta,double *center,double Sigma_crit);
KappaType kappaPowLaw(double *x,double R,double mass,double beta,double *center,double Sigma_crit);
void gammaPowLaw(KappaType *gamma,double *x,double R,double mass,double beta,double *center,double Sigma_crit);
KappaType phiPowLaw(double *x,double R,double mass,double beta,double *center,double Sigma_crit);

// in nfw_lens.c
void alphaNFW(double *alpha,double *x,double Rtrunc,double mass,double r_scale
                ,double *center,double Sigma_crit);
KappaType kappaNFW(double *x,double Rtrunc,double mass,double r_scale
                ,double *center,double Sigma_crit);
void gammaNFW(KappaType *gamma,double *x,double Rtrunc,double mass,double r_scale
                ,double *center,double Sigma_crit);

// in lens_expand.c

double lens_expand(double beta,double *mod,int Nmodes,double *x,double *alpha,KappaType *gamma,KappaType *phi);

// in FullRange/implant_stars.c

//void implant_stars(Point *images,unsigned long Nimages,long *seed);
//void substract_stars_disks(LensHaloBaseNSIE *lens,PosType *ray,PosType *alpha
 //               ,PosType *kappa,PosType *gamma);


#endif /* BASE_ANALENS_H_ */
