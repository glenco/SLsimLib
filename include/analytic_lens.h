/*
 * nsie.h
 *
 *  Created on: Dec 8, 2009
 *      Author: R.B. Metcalf
 */

#ifndef analens_declare
#define analens_declare

#include <lens.h>
#include <Tree.h>
#include <forceTree.h>
#include <quadTree.h>
#include <source.h>

/**
 * \brief An "analytic" model to represent a lens on a single plane.
 *
 * The lens consists of a "host" lens which is a non-singular isothermal ellipsoid (NSIE) plus axial distortion
 * modes, substructures and stars.
 *
 *<pre>
 * Input Paramaters:
 *
 *  **** NSIE parameters
 * 	sigma                  Velocity dispersion of host NSIE.
 *	core                   Core size
 *	axis_ratio             Axis ratio of mass
 *	pos_angle              Position angle
 *	z_lens                 Redshift of lens
 *
 *  **** Distortion parameters
 *	NDistortionModes       Number of distortion modes to be used.  If zero the other distortion parameters are not needed.
 *	beta_perturb
 *	kappa_peturb
 *	gamma_peturb
 *	monopole_peturb
 *	quadrapole_peturb
 *	hexopole_peturb
 *	octopole_peturb
 *
 *  **** Substructure parameters
 *	NdensitySubstruct      Number density of substructures.  They are distributed uniformly.  If zero the other substructure parameters are not needed.
 *	beta_sub               Logorithmic slope of the internal clump profile.  Used if sub_type == powerlaw
 *	alpha_sub              Logorithmic slope of the mass function.
 *	R_submax               Maximum radius of most massive substructure (see Metcalf & Amara 2012)
 *	mass_max               Maximum mass
 *	mass_min               Minimum mass
 *	sub_type               Mass profile of clumps - 0 or nfw,1 or powerlaw, 2 or pointmass
 *
 *  **** Stars parameters
 *	Nstars                 Total number of stars that will be used in the simulation.  If zero the other star parameters are not needed.
 *	fstars                 Fraction of surface denity in stars.
 *	stars_mass             Mass of stars.
 *
 * The stars are not initially present.  They must be implanted later.
 *</pre>
 *
 * TODO BEN finish this.
 */
class AnaLens : public Lens{
public:
	//AnaLens(std::string);
	AnaLens(InputParams& params);
	~AnaLens();

  /// private: Einstein radius of host
  double host_ro;
  double host_sigma;

  /// redshift of lens
  double zlens;
  // private derived quantities
  /// private: conversion factor between Mpc on the lens plane and arcseconds
  double MpcToAsec;
  /// private: critical surface density
  double Sigma_crit;
  /// private: the time delay scale in days/Mpc^2
  double to;

  /// Angular size distance to lens plane
  double Dl;

  // host elliptical
  double *host_x;    /// not used yet
  double host_core;
  double host_axis_ratio;
  double host_pos_angle;    /// position angle

  // perturbations to host
  int perturb_Nmodes;    /// this includes two for external shear
  double perturb_beta;
  double *perturb_rms;
  double *perturb_modes;  ///first two are shear

  // private derived quantities

  /// substructures
  bool substruct_implanted;
  double sub_sigmaScale;
  double sub_Ndensity;
  /// actual number of substructures
  int sub_N;
  double **sub_x;
  float *sub_Rcut;
  float *sub_mass;
  /// slope of mass profile
  double sub_beta;
  /// slope of mass function
  double sub_alpha;
  /// radius of largest mass substructures
  double sub_Rmax;
  double sub_Mmax;
  double sub_Mmin;
  double sub_theta_force;
  ForceTree *sub_tree;
  IndexType *sub_substructures;
  ClumpInternal sub_type;

  /// pointer to function for calculating the deflection caused by a subclump
  void (*sub_alpha_func)(double *alpha,double *x,double Rtrunc,double mass,double r_scale
			,double *center,double Sigma_crit);
  /// pointer to function for calculating the convergence caused by a subclump
  float (*sub_kappa_func)(double *x,double Rtrunc,double mass,double r_scale
  		,double *center,double Sigma_crit);
  /// pointer to function for calculating the shear caused by a subclump
  void (*sub_gamma_func)(float *gamma,double *x,double Rtrunc,double mass,double r_scale
    		,double *center,double Sigma_crit);
  /// pointer to function for calculating the surface potential caused by a subclump
  float (*sub_phi_func)(double *x,double Rtrunc,double mass,double r_scale
    		,double *center,double Sigma_crit);

  /// stars

  int stars_N;
  IndexType *stars;
  PosType **stars_xp;
  //ForceTree *star_tree;
  QuadTree *star_tree;
  double star_massscale;
  /// star masses relative to star_massscles
  float *star_masses;
  double star_fstars;
  double star_theta_force;
  bool stars_implanted;
  int star_Nregions;

 /// Number of regions to be subtracted to compensate for the mass in stars

  double *star_region;
  double *star_kappa;
  double **star_xdisk;

  double getZlens();
  void setZlens(double zlens);
  void setInternalParams(CosmoHndl,SourceHndl);
  void assignParams(InputParams& params);
  void error_message1(std::string name,std::string filename);
  void PrintAnaLens(bool show_substruct,bool show_stars);
  void rayshooterInternal(double *ray, double *alpha, float *gamma, float *kappa, bool kappa_off);

  // in randoimize_lens.c
  void RandomizeHost(long *seed,bool tables);
  void RandomizeSigma(long *seed,bool tables);
  void RandomlyDistortLens(long *seed,int Nmodes);
  void AlignedRandomlyDistortLens(long *seed,double theta,int n);
  void RandomizeSubstructure(double rangeInRei,long *seed);
  void RandomizeSubstructure2(double rangeInRei,long *seed);
  void RandomizeSubstructure3(double rangeInRei,long *seed);
  double FractionWithinRe(double rangeInRei);
  double averageSubMass();

  // in readlens_ana.c
  void reNormSubstructure(double kappa_sub);
  void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off);
  void substract_stars_disks(PosType *ray,PosType *alpha
                  ,float *kappa,float *gamma);
};

double RandomFromTable(double *table,unsigned long Ntable,long *seed);
void setStars(AnaLens *lens, bool implanted);
void implant_stars(AnaLens *lens,Point *images,unsigned long Nimages,long *seed);

void alphaNSIE(double *alpha,double *xt,double f,double bc,double theta);
float kappaNSIE(double *xt,double f,double bc,double theta);
void gammaNSIE(float *gam,double *xt,double f,double bc,double theta);
double invmagNSIE(double *x,double f,double bc,double theta
                     ,float *gam,float kap);
double rmaxNSIE(double sigma,double mass,double f,double rc );
double ellipticRadiusNSIE(double *x,double f,double pa);

void rotation(float *xout,float *xin,double theta);
void rotation(double *xout,double *xin,double theta);

//  in powerlow.c

void alphaPowLaw(double *alpha,double *x,double R,double mass,double beta,double *center,double Sigma_crit);
float kappaPowLaw(double *x,double R,double mass,double beta,double *center,double Sigma_crit);
void gammaPowLaw(float *gamma,double *x,double R,double mass,double beta,double *center,double Sigma_crit);
float phiPowLaw(double *x,double R,double mass,double beta,double *center,double Sigma_crit);

// in nfw_lens.c
void alphaNFW(double *alpha,double *x,double Rtrunc,double mass,double r_scale
                ,double *center,double Sigma_crit);
float kappaNFW(double *x,double Rtrunc,double mass,double r_scale
                ,double *center,double Sigma_crit);
void gammaNFW(float *gamma,double *x,double Rtrunc,double mass,double r_scale
                ,double *center,double Sigma_crit);

// in lens_expand.c

double lens_expand(double beta,double *mod,int Nmodes,double *x,double *alpha,float *gamma,float *phi);

// in FullRange/implant_stars.c

//void implant_stars(Point *images,unsigned long Nimages,long *seed);
//void substract_stars_disks(AnaLens *lens,PosType *ray,PosType *alpha
 //               ,PosType *kappa,PosType *gamma);

// in mark_points.c
void MarkPoints(TreeHndl s_tree,AnaLens *lens,bool sb_cut,short invert);
void _MarkPoints(TreeHndl s_tree,AnaLens *lens,bool *sbcut);
bool InSource(double *ray,AnaLens *lens,bool surfacebright);

#endif
