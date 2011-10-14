/*
 * nsie.h
 *
 *  Created on: Dec 8, 2009
 *      Author: R.B. Metcalf
 */
#include <cosmo.h>
#include <Tree.h>
#include <TreeNB.h>

#ifndef pi
#define pi  3.141593
#endif

#ifndef analens_declare
#define analens_declare

typedef enum {NFW,powerlaw,pointmass} ClumpInternal;
typedef enum {Uniform,Gaussian,BLR_Disk,BLR_Sph1,BLR_Sph2} SBModel;

/**\brief Analytic model for single plane lens and source.
 *
 */
typedef struct analytic_lens{
  /// output file, not always used.
  char outputfile[100];
  /// marks if the lens has been setup.
  bool set;

  /// redshift of lens
  double zlens;
  /// redshift of source
  double zsource;

  // source parameters
  /// lag time
  double source_tau;
  /// frequency
  double source_nu;
  /// internal scale parameter
  double source_gauss_r2;
  /// total source size, ie no flux outside this radius
  double source_r;
  /// center of source
  double source_x[2];

  float source_nuo;
  /// inner radius of BLR
  float source_r_in;
  /// outer radius of BLR
  float source_r_out;
  ///inclination of BLR in radians, face on is
  float source_inclination;
  float source_opening_angle;
  float source_gamma;
  float source_BHmass;
  /// fraction of Keplerian velocity in random motions
  float source_fK;
  /// set to true to integrate over frequency
  bool source_monocrome;

  /// pointer to surface brightness function
  double (*source_sb_func)(double *y);
  SBModel source_sb_type;

  // host elliptical
  double *host_x;    /// not used yet
  double host_core;
  double host_axis_ratio;
  double host_pos_angle;    /// position angle
  double host_sigma;

  // perturbations to host
  long perturb_Nmodes;    /// this includes two for external shear
  double perturb_beta;
  double *perturb_rms;
  double *perturb_modes;  ///first two are shear

  // private derived quantities
  /// private: Einstein radius of host
  double host_ro;
  /// private: conversion factor between Mpc on the lens plane and arcseconds
  double MpcToAsec;
  /// private: critical surface density
  double Sigma_crit;
  /// private: the time delay scale in days/Mpc^2
  double to;

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
  TreeNBHndl sub_tree;
  IndexType *sub_substructures;
  ClumpInternal sub_type;

  /// pointer to function for calculating the deflection caused by a subclump
  void (*sub_alpha_func)(double *alpha,double *x,double Rtrunc,double mass,double r_scale
			,double *center,double Sigma_crit);
  /// pointer to function for calculating the convergence caused by a subclump
  double (*sub_kappa_func)(double *x,double Rtrunc,double mass,double r_scale
  		,double *center,double Sigma_crit);
  /// pointer to function for calculating the shear caused by a subclump
  void (*sub_gamma_func)(double *gamma,double *x,double Rtrunc,double mass,double r_scale
    		,double *center,double Sigma_crit);
  /// pointer to function for calculating the surface potential caused by a subclump
  double (*sub_phi_func)(double *x,double Rtrunc,double mass,double r_scale
    		,double *center,double Sigma_crit);

  /// stars
  bool stars_implanted;
  IndexType stars_N;
  IndexType *stars;
  PosType **stars_xp;
  TreeNBHndl star_tree;
  double star_massscale;
  /// star masses relative to star_massscles
  float *star_masses;
  double star_fstars;
  double star_theta_force;

 /// Number of regions to be subtracted to compensate for the mass in stars
  int star_Nregions;
  double *star_region;
  double *star_kappa;
  double **star_xdisk;

} AnaLens;

#endif

void alphaNSIE(double *alpha,double *xt,double f,double bc,double theta);
double kappaNSIE(double *xt,double f,double bc,double theta);
void gammaNSIE(double gam[2],double *xt,double f,double bc,double theta);
double invmagNSIE(double *x,double f,double bc,double theta
		     ,double *gam,double kap);
void rotation(double *xout,double *xin,double theta);
void ReadParams_AnaLens(char *filename,struct cosmology *cosmo,AnaLens *lens);

//  in powerlow.c

void alphaPowLaw(double *alpha,double *x,double R,double mass,double beta,double *center,double Sigma_crit);
double kappaPowLaw(double *x,double R,double mass,double beta,double *center,double Sigma_crit);
void gammaPowLaw(double *gamma,double *x,double R,double mass,double beta,double *center,double Sigma_crit);
double phiPowLaw(double *x,double R,double mass,double beta,double *center,double Sigma_crit);

// in nfw_lens.c
void alphaNFW(double *alpha,double *x,double Rtrunc,double mass,double r_scale
		,double *center,double Sigma_crit);
double kappaNFW(double *x,double Rtrunc,double mass,double r_scale
		,double *center,double Sigma_crit);
void gammaNFW(double *gamma,double *x,double Rtrunc,double mass,double r_scale
		,double *center,double Sigma_crit);

// in lens_expand.c

double lens_expand(double beta,double *mod,int Nmodes,double *x,double *alpha,double *gamma,double *phi);
void free_AnaLens(AnaLens *lens);
void PrintAnaLens(AnaLens *lens,bool show_substruct,bool show_stars);

// in randoimize_lens.c

void RandomizeHost(AnaLens *lens,double r_source_physical,long *seed,bool tables
		,CosmoHndl cosmo);
void RandomlyDistortLens(AnaLens *lens,long *seed,int Nmodes);
void AlignedRandomlyDistortLens(AnaLens *lens,long *seed,double theta,int n);
double RandomFromTable(double *table,unsigned long Ntable,long *seed);
void RandomizeSubstructure(AnaLens *lens,double rangeInRei,long *seed);
void RandomizeSubstructure2(AnaLens *lens,double rangeInRei,long *seed);
void RandomizeSubstructure3(AnaLens *lens,double rangeInRei,long *seed);
double FractionWithinRe(AnaLens *lens,double rangeInRei);
double averageSubMass(AnaLens *lens);

// in FullRange/implant_stars.c

void implant_stars(AnaLens *lens,Point *images,unsigned long Nimages,long *seed);
void substract_stars_disks(AnaLens *lens,PosType *ray,PosType *alpha
		,PosType *kappa,PosType *gamma);

// in readlens_ana.c
void reNormSubstructure(AnaLens *lens,double kappa_sub);

// in internal_rayshooter.c
double uniform_SB(double *y);
double gaussian_SB(double *y);
double BLR_Disk_SB(double *y);
double BLR_Sph1_SB(double *y);
double BLR_Sph2_SB(double *y);

// in mark_points.c
void MarkPoints(TreeHndl s_tree,AnaLens *lens,bool sb_cut,short invert);
void _MarkPoints(TreeHndl s_tree,AnaLens *lens,bool *sbcut);
bool InSource(double *ray,AnaLens *lens,bool surfacebright);

