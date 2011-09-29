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

typedef struct analytic_lens{
  char outputfile[100];
  Boolean set;   // marks if the lens has been setup

  double zlens;
  double zsource;

  // source parameters
  double source_tau;
  double source_nu;  // frequency
  double source_gauss_r2;  // internal scale parameter
  double source_r;   // total source size, ie no flux outside this radius
  double source_x[2];

  float source_nuo;
  float source_r_in;          //inner radius of BLR
  float source_r_out;         // outer radius of BLR
  float source_inclination;   //inclination of BLR in radians, face on is
  float source_opening_angle;
  float source_gamma;
  float source_BHmass;
  float source_fK;            // fraction of Keplerian velocity in random motions
  Boolean source_monocrome;   // set to true to integrate over frequency

  double (*source_sb_func)(double *y);  // surface brightness function
  SBModel source_sb_type;

  // host elliptical
  double *host_x;    // not used yet
  double host_core;
  double host_axis_ratio;
  double host_pos_angle;    // position angle
  double host_sigma;

  // perturbations to host
  long perturb_Nmodes;    // this includes two for external shear
  double perturb_beta;
  double *perturb_rms;
  double *perturb_modes;  //first two are shear

  // private derived quantities
  double host_ro;
  double MpcToAsec;    // conversion factor between Mpc on the lens plane and arcseconds
  double Sigma_crit;   // critical surface density
  double to;           // the time delay scale in days/Mpc^2

  // substructures
  Boolean substruct_implanted;
  double sub_sigmaScale;
  double sub_Ndensity;
  int sub_N;          // actual number of substructures
  double **sub_x;
  float *sub_Rcut;
  float *sub_mass;
  double sub_beta;   // slope of mass profile
  double sub_alpha;  // slope of mass function
  double sub_Rmax;   // radius of largest mass substructures
  double sub_Mmax;
  double sub_Mmin;
  double sub_theta_force;
  TreeNBHndl sub_tree;
  IndexType *sub_substructures;
  ClumpInternal sub_type;

  void (*sub_alpha_func)(double *alpha,double *x,double Rtrunc,double mass,double r_scale
			,double *center,double Sigma_crit);
  double (*sub_kappa_func)(double *x,double Rtrunc,double mass,double r_scale
  		,double *center,double Sigma_crit);
  void (*sub_gamma_func)(double *gamma,double *x,double Rtrunc,double mass,double r_scale
    		,double *center,double Sigma_crit);
  double (*sub_phi_func)(double *x,double Rtrunc,double mass,double r_scale
    		,double *center,double Sigma_crit);

  // stars
  Boolean stars_implanted;
  IndexType stars_N;
  IndexType *stars;
  PosType **stars_xp;
  TreeNBHndl star_tree;
  double star_massscale;
  float *star_masses;    // star masses relative to star_massscles
  double star_fstars;
  double star_theta_force;

 // regions to be subtracted to compensate for the mass in stars
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
void readparams_ana(char *filename,struct cosmology *cosmo,AnaLens *lens);

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
void PrintAnaLens(AnaLens *lens,Boolean show_substruct,Boolean show_stars);

// in randoimize_lens.c

void RandomizeHost(AnaLens *lens,double r_source_physical,long *seed,Boolean tables
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
void MarkPoints(TreeHndl s_tree,AnaLens *lens,Boolean sb_cut,short invert);
void _MarkPoints(TreeHndl s_tree,AnaLens *lens,Boolean *sbcut);
Boolean InSource(double *ray,AnaLens *lens,Boolean surfacebright);

