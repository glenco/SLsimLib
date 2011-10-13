/*
 * nsie.h
 *
 *  Created on: Dec 8, 2009
 *      Author: R.B. Metcalf
 */
#include <cosmo.h>

#ifndef analens_declare
#define analens_declare

typedef struct analytic_lens{
  char outputfile[40];
  double zlens;
  double zsource;

  // host elliptical
  double *xHost;    // not used yet
  double core;
  double axis_ratio;
  double theta;    // position angle
  double sigma;

  // perturbations to host
  long Nmodes;    // this includes two for external shear
  double *modes;  //first two are shear

  // substructures
  long NSubstruct;
  double **xSubstruct;
  double *RcutSubstruct;
  double *massSubstruct;
  double betaSubstruct;

  // private derived quantities
  double ro;
  double MpcToAsec;    // conversion factor between Mpc on the lens plane and arcseconds
  double Sigma_crit;   // critical surface density

} AnaLens;

#endif

void alphaNSIE(double *alpha,double *xt,double f,double bc,double theta);
double kappaNSIE(double *xt,double f,double bc,double theta);
void gammaNSIE(double gam[2],double *xt,double f,double bc,double theta);
double invmagNSIE(double *x,double f,double bc,double theta
		     ,double *gam,double kap);
double phiNSIE(double *xt,double f,double bc,double theta);
void rotation(double *xout,double *xin,double theta);
void ReadParams_AnaLens(char *filename,struct cosmology *cosmo,AnaLens *lens);

//  in powerlow.c

void alphaPowLaw(double *alpha,double *x,double R,double mass,double beta,double *center,double Sigma_crit);
double kappaPowLaw(double *x,double R,double mass,double beta,double *center,double Sigma_crit);
void gammaPowLaw(double *gamma,double *x,double R,double mass,double beta,double *center,double Sigma_crit);

// in lens_expand.c

double lens_expand(double beta,double *mod,int Nmodes,double *x,double *alpha,double *gamma,double *dt);
