/*
 * nsie.h
 *
 *  Created on: Dec 8, 2009
 *      Author: R.B. Metcalf
 */
#ifndef simple_analens_declare
#define simple_analens_declare

/** \brief A simplified version of LensHaloAnaNSIE */
typedef struct simple_analytic_lens{
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

} SimpAnaLens;

/// deleted repeated function declarations

float phiNSIE(double *xt,double f,double bc,double theta);

#endif


