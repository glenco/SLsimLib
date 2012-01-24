/*
 * lens.h
 *
 *      Author: mpetkova
 */

#ifndef LENS_H_
#define LENS_H_

#include <cosmo.h>
#include <forceTree.h>
#include <point.h>

#ifndef lens_declare
#define lens_declare

class Source{
public:
	 /// names of clump and sb models
	  typedef enum {Uniform,Gaussian,BLR_Disk,BLR_Sph1,BLR_Sph2} SBModel;

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

	  /// redshift of source
	  double zsource;
};

class Lens : public Source {
protected:
	int Nplanes;

public:
	  /// output file, not always used.
	  char outputfile[100];
	  /// marks if the lens has been setup.
	  bool set;

	  /// redshift of lens
	  double zlens;

	  // private derived quantities
	  /// private: conversion factor between Mpc on the lens plane and arcseconds
	  double MpcToAsec;
	  /// private: critical surface density
	  double Sigma_crit;
	  /// private: the time delay scale in days/Mpc^2
	  double to;

	Lens(char*);
	~Lens();

	void readParamfile(char*);
	void setInternal(CosmoHndl cosmo);
	int getNplanes();

	virtual void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off){};
};

typedef Lens *LensHndl;

#endif

#endif /* LENS_H_ */
