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

class Lens {
protected:
	int Nplanes;

public:
	  /// output file, not always used.
	  string outputfile;
	  /// marks if the lens has been setup.
	  bool set;

	  /// redshift of lens
	  double zlens;
	  /// private: Einstein radius of host
	  double host_ro;
	  double host_sigma;

	  // private derived quantities
	  /// private: conversion factor between Mpc on the lens plane and arcseconds
	  double MpcToAsec;
	  /// private: critical surface density
	  double Sigma_crit;
	  /// private: the time delay scale in days/Mpc^2
	  double to;

	Lens();
	~Lens();

	int getNplanes();

	virtual void setInternalParams(CosmoHndl,double){};
	virtual void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off) = 0;
	virtual void RandomizeHost(long *seed,bool tables){};
};

typedef Lens *LensHndl;

#endif /* LENS_H_ */
