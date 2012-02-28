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

	  Lens();
	  ~Lens();

	  int getNplanes();

	  virtual void setInternalParams(CosmoHndl,double) = 0;
	  virtual void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off) = 0;
	  virtual void RandomizeHost(long *seed,bool tables){};
	  virtual void RandomizeSigma(long *seed,bool tables){};
};

typedef Lens *LensHndl;

#endif /* LENS_H_ */
