/*
 * lens.h
 *
 *      Author: mpetkova
 */

#include <cosmo.h>
#include <forceTree.h>
#include <point.h>

#ifndef LENS_H_
#define LENS_H_

/// An abstract base class to represent a gravitational lens.
class Lens {
protected:
	int Nplanes;

public:
	  /// output file, not always used.
	  std::string outputfile;
	  /// marks if the lens has been setup.
	  bool set;

	  Lens();
	  ~Lens();

	  int getNplanes();

	virtual void setInternalParams(CosmoHndl,double) = 0;
	virtual void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off) = 0;
	virtual void RandomizeHost(long *seed,bool tables){};
	virtual void RandomizeSigma(long *seed,bool tables){};
	virtual double getZlens() = 0;
	virtual void setZlens(double zlens) = 0;

};

typedef Lens *LensHndl;

#endif /* LENS_H_ */
