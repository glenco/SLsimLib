/*
 * lens.h
 *
 *      Author: mpetkova
 */

#ifndef LENS_H_
#define LENS_H_

#include <cosmo.h>
#include <point.h>
#include <source.h>

/// An abstract base class to represent a gravitational lens.
class Lens {
protected:
	int Nplanes;

public:
	double zlens;
	/// output file, not always used.
	std::string outputfile;
	/// marks if the lens has been setup.
	bool set;

	Lens();
	virtual ~Lens();

	int getNplanes();

	virtual void setInternalParams(CosmoHndl,SourceHndl) = 0;
	virtual void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off){};
	virtual void rayshooterInternal(double *ray, double *alpha, double *gamma, double *kappa, bool kappa_off){};
	virtual void RandomizeHost(long *seed,bool tables){};
	virtual void RandomizeSigma(long *seed,bool tables){};
	virtual double getZlens() = 0;
	virtual void setZlens(double zlens) = 0;

};

typedef Lens *LensHndl;

#endif /* LENS_H_ */
