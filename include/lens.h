/*
 * lens.h
 *
 *      Author: mpetkova
 */

#ifndef LENS_H_
#define LENS_H_

#include <point.h>
#include <source.h>

/// An abstract base class to represent a gravitational lens.
class Lens {
public:
	int Nplanes;

	double zlens;
	/// output file, not always used.
	//std::string outputfile;
	/// marks if the lens has been setup.
	bool set;

	Lens();
	virtual ~Lens();

	int getNplanes(){return Nplanes;};

	virtual void resetNplanes(CosmoHndl cosmo, int Np){};

	virtual void setInternalParams(CosmoHndl,SourceHndl) = 0;

	virtual void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off){};
	virtual void rayshooterInternal(double *ray, double *alpha, KappaType *gamma, KappaType *kappa, bool kappa_off){};

	virtual void RandomizeHost(long *seed,bool tables){};
	virtual void RandomizeSigma(long *seed,bool tables){};
	virtual double getZlens() = 0;
	virtual void setZlens(CosmoHndl cosmo,double zlens,double zsource = 1000) = 0;
};

/**
 * \brief This is a lens that does no lensing.  It is useful for testing and for running refinement code on sources.
 */
class DummyLens: public Lens{
public:
	DummyLens(): Lens(){};
	~DummyLens(){};

	virtual void setInternalParams(CosmoHndl,SourceHndl){}
	void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off){
		for(unsigned long i=0;i<Npoints;++i){
			i_points[i].image->x[0] = i_points[i].x[0];
			i_points[i].image->x[1] = i_points[i].x[1];
			i_points[i].image->kappa = i_points[i].kappa = 0.0;
			i_points[i].image->gamma[0] = i_points[i].gamma[0] = 0.0;
			i_points[i].image->gamma[1] = i_points[i].gamma[1] = 0.0;
			i_points[i].image->gamma[2] = i_points[i].gamma[2] = 0.0;
		}
	}
	void rayshooterInternal(double *ray, double *alpha, KappaType *gamma, KappaType *kappa, bool kappa_off){
		alpha[0] = alpha[1] = 0.0;
		gamma[0] = gamma[1] = 0.0;
		*kappa = 0.0;
	}
	virtual double getZlens(){return 0.0;}
	virtual void setZlens(CosmoHndl, double, double = 1000){std::cout << "Why would you want to change the lens redshift in a DummyLens?" << std::endl; exit(1);}
};

typedef Lens *LensHndl;

#endif /* LENS_H_ */
