/*
 * MOKAlens.h
 *
 *  Created on: Jun 8, 2012
 *      Author: mpetkova
 */

#include <lens.h>
#include <cosmo.h>
#include <Tree.h>

#ifndef MOKALENS_H_
#define MOKALENS_H_


/// A class to represents the MOKA lens map
class MOKALens : public Lens{
public:

	MOKALens(std::string);
	~MOKALens();

	bool set;

	void readParamfile(std::string);
	void rayshooterInternal(double *ray, double *alpha, double *gamma, double *kappa, bool kappa_off);
	void setParams(float *a1,float *a2,float *g1,float *g2, float *k, double *center,double range,long Np);
	void setZlens(double zlens);
	double getZlens();
	void setInternalParams(CosmoHndl,double);

	/// values for the map
	float* alpha1;
	float* alpha2;
	float* gamma1;
	float* gamma2;
	float* kappa1;
	double* center;
	double range;
	long Npixels;
};

#endif /* MOKALENS_H_ */
