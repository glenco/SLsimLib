/*
 * multiplane.h
 *
 *  Created on: Jan 12, 2012
 *      Author: mpetkova
 */

#ifndef MULTIPLANE_H_
#define MULTIPLANE_H_

#include <forceTree.h>
#include <lens.h>
#include <cosmo.h>
#include <halo.h>
#include <utilities.h>

const int Nmassbin=32;
const double MaxLogm=16.;

class MultiLens : public Lens{
public:
	/// Redshifts of lens planes, 0...Nplanes.  Last one is the source redshift.
	double *redshift;
	double *Dl, *dDl; /// angular diameter distances
	double charge;
	/// mass scaleå
	double mass_scale;
	unsigned long *NhalosinPlane;
	ForceTreeHndl *halo_tree;
	/// min mass for the halo model
	double min_mass;
	AnaLens *analens;
	int flag_analens;

	MultiLens(string);
	~MultiLens();

	void setRedshift(double);
	void printMultiLens();
	void readParamfile(string);

	void setInternalParams(CosmoHndl,double);
	void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off);
};

class haloM{
public:
	PosType **pos;
	IndexType N;

    float *masses,*sizes,*redshifts;
    IndexType *index;
	haloM(double,CosmoHndl, MultiLens*, double,int mfty=1);
	~haloM();
private:
	std:: vector<double> Logm,Nhaloes;
};

typedef haloM *haloHndl;

haloHndl buildHaloTree(MultiLens *, CosmoHndl,double,double);

void swap(float *a,float *b);
void swap(PosType *a,PosType *b);
void swap(IndexType a,IndexType b);
void swap(IndexType *a,IndexType *b);
void quicksort(IndexType *particles,float *redshifts,PosType **pos,float *sizes,float *masses,IndexType N);

#endif /* MULTIPLANE_H_ */
