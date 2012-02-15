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

const int Nmassbin=256;
const double MaxLogm=16.;

class MultiLens : public Lens{
public:
	/// Redshifts of lens planes, 0...Nplanes.  Last one is the source redshift.
	double *redshift;
	double *Dl, *dDl; /// angular diameter distances
	double charge;
	double mass_scale;
	unsigned long *NhalosinPlane;
	//haloHndl halo;
	ForceTreeHndl *halo_tree;
	double mass_resolution;

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

void buildHaloTree(MultiLens *, CosmoHndl,double,double);

void swap(float *a,float *b);
void swap(PosType *a,PosType *b);
void swap(IndexType a,IndexType b);
void swap(IndexType *a,IndexType *b);
void quicksort(IndexType *particles,float *redshifts,PosType **pos,float *sizes,float *masses,IndexType N);

#endif /* MULTIPLANE_H_ */
