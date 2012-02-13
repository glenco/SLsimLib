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

class haloM{
public:
	PosType **pos;
	IndexType N;
	IndexType *index;
	float *masses;
	float *sizes;
	float *redshifts;

	haloM(double,double);
	~haloM();
};

typedef haloM *haloHndl;

class MultiLens : public Lens{
public:
	/// Redshifts of lens planes, 0...Nplanes.  Last one is the source redshift.
	double *redshift;
	double *Dl, *dDl; /// angular diameter distances
	double charge;
	double mass_scale;
	int Nplanes;

	haloHndl halo;
	ForceTreeHndl *halo_tree;

	MultiLens(string);
	~MultiLens();

	void setRedshift(double);
	void printMultiLens();
	void readParamfile(string);

	void setInternalParams(CosmoHndl,double);
	void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off);
};

void buildHaloTrees(MultiLens *);

void swap(float *a,float *b);
void swap(PosType *a,PosType *b);
void swap(IndexType a,IndexType b);
void swap(IndexType *a,IndexType *b);
void quicksort(IndexType *particles,float *redshifts,PosType **pos,float *sizes,float *masses,IndexType N);

#endif /* MULTIPLANE_H_ */
