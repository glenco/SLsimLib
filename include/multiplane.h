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

class multiLens : public Lens{
public:
	double *redshift;
	double *Ds, *Dl, *Dls; //sngular diameter distances
	double *Sigma_crit;
	double mass_scale;

	haloHndl halo;
	ForceTreeHndl *halo_tree;

	multiLens(string);
	~multiLens();

	void setRedshift(double);
	void printMultiLens();
	void readParamfile(string);
	void setInternalParams(CosmoHndl,double);
	void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off);
};

void buildHaloTree(multiLens *);

void swap(float *a,float *b);
void swap(PosType *a,PosType *b);
void swap(IndexType a,IndexType b);
void swap(IndexType *a,IndexType *b);
void quicksort(IndexType *particles,float *redshifts,PosType **pos,float *sizes,float *masses,IndexType N);

#endif /* MULTIPLANE_H_ */
