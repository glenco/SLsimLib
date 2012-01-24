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

class haloM{
public:
	PosType **pos;
	IndexType N;
	float *masses;
	float *sizes;
	float *redshifts;

	haloM(double,double);
	~haloM();
};

typedef struct haloM *haloHndl;

class multiLens : public Lens{
public:
	double *redshift;
	double *Ds, *Dl, *Dls; //sngular diameter distances
	double *Sigma_crit;
	double mass_scale;

	ForceTreeHndl *halo_tree;

	multiLens(char*);
	~multiLens();

	void buildHaloTree();
	void setInternalParams();
	void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off);
};

#endif /* MULTIPLANE_H_ */
