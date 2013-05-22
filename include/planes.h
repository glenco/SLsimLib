/*
 * planes.h
 *
 *  Created on: May 21, 2013
 *      Author: mpetkova
 */

#ifndef PLANES_H_
#define PLANES_H_

#include "quadTree.h"

class LensPlane{
public:
	LensPlane(){};
	~LensPlane(){};

	virtual void force(double *alpha,KappaType *kappa,KappaType *gamma,double *xx,bool kappa_off) = 0;
};

class TreeLensPlane : public LensPlane{
public:
	TreeLensPlane(PosType **xpt,LensHaloHndl *my_halos,IndexType Nhalos,double my_sigma_background);
	~TreeLensPlane();

	void force(double *alpha,KappaType *kappa,KappaType *gamma,double *xx,bool kappa_off);
protected:
	QuadTree *halo_tree;
};

class SingularLensPlane : public LensPlane{
public:
	SingularLensPlane(LensHaloHndl *my_halos, IndexType Nhalos);
	~SingularLensPlane();

	void force(double *alpha,KappaType *kappa,KappaType *gamma,double *xx,bool kappa_off);
protected:
	LensHaloHndl *halos;
	IndexType Nhalos;
};

#endif /* PLANES_H_ */
