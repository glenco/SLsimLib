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
	LensPlane() {}
	virtual ~LensPlane() {} 
	
	virtual void force(double *alpha,KappaType *kappa,KappaType *gamma,double *xx,bool kappa_off) = 0;
	
	virtual void add(LensHalo* halo) = 0;
	virtual void remove(LensHalo* halo) = 0;
};

class LensPlaneTree : public LensPlane{
public:
	LensPlaneTree(PosType **xpt,LensHaloHndl *my_halos,IndexType Nhalos,double my_sigma_background);
	~LensPlaneTree();

	void force(double *alpha,KappaType *kappa,KappaType *gamma,double *xx,bool kappa_off);
	
	void add(LensHalo* halo);
	void remove(LensHalo* halo);
	
private:
	TreeQuad* halo_tree;
};

class LensPlaneSingular : public LensPlane{
public:
	LensPlaneSingular(LensHaloHndl *my_halos, IndexType Nhalos);
	~LensPlaneSingular();

	void force(double *alpha,KappaType *kappa,KappaType *gamma,double *xx,bool kappa_off);
	
	void add(LensHalo* halo);
	void remove(LensHalo* halo);
	
private:
	std::vector<LensHalo*> halos;
};

#endif /* PLANES_H_ */
