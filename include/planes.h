/*
 * planes.h
 *
 *  Created on: May 21, 2013
 *      Author: mpetkova
 */

#ifndef PLANES_H_
#define PLANES_H_

#include "quadTree.h"

/// Base class representing a plane in redshift onto which lenses are placed.
class LensPlane{
public:
	LensPlane() {}
	virtual ~LensPlane() {} 
	
	virtual void force(double *alpha,KappaType *kappa,KappaType *gamma,double *xx,bool kappa_off) = 0;
	
	virtual void add(LensHalo* halo) = 0;
	virtual void remove(LensHalo* halo) = 0;
};

/// A LensPlane with a TreeQuad on it to calculate the deflection caused by field lenses
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

/** \brief A LensPlane with a list of LensHalo's in it.  
 *
 * The deflection is calculated by direct summation which can be slow for large numbers of LensHalo's.  
 * Main lenses are put onto these planes when they are added to the Lens.
*/
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
