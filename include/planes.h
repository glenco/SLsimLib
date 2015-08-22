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
	
  virtual void force(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType *xx) = 0;
	
	virtual void addHalo(LensHalo* halo) = 0;
	virtual void removeHalo(LensHalo* halo) = 0;
	
	virtual std::vector<LensHalo*> getHalos() = 0;
	virtual std::vector<const LensHalo*> getHalos() const = 0;
  virtual void getNeighborHalos(PosType ray[],PosType rmax,std::vector<LensHalo*> &neighbors) const{};
};

/// A LensPlane with a TreeQuad on it to calculate the deflection caused by field lenses
class LensPlaneTree : public LensPlane{
public:
	LensPlaneTree(LensHaloHndl *my_halos,IndexType Nhalos,PosType my_sigma_background,PosType my_inv_screening_scale = 0);
	~LensPlaneTree();

	void force(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType *xx);
  void closestHalos(PosType *xx,std::list<LensHalo *> neighbors){
    
  };

	void addHalo(LensHalo* halo);
	void removeHalo(LensHalo* halo);
	
	std::vector<LensHalo*> getHalos();
	std::vector<const LensHalo*> getHalos() const;
  /// Get the halos on this plane that are wthin rmax of ray[]
  void getNeighborHalos(PosType ray[],PosType rmax,std::vector<LensHalo*> &neighbors) const;
	
private:
	std::vector<LensHalo*> halos;
	
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

  void force(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType *xx);
    
	void addHalo(LensHalo* halo);
	void removeHalo(LensHalo* halo);
	
	std::vector<LensHalo*> getHalos();
	std::vector<const LensHalo*> getHalos() const;
  void getNeighborHalos(PosType ray[],PosType rmax,std::vector<LensHalo*> &neighbors) const;
  
private:
	std::vector<LensHalo*> halos;  
};

#endif /* PLANES_H_ */
