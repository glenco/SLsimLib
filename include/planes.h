/*
 * planes.h
 *
 *  Created on: May 21, 2013
 *      Author: mpetkova
 */

#ifndef PLANES_H_
#define PLANES_H_

#include "quadTree.h"
#include "quadTreeHalos.h"
#include "lens_halos.h"

/// Base class representing a plane in redshift onto which lenses are placed.

class LensPlane{
public:
	LensPlane(float redshift):z(redshift) {}
	virtual ~LensPlane() {} 
	
  virtual void force(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType *xx) = 0;
	
	virtual void addHalo(LensHalo* halo) = 0;
	virtual void removeHalo(LensHalo* halo) = 0;
	
	virtual std::vector<LensHalo*> getHalos() = 0;
	virtual std::vector<const LensHalo*> getHalos() const = 0;
  virtual void getNeighborHalos(PosType ray[],PosType rmax,std::vector<LensHalo*> &neighbors) const{};
  
  float z;
};

/// A LensPlane with a TreeQuad on it to calculate the deflection caused by field lenses
class LensPlaneTree : public LensPlane{
public:

	LensPlaneTree(float z,LensHalo **my_halos,IndexType Nhalos,PosType my_sigma_background,PosType my_inv_screening_scale = 0);

  LensPlaneTree(const LensPlaneTree &p);
  LensPlaneTree(LensPlaneTree &&p):LensPlane(p.z){
    std::swap(*this,p);
  }

	~LensPlaneTree();

  LensPlaneTree & operator=(const LensPlaneTree &p);
  LensPlaneTree & operator=(LensPlaneTree &&p);
  
	void force(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType *xx);
	void addHalo(LensHalo* halo);
	void removeHalo(LensHalo* halo);
	
	std::vector<LensHalo*> getHalos();
	std::vector<const LensHalo*> getHalos() const;
  /// Get the halos on this plane that are wthin rmax of ray[]
  void getNeighborHalos(PosType ray[],PosType rmax,std::vector<LensHalo*> &neighbors) const;
	
private:

	std::vector<LensHalo *> halos;
  TreeQuadHalos<LensHalo> * halo_tree;
};

/** \brief A LensPlane with a list of LensHalo's in it.  
 *
 * The deflection is calculated by direct summation which can be slow for large numbers of LensHalo's.  
 * Main lenses are put onto these planes when they are added to the Lens.
*/
class LensPlaneSingular : public LensPlane{
public:
	LensPlaneSingular(float z,LensHaloHndl *my_halos, IndexType Nhalos);
  LensPlaneSingular(const LensPlaneSingular &p):LensPlane(p.z){
    halos = p.halos;
  }
  LensPlaneSingular(LensPlaneSingular &&p):LensPlane(p.z){
    std::swap(halos,p.halos);
  }
  ~LensPlaneSingular();

  LensPlaneSingular & operator=(const LensPlaneSingular &p){
    if(&p != this){
      z = p.z;
      halos = p.halos;
    }
    return *this;
  }
  LensPlaneSingular & operator=(LensPlaneSingular &&p){
    if(&p != this){
      z = p.z;
      std::swap(halos,p.halos);
    }
    return *this;
  }

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
