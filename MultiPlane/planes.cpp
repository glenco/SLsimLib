/*
 * planes.cpp
 *
 *  Created on: May 21, 2013
 *      Author: mpetkova
 */

#include "planes.h"

#include <iterator>

LensPlaneTree::LensPlaneTree(LensHaloHndl *my_halos,IndexType Nhalos
                             ,PosType my_sigma_background,PosType inv_screening_scale)
: LensPlane(), halos(my_halos, my_halos + Nhalos)
{
	if(inv_screening_scale != 0) halo_tree = new TreeQuad(my_halos,Nhalos,my_sigma_background,5,0.1,true,inv_screening_scale);
	else halo_tree = new TreeQuad(my_halos,Nhalos,my_sigma_background);
}

LensPlaneTree::~LensPlaneTree(){
	delete halo_tree;
}

void LensPlaneTree::force(PosType *alpha,KappaType *kappa,KappaType *gamma
                          ,KappaType *phi,PosType *xx){
	halo_tree->force2D_recur(xx,alpha,kappa,gamma,phi);
}


void LensPlaneTree::addHalo(LensHalo* halo)
{
	bool not_yet_implemented = false;
	assert(not_yet_implemented);
}

void LensPlaneTree::removeHalo(LensHalo* halo)
{
	bool not_yet_implemented = false;
	assert(not_yet_implemented);
}

std::vector<LensHalo*> LensPlaneTree::getHalos()
{
	return halos;
}

void LensPlaneTree::getNeighborHalos(PosType ray[],PosType rmax
                                     ,std::vector<LensHalo*> &neighbors) const
{
  halo_tree->neighbors(ray,rmax,neighbors);
}

std::vector<const LensHalo*> LensPlaneTree::getHalos() const
{
	return std::vector<const LensHalo*>(halos.begin(), halos.end());
}

LensPlaneSingular::LensPlaneSingular(LensHalo** my_halos, std::size_t my_Nhalos)
: LensPlane(), halos(my_halos, my_halos + my_Nhalos)
{
}

LensPlaneSingular::~LensPlaneSingular()
{
}


/** \brief returns the lensing quantities of a ray in physical coordinates
 *
 *  Warning : Be careful, the sign of alpha is changed after the call of force_halo !
 *
 */
void LensPlaneSingular::force(PosType *alpha
                              ,KappaType *kappa
                              ,KappaType *gamma
                              ,KappaType *phi       /// in mass
                              ,PosType *xx          // Position in PhysMpc
                              )
{

	PosType alpha_tmp[2],x_tmp[2];
	KappaType kappa_tmp, gamma_tmp[3];
    KappaType phi_tmp;
    
	alpha[0] = alpha[1] = 0.0;
    x_tmp[0] = x_tmp[1] = 0.0;
	*kappa = 0.0;
	gamma[0] = gamma[1] = gamma[2] = 0.0;
    *phi = 0.0;
  
    // Loop over the different halos present in a given lens plane.
	for(std::size_t i = 0, n = halos.size(); i < n; ++i)
	{
        alpha_tmp[0] = alpha_tmp[1] = 0.0;
        kappa_tmp = 0.0;
        gamma_tmp[0] = gamma_tmp[1] = gamma_tmp[2] = 0.0;
        phi_tmp = 0.0;
            
        // Getting the halo position (in PhysMpc) :
        halos[i]->getX(x_tmp);

        // Taking the shift into account :
        x_tmp[0] = xx[0] - x_tmp[0]; // in PhysMpc
        x_tmp[1] = xx[1] - x_tmp[1];

        halos[i]->force_halo(alpha_tmp,&kappa_tmp,gamma_tmp,&phi_tmp,x_tmp,false);

        // Adding the temporary values to the different quantities :
        alpha[0] -= alpha_tmp[0];
        alpha[1] -= alpha_tmp[1];
        *kappa += kappa_tmp;
        gamma[0] += gamma_tmp[0];
        gamma[1] += gamma_tmp[1];
        gamma[2] += gamma_tmp[2];
        *phi += phi_tmp;
	}
}

/// It is assumed that the position of halo is in physical Mpc
void LensPlaneSingular::addHalo(LensHalo* halo)
{
	halos.push_back(halo);
}

void LensPlaneSingular::removeHalo(LensHalo* halo)
{
	std::vector<LensHalo*>::iterator it = std::find(halos.begin(), halos.end(), halo);
	if(it != halos.end())
		halos.erase(it);
}

std::vector<LensHalo*> LensPlaneSingular::getHalos()
{
	return halos;
}

std::vector<const LensHalo*> LensPlaneSingular::getHalos() const
{
	return std::vector<const LensHalo*>(halos.begin(), halos.end());
}

void LensPlaneSingular::getNeighborHalos(PosType ray[],PosType rmax
                                         ,std::vector<LensHalo*> &neighbors) const{
  PosType x[2],r2 = rmax*rmax;
  for(int i=0;i<halos.size();i++){
     halos[i]->getX(x);
    if(r2 > ( (ray[0]-x[0])*(ray[0]-x[0]) + (ray[1]-x[1])*(ray[1]-x[1]) ))
      neighbors.push_back(halos[i]);
  }
};

