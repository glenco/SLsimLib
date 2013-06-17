/*
 * planes.cpp
 *
 *  Created on: May 21, 2013
 *      Author: mpetkova
 */

#include "planes.h"

#include <iterator>

LensPlaneTree::LensPlaneTree(PosType **xpt,LensHaloHndl *my_halos,IndexType Nhalos,double my_sigma_background) : LensPlane(){
	halo_tree = new TreeQuad(xpt,my_halos,Nhalos,my_sigma_background);
}

LensPlaneTree::~LensPlaneTree(){
	delete halo_tree;
}

void LensPlaneTree::force(double *alpha,KappaType *kappa,KappaType *gamma,double *xx,bool kappa_off){
	halo_tree->force2D_recur(xx,alpha,kappa,gamma,kappa_off);
}

void LensPlaneTree::add(LensHalo* halo)
{
	bool not_yet_implemented = false;
	assert(not_yet_implemented);
}

void LensPlaneTree::remove(LensHalo* halo)
{
	bool not_yet_implemented = false;
	assert(not_yet_implemented);
}

LensPlaneSingular::LensPlaneSingular(LensHalo** my_halos, std::size_t my_Nhalos)
: LensPlane(), halos(my_halos, my_halos + my_Nhalos)
{
}

LensPlaneSingular::~LensPlaneSingular()
{
}

void LensPlaneSingular::force(double *alpha,KappaType *kappa,KappaType *gamma,double *xx,bool kappa_off){
	double alpha_tmp[2];
	KappaType kappa_tmp, gamma_tmp[3];

	alpha[0] = alpha[1] = 0.0;
	*kappa = 0.0;
	gamma[0] = gamma[1] = gamma[2] = 0.0;

	for(std::size_t i = 0, n = halos.size(); i < n; ++i)
	{
		alpha_tmp[0] = alpha_tmp[1] = 0.0;
		kappa_tmp = 0.0;
		gamma_tmp[0] = gamma_tmp[1] = gamma_tmp[2] = 0.0;

		halos[i]->force_halo(alpha_tmp,&kappa_tmp,gamma_tmp,xx,kappa_off,false);

		alpha[0] -= alpha_tmp[0];
		alpha[1] -= alpha_tmp[1];
		*kappa += kappa_tmp;
		gamma[0] += gamma_tmp[0];
		gamma[1] += gamma_tmp[1];
		gamma[2] += gamma_tmp[2];
	}
}

void LensPlaneSingular::add(LensHalo* halo)
{
	halos.push_back(halo);
}

void LensPlaneSingular::remove(LensHalo* halo)
{
	std::vector<LensHalo*>::iterator it = std::find(halos.begin(), halos.end(), halo);
	if(it != halos.end())
		halos.erase(it);
}
