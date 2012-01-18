/*
 * halo_model.cpp
 *
 *  Created on: Jan 14, 2012
 *      Author: mpetkova
 */

#include <slsimlib.h>

extern CosmoHndl cosmo;
extern AnaLens *lens;

void set_HaloModel(haloHndl *halo, int Nplane, char *filename){
}

lensPlane::lensPlane(){
}

lensPlane::~lensPlane(){
	delete halo_tree;
}

void lensPlane::buildHaloTree(){
	PosType **halo_pos;
	IndexType halo_N;
	float *halo_masses;
	float *halo_sizes;
	double halo_theta_force = 0.1;

	halo_tree = new ForceTree(halo_pos,halo_N,halo_masses,halo_sizes,
			true,true,5,2,false,halo_theta_force);
}

void lensPlane::setInternalParams(){
	Ds = cosmo->angDist(0,lens->zsource);
	Dl = cosmo->angDist(0,redshift);
	Dls = cosmo->angDist(redshift,lens->zsource);

	Sigma_crit=cosmo->angDist(0,lens->zsource)/cosmo->angDist(redshift,lens->zsource)
			/cosmo->angDist(0,lens->zlens)/4/pi/Grav;
}
