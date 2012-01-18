/*
 * multiplane.h
 *
 *  Created on: Jan 12, 2012
 *      Author: mpetkova
 */

#ifndef MULTIPLANE_H_
#define MULTIPLANE_H_

#include <forceTree.h>

struct halo_model{
	PosType **halo_pos;
	IndexType halo_N;
	float *halo_masses;
	float *halo_sizes;
	double halo_theta_force;
};

typedef struct halo_model *haloHndl;

void set_HaloModel(haloHndl);

class lensPlane{
public:
	double redshift;
	double Ds, Dl, Dls; //sngular diameter distances
	double Sigma_crit;
	double mass_scale;

	ForceTree *halo_tree;

	lensPlane();
	~lensPlane();

	void buildHaloTree();
	void setInternalParams();
};

void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off, int Nplanes);

#endif /* MULTIPLANE_H_ */
