/*
 * halo_model.cpp
 *
 *  Created on: Jan 14, 2012
 *      Author: mpetkova
 */

#include <slsimlib.h>

extern CosmoHndl cosmo;
extern multiLens *lens;

haloM::haloM(double maxr, double zmax){
	int i;
	long seed;
	double r, theta;

	seed = 2203;

	N = 100;

	masses = new float[N];
	sizes = new float[N];
	redshifts = new float[N];

	pos = PosTypeMatrix(0,N-1,0,2);

	for(i = 0; i < N; i++){
		r = maxr*ran2(&seed);
		theta=2*pi*ran2(&seed);
		pos[i][0] = r*cos(theta);
		pos[i][1] = r*sin(theta);
		pos[i][2] = 0.0;

		redshifts[i] = ran2(&seed) * zmax;

		sizes[i] = ran2(&seed) * 0.3;

		masses[i] = ran2(&seed);
	}
}

haloM::~haloM(){
	free_PosTypeMatrix(pos,0,N-1,0,2);
	delete[] redshifts;
	delete[] sizes;
	delete[] masses;
}

multiLens::multiLens(char filename[]) : Lens(filename){
	redshift = new double[Nplanes];
	Ds = new double[Nplanes];
	Dl = new double[Nplanes];
	Dls = new double[Nplanes];
	Sigma_crit = new double[Nplanes];
	halo_tree = new ForceTreeHndl[Nplanes];
}

multiLens::~multiLens(){
	int i;

	delete[] halo_tree;
	delete[] Sigma_crit;
	delete[] Dls;
	delete[] Dl;
	delete[] Ds;
	delete[] redshift;
}

void multiLens::buildHaloTree(){
	PosType **pos;
	IndexType N;
	float *masses;
	float *sizes;
	double dz;
	int i, j, n;

    haloM *halo = new haloM(0.4, zsource);

    dz = redshift[1] - redshift[0];

	for(j = 0; j < Nplanes; j++){

		for(i = 0, N = 0; i < halo->N; i++)
			if(halo->redshifts[i] >= redshift[j] && halo->redshifts[i] < (redshift[j]+dz))
				N++;

		masses = new float[N];
		sizes = new float[N];
		pos = PosTypeMatrix(0,N-1,0,2);

		for(i = 0, n = 0; i < halo->N; i++)
			if(halo->redshifts[i] >= redshift[j] && halo->redshifts[i] < (redshift[j]+dz)){
				masses[n] = halo->masses[i];
				sizes[n] = halo->sizes[i];
				pos[n][0] = halo->pos[i][0];
				pos[n][1] = halo->pos[i][1];
				pos[n][2] = halo->pos[i][2];
			}

		halo_tree[j] = new ForceTree(pos,N,masses,sizes,true,true,5,2,true,0.1);

		free_PosTypeMatrix(pos,0,N-1,0,2);
		delete[] sizes;
		delete[] masses;
	}

	delete halo;
}

void multiLens::setInternalParams(){
	int j;

	mass_scale = 1.0;

	for(j = 0; j < Nplanes; j++){
		Ds[j] = cosmo->angDist(0,zsource);
		Dl[j] = cosmo->angDist(0,redshift[j]);
		Dls[j] = cosmo->angDist(redshift[j],zsource);

		Sigma_crit[j]=cosmo->angDist(0,zsource)/cosmo->angDist(redshift[j],zsource)
					/cosmo->angDist(0,redshift[j])/4/pi/Grav;
	}
}
