/*
 * multiplane.h
 *
 *  Created on: Jan 12, 2012
 *      Author: mpetkova
 */

#ifndef MULTIPLANE_H_
#define MULTIPLANE_H_

#include <analytic_lens.h>
#include <cosmo.h>
#include <halo.h>
#include <utilities.h>


const int Nmassbin=32;
const double MaxLogm=16.;

class haloM;  // forward declaration
typedef haloM *haloHndl;

/// A class to represents a lens with multiple planes.
class MultiLens : public Lens{
public:
	/// Redshifts of lens planes, 0...Nplanes.  Last one is the source redshift.
	double *redshift;
	double *Dl, *dDl; /// angular diameter distances
	double charge;
	/// mass scale
	double mass_scale;
	IndexType *NhalosinPlane;
	ForceTreeHndl *halo_tree;
	/// min mass for the halo model
	double min_mass;
	AnaLens *analens;
	int flag_analens;
	haloM *haloModel; ///halo model

	// to be read in from the parameter file
	double fieldofview; /// field of view in square degrees
	int mass_func_type; /// type of mass function PS (0) and ST (1) default is ST


	MultiLens(string);
	~MultiLens();

	void buildHaloTree(CosmoHndl cosmo,double zsource);
	double getZlens();
	void setZlens(double zlens);
	void setRedshift(double zsource);
	void printMultiLens();
	void readParamfile(string);

	void setInternalParams(CosmoHndl,double zsource);
	void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off);
};

/// Class to hold information about the halos' number, positions, masses, etc.
class haloM{
public:
	PosType **pos;
	IndexType N;

    float *masses,*sizes,*redshifts;
    IndexType *index;
	haloM(double zsource,CosmoHndl cosmo,MultiLens* lens);
	~haloM();
private:
	std:: vector<double> Logm,Nhaloes;
};

void swap(float *a,float *b);
void swap(PosType *a,PosType *b);
void swap(IndexType a,IndexType b);
void swap(IndexType *a,IndexType *b);
void quicksort(IndexType *particles,float *redshifts,PosType **pos,float *sizes,float *masses,IndexType N);

#endif /* MULTIPLANE_H_ */
