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
typedef haloM *haloMHndl;

/// A class to represents a lens with multiple planes.
class MultiLens : public Lens{
public:
	/// Redshifts of lens planes, 0...Nplanes.  Last one is the source redshift.
	double *redshift;
	/// angular diameter distances
	double *Dl, *dDl;
	/// charge for the tree force solver (4*pi*G*mass_scale)
	double charge;
	/// an array of pointers to the halo models on each plane
	haloMHndl *haloModel;
	/// an array of pointers to halo trees on each plane, uses the haloModel in the construction
	ForceTreeHndl *halo_tree;
	/// a poiner to the analytical lens
	AnaLens *analens;

	/* the following parameters are read in from the parameter file */

	/// field of view in square degrees
	double fieldofview;
	/// type of mass function PS (0) and ST (1) default is ST
	int mass_func_type;
	/// mass scale
	double mass_scale;
	/// min mass for the halo model
	double min_mass;
	/// if = 1 there is an analytical lens, if = 0 there is no analytic lens
	int flag_analens;


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

/// structure to hold information about the halos' positions, masses, etc.
struct HaloStructure{
	/// internal halo parameters
    float mass,Rmax,rscale;
};

class haloM{
public:
	/// halo positions
	PosType **pos;
	/// halo structure with internal halo parameters such as mass, size, etc.
	HaloStructure *halos;
	/// number of halos in the halo model on the plane
	IndexType Nhalos;

	haloM(int jplane,double zsource,CosmoHndl cosmo,MultiLens* lens);
	~haloM();
private:
	/// variables for internal calculations
	std:: vector<double> Logm,Nhalosbin;
};

void swap(float *a,float *b);
void swap(PosType *a,PosType *b);
void swap(IndexType a,IndexType b);
void swap(IndexType *a,IndexType *b);
void quicksort(IndexType *particles,float *redshifts,PosType **pos,float *sizes,float *masses,IndexType N);

#endif /* MULTIPLANE_H_ */
