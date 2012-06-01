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
#include <valarray>

const int Nmassbin=32;
const double MaxLogm=16.;

/// structure to hold information about the halos' positions, masses, etc.
struct HaloStructure{
	/// internal halo parameters
    float mass,Rmax,rscale;
};

/// Class that holds all the information about the halos' positions and their internal parameters on one plane.
class HaloData{
public:
	/// halo positions
	PosType **pos;
	/// halo structure with internal halo parameters such as mass, size, etc.
	HaloStructure *halos;
	/// number of halos in the halo model on the plane
	IndexType Nhalos;
	/// mean mass density
	double kappa_background;

	HaloData(double fov,double min_mass,double mass_scale,double z1,double z2,int mass_func_type,CosmoHndl cosmo,long *seed);
	HaloData(HaloStructure *halostrucs,double **positions,unsigned long Nhaloss);
	~HaloData();


private:
	/// flag which is set to indicate which constructor is used and thus how the the destructor should work.
	bool allocation_flag;
};

typedef HaloData *HaloDataHndl;

/// A class to represents a lens with multiple planes.
class MultiLens : public Lens{
public:

	MultiLens(std::string,long *);
	~MultiLens();

	void buildHaloTrees(CosmoHndl cosmo,double zsource);
	void RandomizeHost(long *seed,bool tables);
	void RandomizeSigma(long *seed,bool tables);
	double getZlens();
	void setZlens(double zlens);
	void setRedshift(double zsource);
	void printMultiLens();

	void setInternalParams(CosmoHndl,double zsource);
	void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off, double zsource=-1);
	void rayshooterInternal(unsigned long Npixels, Point *i_points, bool kappa_off, float*,float*,float*,float*,float*, double* center,double range);

	/// a poiner to the analytical lens
	AnaLens *analens;
	/// field of view in square degrees
	double fieldofview;
	/// if = 1 there is an analytical lens, if = 0 there is no analytic lens
	int flag_analens;
	/// the name of the MOKA input file
	std::string MOKA_input_file;

private:

	void readParamfile(std::string);
	/// Redshifts of lens planes, 0...Nplanes.  Last one is the source redshift.
	double *plane_redshifts;
	/// angular diameter distances
	double *Dl, *dDl;
	/// charge for the tree force solver (4*pi*G*mass_scale)
	double charge;
	/// an array of smart pointers to the halo models on each plane
	//HaloDataHndl *halodata;
	std::auto_ptr<HaloData> *halodata;
	/// an array of smart pointers to halo trees on each plane, uses the haloModel in the construction
	//ForceTreeHndl *halo_tree;
	std::auto_ptr<ForceTree> *halo_tree;

	/* the following parameters are read in from the parameter file */

	/// type of mass function PS (0) and ST (1) default is ST
	int mass_func_type;
	/// mass scale
	double mass_scale;
	/// min mass for the halo model
	double min_mass;
	/// internal profile type, 0=Gauss,1=powerlaw,2=nfw,3=pseudoNfw
	int internal_profile;

	/// read particle/halo data in from a file
	void readInputSimFile(CosmoHndl cosmo);

	std::string input_sim_file;
	bool sim_input_flag;

	/// pointer to first of all the halo internal structures
	HaloStructure *halos;
	/// number of halos on all the planes
	IndexType Nhalos;
	double *halo_zs;
	double **halo_pos;

	long *seed;

	void quicksort(HaloStructure *halos,double **brr,double *arr,unsigned long N);
};


void swap(float *a,float *b);
void swap(PosType *a,PosType *b);
void swap(IndexType a,IndexType b);
void swap(IndexType *a,IndexType *b);
void quicksort(IndexType *particles,float *redshifts,PosType **pos,float *sizes,float *masses,IndexType N);

#endif /* MULTIPLANE_H_ */
