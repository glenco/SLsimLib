/*
 * multiplane.h
 *
 *  Created on: Jan 12, 2012
 *      Author: mpetkova
 */

#ifndef MULTIPLANE_H_
#define MULTIPLANE_H_

#include <analytic_lens.h>
#include <MOKAlens.h>
#include <cosmo.h>
#include <halo.h>
#include <utilities.h>
#include <quadTree.h>
#include <forceTree.h>
#include <sourceAnaGalaxy.h>

const int Nmassbin=32;
const double MaxLogm=16.;

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

	HaloData(double fov,double min_mass,double mass_scale,double z1,double z2,int mass_func_type,double alpha,CosmoHndl cosmo,long *seed);
	HaloData(HaloStructure *halostrucs,double **positions,unsigned long Nhaloss);
	HaloData(CosmoHndl cosmo,double mass,double zlens);
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
	void printMultiLens();

	void setInternalParams(CosmoHndl,SourceHndl);
	void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off);
	/// a poiner to the analytical lens
	LensHndl input_lens;
	AnaLens *analens;
#ifdef WITH_MOKA
	MOKALens *mokalens;
#endif
	/// field of view in square degrees
	double fieldofview;

	// methods used for use with implanted sources
	void ImplantSource(CosmoHndl cosmo,double z,double theta[],OverGalaxy *ana_source);
	void ImplantSource(unsigned long index,CosmoHndl cosmo);
	double *getImplantedSourceX(){
		assert(gal_input_flag);
		return anasource->get_theta();
	}

	unsigned long getImplantedNsources(){
		assert(gal_input_flag);
		return anasource->getNumberOfGalaxies();
	}

private:

	typedef enum {PS, ST, PL} MassFuncType;
	typedef enum {PowerLaw, NFW, PseudoNFW} IntProfType;
	typedef enum {null, ana_lens, moka_lens} InputLens;

	void setRedshifts();
	void setCoorDist(CosmoHndl cosmo, double zsource);
	double QuickFindFromTable(double y);

	static double *coorDist_table;
	static long ob_count;
	void make_table(CosmoHndl cosmo);

	long *seed;
	/// if = 0 there is no input lens, if = 1 there is an analytic lens, if = 2 there is a MOKA lens
	InputLens flag_input_lens;

	void readParamfile(std::string);
	/// Redshifts of lens planes, 0...Nplanes.  Last one is the source redshift.
	double *plane_redshifts;
	/// angular diameter distances
	double *Dl;
	/// dDl[j] is the distance between plane j-1 and j
	double *dDl;
	/// charge for the tree force solver (4*pi*G*mass_scale)
	double charge;
	/// an array of smart pointers to the halo models on each plane
	std::auto_ptr<HaloData> *halodata;
	/// an array of smart pointers to halo trees on each plane, uses the haloModel in the construction
	std::auto_ptr<QuadTree> *halo_tree;

	/* the following parameters are read in from the parameter file */

	/// type of mass function PS (0), ST (1), and power law (2) default is ST
	MassFuncType mass_func_type;
	/// slope of the mass function is mass_func_type == 2
	double pw_alpha;
	/// mass scale
	double mass_scale;
	/// min mass for the halo model
	double min_mass;
	/// internal profile type, 0=powerlaw,1=nfw,2=pseudoNfw
	IntProfType internal_profile;
	/// power law internal profile slope, need to be <= 0
	double pw_beta;
	/// pseudo NFW internal profile slope, needs to be a whole number and > 0
	double pnfw_beta;

	/// read particle/halo data in from a file
	void readInputSimFile(CosmoHndl cosmo);

	std::string input_sim_file;
	bool sim_input_flag;
	std::string input_gal_file;
	bool gal_input_flag;

	/// pointer to first of all the halo internal structures
	HaloStructure *halos;
	/// number of halos on all the planes
	IndexType Nhalos;
	double *halo_zs;
	double **halo_pos;

	// Variables for implanted source
	std::auto_ptr<SourceAnaGalaxy> anasource;
	/// the distance from the source to the next plane
	double dDs_implant;
	double zs_implant,ys_implant[2],Ds_implant;
	/// This is the index of the plane one ahead of the source
	int flag_implanted_source;

	/// nfw tables
	//bool tables_set;

	/// if >= 1, deflection in the rayshooting is wtiched if
	int flag_switch_deflection_off;

	void quicksort(HaloStructure *halos,double **brr,double *arr,unsigned long N);
};

typedef  MultiLens* MultiLensHndl;

void swap(float *a,float *b);
void swap(PosType *a,PosType *b);
void swap(IndexType a,IndexType b);
void swap(IndexType *a,IndexType *b);
void quicksort(IndexType *particles,float *redshifts,PosType **pos,float *sizes,float *masses,IndexType N);

#endif /* MULTIPLANE_H_ */
