/*
 * singlelens.h
 *
 *  Created on: Aug 8, 2012
 *      Author: mpetkova
 */

#ifndef SINGLELENS_H_
#define SINGLELENS_H_

#include <analytic_lens.h>
#include <MOKAlens.h>
#include <cosmo.h>
#include <halo.h>
#include <utilities.h>
#include <quadTree.h>
#include <forceTree.h>
#include <multiplane.h>
/**
 * \brief A test class. Produces one halo (NFW, PeudoNFW, or PowerLaw) centered at {0,0}.
 *
 * This class is or testing the convergence, deflection, and shear of the three halo profiles that
 * we work with. On a single plane a halo is placed in the center.
 *
 * A grid can be created in the executable and then use the function saveProfiles to save the convergence
 * and so on, as a function of radius.
 *
 * */
class SingleLens : public Lens{
public:

	SingleLens(std::string,long *);
	~SingleLens();

	void buildHaloTrees(CosmoHndl cosmo,double zsource);
	void RandomizeHost(long *seed,bool tables);
	void RandomizeSigma(long *seed,bool tables);
	double getZlens();
	void setZlens(double zlens);

	void printSingleLens();

	void setInternalParams(CosmoHndl,SourceHndl);
	void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off);
	void rayshooterInternal(double*, double*, float*, float*, bool){};

private:
	typedef enum {PowerLaw, NFW, PseudoNFW} IntProfType;
	typedef enum {null, ana_lens, moka_lens} InputLens;

	long *seed;

	void readParamfile(std::string);

	/// charge for the tree force solver (4*pi*G*mass_scale)
	double charge;
	/// an array of smart pointers to the halo models on each plane
	HaloData *halo_data;
	/// an array of smart pointers to halo trees on each plane, uses the haloModel in the construction
	QuadTree *halo_tree;

	/* the following parameters are read in from the parameter file */
	///mass of the galaxy or cluster
	double mass;
	/// mass scale
	double mass_scale;
	/// min mass for the halo model
	IntProfType internal_profile;
	/// power law internal profile slope, need to be <= 0
	double pw_beta;
	/// pseudo NFW internal profile slope, needs to be a whole number and > 0
	double pnfw_beta;

	/// pointer to first of all the halo internal structures
	HaloStructure *halos;
	/// number of halos on all the planes
	IndexType Nhalos;

	//bool tables_set;

	void quicksort(HaloStructure *halos,double **brr,double *arr,unsigned long N);
};


void saveProfiles(PointList *points, double boxlMpc,int nx, int ny);

#endif /* SINGLELENS_H_ */
