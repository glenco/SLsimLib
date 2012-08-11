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
//TODO Margarita, This comment does not fully explain the classes use according to the repository comment.  How is it different than MultiLens, etc.
/**
 * \brief A test class. Produces one halo centered at {0,0}.
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
	void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off, double zsource=-1);

private:
	typedef enum {PowerLaw, NFW, PseudoNFW} IntProfType;
	typedef enum {null, ana_lens, moka_lens} InputLens;

	long *seed;

	void readParamfile(std::string);

	/// charge for the tree force solver (4*pi*G*mass_scale)
	double charge;
	/// an array of smart pointers to the halo models on each plane
	HaloData *halodata;
	/// an array of smart pointers to halo trees on each plane, uses the haloModel in the construction
	ForceTree *halo_tree;

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

#endif /* SINGLELENS_H_ */
