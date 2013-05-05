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
#include <quadTree.h>
#include <sourceAnaGalaxy.h>


/** \brief Class that holds all the information about the halos' positions and their internal parameters on one plane.
 *
 */
 class HaloData{
public:
	/// number of halos in the halo model on the plane
	IndexType Nhalos;
	/// halo positions
	PosType **pos;
	/// halo structure with internal halo parameters such as mass, size, etc.
	HaloStructure *halos;
	/// mean mass density
	double sigma_background;
	/// redshifts
	double *z;
	/// ids
	unsigned long *haloID;

	HaloData(HaloStructure *halostrucs,double sb,double **positions,double *z, unsigned long *haloID,unsigned long Nhaloss,double Dl);
	~HaloData();
};

/**
 * \brief A class to represents a lens with multiple planes.
 *
 *<pre>
 *	The rays are traced through multiple deflections.  On each plane there is a deflection
 *	solver.  An AnaLens or MOKALens can be put on one of the planes.  The other planes can be
 *	populated with random halos drawn from a mass function or they can be retrieved from an
 *	external catalog.
 *
 *
 *	Lens plane indexing scheme
 *
 *              --------------------------------  i = Nplanes-1 = source plane, No mass
 *
 *              --------------------------------  i = Nplanes-2 last plane with mass on it
 *
 *
 *              --------------------------------  i = j == (flag_input_lens % Nplanes)
 *
 *
 *              --------------------------------  i = 0 first plane with mass on it at finite distance from observer
 *
 *   Input Parameters:
 *
 *	Nplanes                     Number of lens planes
 *	flag_input_lens             Implanted lens - 0: no lens, 1: AnaLens, 2: MOKALens, The redshifts and internal parameters are the same as for constructing these lenses separately, see classes for each type of lens
 *	fov                         Field of view for the simulation region (not nessisarily the grided region)
 *	internal_profile            The internal profile type for the halos, 0 or PowerLaw,1 or NFW,2 or PseudoNFW, 3 or NSIE, 4 or NFW_NSIE .
 *  halo_to_galaxy_ratio        If NFW_NSIE is chosen this must be set to the ratio of the mass put in each component.
 *	z_source                    The "source" redshift, but actually the redshift of the last plane in the simulation.  The source can be at higher or lower redshift (see ResetSourcePlane)
 *	input_simulation_file       File that contains the input catalog of halos.  If it is missing a random set of halos will be generated.
 *	mass_func_type              The mass function used to generate random halos 0 through 2 or PS (Press & Schechter), ST (Sheth & Torman) or PowLaw (Power-law).  Not needed if input_simulation_file is provided.
 *	min_mass                    Minimum mass of halo when mass function is used (solar masses).  Not used when catalog is used.
 *	mass_scale                  The conversion between the mass units used and solar masses.  Usually 1.
 *	field_buffer                Field of view buffer in physical, rest frame Mpc.  Default is 0. Set to provide a buffer to the field of view so that halos that are centered outside the conical field of view but overlap it will be included.
 *	deflection_off              If true turns deflection off for testing purposes, default if false.
 *  background_off              If true turns deflection caused by background surface density off for testing purposes, default if false
 *
 * </pre>
 */

class MultiLens : public Lens{
public:

  void unusedHalos();

	MultiLens(InputParams& params,long *seed);
	~MultiLens();

	std::string outputfile;

	void resetNplanes(CosmoHndl cosmo, int Np);
	void resetHalos(CosmoHndl cosmo);
	
	void calc_error_test(unsigned long Npoints,Point *point,bool kappa_off);
	void calc_error_test_multi(unsigned long Npoints,Point *i_points,bool kappa_off,CosmoHndl cosmo);

	void buildHaloTrees(CosmoHndl cosmo);
	void buildHaloTrees_test(CosmoHndl cosmo);
	void createHaloData(CosmoHndl cosmo,long *seed);
	void createHaloData_buffered(CosmoHndl cosmo,long *seed);
	void createHaloData_test(CosmoHndl cosmo,long *seed);
	void RandomizeHost(long *seed,bool tables);
	void RandomizeSigma(long *seed,bool tables);
	double getZlens();
	void setZlens(CosmoHndl cosmo,double zlens,double zsource);
	void printMultiLens();

	void setInternalParams(CosmoHndl,SourceHndl);
	void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off);
	void rayshooterInternal(double *ray, double *alpha, KappaType *gamma, KappaType *kappa, bool kappa_off);
	void rayshooterInternal_halos(unsigned long Npoints, Point *i_points, bool kappa_off, double *Dl_halos, double *dDl_halos);
	
	LensHndl input_lens;
	AnaLens *analens;
	MOKALens *mokalens;
	/// field of view in square degrees
	double fieldofview;

	// methods used for use with implanted sources

	short ResetSourcePlane(CosmoHndl cosmo,double z,bool nearest, unsigned long GalID=0, double *xx=NULL);

	/// Revert the source redshift to the value it was when the MultiLens was created.
	void RevertSourcePlane(){ toggle_source_plane = false;}
	//void ImplantSource(unsigned long index,CosmoHndl cosmo);
	double getSourceZ(){
		if(toggle_source_plane){
			return zs_implant;
		}else{
			return plane_redshifts[Nplanes-1];
		}
	}
	double getZmax(){return plane_redshifts[Nplanes-1];}

	/// if = 0 there is no input lens, if = 1 there is an analytic lens, if = 2 there is a MOKA lens
	InputLens flag_input_lens;
	/// Dl[j = 0...]angular diameter distances
	double *Dl;
	/// Redshifts of lens planes, 0...Nplanes.  Last one is the source redshift.
	double *plane_redshifts;
	/// dDl[j] is the distance between plane j-1 and j plane
	double *dDl;
	/// charge for the tree force solver (4*pi*G*mass_scale)
	double charge;
	/// an array of smart pointers to halo trees on each plane, uses the haloModel in the construction
	std::auto_ptr<QuadTree> *halo_tree;
	/// if >= 1, deflection in the rayshooting is switched if
	bool flag_switch_deflection_off;
	/// if >= 1, the background is switched of and only the main lens is present
	bool flag_switch_background_off;

private:

	void setCoorDist(CosmoHndl cosmo);
	
	double *coorDist_table;
	double *redshift_table;
	unsigned long NTABLE;
	bool table_set;
	void make_table(CosmoHndl cosmo);

	long *seed;

	void assignParams(InputParams& params);
	/// an array of smart pointers to the halo models on each plane
	std::auto_ptr<HaloData> *halo_data;

	double r_print_halos;

	/* the following parameters are read in from the parameter file */

	/// type of mass function PS (0), ST (1), and power law (2) default is ST
	MassFuncType mass_func_type;
	/// slope of the mass function is mass_func_type == 2
	double pw_alpha;
	/// mass scale
	double mass_scale;
	/// min mass for the halo model
	double min_mass;
	/// internal profile type, 0=powerlaw,1=nfw,2=pseudoNfw, 3=NSIE
	IntProfType internal_profile;
	/// power law internal profile slope, need to be <= 0
	double pw_beta;
	/// pseudo NFW internal profile slope, needs to be a whole number and > 0
	double pnfw_beta;
	double galaxy_mass_fraction;


	/// read particle/halo data in from a file
	void readInputSimFile(CosmoHndl cosmo);

	std::string input_sim_file;
	bool sim_input_flag;
	//std::string input_gal_file;
	//bool gal_input_flag;
	bool read_sim_file;
    // TODO Margarita:  This parameter is undocumented!!!
	bool partial_cone;

	/// enables to two plane test
	bool flag_run_twop_test;
	/// enables the multi planes halos test
	bool flag_run_multip_test;

	/// pointer to first of all the halo internal structures
	HaloStructure *halos;
	/// number of halos on all the planes
	IndexType Nhalos;
	double *halo_zs;
	double **halo_pos_Mpc;
	double **halo_pos;
	unsigned long *halo_id;

	// Variables for implanted source
	//std::auto_ptr<MultiSourceAnaGalaxy> anasource;
	/// turns source plane on and off
	bool toggle_source_plane;
	/// the distance from the source to the next plane
	double dDs_implant;
	double zs_implant,Ds_implant;
	/// This is the index of the plane at one larger distance than the new source distance
	int index_of_new_sourceplane;

	/// This is the source redshift that is read in from the parameter file and becomes the maximum redshift
	double zsource;

	/// nfw tables
	//bool tables_set;
	double field_buffer;

	void quicksort(HaloStructure *halos,double **brr,double *arr,unsigned long *id,unsigned long N);
};

typedef  MultiLens* MultiLensHndl;

#endif /* MULTIPLANE_H_ */
