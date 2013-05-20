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
#include "tables.h"

/**
 * \brief A class to represents a lens with multiple planes.
 *
 *<pre>
 *	The rays are traced through multiple deflections.  On each plane there is a deflection
 *	solver.  An AnaNSIELensHalo or MOKALensHalo can be put on one of the planes.  The other planes can be
 *	populated with random field_halos drawn from a mass function or they can be retrieved from an
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
 *	flag_input_lens             Implanted lens - 0: no lens, 1: AnaNSIELensHalo, 2: MOKALensHalo, The redshifts and internal parameters are the same as for constructing these lenses separately, see classes for each type of lens
 *	fov                         Field of view for the simulation region (not nessisarily the grided region)
 *	int_prof_type            The internal profile type for the field_halos, 0 or PowerLaw,1 or NFW,2 or PseudoNFW, 3 or NSIE, 4 or NFW_NSIE .
 *  halo_to_galaxy_ratio        If NFW_NSIE is chosen this must be set to the ratio of the mass put in each component.
 *	z_source                    The "source" redshift, but actually the redshift of the last plane in the simulation.  The source can be at higher or lower redshift (see ResetSourcePlane)
 *	input_simulation_file       File that contains the input catalog of field_halos.  If it is missing a random set of field_halos will be generated.
 *	mass_func_type              The mass function used to generate random field_halos 0 through 2 or PS (Press & Schechter), ST (Sheth & Torman) or PowLaw (Power-law).  Not needed if input_simulation_file is provided.
 *	min_mass                    Minimum mass of halo when mass function is used (solar masses).  Not used when catalog is used.
 *	mass_scale                  The conversion between the mass units used and solar masses.  Usually 1.
 *	field_buffer                Field of view buffer in physical, rest frame Mpc.  Default is 0. Set to provide a buffer to the field of view so that field_halos that are centered outside the conical field of view but overlap it will be included.
 *	deflection_off              If true turns deflection off for testing purposes, default if false.
 *  background_off              If true turns deflection caused by background surface density off for testing purposes, default if false
 *
 * </pre>
 */

class MultiLens{
public:
	//MultiLens(InputParams& params,long *seed);
	MultiLens(InputParams& params, CosmoHndl cosmo, SourceHndl source, long *my_seed);
	~MultiLens();

	/// marks if the lens has been setup.
	bool set;

	int getNplanes(){return Nplanes;};

	std::string outputfile;

	void resetNplanes(CosmoHndl cosmo, int Np);
	void resetFieldHalos(CosmoHndl cosmo);
	
	void calc_error_test(unsigned long Npoints,Point *point,bool kappa_off);
	void calc_error_test_multi(unsigned long Npoints,Point *i_points,bool kappa_off,CosmoHndl cosmo);

	void buildHaloTrees(CosmoHndl cosmo);
	void buildHaloTrees_test(CosmoHndl cosmo);
	void createHaloData(CosmoHndl cosmo,long *seed);
	void createHaloData_test(CosmoHndl cosmo,long *seed);
	void RandomizeHost(long *seed,bool tables);
	void RandomizeSigma(long *seed,bool tables);
	double getZlens();
	void setZlens(CosmoHndl cosmo,double zlens,double zsource);
	void printMultiLens();

	void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off);
	void rayshooterInternal(double *ray, double *alpha, KappaType *gamma, KappaType *kappa, bool kappa_off);
	void rayshooterInternal_halos(unsigned long Npoints, Point *i_points, bool kappa_off, double *Dl_halos, double *dDl_halos);

	/// needs to be 0 or nolens, 1 or NFW, 2 or PseudoNFW, 3 or PowerLaw, 4 or NSIE, 5 or AnaLens, 6 or UniLens, 7 or MOKALens
	InputLensType DM_halo_type;
	InputLensType galaxy_halo_type;
	/// main lensing halo in the simulation
	Utilities::MixedVector<LensHaloHndl> main_halos;

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

	/// Dl[j = 0...] angular diameter distances, comoving
	double *Dl;
	/// dDl[j] is the distance between plane j-1 and j plane, comoving
	double *dDl;
	/// Redshifts of lens planes, 0...Nplanes.  Last one is the source redshift.
	double *plane_redshifts;
	/// charge for the tree force solver (4*pi*G*mass_scale)
	double charge;

	/// an array of smart pointers to halo trees on each plane, uses the haloModel in the construction
	std::auto_ptr<QuadTree> *halo_tree;

	/// if >= 1, deflection in the rayshooting is switched off
	bool flag_switch_deflection_off;
	/// if >= 1, the background is switched off and only the main lens is present
	bool flag_switch_background_off;

private:

	int Nplanes;

	void setCoorDist(CosmoHndl cosmo);
	
	/// tables with angular distances (comoving) and corresponding redshifts
	double *coorDist_table;
	double *redshift_table;
	unsigned long NTABLE;
	bool table_set;
	void make_table(CosmoHndl cosmo);

	long *seed;

	void assignParams(InputParams& params);

	double r_print_halos;

	/* the following parameters are read in from the parameter file */

	/// type of mass function PS (0), ST (1), and power law (2) default is ST
	MassFuncType mass_func_type;
	/// slope of the mass function is mass_func_type == 2
	double mass_func_PL_slope;
	/// mass scale
	double mass_scale;
	/// min mass for the halo model
	double min_mass;
	/// internal halo profile type; needs to be 0 or PowerLaw, 1 or NFW, 2 or PseudoNFW, 3 or NSIE, 4 or PointMass
	IntProfType int_prof_type;
	/// internal galaxy profile type; needs to be 0 or PowerLaw, 1 or NFW, 2 or PseudoNFW, 3 or NSIE, 4 or PointMass
	IntProfType int_prof_gal_type;
	/// power law or pseudo NFW internal profile slope
	double halo_slope;
	/// mass fraction in the host galaxy
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
	/// enables the multi planes field_halos test
	bool flag_run_multip_test;

	/// vector of all lens field_halos in the light cone
	std::vector<LensHaloHndl> field_halos;
	/// number of field_halos on all the planes
	IndexType Nhalos;
	double *halo_zs;
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

	double field_buffer;

	bool second_halo;

	void quicksort(LensHaloHndl *halo,double **brr,double *arr,unsigned long *id,unsigned long N);
};

typedef  MultiLens* MultiLensHndl;


#endif /* MULTIPLANE_H_ */
