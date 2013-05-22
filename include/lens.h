/*
 * multiplane.h
 *
 *  Created on: Jan 12, 2012
 *      Author: mpetkova
 */

#ifndef MULTIPLANE_H_
#define MULTIPLANE_H_

#include "analytic_lens.h"
#include "uniform_lens.h"
#include "MOKAlens.h"
#include "quadTree.h"
#include "sourceAnaGalaxy.h"
#include "tables.h"
#include "utilities_slsim.h"
#include "planes.h"

#include <map>

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

class Lens{
public:
	//Lens(InputParams& params,long *seed);
	Lens(InputParams& params, CosmoHndl cosmo, SourceHndl source, long *my_seed);
	~Lens();

	/// marks if the lens has been setup.
	bool set;

	int getNplanes(){return Nplanes;};
	double getfov(){return fieldofview;};
	void setfov(double fov){fieldofview=fov;};

	/// output filename, to be usually used in the executable
	std::string outputfile;

	/// reset te number of planes, but keep the field halos and main lens
	void resetNplanes(CosmoHndl cosmo, int Np);
	/// keep the main lens and the number of planes constant, but generate new field halos
	void resetFieldHalos(CosmoHndl cosmo);

	/// build the lensing planes
	void buildLensPlanes(CosmoHndl cosmo);
	/// generate main halos from the parameter file
	void createMainHalos(InputParams& params, CosmoHndl cosmo, SourceHndl source);
	/// generate field halos from a mass function
	void createFieldHalos(CosmoHndl cosmo,long *seed);
	/// read field halo data in from a file
	void readInputSimFile(CosmoHndl cosmo);

	/// print the main parameters of the lens
	void printMultiLens();
	double getZlens(){
		if(flag_input_lens)
			return main_halos[0]->getZlens();
		else{
			ERROR_MESSAGE();
			std::cout << "error, no main lens present" << std::endl;
			exit(1);
		}
	}

	/// compute the dflection, convergence, and shear for each point on the grid
	void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off);
	/// compute the dflection, convergence, and shear for a single ray
	void rayshooterInternal(double *ray, double *alpha, KappaType *gamma, KappaType *kappa, bool kappa_off);

	// methods used for use with implanted sources

	short ResetSourcePlane(CosmoHndl cosmo,double z,bool nearest, unsigned long GalID=0, double *xx=NULL);

	/// Revert the source redshift to the value it was when the Lens was created.
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

	int flag_input_lens;
	/// the lensing planes
	std::vector<LensPlane *> lensing_planes;
	/// Dl[j = 0...] angular diameter distances, comoving
	std::vector<double> Dl;
	/// dDl[j] is the distance between plane j-1 and j plane, comoving
	std::vector<double> dDl;
	/// Redshifts of lens planes, 0...Nplanes.  Last one is the source redshift.
	std::vector<double> plane_redshifts;
	/// charge for the tree force solver (4*pi*G*mass_scale)
	double charge;

	/// if >= 1, deflection in the rayshooting is switched off
	bool flag_switch_deflection_off;
	/// if >= 1, the background is switched off and only the main lens is present
	bool flag_switch_field_off;

	/* MAIN HALOS */
	/// main lens type: 0 or nolens, 1 or NFW, 2 or PseudoNFW, 3 or PowerLaw, 4 or NSIE, 5 or AnaLens, 6 or UniLens, 7 or MOKALens, 8 or DummyLens
	LensHaloType main_halo_type;
	/// galaxy lens type: 0 or none, 1 or NSIE
	GalaxyLensHaloType galaxy_halo_type;
	/// main lensing halo in the simulation
	Utilities::MixedVector<LensHaloHndl> main_halos;
	/// number of main halo profiles (or main halos)
	IndexType NmainHalos;

private:
	/// number of lensing planes + 1 in the simulation, the last plant is the source plane
	int Nplanes;
	/// field of view in square degrees
	double fieldofview;
	
	/// tables with angular distances (comoving) and corresponding redshifts
	std::map<double,double> coorDist_table;
	const static long NTABLE = 1000;
	bool table_set;
	void make_table(CosmoHndl cosmo);

	/// sets the distances and redshifts of the lensing planes
	void setCoorDist(CosmoHndl cosmo);

	long *seed;

	void assignParams(InputParams& params);

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
	LensHaloType int_prof_type;
	/// power law or pseudo NFW internal profile slope
	double halo_slope;

	/// if true, each field halo contains an NSIE galaxy inside it
	bool flag_galaxy_subhalo;
	/// galaxy subhalo profile type; needs to be 0 or PowerLaw, 1 or NFW, 2 or PseudoNFW, 3 or NSIE, 4 or PointMass
	GalaxyLensHaloType int_prof_gal_type;
	/// mass fraction in the host galaxy
	double galaxy_mass_fraction;


	std::string input_sim_file;
	bool sim_input_flag;
	//std::string input_gal_file;
	//bool gal_input_flag;
	bool read_sim_file;



	/* FIELD HALOS */
	/// vector of all lens field_halos in the light cone
	std::vector<LensHaloHndl> field_halos;
	/// number of field_halos on all the planes
	IndexType Nhalos;
	double **halo_pos;

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

	/// increases are for cosmological mean number density of halos calculation
	double field_buffer;

	void quicksort(LensHaloHndl *halo,double **pos,unsigned long N);
};

typedef  Lens* LensHndl;


#endif /* MULTIPLANE_H_ */
