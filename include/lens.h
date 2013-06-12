/*
 * multiplane.h
 *
 *  Created on: Jan 12, 2012
 *      Author: mpetkova
 */

#ifndef MULTIPLANE_H_
#define MULTIPLANE_H_

#include "quadTree.h"
#include "tables.h"
#include "utilities_slsim.h"
#include "planes.h"

#include <map>

/**
 * \brief A class to represents a lens with multiple planes.
 *
 *<pre>
 *	The rays are traced through multiple deflections.  On each plane there is a deflection
 *	solver.  An LensHaloAnaNSIE or LensHaloMOKA can be put on one of the planes.  The other planes can be
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
 *              --------------------------------  i = j == (main_halo_on % Nplanes)
 *
 *
 *              --------------------------------  i = 0 first plane with mass on it at finite distance from observer
 *
 *   Input Parameters (variable names):
 *
 *   outputfile -- filename for simulation output, usually in the main()
 *   Nplanes -- number of lensing planes
 *   main_halo_on -- 0: no major lens present; 1: there is a major lens present
 *   main_halo_type -- profile type for the main DM lens halo
 *   	0 or nolens, 1 or NFW, 2 or PseudoNFW, 3 or PowerLaw, 4 or NSIE, 5 or AnaLens, 6 or UniLens, 7 or MOKALens, 8 or DummyLens
 *   main_galaxy_halo_type -- profile typ for the main galaxy lens halo 0 or none, 1 or NSIE
 *   redshift_planes_file -- asci file with the redshifts of the lensing planes, if not set then created internaly
 *   flag_switch_field_off -- false: field halos are created, true: no field halos are created; default is false
 *
 *   if flag_switch_field_off == false, i.e. there are field halos then also the following are used:
 *   fieldofview -- field of view of the light cone, filled with field halos
 *   field_internal_profile -- profile type of the DM field lens halos
 *   	0 or nolens, 1 or NFW, 2 or PseudoNFW, 3 or PowerLaw, 4 or NSIE, 5 or AnaLens, 6 or UniLens, 7 or MOKALens, 8 or DummyLens, 9 or Hernquist
 *   field_prof_internal_slope -- slope of the surface density for power law or pseudo nfw profiles
 *   field_internal_profile_galaxy -- profile type of the galaxy field halos; if not set, no galaxies are used
 *   	0 or none, 1 or NSIE
 *   field_galaxy_mass_fraction -- if int_prof_gal_type is set, then this is the galaxy mass fraction
 *
 *   field_input_sim_file -- filename of the Millennium simulation data to be read in and used to populate the light cone with field halos
 *
 *   if field_input_sim_file is  _not_ set, then the field halos are generated from a mass function and the following are used:
 *   field_mass_func_type -- type of the halo mass function
 *   	PS (0), ST (1), and power law (2)
 *   mass_func_PL_slope -- slope of the mass function in the power law case; default is -1/6
 *   field_min_mass -- minimum mass for the generated field halos
 *   field_buffer -- a constant physical size buffer, padding every lens plane to increase its surface
 *
 *   zsource -- source redshift
 *   flag_switch_deflection_off -- false: deflection is on, true: deflection is off; default is false
 *
 *
 * </pre>
 */

class Lens{
public:
	Lens(InputParams& params,CosmoHndl cosmo, long *seed);
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
	/// updates the lensing plane for the main halos
	void updateMainHaloLensPlane();
	/// generate main halos from the parameter file
	void createMainHalos(InputParams& params, CosmoHndl cosmo, SourceHndl source);
	/// generate field halos from a mass function
	void createFieldHalos(CosmoHndl cosmo,long *seed);
	/// read field halo data in from a file
	void readInputSimFile(CosmoHndl cosmo);

	/// print the main parameters of the lens
	void printMultiLens();
	double getZlens(){
		if(main_halo_on)
			return main_halos[0]->getZlens();
		else{
			ERROR_MESSAGE();
			std::cout << "error, no main lens present" << std::endl;
			exit(1);
		}
	}

	/// inserts a single main lens halo and adds it to the existing ones
	void insertSingleMainHalo(CosmoHndl cosmo, SourceHndl source,LensHalo *halo);
	/// inserts a single main lens halo and deletes all previously existing ones
	void insertNewSingleMainHalo(CosmoHndl cosmo, SourceHndl source,LensHalo *halo);
	/// inserts a sequence of main lens halos and adds them to the existing ones
	void insertMainHalos(CosmoHndl cosmo, SourceHndl source,LensHaloHndl *halo, IndexType nhalos);
	/// inserts a sequence of main lens halos and erases all previously existing ones
	void insertNewMainHalos(CosmoHndl cosmo, SourceHndl source,LensHaloHndl *halo, IndexType nhalos);

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

	int main_halo_on;
	/// the lensing planes
	std::vector<LensPlane *> lensing_planes;
	/// Dl[j = 0...] angular diameter distances, comoving
	std::vector<double> Dl;
	/// dDl[j] is the distance between plane j-1 and j plane, comoving
	std::vector<double> dDl;
	/// Redshifts of lens planes, 0...Nplanes.  Last one is the source redshift.
	std::vector<double> plane_redshifts;
	/// charge for the tree force solver (4*pi*G)
	double charge;

	std::string redshift_planes_file;
	bool read_redshift_planes;

	/// if >= 1, deflection in the rayshooting is switched off
	bool flag_switch_deflection_off;
	/// if >= 1, the background is switched off and only the main lens is present
	bool flag_switch_field_off;

	/* MAIN HALOS */
	/// main lens type: 0 or nolens, 1 or NFW, 2 or PseudoNFW, 3 or PowerLaw, 4 or NSIE, 5 or AnaLens, 6 or UniLens, 7 or MOKALens, 8 or DummyLens, 9 or Hernquist
	LensHaloType main_halo_type;
	/// galaxy lens type: 0 or none, 1 or NSIE
	GalaxyLensHaloType main_galaxy_halo_type;
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
	void setCoorDistFromFile(CosmoHndl cosmo);

	long *seed;

	void assignParams(InputParams& params);

	/* the following parameters are read in from the parameter file */

	/// type of mass function PS (0), ST (1), and power law (2) default is ST
	MassFuncType field_mass_func_type;
	/// slope of the mass function is field_mass_func_type == 2
	double mass_func_PL_slope;
	/// min mass for the halo model
	double field_min_mass;
	/// internal halo profile type; needs to be 0 or PowerLaw, 1 or NFW, 2 or PseudoNFW, 3 or NSIE, 4 or PointMass
	LensHaloType field_int_prof_type;
	/// power law or pseudo NFW internal profile slope
	double field_prof_internal_slope;

	/// if true, each field halo contains an NSIE galaxy inside it
	bool flag_field_gal_on;
	/// galaxy subhalo profile type; needs to be 0 or PowerLaw, 1 or NFW, 2 or PseudoNFW, 3 or NSIE, 4 or PointMass
	GalaxyLensHaloType field_int_prof_gal_type;
	/// mass fraction in the host galaxy
	double field_galaxy_mass_fraction;


	std::string field_input_sim_file;
	bool sim_input_flag;
	//std::string input_gal_file;
	//bool gal_input_flag;
	bool read_sim_file;



	/* FIELD HALOS */
	/// vector of all lens field_halos in the light cone
	std::vector<LensHaloHndl> field_halos;
	/// number of field_halos on all the planes
	IndexType Nhalos;
  /// Perpendicular position of halo TODO (In proper distance?)
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
