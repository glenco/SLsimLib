/*
 * multiplane.h
 *
 *  Created on: Jan 12, 2012
 *      Author: mpetkova
 */

#ifndef MULTIPLANE_H_
#define MULTIPLANE_H_

#include "quadTree.h"
#include "utilities_slsim.h"
#include "planes.h"

#include <map>

GLAMER_TEST_USES(LensTest)

/**
 * \brief A class to represents a lens with multiple planes.
 *
 *<pre>
 *	The rays are traced through multiple deflections.  On each plane there is a deflection
 *	solver.  An LensHaloAnaNSIE or LensHaloMOKA can be put on one of the planes.  The other planes can be
 *	populated with random field_halos drawn from a mass function or they can be retrieved from an
 *	external catalog.
 *
 *   Input Parameters (variable names):
 *
 *   main_halo_on -- 0: no major lens present; 1: there is a major lens present
 *   main_halo_type -- profile type for the main DM lens halo
 *   	0 or nolens, 1 or NFW, 2 or PseudoNFW, 3 or PowerLaw, 4 or NSIE, 5 or AnaLens, 6 or UniLens, 7 or MOKALens, 8 or DummyLens
 *   main_galaxy_halo_type -- profile typ for the main galaxy lens halo 0 or none, 1 or NSIE
 *   redshift_planes_file -- asci file with the redshifts of the lensing planes, if not set then created internaly
 *   flag_switch_field_off -- false: field halos are created, true: no field halos are created; default is false
 *
 *   if flag_switch_field_off == false, i.e. there are field halos then also the following are used:
 *   field_Nplanes -- number of field planes
 *   fieldofview -- field of view of the light cone, filled with field halos
 *   field_internal_profile -- profile type of the DM field lens halos
 *   	0 or nolens, 1 or NFW, 2 or PseudoNFW, 3 or PowerLaw, 4 or NSIE, 5 or AnaLens, 6 or UniLens, 7 or MOKALens, 8 or DummyLens, 9 or Hernquist
 *   field_prof_internal_slope -- slope of the surface density for power law or pseudo nfw profiles
 *   field_internal_profile_galaxy -- profile type of the galaxy field halos; if not set, no galaxies are used
 *   	0 or none, 1 or NSIE
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
 *   # Cosmology - Any cosmological parameters that are not set will have default values
 *
 *  Omega_matter -- Total mass (baryons + dark matter) in the units of the critical density, optional
 *  Omega_lambda -- Density in a cosmological constant, if not set it will be = 1 - Omega_matter
 *  Omega_baryon -- Density in baryons
 *  Omega_neutrino -- Density in neutrinos
 *  hubble -- Hubble parameter in units of 100 km/s/Mpc
 *  sigm_8 -- normalization of power spectrum
 *
 * </pre>
 */

class Lens
{
public:
	Lens(long *seed);
	Lens(InputParams& params, Source* source, long *my_seed,CosmoParamSet cosmoset = Planck1yr);
	~Lens();

	/// marks if the lens has been setup.
	bool set;

	/// the total number of lens planes
	int getNplanes(){return lensing_planes.size();}

	/// field of view in square degrees
	double getfov(){return fieldofview;};
	void setfov(double fov){fieldofview=fov;};

	/// reset te number of planes, but keep the field halos and main lens
	void resetFieldNplanes(std::size_t field_Nplanes);
	/// keep the main lens and the number of planes constant, but generate new field halos
	void resetFieldHalos();

	/// print the main parameters of the lens
	void printMultiLens();
	double getZlens(){
		if(flag_switch_main_halo_on)
			return main_halos[0]->getZlens();
		else{
			ERROR_MESSAGE();
			std::cout << "error, no main lens present" << std::endl;
			exit(1);
		}
	}

	/// remove all main halos
	void clearMainHalos();

	/// inserts a single main lens halo and adds it to the existing ones
	void insertMainHalo(Source* source, LensHalo* halo);
	/// inserts a sequence of main lens halos and adds them to the existing ones
	void insertMainHalos(Source* source, LensHalo** halos, std::size_t Nhalos);

	/// replaces existing main halos with a single main halo
	void replaceMainHalos(Source* source, LensHalo* halo);
	/// replaces existing main halos with a sequence of main halos
	void replaceMainHalos(Source* source, LensHalo** halos, std::size_t Nhalos);

	/// compute the dflection, convergence, and shear for each point on the grid
	void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off);
	/// compute the dflection, convergence, and shear for a single ray
	void rayshooterInternal(double *ray, double *alpha, KappaType *gamma, KappaType *kappa, bool kappa_off);

	// methods used for use with implanted sources

	short ResetSourcePlane(double z,bool nearest, unsigned long GalID=0, double *xx=NULL);

	/// Revert the source redshift to the value it was when the Lens was created.
	void RevertSourcePlane(){ toggle_source_plane = false;}
	//void ImplantSource(unsigned long index,CosmoHndl cosmo);
	double getSourceZ(){
		if(toggle_source_plane){
			return zs_implant;
		}else{
			return plane_redshifts.back();
		}
	}

	double getZmax(){return plane_redshifts.back();}

	/// print the cosmological parameters
	void PrintCosmology(){cosmo->PrintCosmology();}

private:
	GLAMER_TEST_FRIEND(LensTest)
	
	COSMOLOGY *cosmo;
	
	long *seed;
	
	/// field of view in square degrees
	double fieldofview;
	
	void readCosmology(InputParams& params);
	void assignParams(InputParams& params);
	
	/// turns source plane on and off
	bool toggle_source_plane;
	/// the distance from the source to the next plane
	double dDs_implant;
	double zs_implant,Ds_implant;
	/// This is the index of the plane at one larger distance than the new source distance
	int index_of_new_sourceplane;
	
	/// This is the source redshift that is read in from the parameter file and becomes the maximum redshift
	double zsource;
	
	void quicksort(LensHaloHndl *halo,double **pos,unsigned long N);
	
private: /* generation */
	/// create the lens planes
	void buildPlanes(InputParams& params, Source* source);
	
	/// sets the distances and redshifts of the field planes equidistant
	void setFieldDist();
	/// load the redshifts of the field planes from a file
	void setFieldDistFromFile();
	/// setup the field plane distances
	void setupFieldPlanes();
	/// create field halos as specified in the parameter file
	void createFieldHalos();
	/// read field halo data in from a file
	void readInputSimFile();
	/// build the field planes and sort halos onto them
	void createFieldPlanes();
	
	/// generate main halo from the parameter file
	void createMainHalos(InputParams& params, Source* source);
	/// generate main halo from the parameter file
	void createMainPlanes();
	/// add a main halo to an existing plane, or create a new plane
	void addMainHaloToPlane(LensHalo* halo);
	
	/// combine field and main planes
	void combinePlanes();
	
private: /* force calculation */
	/// if >= 1, deflection in the rayshooting is switched off
	bool flag_switch_deflection_off;
	
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
	
private: /* field */
	/// if true, the background is switched off and only the main lens is present
	bool flag_switch_field_off;
	
	/// vector of all field halos
	std::vector<LensHalo*> field_halos;
	/// number of field planes
	std::size_t field_Nplanes;
	/// vector of all field planes
	std::vector<LensPlane*> field_planes;
	/// vector of field plane redshifts
	std::vector<double> field_plane_redshifts;
	/// vector of field plane distances
	std::vector<double> field_Dl;
	
	/// Perpendicular position of halo TODO: (In proper distance?)
	double **halo_pos;
	
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
	// mass fraction in the host galaxy
	//double field_galaxy_mass_fraction;
	
	std::string redshift_planes_file;
	bool read_redshift_planes;
	
	std::string field_input_sim_file;
	bool sim_input_flag;
	//std::string input_gal_file;
	//bool gal_input_flag;
	bool read_sim_file;
	
	/// increases are for cosmological mean number density of halos calculation
	double field_buffer;
	
private: /* main */
	/// having a main halo in the paramfile
	bool flag_switch_main_halo_on;
	
	/// vector of all main halos
	std::vector<LensHalo*> main_halos;
	/// vector of own main halos that will be deleted
	std::vector<LensHalo*> main_halos_created;
	/// vector of all main planes
	std::vector<LensPlane*> main_planes;
	/// vector of main plane redshifts
	std::vector<double> main_plane_redshifts;
	/// vector of main plane distances
	std::vector<double> main_Dl;
	
	/// main lens type: 0 or nolens, 1 or NFW, 2 or PseudoNFW, 3 or PowerLaw, 4 or NSIE, 5 or AnaLens, 6 or UniLens, 7 or MOKALens, 8 or DummyLens, 9 or Hernquist
	LensHaloType main_halo_type;
	/// galaxy lens type: 0 or none, 1 or NSIE
	GalaxyLensHaloType main_galaxy_halo_type;
};

typedef Lens* LensHndl;

#endif /* MULTIPLANE_H_ */
