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
#include "geometry.h"

#include <map>

GLAMER_TEST_USES(LensTest)

/**
 * \brief A class to represents a lens with multiple planes.
 *
 *<pre>
 *	The rays are traced through multiple deflections.  On each plane there is a deflection
 *	solver.  An LensHaloAnaNSIE or LensHaloMassMap can be put on one of the planes.  The other planes can be
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
 *   if field_off == false, i.e. there are field halos then also the following are used:
 *   field_Nplanes -- number of field planes
 *   field_fov -- field of view of the light cone, filled with field halos
 *   field_internal_profile -- profile type of the DM field lens halos
 *   	0 or nolens, 1 or NFW, 2 or PseudoNFW, 3 or PowerLaw, 4 or NSIE, 5 or AnaLens, 6 or UniLens, 7 or MOKALens, 8 or DummyLens, 9 or Hernquist, 10 or Jaffe
 *   field_prof_internal_slope_pl -- slope of the surface density for power law
 *   field_prof_internal_slope_pnfw -- slope of the surface density for pseudo nfw profiles
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
 *   flag_switch_deflection_off -- false: deflection is on, but kappa and gamma may be calculated, true: deflection is off; default is false
 *   flag_switch_lensing_off -- false: lensing is on, true: lensing is off (no alpha, kappa or gamma); default is false
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
	Lens(long* seed, PosType z_source,CosmoParamSet cosmoset = WMAP5yr, bool verbose = false);
	Lens(InputParams& params, long* my_seed, CosmoParamSet cosmoset = WMAP5yr, bool verbose = false);
  //Lens(Lens &lens);
  
	~Lens();

	/// marks if the lens has been setup.
	bool set;

	/// the total number of lens planes
	int getNplanes(){return lensing_planes.size();}

	/// field of view in square degrees
	PosType getfov(){return fieldofview;};
	void setfov(PosType fov){fieldofview=fov;};

	/// reset the number of planes, but keep the field halos and main lens
	void resetFieldNplanes(std::size_t field_Nplanes, bool verbose = false);

	/** keep the main lens and the number of planes constant, but generate new field halos.
    This function will also erase the substructure halos so they need to be regenerated using 
   resetSubstructure() if they are desired .
   */
	void resetFieldHalos(bool verbose = false);

	/// print the main parameters of the lens
	void printMultiLens();

  
  /// Redshift of first main lens plane
	PosType getZlens(){
		if(flag_switch_main_halo_on)
			return main_halos[0]->getZlens();
		else{
			ERROR_MESSAGE();
			std::cout << "error, no main lens present" << std::endl;
			exit(1);
		}
	}
  /// Angular size distance (Mpc) to first main lens plane
	PosType getAngDistLens(){
		if(flag_switch_main_halo_on)
			return cosmo.angDist( main_halos[0]->getZlens());
		else{
			ERROR_MESSAGE();
			std::cout << "error, no main lens present" << std::endl;
			exit(1);
		}
	}

	/// remove all main halos
	void clearMainHalos(bool verbose);

	/// inserts a single main lens halo and adds it to the existing ones
	void insertMainHalo(LensHalo* halo,bool verbose = false);

	/// inserts a sequence of main lens halos and adds them to the existing ones
	void insertMainHalos(LensHalo** halos, std::size_t Nhalos,bool verbose = false);

	/// replaces existing main halos with a single main halo
	void replaceMainHalos(LensHalo* halo,bool verbose = false);
	/// replaces existing main halos with a sequence of main halos
	void replaceMainHalos(LensHalo** halos, std::size_t Nhalos,bool verbose = false);

  /** \brief Add substructures to the lens.
   
   This method is meant for inserting substructure to a main lens.  All the substructure will be at 
   one redshift.  The mass function follows a power law.  The density of substructures is constant within 
   a circular region.  The tidal truncation is controlled through the parameter density_contrast which is
   the average density within the substructures orbit in units of the average density to the universe at 
   the redshift where they are places.  For example density_contrast=200 would give them the truncation radius appropriate at R_200. 
   */
  void insertSubstructures(
        PosType Rregion            /// radius of region in which substructures are inserted (radians)
        ,PosType center[]          /// center of region in which the substructures are inserted (radians)
        ,PosType NumberDensity     /// number density per radian^2 of all substructures
        ,PosType Mass_min          /// minimum mass of substructures
        ,PosType Mass_max          /// maximum mass of substructures
        ,PosType redshift          /// redshift of substructures
        ,PosType alpha             /// index of mass function (dN/dm \propto m^alpha)
        ,PosType density_contrast  ///
        ,bool verbose
  );
  /** \brief This function will randomize the substructure without changing the region, mass function, etc.
   
   The Lens::insertSubstructures() function must have been called on this instance of the Lens before.
   */
  void resetSubstructure(bool verbose = false);
  
	/// get number of main halos
	std::size_t getNMainHalos() const;
	/// get number of main halos of given type
	template<typename HaloType>
	std::size_t getNMainHalos() const;
	
	/// get single main halo
	LensHalo* getMainHalo(std::size_t i);
	/// get single main halo of given type
	template<typename HaloType>
	HaloType* getMainHalo(std::size_t i);
	
	void rayshooterInternal(unsigned long Npoints, Point *i_points);
  void info_rayshooter(Point *i_point
                      ,std::vector<std::vector<double>> & ang_positions
                      ,std::vector<KappaType> & kappa_on_planes
                      ,std::vector<std::vector<LensHalo*>> & halo_neighbors
                      ,LensHalo *halo_max
                      ,KappaType &kappa_max
                      ,KappaType gamma_max[]
                      ,PosType rmax,short mode = 0  /// 0:physical distance, 1: comoving distance, 2: angular distance
                                       );

	// methods used for use with implanted sources

	short ResetSourcePlane(PosType z,bool nearest, unsigned long GalID=0, PosType *xx=NULL,bool verbose = false);

	/// Revert the source redshift to the value it was when the Lens was created.
	void RevertSourcePlane(){ toggle_source_plane = false;}
	//void ImplantSource(unsigned long index,CosmoHndl cosmo);
	PosType getSourceZ(){
		if(toggle_source_plane){
			return zs_implant;
		}else{
			return plane_redshifts.back();
		}
	}

	PosType getZmax(){return plane_redshifts.back();}

	/// print the cosmological parameters
	void PrintCosmology() { cosmo.PrintCosmology(); }
	
	/// returns the critical density at the main lens in Msun/ Mpc^2 for a source at zsource
	PosType getSigmaCrit(PosType zsource) { return cosmo.SigmaCrit(getZlens(), zsource); }
	

  /// returns a const reference to the cosmology so that constant functions can be used, but the cosmological parameters cannot be changed.
  const COSMOLOGY & getCosmo(){return cosmo;}
  
  /// set flag_switch_field_off, turn the field On/Off :
  void TurnFieldOff() { flag_switch_field_off = true ; }
  void TurnFieldOn() { flag_switch_field_off = false ; }
  
  /// get the field min mass :
  PosType getFieldMinMass() { return field_min_mass ; }
 
  // get the field_Off value :
  bool getfieldOff() {return flag_switch_field_off ;}
  
private:
	GLAMER_TEST_FRIEND(LensTest)
	
	// seed for random field generation
	long *seed;
  
  long init_seed;
  InputParams init_params;
  
  // the cosmology
	COSMOLOGY cosmo;

	/// field of view in square degrees
	PosType fieldofview;
	
	void readCosmology(InputParams& params);
	void assignParams(InputParams& params,bool verbose = false);
	
	/// turns source plane on and off
	bool toggle_source_plane;
	/// the distance from the source to the next plane
	PosType dDs_implant;
	PosType zs_implant,Ds_implant;
	/// This is the index of the plane at one larger distance than the new source distance
	int index_of_new_sourceplane;
	
	/// This is the source redshift that is read in from the parameter file and becomes the maximum redshift
	PosType zsource;
	
	void quicksort(LensHaloHndl *halo,PosType **pos,unsigned long N);
	
private: /* generation */
	/// create the lens planes
	void buildPlanes(InputParams& params, bool verbose);
	
	/// sets the distances and redshifts of the field planes equidistant
	void setFieldDist();
	/// load the redshifts of the field planes from a file
	void setFieldDistFromFile();
	/// setup the field plane distances
	void setupFieldPlanes();
	/// computes the distribution variables for field halos as specified in the parameter file
  /// this material was before computed in createFieldHalos
  void ComputeHalosDistributionVariables ();
	void createFieldHalos(bool verbose);
  
	/// read field halo data in from a file in Millennium output format
	void readInputSimFileMillennium(bool verbose);
	/// read field halo data in from a file in MultiDarkHalos output format
	void readInputSimFileMultiDarkHalos(bool verbose);
  /// read field halo data in from a file in Cabriel Caminha's input format 
  void readInputSimFileObservedGalaxies(bool verbose);
  
	/// build the field planes and sort halos onto them
	void createFieldPlanes(bool verbose);
	
	/// generate main halo from the parameter file
	void createMainHalos(InputParams& params);
	/// generate main halo from the parameter file
	void createMainPlanes();
	/// add a main halo to an existing plane, or create a new plane
	void addMainHaloToPlane(LensHalo* halo);
	
	/// combine field and main planes
	void combinePlanes(bool verbose);
	
  /* Variables used by buildPlanes, createFieldHalos, and createFieldPlanes */
  /// TO DO : should check the correctness of the following descriptions.
  
  /// number of redshift bins for mass function
  const int Nzbins = 64 ;
  /// number of mass bins for mass function
  const int Nmassbin=64;
  /// number of redshifts sampled for each bin
  int NZSamples = 50;
  // table for redshift bins for mass function
  std::vector<PosType> zbins ;
  // table of number of halos per redshift bins for mass function
  std::vector<PosType> NhalosbinZ ;
  /// same for the cumulative number density in one square degree
  std::vector<PosType> Nhaloestot_Tab ;
  /// averaged number of halos
  PosType aveNhalosField ;
  /// Log(mass) vector
  std::vector<PosType> Logm;
  /// Number of halos  field + substructure
  //std::size_t Nhalos ;
  /// table of halos bins for each sampled redshifts
  std::vector<std::vector<PosType>> NhalosbinMass;
  /// table for sigma_back in createFieldPlanes
  std::vector<PosType> sigma_back_Tab;
  
  /* ----- */
  
  // get the adress of field_plane_redshifts
  std::vector<PosType> & get_field_plane_redshifts () { return field_plane_redshifts ; }
  
  size_t getNFieldHalos() const {return field_halos.size();}
  size_t getNSubHalos() const {return substructure.halos.size();}
  
private: /* force calculation */
	/// if >= 1, deflection in the rayshooting is switched off
	bool flag_switch_deflection_off;
	bool flag_switch_lensing_off;
	
	/// the lensing planes
	std::vector<LensPlane *> lensing_planes;
	/// Dl[j = 0...] angular diameter distances, comoving
	std::vector<PosType> Dl;
	/// dDl[j] is the distance between plane j-1 and j plane, comoving
	std::vector<PosType> dDl;
	/// Redshifts of lens planes, 0...Nplanes.  Last one is the source redshift.
	std::vector<PosType> plane_redshifts;
	/// charge for the tree force solver (4*pi*G)
	PosType charge;
	
private: /* field */
	/// if true, the background is switched off and only the main lens is present
	bool flag_switch_field_off;
	
  /// vector of all field halos
  std::vector<LensHalo*> field_halos;
	/// original number of field planes
  std::size_t field_Nplanes_original;
  /// current number of field planes which may include substructure plane
  std::size_t field_Nplanes_current;
  
	/// vector of all field planes
	std::vector<LensPlane*> field_planes;
	/// vector of field plane redshifts
  std::vector<PosType> field_plane_redshifts;
  /// original field plane redshift
  std::vector<PosType> field_plane_redshifts_original;
  /// vector of field plane distances
  std::vector<PosType> field_Dl;
  /// original vector of field plane distances
  std::vector<PosType> field_Dl_original;
  
  struct SubStructureInfo{
    // things for substructures
    /// vector of all substructure halos
    std::vector<LensHalo*> halos;
    LensPlane *plane;
    PosType Rregion = 0;
    PosType Mmax = 0;
    PosType Mmin = 0;
    PosType alpha = 0;
    PosType Ndensity = 0;
    Point_2d center;
    PosType rho_tidal = 0;
    // Added quantities for the resetting of the substructure
    // (when WasInsertSubStructuresCalled = MAYBE) :
    PosType redshift = 0;
    bool verbose = false;
  };
  
  SubStructureInfo substructure;
  
  /// Flag to know if InsertSubStructures was called
  Boo WasInsertSubStructuresCalled = NO ;
  
	/// Perpendicular position of halo TODO: (In proper distance?)
	//PosType **halo_pos;
	
	/// type of mass function PS (0), ST (1), and power law (2) default is ST
	MassFuncType field_mass_func_type;
	/// slope of the mass function is field_mass_func_type == 2
	PosType mass_func_PL_slope;
	/// min mass for the halo model
	PosType field_min_mass;
	/// internal halo profile type; needs to be 0 or PowerLaw, 1 or NFW, 2 or PseudoNFW, 3 or NSIE, 4 or PointMass
	LensHaloType field_int_prof_type;
	/// power law or pseudo NFW internal profile slope
	PosType field_prof_internal_slope;
	
	/// if true, each field halo contains an NSIE galaxy inside it
	bool flag_field_gal_on;
	/// galaxy subhalo profile type; needs to be 0 or PowerLaw, 1 or NFW, 2 or PseudoNFW, 3 or NSIE, 4 or PointMass
	GalaxyLensHaloType field_int_prof_gal_type;
	// mass fraction in the host galaxy
	//PosType field_galaxy_mass_fraction;
	
	std::string redshift_planes_file;
	bool read_redshift_planes;
	
	std::string field_input_sim_file;
  HaloCatFormats field_input_sim_format;
  
	bool sim_input_flag;
	//std::string input_gal_file;
	//bool gal_input_flag;
	bool read_sim_file;
	
	/// increases are for cosmological mean number density of halos calculation
	PosType field_buffer;
	
private: /* main */
	/// having a main halo in the paramfile
	bool flag_switch_main_halo_on;
	
	/// vector of all main halos
	Utilities::MixedVector<LensHalo*> main_halos;
	/// vector of own main halos that will be deleted
	std::vector<LensHalo*> main_halos_created;
	/// vector of all main planes
	std::vector<LensPlane*> main_planes;
	/// vector of main plane redshifts
	std::vector<PosType> main_plane_redshifts;
	/// vector of main plane distances
	std::vector<PosType> main_Dl;
	
	/// main lens type
	LensHaloType main_halo_type;
	/// galaxy lens type
	GalaxyLensHaloType main_galaxy_halo_type;
	
private: /* input */
	/// file for multiple main halo input
	std::string pixel_map_input_file;
	
 	/// read main halos from a pixelized density map
  short pixel_map_on;
  /// zero padding for FFTs with pixelized density maps
  int pixel_map_zeropad;
  bool pixel_map_zeromean;
	void readPixelizedDensity();
  
  /// the center of the lens in spherical coordinates
  Utilities::Geometry::SphericalPoint central_point_sphere;
  /// optional angular radius of simulation cone that will be included
  PosType sim_angular_radius;
  /// inverse of the angular screening scale in the tree force calculation
  PosType inv_ang_screening_scale;
  
};

inline std::size_t Lens::getNMainHalos() const
{
	return main_halos.size();
}

template<typename HaloType>
inline std::size_t Lens::getNMainHalos() const
{
	return main_halos.size<HaloType>();
}

inline LensHalo* Lens::getMainHalo(std::size_t i)
{
	return main_halos.at(i);
}

template<typename HaloType>
inline HaloType* Lens::getMainHalo(std::size_t i)
{
	return main_halos.at<HaloType>(i);
}

typedef Lens* LensHndl;

#endif /* MULTIPLANE_H_ */
