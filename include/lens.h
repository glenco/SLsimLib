/*
 * multiplane.h
 *
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

class Lens{
public:
  Lens(long* seed, PosType z_source,CosmoParamSet cosmoset, bool verbose = false);
  //Lens(InputParams& params, long* my_seed, CosmoParamSet cosmoset, bool verbose = false);
  Lens(long* seed, PosType z_source,const COSMOLOGY &cosmo, bool verbose = false);
  //Lens(InputParams& params, long* my_seed, const COSMOLOGY &cosmo, bool verbose = false);
  
	~Lens();

	/// marks if the lens has been setup.
	bool set;

	/// the total number of lens planes
	int getNplanes() const {return lensing_planes.size();}
  
	/// field of view in square degrees
	PosType getfov() const {return fieldofview;};
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
	PosType getZlens() const{
		if(flag_switch_main_halo_on)
			return main_halos[0]->getZlens();
		else{
			ERROR_MESSAGE();
			std::cerr << "error, no main lens present" << std::endl;
			exit(1);
		}
	}
  /// Angular size distance (Mpc) to first main lens plane
	PosType getAngDistLens() const{
		if(flag_switch_main_halo_on)
			return cosmo.angDist( main_halos[0]->getZlens());
		else{
			ERROR_MESSAGE();
			std::cerr << "error, no main lens present" << std::endl;
			exit(1);
		}
	}
  
  Utilities::Geometry::SphericalPoint<> getCenter() const {return central_point_sphere;}

	/// remove all main halos
	void clearMainHalos(bool verbose=false);
  /// remaove all main halo of given type
  template<typename HaloType>
  void clearMainHalo(bool verbose=false);


	/// inserts a single main lens halo and adds it to the existing ones
  //void insertMainHalo(LensHalo* halo,PosType zlens, bool addplanes,bool verbose = false);
  
/*  void insertMainHalo(LensHalo *halo, bool addplanes,bool verbose)
  {
//    LensHaloNFW * halo = new LensHaloNFW(halo_in);
    halo->setCosmology(cosmo);
    main_halos.push_back(halo);
    
    flag_switch_main_halo_on = true;
    
    if(addplanes) addMainHaloToPlane(halo);
    else addMainHaloToNearestPlane(halo);
    
    combinePlanes(verbose);
  }
*/
  
  template <typename T>
  void insertMainHalo(const T &halo_in, bool addplanes,bool verbose=false)
  {
    
    T * halo = new T(halo_in);
    halo->setCosmology(cosmo);
    main_halos.push_back(halo);
    
    flag_switch_main_halo_on = true;
    
    if(addplanes) addMainHaloToPlane(halo);
    else addMainHaloToNearestPlane(halo);
    
    combinePlanes(verbose);
  }
  
  /** \brief This has the same effect as insertMainHalo(), but the halo is not
   copied, it is moved.
  
   The Lens will take possession of the halo and will destroy it when it is destroyed.
   This is to avoid copying halos that take up a lot of memory and require
   a lot of time to copy like LensHaloParticles().
  **/
  template <typename T>
  void moveinMainHalo(T &halo_in, bool addplanes,bool verbose=false)
  {
    T * halo = new T(std::move(halo_in));
    halo->setCosmology(cosmo);
    main_halos.push_back(halo);
    
    flag_switch_main_halo_on = true;
    
    if(addplanes) addMainHaloToPlane(halo);
    else addMainHaloToNearestPlane(halo);
    
    combinePlanes(verbose);
  }
  
  /**
   * \brief Inserts a single main lens halo and deletes all previous ones.
   * Then all lensing planes are updated accordingly.
   *
   * Note that this does delete all the halos that were there.
   */
  template <typename T>
  void replaceMainHalo(const T &halo_in,bool addplanes,bool verbose=false)
  {
    //Utilities::delete_container(main_halos);   // ????
    /*while(main_halos.size() > 0){
      delete main_halos.back();
      main_halos.pop_back();
    }*/
    
    main_halos.clear();  // ???? is this a memory leak ????
    
    T * halo = new T(halo_in);
    halo->setCosmology(cosmo);
    main_halos.push_back(halo);
    
    flag_switch_main_halo_on = true;
    
    Utilities::delete_container(main_planes);
    createMainPlanes();
    combinePlanes(verbose);
  }

  /**
   * \brief Inserts a sequense of main lens halos and adds them to the existing ones.
   * Then all lensing planes are updated accordingly.
   * If addplanes is true new planes will be added otherwise
   * the halo is added to the nearest plane and a plane is added only
   * if none exited on entry.
   *
   *  The angular position of the halo should be preserved, but the x coordinates may change
   *  The halos are copied so the input halos can be destoyed without affecting the Lens.
   */
  template <typename T>
  void insertMainHalos(std::vector<T> &my_halos,bool addplanes, bool verbose=false)
  {
    T* ptr;
    //for(std::size_t i = 0; i < my_halos.size() ; ++i)
    for(T &h : my_halos){
      ptr = new T(h);
      ptr->setCosmology(cosmo);
      ptr->setDist(cosmo);
      main_halos.push_back( ptr );
      if(addplanes) addMainHaloToPlane( ptr );
      else addMainHaloToNearestPlane( ptr );
    }
    
    flag_switch_main_halo_on = true;
    
    combinePlanes(verbose);
  }
  /**
   * \brief Inserts a sequense of main lens halos and remove all previous ones.
   *
   * Note that this does delete the halos that were there.
   * Then all lensing planes are updated accordingly.
   */
  template <typename T>
  void replaceMainHalos(std::vector<T> &my_halos,bool verbose)
  {
    Utilities::delete_container(main_halos);   // ????

    T* ptr;
    //for(std::size_t i = 0; i < my_halos.size() ; ++i)
    for(T &h : my_halos){
      ptr = new T(h);
      ptr->setCosmology(cosmo);
      ptr->setDist(cosmo);
      main_halos.push_back( ptr );
    }
    
    flag_switch_main_halo_on = true;
    
    Utilities::delete_container(main_planes);
    createMainPlanes();
    combinePlanes(verbose);
  }

	/// inserts a sequence of main lens halos and adds them to the existing ones
	//void insertMainHalos(LensHalo** halos, std::size_t Nhalos,bool addplanes,bool verbose = false);

	/// replaces existing main halos with a single main halo
  //void replaceMainHalo(LensHalo* halo,PosType zlens, bool addplanes,bool verbose = false);

	/// replaces existing main halos with a sequence of main halos
	// replaceMainHalos(LensHalo** halos, std::size_t Nhalos,bool verbose = false);

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

  /**\brief Using to shoot a single ray
   
   ray.x should be set to the image position.
   The kappa,gamma,deflection, time-delay and
   source position will be calculated at that
   image point.
   */
  void rayshooter(RAY &ray);
	void rayshooterInternal(unsigned long Npoints, Point *i_points, bool RSIverbose = false);
  void info_rayshooter(Point *i_point
                      ,std::vector<Point_2d> & ang_positions
                      ,std::vector<KappaType> & kappa_on_planes
                      ,std::vector<std::vector<LensHalo*>> & halo_neighbors
                      ,LensHalo **halo_max
                      ,KappaType &kappa_max
                      ,KappaType gamma_max[]
                      ,PosType rmax,short mode = 0  /// 0:physical distance, 1: comoving distance, 2: angular distance
                      ,bool verbose = false
                                       );

	// methods used for use with implanted sources

  ///  reset the redshift of the source plane
	short ResetSourcePlane(PosType z,bool nearest=false, unsigned long GalID=0, PosType *xx=NULL,bool verbose = false);

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

	PosType getZmax() const{return plane_redshifts.back();}

	/// print the cosmological parameters
	void PrintCosmology() { cosmo.PrintCosmology(); }
	
	/// returns the critical density at the main lens in Msun/ Mpc^2 for a source at zsource
	PosType getSigmaCrit(PosType zsource) const{ return cosmo.SigmaCrit(getZlens(), zsource); }
	

  /// returns a const reference to the cosmology so that constant functions can be used, but the cosmological parameters cannot be changed.
  const COSMOLOGY & getCosmo(){return cosmo;}
  
  /// set flag_switch_field_off, turn the field On/Off :
  void TurnFieldOff() { flag_switch_field_off = true ; }
  void TurnFieldOn() { flag_switch_field_off = false ; }
  
  /// get the field min mass :
  PosType getFieldMinMass() const { return field_min_mass ; }
 
  // get the field_Off value :
  bool getfieldOff() const {return flag_switch_field_off ;}
  
  /** \brief Add random halos to the light cone according to standard structure formation theory.  A new realization of the light-cone can be made with Lens::resetFieldHalos() after this function is called once.
   
   The cone is filled up until the redshift of the current zsource that is stored in the Lens class.  The field is a circular on the sky.  There is no clustering of the halos.
   */
  void GenerateFieldHalos(double min_mass /// minimum mass of halos
                          ,MassFuncType mass_function /// type of mass function
                          ,double field_of_view  /// in square degrees
                          ,int Nplanes           /// number of lens planes
                          ,LensHaloType halo_type = nfw_lens  /// type of halo
                          ,GalaxyLensHaloType galaxy_type = null_gal  /// type of galaxy, if null_gal no galaxy
                          ,double buffer = 1.0 /// buffer in Mpc for cone
                          ,bool verbose = false
                );
  
  Lens & operator=(Lens &&lens){
    
    fieldofview = lens.fieldofview;
    seed = lens.seed;
    init_seed = lens.init_seed;
    cosmo = lens.cosmo;
    toggle_source_plane = lens.toggle_source_plane;
    dDs_implant = lens.dDs_implant;
    dTs_implant = lens.dTs_implant;
    zs_implant = lens.zs_implant;
    Ds_implant = lens.Ds_implant;
    index_of_new_sourceplane = lens.index_of_new_sourceplane;
    zsource = lens.zsource;

    NZSamples = lens.NZSamples;
    zbins = lens.zbins;
    NhalosbinZ = lens.NhalosbinZ;
    Nhaloestot_Tab = lens.Nhaloestot_Tab;
    aveNhalosField = lens.aveNhalosField;
    Logm = lens.Logm;
    NhalosbinMass = lens.NhalosbinMass;
    sigma_back_Tab = lens.sigma_back_Tab;
    
    flag_switch_deflection_off = lens.flag_switch_deflection_off;
    flag_switch_lensing_off = lens.flag_switch_lensing_off;
    
    Dl = lens.Dl;
    dDl = lens.dDl;
    dTl = lens.dTl;
    plane_redshifts = lens.plane_redshifts;
    charge = lens.charge;
    flag_switch_field_off = lens.flag_switch_field_off;
    
    field_halos = lens.field_halos;

    field_Nplanes_original = lens.field_Nplanes_original;
    field_Nplanes_current = lens.field_Nplanes_current;
    
    field_plane_redshifts = lens.field_plane_redshifts;
    field_plane_redshifts_original = lens.field_plane_redshifts_original;
    field_Dl = lens.field_Dl;
    field_Dl_original = lens.field_Dl_original;
    
    substructure = lens.substructure;
    
    WasInsertSubStructuresCalled = lens.WasInsertSubStructuresCalled;
    field_mass_func_type = lens.field_mass_func_type;
    mass_func_PL_slope = lens.mass_func_PL_slope;
    field_min_mass = lens.field_min_mass;
    field_int_prof_type = lens.field_int_prof_type;
    field_prof_internal_slope = lens.field_prof_internal_slope;
    
    flag_field_gal_on = lens.flag_field_gal_on;
    field_int_prof_gal_type = lens.field_int_prof_gal_type;
    field_int_prof_gal_slope = lens.field_int_prof_gal_slope;
    
    redshift_planes_file = lens.redshift_planes_file;
    read_redshift_planes = lens.read_redshift_planes;
    
    field_input_sim_file = lens.field_input_sim_file;
    field_input_sim_format = lens.field_input_sim_format;
    
    sim_input_flag = lens.sim_input_flag;
    read_sim_file = lens.read_sim_file;
    field_buffer = lens.field_buffer;
    
    flag_switch_main_halo_on = lens.flag_switch_main_halo_on;
    
    main_plane_redshifts = lens.main_plane_redshifts;
    main_Dl = lens.main_Dl;
    
    main_halo_type = lens.main_halo_type;
    main_galaxy_halo_type = lens.main_galaxy_halo_type;

    pixel_map_input_file = lens.pixel_map_input_file;
    pixel_map_on = lens.pixel_map_on;
    pixel_map_zeropad = lens.pixel_map_zeropad;
    pixel_map_zeromean = lens.pixel_map_zeromean;
    
    central_point_sphere = lens.central_point_sphere;
    sim_angular_radius = lens.sim_angular_radius;
    inv_ang_screening_scale = lens.inv_ang_screening_scale;
    
    std::swap(lensing_planes,lens.lensing_planes);
    std::swap(field_planes,lens.field_planes);
    swap(main_halos,lens.main_halos);  /// MixedVector cannot be copyed
    std::swap(main_planes,lens.main_planes);
    
    return *this;
  }
  Lens(Lens &&lens){
    *this = std::move(lens);
  }
  
private:
  Lens & operator=(const Lens &lens);   // block copy
  Lens(const Lens &lens);
  
protected:
  /// field of view in square degrees
  PosType fieldofview;

private:
	GLAMER_TEST_FRIEND(LensTest)
	
	// seed for random field generation
	long *seed;
  
  long init_seed;
  //InputParams init_params;
  
  // the cosmology
	COSMOLOGY cosmo;

	
	//void readCosmology(InputParams& params);
	//void assignParams(InputParams& params,bool verbose = false);
  void defaultParams(PosType zsource,bool verbose = true);
	
	/// turns source plane on and off
	bool toggle_source_plane;
  /// the distance from the source to the next plane
  PosType dDs_implant;
  /// the distance from the source to the next plane
  PosType dTs_implant;
	PosType zs_implant,Ds_implant;
	/// This is the index of the plane at one larger distance than the new source distance
	int index_of_new_sourceplane;
	
	/// This is the source redshift that is read in from the parameter file and becomes the maximum redshift
	PosType zsource;
	
	void quicksort(LensHaloHndl *halo,PosType **pos,unsigned long N);
	
	// create the lens planes
	//void buildPlanes(InputParams& params, bool verbose);
	
	/// sets the distances and redshifts of the field planes equidistant
	void setFieldDist();
	/// load the redshifts of the field planes from a file
	void setFieldDistFromFile();
	/// setup the field plane distances
	void setupFieldPlanes();
	/// computes the distribution variables for field halos as specified in the parameter file
  /// this material was before computed in createFieldHalos
  void ComputeHalosDistributionVariables ();
  
  enum DM_Light_Division {All_DM,Moster};

	void createFieldHalos(bool verbose,DM_Light_Division division = Moster);
  
	/// read field halo data in from a file in Millennium output format
	void readInputSimFileMillennium(bool verbose,DM_Light_Division division = Moster);
  
	/// read field halo data in from a file in MultiDarkHalos output format
	void readInputSimFileMultiDarkHalos(bool verbose,DM_Light_Division division = Moster);
  
  /// read field halo data in from a file in Cabriel Caminha's input format 
  void readInputSimFileObservedGalaxies(bool verbose);
  
	/// build the field planes and sort halos onto them
	void createFieldPlanes(bool verbose);
	
	// generate main halo from the parameter file
	//void createMainHalos(InputParams& params);
	/// generate main halo from the parameter file
	void createMainPlanes();
	/// add a main halo to an existing plane, or create a new one plane if it is not close enough
	void addMainHaloToPlane(LensHalo* halo);
  void addMainHaloToNearestPlane(LensHalo* halo);

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
  
  /// Number of Field Halos
  size_t getNFieldHalos() const {return field_halos.size();}
  /// Number of Sub Halos
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
  /// dTl[j] is the lookback-time between plane j-1 and j plane, comoving
  std::vector<PosType> dTl;
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
  /// if field galaxy PowerLaw is used a slope must be assigned 
  PosType field_int_prof_gal_slope;

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
  Utilities::Geometry::SphericalPoint<PosType> central_point_sphere;
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
  if(main_halos.size<HaloType>() == 0 ) return nullptr;
	return main_halos.at<HaloType>(i);
}

typedef Lens* LensHndl;

#endif /* MULTIPLANE_H_ */

