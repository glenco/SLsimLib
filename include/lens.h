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


//  /* \brief Add substructures to the lens.
//
//   This method is meant for inserting substructure to a main lens.  All the substructure will be at
//   one redshift.  The mass function follows a power law.  The density of substructures is constant within
//   a circular region.  The tidal truncation is controlled through the parameter density_contrast which is
//   the average density within the substructures orbit in units of the average density to the universe at
//   the redshift where they are places.  For example density_contrast=200 would give them the truncation radius appropriate at R_200.
//   */
//  void insertSubstructures(
//        PosType Rregion            /// radius of region in which substructures are inserted (radians)
//        ,PosType center[]          /// center of region in which the substructures are inserted (radians)
//        ,PosType NumberDensity     /// number density per radian^2 of all substructures
//        ,PosType Mass_min          /// minimum mass of substructures
//        ,PosType Mass_max          /// maximum mass of substructures
//        ,PosType redshift          /// redshift of substructures
//        ,PosType alpha             /// index of mass function (dN/dm \propto m^alpha)
//        ,PosType density_contrast  ///
//        ,bool verbose
//  );
//
//  /** \brief This function will randomize the substructure without changing the region, mass function, etc.
//
//   The Lens::insertSubstructures() function must have been called on this instance of the Lens before.
//   */
//  void resetSubstructure(bool verbose = false);
//

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

  /** \brief Main routine for shooting rays in parrallel.
   
   i_points[].x must be set to the image postion in angular radians.
  
   i_points must have linked image points.  `LinkedPoint`s could be used, but for its internal use this is not done.
   */
	void rayshooterInternal(unsigned long Npoints    /// number of points to be shot
                          ,Point *i_points         /// points on the image plane
                          ,bool RSIverbose = false/// verbose option
                          );
  
  /** \brief Routine for shooting rays with differnt source redshifts in parrallel.
   
   i_points[].x must be set to the image postion in angular radians.
   */
  
  void rayshooterInternal(unsigned long Npoints   /// number of points to be shot
                          ,LinkedPoint *i_points        /// poinst on the image plane
                          ,std::vector<double> &source_zs /// source redshifts
                          ,bool RSIverbose = false/// verbose option
                          );

  void info_rayshooter(RAY &i_point
                      ,std::vector<Point_2d> & ang_positions
                      ,std::vector<KappaType> & kappa_on_planes
                      ,std::vector<std::vector<LensHalo*>> & halo_neighbors
                      ,LensHalo &halo_max
                      ,KappaType &kappa_max
                      ,KappaType gamma_max[]
                      ,PosType rmax  /// distance from ray on each plane, units depend on mode parameter
                      ,int tag = 0   /// i f not 0, information on halos with this tag are gathered
                      ,short mode = 0  /// 0:physical distance, 1: comoving distance, 2: angular distance
                      ,bool verbose = false
                                       );
  /** \brief  Find the image position of a source without grid refinement.
  
  This uses Powell's algorithm to minimise the distance between the source point of an image and the desired source point.  No grid is necessary.  This should be fast, but will miss multiple images.  This is useful for finding the position of weakly lensed images or the rough region where a grid should be put down for a strong lens.
  */
  /*RAY find_image(
          Point_2d y_source    /// input position of source (radians)
          ,Point_2d &x_image    /// initial guess for image postion (radians)
          ,PosType z_source     /// redshift of source
          ,PosType ytol2        /// target tolerance in source position squared
          ,PosType &fret        ///
          ,int sign=0             /// sign of magnification, it is found automatically if left out
  );*/
  
  /** \brief  Find the image position of a source without grid refinement.
  
  This finds an image position given a source postion.  No grid is necessary.
   
   If use_image_guess=true the input image position will be used as a first guess and the output image will be guarenteed to have the same pairity.
   
   This is useful for finding the position of weakly lensed images or for refining the image positions that are found on a finite grid.
   
  */
  RAY find_image(
            Point &p              ///  p.image->x should be set to source position
            ,PosType ytol2       /// target tolerance in source position squared
            ,PosType &dy2        /// final value of Delta y ^2
            ,bool use_image_guess /// if true p.x[] will be used as a guess for the image position
  );

  RAY find_image(
            RAY &p       /// p.y[] should be set to source position
            ,PosType ytol2        /// target tolerance in source position squared
            ,PosType &dy2        /// final value of Delta y ^2
            ,bool use_image_guess  // if true p.x[] will be used as a guess for the image position
  );

	// methods used for use with implanted sources

  ///  reset the redshift of the source plane
	short ResetSourcePlane(PosType z,bool nearest=false,bool verbose = false);

  ///  find information on the position of the source plane with respect to the lens planes
  void FindSourcePlane(
                          PosType zs                 /// redshift of implanted source
                          ,long &jmax                /// index of last plane at lower redshift
                          ,double &Dls               /// coordinate distance between last plane and source plane
                          ,double &Ds                /// total coordinate distance to source plane
  );

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
                          ,LensHaloType halo_type = LensHaloType::nfw_lens  /// type of halo
                          ,GalaxyLensHaloType galaxy_type = GalaxyLensHaloType::null_gal  /// type of galaxy, if null_gal no galaxy
                          ,double buffer = 1.0 /// buffer in Mpc for cone
                          ,bool verbose = false
                );
  
  Lens & operator=(Lens &&lens){
    
    fieldofview = lens.fieldofview;
    seed = lens.seed;
    init_seed = lens.init_seed;
    cosmo = lens.cosmo;
    toggle_source_plane = lens.toggle_source_plane;
  
    zs_implant = lens.zs_implant;
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
    
    //substructure = lens.substructure;
    
    //WasInsertSubStructuresCalled = lens.WasInsertSubStructuresCalled;
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

  template<typename P>
  void compute_rays_parallel(int start
                             ,int chunk_size
                             ,P *i_point
                             ,double *source_zs
                             ,bool multiZs
                             ,bool verbose = false);
	
	//void readCosmology(InputParams& params);
	//void assignParams(InputParams& params,bool verbose = false);
  void defaultParams(PosType zsource,bool verbose = true);
	
	/// turns source plane on and off
	bool toggle_source_plane;
 
  // redhsift of true source plane if rest from original
  PosType zs_implant;
	
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
  //size_t getNSubHalos() const {return substructure.halos.size();}
  
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
  std::vector<LensHalo *> field_halos;
	/// original number of field planes
  std::size_t field_Nplanes_original;
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
  
//  struct SubStructureInfo{
//    // things for substructures
//    /// vector of all substructure halos
//    std::vector<LensHaloNFW *> halos;
//    LensPlane *plane;
//    PosType Rregion = 0;
//    PosType Mmax = 0;
//    PosType Mmin = 0;
//    PosType alpha = 0;
//    PosType Ndensity = 0;
//    Point_2d center;
//    PosType rho_tidal = 0;
//    // Added quantities for the resetting of the substructure
//    // (when WasInsertSubStructuresCalled = MAYBE) :
//    PosType redshift = 0;
//    bool verbose = false;
//  };
//
//  SubStructureInfo substructure;
  
  /// Flag to know if InsertSubStructures was called
 // Boo WasInsertSubStructuresCalled = NO ;
  
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
  
  struct MINyFunction{
     MINyFunction(Lens &mylens,Point_2d y,int sign):lens(mylens),y(y),sign(sign),r2max(0){}

    double operator()(double *x){
      point.x[0] = x[1];
      point.x[1] = x[2];
      lens.rayshooterInternal(1,&point);
      double r2 = (y[0]-point.image->x[0])*(y[0]-point.image->x[0])
      + (y[1]-point.image->x[1])*(y[1]-point.image->x[1]);
      
      r2max = MAX(r2,r2max);
      return r2 + r2max*abs(sign - sgn(point.invmag()));
    }
    
    Lens &lens;
    Point_2d y;
    int sign;
    double r2max;
    LinkedPoint point;
  };
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


template <typename P>
void Lens::compute_rays_parallel(int start
                                 ,int chunk_size
                                 ,P *i_points
                                 ,double *source_zs
                                 ,bool multiZs
                                 ,bool verbose)
{
  int end        = start + chunk_size;

  int i, j;
  
  PosType xx[2];
  PosType aa,bb;
  PosType alpha[2];
  
  KappaType kappa,gamma[3];
  KappaType phi;
  
  Matrix2x2<PosType> G;

  PosType SumPrevAlphas[2];
  Matrix2x2<PosType> SumPrevAG;
  
  PosType *theta;
  
  long jmax = lensing_planes.size();
  double Dls_Ds; // this is the ratio between of the distance between the last lens plane and the source to the distance to the source
  double D_Ds; // this is the ratio between of the distance to the last lens plane and the source to the distance to the source

  if(!multiZs){
    if(source_zs[0] == plane_redshifts.back() ){
      Dls_Ds = dDl.back() / Dl.back();
      D_Ds = Dl[Dl.size() - 2] / Dl.back();
    }else{
      PosType Dls,Ds;
      FindSourcePlane(source_zs[0],jmax,Dls,Ds);
      Dls_Ds = Dls / Ds;
      if(jmax > 0) D_Ds = Dl[jmax-1] / Ds;
    }
  }
  
  // Main loop : loop over the points of the image
  for(i = start; i < end; i++)
  {
    // In case e.g. a temporary point is outside of the grid.
    if(i_points[i].in_image == MAYBE) continue;
    
    //theta = i_points[i].image->x;
    theta = i_points[i].ptr_y();
    
    theta[0] = i_points[i].x[0];
    theta[1] = i_points[i].x[1];

    // Initializing SumPrevAlphas :
    SumPrevAlphas[0] = theta[0];
    SumPrevAlphas[1] = theta[1];

    // Initializing SumPrevAG :
    SumPrevAG.setToI();
    
    // Setting phi on the first plane.
    phi = 0.0;
    
    // Default values :
    i_points[i].A.setToI();
    i_points[i].dt = 0;
    
    // In case we don't want to compute the values :
    if(flag_switch_lensing_off)
    {
      i_points[i].image->A.setToI();
      continue;
    }
    
    // Time delay at first plane : position on the observer plane is (0,0) => no need to take difference of positions.
    i_points[i].dt = 0;
    
    //0.5*( p->i_points[i].image->x[0]*p->i_points[i].image->x[0] + p->i_points[i].image->x[1]*p->i_points[i].image->x[1] )/ p->dDl[0] ;
    
    if(multiZs){
      PosType Dls,Ds;
      FindSourcePlane(source_zs[i],jmax,Dls,Ds);
      Dls_Ds = Dls / Ds;
      if(jmax > 0) D_Ds = Dl[jmax-1]/Ds;
    }
    
    // Begining of the loop through the planes :
    // Each iteration leaves i_point[i].image on plane (j+1)

    for(j = 0; j < jmax ; ++j)
      {
      
      double Dphysical = Dl[j]/(1 + plane_redshifts[j]);
      // convert to physical coordinates on the plane j, just for force calculation
      xx[0] = theta[0] *  Dphysical;
      xx[1] = theta[1] *  Dphysical;
      // PhysMpc = ComMpc / (1+z)
      
      assert(xx[0] == xx[0] && xx[1] == xx[1]);
      
      ////////////////////////////////////////////////////////
      
      lensing_planes[j]->force(alpha,&kappa,gamma,&phi,xx);
      // Computed in physical coordinates, xx is in PhysMpc.
      
      ////////////////////////////////////////////////////////
      
      assert(alpha[0] == alpha[0] && alpha[1] == alpha[1]);
      assert(gamma[0] == gamma[0] && gamma[1] == gamma[1]);
      assert(kappa == kappa);
      if(std::isinf(kappa)) { std::cout << "xx = " << xx[0] << " " << xx[1] << std::endl ;}
      assert(!std::isinf(kappa));
      
      G[0] = kappa + gamma[0];    G[1] = gamma[1];
      G[2] = gamma[1]; G[3] = kappa - gamma[0];
  
      /* multiply by fac to obtain 1/comoving_distance/physical_distance
       * such that a multiplication with the charge (in units of physical distance)
       * will result in a 1/comoving_distance quantity */
      
      G *= charge * Dl[j] / (1 + plane_redshifts[j]);
        
      assert(gamma[0] == gamma[0] && gamma[1] == gamma[1]);
      assert(kappa == kappa);
      assert(phi == phi);
      
      // This computes \vec{x}^{j+1} in terms of \vec{x}^{j}
      // according to the corrected Eq. (18) of paper GLAMER II ---------------------------------
      
      // Adding the j-plane alpha contribution to the sum \Sum_{k=1}^{j} \vec{alpha_j} :
      SumPrevAlphas[0] -= charge * alpha[0] ;
      SumPrevAlphas[1] -= charge * alpha[1] ;
      
      if(j < jmax-1 ){
        aa = dDl[j+1] / Dl[j+1];
        bb = Dl[j] / Dl[j+1];
      }else{
        aa = Dls_Ds;
        bb = D_Ds;
      }
      
        if(!flag_switch_deflection_off){
          theta[0] = bb * theta[0] + aa * SumPrevAlphas[0];
          theta[1] = bb * theta[1] + aa * SumPrevAlphas[1];
        }
   
      // ----------------------------------------------------------------------------------------
            
      // Sum_{k=1}^{j} Dl[k] A^k.G^k
      SumPrevAG -= (G * (i_points[i].A)) ;
      
      // Computation of the "plus quantities", i.e. the  next plane quantities :
      i_points[i].A = i_points[i].A * bb + SumPrevAG * aa;
      
      // ----------------------------------------------
      
      // Geometric time delay with added potential
      //p->i_points[i].dt += 0.5*( (xplus[0] - xminus[0])*(xplus[0] - xminus[0]) + (xplus[1] - xminus[1])*(xplus[1] - xminus[1]) ) * p->dTl[j+1] /p->dDl[j+1] /p->dDl[j+1] - phi * p->charge ; /// in Mpc  ???
      
      // Check that the 1+z factor must indeed be there (because the x positions have been rescaled, so it may be different compared to the draft).
      // Remark : Here the true lensing potential is not "phi" but "phi * p->charge = phi * 4 pi G".
      
      
    } // End of the loop going through the planes
    
    if(flag_switch_deflection_off){
      i_points[i].A = Matrix2x2<KappaType>::I() - SumPrevAG;
    }
    
    // Subtracting off a term that makes the unperturbed ray to have zero time delay
    //p->i_points[i].dt -= 0.5*( p->i_points[i].image->x[0]*p->i_points[i].image->x[0] + p->i_points[i].image->x[1]*p->i_points[i].image->x[1] ) / p->Dl[NLastPlane];
    
    // Conversion of dt from Mpc (physical Mpc) to Years -----------------
    i_points[i].dt *= MpcToSeconds * SecondToYears ;
    
    // ---------------------------------------------------------------------------------------------
    
    // Putting the final values of the quantities in the image point -----
    i_points[i].image->A = i_points[i].A;
    i_points[i].image->dt = i_points[i].dt;
    // ------------------------------------------------------------------------
    
    
    // TEST : showing final quantities
    // ------------------------------=
    if(verbose)  std::cout << "RSI final : X X | " << i << "  " << source_zs[i] << " | " << i_points[i].kappa() << " " << i_points[i].gamma1() << " " << i_points[i].gamma2() << " " << i_points[i].gamma3() << " " << i_points[i].invmag() << " | " << i_points[i].dt << std::endl ;
    
  } // End of the main loop.

}

#endif /* MULTIPLANE_H_ */

