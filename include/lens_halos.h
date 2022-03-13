/*
 * lens_halos.h
 *
 *  Created on: 06.05.2013
 */

#ifndef LENS_HALOS_H_
#define LENS_HALOS_H_

//#include "standard.h"
#include "InputParams.h"
//#include "source.h"
//#include "point.h"
#include "quadTree.h"
#include "particle_types.h"
#include "image_processing.h"

#include <complex>
#include <complex.h>
#ifdef ENABLE_CERF
#include <cerf.h>
#endif

/**
 * \brief A base class for all types of lensing "halos" which are any mass distribution that cause lensing.
 *
 * It contains the mass, maximum radius (Rsize), and scale radius (rscale) of the halo,
 * as well as the redshift (zlens).
 *
 * It has get and set functions for the members as well as virtual functions like:
 * force_halo
 * that compute the lensing properties -- deflection, convergence, and shear.
 *
 * Along with the simple set function, there are two general initialization functions,
 * that calculate the rest of the properties based on some input lens halo parameter
 * (e.g. mass).
 *
 * initFromFile
 * is intended to be used when the halo data is read in from a simulation
 * file. In this case the general halo is assumed to be an NFW halo and therefore the
 * maximum velocity and the half-mass radius need to be set. This function is overridden
 * in derived classes and in cases where it is not applicable only the mass is taken
 * into initializing the lensing halo.
 *
 * initFromMassFunc
 * is intended for the cases where the simulation is populated by lensing halos from
 * a mass function. Then one needs all parameters of the halo -- mass, Rsize, and rscale.
 */

class LensHalo{
public:
	LensHalo();
  LensHalo(PosType z,const COSMOLOGY &cosmo);
  //LensHalo(InputParams& params,COSMOLOGY &cosmo,bool needRsize = true);
  //LensHalo(InputParams& params,bool needRsize = true);
  LensHalo(const LensHalo &h);
  
  LensHalo(LensHalo &&h){
//    LensHalo(LensHalo &&h):star_tree(nullptr){
    *this = std::move(h);
  }

	virtual ~LensHalo();
  
  LensHalo & operator=(const LensHalo &h);
  LensHalo & operator=(LensHalo &&h);

  
  /** get the Rmax which is larger than Rsize in Mpc.  This is the region exterior to which the
   halo will be considered a point mass.  Between Rsize and Rmax the deflection and shear are interpolated.
   */
  inline float get_Rmax() const { return Rmax; }
  /// get the Rsize which is the size of the halo in Mpc
  inline float getRsize() const { return Rsize; }
	/// get the mass solar units
  inline float get_mass() const { return mass; }
  /// get the scale radius in Mpc
	inline float get_rscale() const { return rscale; }
	/// get the redshift
	inline PosType getZlens() const { return zlens; }
    
  // set the position of the Halo in physical Mpc on the lens plane
  //void setX(PosType PosX, PosType PosY) { posHalo[0] = PosX ; posHalo[1] = PosY ; }
  // set the position of the Halo in physical Mpc on the lens plane
  //void setX(PosType *PosXY) { posHalo[0] = PosXY[0] ; posHalo[1] = PosXY[1] ; }
  
  /// get the position of the Halo in physical Mpc on the lens plane
  void getX(PosType * MyPosHalo) const {
    assert(Dist != -1);
    MyPosHalo[0] = posHalo[0]*Dist;
    MyPosHalo[1] = posHalo[1]*Dist;
  }

  /// returns position of the Halo in physical Mpc on the lens plane
  PosType operator[](int i) const{return posHalo[i]*Dist;}
  
  /// set the position of the Halo in radians
  void setTheta(PosType PosX, PosType PosY) { posHalo[0] = PosX ; posHalo[1] = PosY ; }
  /// set the position of the Halo in radians
  void setTheta(PosType *PosXY) { posHalo[0] = PosXY[0] ; posHalo[1] = PosXY[1] ; }
  /// set the position of the Halo in radians
  void setTheta(const Point_2d &p) { posHalo[0] = p[0] ; posHalo[1] = p[1] ; }
  /// get the position of the Halo in radians
  void getTheta(PosType * MyPosHalo) const { MyPosHalo[0] = posHalo[0] ; MyPosHalo[1] = posHalo[1]; }
  
  /// Set the angular size distance to the halo.  This should be the distance to the lens plane.
  void setDist(COSMOLOGY &co){Dist = co.angDist(zlens);}

  /// return current angular size distance, ie conversion between angular and special coordinates.  This may not agree with
  /// the getZ() value because of the projection onto the lens plane.
  PosType getDist() const {return Dist;}

  /// get the position of the Halo in physical Mpc
  /// display the position of the halo
  void displayPos() { std::cout << "Halo PosX = " << posHalo[0] << " ; Halo PosY = " << posHalo[1] << std::endl; }
  
	/// initialize from a simulation file
	virtual void initFromFile(float my_mass, long *seed, float vmax, float r_halfmass){};
	/// initialize from a mass function
	virtual void initFromMassFunc(float my_mass, float my_Rsize, float my_rscale, PosType my_slope, long *seed);
  
  /// set Rsize (in Mpc) and reset Rmax
  virtual void set_RsizeRmax(float my_Rsize){Rmax = Rmax*my_Rsize/Rsize; Rsize = my_Rsize; xmax = Rsize/rscale;};
	/// set mass (in solar masses)
	void set_mass(float my_mass){mass=my_mass;};
	/// set scale radius (in Mpc)
	virtual void set_rscale(float my_rscale){rscale=my_rscale; xmax = Rsize/rscale;};
	/// set redshift
  //void setZlens(PosType my_zlens){zlens=my_zlens; Dist=-1;}
  void setZlens(PosType my_zlens,const COSMOLOGY &cosmo){
    zlens=my_zlens;
    Dist=cosmo.angDist(my_zlens);
  }
  void setRsize(PosType R){Rsize = R;}
  
  // ste redshift and distance
  void setZlensDist(PosType my_zlens,const COSMOLOGY &cos){
    zlens=my_zlens;
    Dist = cos.angDist(zlens);
  }
  void setMass(PosType m){mass = m;}
  
  
  /// set slope
	virtual void set_slope(PosType my_slope){beta=my_slope;};
  /// get slope
  virtual PosType get_slope(){return beta;};
  /// flag=True if halo elliptical
  bool get_flag_elliptical(){return elliptical_flag;};
  void set_flag_elliptical(bool ell){elliptical_flag=ell;};
  bool get_switch_flag(){return switch_flag;}; /// flag permits case distinction in force_halo_asym for elliptical NFWs only (get_switch_flag==true), in latter case the mass_norm_factor^2 is used instead of mass_norm_factor.
  void set_switch_flag(bool swt){switch_flag=swt;}; /// used for elliptical NFWs only, in that case get_switch_flag==true
  
  
  /// set cosmology for halo
	virtual void setCosmology(const COSMOLOGY& cosmo) {}
	
	/* calculate the lensing properties -- deflection, convergence, and shear
   if not overwritten by derived class it uses alpha_h(), gamma_h(), etc. of the
   derived case or for a point source in this class
   
   xcm - the physical position on the lens plane relative to the center of the LensHalo in Mpc
  Units :
   
  ALPHA    -    mass/PhysMpc - ALPHA / Sig_crit / Dl is the deflection in radians
  KAPPA    -    surface mass density , mass / /PhysMpc/PhysMpc - KAPPA / Sig_crit is the convergence
  GAMMA    -    mass / /PhysMpc/PhysMpc - GAMMA / Sig_crit is the shear
  PHI      -    mass - PHI / Sig_crit is the lensing potential
   
* returns the lensing quantities of a ray in center of mass coordinates.
   *
   *  Warning: This adds to input value of alpha, kappa, gamma, and phi.  They need
   *  to be zeroed out if the contribution of just this halo is wanted.
   */
  virtual void force_halo(
      PosType *alpha          /// deflection solar mass/Mpc
      ,KappaType *kappa     /// surface density in units of Sigma crit
      ,KappaType *gamma     /// three components of shear
      ,KappaType *phi       /// potential in solar masses
      ,PosType const *xcm   /// position relative to center (in physical Mpc)
      ,bool subtract_point=false /// if true contribution from a point mass is subtracted
      ,PosType screening=1.0   /// the factor by which to scale the mass for screening of the point mass subtraction
  );

	/// force tree calculation for stars
	//void force_stars(PosType *alpha,KappaType *kappa,KappaType *gamma,PosType const *xcm);
  
	/// internal compare redshift function
	bool compareZ(PosType z){return z > zlens;};
  
  /// stars
//  bool AreStarsImplanted() const {return stars_implanted;}
//  void implant_stars(PosType **centers,int Nregions,long *seed, IMFtype type=One);
//  void implant_stars(PosType *center,long *seed,IMFtype type = One);
//  //void implant_stars(PosType *x,PosType *y,int Nregions,long *seed,IMFtype type=One);
//  void remove_stars();
//  IMFtype getStarIMF_type() const {return main_stars_imf_type;}
//  /// Fraction of surface density in stars
//  PosType getFstars() const {return star_fstars;}
//  /// The mass of the stars if they are all the same mass
//  PosType getStarMass() const {if(stars_implanted)return stars_xp[0].mass(); else return 0.0;}
  
  /// the method used to ellipticize a halo if fratio!=1 and halo is not NSIE
  EllipMethod getEllipMethod() const {return main_ellip_method;}
  /// get vector of Fourier modes, which are calculated in the constructors of the LensHaloes when main_ellip_method is set to 'Fourier'
  std::vector<double> get_mod() { std::vector<double> fmodes(Nmod); for(int i=0;i<Nmod;i++){fmodes[i]= mod[i] ;}  ;return fmodes;}
  /// get length of mod array, which is Nmod. Not to be confused with getNmodes in the class LensHaloFit
  const static int get_Nmod() {return Nmod;}
  
	/// get the number of halo parameters
	virtual std::size_t Nparams() const;
	/// get the value of a scaled halo parameter by index
	virtual PosType getParam(std::size_t p) const;
	/// set the value of a scaled halo parameter by index
	virtual PosType setParam(std::size_t p, PosType value);
	
	/// print the halo parameters in CSV format
	virtual void printCSV(std::ostream&, bool header = false) const;
  
	/// Prints star parameters; if show_stars is true, prints data for single stars
	//void PrintStars(bool show_stars);
  
  PosType MassBy2DIntegation(PosType R);
  PosType MassBy1DIntegation(PosType R);
  PosType test_average_gt(PosType R);
  PosType test_average_kappa(PosType R);
  
  // renomalize to make mass match
  
  void set_norm_factor(){mass_norm_factor=1;mass_norm_factor=mass/MassBy1DIntegation(Rsize);}
  
  /// set radius rsize beyond which interpolation values between alpha_ellip and alpha_iso are computed
  void set_rsize(float my_rsize){ Rsize = my_rsize;};
	float get_rsize(){return Rsize;};

  // all of the following functions were used for Ansatz III w derivatives of the Fourier modes
  
  /// perform some basic consistancy checks for halo
  bool test();
  
  size_t getID() const {return idnumber;}
  void setID(size_t id){idnumber = id;}
  
  PosType renormalization(PosType r_max);
 
  /** \brief Map a PixelMap of the surface, density, potential and potential gradient
   centred on (0,0) in LensHalo coordinates

   Units :
   ALPHA    -    mass/PhysMpc - ALPHA / Sig_crit / Dl is the deflection in radians
   KAPPA    -    surface mass density , mass / /PhysMpc/PhysMpc - KAPPA / Sig_crit is the convergence
   GAMMA    -    mass / /PhysMpc/PhysMpc - GAMMA / Sig_crit is the shear
   PHI      -    mass - PHI / Sig_crit / Dl / Dl  is the lensing potential whose angular gradient is the deflection and angular Laplacian is 2 times the convergence
 
   centred on (0,0) in LensHalo coordinates
   
   */

  PixelMap map_variables(
                         LensingVariable lensvar /// lensing variable - KAPPA, ALPHA1, ALPHA2, GAMMA1, GAMMA2 or PHI
                         ,size_t Nx
                         ,size_t Ny
                         ,double res             /// angular resolution in radians
  );
  
private:
  size_t idnumber; /// Identification number of halo.  It is not always used.
  /// Position of the Halo in angle
  PosType posHalo[2];
  PosType zlens;

protected:
  // This is the size of the halo beyond which it does not have the expected profile.
  float Rsize = 0;

  // total mass in Msun
  float mass;
  // angular size distance to this lens
  PosType Dist;
  PosType mnorm;

  // Beyond Rmax the halo will be treated as a point mass.  Between Rsize and Rmax the deflection
  // and shear are interpolated.  For circularly symmetric lenses Rmax should be equal to Rsize
  float Rmax;

 
  PosType alpha_int(PosType x) const;
  PosType norm_int(PosType r_max);

  void force_halo_sym(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType const *xcm,bool subtract_point=false,PosType screening = 1.0);
  void force_halo_asym(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType const *xcm,bool subtract_point=false,PosType screening = 1.0);
  
  bool force_point(PosType *alpha,KappaType *kappa,KappaType *gamma
                   ,KappaType *phi,PosType const *xcm,PosType rcm2
                   ,bool subtract_point,PosType screening);
  
  struct norm_func{
    norm_func(LensHalo& halo, PosType my_r_max): halo(halo), r_max(my_r_max){};
    LensHalo& halo;
    PosType r_max;
    //PosType operator ()(PosType theta) {halo.alpha_asym(r_max, theta, alpha_arr); return alpha_arr[0]*cos(theta)*cos(theta)+alpha_arr[1]*sin(theta)*sin(theta);}
    PosType operator ()(PosType theta) {return halo.alpha_ell(r_max, theta);}
  };
  
  struct Ialpha_func{
    Ialpha_func(LensHalo& halo): halo(halo){};
    LensHalo& halo;
    PosType operator ()(PosType x) {return halo.alpha_h(x)/x ;}
  };
  
  struct Ig_func{
    Ig_func(const LensHalo& halo): halo(halo){};
    const LensHalo& halo;
    PosType operator ()(PosType x) {return halo.gfunction(x)/x ;}
  };
  
//  std::vector<IndexType> stars_index;
//  std::vector<StarType> stars_xp;
//  TreeQuadParticles<StarType> *star_tree;
//
//  int stars_N;
//  PosType star_massscale;
//  /// star masses relative to star_massscles
//  //float *star_masses;
//  PosType star_fstars;
//  PosType star_theta_force;
//  int star_Nregions;
//  std::vector<PosType> star_region;
  PosType beta;
//  void substract_stars_disks(PosType const *ray,PosType *alpha
//                             ,KappaType *kappa,KappaType *gamma);
//  std::vector<float> stellar_mass_function(IMFtype type, unsigned long Nstars, long *seed, PosType minmass=0.0, PosType maxmass=0.0
//                               ,PosType bendmass=0.0, PosType powerlo=0.0, PosType powerhi=0.0);
  
  
  /// read in parameters from a parameterfile in InputParams params
  void assignParams(InputParams& params,bool needRsize);
  /// read in star parameters. This is valid for all halos and not overloaded.
  //void assignParams_stars(InputParams& params);
  
  /// error message printout
  void error_message1(std::string name,std::string filename);
  
  
  /// The factor by which Rmax is larger than Rsize
  float Rmax_to_Rsize_ratio = 1.2;
  
  /// scale length or core size.  Different meaning in different cases.  Not used in NSIE case.
  float rscale;
  // redshift
  //PosType zlens;

  EllipMethod main_ellip_method;

//  bool stars_implanted;
//  /// Number of regions to be subtracted to compensate for the mass in stars
//  IMFtype main_stars_imf_type;
//  PosType main_stars_min_mass;
//  PosType main_stars_max_mass;
//  PosType bend_mstar;
//  PosType lo_mass_slope;
//  PosType hi_mass_slope;
//  /// parameters for stellar mass function: minimal and maximal stellar mass, bending point for a broken power law IMF
//  std::vector<PosType> star_Sigma;
//  std::vector<Point_2d> star_xdisk;
//
  
  /// point mass case
  /// r |alpha(r)| pi Sigma_crit / Mass
  virtual PosType inline alpha_h(PosType x) const {return -1;};
  virtual KappaType inline kappa_h(PosType x) const {return 0;};
  virtual KappaType inline gamma_h(PosType x) const {return -2;};
  virtual KappaType inline phi_h(PosType x) const {return 1;};
  virtual KappaType inline phi_int(PosType x) const {return 1;};
  virtual PosType inline ffunction(PosType x)const {return 0;};
  virtual PosType inline gfunction(PosType x) const {return -1;};
  virtual PosType inline dgfunctiondx(PosType x){return 0;};
  virtual PosType inline bfunction(PosType x){return -1;};
  virtual PosType inline dhfunction(PosType x) const {return 1;};
  virtual PosType inline ddhfunction(PosType x, bool numerical){return 0;};
  virtual PosType inline dddhfunction(PosType x, bool numerical){return 0;};
  virtual PosType inline bnumfunction(PosType x){return -1;};
  virtual PosType inline dbfunction(PosType x){return 0;};
  virtual PosType inline ddbfunction(PosType x){return 0;};
  virtual PosType inline dmoddb(int whichmod, PosType q, PosType b){return 0;};
  virtual PosType inline ddmoddb(int whichmod, PosType q, PosType b){return 0;};
  virtual PosType inline dmoddq(int whichmod, PosType q, PosType b){return 0;};
  virtual PosType inline ddmoddq(int whichmod, PosType q, PosType b){return 0;};
  
  PosType xmax;  /// This is Rsize/rscale !!
  
  PosType mass_norm_factor=1;
  
  // Functions for calculating axial dependence
  float pa;
  float fratio=1.0;
  bool elliptical_flag = false;
  bool switch_flag = false; /// If set to true the correct normalization is applied for asymmetric NFW profiles, the mass_norm_factor is different for the other halos. 
  
  
  void faxial(PosType x,PosType theta,PosType f[]);
  void faxial0(PosType theta,PosType f0[]);
  void faxial1(PosType theta,PosType f1[]);
  void faxial2(PosType theta,PosType f2[]);
  void gradial(PosType r,PosType g[]);
  void gradial2(PosType r,PosType mu, PosType sigma,PosType g[]);
  
  void felliptical(PosType x, PosType q, PosType theta, PosType f[], PosType g[]);
  
  virtual void gamma_asym(PosType x,PosType theta, PosType gamma[]);
  virtual PosType kappa_asym(PosType x,PosType theta);
  virtual void alphakappagamma_asym(PosType x,PosType theta, PosType alpha[]
                                    ,PosType *kappa,PosType gamma[],PosType *phi);
  virtual void alphakappagamma1asym(PosType x,PosType theta, PosType alpha[2]
                                    ,PosType *kappa,PosType gamma[],PosType *phi);
  virtual void alphakappagamma2asym(PosType x,PosType theta, PosType alpha[2]
                                    ,PosType *kappa,PosType gamma[],PosType *phi);
  virtual void alphakappagamma3asym(PosType x,PosType theta, PosType alpha[2]
                                    ,PosType *kappa,PosType gamma[],PosType *phi);
  
  virtual PosType alpha_ell(PosType x,PosType theta);
  
  double fourier_coeff(double n, double q, double beta);
  double IDAXDM(double lambda, double a2, double b2, double x[], double rmax, double mo);
  double IDAYDM(double lambda, double a2, double b2, double x[], double rmax, double mo);
  double SCHRAMMKN(double n, double x[], double rmax);
  double SCHRAMMJN(double n, double x[], double rmax);
  double SCHRAMMI(double x[], double rmax);

  
  void calcModes(double q, double beta, double rottheta, PosType newmod[]);
  void calcModesB(PosType x, double q, double beta, double rottheta, PosType newmod[]);
  void calcModesC(PosType beta_r, double q, double rottheta, PosType newmod[]);
  
  virtual PosType inline InterpolateModes(int whichmod, PosType q, PosType b){return 0;};
  
  void analModes(int modnumber, PosType my_beta, PosType q, PosType amod[3]);
  
  struct fourier_func{
    fourier_func(double my_n, double my_q, double my_beta): n(my_n),q(my_q),beta(my_beta){};
    double n;
    double q;
    double beta;
    double operator ()(double theta) {return cos(n*theta)/pow(cos(theta)*cos(theta) + 1/q/q*sin(theta)*sin(theta),beta/2) ;}
  };
 
  struct SCHRAMMJ_func{
    SCHRAMMJ_func(double my_n, double my_x[], double my_rmax, LensHalo *my_halo): n(my_n), x(my_x), rmx(my_rmax), halo(my_halo){};
    double n, *x, rmx;
    double operator ()(double u) {
      
      // PosType xisq=sqrt(u*(x[0]*x[0]+x[1]*x[1]/(1-(1-halo->fratio*halo->fratio)*u)));
      PosType xisq=sqrt((x[0]*x[0]+x[1]*x[1]/(1-u*(1-halo->fratio*halo->fratio))));
      
      if(xisq*xisq < 1e-20) xisq = 1.0e-10;
      KappaType kappa=halo->kappa_h(xisq)/PI/xisq/xisq;
      //std::cout << "Schramm: n=" << n << " " << u << " " << xisq << " " << halo->fratio << " " << rmx << " " << halo->Rmax << " " << halo->rscale << std::endl;
      return kappa*halo->mass/pow((1.-(1.-halo->fratio*halo->fratio)*u),n+0.5);
    }
  private:
    LensHalo *halo;
  };
  
  struct SCHRAMMK_func{
    SCHRAMMK_func(double my_n, double my_x[], double my_rmax, LensHalo *my_halo): n(my_n), x(my_x), rmx(my_rmax), halo(my_halo){};
    double n, *x, rmx;
    double operator ()(double u) {
      
     // PosType xisq=sqrt(u*(x[0]*x[0]+x[1]*x[1]/(1-(1-halo->fratio*halo->fratio)*u)));
      PosType xisq=sqrt((x[0]*x[0]+x[1]*x[1]/(1-u*(1-halo->fratio*halo->fratio))));
      if(xisq*xisq < 1e-20) xisq = 1.0e-10;
      PosType h=0.0001;
      PosType kp1=halo->kappa_h(xisq+h)/((xisq+h)*(xisq+h));
      PosType km1=halo->kappa_h(xisq-h)/((xisq-h)*(xisq-h));
      PosType kp2=halo->kappa_h(xisq+2*h)/((xisq+2*h)*(xisq+2*h));
      PosType km2=halo->kappa_h(xisq-2*h)/((xisq-2*h)*(xisq-2*h));
      PosType dkdxi=(8*kp1-8*km1-kp2+km2)/12/h/PI;
      


      //std::cout << "Schramm: n=" << n << " " << u << " " << xisq << " " << halo->fratio << " " << rmx << " " << halo->Rmax << " " << halo->rscale << std::endl;
      return u*dkdxi/pow((1.-(1.-halo->fratio*halo->fratio)*u),n+0.5);
      }
    private:
      LensHalo *halo;
    };

  
  struct SCHRAMMI_func{
    SCHRAMMI_func(double my_x[], double my_rmax, LensHalo *my_halo): x(my_x), rmx(my_rmax), halo(my_halo){};
    double *x, rmx;
    double operator ()(double u) {
      
      // PosType xisq=sqrt(u*(x[0]*x[0]+x[1]*x[1]/(1-(1-halo->fratio*halo->fratio)*u)));
      PosType xisq=sqrt((x[0]*x[0]+x[1]*x[1]/(1-u*(1-halo->fratio*halo->fratio))));
      
      if(xisq*xisq < 1e-20) xisq = 1.0e-10;
      PosType alpha=halo->alpha_h(xisq)/PI/xisq;
      //std::cout << "Schramm: n=" << n << " " << u << " " << xisq << " " << halo->fratio << " " << rmx << " " << halo->Rmax << " " << halo->rscale << std::endl;
      return xisq*alpha*halo->mass/pow((1.-(1.-halo->fratio*halo->fratio)*u),0.5);
    }
  private:
    LensHalo *halo;
  };

  
  struct IDAXDM_func{
    IDAXDM_func(double my_lambda,double my_a2, double my_b2, double my_x[], double my_rmax, LensHalo *my_halo): lambda(my_lambda), a2(my_a2), b2(my_b2), x(my_x), rmx(my_rmax), halo(my_halo){};
    double lambda, a2,b2, *x, rmx;
    double operator ()(double m) {
      double ap = m*m*a2 + lambda,bp = m*m*b2 + lambda;
      double p2 = x[0]*x[0]/ap/ap/ap/ap + x[1]*x[1]/bp/bp/bp/bp;  // actually the inverse of equation (5) in Schramm 1990
      PosType tmp = m*rmx;
      if(tmp*tmp < 1e-20) tmp = 1.0e-10;
      PosType xiso=tmp/halo->rscale;
      KappaType kappa=halo->kappa_h(xiso)/PI/xiso/xiso;
      
      PosType alpha[2]={0,0},tm[2] = {m*rmx,0};
      KappaType kappam=0,gamma[2]={0,0},phi=0;
      halo->force_halo_sym(alpha,&kappam,gamma,&phi,tm);
      
      
      //std::cout << "output x: " << m << " " << xiso << " " << m*kappa/(ap*ap*ap*bp*p2)*halo->mass << " " << m*kappam/(ap*ap*ap*bp*p2)<< std::endl;
      double integrandA=m*kappa/(ap*ap*ap*bp*p2)*halo->mass;
      double integrandB=m*kappam/(ap*ap*ap*bp*p2);
      //std::cout << integrandA-integrandB << std::endl;
      assert( std::abs((integrandA - integrandB)/integrandA)<1E-5);

      assert(kappa >= 0.0);

      //return m*kappa/(ap*ap*ap*bp*p2)*halo->mass; // integrand of equation (28) in Schramm 1990
      return integrandB;
    }
  private:
    LensHalo *halo;
  };
  
  struct IDAYDM_func{
    IDAYDM_func(double my_lambda,double my_a2, double my_b2, double my_x[], double my_rmax, LensHalo *my_halo): lambda(my_lambda), a2(my_a2), b2(my_b2), x(my_x), rmx(my_rmax), halo(my_halo){};
    double lambda, a2,b2, *x, rmx;
    double operator ()(double m) {
      double ap = m*m*a2 + lambda,bp = m*m*b2 + lambda;
      double p2 = x[0]*x[0]/ap/ap/ap/ap + x[1]*x[1]/bp/bp/bp/bp;  // actually the inverse of equation (5) in Schramm 1990
      PosType tmp = m*rmx;
      if(tmp*tmp < 1e-20) tmp = 1.0e-10;
      PosType xiso=tmp/halo->rscale;
      KappaType kappa=halo->kappa_h(xiso)/PI/xiso/xiso;
      
      PosType alpha[2]={0,0},tm[2] = {m*rmx,0};
      KappaType kappam=0,gamma[2]={0,0},phi=0;
      halo->force_halo_sym(alpha,&kappam,gamma,&phi,tm);
      
      double integrandA=m*kappa/(ap*ap*ap*bp*p2)*halo->mass;
      double integrandB=m*kappam/(ap*bp*bp*bp*p2);
      assert( std::abs((integrandA - integrandB)/integrandA)<1E-5);

      assert(kappa >= 0.0);
      //return m*kappa/(ap*bp*bp*bp*p2)*halo->mass; // integrand of equation (28) in Schramm 1990
      return integrandB;
    }
  private:
    LensHalo *halo;
  };

  
  const static int Nmod = 32;
  
  // Analytic description of Fourier modes
  
  PosType mod[Nmod];
  PosType mod1[Nmod];
  PosType mod2[Nmod];
  PosType r_eps;
  
  
  // These are stucts used in doing tests
  
  /*struct test_gt_func{
    test_gt_func(LensHalo& halo,PosType my_r): halo(halo),r(my_r){};
    LensHalo& halo;
    PosType r;
    PosType a[2] = {0,0},x[2] = {0,0};
    KappaType k = 0,g[3] = {0,0,0} ,p=0;
    //double operator ()(PosType t) {x[0]=r*cos(t); x[1]=r*sin(t);  halo.force_halo(a,&k,g,&p,x); if(r>1){std::cout << r << " " << g[0]*cos(2*t)+g[1]*sin(2*t) << std::endl;} return (g[0]*cos(2*t)+g[1]*sin(2*t));}
    double operator ()(PosType t) {x[0]=r*cos(t); x[1]=r*sin(t);  halo.force_halo(a,&k,g,&p,x); return (g[0]*cos(2*t)+g[1]*sin(2*t));}
    
  };*/
  
  
  struct test_gt_func{
    test_gt_func(PosType my_r,LensHalo *halo): r(my_r),halo(halo){};
    double operator ()(PosType t) {
      PosType alpha[2] = {0,0},x[2] = {0,0};
      KappaType kappa = 0,gamma[3] = {0,0,0},phi =0 ;
      x[0]=r*cos(t);
      x[1]=r*sin(t);
      halo->force_halo(alpha,&kappa,gamma,&phi,x);
      assert(gamma[0]==gamma[0]);
      assert(gamma[1]==gamma[1]);
       return (gamma[0]*cos(2*t)+gamma[1]*sin(2*t));}
    //double operator ()(PosType t) {x[0]=r*cos(t); x[1]=r*sin(t);  halo.force_halo(a,&k,g,&p,x); if(r>1){std::cout << r << " " << g[0]*cos(2*t)+g[1]*sin(2*t) << std::endl;} return (g[0]*cos(2*t)+g[1]*sin(2*t));}
  private:
    PosType r;
    LensHalo *halo;
  };

  
 /* struct test_kappa_func{
    test_kappa_func(LensHalo& halo,PosType my_r): halo(halo),r(my_r){};
    LensHalo& halo;
    PosType r;
    PosType a[2] = {0,0},x[2] = {0,0};
    KappaType k = 0,g[3] = {0,0,0} ,p=0;
    double operator ()(PosType t) {x[0]=r*cos(t); x[1]=r*sin(t);  halo.force_halo(a,&k,g,&p,x);return 2*PI*k*r*r; }
  };
*/
  
  struct test_kappa_func{
    test_kappa_func(PosType my_r,LensHalo *halo): r(my_r),halo(halo){};
    double operator ()(PosType t) {
      PosType alpha[2] = {0,0},x[2] = {0,0};
      KappaType kappa = 0,gamma[3] = {0,0,0},phi =0 ;
      x[0]=r*cos(t);
      x[1]=r*sin(t);
      halo->force_halo(alpha,&kappa,gamma,&phi,x);
      assert(kappa==kappa);
      return kappa; }
  private:
    PosType r;
    LensHalo *halo;
  };
  
  struct DMDTHETA{
    DMDTHETA(PosType R,LensHalo *halo): R(R),halo(halo){};
    PosType operator()(PosType theta){
      PosType alpha[2] = {0,0},x[2] = {0,0};
      KappaType kappa = 0,gamma[3] = {0,0,0},phi =0 ;
      
      x[0] = R*cos(theta);
      x[1] = R*sin(theta);
      
      halo->force_halo(alpha,&kappa,gamma,&phi,x);
      
      PosType alpha_r = -alpha[0]*cos(theta) - alpha[1]*sin(theta);
      //std::cout << alpha[0] << "  " << alpha[1] << std::endl;
      assert( alpha_r == alpha_r );
      //std::cout << theta << "  " << alpha_r << std::endl;
      return alpha_r;
    }
  private:
    PosType R;
    LensHalo *halo;
  };
  
  struct DMDRDTHETA{
    DMDRDTHETA(PosType R,LensHalo *halo): R(R),halo(halo){};
    PosType operator()(PosType theta){
      
      PosType alpha[2] = {0,0},x[2] = {0,0};
      KappaType kappa = 0,gamma[3] = {0,0,0} ,phi=0;
      
      x[0] = R*cos(theta);
      x[1] = R*sin(theta);
      
      halo->force_halo(alpha,&kappa,gamma,&phi,x);
      assert(kappa == kappa);
      assert(kappa != INFINITY);
      return kappa;
    }
  protected:
    PosType R;
    LensHalo *halo;
  };
  
  struct DMDR{
    DMDR(LensHalo *halo): halo(halo){};
    PosType operator()(PosType logR){
       if(halo->get_flag_elliptical()){
        LensHalo::DMDRDTHETA dmdrdtheta(exp(logR),halo);
        //std::cout << " R = " << exp(logR) << std::endl;
        
        if(exp(2*logR) == 0.0) return 0.0;
         return Utilities::nintegrate<LensHalo::DMDRDTHETA,PosType>(dmdrdtheta,0,2*PI,1.0e-7)
         *exp(2*logR);
      }else{
        PosType alpha[2] = {0,0},x[2] = {0,0};
        KappaType kappa = 0,gamma[3] = {0,0,0} ,phi=0;
        
        x[0] = exp(logR);
        x[1] = 0;
        
        halo->force_halo(alpha,&kappa,gamma,&phi,x);
        return 2*PI*kappa*exp(2*logR);
      }
    }
  protected:
    LensHalo *halo;
    
  };

  struct DALPHAXDM{
    DALPHAXDM(PosType lambda,PosType a2,PosType b2,PosType X[],LensHalo* myisohalo)
    :lambda(lambda),a2(a2),b2(b2),x(X),isohalo(myisohalo){};
    
    PosType operator()(PosType m);
  private:
    PosType lambda;
    PosType a2;
    PosType b2;
    PosType *x;
    LensHalo* isohalo;
  };
  struct DALPHAYDM{
    DALPHAYDM(PosType lambda,PosType a2,PosType b2,PosType X[],LensHalo* isohalo)
    :lambda(lambda),a2(a2),b2(b2),x(X),isohalo(isohalo){};
    
    PosType operator()(PosType m);
  private:
    PosType lambda;
    PosType a2;
    PosType b2;
    PosType *x;
    LensHalo* isohalo;
  };
};

/**
 * \brief A class for calculating the deflection, kappa and gamma caused by an NFW
 * halos.
 *
 * This class uses the true expressions for the NFW profile.  This is
 * time consuming and not usually necessary. See TreeQuadPseudoNFW for a faster alternative.
 *
 * The default value of theta = 0.1 generally gives better than 1% accuracy on alpha.
 * The shear and kappa is always more accurate than the deflection.
 *
 */
class LensHaloNFW: public LensHalo{
public:
  /// Shell constructor.  Sets the halo to zero mass
	LensHaloNFW();
  LensHaloNFW(float my_mass   /// in solar masses
              ,float my_Rsize  /// in Mpc
              ,PosType my_zlens   /// redshift
              ,float my_concentration
              ,float my_fratio    /// axis ratio
              ,float my_pa        /// position angle, it is off by PI/2 and orientation from some others
              ,const COSMOLOGY &cosmo
              ,EllipMethod my_ellip_method=EllipMethod::Pseudo
              );
	//LensHaloNFW(InputParams& params);
  LensHaloNFW(const LensHaloNFW &h):LensHalo(h){
    ++count;
    gmax = h.gmax;
  }
  LensHaloNFW(const LensHaloNFW &&h):LensHalo(std::move(h)){
    ++count;
    gmax = h.gmax;
  }

  LensHaloNFW & operator=(const LensHaloNFW &h){
    if(this != &h){
      LensHalo::operator=(h);
      gmax = h.gmax;
    }
    return *this;
  }
  LensHaloNFW & operator=(const LensHaloNFW &&h){
    if(this != &h){
      LensHalo::operator=(std::move(h));
      gmax = h.gmax;
    }
    return *this;
  }
  
	virtual ~LensHaloNFW();
  
	PosType ffunction(PosType x) const;
	PosType gfunction(PosType x) const;
  PosType dgfunctiondx(PosType x);
	PosType g2function(PosType x) const;
	PosType hfunction(PosType x) const;
  PosType dhfunction(PosType x) const;
  PosType ddhfunction(PosType x, bool numerical);
  PosType dddhfunction(PosType x, bool numerical);
  PosType bfunction(PosType x);
  PosType dbfunction(PosType x);
  PosType ddbfunction(PosType x);
  PosType dmoddb(int whichmod, PosType q, PosType b);
  PosType ddmoddb(int whichmod, PosType q, PosType b);
  PosType dmoddq(int whichmod, PosType q, PosType b);
  PosType ddmoddq(int whichmod, PosType q, PosType b);
  
  //PosType dmod(PosType x, int modnumber, PosType my_slope, PosType my_fratio);    // was used for Ansatz III w derivatives of the Fourier modes
  //PosType ddmod(PosType x, int modnumber, PosType my_slope, PosType my_fratio);   // was used for Ansatz III w derivatives of the Fourier modes
  
  
	// TODO: BEN: the below functions alphaNFW, kappaNFW and gammaNFW are obsolete and better to be deleted to avoid confusion
	void alphaNFW(PosType *alpha,PosType *x,PosType Rtrunc,PosType mass,PosType r_scale
                ,PosType *center,PosType Sigma_crit);
	KappaType kappaNFW(PosType *x,PosType Rtrunc,PosType mass,PosType r_scale
                     ,PosType *center,PosType Sigma_crit);
	void gammaNFW(KappaType *gamma,PosType *x,PosType Rtrunc,PosType mass,PosType r_scale
                ,PosType *center,PosType Sigma_crit);
  
	void initFromFile(float my_mass, long *seed, float vmax, float r_halfmass);
	void initFromMassFunc(float my_mass, float my_Rsize, float my_rscale, PosType my_slope, long *seed);
  /// set Rsize, xmax and gmax
  void set_RsizeXmax(float my_Rsize){LensHalo::setRsize(my_Rsize); xmax = LensHalo::getRsize()/rscale; gmax = InterpolateFromTable(gtable,xmax);};
  /// set scale radius
  /// set rscale, xmax and gmax
	void set_rscaleXmax(float my_rscale){rscale=my_rscale; xmax = LensHalo::getRsize()/rscale; gmax = InterpolateFromTable(gtable,xmax);};
  
protected:
	/// table size
	static const long NTABLE;
	/// maximum Rsize/rscale
	static const PosType maxrm;
	/// keeps track of how many time the tables are created, default is just once
	static int count;
  
	/// tables for lensing properties specific functions
	static PosType *ftable,*gtable,*g2table,*htable,*xtable,*xgtable,***modtable; // modtable was used for Ansatz IV and worked well
  /// make the specific tables
	void make_tables();
	/// interpolates from the specific tables
	PosType InterpolateFromTable(PosType *table, PosType y) const;
  PosType InterpolateModes(int whichmod, PosType q, PosType b);
  
	/// read in parameters from a parameterfile in InputParams params
	void assignParams(InputParams& params);
  
	// Override internal structure of halos
  /// r |alpha(r = x*rscale)| PI Sigma_crit / Mass
	inline PosType alpha_h(PosType x) const {
		//return -1.0*InterpolateFromTable(gtable,x)/InterpolateFromTable(gtable,xmax);
		return -1.0*InterpolateFromTable(gtable,x)/gmax;
    // return -0.5/x*InterpolateFromTable(gtable,x)/gmax;
	}
	inline KappaType kappa_h(PosType x) const{
		return 0.5*x*x*InterpolateFromTable(ftable,x)/gmax;
	}
	inline KappaType gamma_h(PosType x) const{
		//return -0.25*x*x*InterpolateFromTable(g2table,x)/gmax;
		return -0.5*x*x*InterpolateFromTable(g2table,x)/gmax;
	}
	inline KappaType phi_h(PosType x) const{
    return 0.25*(InterpolateFromTable(htable,x) - InterpolateFromTable(htable,LensHalo::getRsize()/rscale))/gmax + log(LensHalo::getRsize()) ;
    // The constant contribution is made to match with the point mass at x = Rsize/rscale.
	}
  inline KappaType phi_int(PosType x) const{
    return -1.0*InterpolateFromTable(xgtable,x)/gmax; //alpha_int(x);
  }
  
private:
  PosType gmax;
};
// ********************


/**
 *
 * \brief A class for calculating the deflection, kappa and gamma caused by a collection of
 * halos with a double power-law mass profile.
 *
 * Spherical halos with \f$ \Sigma \propto 1/(1 + r/r_s )^\beta \f$ so beta would usually be positive.
 *
 * An NFW profile is approximated beta = 2 .
 *
 * The default value of theta = 0.1 generally gives better than 1% accuracy on alpha.
 * The shear and kappa is always more accurate than the deflection.
 */
class LensHaloPseudoNFW: public LensHalo{
public:
  /// shell constructor, should be avoided
	LensHaloPseudoNFW();
  LensHaloPseudoNFW(float my_mass,float my_Rsize,PosType my_zlens,float my_rscale,PosType my_beta,float my_fratio,float my_pa,const COSMOLOGY &cosmo, EllipMethod my_ellip_method=EllipMethod::Pseudo);
	//LensHaloPseudoNFW(InputParams& params);
	~LensHaloPseudoNFW();
  
	PosType mhat(PosType y, PosType beta) const;
  PosType gfunction(PosType x) const;
  
	/// set the slope of the surface density profile
	void set_slope(PosType my_slope){beta=my_slope; make_tables();};
	/// initialize from a mass function
  PosType get_slope(){return beta;};
  /// set Rsize
  //void set_Rsize(float my_Rsize){Rsize = my_Rsize; xmax = Rsize/rscale;};
	
	void initFromMassFunc(float my_mass, float my_Rsize, float my_rscale, PosType my_slope, long *seed);
  
private:
	/// table size
	static const long NTABLE;
	/// maximum Rsize/rscale
	static const PosType maxrm;
	/// keeps track of how many time the tables are created, default is just once
	static int count;
  
	/// tables for lensing properties specific functions
	static PosType *mhattable,*xtable;
	/// make the specific tables
	void make_tables();
	/// interpolates from the specific tables
	PosType InterpolateFromTable(PosType y) const;
  
	/// read in parameters from a parameterfile in InputParams params
	void assignParams(InputParams& params);
  
	/// slope of the surface density profile
	PosType beta;
  
	// Override internal structure of halos
  /// r |alpha(r)| PI Sigma_crit / Mass
	inline PosType alpha_h(PosType x) const {
		return -1.0*InterpolateFromTable(x)/InterpolateFromTable(xmax);
	}
	inline KappaType kappa_h(PosType x) const {
		return 0.5*x*x/InterpolateFromTable(xmax)/pow(1+x,beta);
	}
	inline KappaType gamma_h(PosType x) const {
		//return (0.5*x*x/pow(1+x,beta) - InterpolateFromTable(x))/InterpolateFromTable(xmax);
		return (x*x/pow(1+x,beta) - 2*InterpolateFromTable(x))/InterpolateFromTable(xmax);
	}
	inline KappaType phi_h(PosType x) const {
		return -1.0*alpha_int(x)/InterpolateFromTable(xmax) ;
    //ERROR_MESSAGE();
		//std::cout << "time delay has not been fixed for PseudoNFW profile yet." << std::endl;
		//exit(1);
		//return 0.0;
	}
    inline KappaType phi_int(PosType x) const{
        return -1.0*alpha_int(x)/InterpolateFromTable(xmax) ;
    }
};


/** \brief A class for calculating the deflection, kappa and gamma caused by a collection of halos
 * with truncated power-law mass profiles.
 *
 *The truncation is in 2d not 3d. \f$ \Sigma \propto r^{-\beta} \f$ so beta would usually be positive.
 *
 * The default value of theta = 0.1 generally gives better than 1% accuracy on alpha.
 * The shear and kappa is always more accurate than the deflection.
 */
class LensHaloPowerLaw: public LensHalo{
public:
	LensHaloPowerLaw();
  LensHaloPowerLaw(float my_mass,float my_Rsize,PosType my_zlens,PosType my_beta
                   ,float my_fratio,float my_pa,const COSMOLOGY &cosmo
                   ,EllipMethod my_ellip_method=EllipMethod::Pseudo);
	//LensHaloPowerLaw(InputParams& params);
	~LensHaloPowerLaw();
  
	/// set the slope of the surface density profile
	//void set_slope(PosType my_slope){beta=my_slope;};
  
  /// get slope
  //PosType get_slope(){return beta;};
  
	/// initialize from a mass function
	void initFromMassFunc(float my_mass, float my_Rsize, float my_rscale, PosType my_slope, long *seed);
  
private:
	/// read-in parameters from the parameter file
	void assignParams(InputParams& params);
  
	///	read in parameters from a parameterfile in InputParams params
	// PosType beta;
  
	// Override internal structure of halos
  /// r |alpha(r)| PI Sigma_crit / Mass
	inline PosType alpha_h(
                         PosType x  /// r/rscale
                         ) const{
		if(x==0) x=1e-6*xmax;
		return -1.0*pow(x/xmax,-beta+2);
	}
  /// this is kappa Sigma_crit PI (r/rscale)^2 / mass
	inline KappaType kappa_h(
                           PosType x   /// r/rscale
                           ) const {
		if(x==0) x=1e-6*xmax;
    double tmp = x/xmax;
		return 0.5*(-beta+2)*pow(tmp,2-beta);
	}
  /// this is |gamma| Sigma_crit PI (r/rscale)^2 / mass
	inline KappaType gamma_h(PosType x) const {
		if(x==0) x=1e-6*xmax;
    return -beta*pow(x/xmax,-beta+2);
	}
  /// this is phi Sigma_crit PI / mass, the constants are added so that it is continous over the bourder at Rsize
 	inline KappaType phi_h(PosType x) const {
		if(x==0) x=1e-6*xmax;
    return ( pow(x/xmax,2-beta) - 1 )/(2-beta) + log(LensHalo::getRsize()) ;
	}
  inline KappaType phi_int(PosType x) const{
		//return alpha_int(x);
    return -1.0*pow(x/xmax,-beta+2)/(2-beta);
  }
};


/** \brief Represents a non-singular isothermal elliptical lens
 *
 * This is a true NSIE lens rather than an expansion that approximates one.
 *
 * The maximum radius is set by requireing the total mass to match the input mass.
 *  At radii larger than this radius the halo is treated as a point mass.  In the case of
 *  an ellipitcal halo there is a transition region between these to match the solutions
 *  without a discontinuity.  In this region the lensing quanties have small corrections
 *  which do not correspond to a realistic mass distribution.  If the halo is expected to
 *  be sampled at or beyond the maxiumum radius you should consider using a Truncated
 *  NonSingular Isotherma Ellipsoid (`LensHaloTNSIE') which more naturally deals with the
 *   finite truncation.
*/
class LensHaloRealNSIE : public LensHalo{
public:

  /// explicit constructor, Warning: If my_rcore > 0.0 and my_fratio < 1 then the mass will be somewhat less than my_mass.
  LensHaloRealNSIE(float my_mass,PosType my_zlens,float my_sigma
                   ,float my_rcore,float my_fratio,float my_pa,const COSMOLOGY &cosmo);
  
  // Warning: If my_rcore > 0.0 and my_fratio < 1 then the mass will be somewhat less than my_mass.
	//LensHaloRealNSIE(InputParams& params);
  
  LensHaloRealNSIE(const LensHaloRealNSIE &h):
  LensHalo(h)
  {
    ++objectCount;
    sigma = h.sigma;
    fratio = h.fratio;
    pa = h.pa;
    rcore = h.rcore;
    units = h.units;
  }
  
  LensHaloRealNSIE &operator=(const LensHaloRealNSIE &h){
    if(&h == this) return *this;
    LensHalo::operator=(h);
    sigma = h.sigma;
    fratio = h.fratio;
    pa = h.pa;
    rcore = h.rcore;
    units = h.units;
    return *this;
  }
  
	~LensHaloRealNSIE();
  
	/// overridden function to calculate the lensing properties
	void force_halo(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType const *xcm,bool subtract_point=false,PosType screening = 1.0);
	
	/// get the velocity dispersion
	float get_sigma(){return sigma;};
	// get the NSIE radius
	//float get_Rsize(){return Rsize;};
	/// get the axis ratio
	float get_fratio(){return fratio;};
	/// get the position angle
	float get_pa(){return pa;};
	/// get the core radius
	float get_rcore(){return rcore;};
  
	/// set the velocity dispersion
	void set_sigma(float my_sigma){sigma=my_sigma;};
	// set the NSIE radius
	//void set_Rsize(float my_Rsize){Rsize=my_Rsize;};
	///set the axis ratio
	void set_fratio(float my_fratio){fratio=my_fratio;};
	/// set the position angle
	void set_pa(float my_pa){pa=my_pa;};
	/// set the core radius Einstein radius
	void set_rcore(float my_rcore){rcore=my_rcore;};
  
  void setZlens(PosType my_zlens,const COSMOLOGY &cosmo){
    LensHalo::setZlens(my_zlens,cosmo);
  }

  
protected:
  
  float units;
  
  static size_t objectCount;
  static std::vector<double> q_table;
  static std::vector<double> Fofq_table;
  
	/// initialize from a simulation file
	//void initFromFile(float my_mass, long *seed, float vmax, float r_halfmass);
	/// initialize from a mass function
	//void initFromMassFunc(float my_mass, float my_Rsize, float my_rscale, PosType my_slope, long *seed);
	/// simple initialize from mass while setting a random position angle and ellipticity
	//void initFromMass(float my_mass, long *seed);
  
	/// read-in parameters from a parameter file
	void assignParams(InputParams& params);
  
	/// velocity dispersion of NSIE
	float sigma;
	/// axis ratio of surface mass distribution
	float fratio;
	/// position angle on sky, radians
	float pa;
	/// core size of NSIE
	float rcore;
  
  /// for the set fratio, sigma and rcore calculate the radius that contains the correct mass
  PosType rmax_calc();
  void construct_ellip_tables();
  
  struct NormFuncer{
    NormFuncer(double my_q):q(my_q){};
    
    double operator()(double t){
      return 1.0/sqrt( 1 + (q*q-1)*sin(t)*sin(t));
    };
    
  private:
    double q;
  };
  
};

/** \brief Truncated non-singular isothermal ellipsoid

This is a true TNSIE lens rather than an expansion that approximates one.
*/
class LensHaloTNSIE : public LensHalo{
public:

  LensHaloTNSIE(float my_mass  /// total mass in Msun
                ,PosType my_zlens /// redshift
                ,float my_sigma  /// vleocity dispertion in km/s
                ,float my_rcore  ///  core size Mpc
                ,float my_fratio /// axis ratio
                ,float my_pa     /// position angle
                ,const COSMOLOGY &cosmo  /// cosmology
                ,float f=20 /// cuttoff radius in units of truncation radius
                );
  
  LensHaloTNSIE(const LensHaloTNSIE &h):
  LensHalo(h)
  {
    sigma = h.sigma;
    fratio = h.fratio;
    pa = h.pa;
    rcore = h.rcore;
    rtrunc = h.rtrunc;
    units = h.units;
  }
  
  LensHaloTNSIE &operator=(const LensHaloTNSIE &h){
    if(&h == this) return *this;
    LensHalo::operator=(h);
    sigma = h.sigma;
    fratio = h.fratio;
    pa = h.pa;
    rcore = h.rcore;
    rtrunc = h.rtrunc;
    units = h.units;
    return *this;
  }
  
  ~LensHaloTNSIE(){};
  
  /// returns the sigma in km / s for the other parameters fixed
  static double calc_sigma(float mass,float Rtrunc,float Rcore,float fratio){
    return lightspeed*sqrt( Grav*mass*sqrt(fratio) / (Rtrunc-Rcore) / PI);
  }

  /// overridden function to calculate the lensing properties
  void force_halo(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType const *xcm,bool subtract_point=false,PosType screening = 1.0);
  
  /// get the velocity dispersion
  float get_sigma(){return sigma;};
  // get the NSIE radius
  //float get_Rsize(){return Rsize;};
  /// get the axis ratio
  float get_fratio(){return fratio;};
  /// get the position angle
  float get_pa(){return pa;};
  /// get the core radius
  float get_rcore(){return rcore;};
  /// get the truncation radius
  float get_rtrunc(){return rtrunc;};
  
  void setZlens(PosType my_zlens,const COSMOLOGY &cosmo){
    LensHalo::setZlens(my_zlens,cosmo);
  }

  /// set the position angle
  void set_pa(float my_pa){pa=my_pa;};
  
protected:
  
  float units;
  
  /// read-in parameters from a parameter file
  void assignParams(InputParams& params);
  /// velocity dispersion of TNSIE
  float sigma;
  /// axis ratio of surface mass distribution
  float fratio;
  /// position angle on sky, radians
  float pa;
  /// core size of NSIE
  float rcore;
  /// core size of NSIE
  float rtrunc;
};

/** \brief A truncated elliptical power-law profile

 This is an implementation of O’Riordan, Warren, & Mortlock (2020)
*/

class LensHaloTEPL : public LensHalo{
public:

  LensHaloTEPL(float my_mass  /// total mass in Msun
                ,PosType my_zlens /// redshift
                ,PosType r_trunc  /// elliptical truncation radius in Mpc
                ,PosType gamma    /// power-law index
                ,float my_fratio /// axis ratio
                ,float my_pa     /// position angle, 0 has long axis along the veritical axis and goes clockwise
                ,const COSMOLOGY &cosmo  /// cosmology
                ,float f=10 /// cuttoff radius in units of truncation radius
  );
  
  LensHaloTEPL(const LensHaloTEPL &h):
  LensHalo(h)
  {
    tt = h.tt;
    x_T = h.x_T;
    q = h.q;
    q_prime = h.q_prime;
    SigmaT = h.SigmaT;
    pa = h.pa;
    mass_pi = h.mass_pi;
    R = h.R;
  }
  
  LensHaloTEPL &operator=(const LensHaloTEPL &h){
    if(&h == this) return *this;
    LensHalo::operator=(h);

    tt = h.tt;
    x_T = h.x_T;
    q = h.q;
    q_prime = h.q_prime;
    SigmaT = h.SigmaT;
    pa = h.pa;
    mass_pi = h.mass_pi;
    R = h.R;

    return *this;
  }
  
  ~LensHaloTEPL(){};


  /// overridden function to calculate the lensing properties
  void force_halo(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType const *xcm,bool subtract_point=false,PosType screening = 1.0);
  
  /// get the axis ratio
  float get_fratio(){return q;};
  /// get the position angle
  float get_pa(){return pa;};
  /// get the truncation radius
  float get_rtrunc(){return x_T;};
  /// get the pwer-law index
  float get_t(){return tt;}
  
  void set_pa(double p){pa = p;}
  
  void deflection(std::complex<double> &z
                  ,std::complex<double> &a
                  ,std::complex<double> &g
                  ,KappaType &sigma) const;
  
  
protected:
  
  //float units;
  double tt;
  double x_T;       // truncation radius
  double q;
  double q_prime;
  double SigmaT;   // SigmaT - surface density at trunction radius
  double pa;
  double mass_pi;
  std::complex<double> R; // rotation
  
  std::complex<double> F(double r_e,double t,std::complex<double> z) const;
};

class LensHaloTEBPL : public LensHalo{
public:
  
  LensHaloTEBPL(float my_mass  /// total mass in Msun
                ,PosType my_zlens /// redshift
                ,PosType r_break  /// elliptical truncation radius in Mpc
                ,PosType r_trunc  /// elliptical truncation radius in Mpc
                ,PosType t1    /// inner power-law index
                ,PosType t2    /// outer power-law index
                ,float my_fratio /// axis ratio
                ,float my_pa     /// position angle, 0 has long axis along the veritical axis and goes clockwise
                ,const COSMOLOGY &cosmo  /// cosmology
                ,float f=10 /// cuttoff radius in units of truncation radius
  );

  ~LensHaloTEBPL(){}
  
  LensHaloTEBPL(const LensHaloTEBPL &h):
  LensHalo(h),h1(h.h1),h2(h.h2),h3(h.h3)
  {
 
    m2 = h.m2;
    m3 = h.m3;
    m1 = h.m1;
    
    rb = h.rb;
    rt = h.rt;
    q = h.q;

    R = h.R;
  }
  
  LensHaloTEBPL &operator=(const LensHaloTEBPL &h){
    if(&h == this) return *this;
    LensHalo::operator=(h);

    h1 = h.h1;
    h2 = h.h2;
    h3 = h.h3;
    
    m2 = h.m2;
    m3 = h.m3;
    m1 = h.m1;
    
    rb = h.rb;
    rt = h.rt;
    q = h.q;

    R = h.R;

    return *this;
  }
  
  /// get the axis ratio
  float get_fratio(){return q;};
  /// get the position angle
  float get_pa(){return pa;};
  /// get the truncation radius
  float get_rtrunc(){return rt;};
  /// get the break radius
  float get_rbreak(){return rb;};
  /// get inner pwer-law index
  float get_t1(){return h1.get_t();}
  /// get outer pwer-law index
  float get_t2(){return h2.get_t();}

  void set_pa(double p){pa = p;}

  void force_halo(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType const *xcm,bool subtract_point=false,PosType screening = 1.0);
  
private:

  double q;
  double rb;
  double rt;

  double m2 ;
  double m3 ;
  double m1 ;

  LensHaloTEPL h1;
  LensHaloTEPL h2;
  LensHaloTEPL h3;
  
  std::complex<PosType> R;
  
};

#ifdef ENABLE_CERF
/**
 
 This class uses the libcerf library (https://jugit.fz-juelich.de/mlz/libcerf).
 It can be installed with homebreww.
 */
class LensHaloGaussian : public LensHalo{
public:

  LensHaloGaussian(float my_mass  /// total mass in Msun
                ,PosType my_zlens /// redshift
                ,PosType r_scale  /// scale hight along the largest dimension
                ,float my_fratio /// axis ratio
                ,float my_pa     /// position angle, 0 has long axis along the veritical axis and goes clockwise
                ,const COSMOLOGY &cosmo  /// cosmology
                ,float f=100 /// cuttoff radius in units of truncation radius
  );
  
  LensHaloGaussian(const LensHaloGaussian &h):
  LensHalo(h)
  {
    Rhight = h.Rhight;
    q = h.q;
    q_prime = h.q_prime;
    SigmaO = h.SigmaO;
    pa = h.pa;
    R = h.R;
    
    norm = h.norm;
    norm_g = h.norm_g;
    ss = h.ss;
    I = h.I;
    I_sqpi = h. I_sqpi;
  }
  
  LensHaloGaussian &operator=(const LensHaloGaussian &h){
    if(&h == this) return *this;
    LensHalo::operator=(h);
    
    Rhight = h.Rhight;
    q = h.q;
    q_prime = h.q_prime;
    SigmaO = h.SigmaO;
    pa = h.pa;
    R = h.R;
    
    norm = h.norm;
    norm_g = h.norm_g;
    ss = h.ss;
    I = h.I;
    I_sqpi = h. I_sqpi;
    
    return *this;
  }
  
  ~LensHaloGaussian(){};


  /// overridden function to calculate the lensing properties
  void force_halo(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType const *xcm,bool subtract_point=false,PosType screening = 1.0);
  
  /// get the axis ratio
  float get_fratio(){return q;};
  /// get the position angle
  float get_pa(){return pa;};
  /// get the truncation radius
  float get_scalehight(){return Rhight;};

  /// central surface density
  float get_SigmaCentral(){return SigmaO;};

  void set_pa(double p){pa = p;}
  
  void deflection(std::complex<double> &z
                  ,std::complex<double> &a
                  ,std::complex<double> &g
                  ,KappaType &sigma) const;
  
  
protected:
  
  std::complex<double> dwdz(std::complex<double> z) const{
    return 2.0*(I_sqpi - z*wF(z));
  };
  
  std::complex<double> wF(std::complex<double> z) const;
  std::complex<double> my_erfc(std::complex<double> z) const;

  // https://arxiv.org/pdf/1407.0748.pdf
//  std::complex<double> erfc(std::complex<double> z) const{
//    if(z.real() >= 0){
//    std::complex<double> ans = wf(I*z);
//    ans *= exp(-z*z);
    //return ans;
//      return exp(-z*z)*wf(I*z);
//    }else{
//      return 2. - exp(-z*z)*wf(-I*z);
//    }
//  }
  
  double Rhight;  // scale hight in larges dimension
  double q;
  double q_prime;
  double SigmaO;   // SigmaO - central surface density
  double pa;
  std::complex<double> R; // rotation
  
  double norm;
  double norm_g;
  double ss;
  
  //std::complex<double> F(double r_e,double t,std::complex<double> z) const;
  std::complex<double> erfi(std::complex<double> z) const{
    return I*(1. - exp(z*z)*wF(z));
  }
  
//  std::complex<double> rho(std::complex<double> z) const{
//
//    //std::complex<double> y=(q*q*z.real() + I*z.imag())/ss/sqrt(1-q*q) ;
//    //return -I * ( wF( z/ss/sqrt(q_prime) )
//    //                        -  exp( y*y -z*z/ss/ss/q_prime ) * wF(y));
//
//    return wbar(z/ss/ss/sqrt(q_prime),1) - wbar(z/ss/ss/sqrt(q_prime),q);
//  }
//
//  std::complex<double> wbar(std::complex<double> z, double p) const{
//    double x = z.real();
//    double y = z.imag();
//    std::cout << exp(-z*z) << "," <<
//    exp(-x*x*(1-p*p) - y*y*(1./p/p -1)) << "," << -I* wF(p*x + I*y/p) << std::endl;
//
//    double tmp = exp(-x*x*(1-p*p) - y*y*(1./p/p -1));
//    if(tmp==0) return 0;
//
//    return  -(exp(-z*z)*I*tmp)*wF(p*x + I*y/p);
//  }
 
  // Faddeeva function according to Zaghloul (2017)
  // It oesn't seem right near the real axis for small |z|!
  // not used
//  std::complex<double> wf(std::complex<double> z) const{
//    //std::cout << z << std::endl;
//    double y2 = z.imag()*z.imag();
//    double z2 = std::norm(z);
//    if( z2 >= 3.8e4 ){
//      return I_sqpi / z;
//    }
//    if(z2 >= 256){
//      //std::cout << z << "," << (z*z -0.5) << "," << I * PI << std::endl;
//      return z/(z*z -0.5) * I_sqpi;
//    }
//    else if(z2 >= 62 && z2 < 256){
//      return (z*z - 1.) / (z * (z*z-1.5)) * I_sqpi;
//    }
//    else if(z2 < 62 && z2 >= 30 && y2 >= 1.0e-13){
//      return z*(z*z -2.5) / (z*z*(z*z-3.) + 0.75) * I_sqpi;
//    }
//    else if ( (z2 < 62 && z2 >= 30 && y2 < 1.0e-13) ||  (z2 < 30 && z2 >= 2.5 && y2 < 0.072) ){
//      return exp(-z*z) + I * z * (U1[5] + z*z*(U1[4] + z*z*(U1[3] + z*z*(U1[2]+                     z*z*(U1[1] + z*z*(U1[0] + z*z*sqrt(PI) ))))))
//      /(V1[6] + z*z*(V1[5]+ z*z*(V1[4] + z*z*(V1[3] + z*z*(V1[2] + z*z*(V1[1] + z*z*(V1[0]+z*z)))))));
//    }
//
//    return (U2[5]-I*z*(U2[4]-I*z*(U2[3]-I*z*(U2[2]-I*z*(U2[1]-I*z*(U2[0]-I*z*sqrt(PI) ))))))
//      /(V2[6] - I*z*(V2[5] - I*z*(V2[4] - I*z*(V2[3] - I*z*(V2[2]-I*z*(V2[1] - I*z*(V2[0] - I*z)))))));
//  }
  
  std::complex<double> I;
  std::complex<double> I_sqpi;
//  double U1[6] = {1.320522,35.7668,219.031,1540.787,3321.990,36183.31};
//  double V1[7] = {1.841439,61.57037,364.2191,2186.181,9022.228,24322.84,32066.6};
//
//  double U2[6] = {5.9126262,30.180142,93.15558,181.92853,214.38239,122.60793};
//  double V2[7] = {10.479857,53.992907,170.35400,348.70392,457.33448,352.73063,122.60793};
};

#endif


/**
 *
 * \brief A class for calculating the deflection, kappa and gamma caused by a halos
 * with truncated Hernquist mass profiles.
 *
 * The profile is \f$ \rho \propto \left( \frac{r}{r_s} \right)^{-1} \left( 1 + \frac{r}{r_s} \right)^{-3} \f$.
 *
 *
 */

class LensHaloHernquist: public LensHalo{
public:
	//LensHaloHernquist();
  LensHaloHernquist(float my_mass,float my_Rsize,PosType my_zlens,float my_rscale,float my_fratio,float my_pa,const COSMOLOGY &cosmo, EllipMethod my_ellip_method=EllipMethod::Pseudo);
	//LensHaloHernquist(InputParams& params);
	virtual ~LensHaloHernquist();
  
  PosType ffunction(PosType x) const;
	PosType gfunction(PosType x) const;
	PosType hfunction(PosType x) const;
	PosType g2function(PosType x) const;
  
	/* the below functions alphaHern, kappaHern and gammaHern are obsolete and better to be deleted to avoid confusion
   void alphaHern(PosType *alpha,PosType *x,PosType Rtrunc,PosType mass,PosType r_scale
   ,PosType *center,PosType Sigma_crit);
   KappaType kappaHern(PosType *x,PosType Rtrunc,PosType mass,PosType r_scale
   ,PosType *center,PosType Sigma_crit);
   void gammaHern(KappaType *gamma,PosType *x,PosType Rtrunc,PosType mass,PosType r_scale
   ,PosType *center,PosType Sigma_crit);
   */
	//void initFromFile(float my_mass, long *seed, float vmax, float r_halfmass);
  
	/// set scale radius
  void set_rscale(float my_rscale){rscale=my_rscale; xmax = LensHalo::getRsize()/rscale; gmax = InterpolateFromTable(gtable,xmax);};
  // friend struct Ialpha_func;
  
protected:
	/// table size
	static const long NTABLE;
	/// maximum Rsize/rscale
	static const PosType maxrm;
	/// keeps track of how many time the tables are created, default is just once
	static int count;
  
	/// tables for lensing properties specific functions
	static PosType *ftable,*gtable,*g2table,*htable,*xtable,*xgtable;
	/// make the specific tables
	void make_tables();
	/// interpolates from the specific tables
	PosType InterpolateFromTable(PosType *table, PosType y) const;
  
	/// read in parameters from a parameterfile in InputParams params
	void assignParams(InputParams& params);
  
	// Override internal structure of halos
  /// r |alpha(r)| PI Sigma_crit / Mass
	inline PosType alpha_h(PosType x) const {
		return -0.25*InterpolateFromTable(gtable,x)/gmax;
	}
	inline KappaType kappa_h(PosType x) const {
		return 0.25*x*x*InterpolateFromTable(ftable,x)/gmax; // 0.5*
	}
	inline KappaType gamma_h(PosType x) const {
		return -0.5*x*x*InterpolateFromTable(g2table,x)/gmax;
	}
	inline KappaType phi_h(PosType x) const {
		return -0.25*InterpolateFromTable(htable,x)/gmax;
		//return -1.0*InterpolateFromTable(htable,x)/gmax;
	}
    inline KappaType phi_int(PosType x) const{
		return -0.25*InterpolateFromTable(xgtable,x)/gmax;
  }
  
  
private:
  PosType gmax;
};

/** 
 *
 * \brief A class for calculating the deflection, kappa and gamma caused by a collection of halos
 * with truncated Jaffe mass profiles.
 *
 * The profile is \f$ \rho \propto \left( \frac{r}{r_s} \right)^{-2} \left( 1 + \frac{r}{r_s} \right)^{-2} \f$ so beta would usually be negative.
 *
 *
 */

class LensHaloJaffe: public LensHalo{
public:
	//LensHaloJaffe();
  LensHaloJaffe(float my_mass,float my_Rsize,PosType my_zlens,float my_rscale,float my_fratio
                ,float my_pa,const COSMOLOGY &cosmo
                , EllipMethod my_ellip_method=EllipMethod::Pseudo);
	//LensHaloJaffe(InputParams& params);
	virtual ~LensHaloJaffe();
  
	/// set scale radius
  void set_rscale(float my_rscale){rscale=my_rscale; xmax = LensHalo::getRsize()/rscale; gmax = InterpolateFromTable(gtable,xmax);};
  
  PosType ffunction(PosType x) const;
	PosType gfunction(PosType x) const;
	PosType hfunction(PosType x) const;
	PosType g2function(PosType x) const;
  PosType bfunction(PosType x);
  PosType dbfunction(PosType x);
  PosType ddbfunction(PosType x);
  
protected:
  
	/// table size
	static const long NTABLE;
	/// maximum Rsize/rscale
	static const PosType maxrm;
	/// keeps track of how many time the tables are created, default is just once
	static int count;
  
	/// tables for lensing properties specific functions
	static PosType *ftable,*gtable,*g2table,*htable,*xtable,*xgtable;
	/// make the specific tables
	void make_tables();
	/// interpolates from the specific tables
	PosType InterpolateFromTable(PosType *table, PosType y) const;
  
	/// read in parameters from a parameterfile in InputParams params
	void assignParams(InputParams& params);
  
	// Override internal structure of halos
  /// r |alpha(r)| PI Sigma_crit / Mass
	inline PosType alpha_h(PosType x) const{
		return -0.25*InterpolateFromTable(gtable,x)/gmax;
	}
	inline KappaType kappa_h(PosType x) const {
		return 0.125*x*x*InterpolateFromTable(ftable,x)/gmax;
	}
	inline KappaType gamma_h(PosType x) const {
		//return -0.125*x*x*InterpolateFromTable(g2table,x)/gmax;
		return -0.25*x*x*InterpolateFromTable(g2table,x)/gmax;
	}
  ///std::cout << "no analytic expression defined yet for Jaffe profile, take numerical" << std::endl;
	inline KappaType phi_h(PosType x) const {
        return -0.25*InterpolateFromTable(xgtable,x)/gmax;
    }
    inline KappaType phi_int(PosType x) const{
		return -0.25*InterpolateFromTable(xgtable,x)/gmax;
  }
  
private:
  PosType gmax;
  
  // I have temporarily set these functions to 0 to make the code compile, Ben
  //  PosType ffunction(PosType x){throw std::runtime_error("Set to temporary invalid value"); return 0;}
  //	PosType gfunction(PosType x){throw std::runtime_error("Set to temporary invalid value"); return 0;}
  //	PosType hfunction(PosType x){throw std::runtime_error("Set to temporary invalid value"); return 0;}
  //	PosType g2function(PosType x){throw std::runtime_error("Set to temporary invalid value"); return 0;}
  
};




/**
 * \brief This is a lens that does no lensing.  It is useful for testing and for running refinement code on sources.
 */
class LensHaloDummy: public LensHalo{
public:
	LensHaloDummy();
  LensHaloDummy(float my_mass,float my_Rsize,PosType my_zlens,float my_rscale,const COSMOLOGY &cosmo);
	//LensHaloDummy(InputParams& params);
	~LensHaloDummy(){};
	
	/// overridden function to calculate the lensing properties
	void force_halo(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType const *xcm,bool subtract_point=false,PosType screening = 1.0);
  // void force_halo(PosType *alpha,KappaType *kappa,KappaType *gamma,PosType *xcm,bool subtract_point=false,PosType screening = 1.0);
  
	/// initialize from a mass function
	void initFromMassFunc(float my_mass, float my_Rsize, float my_rscale, PosType my_slope, long *seed);
  
	
private:
	/// read-in parameters from a parameter file
	void assignParams(InputParams& params);
	inline PosType alpha_h(PosType x) const {return  0.;}
	inline KappaType kappa_h(PosType x) const {return  0.;}
	inline KappaType gamma_h(PosType x) const {return  0.;}
};


typedef LensHalo* LensHaloHndl;

//bool LensHaloZcompare(LensHalo *lh1,LensHalo *lh2);//{return (lh1->getZlens() < lh2->getZlens());}
//bool LensHaloZcompare(LensHalo lh1,LensHalo lh2){return (lh1.getZlens() < lh2.getZlens());}


#endif /* LENS_HALOS_H_ */
