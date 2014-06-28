/*
 * lens_halos.h
 *
 *  Created on: 06.05.2013
 *      Author: mpetkova
 */

#ifndef LENS_HALOS_H_
#define LENS_HALOS_H_

#include "standard.h"
#include "InputParams.h"
#include "source.h"
#include "point.h"

//#include "quadTree.h"

class TreeQuad;

/**
 * \brief A base class for all types of lensing halos.
 *
 * It contains the mass, maximum radius (Rmax), and scale radius (rscale) of the halo,
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
 * a mass function. Then one needs all parameters of the halo -- mass, Rmax, and rscale.
 */

class LensHalo{
public:
	LensHalo();
	LensHalo(InputParams& params);
	virtual ~LensHalo();

	/// get the Rmax
	float get_Rmax() const { return Rmax; }
	/// get the mass
	float get_mass() const { return mass; }
	/// get the scale radius
	float get_rscale() const { return rscale; }
	/// get the redshift

	PosType getZlens() const { return zlens; }

    /// set the position of the Halo :
    /// taking as argument the position in comoving Mpc
    /// storing the position in physical Mpc
    void setX(PosType PosX, PosType PosY) { posHalo[0] = PosX / (1 + zlens) ; posHalo[1] = PosY / (1 + zlens) ; }
    /// get the position of the Halo in physical Mpc
    void getX(PosType * MyPosHalo) const { MyPosHalo[0] = posHalo[0] ; MyPosHalo[1] = posHalo[1]; }
    /// display the position of the halo
    void displayPos() { std::cout << "Halo PosX = " << posHalo[0] << " ; Halo PosY = " << posHalo[1] << std::endl; }
 
	/// initialize from a simulation file
	virtual void initFromFile(float my_mass, long *seed, float vmax, float r_halfmass){};
	/// initialize from a mass function
	virtual void initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, PosType my_slope, long *seed);

	/// set Rmax (in Mpc)
	virtual void set_Rmax(float my_Rmax){Rmax=my_Rmax; xmax = Rmax/rscale;};
	/// set mass (in solar masses)
	void set_mass(float my_mass){mass=my_mass;};
	/// set scale radius (in Mpc)
	virtual void set_rscale(float my_rscale){rscale=my_rscale; xmax = Rmax/rscale;};
	/// set redshift
	void setZlens(PosType my_zlens){ zlens=my_zlens; };

	/// set slope
	virtual void set_slope(PosType my_slope){beta=my_slope;};
    /// get slope
    virtual PosType get_slope(){return beta;};
    /// flag=True if halo elliptical
    bool get_flag_elliptical(){return elliptical_flag;};
    void set_flag_elliptical(bool ell){elliptical_flag=ell;};
    
    /// set cosmology for halo
	virtual void setCosmology(const COSMOLOGY& cosmo) {}
	
	/// calculate the lensing properties -- deflection, convergence, and shear

	virtual void force_halo(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType const *xcm,bool subtract_point=false);

	/// force tree calculation for stars
	void force_stars(PosType *alpha,KappaType *kappa,KappaType *gamma,PosType const *xcm);

	/// internal compare redshift function
	bool compareZ(PosType z){return z > zlens;};

  /// stars
  bool AreStarsImaplated() const {return stars_implanted;}
  void implant_stars(PosType **centers,int Nregions,long *seed, IMFtype type=One);
  void implant_stars(PosType *center,long *seed,IMFtype type = One);
  //void implant_stars(PosType *x,PosType *y,int Nregions,long *seed,IMFtype type=One);
  void remove_stars();
  IMFtype getStarIMF_type() const {return main_stars_imf_type;}
  /// Fraction of surface density in stars
  PosType getFstars() const {return star_fstars;}
  /// The mass of the stars if they are all the same mass
  PosType getStarMass() const {if(stars_implanted)return star_masses[0]; else return 0.0;}

	/// get the number of halo parameters
	virtual std::size_t Nparams() const;
	/// get the value of a scaled halo parameter by index
	virtual PosType getParam(std::size_t p) const;
	/// set the value of a scaled halo parameter by index
	virtual PosType setParam(std::size_t p, PosType value);
	
	/// print the halo parameters in CSV format
	virtual void printCSV(std::ostream&, bool header = false) const;

	/// Prints star parameters; if show_stars is true, prints data for single stars
	void PrintStars(bool show_stars) const;
 
  
protected:
  PosType alpha_int(PosType x);
  
  void force_halo_sym(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType const *xcm,bool subtract_point=false);
	void force_halo_asym(PosType *alpha,KappaType *kappa,KappaType *gamma,PosType const *xcm,bool subtract_point=false);

  //friend struct Ig_func;
  
  struct Ialpha_func{
    Ialpha_func(LensHalo& halo): halo(halo){};
    LensHalo& halo;
    PosType operator ()(PosType x) {return halo.alpha_h(x)/x ;}
  };
  
  struct Ig_func{
    Ig_func(LensHalo& halo): halo(halo){};
    LensHalo& halo;
    PosType operator ()(PosType x) {return halo.gfunction(x)/x ;}
  };

  IndexType *stars;
  PosType **stars_xp;
  TreeQuad *star_tree;
  int stars_N;
  PosType star_massscale;
  /// star masses relative to star_massscles
  float *star_masses;
  PosType star_fstars;
  PosType star_theta_force;
  int star_Nregions;
  PosType *star_region;
  PosType beta;
  void substract_stars_disks(PosType const *ray,PosType *alpha
                             ,KappaType *kappa,KappaType *gamma);
  float* stellar_mass_function(IMFtype type, unsigned long Nstars, long *seed, PosType minmass=0.0, PosType maxmass=0.0
                               ,PosType bendmass=0.0, PosType powerlo=0.0, PosType powerhi=0.0);

  
	/// read in parameters from a parameterfile in InputParams params
	void assignParams(InputParams& params);
    /// read in star parameters. This is valid for all halos and not overloaded.
    void assignParams_stars(InputParams& params);

	/// error message printout
	void error_message1(std::string name,std::string filename);


    float mass;
    /// Radius of halo and NSIE if it exists,  This is the radius used in the tree force solver
    /// to determine when a ray intersects an object.
    float Rmax;
    /// scale length or core size.  Different meaning in different cases.  Not used in NSIE case.
    float rscale;
    /// redshift
    //PosType zlens;

    bool stars_implanted;
    /// Number of regions to be subtracted to compensate for the mass in stars
    IMFtype main_stars_imf_type;
    PosType main_stars_min_mass;
    PosType main_stars_max_mass;
    PosType bend_mstar;
    PosType lo_mass_slope;
    PosType hi_mass_slope;
    /// parameters for stellar mass function: minimal and maximal stellar mass, bending point for a broken power law IMF
    PosType *star_Sigma;
    PosType **star_xdisk;


  /// point mass case
	virtual PosType inline alpha_h(PosType x){return -1;};
	virtual KappaType inline kappa_h(PosType x){return 0;};
	virtual KappaType inline gamma_h(PosType x){return -2;};
	virtual KappaType inline phi_h(PosType x){return 1;};
    virtual KappaType inline phi_int(PosType x){return 1;};
    
    virtual PosType inline gfunction(PosType x){return -1;};
	
    PosType xmax;

  // Functions for calculating axial dependence
    float pa;
    float fratio=1;
    bool elliptical_flag = false;

    void faxial(PosType theta,PosType f[]);
    void gradial(PosType r,PosType g[]);
    void gradial2(PosType r,PosType g[]);
    
    //void felliptical(PosType x, PosType q, PosType theta, PosType f[], PosType g[]);

	virtual void gamma_asym(PosType x,PosType theta, PosType gamma[2]);
	virtual PosType kappa_asym(PosType x,PosType theta);
	virtual void alpha_asym(PosType x,PosType theta, PosType alpha[2]);
    double fourier_coeff(double n, double q, double beta);
    
    void calcModes(double q, double beta, double rottheta, PosType newmod[]);
    
    struct fourier_func{
        fourier_func(double my_n, double my_q, double my_beta): n(my_n),q(my_q),beta(my_beta){};
        double n;
        double q;
        double beta;
        double operator ()(double theta) {return cos(n*theta)/pow(cos(theta)*cos(theta) + (1/q/q)*sin(theta)*sin(theta),beta/2) ;}
    };


  const static int Nmod = 32;
  
  PosType mod[Nmod];
  PosType r_eps;
  

    PosType zlens;
    
    /// Position of the Halo
    PosType posHalo[2];

};

/** \ingroup DeflectionL2
 *
 * \brief A class for calculating the deflection, kappa and gamma caused by a collection of NFW
 * halos.
 *
 * Derived from the TreeQuad class.  The "particles" are replaced with spherical NFW halos.
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
  /// Shell constructor that should be avoided
	LensHaloNFW();
  LensHaloNFW(float my_mass,float my_Rmax,PosType my_zlens,float my_rscale,float my_fratio,float my_pa,int my_stars_N);
	LensHaloNFW(InputParams& params);
	virtual ~LensHaloNFW();

	PosType ffunction(PosType x);
	PosType gfunction(PosType x);
	PosType g2function(PosType x);
	PosType hfunction(PosType x);
    

	// TODO: BEN: the below functions alphaNFW, kappaNFW and gammaNFW are obsolete and better to be deleted to avoid confusion
	void alphaNFW(PosType *alpha,PosType *x,PosType Rtrunc,PosType mass,PosType r_scale
			,PosType *center,PosType Sigma_crit);
	KappaType kappaNFW(PosType *x,PosType Rtrunc,PosType mass,PosType r_scale
			,PosType *center,PosType Sigma_crit);
	void gammaNFW(KappaType *gamma,PosType *x,PosType Rtrunc,PosType mass,PosType r_scale
			,PosType *center,PosType Sigma_crit);

	void initFromFile(float my_mass, long *seed, float vmax, float r_halfmass);
	void initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, PosType my_slope, long *seed);
  /// set Rmax
    void set_Rmax(float my_Rmax){Rmax=my_Rmax; xmax = Rmax/rscale; gmax = InterpolateFromTable(gtable,xmax);};
  /// set scale radius
	void set_rscale(float my_rscale){rscale=my_rscale; xmax = Rmax/rscale; gmax = InterpolateFromTable(gtable,xmax);};
    

protected:
	/// table size
	static const long NTABLE;
	/// maximum Rmax/rscale
	static const PosType maxrm;
	/// keeps track of how many time the tables are created, default is just once
	static int count;

	/// tables for lensing properties specific functions
	static PosType *ftable,*gtable,*g2table,*htable,*xtable,*xgtable;
	/// make the specific tables
	void make_tables();
	/// interpolates from the specific tables
	PosType InterpolateFromTable(PosType *table, PosType y);

	/// read in parameters from a parameterfile in InputParams params
	void assignParams(InputParams& params);

	/// Override internal structure of halos
	inline PosType alpha_h(PosType x){
		//return -1.0*InterpolateFromTable(gtable,x)/InterpolateFromTable(gtable,xmax);
		return -1.0*InterpolateFromTable(gtable,x)/gmax;
	// return -0.5/x*InterpolateFromTable(gtable,x)/gmax;
	}
	inline KappaType kappa_h(PosType x){
		return 0.5*x*x*InterpolateFromTable(ftable,x)/gmax;
	}
	inline KappaType gamma_h(PosType x){
		return -0.25*x*x*InterpolateFromTable(g2table,x)/gmax;
	}
	inline KappaType phi_h(PosType x){
		//ERROR_MESSAGE();
		//std::cout << "time delay has not been fixed for NFW profile yet." << std::endl;
		//exit(1);
		return -0.25*InterpolateFromTable(htable,x)/gmax;
	}
    inline KappaType phi_int(PosType x){
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
 * Derived from the TreeQuad class.  The "particles" are replaced with spherical halos
 * with \f$ \Sigma \propto 1/(1 + r/r_s )^\beta \f$ so beta would usually be positive.
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
  LensHaloPseudoNFW(float my_mass,float my_Rmax,PosType my_zlens,float my_rscale,PosType my_beta,float my_fratio,float my_pa,int my_stars_N);
	LensHaloPseudoNFW(InputParams& params);
	~LensHaloPseudoNFW();

	PosType mhat(PosType y, PosType beta);
    PosType gfunction(PosType x);
    
	/// set the slope of the surface density profile
	void set_slope(PosType my_slope){beta=my_slope; make_tables();};
	/// initialize from a mass function
    PosType get_slope(){return beta;};
    /// set Rmax
    void set_Rmax(float my_Rmax){Rmax=my_Rmax; xmax = Rmax/rscale;};
	
	void initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, PosType my_slope, long *seed);

private:
	/// table size
	static const long NTABLE;
	/// maximum Rmax/rscale
	static const PosType maxrm;
	/// keeps track of how many time the tables are created, default is just once
	static int count;

	/// tables for lensing properties specific functions
	static PosType *mhattable,*xtable;
	/// make the specific tables
	void make_tables();
	/// interpolates from the specific tables
	PosType InterpolateFromTable(PosType y);

	/// read in parameters from a parameterfile in InputParams params
	void assignParams(InputParams& params);

	/// slope of the surface density profile
	PosType beta;

	// Override internal structure of halos
	inline PosType alpha_h(PosType x){
		return -1.0*InterpolateFromTable(x)/InterpolateFromTable(xmax);
	}
	inline KappaType kappa_h(PosType x){
		return 0.5*x*x/InterpolateFromTable(xmax)/pow(1+x,beta);
	}
	inline KappaType gamma_h(PosType x){
		return (0.5*x*x/pow(1+x,beta) - InterpolateFromTable(x))/InterpolateFromTable(xmax);
	}
	inline KappaType phi_h(PosType x){
        return -1.0*alpha_int(x)/InterpolateFromTable(xmax) ;
		//ERROR_MESSAGE();
		//std::cout << "time delay has not been fixed for PseudoNFW profile yet." << std::endl;
		//exit(1);
		//return 0.0;
	}
    inline KappaType phi_int(PosType x){
        return -1.0*alpha_int(x)/InterpolateFromTable(xmax) ;
    }
};





/** \ingroup DeflectionL2
 *
 * \brief A class for calculating the deflection, kappa and gamma caused by a collection of halos
 * with truncated power-law mass profiles.
 *
 * Derived from the TreeQuad class.  The "particles" are replaced with spherical halos.
 *The truncation is in 2d not 3d. \f$ \Sigma \propto r^{-\beta} \f$ so beta would usually be positive.
 *
 *
 * The default value of theta = 0.1 generally gives better than 1% accuracy on alpha.
 * The shear and kappa is always more accurate than the deflection.
 */
class LensHaloPowerLaw: public LensHalo{
public:
	LensHaloPowerLaw();
    LensHaloPowerLaw(float my_mass,float my_Rmax,PosType my_zlens,float my_rscale,PosType my_beta,float my_fratio,float my_pa,int my_stars_N);
	LensHaloPowerLaw(InputParams& params);
	~LensHaloPowerLaw();

	/// set the slope of the surface density profile
	void set_slope(PosType my_slope){beta=my_slope;};
    
  /// get slope
  PosType get_slope(){return beta;};
    
	/// initialize from a mass function
	void initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, PosType my_slope, long *seed);
    
    

private:
	/// read-in parameters from the parameter file
	void assignParams(InputParams& params);

	///	read in parameters from a parameterfile in InputParams params
	PosType beta;
  PosType fratio;
  PosType pa;
    

	// Override internal structure of halos
	inline PosType alpha_h(PosType x){
		if(x==0) x=1e-6*xmax;
		//assert(beta==2);
		//assert(-1.0*pow(x/xmax,beta+2) != 0.0);
        //std::cout << x << " " << beta << "  " << -1.0*pow(x/xmax,beta+2) << std::endl;
		return -1.0*pow(x/xmax,-beta+2);
	}
    /// this is not kappa to be overridden and there is an extra factor of M/pi/x^2 in force_halo
	inline KappaType kappa_h(PosType x){
		if(x==0) x=1e-6*xmax;
		//assert(0.5*(beta+2)*pow(x/xmax,beta)*x*x/(xmax*xmax) != 0);
        //std::cout << x << " " << beta << "  " << -1.0*pow(x/xmax,beta+2) << std::endl;
        
		return 0.5*(-beta+2)*pow(x/xmax,-beta)*x*x; //(xmax*xmax); 
	}
	inline KappaType gamma_h(PosType x){
		if(x==0) x=1e-6*xmax;
		//assert(0.5*beta*pow(x/xmax,beta+2) != 0);
        //std::cout << "gamma_h(" << 0.5*beta*pow(x/xmax,beta+2) << ")" << std::endl;
		return -0.5*beta*pow(x/xmax,-beta+2);
	}
	inline KappaType phi_h(PosType x){
		if(x==0) x=1e-6*xmax;
		return -1.0*pow(x/xmax,-beta+2)/(-beta+2);
	}
  inline KappaType phi_int(PosType x){
		//return alpha_int(x);
        return -1.0*pow(x/xmax,-beta+2)/(-beta+2);
    }
};






class LensHaloSimpleNSIE : public LensHalo{
public:
  /*LensHaloSimpleNSIE(){
    sigma = zlens = fratio = pa = rcore = 0.;
  }*/

  //LensHaloSimpleNSIE(float my_mass,float my_Rmax,PosType my_zlens,float my_rscale,float my_sigma, float my_rcore,float my_fratio,float my_pa,int my_stars_N);
  /// explicit constructor
  LensHaloSimpleNSIE(float my_mass,PosType my_zlens,float my_sigma
              ,float my_rcore,float my_fratio,float my_pa,int my_stars_N);
	LensHaloSimpleNSIE(InputParams& params);
	~LensHaloSimpleNSIE();

	/// overridden function to calculate the lensing properties
	void force_halo(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType const *xcm,bool subtract_point=false);
	// void force_halo(PosType *alpha,KappaType *kappa,KappaType *gamma,PosType *xcm,bool subtract_point=false);
    
	/// get the velocity dispersion
	float get_sigma(){return sigma;};
	/// get the NSIE radius
	float get_Rsize(){return Rsize;};
	/// get the axis ratio
	float get_fratio(){return fratio;};
	/// get the position angle
	float get_pa(){return pa;};
	/// get the core radius
	float get_rcore(){return rcore;};

	/// set the velocity dispersion
	void set_sigma(float my_sigma){sigma=my_sigma;};
	/// set the NSIE radius
	void set_Rsize(float my_Rsize){Rsize=my_Rsize;};
	///set the axis ratio
	void set_fratio(float my_fratio){fratio=my_fratio;};
	/// set the position angle
	void set_pa(float my_pa){pa=my_pa;};
	/// set the core radius
	void set_rcore(float my_rcore){rcore=my_rcore;};

protected:

	/// initialize from a simulation file
	//void initFromFile(float my_mass, long *seed, float vmax, float r_halfmass);
	/// initialize from a mass function
	//void initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, PosType my_slope, long *seed);
	/// simple initialize from mass while setting a random position angle and ellipticity
	//void initFromMass(float my_mass, long *seed);

	/// read-in parameters from a parameter file
	void assignParams(InputParams& params);

	/// velocity dispersion of NSIE
	float sigma;
	/// Actual edge of mass distribution in elliptical radius, Rmax is the range beyond which the halo is a point mass
	float Rsize;
	/// axis ratio of surface mass distribution
	float fratio;
	/// position angle on sky, radians
	float pa;
	/// core size of NSIE
	float rcore;

};

/** \ingroup DeflectionL2
 *
 * \brief A class for calculating the deflection, kappa and gamma caused by a collection of halos
 * with truncated Hernquist mass profiles.
 *
 * Derived from the TreeQuad class.  The "particles" are replaced with spherical halos.
 *The truncation is in 2d not 3d. \f$ \Sigma \propto r^\beta \f$ so beta would usually be negative.
 *
 *
 * The default value of theta = 0.1 generally gives better than 1% accuracy on alpha.
 * The shear and kappa is always more accurate than the deflection.
 */

class LensHaloHernquist: public LensHalo{
public:
	//LensHaloHernquist();
  LensHaloHernquist(float my_mass,float my_Rmax,PosType my_zlens,float my_rscale,float my_fratio,float my_pa,int my_stars_N);
	LensHaloHernquist(InputParams& params);
	virtual ~LensHaloHernquist();

    PosType ffunction(PosType x);
	PosType gfunction(PosType x);
	PosType hfunction(PosType x);
	PosType g2function(PosType x);

	/* the below functions alphaHern, kappaHern and gammaHern are obsolete and better to be deleted to avoid confusion
	void alphaHern(PosType *alpha,PosType *x,PosType Rtrunc,PosType mass,PosType r_scale
			,PosType *center,PosType Sigma_crit);
	KappaType kappaHern(PosType *x,PosType Rtrunc,PosType mass,PosType r_scale
			,PosType *center,PosType Sigma_crit);
	void gammaHern(KappaType *gamma,PosType *x,PosType Rtrunc,PosType mass,PosType r_scale
			,PosType *center,PosType Sigma_crit);
  */
	//void initFromFile(float my_mass, long *seed, float vmax, float r_halfmass);

	/// set Rmax
	void set_Rmax(float my_Rmax){Rmax=my_Rmax; xmax = Rmax/rscale; gmax = InterpolateFromTable(gtable,xmax);};
	/// set scale radius
	void set_rscale(float my_rscale){rscale=my_rscale; xmax = Rmax/rscale; gmax = InterpolateFromTable(gtable,xmax);};
   // friend struct Ialpha_func;
    


    
protected:
	/// table size
	static const long NTABLE;
	/// maximum Rmax/rscale
	static const PosType maxrm;
	/// keeps track of how many time the tables are created, default is just once
	static int count;

	/// tables for lensing properties specific functions
	static PosType *ftable,*gtable,*g2table,*htable,*xtable,*xgtable;
	/// make the specific tables
	void make_tables();
	/// interpolates from the specific tables
	PosType InterpolateFromTable(PosType *table, PosType y);

	/// read in parameters from a parameterfile in InputParams params
	void assignParams(InputParams& params);

	/// Override internal structure of halos
	inline PosType alpha_h(PosType x){
		return -0.25*InterpolateFromTable(gtable,x)/gmax;
	}
	inline KappaType kappa_h(PosType x){
		return 0.25*x*x*InterpolateFromTable(ftable,x)/gmax; // 0.5*
	}
	inline KappaType gamma_h(PosType x){
		return -0.25*x*x*InterpolateFromTable(g2table,x)/gmax;
	}
	inline KappaType phi_h(PosType x){
		return -0.25*InterpolateFromTable(htable,x)/gmax;
		//return -1.0*InterpolateFromTable(htable,x)/gmax;
	}
    inline KappaType phi_int(PosType x){
		return -0.25*InterpolateFromTable(xgtable,x)/gmax;
    }
    
    
private:
  PosType gmax;
};

/** \ingroup DeflectionL2
 *
 * \brief A class for calculating the deflection, kappa and gamma caused by a collection of halos
 * with truncated Jaffe mass profiles.
 *
 * Derived from the TreeQuad class.  The "particles" are replaced with spherical halos.
 *The truncation is in 2d not 3d. \f$ \Sigma \propto r^\beta \f$ so beta would usually be negative.
 *
 *
 * The default value of theta = 0.1 generally gives better than 1% accuracy on alpha.
 * The shear and kappa is always more accurate than the deflection.
 */

class LensHaloJaffe: public LensHalo{
public:
	//LensHaloJaffe();
  LensHaloJaffe(float my_mass,float my_Rmax,PosType my_zlens,float my_rscale,float my_fratio,float my_pa,int my_stars_N);
	LensHaloJaffe(InputParams& params);
	virtual ~LensHaloJaffe();

	/// set Rmax
	void set_Rmax(float my_Rmax){Rmax=my_Rmax; xmax = Rmax/rscale; gmax = InterpolateFromTable(gtable,xmax);};
	/// set scale radius
	void set_rscale(float my_rscale){rscale=my_rscale; xmax = Rmax/rscale; gmax = InterpolateFromTable(gtable,xmax);};

    PosType ffunction(PosType x);
	PosType gfunction(PosType x);
	PosType hfunction(PosType x);
	PosType g2function(PosType x);
    
protected:
    
	/// table size
	static const long NTABLE;
	/// maximum Rmax/rscale
	static const PosType maxrm;
	/// keeps track of how many time the tables are created, default is just once
	static int count;

	/// tables for lensing properties specific functions
	static PosType *ftable,*gtable,*g2table,*htable,*xtable,*xgtable;
	/// make the specific tables
	void make_tables();
	/// interpolates from the specific tables
	PosType InterpolateFromTable(PosType *table, PosType y);

	/// read in parameters from a parameterfile in InputParams params
	void assignParams(InputParams& params);

	/// Override internal structure of halos
	inline PosType alpha_h(PosType x){
		return -0.25*InterpolateFromTable(gtable,x)/gmax;
	}
	inline KappaType kappa_h(PosType x){
		return 0.125*x*x*InterpolateFromTable(ftable,x)/gmax;
	}
	inline KappaType gamma_h(PosType x){
		return -0.125*x*x*InterpolateFromTable(g2table,x)/gmax;
	}
    ///std::cout << "no analytic expression defined yet for Jaffe profile, take numerical" << std::endl;
	inline KappaType phi_h(PosType x){
        return -0.25*InterpolateFromTable(xgtable,x)/gmax;
    }
    inline KappaType phi_int(PosType x){
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
    LensHaloDummy(float my_mass,float my_Rmax,PosType my_zlens,float my_rscale, int my_stars_N);
	LensHaloDummy(InputParams& params);
	~LensHaloDummy(){};
	
	/// overridden function to calculate the lensing properties
	void force_halo(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType const *xcm,bool subtract_point=false);
    // void force_halo(PosType *alpha,KappaType *kappa,KappaType *gamma,PosType *xcm,bool subtract_point=false);
    
	/// initialize from a mass function
	void initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, PosType my_slope, long *seed);

	
private:
	/// read-in parameters from a parameter file
	void assignParams(InputParams& params);
	inline PosType alpha_h(PosType x){return  0.;}
	inline KappaType kappa_h(PosType x){return  0.;}
	inline KappaType gamma_h(PosType x){return  0.;}
};


typedef LensHalo* LensHaloHndl;


#endif /* LENS_HALOS_H_ */
