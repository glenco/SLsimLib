/*
 * lens_halos.h
 *
 *  Created on: 06.05.2013
 *      Author: mpetkova
 */

#ifndef LENS_HALOS_H_
#define LENS_HALOS_H_

#include "standard.h"
#include "tables.h"
#include "InputParams.h"
#include "source.h"

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
	~LensHalo();

	/// get the Rmax
	float get_Rmax(){return Rmax;};
	/// get the mass
	float get_mass(){return mass;};
	/// get the scale radius
	float get_rscale(){return rscale;};
	/// get the redshift
	double getZlens() const {return zlens;};

	/// initialize from a simulation file
	virtual void initFromFile(float my_mass, long *seed, float vmax, float r_halfmass){};
	/// initialize from a mass function
	virtual void initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, double my_slope, long *seed);

	/// set Rmax
	void set_Rmax(float my_Rmax){Rmax=my_Rmax; xmax = Rmax/rscale; };
	/// set mass
	void set_mass(float my_mass){mass=my_mass;};
	/// set scale radius
	void set_rscale(float my_rscale){rscale=my_rscale; xmax = Rmax/rscale;};
	/// set redshift
	void setZlens(double my_zlens){zlens=my_zlens;};
	/// set redshift, where the cosmology and the source redshift are needed (BaseNSIELensHalo)
	virtual void setZlens(CosmoHndl cosmo,double z,double dummy){zlens=z;};
	/// set slope
	virtual void set_slope(double my_slope){};

	/// set internal params that need either the cosmology or the source
	virtual void setInternalParams(CosmoHndl cosmo, SourceHndl source){};

	/// calculate the lensing properties -- deflection, convergence, and shear
	virtual void force_halo(double *alpha,KappaType *kappa,KappaType *gamma,double *xcm,bool no_kappa,bool subtract_point=false);

	/// internal compare redshift function
	bool compare(double z){return z > zlens;};

protected:
	/// read in parameters from a parameterfile in InputParams params
	void assignParams(InputParams& params);
	/// error message printout
	void error_message1(std::string name,std::string filename);

    float mass;
    /// Radius of halo and NSIE if it exists,  This is the radius used in the tree force solver
    /// to determine when a ray intersects an object.
    float Rmax;
    /// scale length or core size.  Different meaning in different cases.  Not used in NSIE case.
    float rscale;
    /// redshift
    double zlens;

    /// point mass case
	virtual double inline alpha_h(double x){return -1;};
	virtual KappaType inline kappa_h(double x){return 0;};
	virtual KappaType inline gamma_h(double x){return -2;};
	virtual KappaType inline phi_h(double x){return 0;};
  double xmax;
  
  // Functions for calculating axial dependence
  void setModesToEllip(double q,double theta);
  void faxial(double theta,double f[]);
  void gradial(double r,double g[]);
  void desymmeterize(double r,double theta,double *alpha,double *kappa,double *gamma);
  const static int Nmod = 18;
  double mod[18];
  double r_eps;
};

/** \ingroup DeflectionL2
 *
 * \brief A class for calculating the deflection, kappa and gamma caused by a collection of NFW
 * halos.
 *
 * Derived from the QuadTree class.  The "particles" are replaced with spherical NFW halos.
 *
 * This class uses the true expressions for the NFW profile.  This is
 * time consuming and not usually necessary. See QuadTreePseudoNFW for a faster alternative.
 *
* The default value of theta = 0.1 generally gives better than 1% accuracy on alpha.
 * The shear and kappa is always more accurate than the deflection.
 *
 */
class NFWLensHalo: public LensHalo{
public:
	NFWLensHalo();
	NFWLensHalo(InputParams& params);
	virtual ~NFWLensHalo();

	void initFromFile(float my_mass, long *seed, float vmax, float r_halfmass);

protected:
	/// table size
	const static long NTABLE = 1000;
	/// maximum Rmax/rscale
	const static double maxrm = 100.0;
	/// keeps track of how many time the tables are created, default is just once
	static int count;

	/// tables for lensing properties specific functions
	static double *ftable,*gtable,*g2table,*xtable;
	/// make the specific tables
	void make_tables();
	/// interpolates from the specific tables
	double InterpolateFromTable(double *table, double y);

	/// read in parameters from a parameterfile in InputParams params
	void assignParams(InputParams& params);

	/// Override internal structure of halos
	inline double alpha_h(double x){
		//return -1.0*InterpolateFromTable(gtable,x)/InterpolateFromTable(gtable,xmax);
		return -1.0*InterpolateFromTable(gtable,x)/gmax;
	}
	inline KappaType kappa_h(double x){
		return 0.5*x*x*InterpolateFromTable(ftable,x)/gmax;
	}
	inline KappaType gamma_h(double x){
		return -0.25*x*x*InterpolateFromTable(g2table,x)/gmax;
	}
	inline KappaType phi_h(double x){
		ERROR_MESSAGE();
		std::cout << "time delay has not been fixed for NFW profile yet." << std::endl;
		exit(1);
		return 0.0;
	}
  
private:
  double gmax;
};

/** \ingroup DeflectionL2
 *
 * \brief A class for calculating the deflection, kappa and gamma caused by a collection of
 * halos with a double power-law mass profile.
 *
 * Derived from the QuadTree class.  The "particles" are replaced with spherical halos
 * with \f$ \Sigma \propto 1/(1 + r/r_s )^\beta \f$ so beta would usually be positive.
 *
 * An NFW profile is approximated beta = 2 .
 *
 * The default value of theta = 0.1 generally gives better than 1% accuracy on alpha.
 * The shear and kappa is always more accurate than the deflection.
 */
class PseudoNFWLensHalo: public LensHalo{
public:
	PseudoNFWLensHalo();
	PseudoNFWLensHalo(InputParams& params);
	~PseudoNFWLensHalo();

	/// set the slope of the surface density profile
	void set_slope(double my_slope){beta=my_slope; make_tables();};
	/// initialize from a mass function
	void initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, double my_slope, long *seed);

private:
	/// table size
	const static long NTABLE = 1000;
	/// maximum Rmax/rscale
	const static double maxrm = 100.0;
	/// keeps track of how many time the tables are created, default is just once
	static int count;

	/// tables for lensing properties specific functions
	static double *mhattable,*xtable;
	/// make the specific tables
	void make_tables();
	/// interpolates from the specific tables
	double InterpolateFromTable(double y);

	/// read in parameters from a parameterfile in InputParams params
	void assignParams(InputParams& params);

	/// slope of the surface density profile
	double beta;

	// Override internal structure of halos
	inline double alpha_h(double x){
		return -1.0*InterpolateFromTable(x)/InterpolateFromTable(xmax);
	}
	inline KappaType kappa_h(double x){
		return 0.5*x*x/InterpolateFromTable(xmax)/pow(1+x,beta);
	}
	inline KappaType gamma_h(double x){
		return (0.5*x*x/pow(1+x,beta) - InterpolateFromTable(x))/InterpolateFromTable(xmax);
	}
	inline KappaType phi_h(double x){
		ERROR_MESSAGE();
		std::cout << "time delay has not been fixed for PseudoNFW profile yet." << std::endl;
		exit(1);
		return 0.0;
	}
};

/** \ingroup DeflectionL2
 *
 * \brief A class for calculating the deflection, kappa and gamma caused by a collection of halos
 * with truncated power-law mass profiles.
 *
 * Derived from the QuadTree class.  The "particles" are replaced with spherical halos.
 *The truncation is in 2d not 3d. \f$ \Sigma \propto r^\beta \f$ so beta would usually be negative.
 *
 *
 * The default value of theta = 0.1 generally gives better than 1% accuracy on alpha.
 * The shear and kappa is always more accurate than the deflection.
 */
class PowerLawLensHalo: public LensHalo{
public:
	PowerLawLensHalo();
	PowerLawLensHalo(InputParams& params);
	~PowerLawLensHalo();

	/// set the slope of the surface density profile
	void set_slope(double my_slope){beta=my_slope;};
	/// initialize from a mass function
	void initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, double my_slope, long *seed);

private:
	/// read-in parameters from the parameter file
	void assignParams(InputParams& params);

	///	read in parameters from a parameterfile in InputParams params
	double beta;

	// Override internal structure of halos
	inline double alpha_h(double x){
		if(x==0) x=1e-6*xmax;
		return -1.0*pow(x/xmax,beta+2);
	}
	inline KappaType kappa_h(double x){
		if(x==0) x=1e-6*xmax;
		return 0.5*(beta+2)*pow(x/xmax,beta)*x*x/(xmax*xmax);
	}
	inline KappaType gamma_h(double x){
		if(x==0) x=1e-6*xmax;
		return 0.5*beta*pow(x/xmax,beta+2);
	}
	inline KappaType phi_h(double x){
		ERROR_MESSAGE();
		std::cout << "time delay has not been fixed for PowerLaw profile yet." << std::endl;
		exit(1);
		return 0.0;
	}
};

class SimpleNSIELensHalo : public LensHalo{
public:
	SimpleNSIELensHalo();
	SimpleNSIELensHalo(InputParams& params);
	~SimpleNSIELensHalo();

	/// overridden function to calculate the lensing properties
	void force_halo(double *alpha,KappaType *kappa,KappaType *gamma,double *xcm,bool no_kappa,bool subtract_point=false);

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


	/// initialize from a simulation file
	void initFromFile(float my_mass, long *seed, float vmax, float r_halfmass);
	/// initialize from a mass function
	void initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, double my_slope, long *seed);
	/// simple initialize from mass
	void initFromMass(float my_mass, long *seed);

protected:
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


/**
 * \brief This is a lens that does no lensing.  It is useful for testing and for running refinement code on sources.
 */
class DummyLensHalo: public LensHalo{
public:
	DummyLensHalo();
	DummyLensHalo(InputParams& params);
	~DummyLensHalo(){};
	
	/// overridden function to calculate the lensing properties
	void force_halo(double *alpha,KappaType *kappa,KappaType *gamma,double *xcm,bool no_kappa,bool subtract_point=false);
	
private:
	/// read-in parameters from a parameter file
	void assignParams(InputParams& params);
};


typedef LensHalo* LensHaloHndl;


#endif /* LENS_HALOS_H_ */
