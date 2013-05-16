/*
 * lens_halos.h
 *
 *  Created on: 06.05.2013
 *      Author: mpetkova
 */

#ifndef LENS_HALOS_H_
#define LENS_HALOS_H_

#include "standard.h"
#include "simpleTree.h"
#include "tables.h"
#include "InputParams.h"

class LensHalo{
public:
	LensHalo();
	LensHalo(InputParams& params);
	virtual ~LensHalo();

	virtual float get_Rmax(){return Rmax;};
	virtual float get_mass(){return mass;};
	virtual float get_rscale(){return rscale;};

	double get_zlens(){return zlens;};

	virtual void set_internal(long*,float,float){};

	void set_Rmax(float my_Rmax){Rmax=my_Rmax;};
	void set_mass(float my_mass){mass=my_mass;};
	void set_rscale(float my_rscale){rscale=my_rscale;};
	void set_zlens(float my_zlens){zlens=my_zlens;};

	virtual void set_slope(double my_slope){};

	virtual void setInternalParams(CosmoHndl cosmo, double zsource){};

	virtual void force_halo(double *alpha,KappaType *kappa,KappaType *gamma,double *xcm,bool no_kappa);

protected:

	virtual void error_message1(std::string name,std::string filename);
	virtual void assignParams(InputParams& params);

    float mass;
    /// Radius of halo and NSIE if it exists,  This is the radius used in the tree force solver
    /// to determine when a ray intersects an object.
    float Rmax;
    /// scale length or core size.  Different meaning in different cases.  Not used in NSIE case.
    float rscale;
    /// redhsift
    double zlens;

    /// point mass case
	virtual double inline alpha_h(double x, double xmax){return -1;};
	virtual KappaType inline kappa_h(double x, double xmax){return 0;};
	virtual KappaType inline gamma_h(double x, double xmax){return -2;};
	virtual KappaType inline phi_h(double x, double xmax){return 0;};
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

	void set_internal(long*,float,float);

protected:
	void assignParams(InputParams& params);

	// Override internal structure of halos
	inline double alpha_h(double x,double xmax){
		return -1.0*InterpolateFromTable(gtable,xtable,x)/InterpolateFromTable(gtable,xtable,xmax);
	}
	inline KappaType kappa_h(double x,double xmax){
		return 0.5*x*x*InterpolateFromTable(ftable,xtable,x)/InterpolateFromTable(gtable,xtable,xmax);
	}
	inline KappaType gamma_h(double x,double xmax){
		return -0.25*x*x*InterpolateFromTable(g2table,xtable,x)/InterpolateFromTable(gtable,xtable,xmax);
	}
	inline KappaType phi_h(double x,double xmax){
		ERROR_MESSAGE();
		std::cout << "time delay has not been fixed for NFW profile yet." << std::endl;
		exit(1);
		return 0.0;
	}
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

	void set_slope(double my_slope){beta=my_slope;};

private:
	void assignParams(InputParams& params);

	double beta;

	// Override internal structure of halos
	inline double alpha_h(double x,double xmax){
		return -1.0*InterpolateFromTable(mhattable,xtable,x)/InterpolateFromTable(mhattable,xtable,xmax);
	}
	inline KappaType kappa_h(double x,double xmax){
		return 0.5*x*x/InterpolateFromTable(mhattable,xtable,xmax)/pow(1+x,beta);
	}
	inline KappaType gamma_h(double x,double xmax){
		return (0.5*x*x/pow(1+x,beta) - InterpolateFromTable(mhattable,xtable,x))/InterpolateFromTable(mhattable,xtable,xmax);
	}
	inline KappaType phi_h(double x,double xmax){
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

	void set_slope(double my_slope){beta=my_slope;};

private:
	void assignParams(InputParams& params);

	double beta; // logarithmic slope of 2d mass profile

	// Override internal structure of halos
	inline double alpha_h(double x,double xmax){
		if(x==0) x=1e-6*xmax;
		return -1.0*pow(x/xmax,beta+2);
	}
	inline KappaType kappa_h(double x,double xmax){
		if(x==0) x=1e-6*xmax;
		return 0.5*(beta+2)*pow(x/xmax,beta)*x*x/(xmax*xmax);
	}
	inline KappaType gamma_h(double x,double xmax){
		if(x==0) x=1e-6*xmax;
		return 0.5*beta*pow(x/xmax,beta+2);
	}
	inline KappaType phi_h(double x,double xmax){
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
	virtual ~SimpleNSIELensHalo();

	void force_halo(double *alpha,KappaType *kappa,KappaType *gamma,double *xcm,bool no_kappa);

	float get_sigma(){return sigma;};
	float get_Rsize(){return Rsize;};
	float get_fratio(){return fratio;};
	float get_pa(){return pa;};
	float get_rcore(){return rcore;};

	void set_sigma(float my_sigma){sigma=my_sigma;};
	void set_Rsize(float my_Rsize){Rsize=my_Rsize;};
	void set_fratio(float my_fratio){fratio=my_fratio;};
	void set_pa(float my_pa){pa=my_pa;};
	void set_rcore(float my_rcore){rcore=my_rcore;};

	float get_Rmax(){return Rmax*MAX(1.0,1.0/fratio);};

	void set_internal(long*,float,float);

protected:
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


typedef LensHalo* LensHaloHndl;


#endif /* LENS_HALOS_H_ */
