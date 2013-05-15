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

class BaseNSIELensHalo : public LensHalo{
public:
	BaseNSIELensHalo();
	BaseNSIELensHalo(InputParams& params);
	virtual ~BaseNSIELensHalo();

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

class QuadTree;

class NSIELensHalo : public BaseNSIELensHalo{
public:
	NSIELensHalo(InputParams& params);
	virtual ~NSIELensHalo();

	bool AreSubStructImaplated(){return substruct_implanted;};
	bool AreStarsImaplated(){return stars_implanted;};

	void force_halo(double *alpha,KappaType *kappa,KappaType *gamma,double *xcm,bool no_kappa);

	void setInternalParams(CosmoHndl cosmo, double zsource);

	// These need to be in the base class because they are used in the rayshooter function which because of parrelleization is not a member function
	double getEinstein_ro(){return Einstein_ro;}

	double getPerturb_beta(){return perturb_beta;}

	int getPerturb_Nmodes(){return perturb_Nmodes;}    /// this includes two for external shear

	double averageSubMass();
	void reNormSubstructure(double kappa_sub);
	void substract_stars_disks(PosType *ray,PosType *alpha
			,KappaType *kappa,KappaType *gamma);
	void implant_stars(Point *centers,unsigned long Nregions,long *seed,int mftype=0);
	float* stellar_mass_function(int mftype, unsigned long Nstars, long *seed, double minmass=0.1, double maxmass=100
			,double bendmass=1.0, double powerlo=-0.3, double powerhi=-2.35);

protected:
	/// substructures
	double sub_sigmaScale;
	double sub_Ndensity;
	/// actual number of substructures
	int sub_N;
	double **sub_x;
	LensHalo *subs;
	/// slope of mass profile
	double sub_beta;
	/// slope of mass function
	double sub_alpha;
	/// radius of largest mass substructures
	double sub_Rmax;
	double sub_Mmax;
	double sub_Mmin;
	double sub_theta_force;
	QuadTree *sub_tree;
	IndexType *sub_substructures;
	ClumpInternal sub_type;

	/// stars
	int stars_N;
	IndexType *stars;
	PosType **stars_xp;
	//ForceTree *star_tree;
	QuadTree *star_tree;
	double star_massscale;
	/// star masses relative to star_massscles
	float *star_masses;
	double star_fstars;
	double star_theta_force;
	int star_Nregions;
	double *star_region;

	void assignParams(InputParams& params);

	void error_message1(std::string name,std::string filename);

	// in readlens_ana.c

	double *perturb_modes;  ///first two are shear

	// perturbations to host.  These are protected so that in some derived classes they can or cann't be changed.
	int perturb_Nmodes;    /// this includes two for external shear
	double perturb_beta;
	double *perturb_rms;

	bool substruct_implanted;

	bool stars_implanted;
	/// Number of regions to be subtracted to compensate for the mass in stars

	double *star_kappa;
	double **star_xdisk;

	// private derived quantities
	double Sigma_crit;
	double Einstein_ro;
	/// private: conversion factor between Mpc on the lens plane and arcseconds
	double MpcToAsec;
	double to;
};

class AnaNSIELensHalo : public NSIELensHalo{
public:
	AnaNSIELensHalo(InputParams& params);
	~AnaNSIELensHalo();

	double FractionWithinRe(double rangeInRei);

	// in randoimize_lens.c
	void RandomizeHost(long *seed,bool tables);
	void RandomizeSigma(long *seed,bool tables);
	void RandomlyDistortLens(long *seed,int Nmodes);
	void AlignedRandomlyDistortLens(long *seed,double theta,int n);
	void RandomizeSubstructure(double rangeInRei,long *seed);
	void RandomizeSubstructure2(double rangeInRei,long *seed);
	void RandomizeSubstructure3(double rangeInRei,long *seed);

	void FindLensSimple(int Nimages,Point *image_positions,double *y,double **dx_sub);
	void FindLensSimple(ImageInfo *imageinfo ,int Nimages ,double *y,double **dx_sub);

protected:

	void assignParams(InputParams& params);

	// Things added to manipulate and fit lenses.
	int check_model(int Nimages,int Nsources,int Nlenses,int *pairing,double **xob,double *x_center
			,int Nmod,double *mod,double **xg,double Re2,double **dx_sub,double **Amag,double ytol);
	double find_axis(double *mod,int Nmod);
	double deflect_translated(double beta,double *mod,double *x,double *y,double *mag,int N
			,int Nlenses,double Re2,double *x2);
	double ElliptisizeLens(int Nimages,int Nsources,int Nlenses,int *pairing,double **xob
			,double *xc,double **xg,double sigG,double beta,int Nmod
			,double *mod,double **dx,double *re2,double *q);
};

class UniNSIELensHalo : public NSIELensHalo{
public:
	UniNSIELensHalo(InputParams& params);
	~UniNSIELensHalo();

	void implant_stars(double x,double y,unsigned long Nregions,long *seed,int mftype=0);
	float getKappa_uniform(){return kappa_uniform;}
	float* getGamma_uniform(){return gamma_uniform;}
	double getAveMag(){ return 1.0/( pow(1-kappa_uniform,2) - gamma_uniform[0]*gamma_uniform[0] - gamma_uniform[1]*gamma_uniform[1]);}

protected:
	void assignParams(InputParams& params);

	float kappa_uniform;
	float gamma_uniform[3];
};

typedef LensHalo* LensHaloHndl;


namespace Utilities{
	double RandomFromTable(double *table,unsigned long Ntable,long *seed);
	void rotation(float *xout,float *xin,double theta);
	void rotation(double *xout,double *xin,double theta);
}


double lens_expand(double beta,double *mod,int Nmodes,double *x,double *alpha,KappaType *gamma,KappaType *phi);

void alphaNSIE(double *alpha,double *xt,double f,double bc,double theta);
KappaType kappaNSIE(double *xt,double f,double bc,double theta);
void gammaNSIE(KappaType *gam,double *xt,double f,double bc,double theta);
KappaType invmagNSIE(double *x,double f,double bc,double theta
                     ,float *gam,float kap);
double rmaxNSIE(double sigma,double mass,double f,double rc );
double ellipticRadiusNSIE(double *x,double f,double pa);
void quadMomNSIE(float mass,float Rmax,float f,float rc,float theta,double *quad);

//  in powerlow.c

void alphaPowLaw(double *alpha,double *x,double R,double mass,double beta,double *center,double Sigma_crit);
KappaType kappaPowLaw(double *x,double R,double mass,double beta,double *center,double Sigma_crit);
void gammaPowLaw(KappaType *gamma,double *x,double R,double mass,double beta,double *center,double Sigma_crit);
KappaType phiPowLaw(double *x,double R,double mass,double beta,double *center,double Sigma_crit);

// in nfw_lens.c
void alphaNFW(double *alpha,double *x,double Rtrunc,double mass,double r_scale
                ,double *center,double Sigma_crit);
KappaType kappaNFW(double *x,double Rtrunc,double mass,double r_scale
                ,double *center,double Sigma_crit);
void gammaNFW(KappaType *gamma,double *x,double Rtrunc,double mass,double r_scale
                ,double *center,double Sigma_crit);
#endif /* LENS_HALOS_H_ */
