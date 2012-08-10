/*
 * forceTree.h
 *
 *  Created on: Oct 14, 2011
 *      Author: bmetcalf
 */

#ifndef FORCE_TREE_H_
#define FORCE_TREE_H_

#include <simpleTree.h>

//enum PartProf {gaussian,powerlaw};

/** \ingroup DeflectionL2
 *
 * \brief Object used to calculate the force or deflection caused by a collection
 * of "particles" by the tree method.
 *
 * The particles can be point masses or have multiple sizes in which case they have a Gaussian profile.
 * They can also have the same mass or multiple masses.
 *
 * xp[][], masses[] and rsph[] need to be allocated before constructing a ForceTree is
 * constructed and de-allocated after it is destruction.  If the boolian flags are set to
 * false these arrays need only be one element long.  Multiple ForceTrees can be
 * made from the same particles.  Do not rotate the particles without reconstructing
 * a ForceTree.
 *
 * Most of the code in the earlier TreeNBForce.c is duplicated here as private methods and
 * a few public ones.
 *
 * The default value of theta = 0.1 generally gives better than 1% accuracy on alpha.
 * The shear and kappa is always more accurate than the deflection.
 *
 */
class ForceTree : public SimpleTree {
public:
	ForceTree(PosType **xp,IndexType Npoints,float *masses,float *rsph,bool Multimass,bool Multisize
			,double my_kappa_background = 0,int bucket = 5,int dimensions = 2,bool median = false,PosType theta = 0.1
			);

	virtual ~ForceTree();

	/// calculated sph smoothing and store them in the tree, also provide pointer to them
	float * CalculateSPHsmoothing(int N);
	/// calculate the deflection and lensing properties
	void force2D(PosType *ray,PosType *alpha,PosType *kappa,PosType *gamma,bool no_kappa = true);
	/// provides a way to change the profiles of the particles, by default Gaussian
	//void ChangeParticleProfile(PartProf partprof);


protected:

	bool init;

	/// true if particles have different masses.
	bool MultiMass;
	/// true if particles have different sizes.
	bool MultiRadius;
	/// Array of particle masses
	float *masses;
	/// Array of particle sizes
	float *rsph;
	/// A uniform mass sheet in units of mass_scale/Mpc^2 used to subtract of the contribution
	/// of the particles to the mean density of the universe
	double kappa_background;

	PosType force_theta;

	//PosType (*alpha_particle)(PosType r,float rmax);
	//PosType (*kappa_particle)(PosType r,float rmax);
	//PosType (*gamma_particle)(PosType r,float rmax);

	//PosType (*alpha_halo)(PosType r,HaloInternal &par);
	//PosType (*kappa_halo)(PosType r,HaloInternal &par);
	//PosType (*gamma_halo)(PosType r,HaloInternal &par);

	void CalcMoments();
	void rotate_coordinates(PosType **coord);
	//void spread_particles();

	// Internal profiles for a Gaussian particle
	inline virtual double alpha_h(double r2s2,double sigma){
	  return (sigma > 0.0 ) ? ( exp(-0.5*r2s2) - 1.0 ) : -1.0;
	}
	inline virtual double kappa_h(double r2s2,double sigma){
	  return 0.5*r2s2*exp(-0.5*r2s2);
	}
	inline virtual double gamma_h(double r2s2,double sigma){
	  return (sigma > 0.0 ) ? (-2.0 + (2.0 + r2s2)*exp(-0.5*r2s2) ) : -2.0;
	}
	inline virtual double phi_o(double r2,double sigma){
		ERROR_MESSAGE();  // not yet written
		std::cout << "time delay has not been fixed fot this profile yet." << std::endl;
		exit(1);
		return 0;
	}

	bool haloON;
	HaloStructure *halo_params;
};

typedef ForceTree *ForceTreeHndl;

/** \ingroup DeflectionL2
 *
 * \brief A class for calculating the deflection, kappa and gamma caused by a collection of halos
 * with truncated power-law mass profiles.
 *
 * Derived from the ForceTree class.  The "particles" are replaced with spherical halos.
 *The truncation is in 2d not 3d. \f$ \Sigma \propto r^\beta \f$ so beta would usually be negative.
 *
 * The default value of theta = 0.1 generally gives better than 1% accuracy on alpha, but
 * only ~ 20% accuracy on gamma.  For high accuracy gamma use theta <~ 0.01
 */
class ForceTreePowerLaw : public ForceTree{

public:
	ForceTreePowerLaw(float beta,PosType **xp,IndexType Npoints,HaloStructure *par_internals
			,bool Multisize = true,double my_kappa_bk=0.0,int bucket = 5,int dimensions = 2
			,bool median = false,PosType theta = 0.1
			);
	~ForceTreePowerLaw();

private:
	float beta; // logorithmic slop of 2d mass profile

	// Override internal structure of halos
	inline double alpha_h(double r,double rm){
		return (r < rm) ? -1.0*pow(r/rm,beta+2) : -1.0;
	}
	inline double kappa_h(double r,double rm){
		return (r < rm) ? (beta+2)*pow(r/rm,beta)*r*r/(rm*rm) : 0.0;
	}
	inline double gamma_h(double r,double rm){
		return (r < rm) ? -0.5*beta*pow(r/rm,beta+2) : -2.0;
	}
	inline double phi_h(double r,double rm){
		ERROR_MESSAGE();
		std::cout << "time delay has not been fixed for PowerLaw profile yet." << std::endl;
		exit(1);
		return 0.0;
	}
};

/** \ingroup DeflectionL2
 *
 * \brief A class for calculating the deflection, kappa and gamma caused by a collection of NFW
 * halos.
 *
 * Derived from the ForceTree class.  The "particles" are replaced with spherical NFW halos.
 *
 * This class uses the true expressions for the NFW profile.  This is
 * time consuming and not usually necessary. See ForceTreePseudoNFW for a faster alternative.
 *
 * The default value of theta = 0.1 generally gives better than 1% accuracy on alpha, but
 * only ~ 20% accuracy on gamma.  For high accuracy gamma use theta <~ 0.01
 */
class ForceTreeNFW : public ForceTree{

public:
	ForceTreeNFW(PosType **xp,IndexType Npoints,HaloStructure *par_internals
			,bool Multisize = true,double my_kappa_bk = 0.0
			,int bucket = 5,int dimensions = 2,bool median = false,PosType theta = 0.1
			);
	~ForceTreeNFW();

private:

	double *ft, *gt, *g2t;

	// Override internal structure of halos
	inline double alpha_h(double r,double rm){
		return (r < rm) ? -1.0*InterpolateFromTable(gt,r)/InterpolateFromTable(gt,rm) : -1.0;
	}
	inline double kappa_h(double r,double rm){
		return (r < rm) ? 0.5*r*r*InterpolateFromTable(ft,r)/InterpolateFromTable(gt,rm) : 0.0;
	}
	inline double gamma_h(double r,double rm){
		return (r < rm) ? -0.25*r*r*InterpolateFromTable(g2t,r)/InterpolateFromTable(gt,rm) : -2.0;
	}
	inline double phi_h(double r,double rm){
		ERROR_MESSAGE();
		std::cout << "time delay has not been fixed for NFW profile yet." << std::endl;
		exit(1);
		return 0.0;
	}
	void point_tables();

};

double gfunction(double x);
double ffunction(double x);
double g2function(double x);
double rhos(double x);

/** \ingroup DeflectionL2
 *
 * \brief A class for calculating the deflection, kappa and gamma caused by a collection of
 * halos with a double power-law mass profile.
 *
 * Derived from the ForceTree class.  The "particles" are replaced with spherical halos
 * with \f$ \Sigma \propto 1/(1 + r/r_s )^\beta \f$ so beta would usually be positive.
 *
 * An NFW profile is approximated beta = 2 .
 *
 * The default value of theta = 0.1 generally gives better than 1% accuracy on alpha, but
 * only ~ 20% accuracy on gamma.  For high accuracy gamma use theta <~ 0.01
 */
class ForceTreePseudoNFW : public ForceTree{

public:
	ForceTreePseudoNFW(double beta,PosType **xp,IndexType Npoints,HaloStructure *par_internals
			,bool Multisize = true,double my_kappa_bk = 0.0,int bucket = 5,int dimensions = 2
			,bool median = false,PosType theta = 0.1
			);
	~ForceTreePseudoNFW();

private:

	double beta;
	double *mhatt;

	// Override internal structure of halos
	inline double alpha_h(double r,double rm){
		return (r < rm) ? -1.0*InterpolateFromTable(mhatt,r)/InterpolateFromTable(mhatt,rm) : -1.0;
	}
	inline double kappa_h(double r,double rm){
		return (r < rm) ? 0.5*r*r/InterpolateFromTable(mhatt,rm)/pow(1+r,beta) : 0.0;
	}
	inline double gamma_h(double r,double rm){
		return (r < rm) ? (0.5*r*r/pow(1+r,beta) - InterpolateFromTable(mhatt,r))/InterpolateFromTable(mhatt,rm) : -2.0;
	}
	inline double phi_h(double r,double rm){
		ERROR_MESSAGE();
		std::cout << "time delay has not been fixed for PseudoNFW profile yet." << std::endl;
		exit(1);
		return 0.0;
	}
	void point_tables();
};

double mhat(double y, double beta);

#endif /* FORCE_TREE_H_ */
