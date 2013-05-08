/*
 * forceTree.h
 *
 *  Created on: Oct 14, 2011
 *      Author: bmetcalf
 */

#ifndef FORCE_TREE_H_
#define FORCE_TREE_H_

#include "simpleTree.h"

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

	ForceTree(PosType **xp,IndexType Npoints,LensHalo *my_halos
			,bool Multisize = true,double my_kappa_bk=0.0,int bucket = 5,int dimensions = 2
			,bool median = false,PosType theta = 0.1
			);

	~ForceTree();

	/// calculated sph smoothing and store them in the tree, also provide pointer to them
	float * CalculateSPHsmoothing(int N);
	/// calculate the deflection and lensing properties
	void force2D(double *ray,double *alpha,KappaType *kappa,KappaType *gamma,bool no_kappa = true);
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

	//PosType (*alpha_particle)(PosType r,float xmax);
	//PosType (*kappa_particle)(PosType r,float xmax);
	//PosType (*gamma_particle)(PosType r,float xmax);

	//PosType (*alpha_halo)(PosType r,HaloInternal &par);
	//PosType (*kappa_halo)(PosType r,HaloInternal &par);
	//PosType (*gamma_halo)(PosType r,HaloInternal &par);

	void CalcMoments();
	void rotate_coordinates(PosType **coord);
	//void spread_particles();

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
	LensHalo *halos;
};

typedef ForceTree *ForceTreeHndl;

#endif /* FORCE_TREE_H_ */
