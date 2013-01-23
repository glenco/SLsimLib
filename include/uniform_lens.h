/*
 * uniform_lens.h
 *
 *  Created on: Jan 21, 2013
 *      Author: D. Leier
 */

#ifndef UNIFORM_LENS_H_
#define UNIFORM_LENS_H_

#include <lens.h>
#include <Tree.h>
#include <forceTree.h>
#include <quadTree.h>
#include <source.h>
#include <base_analens.h>

/**
 * \brief An uniform model to represent a lens on a single plane.
 *
 * The lens consists of a "host" lens which is a non-singular isothermal ellipsoid (NSIE) plus axial distortion
 * modes, substructures and stars.
 *
 *<pre>
 * Input Parameters:
 *
 *  ****  Uniform model parameters
 *	kappa_uniform
 *  gamma_uniform_1
 *  gamma_uniform_2
 *
 *  **** Distortion parameters
 *	NDistortionModes       Number of distortion modes to be used.  If zero the other distortion parameters are not needed.
 *	beta_perturb
 *	kappa_peturb
 *	gamma_peturb
 *	monopole_peturb
 *	quadrapole_peturb
 *	hexopole_peturb
 *	octopole_peturb
 *
 *  **** Substructure parameters
 *	NdensitySubstruct      Number density of substructures.  They are distributed uniformly.  If zero the other substructure parameters are not needed.
 *	beta_sub               Logorithmic slope of the internal clump profile.  Used if sub_type == powerlaw
 *	alpha_sub              Logorithmic slope of the mass function.
 *	R_submax               Maximum radius of most massive substructure (see Metcalf & Amara 2012)
 *	mass_max               Maximum mass
 *	mass_min               Minimum mass
 *	sub_type               Mass profile of clumps - 0 or nfw,1 or powerlaw, 2 or pointmass
 *
 *  **** Stars parameters
 *	Nstars                 Total number of stars that will be used in the simulation.  If zero the other star parameters are not needed.
 *	fstars                 Fraction of surface denity in stars.
 *	stars_mass             Mass of stars.
 *
 * The stars are not initially present.  They must be implanted later.
 *</pre>
 *
 */

class UniLens: public BaseAnaLens{
public:
	UniLens(InputParams& params);
	~UniLens();

	virtual void assignParams(InputParams& params);
	void PrintLens(bool show_substruct,bool show_stars);
/*	void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off){
		for(unsigned long i=0;i<Npoints;++i){
			i_points[i].image->kappa = i_points[i].kappa = kappa;
			i_points[i].image->gamma[0] = i_points[i].gamma[0] = gamma[0];
			i_points[i].image->gamma[1] = i_points[i].gamma[1] = gamma[1];
			//i_points[i].image->gamma[2] = i_points[i].gamma[2] = 0.0;
			i_points[i].image->x[0] = (1-kappa-gamma[0])*i_points[i].x[0]-gamma[1]*i_points[i].x[1];
			i_points[i].image->x[1] = (1-kappa+gamma[0])*i_points[i].x[1]-gamma[1]*i_points[i].x[0];
		}
	}
	void rayshooterInternal(double *ray, double *alpha, float *gamma, float *kappa, bool kappa_off){
		alpha[0] = alpha[1] = 0.0;
		gamma[0] = gamma[1] = 0.0;
		*kappa = 0.0;
	}*/


private:
   float kappa_uniform;
   float gamma_uniform[2];

};


#endif /* UNIFORM_LENS_H_ */
