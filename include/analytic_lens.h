/*
 * nsie.h
 *
 *  Created on: Dec 8, 2009
 *      Author: R.B. Metcalf
 */

#ifndef analens_declare
#define analens_declare

#include <lens.h>
#include <Tree.h>
#include <forceTree.h>
#include <quadTree.h>
#include <source.h>
#include <base_analens.h>

/**
 * \brief An "analytic" model to represent a lens on a single plane.
 *
 * The lens consists of a "host" lens which is a non-singular isothermal ellipsoid (NSIE) plus axial distortion
 * modes, substructures and stars.
 *
 *<pre>
 * Input Parameters:
 *
 *  **** NSIE parameters
 * 	sigma                  Velocity dispersion of host NSIE.
 *	core                   Core size
 *	axis_ratio             Axis ratio of mass
 *	pos_angle              Position angle
 *	z_lens                 Redshift of lens
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
 *	beta_sub               Logarithmic slope of the internal clump profile.  Used if sub_type == powerlaw
 *	alpha_sub              Logarithmic slope of the mass function.
 *	R_submax               Maximum radius of most massive substructure (see Metcalf & Amara 2012)
 *	mass_max               Maximum mass
 *	mass_min               Minimum mass
 *	sub_type               Mass profile of clumps - 0 or nfw,1 or powerlaw, 2 or pointmass
 *
 *  **** Stars parameters
 *	Nstars                 Total number of stars that will be used in the simulation.  If zero the other star parameters are not needed.
 *	fstars                 Fraction of surface density in stars.
 *	stars_mass             Mass of stars.
 *
 * The stars are not initially present.  They must be implanted later.
 *</pre>
 *
 * TODO BEN finish this documentation.
 */
class AnaLens : public BaseAnaLens{
public:
	AnaLens(InputParams& params);
	~AnaLens();

  virtual void assignParams(InputParams& params);
  double FractionWithinRe(double rangeInRei);
  void PrintLens(bool show_substruct,bool show_stars);

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

private:
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


// in mark_points.c
void MarkPoints(TreeHndl s_tree,AnaLens *lens,bool sb_cut,short invert);
void _MarkPoints(TreeHndl s_tree,AnaLens *lens,bool *sbcut);
bool InSource(double *ray,AnaLens *lens,bool surfacebright);

void find_lens(int Nimages,int Nsources,int *pairing,double **xob,double *xg,double beta
   		 ,int N,int *degen,double *mod,double **v,double **dx);

#endif
