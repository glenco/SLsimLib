/*
 * nsie.h
 *
 *  Created on: Dec 8, 2009
 *      Author: R.B. Metcalf
 */

#ifndef analens_declare
#define analens_declare

#include <Tree.h>
#include <source.h>
#include <base_analens.h>

/**
 * \brief An "analytic" model to represent a lens on a single plane.
 *
 * The lens consists of a "host" lens which is a non-singular isothermal ellipsoid (NSIE) plus axial distortion
 * modes, substructures and stars.  LensHaloAnaNSIE differs from a LensHaloBaseNSIE in that there are additional 
 * functions for fitting the lens to image positions and for giving the lens random substructures and distortions.
 *
 *<pre>
 * Input Parameters:
 *
 *  **** NSIE parameters
 * 	main_sigma                  Velocity dispersion of host NSIE.
 *	main_core                   Core size
 *	main_axis_ratio             Axis ratio of mass
 *	main_pos_angle              Position angle in radiants
 *	main_z_lens                 Redshift of lens
 *
 *  **** Distortion parameters
 *	main_NDistortionModes       Number of distortion modes to be used.  If zero the other distortion parameters are not needed.
 *	main_perturb_beta
 *	kappa_peturb
 *	gamma_peturb
 *	monopole_peturb
 *	quadrapole_peturb
 *	hexopole_peturb
 *	octopole_peturb
 *
 *  **** Substructure parameters
 *	main_sub_Ndensity      Number density of substructures.  They are distributed uniformly.  If zero the other substructure parameters are not needed.
 *	main_sub_beta               Logarithmic slope of the internal clump profile.  Used if main_sub_type == powerlaw
 *	main_sub_alpha              Logarithmic slope of the mass function.
 *	main_sub_Rmax               Maximum radius of most massive substructure (see Metcalf & Amara 2012)
 *	main_sub_mass_max               Maximum mass
 *	main_sub_mass_min               Minimum mass
 *	main_sub_type               Mass profile of clumps - 0 or nfw,1 or powerlaw, 2 or pointmass
 *
 *  **** Stars parameters
 *	main_stars_N                 Total number of stars that will be used in the simulation.  If zero the other star parameters are not needed.
 *	main_stars_fraction                 Fraction of surface density in stars.
 *	main_stars_mass             Mass of stars.
 *
 * The stars are not initially present.  They must be implanted later.
 *</pre>
 *
 */
// TODO: BEN finish this documentation for perturbation parameters.
class LensHaloAnaNSIE : public LensHaloBaseNSIE{
public:
	LensHaloAnaNSIE(InputParams& params);
	~LensHaloAnaNSIE();

  void assignParams(InputParams& params);
  double FractionWithinRe(double rangeInRei);
  void PrintLens(bool show_substruct,bool show_stars);
  void setCosmology(const COSMOLOGY& cosmo);

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
void MarkPoints(TreeHndl s_tree,LensHaloAnaNSIE *lens,bool sb_cut,short invert);
void _MarkPoints(TreeHndl s_tree,LensHaloAnaNSIE *lens,bool *sbcut);
bool InSource(double *ray,LensHaloAnaNSIE *lens,bool surfacebright);

void find_lens(int Nimages,int Nsources,int *pairing,double **xob,double *xg,double beta
		 ,int N,int *degen,double *mod,double **v,double **dx);

#endif
