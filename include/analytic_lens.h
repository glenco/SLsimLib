/*
 * nsie.h
 *
 *  Created on: Dec 8, 2009
 *      Author: R.B. Metcalf
 */

#ifndef analens_declare
#define analens_declare

#include "Tree.h"
#include "source.h"
#include "base_analens.h"

/**  \brief LensHalo class primarily used for fitting point source lenses.
 
 The model is a power law with axial modes made to be as elliptical as possible.
 */
// TODO: BEN finish this documentation for perturbation parameters.
class LensHaloFit : public LensHaloBaseNSIE{
public:
  //LensHaloAnaNSIE(InputParams& params,bool verbose = false);
  /// Creates a AnaLens which initially has no mass,  Use FindLensSimple() to give it mass
  LensHaloFit(const COSMOLOGY& cosmo, int MyNmodes, PosType beta,PosType zlensref, PosType zsourceref);
	~LensHaloFit();

  void PrintLens(bool show_substruct,bool show_stars);

  void FindLensSimple(int Nimages,Point *image_positions,double *y,double **dx_sub);
  void FindLensSimple(ImageInfo *imageinfo ,int Nimages ,double *y,double **dx_sub);
  
  // these need to be written so that they translate between modes and these quantities
  /// get the velocity dispersion
  virtual PosType get_sigma(){return 0.0;};
  /// get the axis ratio
  virtual PosType get_fratio(){return 0.0;};
  /// get the position angle
  virtual PosType get_pa(){return 0.0;};
  /// get the core radius
  virtual PosType get_rcore(){return 0.0;};
  

  /// set the number of perturbation modes -- Does the same as LensHaloBaseNSIE::getPerturb_Nmodes().
  void setNmodes(int my_Nmodes){perturb_Nmodes = my_Nmodes;};

  /// get the number of perturbation modes
  int getNmodes(){return perturb_Nmodes;};
  /// set the perturbation modes
  void set_perturbmodes(PosType * ListModes, const int Nmodes);
  /// get the perturbation modes
  void get_perturbmodes(PosType * ListModes, const int Nmodes);
  /// get the ouput of ElliptisizeLens :
  double * getq() { return qpriv; };

  /// set the velocity dispersion
  void set_sigma(PosType my_sigma){sigma = my_sigma; };
  /// set the axis ratio
  void set_fratio(PosType my_fratio){fratio = my_fratio;};
  
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
  void setCosmology(const COSMOLOGY& cosmo);
  
  // output of ElliptisizeLens
  double qpriv[7];

  //void find_lens(int Nimages,int Nsources,int *pairing,double **xob,double *xg,double beta
  //               ,int N,int *degen,double *mod,double **v,double **dx);

  //void assignParams(InputParams& params);
};

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
  LensHaloAnaNSIE(InputParams& params,bool verbose = false);
  /// Creates a AnaLens which initially has no mass,  Use FindLensSimple() to give it mass
  //LensHaloAnaNSIE(const COSMOLOGY& cosmo);
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
  
  /// get the velocity dispersion
  virtual PosType get_sigma(){return sigma;};
  /// get the NSIE radius
  //PosType get_Rsize(){return Rsize;};
  /// get the axis ratio
  virtual PosType get_fratio(){return fratio;};
  /// get the position angle
  virtual PosType get_pa(){return pa;};
  /// get the core radius
  virtual PosType get_rcore(){return rcore;};
  

private:
  
  // Things added to manipulate and fit lenses.
  int check_model(int Nimages,int Nsources,int Nlenses,int *pairing,double **xob,double *x_center
                  ,int Nmod,double *mod,double **xg,double Re2,double **dx_sub,double **Amag,double ytol);
  double find_axis(double *mod,int Nmod);
  double deflect_translated(double beta,double *mod,double *x,double *y,double *mag,int N
                            ,int Nlenses,double Re2,double *x2);
  
};


// in mark_points.c
void MarkPoints(TreeHndl s_tree,LensHaloAnaNSIE *lens,bool sb_cut,short invert);
void _MarkPoints(TreeHndl s_tree,LensHaloAnaNSIE *lens,bool *sbcut);
bool InSource(double *ray,LensHaloAnaNSIE *lens,bool surfacebright);


#endif
