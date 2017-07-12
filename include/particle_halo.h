//
//  particle_halo.h
//  GLAMER
//
//  Created by bmetcalf on 16/06/15.
//
//

#ifndef GLAMER_particle_halo_h
#define GLAMER_particle_halo_h

#include "quadTree.h"
#include "geometry.h"

/**
 *  \brief A class that represents the lensing by a collection of simulation particles.
 
   Smoothing is done according to the density of particles in 3D.  Smoothing sizes are
   either read in from a file (names simulation_filename + "." + Nsmooth + "sizes") or calculated
   if the file does not exist (in which case the file is created).  This can be
   time and memory consuming when there are a large number of particles.
 
   Input format:
     ASCCI - a table of three floats for positions in comoving Mpc (no h factor).
             The lines "# nparticles ...." and "# mass ...." must be in
             header at the top of the file. # is otherwise a comment character.
             Only one type of particle in a single input file.
 
    More input formats will be added in the future.
*/
class LensHaloParticles : public LensHalo
{
public:
  LensHaloParticles(
                    const std::string& simulation_filename  /// name of file containing particles
                    ,PosType redshift     /// redshift of particles
                    ,int Nsmooth          /// number of neighbours for the smoothing
                    ,const COSMOLOGY& lenscosmo /// cosmology
                    ,Point_2d theta_rotate      /// rotation of simulation, x-axis ratation angle and then y-axis rotation angle
                    ,bool recenter = false     /// re-center the coordinates to the center of mass
                    ,bool multimass = false     /** allows the particles to have different masses
                                                 , the masses for each particles must be provided as a 4th column in the particle data file */
                    ,PosType MinPSize = 0.0    /// Minimum smoothing size of particles
                    );
  
  LensHaloParticles(
                    PosType **positions    /// 3d positions in physical coordinates, only first two are used
                    ,size_t Nparticles    /// number of particles
                    ,std::vector<float> &my_sizes  /// smoothing sizes, will be swapped so returns empty, if size 1 the first size will be used for all the particles
                    ,std::vector<float> &my_masses /// masses, will be swapped so returns empty, if size 1 the first mass will be used for all the particles
                    ,PosType redshift     /// redshift of origin
                    ,const COSMOLOGY& cosmo /// cosmology
                    ,PosType sigma_back  /// background mass sheet
                    );
  
  ~LensHaloParticles();
  
  void force_halo(double *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,double const *xcm
                  ,bool subtract_point=false,PosType screening = 1.0);

  size_t getN() const { return Npoints; };
  
  /// rotate the simulation around the origin of the simulation coordinates, (radians)
  void rotate(Point_2d theta);
  
  /// center of mass in input coordinates
  Point_3d CenterOfMass(){return mcenter;}
  /** \brief This is a test class that makes a truncated SIE out of particles and puts it into a file in the right format for constructing a LensHaloParticles.
   
   This is useful for calculating the level of shot noise and finite source size.  The particles are distributed in 3D according to the SIE profile with only the perpendicular coordinates (1st and 2nd) distorted into an elliptical shape. If the halo is rotated from the original orientation it will not be a oblate spheroid.
   */
  static void makeSIE(
                      std::string new_filename  /// file name to store the particles
                      ,PosType redshift     /// redshift of particles
                      ,double particle_mass /// particle mass
                      ,double total_mass  /// total mass of SIE
                      ,double sigma       /// velocity dispersion in km/s
                      ,double q  /// axis ratio
                      ,Utilities::RandomNumbers_NR &ran
                      );
  
  void setSigmaBackground(PosType sigma_back){qtree->SetSigmaBackground(sigma_back); }
  
  static void find_smoothing(PosType **xp,size_t N,std::vector<float> &s,int Nneighbores);

private:

  static void find_smoothing_(TreeSimple *tree3d,PosType **xp,float *sizes,size_t N,int Nsmooth);

  Point_3d mcenter;
  void rotate_particles(PosType theta_x,PosType theta_y);

  //void calculate_smoothing(int Nsmooth);
  //void smooth_(TreeSimple *tree3d,PosType **xp,float *sizes,size_t N,int Nsmooth);

  void readPositionFileASCII(const std::string& filename);
  bool readSizesFile(const std::string& filename,int Nsmooth,PosType min_size);
  void writeSizes(const std::string& filename,int Nsmooth);
  
  void assignParams(InputParams& params);

  PosType **xp;
  std::vector<float> masses;
  std::vector<float> sizes;
  PosType min_size;
  bool multimass;
  
  Utilities::Geometry::SphericalPoint center;
  
  size_t Npoints;
  
  std::string simfile;
  std::string sizefile;
  
  TreeQuad * qtree;
};


#endif
