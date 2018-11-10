//
//  particle_halo.h
//  GLAMER
//
//  Created by bmetcalf on 16/06/15.
//
//

#ifndef GLAMER_particle_halo_h
#define GLAMER_particle_halo_h

#include "geometry.h"
#include "quadTree.h"
#include "simpleTree.h"

/// Atomic data class for simulation particles with individual sizes and masses
template<typename T = float>
struct ParticleData{
  T &operator[](int i){return x[i];}
  T *operator*(){return x;}
  T x[3];
  float Mass;
  float Size;
  int type;
  
  float size(){return Size;}
  float mass(){return Mass;}
};

/// Atomic data class for simulation particles of the same size and mass
struct ParticleDataSimple{
  float &operator[](int i){return x[i];}
  float *operator*(){return x;}
  float x[3];
  
  static float Mass;
  static float Size;

  float size(){return ParticleDataSimple::Size;}
  float mass(){return ParticleDataSimple::Mass;}
};
float ParticleDataSimple::Size = 0;
float ParticleDataSimple::Mass = 0;

enum SimFileFormats {ascii,gadget2};

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

template<typename PType>
class LensHaloParticles : public LensHalo
{
public:
  
  
  LensHaloParticles(const std::string& simulation_filename /// name of data files
                    ,SimFileFormats format   /// format of data file
                    ,PosType redshift        /// redshift of origin
                    ,int Nsmooth             /// number of neighbours for adaptive smoothing
                    ,const COSMOLOGY& cosmo  /// cosmology
                    ,Point_2d theta_rotate   /// rotation of particles around the origin
                    ,bool recenter           /// center on center of mass
                    ,bool my_multimass       /// set to true is particles have different sizes
                    ,PosType MinPSize        /// minimum particle size
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
  
  static void calculate_smoothing(int Nsmooth
                                  ,std::vector<PType> &xxp
                                  );
  
  static void readPositionFileASCII(const std::string &filename
                                    ,bool multimass
                                    ,std::vector<PType> &xxp
                                    );

  static void writeSizes(const std::string &filename
                         ,int Nsmooth
                         ,const std::vector<PType> &xxp
                         );

private:

  
  Point_3d mcenter;
  void rotate_particles(PosType theta_x,PosType theta_y);

  static void smooth_(TreeSimple<PType> *tree3d,PType *xp,size_t N,int Nsmooth);
  
  bool readSizesFile(const std::string& filename,int Nsmooth,PosType min_size);
  
  void assignParams(InputParams& params);

  std::vector<PType> xxp;
  //PosType **xp;
  //std::vector<float> masses;
  //std::vector<float> sizes;
  
  PosType min_size;
  bool multimass;
  
  Utilities::Geometry::SphericalPoint center;
  
  size_t Npoints;
  
  std::string simfile;
  std::string sizefile;
  
  TreeQuadParticles<PType> * qtree;
};


#endif
