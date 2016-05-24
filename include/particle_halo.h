//
//  particle_halo.h
//  GLAMER
//
//  Created by bmetcalf on 16/06/15.
//
//
#include "quadTree.h"

#ifndef GLAMER_particle_halo_h
#define GLAMER_particle_halo_h

/**
 *  \brief A class that represents the lensing by a collection of simulation particles.
 *
 *  Smoothing is done according to the density of particles in 3D.  Smoothing sizes are
 *  either read in from a file (names simulation_filename + "." + Nsmooth + "sizes") or calculated 
 *  if the file does not exist (in which case the file is created).  This can be
 *  time and memory consuming when there are a large number of particles.
 *
 *  Input format:
 *    ASCCI - a table of three floats for positions in comoving Mpc
 *            The lines "# nparticles ...." and "# mass ...." must be in the
 *            near the top of the file. # is otherwise a comment character.  
 *            Only one type of particle in a single input file.
 *
 *   More input formats will be added in the future. 
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
                    );
  
  ~LensHaloParticles();
  
  void force_halo(double *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,double const *xcm
                  ,bool subtract_point=false,PosType screening = 1.0);

  size_t getN() const { return Npoints; };
  
  /// rotate the simulation, radians
  void rotate(Point_2d theta);
  
  /// center of mass in input coordinates
  Point_3d CenterOfMass(){return mcenter;}
  
private:

  Point_3d mcenter;
  void rotate_particles(PosType theta_x,PosType theta_y);

  void calculate_smoothing(int Nsmooth);
  void smooth_(TreeSimple *tree3d,PosType **xp,float *sizes,size_t N,int Nsmooth);

  void readPositionFileASCII(const std::string& filename);
  bool readSizesFile(const std::string& filename,int Nsmooth);
  void writeSizes(const std::string& filename,int Nsmooth);
  
  void assignParams(InputParams& params);

  PosType **xp;
  std::vector<float> masses;
  std::vector<float> sizes;
  bool multimass;
  
  Utilities::Geometry::SphericalPoint center;
  
  size_t Npoints;
  
  std::string simfile;
  std::string sizefile;
  
  TreeQuad * qtree;
};


#endif
