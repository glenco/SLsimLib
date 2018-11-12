//
//  particle_halo.h
//  GLAMER
//
//  Created by bmetcalf on 16/06/15.
//
//

things to do:

protect particle data by storing it in a class that reads data and generates LensHaloParticles

make multi-type halos center on a common cneter of mass

#ifndef GLAMER_particle_halo_h
#define GLAMER_particle_halo_h

#include "geometry.h"
#include "quadTree.h"
#include "simpleTree.h"
#include "particle_types.h"
#include "gadget.hh"

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
 
  LensHaloParticles(PType  *pdata
                    ,float redshift        /// redshift of origin
                    ,const COSMOLOGY& cosmo  /// cosmology
                    ,Point_2d theta_rotate   /// rotation of particles around the origin
                    ,bool recenter           /// center on center of mass
                    ,float MinPSize        /// minimum particle size
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
                                  ,PType *pp
                                  ,size_t Npoints);
  
  static void readPositionFileASCII(const std::string &filename
                                    ,bool multimass
                                    ,PType *pp
                                    );
  
  
  
  static void writeSizes(const std::string &filename
                         ,int Nsmooth
                         ,const PType *pp
                         ,size_t Npoints
                         );

  static bool readSizesFile(const std::string &filename
                            ,PType * pp
                            ,size_t Npoints
                            ,int Nsmooth
                            ,PosType min_size);

private:

  Point_3d mcenter;
  void rotate_particles(PosType theta_x,PosType theta_y);

  static void smooth_(TreeSimple<PType> *tree3d,PType *xp,size_t N,int Nsmooth);
  
  void assignParams(InputParams& params);

  PType *pp;
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

void makeHalosFromGadget(const std::string &filename
                                ,std::vector<ParticleType<float> > &data
                                ,std::vector<LensHaloParticles<ParticleType<float> > > &halos
                                ,int Nsmooth
                                ,float redshift
                                ,Point_2d theta_rotate   /// rotation of particles around the origin
                                ,const COSMOLOGY &cosmo
                                ){
  
  
  GadgetFile<ParticleType<float> > gadget_file(filename,data);
  
  for(int n=0 ; n < gadget_file.numfiles ; ++n){
    gadget_file.openFile();
    gadget_file.readBlock("POS");
    gadget_file.readBlock("MASS");
    gadget_file.closeFile();
  }
  
  // sort by type
  std::sort(data.begin(),data.end(),[](ParticleType<float> &a1,ParticleType<float> &a2){return a1.type < a2.type;});
  
  ParticleType<float> *pp;
  
  size_t skip = 0;
  std::string sizefile;
  for(int i = 0 ; i < 6 ; ++i){  //loop through type
    if(gadget_file.npart[i] > 0){
      
      pp = data.data() + skip;  // pointer to first particle of type
      size_t N = gadget_file.npart[i];
      
      sizefile = filename + "_T" + std::to_string(i) + "S"
      + std::to_string(Nsmooth) + "sizes";
      
      if(!LensHaloParticles<ParticleType<float> >::readSizesFile(sizefile,pp,gadget_file.npart[i],Nsmooth,0)){
        // calculate sizes
        //sizes.resize(Npoints);
        LensHaloParticles<ParticleType<float> >::calculate_smoothing(Nsmooth,pp,N);
        
        // save result to a file for future use
        LensHaloParticles<ParticleType<float> >::writeSizes(sizefile,Nsmooth,pp,N);
      }
      
      halos.emplace_back(pp
                         ,redshift
                         ,cosmo  /// cosmology
                         ,theta_rotate   /// rotation of particles around the origin
                         ,false
                         ,0);
      
    }
    skip += gadget_file.npart[i];
  }
};


#endif
