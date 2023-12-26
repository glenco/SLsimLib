//
//  particle_types.h
//  GLAMER
//
//  Created by Ben Metcalf on 12/11/2018.
//
/*! \file particle_types.h

 */
#ifndef particle_types_h
#define particle_types_h

#include "point.h"

/** Formats for particle data files

 For csv3,csv4,csv5,csv6 - CSV ascii format without header.  The first three columns
 are the positions.  Next columns are used for the other formats being and
 interpreted as (column 4) masses are in Msun/h, (column 5) the paricle smoothing
 size in Mpc/h and (column 6) an integer for type of particle.  There can be more
 columns in the file than are uesed.  In the case of csv6, when there are more then one
 type of halo each type will be in a differeent LensHaloParticles with differnt smoothing.

 ascii - This is the original ascii GLAMER format.
 Three floats for positions in comoving Mpc (no h factor).
 The lines "# nparticles ...." and "# mass ...." must be
 in header at the top of the file. # is otherwise a
 comment character.  Only one type of particle in a single
 input file.
 
 */
enum class SimFileFormat {
  glmb     /**This is a binary format internal to GLAMER used to store
            the positions, masses and sizes of the particles.  If
            GLAMER has generated one, it should be all that is
            needed to recreate the LensHaloParticles.*/
  ,csv3    /// see above
  ,csv4    /// see above
  ,csv5    /// see above
  ,csv6    /// see above
  ,gadget2 /// Gadget 2 output file format
  ,ascii   /// the original ascii GLAMER format.
};

// Atomic data class for simulation particles with individual sizes and masses
template<typename T = float>
struct ParticleType{
  T &operator[](int i){return x[i];}
  T *operator*(){return x;}
  T x[3];
  
  float Mass;
  float Size;
  int type;
  
  float size(){return Size;}
  float mass(){return Mass;}
};

// Atomic data class for simulation particles with individual sizes, masses and velocities
template<typename T = float>
struct ParticleTypeV{
  T &operator[](int i){return x[i];}
  T *operator*(){return x;}
  T x[3];

  T v[3];

  float Mass;
  float Size;
  int type;
  
  float size(){return Size;}
  float mass(){return Mass;}
};

/// Atomic data class for simulation particles of the same size and mass
struct ParticleTypeSimple{
  float &operator[](int i){return x[i];}
  float *operator*(){return x;}
  float x[3];
  
  static float Mass;
  static float Size;
  
  float size(){return ParticleTypeSimple::Size;}
  float mass(){return ParticleTypeSimple::Mass;}
};

/// Atomic data class for simulation particles of the same size and mass
struct ParticleType2D{
  float &operator[](int i){return x[i];}
  float *operator*(){return x;}
  float x[2];
  
  float Mass;
  float Size;
  
  float size(){return Size;}
  float mass(){return Mass;}
};

/// Atomic data class for stars with different masses
struct StarType{
  double &operator[](int i){return x[i];}
  Point_3d<> operator*(){return x;}
  Point_3d<> x;
  float Mass;
  
  float mass() const {return Mass;}
  static float size() {return 0;}
};


#endif /* particle_types_h */
