//
//  particle_types.h
//  GLAMER
//
//  Created by Ben Metcalf on 12/11/2018.
//

#ifndef particle_types_h
#define particle_types_h

enum SimFileFormat {glmb,csv3,csv4,csv5,csv6,gadget2,ascii};

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
//float ParticleTypeSimple::Size = 0;
//float ParticleTypeSimple::Mass = 0;

/// Atomic data class for stars with different masses
struct StarType{
  double &operator[](int i){return x[i];}
  double *operator*(){return x;}
  double x[3];
  float Mass;
  
  float mass() const {return Mass;}
  static float size() {return 0;}
};


#endif /* particle_types_h */
