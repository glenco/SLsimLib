#include "particle_halo.h"

#ifndef GLAMER_disk_halo_h
#define GLAMER_disk_halo_h

using Utilities::Geometry::Quaternion;

template <typename T=float>
class LensHaloDisk:public LensHaloParticles<ParticleType<T> > {
  
public:
  LensHaloDisk(double mass
               ,double disk_scale       /// scale length in Mpc
               ,double Rperp            /// vertical scale length in Mpc
               ,double mass_res
               ,float my_inclination    /// disk inclination, 0 is face on (radians)
               ,float my_PA             /// position angle (radians)
               ,Utilities::RandomNumbers_NR &ran
               ,float redshift
               ,const COSMOLOGY &cosmo
               ,int Nsmooth = 64
               ):
  LensHaloParticles<ParticleType<T> >(redshift,cosmo),qrot_invers(1,0,0,0),
  Rscale(disk_scale),Rhight(Rperp),PA(my_PA),inclination(my_inclination)
  {
    
    /// set up base Lenshalo
    size_t N = (size_t)(mass/mass_res + 1);
    particles.resize(N);
    
    double r,theta;
    
    size_t i = 0;
    for(auto &p : particles){
      r = -Rscale * log(1 - (float)(i) / N );
      theta = 2*PI*ran();
      
      p.x[0] = r*cos(theta);
      p.x[1] = r*sin(theta);
      p.x[2] = -Rhight * log( 1 - ran() );
      p.x[2] *= 2*(int)(2*ran()) - 1;  // random sign
      
      p.Mass = mass_res;
      ++i;
    }
    
    LensHalo::Rmax = -3 * Rscale * log(1 - (float)(N-1) / N );
    LensHalo::setRsize( LensHalo::Rmax );
    
    LensHaloParticles<ParticleType<T> >::calculate_smoothing(Nsmooth,particles.data(),particles.size());
    
    // rotate particles to required inclination and position angle
    Quaternion R = Quaternion::q_z_rotation(PA)*Quaternion::q_x_rotation(inclination);
    rotate_all(R);
    qrot_invers = R.conj();
    
    LensHaloParticles<ParticleType<T> >::mcenter *= 0.0;
    LensHalo::setMass(mass);
    
    LensHaloParticles<ParticleType<T> >::qtree = new TreeQuadParticles<ParticleType<T> >(particles.data(),particles.size(),-1,-1,0,20);
  }
  
  LensHaloDisk(LensHaloDisk &&h):LensHaloParticles<ParticleType<T> >(std::move(h)){
    particles=std::move(h.particles);
    
    //??? This and the move assignment operator depends on the move operator
    //  for the std::vector keeping the pointers valid.  I think this is always
    //  the case but it is not required by the C++ standard.
    
    qrot_invers = h.qrot_invers;  // rotation that brings the disk back to face on
    Rscale = h.Rscale;
    Rhight = h.Rhight;
    PA = h.PA;
    inclination = h.inclination;
  }
  
  LensHaloDisk & operator=(LensHaloDisk &&h){
    LensHaloParticles<ParticleType<T> >::operator=(std::move(h));
    particles=std::move(h.particles);
    
    qrot_invers = h.qrot_invers;  // rotation that brings the disk back to face on
    Rscale = h.Rscale;
    Rhight = h.Rhight;
    PA = h.PA;
    inclination = h.inclination;
    
    return *this;
  }

  
  /// Reorient the disk
  void reorient(float my_inclination,float my_PA){
    inclination = my_inclination;
    PA=my_PA;
    
    Quaternion R = Quaternion::q_z_rotation(PA)*Quaternion::q_x_rotation(inclination)*qrot_invers;
    rotate_all(R);
    qrot_invers = Quaternion::q_x_rotation(-inclination)*Quaternion::q_z_rotation(-PA);
    
    // reconstruct LensHaloParticles base class
    delete LensHaloParticles<ParticleType<T> >::qtree;
    LensHaloParticles<ParticleType<T> >::qtree = new TreeQuadParticles<ParticleType<T> >(particles.data(),particles.size(),-1,-1,0,20);
  }
  
  
private:
  
  void rotate_all(Quaternion &R){
    Quaternion q(1,0,0,0);
    for(auto &p : particles){
      q[0] = 0; q[1] = p[0]; q[2] = p[1]; q[3] = p[2];
      q.Rotate(R);
      p[0] = q[1]; p[1] = q[2]; p[2] = q[3];
    }
  }
  
  
  Utilities::Geometry::Quaternion qrot_invers;  // rotation that brings the disk back to face on
  double Rscale;
  double Rhight;
  double PA;
  double inclination;
  
  std::vector<ParticleType<T> > particles;
};

#endif
