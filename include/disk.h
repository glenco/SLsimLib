#include "particle_halo.h"

#ifndef GLAMER_disk_halo_h
#define GLAMER_disk_halo_h

using Utilities::Geometry::Quaternion;

template <typename T=float>
class LensHaloDisk:public LensHaloParticles<ParticleType<T> > {
  
public:
  LensHaloDisk(
               double mass            /// mass of disk
               ,double disk_scale     /// scale hight of disk
               ,double Rperp          /// vertical scale hight of disk
               ,double mass_res       /// mass resolution, mass of particles
               ,float my_PA           /// position angle in radians
               ,float my_inclination  /// inclination of disk in radians, 0 is face on
               ,Utilities::RandomNumbers_NR &ran  /// random number generator
               ,float redshift
               ,const COSMOLOGY &cosmo
               ,int Nsmooth = 64      /// number of neighbors used in smoothing
               );
  
  LensHaloDisk(LensHaloDisk &&h):
  LensHaloParticles<ParticleType<T> >(std::move(h)),
  particles(LensHaloParticles<ParticleType<T> >::trash_collector)
  {
    qrot_invers = h.qrot_invers;  // rotation that brings the disk back to face on
    Rscale = h.Rscale;
    Rhight = h.Rhight;
    zpa = h.zpa;
    inclination = h.inclination;
  }
  
  LensHaloDisk & operator=(LensHaloDisk &&h){
    LensHaloParticles<ParticleType<T> >::operator=(std::move(h));
    
    qrot_invers = h.qrot_invers;  // rotation that brings the disk back to face on
    Rscale = h.Rscale;
    Rhight = h.Rhight;
    zpa = h.zpa;
    inclination = h.inclination;
    
    return *this;
  }

  
  /// Reorient the disk
  void reorient(float my_inclination,float my_PA){
    inclination = my_inclination;
    zpa = -my_PA;
    
    Quaternion R = Quaternion::q_z_rotation(zpa)*Quaternion::q_x_rotation(inclination)*qrot_invers;
    rotate_all(R);
    qrot_invers = Quaternion::q_x_rotation(-inclination)*Quaternion::q_z_rotation(-zpa);
    
    // reconstruct LensHaloParticles base class
    delete LensHaloParticles<ParticleType<T> >::qtree;
    LensHaloParticles<ParticleType<T> >::qtree = new TreeQuadParticles<ParticleType<T> >(particles.data(),particles.size(),-1,-1,0,20);
  }
  /// inclination in radians, 0 is face on
  float getInclination(){return inclination;}
  /// postion angle in radians,
  float getPA(){return -zpa;}

private:
  
  std::vector<ParticleType<T> > &particles;
  
  void rotate_all(Quaternion &R){
    Quaternion q(1,0,0,0);
    for(auto &p : particles){
      q[0] = 0; q[1] = p[0]; q[2] = p[1]; q[3] = p[2];
      q.RotInplace(R);
      p[0] = q[1]; p[1] = q[2]; p[2] = q[3];
    }
  }
  
  Utilities::Geometry::Quaternion qrot_invers;  // rotation that brings the disk back to face on
  double Rscale;
  double Rhight;
  double zpa;
  double inclination;
};

template<typename T>
LensHaloDisk<T>::LensHaloDisk(
              double mass
             ,double disk_scale
             ,double Rperp
             ,double mass_res
             ,float my_PA
             ,float my_inclination
             ,Utilities::RandomNumbers_NR &ran
             ,float redshift
             ,const COSMOLOGY &cosmo
             ,int Nsmooth
             ):
LensHaloParticles<ParticleType<T> >(redshift,cosmo),
particles(LensHaloParticles<ParticleType<T> >::trash_collector),
qrot_invers(1,0,0,0),Rscale(disk_scale),Rhight(Rperp),zpa(-my_PA),inclination(my_inclination)
{
  
  /// set up base Lenshalo
  size_t N = (size_t)(mass/mass_res + 1);
  particles.resize(N);
  LensHaloParticles<ParticleType<T> >::pp = particles.data();
  LensHaloParticles<ParticleType<T> >::Npoints = particles.size();

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
  Quaternion R = Quaternion::q_z_rotation(zpa)*Quaternion::q_x_rotation(inclination);
  rotate_all(R);
  qrot_invers = R.conj();
  
  LensHaloParticles<ParticleType<T> >::mcenter *= 0.0;
  LensHalo::setMass(mass);
  
  LensHaloParticles<ParticleType<T> >::min_size=0;
  LensHaloParticles<ParticleType<T> >::multimass=false;
  
  Point_2d no_rotation;
  LensHaloParticles<ParticleType<T> >::set_up(redshift,cosmo,no_rotation,false,false);
};

#endif
