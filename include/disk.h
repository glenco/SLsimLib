#include "particle_halo.h"

#ifndef GLAMER_disk_halo_h
#define GLAMER_disk_halo_h

using Utilities::Geometry::Quaternion;

/** \brief Creates a exponential disk out of particles.
 
 The disk is created out of particles and the smoothing done by nearest-N neighbour
 B-spline smoothing as if they came from a simulation, but they are placed more regularly
 so that the surface density is relatively smooth.
 */

template <typename T=float>
class LensHaloDisk:public LensHaloParticles<ParticleType<T> > {
  
public:
  LensHaloDisk(
               double mass            /// mass of disk
               ,double disk_scale     /// scale hight of disk (Mpc)
               ,double Rperp          /// vertical scale hight of disk (Mpc)
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
    assert(&h != this);
    
    LensHaloParticles<ParticleType<T> >::operator=(std::move(h));
    
    particles = LensHaloParticles<ParticleType<T> >::trash_collector;
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
    zpa = my_PA;
    
    Quaternion<> R = Quaternion<>::q_z_rotation(zpa)*Quaternion<>::q_x_rotation(inclination)*qrot_invers;
    rotate_all(R);
    qrot_invers = Quaternion<>::q_x_rotation(-inclination)*Quaternion<>::q_z_rotation(-zpa);
    
    // reconstruct LensHaloParticles base class
    delete LensHaloParticles<ParticleType<T> >::qtree;
    LensHaloParticles<ParticleType<T> >::qtree = new TreeQuadParticles<ParticleType<T> >(particles.data(),particles.size(),-1,-1,0,20);
  }
  /// inclination in radians, 0 is face on
  float getInclination(){return inclination;}
  /// postion angle in radians,
  float getPA(){return zpa;}

private:
  
  std::vector<ParticleType<T> > &particles;
  
  void rotate_all(Quaternion<T> &R){
    Quaternion<T> q(1,0,0,0);
    for(auto &p : particles){
      q[0] = 0; q[1] = p[0]; q[2] = p[1]; q[3] = p[2];
      q.RotInplace(R);
      p[0] = q[1]; p[1] = q[2]; p[2] = q[3];
    }
  }
  
  Utilities::Geometry::Quaternion<T> qrot_invers;  // rotation that brings the disk back to face on
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
qrot_invers(1,0,0,0),Rscale(disk_scale),Rhight(Rperp),zpa(my_PA),
inclination(my_inclination)
{
  
  /// set up base Lenshalo
  size_t N = (size_t)(mass/mass_res + 1);
  particles.resize(N);
  LensHaloParticles<ParticleType<T> >::pp = particles.data();
  LensHaloParticles<ParticleType<T> >::Npoints = particles.size();

  double r,theta=0;
  double dt = PI +  PI*ran()/10.;

  //std::cout << "dt = " << dt/PI << std::endl;
  size_t i = 0;
  double deltaF = 1.0/(N-1);
  double x = 0;
  for(auto &p : particles){
    r = Rscale * x;
    
    // quadratic recursive approximation
    //x = x + 0.5*( sqrt(x + 4*(1-x)*exp(x)*deltaF ) - x)
    ///(1-x);
    
    if(x==0){
      x = sqrt(deltaF);
    }else{
      x = x + deltaF*exp(x)/x;
    }
    
    assert(!isnan(x));
    //r = -Rscale * log(1 - (float)(i) / N );
    //theta = 2*PI*ran();
    theta += dt;
      
    p.x[0] = r*cos(theta);
    p.x[1] = r*sin(theta);
    if(Rhight > 0){
      p.x[2] = -Rhight * log( 1 - ran() );
      p.x[2] *= 2*(int)(2*ran()) - 1;  // random sign
    }else{
      p.x[2] = 0;
    }
    p.Mass = mass_res;
    ++i;
  }
  
 
  LensHaloParticles<ParticleType<T> >::calculate_smoothing(Nsmooth,particles.data(),particles.size());
  
  // rotate particles to required inclination and position angle
  Quaternion<T> R = Quaternion<T>::q_z_rotation(zpa) * Quaternion<T>::q_y_rotation(inclination);
  rotate_all(R);
  qrot_invers = R.conj();
  
  LensHaloParticles<ParticleType<T> >::mcenter *= 0.0;
  LensHalo::setMass(mass);
  
  //LensHaloParticles<ParticleType<T> >::min_size;
  LensHaloParticles<ParticleType<T> >::multimass=false;
  
  Point_2d no_rotation;
  LensHaloParticles<ParticleType<T> >::set_up(redshift,cosmo,no_rotation,-1,false,false);

  LensHalo::Rmax = -8 * Rscale * log(1 - (float)(N-1) / N );
  LensHalo::setRsize( LensHalo::Rmax );
  
};

#endif
