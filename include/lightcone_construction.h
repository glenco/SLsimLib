//
//  lightcone_construction.h
//  GLAMER
//
//  Created by Ben Metcalf on 26/10/16.
//
//

#ifndef __GLAMER__lightcone_construction__
#define __GLAMER__lightcone_construction__

#include <stdio.h>
#include "point.h"
#include "lens_halos.h"
#include "geometry.h"


namespace LightCones{
  
  struct DataRockStar{
    unsigned int id;
    Point_3d x;
    double mass;
    double Rvir;     // in comoving units
    double Rscale;   // in comoving units
    double Vmax;
  };
  
  void ReadLightConeNFW(std::string filename,COSMOLOGY &cosmo
                        ,std::vector<LensHalo* > &lensVec
                        ,PosType &theta_max);
  
  void ReadLightConeParticles(std::string filename,COSMOLOGY &cosmo
                              ,std::vector<LensHaloParticles *> &lensVec
                              ,int Nplanes
                              ,float particle_mass
                              ,float particle_size
                              ,bool angular_sizes = false
                              ,bool verbose = false);
  
  void ReadLightConeParticles(std::string filename,COSMOLOGY &cosmo
                              ,std::vector<LensHaloMassMap *> &lensVec
                              ,int Nplanes
                              ,float particle_mass
                              ,float angular_resolution
                              ,bool verbose = false);
  
  /*** \brief Class for constructing light cones form simulation boxes.
   *
   *
   */
  class LightCone{
  public:
    LightCone(double angular_radius);
    
    void WriteLightCone(std::string filename,std::vector<DataRockStar> &vec);
    void WriteLightCone(std::string filename,std::vector<Point_3d> &vec);
    
    friend class MultiLightCone;
    
    template <typename T>
    void select(Point_3d xo,Point_3d v,double Length,double rlow,double rhigh
                ,T* begin
                ,T* end
                ,Utilities::LockableContainer<std::vector<T> > &incone
                ,bool periodic_boundaries = true);
    
  private:
    
    double r_theta; // angular radius of light cone (radians)
    double sin_theta_sqrt;
  };
  /// select the halos from the box that are within the light cone
  template <typename T>
  void LightCone::select(Point_3d xo,Point_3d v,double Length,double rlow,double rhigh
                         ,T* begin
                         ,T* end
                         ,Utilities::LockableContainer<std::vector<T> > &incone
                         ,bool periodic_boundaries){
    
    if(begin == end) return;
    
    Point_3d dx,x;
    double r2,xp;
    Point_3d n;
    double rhigh2 = rhigh*rhigh;
    double rlow2 = rlow*rlow;
    T b;
    
    v /= v.length();
    
    // standard coordinates so cone is always centered on the x-axis
    Point_3d y_axis,z_axis;
    
    y_axis[0] = -v[1];
    y_axis[1] =  v[0];
    y_axis[2] =  0;
    
    y_axis /= y_axis.length();
    z_axis = v.cross(y_axis);
    
    //std::cout << v*y_axis << " " << v*z_axis << " " << y_axis*z_axis << std::endl;
    
    // we can do better than this !!!!
    int n0[2] = {(int)((rhigh + xo[0])/Length),(int)((xo[0] - rhigh)/Length) - 1 };
    int n1[2] = {(int)((rhigh + xo[1])/Length),(int)((xo[1] - rhigh)/Length) - 1 };
    int n2[2] = {(int)((rhigh + xo[2])/Length),(int)((xo[2] - rhigh)/Length) - 1 };
    if(!periodic_boundaries){
      n0[0] = n0[1] = 0;
      n1[0] = n1[1] = 0;
      n2[0] = n2[1] = 0;
    }
    for(auto a = begin ; a != end ; ++a){
      dx[0] = (*a).x[0] - xo[0];
      dx[1] = (*a).x[1] - xo[1];
      dx[2] = (*a).x[2] - xo[2];
      
      int count = 0;
      
      for(n[0] = n0[1]; n[0] <= n0[0] ; ++n[0]){
        for(n[1] = n1[1]; n[1] <= n1[0] ; ++n[1]){
          for(n[2] = n2[1]; n[2] <= n2[0] ; ++n[2]){
            
            x = dx + n*Length;
            xp = x*v;
            if(xp > 0){
              r2 = x.length_sqr();
              
              if(r2 > rlow2 && r2 < rhigh2){
                if( r2*sin_theta_sqrt > r2 - xp*xp){
                  //std::cout << n << " | " << xo << " | "
                  //          << v << " |  " << x << std::endl;
                  
                  b = *a;
                  //b.x = x;
                  // rotate to standard reference frame
                  b.x[0] = xp;
                  b.x[1] = y_axis*x;
                  b.x[2] = z_axis*x;
                  
                  //std::cout << "b = " << b.x << std::endl;
                  
                  incone.push_back(b);
                  ++count;
                }
              }
            }
          }
        }
      }
      if(count > 1) std::cout << " Warning: Same halo appears " << count << " times in light-cone." << std::endl;
    }
  }
  
  class MultiLightCone{
  public:
    MultiLightCone(
                   double angular_radius
                   ,const std::vector<Point_3d> &observers    /// postion of observers within the simulation box
                   ,const std::vector<Point_3d> &directions   /// direction of light cones
    ):
    xos(observers),vs(directions)
    {
      assert(observers.size() == directions.size());
      for(int i=0;i<observers.size(); ++i) cones.push_back(angular_radius);
    }
    
    /** Read the points in from a snapshot of the halos created by RockStar.
     The output will be the halos in the cone in coordinates where the x-axis is
     the center of the cone.  Changing the observer, xo, and direction of view, V,
     effectively translates and rotates the box, respectively.  Only halos with radii
     between rlow and rhigh are added.
     */
    void ReadBoxRockStar(std::string filename
                         ,double rlow,double rhigh
                         ,std::vector<std::vector<LightCones::DataRockStar> > &conehalos
                         ,bool periodic_boundaries = true
                         ,bool allow_subhalos = false);
    
    void ReadBoxXYZ(std::string filename
                    ,double rlow,double rhigh
                    ,std::vector<std::vector<Point_3d> > &conehalos
                    ,double hubble
                    ,double BoxLength
                    ,bool periodic_boundaries = true
                    );
    
    void ReadBoxToMaps(std::string filename
                       ,double rlow,double rhigh
                       ,std::vector<std::vector<PixelMap>  > &maps
                       ,double hubble
                       ,double BoxLength
                       ,bool periodic_boundaries = true
                       );
    
    
  private:
    
    std::vector<Point_3d> xos;
    std::vector<Point_3d> vs;
    std::vector<LightCone> cones;
  };
  
  /** \brief
   
   The Born approximation is used.
   */
  void FastLightCones(
                      const COSMOLOGY &cosmo
                      ,const std::vector<double> &zsources    // vector of source redshifts desired
                      ,std::vector<std::vector<PixelMap> > &maps
                      ,double range
                      ,double angular_resolution
                      ,std::vector<Point_3d> &observers    /// position of observers within the simulation box
                      ,std::vector<Point_3d> &directions   /// direction of light cones
                      ,const std::vector<std::string> &snap_filenames
                      ,const std::vector<float> &snap_redshifts
                      ,double BoxLength
                      ,double particle_mass
                      ,bool verbose = false
                      ,bool addtocone = false  /// if false the maps will be cleared and new maps made, if false particles are added to the existing maps
  );

  using Utilities::Geometry::Quaternion;
  using Utilities::Geometry::SphericalPoint;

  void _fastplanes_parallel_(Point_3d *begin,Point_3d *end
                 ,const COSMOLOGY &cosmo
                 ,std::vector<std::vector<Point_3d> > &boxes
                 ,std::vector<Point_3d> &observers
                 ,std::vector<Quaternion> &rotationQs
                 ,std::vector<double> &dsources
                 ,std::vector<std::vector<PixelMap> > &maps
                 ,double dmin
                 ,double dmax
                 ,double BoxLength
                 );

  
  /** \brief class for generating positions in proportion to mass in an NFW profiles
   */
  class NFWgenerator{
  public:
    NFWgenerator(Utilities::RandomNumbers_NR &ran_in,double max_cons)
    :ran(ran_in)
    {
      X.resize(N);
      F.resize(N);
      X[0] = F[0] = 0.0;
      dx = max_cons/(N-1);
      for(int i=1;i<N;++i){
        X[i] = i*dx;
        F[i] = log(1+X[i]) - X[i]/(1+X[i]) ;
      }
    }
    
    /// returns a vector of points drawn from a spherical halo
    void drawSpherical(std::vector<Point_3d> &points  /// output points
                       ,double cons                   /// concentration
                       ,double Rvir                   /// maximum elliptical radius
                       ){
      double Fmax = log(1+cons) - cons/(1+cons);
      double rs = Rvir/cons;
      for(auto p : points){
        double f = Fmax*ran();
        size_t i = Utilities::locate(F,f);
        double x = X[i] + dx*(f - F[i])/(F[i+1] - F[i]);
        
        p[0] = ran.gauss();
        p[1] = ran.gauss();
        p[2] = ran.gauss();
        
        p *= x*rs/p.length();
      }
    }
    
    ///  STILL UNDER CONSTRUCTION returns a vector of points drawn from a triaxial halo,
    void drawTriAxial(std::vector<Point_3d> &points  /// output points
                      ,double cons                   /// concentration
                      ,double Rvir                   /// maximum elliptical radius
                      ,double f1                     /// axis ratio 1 to 3
                      ,double f2                     /// axis ratio 2 to 3
                      ,SphericalPoint v              /// direction of axis 3
                      ){
      
      double a3 = 1.0/pow(f1*f2,1.0/3.0);
      double a1 = f1*a3,a2 = f2*a3;
      Quaternion rot = Quaternion::q_z_rotation( v.phi )*Quaternion::q_y_rotation( v.theta );
      
      double Fmax = log(1+cons) - cons/(1+cons);
      double rs = Rvir/cons;
      for(auto p : points){

        p[0] = ran.gauss();
        p[1] = ran.gauss();
        p[2] = ran.gauss();

        double f = Fmax*ran();
        size_t i = Utilities::locate(F,f);
        double x = X[i] + dx*(f - F[i])/(F[i+1] - F[i]);
        
        p *= x*rs/p.length();

        p[0] *= a1;
        p[1] *= a2;
        p[2] *= a3;
        
        p = rot.Rotate(p);
      }
    }
    
  private:
    Utilities::RandomNumbers_NR &ran;
    double dx;
    std::vector<double> X;
    std::vector<double> F;
    const int N = 1000;
  };
  
  
}
#endif /* defined(__GLAMER__lightcone_construction__) */
