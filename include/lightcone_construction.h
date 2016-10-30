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


/*** \brief Class for constructing light cones form boxes
 */
class LightCone{
public:
  LightCone(double angular_radius);
  
  struct DataRockStar{
    unsigned int id;
    Point_3d x;
    double mass;
    double Rvir;     // in comoving units
    double Rscale;   // in comoving units
  };

  void ReadBoxRockStar(std::string filename,Point_3d xo,Point_3d V
                                  ,double rlow,double rhigh
                                  ,std::vector<DataRockStar> &conehalos);
  
  static void ReadLightCone(std::string filename,COSMOLOGY &cosmo
                     ,std::vector<LensHalo* > &lensVec);
  
  static void WriteLightCone(std::string filename,std::vector<DataRockStar> &vec);
  
  /// select the halos from the box that are within the light cone
  template <typename T>
  void select(Point_3d xo,Point_3d v,double Length,double rlow,double rhigh
              ,std::vector<T> &input,std::vector<T> &incone){
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

    std::cout << v*y_axis << " " << v*z_axis << " " << y_axis*z_axis << std::endl;

    // we can do better than this !!!!
    int n0[2] = {(int)((rhigh + xo[0])/Length),(int)((xo[0] - rhigh)/Length) - 1 };
    int n1[2] = {(int)((rhigh + xo[1])/Length),(int)((xo[1] - rhigh)/Length) - 1 };
    int n2[2] = {(int)((rhigh + xo[2])/Length),(int)((xo[2] - rhigh)/Length) - 1 };
    
    for(auto a : input){
      dx = a.x - xo;
      
      for(n[0] = n0[1]; n[0] <= n0[0] ; ++n[0]){
        for(n[1] = n1[1]; n[1] <= n1[0] ; ++n[1]){
          for(n[2] = n2[1]; n[2] <= n2[0] ; ++n[2]){
            
            x = dx + n*Length;
            xp = x*v;
            if(xp > 0){
              r2 = x.length_sqr();
              
              if( r2*sin_theta_sqrt > r2 - xp*xp){
                if(r2 > rlow2 && r2 < rhigh2){
                  std::cout << n << " | " << xo << " | "
                            << v << " |  " << x << std::endl;
                  
                  b = a;
                  b.x = x;
                  // rotate to standard reference frame
                  b.x[0] = xp;
                  b.x[1] = y_axis*x;
                  b.x[2] = z_axis*x;
                  
                  std::cout << "b = " << b.x << std::endl;
                  
                  incone.push_back(b);
                }
              }
            }
          }
        }
      }
    }
  }
  
  
private:
  //Point_3d xo; // position of observer
  //Point_3d v; // direction of center of cone
  double r_theta; // angular radius of light cone (radians)
  double sin_theta_sqrt;
};

#endif /* defined(__GLAMER__lightcone_construction__) */
