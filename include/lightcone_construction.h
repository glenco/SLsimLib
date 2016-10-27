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


/*** \brief Class for constructing light cones form boxes
 */
class LightCone{
  LightCone(Point_3d x_observer,Point_3d v_direction,double angular_radius):
  xo(x_observer),v(v_direction),r_theta(angular_radius)
    {
      v = v/v.length();
    };
  
  
  struct DataRockStar{
    double Vmax;
    Point_3d x;
  };

  void ReadBoxRockStar(std::string filename,std::vector<DataRockStar> &output);
  void PrintConeRockStar(std::string outfile);
  
  
  
  /// select the halos from the box that are within the light cone
  template <typename T>
  void select(double Length,double rlow,double rhigh,std::vector<T> &input
              ,std::vector<T> &incone){
    Point_3d dx,x;
    double vx,r2,xp;
    Point_3d n;
    double tmp = 1+r_theta*r_theta;
    double rhigh2 = rhigh*rhigh;
    double rlow2 = rlow*rlow;
    T b;
    
    // we can do better than this !!!!
    int n0 = (int)(rhigh/Length);
    //int n1 = (int)(v[1]*rhigh/Length);
    //int n2 = (int)(v[2]*rhigh/Length);
    
    for(auto a : input){
      dx = a.x - xo;
      vx = v*dx;
      for(n[0] = -n0; n[0]<=n0 ; ++n[0]){
        for(n[1] = -n0; n[1]<=n0 ; ++n[1]){
          for(n[2] = -n0; n[2]<n0 ; ++n[2]){
            
            x = dx + n*Length;
            r2 = x.length_sqr();
            xp = vx + (v*n)*Length;
            
            if( r2 > tmp*xp*xp ){
              
              if(r2 < rlow2 && r2 < rhigh2){
                b = a;
                b.x = x;
                incone.push_back(b);
              }
            }
          }
        }
      }
    }
  }
  
  
private:
  Point_3d xo; // position of observer
  Point_3d v; // direction of center of cone
  double r_theta; // angular redius of light cone (radians)
};

#endif /* defined(__GLAMER__lightcone_construction__) */
