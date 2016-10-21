//
//  geometry.h
//  GLAMER
//
//  Created by bmetcalf on 7/24/14.
//
//

#ifndef __GLAMER__geometry__
#define __GLAMER__geometry__

#include <iostream>
#include "standard.h"
#include "point.h"
//
//  geometry.h
//  GLAMER
//
//  Created by bmetcalf on 7/23/14.
//
//
namespace Utilities {
  /// Namespace for geometrical functions mostly having to do with spherical coordinates
  namespace Geometry{
    
    /// represents a point in spherical coordinates theta = 0 is equator
    class SphericalPoint{
      
    public:
      SphericalPoint(PosType r,PosType theta,PosType phi):r(r),theta(theta),phi(phi){};
      SphericalPoint():r(0),theta(0),phi(0){};
      
      PosType r;
      PosType theta;
      PosType phi;
      
      /// output Cartesian coordinates of the point
      void sphericalTOcartisian(PosType x[]) const;
      void cartisianTOspherical(PosType const x[]);
      void StereographicProjection(const SphericalPoint &central,PosType x[]) const;
      void OrthographicProjection(const SphericalPoint &central,PosType x[]) const;
      void InverseOrthographicProjection(const SphericalPoint &central,PosType const x[]);
    };
    
    ///  3 dimensional distance between points
    PosType Seporation(const SphericalPoint &p1,const SphericalPoint &p2);
    ///  Angular seporation between points
    PosType AngleSeporation(const SphericalPoint &p1,const SphericalPoint &p2);

    /// Determine if line segments a1a2 and b1b2 intersect.  Sharing an endpoint does not count as intersecting
    bool intersect(const PosType a1[],const PosType a2[],const PosType b1[],const PosType b2[]);
    /// returns the number of times a closed curve intersects itself
    int intersect(const std::vector<Point_2d> &curve);
    
    /** \brief To find orientation of the triangle formed by the ordered triplet (p, q, r).
     
     The function returns following values
     0 --> p, q and r are collinear
     1 --> Clockwise
     2 --> Counterclockwise
     */
    int orientation(const PosType p[],const PosType q[],const PosType r[]);
    /** \brief Given three collinear points p, q, r, the function checks if
     point q lies on line segment 'pr', but not at p or r
     */
    bool onSegment(const PosType p[], const PosType q[], const PosType r[]);
    
    /// returns the angle between two 2 dimensional vectors in radians.
    double AngleBetween2d(double v1[],double v2[]);
    /**   \brief returns 1 if it is in the curve and 0 if it is out.  Borders count as in.
    *
     *  This is faster than the windings() functions which also calulate the area
     */
    int incurve(PosType x[],std::vector<double *> curve);
    
    /// K=3
    class HealPixPoint{
    public:
      HealPixPoint(int h=4):H(h)
      {
        theta_c = asin( 2./3. );
      }
      HealPixPoint(const SphericalPoint &sp,int h=4):H(h)
      {
        theta_c = asin( 2./3. );
        hp2sp(sp);
      }
      HealPixPoint(double x,double y,int h=4):H(h)
      {
        theta_c = asin( 2./3. );
        xx[0] = x;
        xx[1] = y;
      }
    
      /// set a HealPix point coordinates from sphereical coordinates
      void hp2sp(const SphericalPoint &sp){
        
        if(fabs(sp.theta) < theta_c){
          xx[0] = sp.phi;
          xx[1] = 1.5*pi*sin(sp.theta)/H;
        }else{
          
          double sigma = sqrt(3*(1 - fabs(sin(sp.theta)) ) );
          xx[0] = std::copysign(1.0,sp.theta)*pi*(2-sigma)/H;
          
          double phi_c = -pi + ( 2*std::floor( (sp.phi + pi)*H/2/pi ) + 1 )*pi/H;

          xx[1] = phi_c + (sp.phi - phi_c)*sigma;
        }
      }
      
      /// output HealPix point in spherical coordinates
      SphericalPoint sp2hp(){
        SphericalPoint sp;
        if(fabs(xx[1]) < pi/H){
          sp.phi = xx[0];
          sp.theta = asin(xx[1]*H/pi/1.5);
        }else{
          double sigma = 2 - fabs(xx[1]*H)/pi;
          double xc = -pi +(2*std::floor((xx[0]+pi)*H/2/pi) +1)*pi/H;
          sp.phi = xc + (xx[0] -xc)/sigma;
          sp.theta = std::copysign(1,xx[1])*asin(1-sigma*sigma/3);
        }
        
        return sp;
      }
      
      double x(){return xx[0];}
      double y(){return xx[1];}
      
    private:
      Point_2d xx;
      double theta_c;
      int H;
    
    };
}
}

#endif /* defined(__GLAMER__geometry__) */
