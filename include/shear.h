//
//  shear.h
//  GLAMER
//
//  Created by Robert Benton Metcalf on 24/06/2020.
//
// These data structures and functions for storing and
// minupulating shear and correlations between shears
// on the sphere and flat-sky

#ifndef shear_h
#define shear_h

#include "geometry.h"

/// this data class represents a postion inspherical coordinates and a polarization relative to the sphircal coordinate system
template <typename T>
struct Polar : public Utilities::Geometry::SphericalPoint<T>
{
  Point_2d shear;
  T kappa;
};

/// finds the tangential and x components of polarization realtive to the great circle seporating the points of the sphere represent by x0 and p.  This done in spherical coordinates.
template <typename T>
Point_2d tangent_cross_shear(
                             const Utilities::Geometry::SphericalPoint<T> &xo
                             ,const Polar<T> &p
                             ){
  //double d = xo.angular_separation(p);
  //double cd = cos(d);
  //double sd = sin(d);
 
  Point_3d<T> v1 = xo.TOcartisian() /( (xo.r <= 0) ? 1 : xo.r );
  Point_3d<T> v2 = p.TOcartisian() /( (p.r <= 0) ? 1 : p.r );
  double cd = v1*v2;
  
  // tangent vectors to the great circle at the two points
  Point_3d<T> t2 = (v2 * cd - v1);
  //Point_3d<T> t2 = ( v2 - v1*cd );
  t2.unitize();
  
  // unit vector in theta direction
  //Point_3d<T> v_theta_2 = p.theta_hat();
  Point_3d<T> v_theta_2 = p.phi_hat();

  T ct2 =  t2*v_theta_2;
  T st2 = (t2^v2)*v_theta_2;

  double c2theta = ct2*ct2-st2*st2;
  double s2theta = 2*st2*ct2;

  Point_2d g2;
 
  g2[0] = c2theta * p.shear[0] + s2theta * p.shear[1];
  g2[1] = c2theta * p.shear[1] - s2theta * p.shear[0];
  
  return g2;
}

/// flat sky approximation for fthe tangential and x components of polarization realtive to the vector seporating the pointst
Point_2d tangent_cross_shear(
                             const Point_2d &xo
                             ,const Point_2d &x
                             ,const Point_2d &gamma
                             ){
  
  Point_2d t = (x-xo);
  t.unitize();
  
  double c2theta = t[0]*t[0]-t[1]*t[1];
  double s2theta = 2*t[0]*t[1];
  
  Point_2d g2;

  g2[0] = c2theta * gamma[0] - s2theta * gamma[1];
  g2[1] = c2theta * gamma[1] + s2theta * gamma[0];
  
  return g2;
}

/// finds the ++ , xx and x+ products of shear in spherical coordinates
template <typename T>
Point_3d<double> correlate(Polar<T> &p1,Polar<T> &p2){

  Point_2d g1,g2;
  
  g1 = tangent_cross_shear(p2, p1);
  g2 = tangent_cross_shear(p1, p2);

  return Point_3d<double>(g1[0]*g2[0],g1[1]*g2[1],g1[0]*g2[1]);
}

/// finds the ++ , xx and x+ products of shear with a plat sky approximation
Point_3d<double> correlate(Point_2d &x1,Point_2d &g1,Point_2d &x2,Point_2d &g2){

  Point_2d gt1,gt2;
  
  gt1 = tangent_cross_shear(x1,x2,g2);
  gt2 = tangent_cross_shear(x2,x1,g1);

  return Point_3d<double>(gt1[0]*gt2[0],gt1[1]*gt2[1],gt1[0]*gt2[1]);
}


#endif /* shear_h */
