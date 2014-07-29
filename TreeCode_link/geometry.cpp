//
//  geometry.cpp
//  GLAMER
//
//  Created by bmetcalf on 7/24/14.
//
//

#include "geometry.h"

/// output cartisian coordinates of the point
void Utilities::Geometry::SphericalPoint::sphericalTOcartisian(PosType x[]) const{
  x[0] = r*cos(theta)*cos(phi);
  x[1] = r*cos(theta)*sin(phi);
  x[2] = r*sin(theta);
}

/// set the spherical coordinates of the point from the cartisian coordinates
void Utilities::Geometry::SphericalPoint::cartisianTOspherical(PosType const x[]){
  r = sqrt( x[0]*x[0] + x[1]*x[1] +x[2]*x[2]);
  theta = asin(x[2]/r);
  phi = atan2(x[1],x[0]);
}
      
      /** \brief Calculates the stereographic projection of the point onto a plane.
       *
       * The result is in radian units.  Near the central point this is a rectolinear projection
       * onto a tangent plane.
       */
void Utilities::Geometry::SphericalPoint::StereographicProjection(
      const SphericalPoint &central   /// point on the sphere where the tangent plane touches
      ,PosType x[]             /// 2D output coordinate on projection
  ) const{

  PosType k = 2/( 1 + sin(central.theta)*sin(theta) + cos(central.theta)*cos(theta)*cos(phi - central.phi) );
  
  x[0] = k*(cos(theta)*sin(phi - central.phi));
  x[1] = k*(cos(central.theta)*sin(theta) - sin(central.theta)*cos(theta)*cos(phi - central.phi));
}
      
/** \brief Calculates the orthographic projection of the point onto a plane.
*
* The result is in radian units.  Near the central point this is a rectolinear projection
* onto a tangent plane.
*/
void Utilities::Geometry::SphericalPoint::OrthographicProjection(
                            const SphericalPoint &central   /// point on the sphere where the tangent plane touches
                            ,PosType x[]             /// 2D output coordinate on projection
                            ) const{
  x[0] = cos(theta)*sin(phi - central.phi);
  x[1] = cos(central.theta)*sin(theta) - sin(central.theta)*cos(theta)*cos(phi - central.phi);
}
      
/** \brief Convert from an orthographic projection of the plane onto the unit sphere
*/
void Utilities::Geometry::SphericalPoint::InverseOrthographicProjection(
                                   const SphericalPoint &central   /// point on the sphere where the tangent plane touches
                                    ,PosType const x[]             /// 2D output coordinate on projection
                                    ){
  PosType rho = sqrt(x[0]*x[0] + x[1]*x[1]);
  PosType c = asin(rho);
  r=1.0;
  theta = asin( cos(c)*sin(central.theta) + x[1]*sin(c)*cos(central.theta)/rho );
  phi = central.phi + atan2(x[0]*sin(c),rho*cos(central.theta)*cos(c)
                                  - x[1]*sin(central.theta)*sin(c) );
}

///  3 dimensional distance between points
PosType Utilities::Geometry::Seporation(const SphericalPoint &p1,const SphericalPoint &p2){
  PosType x1[3],x2[3];

  p1.sphericalTOcartisian(x1);
  p2.sphericalTOcartisian(x2);
  return sqrt( (x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1])
                  + (x1[2]-x2[2])*(x1[2]-x2[2]) );
}
    
///  Angular seporation between points
PosType Utilities::Geometry::AngleSeporation(const SphericalPoint &p1,const SphericalPoint &p2){
  return acos(sin(p1.theta)*sin(p2.theta) + cos(p1.theta)*cos(p2.theta)*cos(p1.phi-p2.phi));
}
