//
//  geometry.cpp
//  GLAMER
//
//  Created by bmetcalf on 7/24/14.
//
//

#include "geometry.h"
#include "Tree.h"

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

bool Utilities::Geometry::intersect(PosType a1[],PosType a2[],PosType b1[],PosType b2[]){
  
  if(a1 == b1 || a1 == b2) return false;
  if(a2 == b1 || a2 == b2) return false;
  
  int o1 = Utilities::Geometry::orientation(a1, a2, b1);
  int o2 = Utilities::Geometry::orientation(a1, a2, b2);
  int o3 = Utilities::Geometry::orientation(b1, b2, a1);
  int o4 = Utilities::Geometry::orientation(b1, b2, a2);
  
  // General case
  if (o1 != o2 && o3 != o4)
    return true;
  
  // Special Cases
  // p1, q1 and p2 are colinear and p2 lies on segment p1q1
  if (o1 == 0 && Utilities::Geometry::onSegment(a1, b1, a2)) return true;
  
  // p1, q1 and p2 are colinear and q2 lies on segment p1q1
  if (o2 == 0 && Utilities::Geometry::onSegment(a1, b2, a2)) return true;
  
  // p2, q2 and p1 are colinear and p1 lies on segment p2q2
  if (o3 == 0 && Utilities::Geometry::onSegment(b1, a1, b2)) return true;
  
  // p2, q2 and q1 are colinear and q1 lies on segment p2q2
  if (o4 == 0 && Utilities::Geometry::onSegment(b1, a2, b2)) return true;
  
  return false; // Doesn't fall in any of the above cases
}

int Utilities::Geometry::orientation(PosType p[],PosType q[],PosType r[])
{
  // See 10th slides from following link for derivation of the formula
  int val = (q[1] - p[1]) * (r[0] - q[0]) -
  (q[0] - p[0]) * (r[1] - q[1]);
  
  if (val == 0) return 0;  // colinear
  
  return (val > 0)? 1: 2; // clock or counterclock wise
}
/** \brief Given three colinear points p, q, r, the function checks if
 point q lies on line segment 'pr', but not at p or r
 */
bool Utilities::Geometry::onSegment(PosType p[], PosType q[], PosType r[])
{
  if (q[0] < MAX(p[0], r[0]) && q[0] > MIN(p[0], r[0]) &&
      q[1] < MAX(p[1], r[1]) && q[1] > MIN(p[1], r[1]))
    return true;
  
  return false;
}
