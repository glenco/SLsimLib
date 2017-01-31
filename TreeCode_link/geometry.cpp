//
//  geometry.cpp
//  GLAMER
//
//  Created by bmetcalf on 7/24/14.
//
//

#include "geometry.h"
#include "Tree.h"
#include "point.h"

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
PosType Utilities::Geometry::SphericalPoint::AngleSeporation(const SphericalPoint &p2){
  return acos(sin(theta)*sin(p2.theta) + cos(theta)*cos(p2.theta)*cos(phi-p2.phi));
}


bool Utilities::Geometry::Cone::intersect_line_segment(const Point_3d &p1,const Point_3d &p2){
  
  {  // check if either of the end points are within cone
    Point_3d pc = p1-p;
    if( pc*v > costheta*pc.length() && pc*v > 0) return true;
    pc = p2-p;
    if( pc*v > costheta*pc.length() && pc*v > 0 ) return true;
  }
  
  
  Point_3d t = p2 - p1;  // tangent to curve
  Point_3d w = p1 - p;
  
  double tv = t*v,wv=w*v,tw=t*w;
  
  double a = tv*tv - costheta*costheta*t.length_sqr();
  double b = 2*( wv*tv - costheta*costheta*tw );
  double c = wv*wv - costheta*costheta*w.length_sqr();
  
  double tmp = b*b - 4*a*c;
  if(tmp < 0){  // no intersection
    return false;
  }else{
    
    double s1 = (-b + sqrt(tmp))/2/a;
    double s2 = (-b - sqrt(tmp))/2/a;
    
    // case where segment intersects the wall of the cone
    if(s1 < 1 && s1 > 0){
      if( (p1 + t*s1)*v > 0) return true;  // check that it is the forward half of the cone
    }
    if(s2 < 1 && s2 > 0){
      if( (p1 + t*s2)*v > 0) return true;
    }
  }
  
  return false;  // line intersects the cone but not the segment
};

bool Utilities::Geometry::Cone::intersect_face(const Point_3d &p1,const Point_3d &p2,const Point_3d &p3){
  
  Point_3d p4 = p1 + p3 - p2;
  if(intersect_line_segment(p1,p2)) return true;
  if(intersect_line_segment(p2,p3)) return true;
  if(intersect_line_segment(p3,p4)) return true;
  if(intersect_line_segment(p4,p1)) return true;
  
  // is the whole cone in the face the center of the cone will be
  
  Point_3d p12 = p1 - p2,p14 = p1 - p4;
  Point_3d norm = p12.cross(p14);  // normal vector to plane of the face
  
  double area = norm.length();
  
  // point on plane of face
  double s = (p1*norm)/(v*norm);
  if(s < 0) return false;
  Point_3d center = p1 - v*s;
  p12 = p12.cross(center);
  p14 = p14.cross(center);
  
  if(p12*p14 < 0 && p12.length() < area) return true;
  
  return false;
};
bool Utilities::Geometry::Cone::intersect_box(const Point_3d &p1,const Point_3d &p2){
  
  Point_3d d11 = p1 + Point_3d(p2(0)-p1(0),0,0);
  Point_3d d12 = p1 + Point_3d(0,p2(1)-p1(1),0);
  Point_3d d13 = p1 + Point_3d(0,0,p2(2)-p1(2));
  
  if(intersect_face(d11,p1,d12)) return true;
  if(intersect_face(d11,p1,d13)) return true;
  if(intersect_face(d12,p1,d13)) return true;
  
  Point_3d d21 = p2 + Point_3d(p1(0)-p2(0),0,0);
  Point_3d d22 = p2 + Point_3d(0,p1(1)-p2(1),0);
  Point_3d d23 = p2 + Point_3d(0,0,p1(2)-p2(2));
  
  if(intersect_face(d21,p2,d22)) return true;
  if(intersect_face(d21,p2,d23)) return true;
  if(intersect_face(d22,p2,d23)) return true;
  
  return false;
}

bool Utilities::Geometry::intersect(const PosType a1[],const PosType a2[],const PosType b1[],const PosType b2[]){
  
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

int Utilities::Geometry::intersect(const std::vector<Point_2d> &curve){
  if(curve.size() < 4) return 0;
  
  int intersections = 0;
  for(size_t i=0;i<curve.size()-2;++i){
    for(size_t j=i+2;j<curve.size()-1;++j){
      intersections += Utilities::Geometry::intersect(curve[i].x,curve[i+1].x,curve[j].x,curve[j+1].x);
    }
    intersections += Utilities::Geometry::intersect(curve[i].x,curve[i+1].x
                                                    ,curve.back().x,curve[0].x);
  }
  
  return intersections;
}

int Utilities::Geometry::orientation(const PosType p[],const PosType q[],const PosType r[])
{
  // See 10th slides from following link for derivation of the formula
  double val = (q[1] - p[1]) * (r[0] - q[0]) -
  (q[0] - p[0]) * (r[1] - q[1]);
  
  if (val == 0.0) return 0;  // colinear
  
  return (val > 0)? 1: 2; // clock or counterclock wise
}
/** \brief Given three colinear points p, q, r, the function checks if
 point q lies on line segment 'pr', but not at p or r
 */
bool Utilities::Geometry::onSegment(const PosType p[], const PosType q[], const PosType r[])
{
  if (q[0] < MAX(p[0], r[0]) && q[0] > MIN(p[0], r[0]) &&
      q[1] < MAX(p[1], r[1]) && q[1] > MIN(p[1], r[1]))
    return true;
  
  return false;
}

double Utilities::Geometry::AngleBetween2d(double v1[],double v2[]){
  double y = (v1[0] * v2[1]) - (v2[0] * v1[1]);
  double x = (v1[0] * v2[0]) + (v1[1] * v2[1]);

  if(y == 0 && x < 0 ) return pi;
  return atan2(y, x);
}

int Utilities::Geometry::incurve(PosType x[],std::vector<double *> curve){
  int number = 0;
  size_t i;
  
  Point point;
  for(i=0;i<curve.size()-1;++i){
    
    if( (x[1] >= curve[i][1])*(x[1] <= curve[i+1][1]) ){
      if(Utilities::Geometry::orientation(curve[i], x, curve[i+1]) <= 1) ++number;
    }else if( (x[1] <= curve[i][1])*(x[1] > curve[i+1][1]) ){
      if(Utilities::Geometry::orientation(curve[i], x, curve[i+1]) == 2) --number;
    }
    
  }
  
  if( (x[1] >= curve[i][1])*(x[1] <= curve[0][1]) ){
    if(Utilities::Geometry::orientation(curve[i], x, curve[0]) <= 1) ++number;
  }else if( (x[1] <= curve[i][1])*(x[1] > curve[0][1]) ){
    if(Utilities::Geometry::orientation(curve[i], x, curve[0]) == 2) --number;
  }
  
  return number == 0 ? 0 : 1;
}

std::ostream &operator<<(std::ostream &os, Utilities::Geometry::SphericalPoint &p){
  return os << "r: " << p.r << " theta: " << p.theta << " phi: " << p.phi;

}

