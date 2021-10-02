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
#include <array>
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

/// represents a point in spherical coordinates, theta = 0 is equator
template <typename T = double>
class SphericalPoint{
  
public:
  SphericalPoint(T r,T theta,T phi):
  r(r),theta(theta),phi(phi){};
  SphericalPoint():r(0),theta(0),phi(0){};
  SphericalPoint(Point_3d<T> &x){
    TOspherical(x);
  }
  
  SphericalPoint & operator=(const SphericalPoint &p){
    if(&p != this){
      r = p.r;
      theta = p.theta;
      phi = p.phi;
    }
    return *this;
  }
  
  SphericalPoint & operator=(Point_3d<T> &x){
    TOspherical(x);
    return *this;
  }
  
  bool operator==(const SphericalPoint &p) const{
    if(&p == this) return true;
    return (r==p.r)*(theta==p.theta)*(phi==p.phi);
  }
  
  
  T r;
  T theta;
  T phi;
  
  /// output Cartesian coordinates of the point
  void TOcartisian(T x[]) const;
  Point_3d<T> TOcartisian() const;
  
  void cartisianTOspherical(T const x[]);
  void TOspherical(Point_3d<T> &x);
  
  void StereographicProjection(const SphericalPoint &central
                               ,T x[]) const;
  void StereographicProjection(const SphericalPoint &central,Point_2d &x) const;
  Point_2d StereographicProjection(const SphericalPoint &central) const;
  
  void OrthographicProjection(const SphericalPoint &central
                              ,T x[]) const;
  Point_2d OrthographicProjection(const SphericalPoint &central) const;
  void InverseOrthographicProjection(const SphericalPoint &central,T const x[]);
  void InverseOrthographicProjection(const SphericalPoint &central,const Point_2d &x);
  SphericalPoint<T> InverseOrthographicProjection(const Point_2d &x);

  
  /// angle between points.  This uses the haversine formula that is more stable for small angles than the more commin formula
  T angular_separation(SphericalPoint &p){
    double s1 = sin((theta - p.theta)/2);
    double s2 = sin((phi - p.phi)/2);
    
    return 2*asin(sqrt( s1*s1 + cos(theta)*cos(p.theta)*s2*s2 ) );
  }
  
  /// unit vector in phi direction
  Point_3d<T> unitPhi();
  /// unit vector in theta direction
  Point_3d<T> unitTheta();

  /// the angle between the orthographic x-axis  and the constant theta curve
   T OrthographicAngleTheta(const SphericalPoint &central);
  /// the angle between the orthographic x-axis  and the constant Phi curve
   T OrthographicAnglePhi(const SphericalPoint &central);
  
  // returns the unit theta vector
  Point_3d<T> theta_hat() const{
    Point_3d<T> p;
    p[0] = -sin(theta)*cos(phi);
    p[1] = -sin(theta)*sin(phi);
    p[2] = cos(theta);
    
    return p;
  }

  // returns the unit phi vector
  Point_3d<T> phi_hat() const{
     Point_3d<T> p;
     p[0] = -sin(phi);
     p[1] =  cos(phi);
     p[2] = 0;
     
     return p;
   }

};
 

/** \brief Quaternion class that is especially useful for rotations.
 
 
 <pre>
 The Quaternion can be easily transformed between classes Point_3d and Utilities::Geomotry::SphericalPoint.
 
 Below is an example of a rotation of a vector initially defined as a Utilities::Geometry::SphericalPoint.
 
 
 #include "geometry.h"
 {
 using Utilities::Geometry::SphericalPoint;
 using Utilities::Geometry::Quaternion;
 
 SphericalPoint sp;
 
 sp.phi = -pi/2;
 sp.theta = pi/3;
 sp.r = 1.0;
 
 cout << sp.r << " " << sp.theta << " " << sp.phi << endl;
 
 // make a rotation Quaternion that will rotate the vector to the x-axis.
 Quaternion R = Quaternion::q_y_rotation(-sp.theta)*Quaternion::q_z_rotation(-sp.phi);
 
 
 // construct a Quaternion from a SphericalPoint
 Quaternion q(sp);
 
 // apply the rotation
 p = q.Rotate(R);
 // this could also be written as
 p = R*q*R.conj();
 
 // rotate in place, same as above but with one less copy
 q.RotInplace(R);
 
 // convert to a Point_3d
 Point_3d p = q.to_point_3d();
 
 cout << q.to_point_3d() << endl;
 
 // convert back to a SphericalPoint
 cout << q.to_SpericalPoint() << endl;
 }
 
 <\pre>
 */
template <typename T = double>
class Quaternion{
public:
  Quaternion(){
    v[0] = v[1] = v[2] = v[3] = 0;
  }
  Quaternion(T s,T x,T y,T z){
    v[0] = s;
    v[1] = x;
    v[2] = y;
    v[3] = z;
  }
  Quaternion(const Quaternion &q){
    v[0] = q.v[0];
    v[1] = q.v[1];
    v[2] = q.v[2];
    v[3] = q.v[3];
  }
  Quaternion(Point_3d<T> p){
    v[0] = 0;
    v[1] = p[0];
    v[2] = p[1];
    v[3] = p[2];
  }
  Quaternion(SphericalPoint<T> sp){
    sp.TOcartisian(v+1);
    v[0] = 0;
  }
  
  Quaternion &operator=(const Quaternion &q){
    v[0] = q.v[0];
    v[1] = q.v[1];
    v[2] = q.v[2];
    v[3] = q.v[3];
    
    return *this;
  }
  
  Quaternion operator*(const Quaternion &q) const{
    return Quaternion(
                      v[0]*q.v[0] - v[1]*q.v[1] - v[2]*q.v[2] - v[3]*q.v[3],
                      v[0]*q.v[1] + v[1]*q.v[0] + v[2]*q.v[3] - v[3]*q.v[2],
                      v[0]*q.v[2] - v[1]*q.v[3] + v[2]*q.v[0] + v[3]*q.v[1],
                      v[0]*q.v[3] + v[1]*q.v[2] - v[2]*q.v[1] + v[3]*q.v[0]);
  }
  
  /// returns the conjugate of the Quaternion
  Quaternion conj() const{
    return Quaternion(v[0],-v[1],-v[2],-v[3]);
  }
  /** \brief returns the reciprocal of the Quaternion
   
   This is the Quaternion whose product with the original Quaternion is
   the real number 1.
   */
  Quaternion inverse() const{
    return this->conj()/(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
  }
  
  /// the norm of the Quaternion
  double norm() const{
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
  }
  
  /// returns the Quaternion rotated around the x-axis by theta in radians
  Quaternion RotateX(T theta){
    Quaternion R = q_x_rotation(theta);
    return R*(*this)*R.conj();
  }
  /// returns the Quaternion rotated around the y-axis by theta in radians
  Quaternion RotateY(T theta){
    Quaternion R = q_y_rotation(theta);
    return R*(*this)*R.conj();
  }
  /// returns the Quaternion rotated around the z-axis by theta in radians
  Quaternion RotateZ(T theta){
    Quaternion R = q_z_rotation(theta);
    return R*(*this)*R.conj();
  }
  /** returns the Quaternion rotated first around the z-axis by phi and then
   around the z-axis by theta, this is more efficient than doing two rotation
   sequentially.  Usefull for recentering a field of points.
   */
  Quaternion RotateYZ(T theta,T phi){
    Quaternion R = q_z_rotation(theta)*q_z_rotation(phi);
    return R*(*this)*R.conj();
  }
  
  /** returns the Quaternion rotated with the rotation Quaternion R.
   It is assumed that R has norm = 1, ie ||R|| = 1 */
  Quaternion Rotate(const Quaternion &R) const{
    return R*(*this)*R.conj();
  }
  
  /** returns the Quaternion rotated with the rotation Quaternion R.
   It is assumed that R has norm = 1, ie ||R|| = 1 */
  void RotInplace(const Quaternion &R){
    *this = R*(*this)*R.conj();
  }
  
  
  /// rotate a Point_3d using a rotation Quaternion
  static Point_3d<T> Rotate(const Point_3d<T> &p
                            ,const Quaternion &R){
    Quaternion q(p);
    q.RotInplace(R);
    return q.to_point_3d();
  }
  
  /// rotate a SpericalPoint using a rotation Quaternion
  
  static SphericalPoint<T> Rotate(const SphericalPoint<T> &p
                                  ,const Quaternion &R){
    Quaternion q(p);
    q.RotInplace(R);
    return q.to_SpericalPoint();
  }
  
  Quaternion operator+(const Quaternion &q) const{
    return Quaternion(v[0] + q.v[0] , v[1] + q.v[1] , v[2] + q.v[2] , v[3] + q.v[3] );
    
  }
  Quaternion operator-(const Quaternion &q) const{
    return Quaternion(v[0] - q.v[0] , v[1] - q.v[1] , v[2] - q.v[2] , v[3] - q.v[3] );
  }
  
  /// division by a scaler
  Quaternion operator/(double s) const{
    return Quaternion(v[0]/s,v[1]/s,v[2]/s,v[3]/s);
  }
  /// multiply by a scaler
  Quaternion operator*(double s) const{
    return Quaternion(v[0]*s,v[1]*s,v[2]*s,v[3]*s);
  }
  
  /// cross (outer) product
  Point_3d<T> cross(const Quaternion &q) const{
    return (*this*q).to_point_3d();
  }
  
  /// scalar product of vector part of the Quaternions
  T scalor(const Quaternion &q) const{
    return v[1]*q.v[1] + v[2]*q.v[2] + v[3]*q.v[3];
  }
  
  /// returns the vector part of the Quaternion as a Point_3d
  Point_3d<T> to_point_3d() const{
    return Point_3d<T>(v[1],v[2],v[3]);
  }
  
  /// returns the vector part of the Quaternion as a SpericalPoint
  SphericalPoint<T> to_SpericalPoint() const{
    SphericalPoint<T> sp;
    sp.cartisianTOspherical(v+1);
    return sp;
  }
  
  
  /// the components of the Quaternion
  T v[4];
  /// components, 0,1,2,3
  T & operator[](int i){return v[i];}
  
  /// returns the rotation Quaternion for a rotation around the x-axis
  static Quaternion q_x_rotation(T theta){
    return Quaternion(cos(theta/2),sin(theta/2),0,0);
  }
  /// returns the rotation Quaternion for a rotation around the y-axis
  static Quaternion q_y_rotation(T theta){
    return Quaternion(cos(theta/2),0,-sin(theta/2),0);
  }
  /// returns the rotation Quaternion for a rotation around the z-axis
  static Quaternion q_z_rotation(T theta){
    return Quaternion(cos(theta/2),0,0,sin(theta/2));
  }
  
};

/// Determine if line segments a1a2 and b1b2 intersect.  Sharing an endpoint does not count as intersecting
bool intersect(const PosType a1[],const PosType a2[]
               ,const PosType b1[],const PosType b2[]);
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

}
}

namespace Utilities {
namespace Geometry{

/// output cartisian coordinates
template <typename T>
void SphericalPoint<T>::TOcartisian(T x[]) const{
  x[0] = r*cos(theta)*cos(phi);
  x[1] = r*cos(theta)*sin(phi);
  x[2] = r*sin(theta);
}

/// output cartisian coordinates of the point
template <typename T>
Point_3d<T> SphericalPoint<T>::TOcartisian() const{
  return Point_3d<T>(r*cos(theta)*cos(phi)
                     ,r*cos(theta)*sin(phi),r*sin(theta) );
}

/// set the spherical coordinates of the point from the cartisian coordinates
template <typename T>
void SphericalPoint<T>::cartisianTOspherical(T const x[]){
  r = sqrt( x[0]*x[0] + x[1]*x[1] +x[2]*x[2]);
  theta = asin(x[2]/r);
  phi = atan2(x[1],x[0]);
}

template <typename T>
void SphericalPoint<T>::TOspherical(Point_3d<T> &x){
  r = sqrt( x[0]*x[0] + x[1]*x[1] +x[2]*x[2]);
  theta = asin(x[2]/r);
  phi = atan2(x[1],x[0]);
}


/** \brief Calculates the stereographic projection of the point onto a plane.
 *
 * The result is in radian units.  Near the central point this is a rectolinear projection
 * onto a tangent plane.
 */
template <typename T>
void SphericalPoint<T>::StereographicProjection(
                                                const SphericalPoint<T> &central   /// point on the sphere where the tangent plane touches
                                                ,T x[]             /// 2D output coordinate on projection
) const{
  double ct = cos(theta),st = sin(theta);
  double cosphi = cos(phi - central.phi);
  double so = sin(central.theta);
  double co = cos(central.theta);

  PosType k = 2/( 1 + so*st + co*ct*cosphi );
  
  x[0] = k*(ct*sin(phi - central.phi));
  x[1] = k*(co*st - so*ct*cosphi);
}
template <typename T>
void SphericalPoint<T>::StereographicProjection(
                                                const SphericalPoint<T> &central   /// point on the sphere where the tangent plane touches
                                                ,Point_2d &x  /// 2D output coordinate on projection
) const{
  double ct = cos(theta),st = sin(theta);
  double cosphi = cos(phi - central.phi);
  double so = sin(central.theta);
  double co = cos(central.theta);

  PosType k = 2/( 1 + so*st + co*ct*cosphi );
  
  x[0] = k*(ct*sin(phi - central.phi));
  x[1] = k*(co*st - so*ct*cosphi);
}

template <typename T>
Point_2d SphericalPoint<T>::StereographicProjection(
                                                    const SphericalPoint<T> &central   /// point on the sphere where the tangent plane touches
) const{
  double ct = cos(theta),st = sin(theta);
  double cosphi = cos(phi - central.phi);
  double so = sin(central.theta);
  double co = cos(central.theta);

  PosType k = 2/( 1 + so*st + co*ct*cosphi );
  
  return Point_2d(k*(ct*sin(phi - central.phi)),k*(co*st - so*ct*cosphi));
}


/** \brief Calculates the orthographic projection of the point onto a plane.
 *
 * The result is in radian units.  Near the central point this is a rectolinear projection
 * onto a tangent plane.
 * Points in the oposite hemosphere are maped to points on the border at |x| = 1
 */
template <typename T>
void SphericalPoint<T>::OrthographicProjection(
                                               const SphericalPoint<T> &central   /// point on the sphere where the tangent plane touches
                                               ,T x[]             /// 2D output coordinate on projection
) const{
  double ct = cos(theta),st = sin(theta);
  double cosphi = cos(phi - central.phi);
  double so = sin(central.theta);
  double co = cos(central.theta);
  
  x[0] = ct*sin(phi - central.phi);
  x[1] = co*st - so*ct*cosphi;
  
  if(st*so + ct*co*cosphi < 0 ){  // point is on the other hemosphere
    double s = sqrt( x[0]*x[0] + x[1]*x[0] );
    x[0] /=s;
    x[1] /=s;
  }
}

template <typename T>
Point_2d SphericalPoint<T>::OrthographicProjection(
                                                   const SphericalPoint<T> &central   /// point on the sphere where the tangent plane touches
) const{
  double c = cos(theta);
  return Point_2d( c*sin(phi - central.phi),
                  cos(central.theta)*sin(theta) - sin(central.theta)*c*cos(phi - central.phi) );
}

/** \brief Convert from an orthographic projection of the plane onto the unit sphere
 */
template <typename T>
void SphericalPoint<T>::InverseOrthographicProjection(
                                                      const SphericalPoint<T> &central   /// point on the sphere where the tangent plane touches
                                                      ,T const x[]             /// 2D output coordinate on projection
){
  PosType rho = sqrt(x[0]*x[0] + x[1]*x[1]);
  PosType c = asin(rho);
  r=1.0;
  theta = asin( cos(c)*sin(central.theta) + x[1]*cos(central.theta) );
  phi = central.phi + atan2(x[0] , cos(central.theta)*cos(c)
                            - x[1]*sin(central.theta) );
}

template <typename T>
void SphericalPoint<T>::InverseOrthographicProjection(
                                                      const SphericalPoint<T> &central   /// point on the sphere where the tangent plane touches
                                                      ,const Point_2d &x             /// 2D output coordinate on projection
){
  PosType rho = sqrt(x[0]*x[0] + x[1]*x[1]);
  PosType c = asin(rho);
  r=1.0;
  double at=cos(c)*sin(central.theta) + x[1]*cos(central.theta);
  at = (at > 1) ? 1 : at;
  at = (at < -1) ? -1 : at;
  theta = asin( at );
  phi = central.phi + atan2(x[0] , cos(central.theta)*cos(c)
                             - x[1]*sin(central.theta) );
}

template <typename T>
T SphericalPoint<T>::OrthographicAngleTheta(
          const SphericalPoint<T> &central   /// point on the sphere where the tangent plane touches
  ){
    return atan2( sin(central.theta) * sin(phi - central.phi)
                 , cos(phi - central.phi) );
  }

  template <typename T>
  T SphericalPoint<T>::OrthographicAnglePhi(
            const SphericalPoint<T> &central   /// point on the sphere where the tangent plane touches
    ){
      return atan2( cos(central.theta) * cos(theta) + sin(central.theta) * sin(central.theta)
                   *cos(phi - central.phi)
          , -sin(phi - central.phi) * sin(phi - central.phi) );
    }

//// deprojection where this is the center of the projection

template <typename T>
SphericalPoint<T> SphericalPoint<T>::InverseOrthographicProjection(
                                                                   const Point_2d &x             /// 2D coordinate on projection
){
  PosType rho = sqrt(x[0]*x[0] + x[1]*x[1]);
  if(rho==0) return *this;
  
  PosType c = asin(rho);
  T new_theta = asin( cos(c)*sin(theta) + x[1]*sin(c)*cos(theta)/rho );
  T new_phi = phi + atan2(x[0]*sin(c),rho*cos(theta)*cos(c)
                          - x[1]*sin(theta)*sin(c) );
  
  return SphericalPoint<T>(1,new_theta,new_phi);
}

template <typename T>
Point_3d<T> SphericalPoint<T>::unitPhi(){
  return Point_3d<T>(-sin(phi),cos(phi),0);
}

template <typename T>
Point_3d<T> SphericalPoint<T>::unitTheta(){
  return Point_3d<T>(-cos(phi)*sin(theta),-sin(theta)*sin(phi),cos(theta));
}

///  3 dimensional distance between points
template <typename T>
T Seporation(const SphericalPoint<T> &p1
             ,const SphericalPoint<T> &p2){
  T x1[3],x2[3];
  
  p1.TOcartisian(x1);
  p2.TOcartisian(x2);
  return sqrt( (x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1])
              + (x1[2]-x2[2])*(x1[2]-x2[2]) );
}

///  Angular seporation between points
template <typename T>
T AngleSeporation(const SphericalPoint<T> &p1
                  ,const SphericalPoint<T> &p2){
  if(p1 == p2) return 0;
  double ans = sin(p1.theta)*sin(p2.theta) + cos(p1.theta)*cos(p2.theta)*cos(p1.phi-p2.phi);
  if(fabs(ans) > 1) return 0;
  return acos(sin(p1.theta)*sin(p2.theta) + cos(p1.theta)*cos(p2.theta)*cos(p1.phi-p2.phi));
}

}
}

template <typename T>
std::ostream &operator<<(std::ostream &os, Utilities::Geometry::SphericalPoint<T> const &p){
  return os << "r: " << p.r << " theta: " << p.theta << " phi: " << p.phi;
}


#endif /* defined(__GLAMER__geometry__) */
