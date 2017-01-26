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
    
    /// represents a point in spherical coordinates theta = 0 is the equator
    class SphericalPoint{
      
    public:
      SphericalPoint(PosType r,PosType theta,PosType phi):r(r),theta(theta),phi(phi){};
      SphericalPoint(const Point_3d &x){
        cartisianTOspherical(x.x);
      }
      SphericalPoint():r(0),theta(0),phi(0){};

      SphericalPoint &operator=(Point_3d &x){
        cartisianTOspherical(x.x);
        return *this;
      }
      
      PosType r;
      PosType theta;
      PosType phi;
      
      /// output Cartesian coordinates of the point
      void sphericalTOcartisian(PosType x[]) const;
      void cartisianTOspherical(PosType const x[]);
      void StereographicProjection(const SphericalPoint &central,PosType x[]) const;
      void OrthographicProjection(const SphericalPoint &central,PosType x[]) const;
      void InverseOrthographicProjection(const SphericalPoint &central,PosType const x[]);
      PosType AngleSeporation(const SphericalPoint &p2);

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
     q = q.Rotate(R);
     
     // convert to a Point_3d
     Point_3d p = q.to_point_3d();
     
     cout << q.to_point_3d() << endl;
     
     // convert back to a SphericalPoint
     cout << q.to_SpericalPoint() << endl;
     }
     
     <\pre>
     */
    
    class Quaternion{
    public:
      Quaternion(){
        v[0] = v[1] = v[2] = v[3] = 0.0;
      }
      Quaternion(double s,double x,double y,double z){
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
      Quaternion(Point_3d p){
        v[0] = 0;
        v[1] = p[0];
        v[2] = p[1];
        v[3] = p[2];
      }
      Quaternion(SphericalPoint sp){
        sp.sphericalTOcartisian(v+1);
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
      Quaternion RotateX(double theta){
        Quaternion R = q_x_rotation(theta);
        return R*(*this)*R.conj();
      }
      /// returns the Quaternion rotated around the y-axis by theta in radians
      Quaternion RotateY(double theta){
        Quaternion R = q_y_rotation(theta);
        return R*(*this)*R.conj();
      }
      /// returns the Quaternion rotated around the z-axis by theta in radians
      Quaternion RotateZ(double theta){
        Quaternion R = q_z_rotation(theta);
        return R*(*this)*R.conj();
      }
      /** returns the Quaternion rotated first around the z-axis by phi and then
      around the z-axis by theta, this is more efficient than doing two rotation 
       sequentially.  Usefull for recentering a field of points.
       */
      Quaternion RotateYZ(double theta,double phi){
        Quaternion R = q_z_rotation(theta)*q_z_rotation(phi);
        return R*(*this)*R.conj();
      }

      /** returns the Quaternion rotated with the rotation Quaternion R.
       It is assumed that R has norm = 1, ie ||R|| = 1 */
      Quaternion Rotate(const Quaternion &R) const{
        return R*(*this)*R.conj();
      }
      
      
      /// rotate a Point_3d using a rotation Quaternion
      Point_3d Rotate(const Point_3d &p){
        Quaternion q(p);
        q.Rotate(*this);
        return q.to_point_3d();
      }

      /// rotate a SpericalPoint using a rotation Quaternion
      SphericalPoint Rotate(const SphericalPoint &p){
        Quaternion q(p);
        q.Rotate(*this);
        return q.to_SpericalPoint();
      }
      
      Quaternion operator+(const Quaternion &q) const{
        Quaternion p = q;
        
        p.v[0] += v[0];
        p.v[1] += v[1];
        p.v[2] += v[2];
        p.v[3] += v[3];
        
        return p;
      }
      Quaternion operator-(const Quaternion &q) const{
        Quaternion p = *this;
        
        p.v[0] -= q.v[0];
        p.v[1] -= q.v[1];
        p.v[2] -= q.v[2];
        p.v[3] -= q.v[3];
        
        return p;
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
      Point_3d cross(const Quaternion &q) const{
        return (*this*q).to_point_3d();
      }

      /// scalar product of vector part of the Quaternions
      double scalor(const Quaternion &q) const{
        return v[1]*q.v[1] + v[2]*q.v[2] + v[3]*q.v[3];
      }
      
      /// returns the vector part of the Quaternion as a Point_3d
      Point_3d to_point_3d() const{
        return Point_3d(v[1],v[2],v[3]);
      }

      /// returns the vector part of the Quaternion as a SpericalPoint
      SphericalPoint to_SpericalPoint() const{
        SphericalPoint sp;
        sp.cartisianTOspherical(v+1);
        return sp;
      }

      
      /// the components of the Quaternion
      double v[4];

      /// returns the rotation Quaternion for a ratation around the x-axis
      static Quaternion q_x_rotation(double theta){
        return Quaternion(cos(theta/2),sin(theta/2),0,0);
      }
      /// returns the rotation Quaternion for a ratation around the y-axis
      static Quaternion q_y_rotation(double theta){
        return Quaternion(cos(theta/2),0,-sin(theta/2),0);
      }
      /// returns the rotation Quaternion for a ratation around the z-axis
      static Quaternion q_z_rotation(double theta){
        return Quaternion(cos(theta/2),0,0,sin(theta/2));
      }
      
    };
    
    
    ///  3 dimensional distance between points
    PosType Seporation(const SphericalPoint &p1,const SphericalPoint &p2);
    ///  Angular seporation between points
    PosType AngleSeporation(const SphericalPoint &p1,const SphericalPoint &p2);

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

std::ostream &operator<<(std::ostream &os,Utilities::Geometry::SphericalPoint &p);

#endif /* defined(__GLAMER__geometry__) */
