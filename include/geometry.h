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
    
    
    class Quaternian{
    public:
      Quaternian(double s,double x,double y,double z){
        v[0] = s;
        v[1] = x;
        v[2] = y;
        v[3] = z;
      }
      Quaternian(const Quaternian &q){
        v[0] = q.v[0];
        v[1] = q.v[1];
        v[2] = q.v[2];
        v[3] = q.v[3];
      }
      Quaternian(Point_3d p){
        v[0] = 0;
        v[1] = p[0];
        v[2] = p[1];
        v[3] = p[2];
      }
      Quaternian(SphericalPoint sp){
        sp.sphericalTOcartisian(v+1);
        v[0] = 0;
      }
      
      Quaternian &operator=(const Quaternian &q){
        v[0] = q.v[0];
        v[1] = q.v[1];
        v[2] = q.v[2];
        v[3] = q.v[3];
        
        return *this;
      }
      
      Quaternian operator*(const Quaternian &q) const{
        return Quaternian(
                          v[0]*q.v[0] - v[1]*q.v[1] - v[2]*q.v[2] - v[3]*q.v[3],
                          v[0]*q.v[1] + v[1]*q.v[0] + v[2]*q.v[3] - v[3]*q.v[2],
                          v[0]*q.v[2] - v[1]*q.v[3] + v[2]*q.v[0] + v[3]*q.v[1],
                          v[0]*q.v[3] + v[1]*q.v[2] - v[2]*q.v[1] + v[3]*q.v[0]);
      }
      
      /// returns the conjugate of the quaternian
      Quaternian conj() const{
        return Quaternian(v[0],-v[1],-v[2],-v[3]);
      }
      /** \brief returns the reciprocal of the quaternian
      
      This is the quaternian whose product with the original quaternian is 
       the real number 1.
      */
      Quaternian inverse() const{
        return this->conj()/(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
      }
      
      /// the norm of the quaternian
      double norm() const{
        return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
      }
      
      /// returns the quaternian rotated around the x-axis by theta in radians
      Quaternian RotateX(double theta){
        Quaternian R = q_x_rotation(theta);
        return R*(*this)*R.conj();
      }
      /// returns the quaternian rotated around the y-axis by theta in radians
      Quaternian RotateY(double theta){
        Quaternian R = q_y_rotation(theta);
        return R*(*this)*R.conj();
      }
      /// returns the quaternian rotated around the z-axis by theta in radians
      Quaternian RotateZ(double theta){
        Quaternian R = q_z_rotation(theta);
        return R*(*this)*R.conj();
      }
      /** returns the quaternian rotated first around the z-axis by phi and then
      around the z-axis by theta, this is more efficient than doing two rotation 
       sequentially.  Usefull for recentering a field of points.
       */
      Quaternian RotateYZ(double theta,double phi){
        Quaternian R = q_z_rotation(theta)*q_z_rotation(phi);
        return R*(*this)*R.conj();
      }

      /** returns the quaternian rotated with the rotation quaternian R.
       It is assumed that R has norm = 1, ie ||R|| = 1 */
      Quaternian Rotate(const Quaternian &R) const{
        return R*(*this)*R.conj();
      }
      
      Quaternian operator+(const Quaternian &q) const{
        Quaternian p = q;
        
        p.v[0] += v[0];
        p.v[1] += v[1];
        p.v[2] += v[2];
        p.v[3] += v[3];
        
        return p;
      }
      Quaternian operator-(const Quaternian &q) const{
        Quaternian p = *this;
        
        p.v[0] -= q.v[0];
        p.v[1] -= q.v[1];
        p.v[2] -= q.v[2];
        p.v[3] -= q.v[3];
        
        return p;
      }
      
      /// division by a scaler
      Quaternian operator/(double s) const{
       return Quaternian(v[0]/s,v[1]/s,v[2]/s,v[3]/s);
      }
      /// division by a scaler

      Quaternian operator*(double s) const{
        return Quaternian(v[0]*s,v[1]*s,v[2]*s,v[3]*s);
      }
      
      /// cross (outer) product
      Quaternian cross(const Quaternian &q) const{
        return (*this*q - this->conj()*q.conj());
      }

      /// scalar product of vector part of the quaternians
      double scalor(const Quaternian &q) const{
        return v[1]*q.v[1] + v[2]*q.v[2] + v[3]*q.v[3];
      }
      
      /// returns the vector part of the quaternian as a Point_3d
      Point_3d to_point_3d() const{
        return Point_3d(v[1],v[2],v[3]);
      }

      /// returns the vector part of the quaternian as a SpericalPoint
      SphericalPoint to_SpericalPoint() const{
        SphericalPoint sp;
        sp.cartisianTOspherical(v+1);
        return sp;
      }

      
      /// the components of the quaternian
      double v[4];

      /// returns the rotation quaternian for a ratation around the x-axis
      static Quaternian q_x_rotation(double theta){
        return Quaternian(cos(theta/2),sin(theta/2),0,0);
      }
      /// returns the rotation quaternian for a ratation around the y-axis
      static Quaternian q_y_rotation(double theta){
        return Quaternian(cos(theta/2),0,sin(theta/2),0);
      }
      /// returns the rotation quaternian for a ratation around the z-axis
      static Quaternian q_z_rotation(double theta){
        return Quaternian(cos(theta/2),0,0,sin(theta/2));
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

#endif /* defined(__GLAMER__geometry__) */
