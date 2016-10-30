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
    
    /// represents a point in spherical coordinates theta = 0 is the equator
    class SphericalPoint{
      
    public:
      SphericalPoint(PosType r,PosType theta,PosType phi):r(r),theta(theta),phi(phi){};
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

  }
}

std::ostream &operator<<(std::ostream &os,const Utilities::Geometry::SphericalPoint &p);
#endif /* defined(__GLAMER__geometry__) */
