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
      
      /// output cartisian coordinates of the point
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
  }
}

#endif /* defined(__GLAMER__geometry__) */
