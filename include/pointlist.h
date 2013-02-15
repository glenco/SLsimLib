/*
 * point.h
 *
 *  Created on: Nov 15, 2010
 *      Author: bmetcalf
 *
 *      Defines Point and Branch type.
 *
 *      The Branch type needs to be defined here so that point.leaf can be defined.
 */


#ifndef pointtypes_declare
#define pointtypes_declare

#include "standard.h"
#include <point.h>

//struct Point;
//struct Branch;

/** \brief A point on the source or image plane that contains a position and the lensing quantities */
/*typedef struct Point{
  struct Point *next;    // pointer to next point in linked list
  struct Point *prev;
  struct Point *image;  // pointer to point on image or source plane
  unsigned long id;
  double *x;         // the position of the point
  unsigned long head;         // marks beginning of allocated array of points for easy deallocation
  Boo in_image; // marks if point is in image

  // redundant information in image and source points
  KappaType kappa;        // surface density
  KappaType gamma[3];    // shear, third component is the rotation quantity that is only non-zero for multi-plane lensing
  KappaType dt;          // time delay
  float invmag;     // inverse of magnification
  double gridsize;   // the size of the most refined grid the point is in
  float surface_brightness;  // the surface brightness at this points

  struct Branch;
  Branch *leaf;

} Point;
/*
struct branchstruct{
  struct Point *points;        // pointer to first points in Branch
  unsigned long npoints;
  double center[2];
  int level;
  unsigned long number;
  double boundary_p1[2];
  double boundary_p2[2];
  struct branchstruct *child1;
  struct branchstruct *child2;
  struct branchstruct *brother;
  struct branchstruct *prev;
  /// Marks point as start of a level of refinement
  bool refined;
} ;

typedef struct branchstruct Branch;
*/



#endif
