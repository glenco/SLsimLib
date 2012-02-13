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

#ifndef pi
#define pi  3.141593
#endif

#ifndef error_message
#define error_message
#define ERROR_MESSAGE() std::printf("ERROR: file: %s line: %i\n",__FILE__,__LINE__)
#endif

#ifndef pointtypes_declare
#define pointtypes_declare

struct branchstruct;

/** \brief A point on the source or image plane that contains a position and the lensing quantities */
typedef struct Point{
  struct Point *next;    // pointer to next point in linked list
  struct Point *prev;
  struct Point *image;  // pointer to point on image or source plane
  unsigned long id;
  double kappa;        // surface density
  double gamma[3];    // shear, third component is the rotation quantity that is only non-zero for multi-plane lensing
  double dt;          // time delay
  double invmag;     // inverse of magnification
  double *x;         // the position of the point
  double gridsize;   // the size of the most refined grid the point is in
  unsigned long head;         // marks beginning of allocated array of points for easy deallocation
  bool in_image; // marks if point is in image
  double surface_brightness;  // the surface brightness at this points

  struct branchstruct *leaf;

} Point;

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
} ;

typedef struct branchstruct Branch;

#endif
