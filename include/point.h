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

#ifndef pi
#define pi  3.141593
#endif

#ifndef error_message
#define error_message
#define ERROR_MESSAGE() std::printf("ERROR: file: %s line: %i\n",__FILE__,__LINE__)
#endif

#ifndef bool_declare
#define bool_declare
typedef enum {FALSE, TRUE, MAYBE} Boo;
#endif

struct branchstruct;

/** \brief A point on the source or image plane that contains a position and the lensing quantities */
typedef struct Point{
  struct Point *next;    // pointer to next point in linked list
  struct Point *prev;
  struct Point *image;  // pointer to point on image or source plane
  unsigned long id;
  double *x;         // the position of the point
  unsigned long head;         // marks beginning of allocated array of points for easy deallocation
  Boo in_image; // marks if point is in image

  // redundant information in image and source points
  double kappa;        // surface density
  double gamma[3];    // shear, third component is the rotation quantity that is only non-zero for multi-plane lensing
  double dt;          // time delay
  double invmag;     // inverse of magnification
  double gridsize;   // the size of the most refined grid the point is in
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
  /// Marks point as start of a level of refinement
  bool refined;
} ;

typedef struct branchstruct Branch;

#endif
