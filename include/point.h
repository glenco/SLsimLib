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
#define ERROR_MESSAGE() printf("ERROR: file: %s line: %i\n",__FILE__,__LINE__)
#endif

#ifndef Boolean_declare
#define Boolean_declare
typedef enum {False, True} Boolean;
#endif

#ifndef pointtypes_declare
#define pointtypes_declare

typedef struct Point{
  struct Point *next;
  struct Point *prev;
  struct Point *image;  // pointer to point on image or source plane
  unsigned long id;
  double kappa;        // surface density
  double gamma[2];    // shear
  double dt;          // time delay
  double invmag;     // inverse of magnification
  double *x;
  double gridsize;
  unsigned long head;         // marks beginning of allocated array of points for easy deallocation
  Boolean in_image; // marks if point is in image
  double surface_brightness;

  struct branchstruct{
    struct Point *points;        // pointer to first points in Branch
    unsigned long npoints;
    double center[2];
    int level;
    unsigned long number;
    double boundery_p1[2];
    double boundery_p2[2];
    struct branchstruct *child1;
    struct branchstruct *child2;
    struct branchstruct *brother;
    struct branchstruct *prev;
  } *leaf;
} Point;

typedef struct branchstruct Branch;

#endif
