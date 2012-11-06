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
  float kappa;        // surface density
  float gamma[3];    // shear, third component is the rotation quantity that is only non-zero for multi-plane lensing
  float dt;          // time delay
  float invmag;     // inverse of magnification
  double gridsize;   // the size of the most refined grid the point is in
  float surface_brightness;  // the surface brightness at this points

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


/** \brief link list for points, uses the linking pointers within the Point type unlike  Kist */
typedef struct PointList{
  Point *top;
  Point *bottom;
  Point *current;
  unsigned long Npoints;
} PointList;

typedef struct PointList *ListHndl;
//bool AtTopList(ListHndl list);
//bool AtBottomList(ListHndl list);

inline bool AtTopList(ListHndl list){
	assert(list);
	if(list->current==list->top) return true;
	else return false;
};
inline bool AtBottomList(ListHndl list){
	assert(list);
	if(list->current==list->bottom) return true;
	else return false;
};
//inline void MoveToTopList(ListHndl list);
//inline void MoveToBottomList(ListHndl list);
inline void MoveToTopList(ListHndl list){
  list->current=list->top;
};
inline void MoveToBottomList(ListHndl list){
  list->current=list->bottom;
};

/***********************************************************
   routines for linked list of points
************************************************************/

ListHndl NewList(void);
Point *NewPoint(double *x,unsigned long id);
void InsertAfterCurrent(ListHndl list,double *x,unsigned long id,Point *image);
void InsertBeforeCurrent(ListHndl list,double *x,unsigned long id,Point *image);
void InsertPointAfterCurrent(ListHndl list,Point *);
void InsertPointBeforeCurrent(ListHndl list,Point *);

void JumpDownList(ListHndl list,int jump);
bool MoveDownList(ListHndl list);
bool MoveUpList(ListHndl list);
void ShiftList(ListHndl list);


void FillList(ListHndl list,double **x,unsigned long N
	      ,unsigned long idmin);
void SwapPointsInList(ListHndl list,Point *p1,Point *p2);
void PrintList(ListHndl list);
Point *sortList(long n, double arr[],ListHndl list,Point *firstpoint);
void MoveCurrentToBottom(ListHndl list);
Point *TakeOutCurrent(ListHndl list);
void MergeLists(ListHndl list1,ListHndl list2);
void InsertListAfterCurrent(ListHndl list1,ListHndl list2);
void InsertListBeforeCurrent(ListHndl list1,ListHndl list2);
void EmptyList(ListHndl list);
void UnionList(ListHndl list1,ListHndl list2);
bool ArePointsUniqueList(ListHndl list);
bool IntersectionList(ListHndl list1,ListHndl list2,ListHndl intersection);


#endif
