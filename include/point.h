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

#include <standard.h>
#include "Kist.h"

#ifndef pi
#define pi  3.141593
#endif

#ifndef error_message
#define error_message
#define ERROR_MESSAGE() std::cout << "ERROR: file: " << __FILE__ << " line: " << __LINE__ << std::endl;
#endif

#ifndef bool_declare
#define bool_declare
typedef enum {FALSE, TRUE, MAYBE} Boo;
#endif


//struct branchstruct;
struct Branch;


/** \brief A point on the source or image plane that contains a position and the lensing quantities */

struct Point{
    
  Point();
  struct Point *next;    // pointer to next point in linked list
  struct Point *prev;
  struct Point *image;  // pointer to point on image or source plane
  unsigned long id;
  double x[2];         // the position of the point
  unsigned long head;         // marks beginning of allocated array of points for easy deallocation
  Boo in_image; // marks if point is in image

    
  // redundant information in image and source points
  KappaType kappa;           // surface density
  KappaType gamma[3];        // shear, third component is the rotation quantity that is only non-zero for multi-plane lensing
  double dt;                 // time delay : double implies permanent precision independently from DOUBLE_PRECISION
  KappaType invmag;          // inverse of magnification
    
  double gridsize;           // the size of the most refined grid the point is in
  float surface_brightness;  // the surface brightness at this points

  Branch *leaf;
  bool flag;

  void print();
  void Print();
  
  /// cross product of points 2d positions
  double cross(Point &p){
    return x[0]*p.x[1] - x[1]*p.x[0];
  }
  /// dot product of points in 2d
  double dot(Point &p){
    return x[0]*p.x[0] + x[1]*p.x[1];
  }

  static bool orderX(Point *p1,Point *p2){
    return (p1->x[0] < p2->x[0]);
  }
  static bool orderXrev(Point *p1,Point *p2){
    return (p1->x[0] > p2->x[0]);
  }
  static bool orderY(Point *p1,Point *p2){
    return (p1->x[1] < p2->x[1]);
  }
  static bool orderYrev(Point *p1,Point *p2){
    return (p1->x[1] > p2->x[1]);
  }

};


/// The box representing a branch of a binary tree structure.  Used specifically in TreeStruct for organizing points in the grid.
struct Branch{
	Branch(Point *my_points,unsigned long my_npoints
			  ,double my_boundary_p1[2],double my_boundary_p2[2]
			  ,double my_center[2],int my_level);
	~Branch();

  struct Point *points;        /// pointer to first points in Branch
  Kist<Point>::iterator pointit;       /// Kist iterator pointing to first point in branch
  
  unsigned long npoints;
  double center[2];
  int level;
  unsigned long number;
  double boundary_p1[2];
  double boundary_p2[2];
  Branch *child1;
  Branch *child2;
  Branch *brother;
  Branch *prev;
  /// Marks point as start of a level of refinement
  bool refined;

  void print();

  PosType area(){return (boundary_p2[0]-boundary_p1[0])*(boundary_p2[1]-boundary_p1[1]);}
    
  std::list<Branch *> neighbors;
private:
  static unsigned long countID;
} ;

//typedef struct branchstruct Branch;

/** \brief link list for points, uses the linking pointers within the Point type unlike  Kist */
struct PointList{
  Point *top;
  Point *bottom;
  Point *current;
  unsigned long Npoints;
};

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

/** \brief Transform all points in kist from image to source plane or vis versus.
 * Data type must have a "image" attribute.
 *
void TranformPlanes(Kist<Point> &kist){
  
	if(kist.Nunits() == 0) return;
	kist.MoveToTop();
	do{
		assert(kist.getCurrent());
		assert(kist.getCurrent()->image);
		kist.getCurrent() = kist.getCurrent()->image;
	}while(kist.Down());
  
	return;
}*/

#endif
