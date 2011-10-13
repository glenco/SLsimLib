/*
 * List.h
 *
 *  Created on: Nov 15, 2010
 *      Author: bmetcalf
 */

#include <point.h>

#ifndef pointlist_declare
#define pointlist_declare

/** \brief link list for points, uses the linking pointers within the Point type unlike  Kist */
typedef struct PointList{
  Point *top;
  Point *bottom;
  Point *current;
  unsigned long Npoints;
} PointList;

typedef struct PointList *ListHndl;

#endif
/***********************************************************
   routines for linked list of points
************************************************************/

ListHndl NewList(void);
Point *NewPoint(double *x,unsigned long id);
void InsertAfterCurrent(ListHndl list,double *x,unsigned long id,Point *image);
void InsertBeforeCurrent(ListHndl list,double *x,unsigned long id,Point *image);
void InsertPointAfterCurrent(ListHndl list,Point *);
void InsertPointBeforeCurrent(ListHndl list,Point *);
Boolean AtTopList(ListHndl list);
Boolean AtBottomList(ListHndl list);
void JumpDownList(ListHndl list,int jump);
Boolean MoveDownList(ListHndl list);
Boolean MoveUpList(ListHndl list);
void ShiftList(ListHndl list);
inline void MoveToTopList(ListHndl list);
inline void MoveToBottomList(ListHndl list);
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
Boolean ArePointsUniqueList(ListHndl list);
Boolean IntersectionList(ListHndl list1,ListHndl list2,ListHndl intersection);
