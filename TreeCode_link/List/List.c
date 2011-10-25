/*
 * List.c
 *
 *  Created on: Nov 15, 2010
 *      Author: bmetcalf
 */
/*#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "List.h"*/

#include <slsimlib.h>

/***********************************************************
   routines for linked list of points
************************************************************/
/** \ingroup ConstructorL2
 *  \brief Constructor for a new list
 */
ListHndl NewList(void){
  PointList *pointlist;

  pointlist=(PointList *) malloc(sizeof(PointList));
  if (!pointlist){
    ERROR_MESSAGE(); fprintf(stderr,"allocation failure in NewList()\n");
    exit(1);
  }
  pointlist->top=NULL;
  pointlist->Npoints=0;
  pointlist->bottom = pointlist->top;
  pointlist->current = pointlist->top;

  return pointlist;
}
/** \ingroup ConstructorL2
 *  \brief Destructor for a new list
 */

void freeList(ListHndl list){

	EmptyList(list);
	free(list);
}

inline bool AtTopList(ListHndl list){
	assert(list);
	if(list->current==list->top) return true;
	else return false;
}
inline bool AtBottomList(ListHndl list){
	assert(list);
	if(list->current==list->bottom) return true;
	else return false;
}


void InsertAfterCurrent(ListHndl list,double *x,unsigned long id,Point *image){
	assert(list);

	Point *point;
  /* leaves current unchanged */

    point=NewPoint(x,id);
    point->image=image;

    if(list->Npoints > 0){
    	assert(list->top);
    	assert(list->bottom);
    	assert(list->current);
      point->prev=list->current;
      point->next=list->current->next;

      if(list->current == list->bottom) list->bottom=point;
      else list->current->next->prev=point;
      list->current->next=point;
    }else{  /* empty list case */
      list->current=point;
      list->top=point;
      list->bottom=point;
    }
    list->Npoints++;
    return;
}

void InsertBeforeCurrent(ListHndl list,double *x,unsigned long id,Point *image){
    Point *point;
  /* leaves current unchanged */

    assert(list);
    point=NewPoint(x,id);
    point->image=image;

    if(list->Npoints > 0){
    	assert(list->top);
    	assert(list->bottom);
    	assert(list->current);

      point->prev=list->current->prev;
      point->next=list->current;
      if(list->current == list->top) list->top=point;
      else list->current->prev->next=point;
      list->current->prev=point;
    }else{  /* empty list case */
      list->current=point;
      list->top=point;
      list->bottom=point;
    }

    list->Npoints++;
    return;
}

void InsertPointAfterCurrent(ListHndl list,Point *point){
 // leaves current unchanged
  // changes only list and links in point
	assert(list);
    if(list->Npoints > 0){
    	assert(list);
    	assert(list->top);
    	assert(list->bottom);
    	assert(list->current);

      point->prev=list->current;
      point->next=list->current->next;

      if(list->current == list->bottom) list->bottom=point;
      else list->current->next->prev=point;
      list->current->next=point;
    }else{  /* empty list case */
      list->current=point;
      list->top=point;
      list->bottom=point;
    }
    list->Npoints++;
    return;
}

void InsertPointBeforeCurrent(ListHndl list,Point *point){
	assert(list);

	/* leaves current unchanged */

    if(list->Npoints > 0){
      point->prev=list->current->prev;
      point->next=list->current;
      if(list->current == list->top) list->top=point;
      else list->current->prev->next=point;
      list->current->prev=point;
    }else{
    	assert(list->top);
    	assert(list->bottom);
    	assert(list->current);

      list->current=point;
      list->top=point;
      list->bottom=point;
      point->prev=NULL;
      point->next=NULL;
    }
    list->Npoints++;
    return;
}

void MoveCurrentToBottom(ListHndl list){
	assert(list);
	assert(list->top);
	assert(list->bottom);
	assert(list->current);

	Point *point,*point_current;
	// leaves current one above former current

	point=TakeOutCurrent(list);
	point_current=list->current;
	MoveToBottomList(list);
	InsertPointAfterCurrent(list,point);
	list->current=point_current;
}

/**
 *  takes out current point and set current to point previous */
/* Except at top where current is set to new top */
/* returns pointer to removed point */
Point *TakeOutCurrent(ListHndl list){

    Point *point;

    if(list == NULL || list->Npoints <= 0) return NULL;
    assert(list->current);

    point=list->current;

    if(list->top == list->bottom){  /* last point */
      list->current=NULL;
      list->top=NULL;
      list->bottom=NULL;
    }else if(list->current==list->top){
      list->top=list->top->next;
      list->current=list->top;
      list->top->prev=NULL;
    } else if(list->current==list->bottom){
      list->bottom=list->bottom->prev;
      list->current=list->bottom;
      list->bottom->next=NULL;
    }else{
      list->current->prev->next=list->current->next;
      list->current->next->prev=list->current->prev;
      list->current=list->current->prev;
    }

    list->Npoints--;

    point->prev=point->next=NULL;

    return point;
}

void MergeLists(ListHndl list1,ListHndl list2){
	/* list 1 is made into the union of lists 1 & 2
	 *   list 2 remains a sublist
	 */

	if(list2->Npoints==0) return;
	if(list1->Npoints == 0){
		list1->Npoints=list2->Npoints;
		list1->top=list2->top;
		list1->bottom=list2->bottom;
		list1->current=list1->current;
		return;
	}

	list1->bottom->next=list2->top;
	list2->top->prev=list1->bottom;
	list1->bottom=list2->bottom;
	list1->Npoints+=list2->Npoints;
	return;
}

void InsertListAfterCurrent(ListHndl list1,ListHndl list2){
	/* inserts list2 into list1
	 *  leaves list2 as a sublist
	 *  leaves currents unchanged
	 */
	Point *point;

	if(list2->Npoints==0) return;
	if(list1->Npoints == 0){
		list1->Npoints=list2->Npoints;
		list1->top=list2->top;
		list1->bottom=list2->bottom;
		list1->current=list1->current;
		return;
	}

	point=list1->current->next;
	list1->current->next=list2->top;
	list2->top->prev=list1->current;

	if(list1->current==list1->bottom){
		list1->bottom=list2->bottom;
	}else{
		point->prev=list2->bottom;
		list2->bottom->next=point;
	}

	list1->Npoints+=list2->Npoints;
	return;
}
void InsertListBeforeCurrent(ListHndl list1,ListHndl list2){
	/* inserts list2 into list1
	 *  leaves list2 as a sublist
	 *  leaves currents unchanged
	 */
	Point *point;

	if(list2->Npoints==0) return;
	if(list1->Npoints == 0){
		list1->Npoints=list2->Npoints;
		list1->top=list2->top;
		list1->bottom=list2->bottom;
		list1->current=list1->current;
		return;
	}

	point=list1->current->prev;
	list1->current->prev=list2->bottom;
	list2->bottom->next=list1->current;

	if(list1->current==list1->top){
		list1->top=list2->top;
	}else{
		point->next=list2->top;
		list2->top->prev=point;
	}

	list1->Npoints+=list2->Npoints;
	return;
}

void EmptyList(ListHndl list){
	/* This function should properly release the memory for all the
	 * points in a list leaving the list with NULL pointers
	 *  points need to have been allocated in blocks of 1
	 */

	if(list == NULL || list->Npoints==0) return;

	assert(list->top);
	assert(list->bottom);
	assert(list->current);

	Point **point;
	unsigned long blocks=0,i=0,Nfreepoints=0;

	MoveToTopList(list);
	do{
	  if(list->current->head > 0){ ++blocks; Nfreepoints += list->current->head;}
	}while(MoveDownList(list));

	assert(blocks <= list->Npoints);
	assert(Nfreepoints == list->Npoints);

	point=(Point **)malloc(blocks*sizeof(Point*));

	MoveToTopList(list);
	do{
	  if(list->current->head > 0){ point[i] = list->current; ++i;}
	}while(MoveDownList(list));

	for(i=0;i<blocks;++i) free(point[i]);
	free(point);

	list->Npoints = 0;
	list->top = NULL;
	list->bottom = NULL;
	list->current = NULL;

  /*
  Point *point;

  while(list->Npoints > 0 ){
    point = TakeOutCurrent(list);
    assert(point->head == 1);

    free(point);
  }
*/
}

void JumpDownList(ListHndl list,int jump){
  int i;

  if(jump > 0) for(i=0;i<jump;++i) MoveDownList(list);
  if(jump < 0) for(i=0;i<abs(jump);++i) MoveUpList(list);
}

bool MoveDownList(ListHndl list){

	if(list->Npoints == 0) return false;
	if(list->current==list->bottom) return false;
	list->current=list->current->next;

	return true;
}

bool MoveUpList(ListHndl list){

	if(list->Npoints == 0) return false;
	if(list->current==list->top) return false;
	list->current=list->current->prev;

	return true;
}

void ShiftList(ListHndl list){
	// cyclic shift of list
	// to make current top

	// link into ring
	list->top->prev=list->bottom;
	list->bottom->next=list->top;

	list->top=list->current;
	list->bottom=list->top->prev;

	// break ring
	list->bottom->next=NULL;
	list->top->prev=NULL;
}

inline void MoveToTopList(ListHndl list){
  list->current=list->top;
}
inline void MoveToBottomList(ListHndl list){
  list->current=list->bottom;
}

void FillList(ListHndl list,double **x,unsigned long N
	      ,unsigned long idmin){
  unsigned long i;
  /* add N points to to end of list */
  /* id numbers are given in order from idmin */
  /* this is used to initialize list */

  printf("Inserting  %li points after %li x= %e %e\n",N,list->current->id
	 ,list->current->x[0],list->current->x[1]);
  MoveToBottomList(list);
  InsertAfterCurrent(list,x[0],idmin,NULL);
  for(i=1;i<N;++i){
    MoveDownList(list);
    InsertAfterCurrent(list,x[i],i+idmin,NULL);
  }
}

void SwapPointsInList(ListHndl list,Point *p1,Point *p2){
  Point pt;

  if(p1==p2) return;

  pt.prev=p1->prev;
  pt.next=p1->next;

  if(p1->prev != NULL) p1->prev->next = p2;
  if(p1->next != NULL) p1->next->prev = p2;

  if(p2->prev != NULL) p2->prev->next = p1;
  if(p2->next != NULL) p2->next->prev = p1;

  p1->prev=p2->prev;
  p1->next=p2->next;

  p2->prev=pt.prev;
  p2->next=pt.next;

  if(list->top == p1) list->top = p2;
  else if(list->top == p2) list->top = p1;
  if(list->bottom == p1) list->bottom = p2;
  else if(list->bottom == p2) list->bottom = p1;
}

void PrintList(ListHndl list){
  unsigned long i;
  Point *placemark;

  placemark=list->current;
  MoveToTopList(list);
//  printf("%i points in list\n",list->Npoints);
  printf("%li\n",list->Npoints);

  for(i=0;i<list->Npoints;++i){
    printf("%li  %li  %e %e %e\n",i,list->current->id,list->current->x[0],list->current->x[1]
                                  ,list->current->gridsize);
    MoveDownList(list);
  }

  list->current=placemark;
}
