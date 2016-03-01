/*
 * List.c
 *
 *  Created on: Nov 15, 2010
 *      Author: bmetcalf
 */
#include "slsimlib.h"

/***********************************************************
   routines for linked list of points
************************************************************/

void PointList::InsertAfterCurrent(iterator &current,PosType *x,unsigned long id,Point *image){

	Point *point;
  /* leaves current unchanged */

  point=NewPoint(x,id);
  point->image=image;

  InsertPointAfterCurrent(current,point);
    return;
}

void PointList::InsertPointAfterCurrent(iterator &current,Point *point){
  // leaves current unchanged
  // changes only list and links in point
  
  if(Npoints > 0){
    assert(top);
    assert(bottom);
    
    point->prev=*current;
    point->next=(*current)->next;
    
    if(*current == bottom) bottom=point;
    else (*current)->next->prev=point;
    (*current)->next=point;
  }else{  /* empty list case */
    current=point;
    top=point;
    bottom=point;
    point->prev = point->next = NULL;
    current = point;
  }
  Npoints++;
  return;
}

void PointList::InsertBeforeCurrent(iterator &current,PosType *x,unsigned long id,Point *image){
    Point *point;
  /* leaves current unchanged */

    point=NewPoint(x,id);
    point->image=image;

  InsertPointBeforeCurrent(current,point);
    return;
}

void PointList::InsertPointBeforeCurrent(iterator &current,Point *point){
  
  if(Npoints > 0){
    assert(top);
    assert(bottom);
    
    point->prev = (*current)->prev;
    point->next = *current;
    if(*current == top) top=point;
    else (*current)->prev->next=point;
    (*current)->prev=point;
  }else{  /* empty list case */
    current=point;
    top=point;
    bottom=point;
    point->prev = point->next = NULL;
    current = point;
  }
  
  Npoints++;

  return;
}

void PointList::MoveCurrentToBottom(iterator &current){

	assert(top);
	assert(bottom);

	Point *point,*point_current;
	// leaves current one above former current

	point=TakeOutCurrent(current);
	point_current=*current;
  current = bottom;
	InsertPointAfterCurrent(current,point);
	current=point_current;
}

/**
 *  takes out current point and set current to point previous */
/* Except at top where current is set to new top */
/* returns pointer to removed point */
Point *PointList::TakeOutCurrent(iterator &current){

    Point *point;

    if(Npoints <= 0) return NULL;
    assert(*current);

    point = *current;

    if(top == bottom){  /* last point */
      current=NULL;
      top=NULL;
      bottom=NULL;
    }else if(*current==top){
      top=top->next;
      current = top;
      top->prev=NULL;
    } else if(*current == bottom){
      bottom = bottom->prev;
      current = bottom;
      bottom->next=NULL;
    }else{
      (*current)->prev->next=(*current)->next;
      (*current)->next->prev=(*current)->prev;
      current = (*current)->prev;
    }

    Npoints--;

    point->prev=point->next=NULL;

    return point;
}

void PointList::MergeLists(ListHndl list2){
	/* list 1 is made into the union of lists 1 & 2
	 *   list 2 remains a sublist
	 */

	if(list2->Npoints==0) return;
	if(Npoints == 0){
		Npoints=list2->Npoints;
		top=list2->top;
		bottom=list2->bottom;
		return;
	}

	bottom->next=list2->top;
	list2->top->prev=bottom;
	bottom=list2->bottom;
	Npoints += list2->Npoints;
	return;
}

void PointList::InsertListAfterCurrent(iterator &current,ListHndl list2){
	/* inserts list2 into list1
	 *  leaves list2 as a sublist
	 *  leaves currents unchanged
	 */
	Point *point;

	if(list2->Npoints==0) return;
  
	if(Npoints == 0){
		Npoints=list2->Npoints;
		top=list2->top;
		bottom=list2->bottom;
		return;
	}

	point = (*current)->next;
	(*current)->next = list2->top;
	list2->top->prev = *current;

	if(*current == bottom){
		bottom=list2->bottom;
	}else{
		point->prev=list2->bottom;
		list2->bottom->next=point;
	}

	Npoints += list2->Npoints;
	return;
}
void PointList::InsertListBeforeCurrent(iterator &current,ListHndl list2){
	/* inserts list2 into list
	 *  leaves list2 as a sublist
	 *  leaves currents unchanged
	 */
	Point *point;

	if(list2->Npoints==0) return;
	if(Npoints == 0){
		Npoints=list2->Npoints;
		top=list2->top;
		bottom=list2->bottom;
		return;
	}

	point = (*current)->prev;
	(*current)->prev=list2->bottom;
	list2->bottom->next = *current;

	if( *current == top){
		top=list2->top;
	}else{
		point->next=list2->top;
		list2->top->prev=point;
	}

	Npoints += list2->Npoints;
	return;
}

/* This function should properly release the memory for all the
 * points in a list leaving the list with NULL pointers
 *  points need to have been allocated in blocks of 1
 */
void PointList::EmptyList(){

	if(Npoints==0) return;

	assert(top);
	assert(bottom);

	Point **point;
	unsigned long blocks=0,i=0,Nfreepoints=0;

  iterator current(top);

	do{
	  if((*current)->head > 0){ ++blocks; Nfreepoints += (*current)->head;}
	}while(current--);

	assert(blocks <= Npoints);
	assert(Nfreepoints == Npoints);

	point=(Point **)malloc(blocks*sizeof(Point*));

  current = top;
	do{
	  if((*current)->head > 0){ point[i] = *current; ++i;}
	}while(current--);

	for(i=0;i<blocks;++i) free(point[i]);
	free(point);

	Npoints = 0;
	top = NULL;
	bottom = NULL;

  /*
  Point *point;

  while(list->Npoints > 0 ){
    point = TakeOutCurrent(list);
    assert(point->head == 1);

    free(point);
  }
*/
}


/*bool MoveDownList(ListHndl list){

	if(list->Npoints == 0) return false;
	if(list->current==list->bottom) return false;
	list->current=list->current->next;

	return true;
}*/


void PointList::ShiftList(iterator &current){
	// cyclic shift of list
	// to make current top

	// link into ring
	top->prev=bottom;
	bottom->next=top;

	top= *current;
	bottom=top->prev;

	// break ring
	bottom->next=NULL;
	top->prev=NULL;
}

void PointList::FillList(PosType **x,unsigned long N
	      ,unsigned long idmin){
  unsigned long i;
  /* add N points to to end of list */
  /* id numbers are given in order from idmin */
  /* this is used to initialize list */


  iterator current(bottom);
  InsertAfterCurrent(current,x[0],idmin,NULL);
  for(i=1;i<N;++i){
    --current;
    InsertAfterCurrent(current,x[i],i+idmin,NULL);
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

  if(list->Top() == p1) list->setTop(p2);
  else if(list->Top() == p2) list->setTop(p1);
  if(list->Bottom() == p1) list->setBottom(p2);
  else if(list->Bottom() == p2) list->setBottom(p1);
}

void PointList::PrintList(){
  unsigned long i;

  iterator current(top);

//  std::printf("%i points in list\n",list->Npoints);
  std::printf("%li\n",Npoints);

  for(i=0;i<Npoints;++i){
    std::printf("%li  %li  %e %e %e\n",i,(*current)->id,(*current)->x[0],(*current)->x[1]
                                  ,(*current)->gridsize);
    --current;
  }

}

std::ostream &operator<<(std::ostream &os, Point_2d const &p) {
  return os << p.x[0] << " " << p.x[1];
}

std::ostream &operator<<(std::ostream &os, Point_3d const &p) {
  return os << p.x[0] << " " << p.x[1] << " " << p.x[2];
}

std::ostream &operator<<(std::ostream &os, CritType const &p) {
  if(p == ND) return os << "not determined";
  if(p == pseudo) return os << "Pseudo";
  if(p == tangential) return os << "Tangential";
  return os << "Radial";
}


