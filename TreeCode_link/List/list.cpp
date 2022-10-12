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
/*
void PointList::InsertAfterCurrent(iterator &current,PosType *x,unsigned long id,Point *image){

	Point *point;
  // leaves current unchanged

  point=NewPoint(x,id);
  point->image=image;

  InsertPointAfterCurrent(current,point);
    return;
}
*/
void PointList::InsertPointAfterCurrent(PointList::iterator &it,Point *point){
  // leaves current unchanged
  // changes only list and links in point
  
  if(Npoints > 0){
    assert(top);
    assert(bottom);
    
    point->prev=*it;
    point->next=(*it)->next;
    
    if(*it == bottom) bottom=point;
    else (*it)->next->prev=point;
    (*it)->next=point;
  }else{  /* empty list case */
    it.current=point;
    top=point;
    bottom=point;
    point->prev = point->next = NULL;
    it.current = point;
  }
  Npoints++;
  return;
}
/*
void PointList::InsertBeforeCurrent(iterator &current,PosType *x,unsigned long id,Point *image){
    Point *point;
  // leaves current unchanged

    point=NewPoint(x,id);
    point->image=image;

    InsertPointBeforeCurrent(current,point);
    return;
}
*/

void PointList::InsertPointBeforeCurrent(iterator &it,Point *point){
  
  if(Npoints > 0){
    assert(top);
    assert(bottom);
    
    point->prev = (*it)->prev;
    point->next = *it;
    if(*it == top) top=point;
    else (*it)->prev->next=point;
    (*it)->prev=point;
  }else{  /* empty list case */
    it.current=point;
    top=point;
    bottom=point;
    point->prev = point->next = NULL;
    it.current = point;
  }
  
  Npoints++;

  return;
}

void PointList::MoveCurrentToBottom(iterator &it){

	assert(top);
	assert(bottom);

	Point *point,*point_current;
	// leaves current one above former current

	point=TakeOutCurrent(it);
	point_current=*it;
  it.current = bottom;
	InsertPointAfterCurrent(it,point);
	it.current=point_current;
}

/**
 *  takes out current point and set current to point previous */
/* Except at top where current is set to new top */
/* returns pointer to removed point */
Point *PointList::TakeOutCurrent(iterator &it){

    Point *point;

    if(Npoints <= 0) return NULL;
    assert(*it);

    point = *it;

    if(top == bottom){  /* last point */
      it.current=NULL;
      top=NULL;
      bottom=NULL;
    }else if(*it==top){
      top=top->next;
      it.current = top;
      top->prev=NULL;
    } else if(*it == bottom){
      bottom = bottom->prev;
      it.current = bottom;
      bottom->next=NULL;
    }else{
      (*it)->prev->next=(*it)->next;
      (*it)->next->prev=(*it)->prev;
      it.current = (*it)->prev;
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

/* This function should empty the PointList without deleting the points
 */
void PointList::EmptyList(){

	if(Npoints==0) return;
/*
	assert(top);
	assert(bottom);

	//Point **point;
  std::vector<Point *> point;
	unsigned long blocks=0,i=0,Nfreepoints=0;

  iterator it;
  it.current = top;

	do{
	  if((*it)->head > 0){ ++blocks; Nfreepoints += (*it)->head;}
	}while(it--);

	assert(blocks <= Npoints);
	assert(Nfreepoints == Npoints);

	//point=(Point **)malloc(blocks*sizeof(Point*));
  point.resize(blocks);
  
  it.current = top;
	do{
	  if((*it)->head > 0){ point[i] = *it; ++i;}
	}while(it--);

	for(i=0;i<blocks;++i) free(point[i]);
	//free(point);
*/
	Npoints = 0;
	top = NULL;
	bottom = NULL;
}


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

/*
void PointList::FillList(PosType **x,unsigned long N
	      ,unsigned long idmin){
  unsigned long i;
  // add N points to to end of list
  // id numbers are given in order from idmin
  // this is used to initialize list


  iterator it;
  it.current = bottom;
  InsertAfterCurrent(it,x[0],idmin,NULL);
  for(i=1;i<N;++i){
    --it;
    InsertAfterCurrent(it,x[i],i+idmin,NULL);
  }
}
*/

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

  iterator it;
  it.current = top;
  
//  std::printf("%i points in list\n",list->Npoints);
  std::printf("%li\n",Npoints);

  for(i=0;i<Npoints;++i){
    std::printf("%li  %li  %e %e %e\n",i,(*it)->id,(*it)->x[0],(*it)->x[1]
                                  ,(*it)->gridsize);
    --it;
  }

}

std::ostream &operator<<(std::ostream &os, Point_2d const &p) {
  return os << p.x[0] << " " << p.x[1];
}


std::ostream &operator<<(std::ostream &os, CritType const &p) {
  if(p == CritType::ND) return os << "not determined";
  if(p == CritType::pseudo) return os << "Pseudo";
  if(p == CritType::tangential) return os << "Tangential";
  return os << "Radial";
}


