/** \file
 *
 * \brief   routines for linked list of Datas
 *
 * Code Name:     kist
 * Programmer:    R Ben Metcalf
 * Last Revised:  Nov, 2005                                   
 * Discription:  replaces linked list programs
 * Comments:                           
 */

/*#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "Kist.h"*/

#include <slsimlib.h>

/***** Structs *****/

/// constructor
Kist::Kist(){
	  top=NULL;
	  Number=0;
	  bottom = top;
	  current = top;

	  return;
}

/// Destructor.  Needs to remove all the units.
Kist::~Kist(){

	Kist::Empty();

	return;
}
/// Copy constructor.  Two lists can point to the same data.
Kist::Kist(Kist &a){
	Unit *current;

	top=NULL;
	Number=0;
	bottom = top;
	current = top;

	current = a.current;
	a.MoveToTop();
	do{
		InsertAfterCurrent(a.getCurrent());
		Down();
	}while(a.Down());
	assert(a.Nunits() == Nunits());
	a.current = current;
}
/*
 * This is the obsolete C style constructor.
KistHndl NewKist(void){
  KistHndl kist;

  kist=(Kist *) malloc(sizeof(Kist));
  if (!kist){
    std::fprintf(stderr,"allocation failure in NewKist()\n");
    exit(1);
  }
  kist->top=NULL;
  kist->Nunits=0;
  kist->bottom = kist->top;
  kist->current = kist->top;

  return kist;
}
*/

/** \ingroup ImageFindingL2
 * \brief Removes all elements from list without destroy the data.
 */
void Kist::Empty(){

	while(Number > 0) TakeOutCurrent();
	return;
}

void EmptyKist(KistHndl kist){

	kist->Empty();

	return;
}

/* \ingroup ConstructorL2
*  deallocate memory for kist
*  Note: Does not destroy data.  A handle must still point to the data.
 * *
void freeKist(KistHndl kist){
	EmptyKist(kist);
	free(kist);
	return;
}*/

/**
 * \brief Deallocate memory for all points in list.
 * Note: Does destroy the data points.  Leaves Kist
 * in empty state.
 */
void Kist::FreeAll(){
	Data *data_t[Number];
	unsigned long Nheads=0,Ndata=0,i,Nunits_t;

	Nunits_t = Number;
	i=0;
	Nheads = 0;
	while(Number > 0){
		data_t[Nheads] = Kist::TakeOutCurrent();
		Ndata += data_t[i]->head;
		if(data_t[Nheads]->head) ++Nheads;
	}

	for(i = 0; i < Nheads ; ++i){
		std::free(data_t[i]->x);
		std::free(data_t[i]);
	}

	if(Number != Ndata){
		ERROR_MESSAGE();
		std::printf("FreeAllKist freed not all of data in kist\n");
		exit(0);
	}
	return;
}
void FreeAllKist(KistHndl kist){
	assert(kist);

	kist->FreeAll();

	return;
}

// Check state
bool Kist::AtTop(){
	if(current==top) return true;
	else return false;
}
bool AtTopKist(KistHndl kist){
	assert(kist);

	return kist->AtTop();
}

bool Kist::AtBottom(){
	if(current==bottom) return true;
	else return false;
}
bool AtBottomKist(KistHndl kist){
	assert(kist);

	return kist->AtBottom();
}

// Insert and remove

/// Insert data into kist.  Leaves current unchanged.
void Kist::InsertAfterCurrent(Data *data){
	assert(data);

    Unit *unit = (Unit *)malloc(sizeof(Unit));

    assert(unit);

    unit->data = data;

    if(Number > 0){
      assert(current);

      unit->prev = current;
      unit->next = current->next;

      if(current == bottom) bottom = unit;
      else current->next->prev=unit;
      current->next=unit;
    }else{  // empty kist case
      current=unit;
      top=unit;
      bottom=unit;

      unit->prev = NULL;
      unit->next = NULL;
    }

    Number++;
    return;
}
void InsertAfterCurrentKist(KistHndl kist,Data *data){
	/* leaves current unchanged */
	assert(kist);

	kist->InsertAfterCurrent(data);
    return;
}
/// Insert data into kist.  Leaves current unchanged.
void Kist::InsertBeforeCurrent(Data *data){
	// leaves current unchanged
	assert(data);

	Unit *unit;

	unit=(Unit *)malloc(sizeof(Unit));
    assert(unit);

	unit->data=data;

    if(Number > 0){
      assert(current);

      unit->prev=current->prev;
      unit->next=current;
      if(current == top) top=unit;
      else current->prev->next=unit;
      current->prev=unit;
    }else{  /* empty kist case */
      current=unit;
      top=unit;
      bottom=unit;

      unit->prev = NULL;
      unit->next = NULL;
    }

    Number++;
    return;
}
void InsertBeforeCurrentKist(KistHndl kist,Data *data){
	// leaves current unchanged
	assert(kist);

	kist->InsertBeforeCurrent(data);
    return;
}
/**
 * Swaps current data with bottom data leaving current one above former current.
 */
void Kist::SwapCurrentWithBottom(){

	Data *data;

	data=current->data;
	current->data=bottom->data;
	bottom->data=data;
	Up();
}
void SwapCurrentWithBottomKist(KistHndl kist){
	assert(kist);

	kist->SwapCurrentWithBottom();
}
/**
 * Moves current to the bottom of the kist.  Current
 * is left at the bottom
 */
void Kist::MoveCurrentToBottom(){

	Data *data = TakeOutCurrent();
	MoveToBottom();
	InsertAfterCurrent(data);
	Down();
}

/**
 * \brief Takes out current data and set current to data previous
* except at top where current is set to new top.
* Returns pointer to removed data.
*/
Data *Kist::TakeOutCurrent(){

    Data *data;
    Unit *unit;

    if( Number <= 0) return NULL;

    assert(current);
    assert(top);
    assert(bottom);

    data = current->data;
    unit = current;

    if(top == bottom){  /* last data */
      current=NULL;
      top=NULL;
      bottom=NULL;
    }else if(current==top){
      top=top->next;
      current=top;
      top->prev=NULL;
    } else if(current==bottom){
      bottom=bottom->prev;
      current=bottom;
      bottom->next=NULL;
    }else{
      current->prev->next=current->next;
      current->next->prev=current->prev;
      current=current->prev;
    }
    
    std::free(unit);
    Number--;

    return data;
}
Data *TakeOutCurrentKist(KistHndl kist){
	assert(kist);

	return kist->TakeOutCurrent();
}


/*
 *  Moving through kist
 */

/// Moves jump elements down from current if possible. Negative jump moves up list.
/// Returns false if and of list is reached before jump is finished.
bool Kist::JumpDown(int jump){
	int i;
	bool ans;

	if(jump > 0) for(i=0;i<jump;++i) ans=Down();
	if(jump < 0) for(i=0;i<abs(jump);++i) ans=Up();
	return ans;
}
void JumpDownKist(KistHndl kist,int jump){
	assert(kist);

	kist->JumpDown(jump);
	return ;
}
/// Move current down the list. Returns false when at the bottom of the list.
bool Kist::Down(){

	if(Number == 0) return false;
	if(current==bottom) return false;
	current=current->next;

	return true;
}
bool MoveDownKist(KistHndl kist){
	assert(kist);

	return kist->Down();
}
/// Move current up the list. Returns false when at the top of the list.
bool Kist::Up(){

	if(Number == 0) return false;
	if(current==top) return false;
	current=current->prev;

	return true;
}
bool MoveUpKist(KistHndl kist){
	assert(kist);

	return kist->Up();
}


void Kist::MoveToTop(){
	current=top;
}
void MoveToTopKist(KistHndl kist){
	assert(kist);
	kist->MoveToTop();
}

void Kist::MoveToBottom(){
	current=bottom;
}
void MoveToBottomKist(KistHndl kist){
	assert(kist);
	kist->SwapCurrentWithBottom();
}

/** \brief Put an array of data into a kist.
 * */
void Kist::Fill(Data *data_array,unsigned long N){
  unsigned long i;

   for(i=1;i<N;++i){
    InsertAfterCurrent(&(data_array[i]));
    Kist::Down();
  }
}
void FillKist(KistHndl kist,Data *data_array,unsigned long N){
	assert(kist);

	kist->Fill(data_array,N);
}

/*
void SwapDataInKist(KistHndl kist,Unit *p1,Unit *p2){
	assert(kist);
	assert(p1);
	assert(p2);

	Data *data;

	data=p1->data;
	p1->data=p2->data;
	p2->data=data;

	return;
 }
*/
/*********************************/
/*  data extraction routines */
/*********************************/
/// Return pointer to data in current element.
Data * Kist::getCurrent(){
	assert(current);
	assert(current->data);

	return current->data;
}
Data *getCurrentKist(KistHndl kist){
	assert(kist);

	return kist->getCurrent();
}

/*********************************
 * specific to image points
 * ******************************/
/**
 *\brief Transform all points in kist from image to source plane or vis versus
 */
void Kist::TranformPlanes(){

	if(Number == 0) return;
	MoveToTop();
	do{
		assert(current->data);
		assert(current->data->image);
		current->data = current->data->image;
	}while(Down());

	return;
}
void TranformPlanesKist(KistHndl kist){
	assert(kist);

	kist->TranformPlanes();
	return;
}

void Kist::Print(){
	cout << Nunits() << endl;
	MoveToTop();
	do{
		cout << getCurrent()->x[0] <<  "  " << getCurrent()->x[1] << "  " << getCurrent()->gridsize << endl;
	}while(Down());
}
/*
bool AreDataUniqueKist(KistHndl kist){
	assert(kist);
	if(kist->Nunits() < 2) return true;

	unsigned long i,j;
	KistHndl tmpkist = NewKist();

	// clone kist
	//  the two lists use the same units as well
	//  as data.
	tmpkist->Nunits() = kist->Nunits();
	tmpkist->top = kist->top;
	tmpkist->bottom = kist->bottom;

	MoveToTopKist(kist);
	for(i=0;i<kist->Nunits();++i,MoveDownKist(kist)){
		MoveToTopKist(tmpkist);
		for(j=0;j<i;++j,MoveDownKist(tmpkist)){
			if( getCurrentKist(kist) == getCurrentKist(tmpkist) ){
				free(tmpkist);
				return false;
			}
		}
	}

	free(tmpkist);
	return true;
}
*/
