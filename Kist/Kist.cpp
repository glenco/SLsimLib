/** \file
 *
 * \brief   routines for linked list of Datas  types
 *
 * Code Name:     kist
 * Programmer:    R Ben Metcalf
 * Last Revised:  Nov, 2005                                   
 * Discription:  replaces linked list programs
 * Comments:                           
 */

#include "slsimlib.h"

/***** Structs *****/

/// constructor
/**Kist<Data>::Kist(
		unsigned long my_blocksize  /// optional argument with size of memory blocks used
		){
	  top=NULL;
	  Number=0;
	  bottom = top;
	  current = top;
	  blocksize = my_blocksize;

	  reserve_top = NULL;
	  Nreserve = 0;

	  return;
}*/


/// Destructor.  Needs to remove all the units.
/*Kist<Data>::~Kist(){

	//Kist<Data>::Empty();

	Unit *units;
	while(heads.size() > 0){
		units = heads.back();
		heads.pop_back();
		delete[] units;
	}
	//units.clear();
	//delete[] units;


	return;
}*/
/// Copy constructor.  Two lists can point to the same data.
/*Kist<Data>::Kist(Kist &a){
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
}*/


void EmptyKist(Kist<Point> * kist){

	kist->Empty();

	return;
}

/**
 * \brief Deallocate memory for all points in list.
 * Note: Does destroy the data points.  Leaves Kist
 * in empty state.
 *
void Kist<Data>::FreeAll(){
	Data *data_t[Number];
	unsigned long Nheads=0,Ndata=0,i,Nunits_t;

	Nunits_t = Number;
	i=0;
	Nheads = 0;
	while(Number > 0){
		data_t[Nheads] = Kist<Data>::TakeOutCurrent();
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
void FreeAllKist(Kist * kist){
	assert(kist);

	kist->FreeAll();

	return;
}*/

// Check state

bool AtTopKist(Kist<Point> * kist){
	assert(kist);

	return kist->AtTop();
}

bool AtBottomKist(Kist<Point> * kist){
	assert(kist);

	return kist->AtBottom();
}

// Insert and remove

void InsertAfterCurrentKist(Kist<Point> * kist,Point *data){
	/* leaves current unchanged */
	assert(kist);

	kist->InsertAfterCurrent(data);
    return;
}
void InsertBeforeCurrentKist(Kist<Point> * kist,Point *data){
	// leaves current unchanged
	assert(kist);

	kist->InsertBeforeCurrent(data);
    return;
}
void SwapCurrentWithBottomKist(Kist<Point> * kist){
	assert(kist);

	kist->SwapCurrentWithBottom();
}


Point *TakeOutCurrentKist(Kist<Point> * kist){
	assert(kist);

	return kist->TakeOutCurrent();
}


/*
 *  Moving through kist
 */

/// Moves jump elements down from current if possible. Negative jump moves up list.
/// Returns false if and of list is reached before jump is finished.

void JumpDownKist(Kist<Point> * kist,int jump){
	assert(kist);

	kist->JumpDown(jump);
	return ;
}
bool MoveDownKist(Kist<Point> * kist){
	assert(kist);

	return kist->Down();
}

bool MoveUpKist(Kist<Point> * kist){
	assert(kist);

	return kist->Up();
}

bool MoveToTopKist(Kist<Point> * kist){
	return kist->MoveToTop();
}

bool MoveToBottomKist(Kist<Point> * kist){
	assert(kist);
	kist->MoveToBottom();
}

void FillKist(Kist<Point> * kist,Point *data_array,unsigned long N){
	assert(kist);
	kist->Fill(data_array,N);
}

/*
void SwapDataInKist(Kist * kist,Unit *p1,Unit *p2){
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

Point *getCurrentKist(Kist<Point> * kist){
	assert(kist);

	return kist->getCurrent();
}

void TranformPlanesKist(Kist<Point> * kist){
	assert(kist);

	kist->TranformPlanes();
	return;
}



