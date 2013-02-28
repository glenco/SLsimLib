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

#include "slsimlib.h"

/***** Structs *****/

/// constructor
Kist::Kist(
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
}


/// Destructor.  Needs to remove all the units.
Kist::~Kist(){

	//Kist::Empty();

	Unit *units;
	while(heads.size() > 0){
		units = heads.back();
		heads.pop_back();
		delete[] units;
	}
	//units.clear();
	//delete[] units;


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
/// Take a Unit out of the reservoir, expand reservoir if necessary
Unit * Kist::pop_from_reserve(){

	if(Nreserve <= 1){
		// expand reservoir of Units for later use
		  Unit *units = new Unit[blocksize];
		  for(unsigned long i=0;i<blocksize-1;++i) units[i].next = &units[i+1];
		  units[blocksize-1].next = NULL;

		  heads.push_back(units);  // keep list of first point in each block of memory
		  if(Nreserve == 1) reserve_top->next = units;
		  else reserve_top = units;
		  Nreserve += blocksize;
	}

	Unit *unit = reserve_top;
	reserve_top = unit->next;
	--Nreserve;

	return unit;
}
/// put a Unit into the reservoir
void Kist::push_to_reserve(Unit *unit){
	assert(reserve_top);
	assert(Nreserve > 0);

	unit->next = reserve_top;
	reserve_top = unit;
	++Nreserve;
}
/** \ingroup ImageFindingL2
 * \brief Removes all elements from list without destroy the data.
 */
void Kist::Empty(){

	while(Number > 0) TakeOutCurrent();
	//Number = 0;
	//top = bottom = current = NULL;
	return;
}

void EmptyKist(KistHndl kist){

	kist->Empty();

	return;
}

/**
 * \brief Deallocate memory for all points in list.
 * Note: Does destroy the data points.  Leaves Kist
 * in empty state.
 *
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
}*/

// Check state
bool Kist::AtTop(){
	if(current==top) return true;
	return false;
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
    assert(data->x);
    assert(data->gridsize >= 0);

    //Unit *unit = (Unit *)malloc(sizeof(Unit));
    Unit *unit = pop_from_reserve();

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

	assert(data);
    assert(data->x);
    assert(data->gridsize >= 0);

    //Unit *unit = (Unit *)malloc(sizeof(Unit));
    Unit *unit = pop_from_reserve();

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
    
    push_to_reserve(unit);
    //std::free(unit);
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
/// Move down the list. Returns false when at the bottom of the list.
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
/// Move up the list. Returns false when at the top of the list.
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

/// Moves to the top.  Returns false if the kist is empty.
bool Kist::MoveToTop(){
	if(!Number)return false;
	current=top;
	return true;
}
bool MoveToTopKist(KistHndl kist){
	return kist->MoveToTop();
}

/// Moves to the top.  Returns false if the kist is empty.
bool Kist::MoveToBottom(){
	if(!Number)return false;
	current=bottom;
	return true;
}
bool MoveToBottomKist(KistHndl kist){
	assert(kist);
	kist->MoveToBottom();
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
	if(Number==0) return NULL;
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
	std::cout << Nunits() << std::endl;
	MoveToTop();
	do{
		std::cout << getCurrent()->x[0] <<  "  " << getCurrent()->x[1] << "  " << getCurrent()->gridsize << std::endl;
	}while(Down());
}
/**
 * \brief return true if all the points in the kist are unique.
 */
bool Kist::AreDataUnique(){

	if(Nunits() < 2) return true;
	Unit *unit = current;

	unsigned long i,j;
	Data **data = new Data*[Nunits()];

	MoveToTop();
	for(i=0;i<Nunits();++i,Down()){
		data[i]=getCurrent();
	}

	MoveToTop();
	for(i=0;i<Nunits();++i,Down()){
		for(j=0;j<i;++j){
			if( getCurrent() == data[j] ){
				current = unit;
				delete[] data;
				return false;
			}
		}
	}

	delete[] data;
	current = unit;
	return true;
}



