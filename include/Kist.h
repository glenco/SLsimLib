/*
 * Kist.h
 *
 *  Created on: Feb 2, 2010
 *      Author: R.B. Metcalf
 */

#ifndef kisttypes_declare
#define kisttypes_declare

//#include "pointlist.h"
#include "point.h"

// Used as internal container in Kist
template <class T>
struct Unit{
	T * data;
	struct Unit<T> *next;
	struct Unit<T> *prev;
};// Unit;

//typedef struct Point Data;  // change this to make a kist of other objects

/** \ingroup ImageFindingL2
 * \brief A Kist is a linked list of Units which each point to a Data type.
 *
 * In this implementation the Data type is set to Point type, but this could
 * be changed for other applications.  Multiple Kists of the same points can be
 * made without copying data.  Memory is allocated in blocks in held in a reservoir
 * to improve efficiency when adding and removing data.
 */
template <class Data = Point>
struct Kist{
public:

	Kist(unsigned long my_blocksize = 10000){
		  top=NULL;
		  Number=0;
		  bottom = top;
		  current = top;
		  blocksize = my_blocksize;

		  reserve_top = NULL;
		  Nreserve = 0;

		  return;
	};
	Kist(Kist &a){
		Unit<Data> *current;

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
	};
	~Kist(){
		Unit<Data> *units;
		while(heads.size() > 0){
			units = heads.back();
			heads.pop_back();
			delete[] units;
		}
	};

	inline Data *getCurrent();

	void InsertAfterCurrent(Data * data);
	void InsertBeforeCurrent(Data * data);
	Data *TakeOutCurrent();
	void Empty();
	//void FreeAll();
	void Fill(Data * data,unsigned long N);
	void SwapCurrentWithBottom();
	void MoveCurrentToBottom();

	// movement
	bool JumpDown(int jump);
	bool Down();
	bool Up();
	bool MoveToTop();
	bool MoveToBottom();

	// status
	/// Number of elements in list.
	unsigned long Nunits(){return Number;}
	bool AtTop();
	bool AtBottom();
	void Print();

	void TranformPlanes();
	bool AreDataUnique();
	void SetInImage(Boo value);


	/// Returns a pointer to the current data
	Data *operator*(){return getCurrent();}
	/// Same as Up()
	bool operator++(){return Up();}
	/// Same as Down()
	bool operator--(){return Down();}
	/// Same as Up()
	bool operator++(int x){return Up();}
	/// Same as Down()
	bool operator--(int x){return Down();}

private:

	Unit<Data> * pop_from_reserve();
	void push_to_reserve(Unit<Data> *unit);

	unsigned long Number;
	unsigned long blocksize;

	std::vector<Unit<Data> *> heads;
	//Unit *units;
	Unit<Data> *reserve_top;
	unsigned long Nreserve;

	Unit<Data> *top;
	Unit<Data> *bottom;
	Unit<Data> *current;

};// Kist;

//typedef struct Kist * Kist *;


/** \ingroup ImageFindingL2
 * \brief Removes all elements from list without destroy the data.
 */
template <class Data> void Kist<Data>::Empty(){

	while(Number > 0) TakeOutCurrent();
	//Number = 0;
	//top = bottom = current = NULL;
	return;
}

/// Take a Unit out of the reservoir, expand reservoir if necessary
template <class Data> Unit<Data> * Kist<Data>::pop_from_reserve(){

	if(Nreserve <= 1){
		// expand reservoir of Units for later use
		  Unit<Data> *units = new Unit<Data>[blocksize];
		  for(unsigned long i=0;i<blocksize-1;++i) units[i].next = &units[i+1];
		  units[blocksize-1].next = NULL;

		  heads.push_back(units);  // keep list of first point in each block of memory
		  if(Nreserve == 1) reserve_top->next = units;
		  else reserve_top = units;
		  Nreserve += blocksize;
	}

	Unit<Data> *unit = reserve_top;
	reserve_top = unit->next;
	--Nreserve;

	return unit;
}

/// put a Unit into the reservoir
template <class Data> void Kist<Data>::push_to_reserve(Unit<Data> *unit){
	assert(reserve_top);
	assert(Nreserve > 0);

	unit->next = reserve_top;
	reserve_top = unit;
	++Nreserve;
}

template <class Data> bool Kist<Data>::AtTop(){
	if(current==top) return true;
	return false;
}

template <class Data> bool Kist<Data>::AtBottom(){
	if(current==bottom) return true;
	else return false;
}

/// Insert data into kist.  Leaves current unchanged.
template <class Data> void Kist<Data>::InsertAfterCurrent(Data *data){

    assert(data);
    assert(data->x);
    assert(data->gridsize >= 0);

    Unit<Data> *unit = pop_from_reserve();

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

/// Insert data into kist.  Leaves current unchanged.
template <class Data> void Kist<Data>::InsertBeforeCurrent(Data *data){

	assert(data);
    assert(data->x);
    assert(data->gridsize >= 0);

    //Unit *unit = (Unit *)malloc(sizeof(Unit));
    Unit<Data> *unit = pop_from_reserve();

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

/**
 * Swaps current data with bottom data leaving current one above former current.
 */
template <class Data> void Kist<Data>::SwapCurrentWithBottom(){

	Data *data;

	data=current->data;
	current->data=bottom->data;
	bottom->data=data;
	Up();
}

/**
 * Moves current to the bottom of the kist.  Current
 * is left at the bottom
 */
template <class Data> void Kist<Data>::MoveCurrentToBottom(){

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
template <class Data> Data *Kist<Data>::TakeOutCurrent(){

    Data *data;
    Unit<Data> *unit;

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

template <class Data> bool Kist<Data>::JumpDown(int jump){
	int i;
	bool ans;

	if(jump > 0) for(i=0;i<jump;++i) ans=Down();
	if(jump < 0) for(i=0;i<abs(jump);++i) ans=Up();
	return ans;
}
/// Move down the list. Returns false when at the bottom of the list.
template <class Data> bool Kist<Data>::Down(){

	if(Number == 0) return false;
	if(current==bottom) return false;
	current=current->next;

	return true;
}
/// Move up the list. Returns false when at the top of the list.
template <class Data> bool Kist<Data>::Up(){

	if(Number == 0) return false;
	if(current==top) return false;
	current=current->prev;

	return true;
}
/// Moves to the top.  Returns false if the kist is empty.
template <class Data> bool Kist<Data>::MoveToTop(){
	if(!Number)return false;
	current=top;
	return true;
}

/// Moves to the top.  Returns false if the kist is empty.
template <class Data> bool Kist<Data>::MoveToBottom(){
	if(!Number)return false;
	current=bottom;
	return true;
}

/** \brief Put an array of data into a kist.
 * */
template <class Data> void Kist<Data>::Fill(Data *data_array,unsigned long N){
  unsigned long i;

   for(i=1;i<N;++i){
    InsertAfterCurrent(&(data_array[i]));
    Kist<Data>::Down();
  }
}
/*********************************/
/*  data extraction routines */
/*********************************/
/// Return pointer to data in current element.
template <class Data> inline Data * Kist<Data>::getCurrent(){
	return Number ? current->data : NULL;
}
/*********************************
 * specific to image points
 * ******************************/
/**
 *\brief Transform all points in kist from image to source plane or vis versus
 */
template <class Data> void Kist<Data>::TranformPlanes(){

	if(Number == 0) return;
	MoveToTop();
	do{
		assert(current->data);
		assert(current->data->image);
		current->data = current->data->image;
	}while(Down());

	return;
}

/// Print positions and gridsizes of all points in Kist to standard out
template <class Data> void Kist<Data>::Print(){
	std::cout << Nunits() << std::endl;
	MoveToTop();
	do{
		std::cout << getCurrent()->x[0] <<  "  " << getCurrent()->x[1] << "  " << getCurrent()->gridsize << std::endl;
	}while(Down());
}
/**
 * \brief return true if all the points in the kist are unique.
 */
template <class Data> bool Kist<Data>::AreDataUnique(){

	if(Nunits() < 2) return true;
	Unit<Data> *unit = current;

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

/**
 * \brief Set the in_image variable in every member of kist to value.
 *
 * Does not return current to previous value.
 */
template <class Data> void Kist<Data>::SetInImage(Boo value){

	if(Nunits() == 0) return;

	MoveToTop();
	do{
		getCurrent()->in_image = value;
		getCurrent()->image->in_image = value;
	}while(Down());
}

/*
 * These functions are provided for backwards compatibility, but the kist member methods
 * should be used in the future.
 */

void InsertAfterCurrentKist(Kist<Point> * kist,Point * data);
void InsertBeforeCurrentKist(Kist<Point> * kist,Point * data);
bool AtTopKist(Kist<Point> * kist);
bool AtBottomKist(Kist<Point> * kist);
void JumpDownKist(Kist<Point> * kist,int jump);
bool MoveDownKist(Kist<Point> * kist);
bool MoveUpKist(Kist<Point> * kist);
bool MoveToTopKist(Kist<Point> * kist);
bool MoveToBottomKist(Kist<Point> * kist);
void FillKist(Kist<Point> * kist,Point * data,unsigned long N);
void PrintKist(Kist<Point> * kist);
void SwapCurrentWithBottomKist(Kist<Point> * kist);
Point *TakeOutCurrentKist(Kist<Point> * kist);
Point *GetCurrentKist(Kist<Point> * kist);
void EmptyKist(Kist<Point> * kist);
//void FreeAllKist(Kist<Point> * kist);
//void UnionKist(Kist<Point> * kist1,Kist<Point> * kist2);
//bool AreDataUniqueKist(Kist<Point> * kist);
bool IntersectionKist(Kist<Point> * kist1,Kist<Point> * kist2,Kist<Point> * intersection);
Point *getCurrentKist(Kist<Point> * kist);
void TranformPlanesKist(Kist<Point> * kist);

#endif
