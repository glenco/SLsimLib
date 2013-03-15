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
struct KistUnit{
	T * data;
	struct KistUnit<T> *next;
	struct KistUnit<T> *prev;
};// Unit;

//typedef struct Point Data;  // change this to make a kist of other objects

/** \ingroup ImageFindingL2
 * \brief A Kist is a class template for a linked list of any data type (default is Point).
 *
 * Multiple Kists of the same data can be made without copying data.  Memory is allocated in
 * blocks and held in a reservoir to improve efficiency when adding and removing data.  This
 * makes it faster then the STL list when data is repeatedly being added and subtracted from the
 * list.
 */
template <class Data = Point>
struct Kist{
public:

	Kist(unsigned long my_blocksize = 10000   /// Number of KistUnits allocated at one time.  Larger it is the less time spent allocating memory, but memory will be waisted for small Kists
			){
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
		KistUnit<Data> *current;

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
		KistUnit<Data> *units;
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


	/// Returns a pointer to the current data.  Same as getCurrent.
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

	KistUnit<Data> * pop_from_reserve();
	void push_to_reserve(KistUnit<Data> *unit);

	unsigned long Number;
	unsigned long blocksize;

	std::vector<KistUnit<Data> *> heads;
	//Unit *units;
	KistUnit<Data> *reserve_top;
	unsigned long Nreserve;

	KistUnit<Data> *top;
	KistUnit<Data> *bottom;
	KistUnit<Data> *current;

};// Kist;

//typedef struct Kist * Kist *;

/**
 * \brief Removes all elements from list without destroy the data.
 */
template <class Data> void Kist<Data>::Empty(){

	while(Number > 0) TakeOutCurrent();
	//Number = 0;
	//top = bottom = current = NULL;
	return;
}

/// Take a Unit out of the reservoir, expand reservoir if necessary
template <class Data> KistUnit<Data> * Kist<Data>::pop_from_reserve(){

	if(Nreserve <= 1){
		// expand reservoir of Units for later use
		  KistUnit<Data> *units = new KistUnit<Data>[blocksize];
		  for(unsigned long i=0;i<blocksize-1;++i) units[i].next = &units[i+1];
		  units[blocksize-1].next = NULL;

		  heads.push_back(units);  // keep list of first point in each block of memory
		  if(Nreserve == 1) reserve_top->next = units;
		  else reserve_top = units;
		  Nreserve += blocksize;
	}

	KistUnit<Data> *unit = reserve_top;
	reserve_top = unit->next;
	--Nreserve;

	return unit;
}

/// put a Unit into the reservoir
template <class Data> void Kist<Data>::push_to_reserve(KistUnit<Data> *unit){
	assert(reserve_top);
	assert(Nreserve > 0);

	unit->next = reserve_top;
	reserve_top = unit;
	++Nreserve;
}
/// True if current is at top of list
template <class Data> bool Kist<Data>::AtTop(){
	if(current==top) return true;
	return false;
}
/// True if current is at bottom of list
template <class Data> bool Kist<Data>::AtBottom(){
	if(current==bottom) return true;
	else return false;
}

/// Insert data into kist.  Leaves current unchanged.
template <class Data> void Kist<Data>::InsertAfterCurrent(Data *data){

    assert(data);
    assert(data->x);
    assert(data->gridsize >= 0);

    KistUnit<Data> *unit = pop_from_reserve();

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
    KistUnit<Data> *unit = pop_from_reserve();

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

/// Swaps current data with bottom data leaving current one above former current.
template <class Data> void Kist<Data>::SwapCurrentWithBottom(){

	Data *data;

	data=current->data;
	current->data=bottom->data;
	bottom->data=data;
	Up();
}


/// Moves data at current location to the bottom of the kist.  Current is left at the bottom
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
    KistUnit<Data> *unit;

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
/// Move down the list jump units from current position
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

/// Put an array of data into a kist.
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
 /** \brief Transform all points in kist from image to source plane or vis versus.
  * Data type must have a "image" attribute.
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
/// \brief Print data to stdout. Data type must have a "print()" public function.
template <class Data> void Kist<Data>::Print(){
	std::cout << Nunits() << std::endl;
	MoveToTop();
	do{
		getCurrent()->print();
	}while(Down());
}

 /// returns true if all the members have unique addresses.  Warning: This can be slow.
template <class Data> bool Kist<Data>::AreDataUnique(){

	if(Nunits() < 2) return true;
	KistUnit<Data> *unit = current;

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
 * Does not return current to previous value.  Data type must have a
 * "in_image" attribute.
 */
/*template <class Data> void Kist<Data>::SetInImage(Boo value){
	return;
}*/
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
inline void EmptyKist(Kist<Point> * kist){
	kist->Empty();
}
inline bool AtTopKist(Kist<Point> * kist){
	return kist->AtTop();
}
inline bool AtBottomKist(Kist<Point> * kist){
	return kist->AtBottom();
}
inline void InsertAfterCurrentKist(Kist<Point> * kist,Point *data){
	kist->InsertAfterCurrent(data);
}
inline void InsertBeforeCurrentKist(Kist<Point> * kist,Point *data){
	kist->InsertBeforeCurrent(data);
}
inline void SwapCurrentWithBottomKist(Kist<Point> * kist){
	kist->SwapCurrentWithBottom();
}
inline Point *TakeOutCurrentKist(Kist<Point> * kist){
	return kist->TakeOutCurrent();
}
inline void JumpDownKist(Kist<Point> * kist,int jump){
	kist->JumpDown(jump);
}
inline bool MoveDownKist(Kist<Point> * kist){
	return kist->Down();
}
inline bool MoveUpKist(Kist<Point> * kist){
	return kist->Up();
}
inline bool MoveToTopKist(Kist<Point> * kist){
	return kist->MoveToTop();
}
inline bool MoveToBottomKist(Kist<Point> * kist){
	return kist->MoveToBottom();
}
inline void FillKist(Kist<Point> * kist,Point *data_array,unsigned long N){
	kist->Fill(data_array,N);
}
inline Point *getCurrentKist(Kist<Point> * kist){
	return kist->getCurrent();
}
inline void TranformPlanesKist(Kist<Point> * kist){
	kist->TranformPlanes();
}

#endif
