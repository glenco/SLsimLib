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

/// Internal container used in Kist container class
template <class T>
struct KistUnit{
	T * data;
	struct KistUnit<T> *next;
	struct KistUnit<T> *prev;
};


//typedef struct Point Data;  // change this to make a kist of other objects

/** \ingroup ImageFindingL2
 * \brief A Kist is a class template for a linked list of any data type (default is Point).
 *
 * Multiple Kists of the same data can be made without copying data.  Memory is allocated in
 * blocks and held in a reservoir to improve efficiency when adding and removing data.  This
 * makes it faster then the STL list when data is repeatedly being added and subtracted from the
 * list.
 *
 <pre>
 example:
 
 Kist<Point> kist;
 Point *p_point = NewPointArray(100);
 int i=0;
 
 for(i = 0; i< 100; ++i){
 p_point[i].x[0] = p_point[i].x[1] = ran2(&seed);
 kist.InsertAfterCurrent(&p_point[i]);
 --kist;
 }
 cout << "Number of points " << kist.Nunits() << endl;
 
 for(i=0,kist.MoveToTop() ; !(kist.OffBottom()) ; --kist ){
 cout << i << " x = " << (*kist)->x[0] << "  " << (*kist)->x[1] << endl;
 ++i;
 }
 
 There is also an iterator Kist<T>::iterator class which works like this:
 
 Kist<Point>::iterator kist;
 
 ....
 
 for(Kist<Point>::iterator it = kist.begin(); !(it.atend()) ;++it){
   cout << i << " x = " << (*it).x[0] << "  " << (*it).x[1] << endl;
 }
 
 or
 
 for(auto it = kist.begin() ; it != kist.end() ; ++it){
  cout  " x = " << (*it).x[0] << "  " << (*it).x[1] << endl;
 }
 
 
 
<\pre>
 */
template <class Data>
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
		KistUnit<Data> *tmp_current;

		top=NULL;
		Number=0;
		bottom = top;
		current = top;

		tmp_current = a.current;
		a.MoveToTop();
		do{
			InsertAfterCurrent(a.getCurrent());
			Down();
		}while(a.Down());
		assert(a.Nunits() == Nunits());
		a.current = tmp_current;
	};
  
	~Kist(){
		KistUnit<Data> *units;
		while(heads.size() > 0){
			units = heads.back();
			heads.pop_back();
			delete[] units;
		}
	};
  
  Kist & operator=(const Kist &a){
    if(this == &a) return *this;
    copy(a);
    return *this;
  }

	inline Data *getCurrent();

	void InsertAfterCurrent(Data * data);
	void InsertBeforeCurrent(Data * data);
	Data *TakeOutCurrent();
	void Empty();
  void Clear();
    
	//void FreeAll();
	void Fill(Data * data,unsigned long N);
	void SwapCurrentWithBottom();
	void MoveCurrentToBottom();
	void MoveCurrentToTop();
	void copy(Kist<Data> *kist);
  void copy(const Kist<Data> &kist);
  void copy(std::vector<Data *> &vector);
  void add(Kist<Data> *kist);
  void add(const Kist<Data> &kist);
  void test_iterator();

	// movement
	bool JumpDown(int jump);
	bool Down();
	bool Up();
	bool MoveToTop();
	bool MoveToBottom();
    
    /// Custom iterator class for Kist container
    class iterator{
      
    public:
      
      iterator(){
        unit = NULL;
      }
 
      iterator(const iterator &my_it ){
        unit = my_it.unit;
      }

        /// Returns a pointer to the current data.  Same as getCurrent.
        Data & operator*(){return *(unit->data);}
      
        bool atend(){return (unit==NULL);}
      
        iterator& operator++(){
          
          if(unit == NULL ) return *this;
          unit = unit->prev;
          return *this;
       }

        iterator operator++(int){
          
          if(unit == NULL) return *this;
             iterator tmp = *this;
            unit = unit->prev;
             return tmp;
        }
        
        /// Same as Down()
        iterator& operator--(){
          
          
          if(unit == NULL ) return *this;
          unit = unit->next;
          return *this;
      }
        
        /// Same as Down()
        iterator operator--(int){
          
          if(unit == NULL) return *this;
             iterator tmp = *this;
             unit = unit->next;
             return tmp;
        }
        
        iterator& operator=(KistUnit<Data> *my_unit){
             unit = my_unit;
            return *this;
        }
        
        iterator& operator=(const iterator &my_it){
            if(this == &my_it) return *this;
          unit = my_it.unit;
            return *this;
        }
      
      bool operator==(const iterator my_it){
        return (unit == my_it.unit);
       }
      
      bool operator!=(const iterator my_it){
        return (unit != my_it.unit);
      }

    private:
      struct KistUnit<Data> *unit;
      //struct KistUnit<Data> offbot2;

    };

//  void SetCurrentIt(iterator it){current = it.getUnit();}
//    void SetCurrentIt(iterator it){current = it.unit;}
    
  /// returns iterator pointing to entry that current it pointing to
  Kist<Data>::iterator CurrentIt(){
      Kist<Data>::iterator it ;
      it = current;
      return it;
  }
  /// returns iterator pointing to first entry with data
  Kist<Data>::iterator TopIt(){
      Kist<Data>::iterator it ;
      it = top;
      return it;
  }
  /// returns iterator pointing to last entry with data
  Kist<Data>::iterator BottomIt(){
    Kist<Data>::iterator it;
    it = bottom;
    return it;
  }
  
  Kist<Data>::iterator begin(){
    Kist<Data>::iterator it;
    it = bottom;
    return it;
  }
  
  Kist<Data>::iterator end(){
    Kist<Data>::iterator it;
    it = NULL;
    return it;
  }

	// status
	/// Number of elements in list.
	unsigned long Nunits(){return Number;}
	bool AtTop();
	bool AtBottom();
	void Print();

	void TranformPlanes();
	bool AreDataUnique();
	void SetInImage(Boo value);

  /// Test if Down() (or --) was last called from last element in list.  Used to stop a for loop. 
  bool OffBottom(){
    return (current == &offbot);
  }
  
	/// Returns a pointer to the current data.  Same as getCurrent.
	Data *operator*(){return getCurrent();}
	/// Same as Up()
	bool operator++(){return Up();}
	/// Same as Down()
	bool operator--(){return Down();}
	/// Same as Up()
	bool operator++(int){return Up();}
	/// Same as Down()
	bool operator--(int){return Down();}
  
  /// make a copy of the data in a vector of pointers form
  std::vector<Data *> copytovector(){
    std::vector<Data *> vec(Number);

    Kist<Data>::iterator it=TopIt();
    for(size_t i=0;!(it.atend()) ;--it,++i){
      vec[i] = &(*it);
    }
    
    return vec;
  }
  

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
  
  KistUnit<Data> offbot;

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

/**
 * \brief Removes all elements from list and destroy the data.
 */
template <class Data> void Kist<Data>::Clear(){
  std::list<Data *> list;
  Data* point;
  
	while(Number > 0){
    point = TakeOutCurrent();
    if(point->head) list.push_back(point);
  }
  while(list.size() > 0){
    point = list.back();
    list.pop_back();
    delete [] point;
  }

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

/// Insert data into kist.  Leaves current unchanged. If off end it acts as if it where at the bottom.
template <class Data> void Kist<Data>::InsertAfterCurrent(Data *data){

    assert(data);
    assert(data->x);
    assert(data->gridsize >= 0);
    //assert(current != &offbot);
  
    KistUnit<Data> *unit = pop_from_reserve();

    unit->data = data;

    if(Number > 0){
      assert(current);

      if(current == &offbot) current = bottom;
      unit->prev = current;
      unit->next = current->next;

      if(current == bottom) bottom = unit;
      else current->next->prev = unit;
      current->next = unit;

      if(unit == bottom) unit->next = NULL;

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

/// Insert data into kist.  Leaves current unchanged. If offend acts as if it where at bottom.
template <class Data> void Kist<Data>::InsertBeforeCurrent(Data *data){

	assert(data);
  assert(data->x);
  assert(data->gridsize >= 0);

  //Unit *unit = (Unit *)malloc(sizeof(Unit));
  KistUnit<Data> *unit = pop_from_reserve();

	unit->data=data;

    if(Number > 0){
      assert(current);

      if(current == &offbot) current=bottom;
      
      unit->prev=current->prev;
      unit->next=current;
      if(current == top) top=unit;
      else current->prev->next=unit;
      current->prev=unit;
      
      if(unit == top) unit->prev = NULL;

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

  if(Number == 0 || current==&offbot) return;
	Data *data;

	data=current->data;
	current->data=bottom->data;
	bottom->data=data;
	Up();
}


/// Moves data at current location to the bottom of the kist.  Current is left at the bottom
template <class Data> void Kist<Data>::MoveCurrentToBottom(){
  
  if(Number == 0 || current==&offbot) return;
	Data *data = TakeOutCurrent();
	MoveToBottom();
	InsertAfterCurrent(data);
	Down();
}
/// Moves data at current location to the top of the kist.  Current is left at the top
template <class Data> void Kist<Data>::MoveCurrentToTop(){
  
  if(Number == 0 || current==&offbot) return;
	Data *data = TakeOutCurrent();
	MoveToTop();
	InsertBeforeCurrent(data);
	Up();
}

/**
 * \brief Takes out current data and set current to data previous
* except at top where current is set to new top.
* Returns pointer to removed data.  If offend acts as if at bottom.
*/
template <class Data> Data *Kist<Data>::TakeOutCurrent(){

    Data *data;
    KistUnit<Data> *unit;

    if( Number <= 0 ) return NULL;
  
    if( current == &offbot ) current = bottom;

    assert(current);
    assert(top);
    assert(bottom);

    data = current->data;
    unit = current;

    if(top == bottom){  /* last data */
      current=&offbot;
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
	bool ans = true;

	if(jump > 0) for(i=0;i<jump;++i) ans=Down();
	if(jump < 0) for(i=0;i<abs(jump);++i) ans=Up();
	return ans;
}
/// Move down the list. Returns false when at the bottom of the list.
template <class Data> bool Kist<Data>::Down(){

	if(Number == 0) return false;
	if(current==bottom || current == &offbot){
    current = &offbot;
    return false;
  }
	current=current->next;

	return true;
}
/// Move up the list. Returns false when at the top of the list. If offend act as if at bottom.
template <class Data> bool Kist<Data>::Up(){

  if(current == &offbot) current=bottom;
	if(Number == 0) return false;
	if(current==top) return false;
	current=current->prev;

	return true;
}
/// Moves to the top.  Returns false if the kist is empty.
template <class Data> bool Kist<Data>::MoveToTop(){
	if(!Number){current = &offbot; return false;}
	current = top;
	return true;
}

/// Moves to the top.  Returns false if the kist is empty.
template <class Data> bool Kist<Data>::MoveToBottom(){
	if(!Number){current=&offbot; return false;}
	current = bottom;
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
/// Return pointer to data in current element.  If offend act as if it is at bottem.
template <class Data> inline Data * Kist<Data>::getCurrent(){
  if(Number == 0 || current == NULL) return NULL;
  if(current == &offbot) current = bottom;
	return current->data;
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
/**
 * \brief copy contents of kist into this kist.  Destroys former content and leaves current of kist in
 * new place.
 */
template <class Data>
void Kist<Data>::copy(Kist<Data> *kist){
  Empty();
  if(kist->Nunits() == 0) return;
  kist->MoveToTop();
  do{
    InsertAfterCurrent(kist->getCurrent());
    Down();
  }while(kist->Down());
}
/**
 * \brief copy contents of kist into this kist.  Destroys former content and leaves current of kist in
 * new place.
 */
template <class Data>
void Kist<Data>::copy(const Kist<Data> &kist){
  Empty();
  if(kist.Nunits() == 0) return;
  Kist<Data>::iterator it = kist.getBottomIt();
  for(;!(it.atend()) ;++it){
    InsertAfterCurrent(*it);
    Down();
  };
}
template <class Data>
void Kist<Data>::copy(std::vector<Data *> &vector){
	Empty();
	if(vector.size() == 0) return;
	for(size_t i=0;i<vector.size();++i){
		InsertAfterCurrent(vector[i]);
		Down();
	}
}
/**
 * \brief copy contents of kist into this kist without destroying what is already there.
 */
template <class Data>
void Kist<Data>::add(Kist<Data> *kist){
  if(kist->Nunits() == 0) return;
  kist->MoveToTop();
  do{
    InsertAfterCurrent(kist->getCurrent());
    Down();
  }while(kist->Down());
}
/**
 * \brief copy contents of kist into this kist without destroying what is already there.
 */
template <class Data>
void Kist<Data>::add(const Kist<Data> &kist){
  if(kist.Nunits() == 0) return;
  Kist<Data>::iterator it = kist.getBottomIt();
  for(;!(it.atend()) ;++it){
    InsertAfterCurrent(*it);
    Down();
  };
}


template <class Data> void Kist<Data>::test_iterator(){
    
    MoveToBottom();
    int i;
    Kist<Data>::iterator it = this->getBottomIt();
    for(i=0;!(it.atend()) ;++it,++i){
        assert(current->data == (*it));
        std::cout << "i = " << i << std::endl;
        Up();
    }
    assert(i == Nunits()); 
}

template <class Data> void test_iterator(Kist<Data> &kist){
  
  kist.MoveToBottom();
  int i=0;
  for(auto p = kist.begin() ; p != kist.end(); ++p ){
    assert( &(*p)== kist.getCurrent() );
    (*p).print();
    std::cout << "i = " << i << std::endl;
    kist.Up();
    ++i;
  }
  assert(i == kist.Nunits());
}

template <class Data> void test_iterator2(Kist<Data> &kist){
  int i=0;
  auto p2 = kist.begin();
  for(auto p : kist ){
    
    assert(p.id == (*p2).id);
    
    p.print();
    std::cout << "i = " << i << std::endl;
    ++p2;
    ++i;
  }
  assert(i == kist.Nunits());
}
/*
 * These functions are provided for backwards compatibility, but the kist member methods
 * should be used in the future.
 *
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
*/

#endif
