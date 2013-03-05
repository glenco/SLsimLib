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

typedef struct Point Data;  // change this to make a kist of other objects

/** \brief Used as internal container in Kist */
typedef struct Unit{
	Data * data;
	struct Unit *next;
	struct Unit *prev;
} Unit;

/** \ingroup ImageFindingL2
 * \brief A Kist is a linked list of Units which each point to a Data type.
 *
 * In this implementation the Data type is set to Point type, but this could
 * be changed for other applications.  Multiple Kists of the same points can be
 * made without copying data.  Memory is allocated in blocks in held in a reservoir
 * to improve efficiency when adding and removing data.
 */
typedef struct Kist{
public:

	Kist(unsigned long blocksize = 10000);
	Kist(Kist &a);
	~Kist();

	Data *getCurrent();

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

	Unit * pop_from_reserve();
	void push_to_reserve(Unit *unit);

	unsigned long Number;
	unsigned long blocksize;

	std::vector<Unit *> heads;
	//Unit *units;
	Unit *reserve_top;
	unsigned long Nreserve;

	Unit *top;
	Unit *bottom;
	Unit *current;

} Kist;

typedef struct Kist *KistHndl;

/*
 * These functions are provided for backwards compatibility, but the kist member methods
 * should be used in the future.
 */

void InsertAfterCurrentKist(KistHndl kist,Data * data);
void InsertBeforeCurrentKist(KistHndl kist,Data * data);
bool AtTopKist(KistHndl kist);
bool AtBottomKist(KistHndl kist);
void JumpDownKist(KistHndl kist,int jump);
bool MoveDownKist(KistHndl kist);
bool MoveUpKist(KistHndl kist);
bool MoveToTopKist(KistHndl kist);
bool MoveToBottomKist(KistHndl kist);
void FillKist(KistHndl kist,Data * data,unsigned long N);
void PrintKist(KistHndl kist);
void SwapCurrentWithBottomKist(KistHndl kist);
Data *TakeOutCurrentKist(KistHndl kist);
Data *GetCurrentKist(KistHndl kist);
void EmptyKist(KistHndl kist);
//void FreeAllKist(KistHndl kist);
//void UnionKist(KistHndl kist1,KistHndl kist2);
//bool AreDataUniqueKist(KistHndl kist);
bool IntersectionKist(KistHndl kist1,KistHndl kist2,KistHndl intersection);
Data *getCurrentKist(KistHndl kist);
void TranformPlanesKist(KistHndl kist);

/*
 * Inline functions.
 */

/*********************************/
/*  data extraction routines */
/*********************************/
/// Return pointer to data in current element.
inline Data * Kist::getCurrent(){
	if(Number==0) return NULL;
	assert(current);
	assert(current->data);
	
	return current->data;
}
inline Data *getCurrentKist(KistHndl kist){
	assert(kist);
	
	return kist->getCurrent();
}

#endif
