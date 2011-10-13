/*
 * Kist.h
 *
 *  Created on: Feb 2, 2010
 *      Author: R.B. Metcalf
 */

#include <point.h>  // only reason is to define Point

#ifndef kisttypes_declare
	#define kisttypes_declare

	typedef struct Point Data;  // change this to make a kist of other objects

	/** \brief Used as internal container in Kist */
	typedef struct Unit{
		Data * data;
		struct Unit *next;
		struct Unit *prev;
	} Unit;

	/** \brief
	 * A Kist is a linked list of Units which each point to a Data type.
	 *
	 * In this implementation the Data type is set to Point type, but this could
	 * be changed for other applications.  Multiple Kists of the same points can be
	 * made without copying data.
	 */
	typedef struct Kist{
		Unit *top;
		Unit *bottom;
		Unit *current;
		unsigned long Nunits;
	} Kist;

	typedef struct Kist *KistHndl;

#endif

KistHndl NewKist(void);
void InsertAfterCurrentKist(KistHndl kist,Data * data);
void InsertBeforeCurrentKist(KistHndl kist,Data * data);
Boolean AtTopKist(KistHndl kist);
Boolean AtBottomKist(KistHndl kist);
void JumpDownKist(KistHndl kist,int jump);
Boolean MoveDownKist(KistHndl kist);
Boolean MoveUpKist(KistHndl kist);
void MoveToTopKist(KistHndl kist);
void MoveToBottomKist(KistHndl kist);
void FillKist(KistHndl kist,Data * data,unsigned long N);
void SwapDataInKist(KistHndl kist,Unit *p1,Unit *p2);
void PrintKist(KistHndl kist);
void MoveCurrentToBottomKist(KistHndl kist);
Data *TakeOutCurrentKist(KistHndl kist);
Data *GetCurrentKist(KistHndl kist);
void EmptyKist(KistHndl kist);
void freeKist(KistHndl kist);
void FreeAllKist(KistHndl kist);
void UnionKist(KistHndl kist1,KistHndl kist2);
Boolean AreDataUniqueKist(KistHndl kist);
Boolean IntersectionKist(KistHndl kist1,KistHndl kist2,KistHndl intersection);
Data *getCurrentKist(KistHndl kist);
void TranformPlanesKist(KistHndl kist);
