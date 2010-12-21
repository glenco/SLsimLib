/*
 * Code Name:     kist
 * Programmer:    R Ben Metcalf
 * Last Revised:  Nov, 2005                                   
 * Discription:  replaces linked list programs
 * Comments:                           
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "Kist.h"

/***** Structs *****/

/***********************************************************
   routines for linked list of Datas
************************************************************/

KistHndl NewKist(void){
  KistHndl kist;

  kist=(Kist *) malloc(sizeof(Kist));
  if (!kist){
    fprintf(stderr,"allocation failure in NewKist()\n");
    exit(1);
  }
  kist->top=NULL;
  kist->Nunits=0;
  kist->bottom = kist->top;
  kist->current = kist->top;

  return kist;
}

Boolean AtTopKist(KistHndl kist){
	assert(kist);

	if(kist->current==kist->top) return True;
	else return False;
}
Boolean AtBottomKist(KistHndl kist){
	assert(kist);

	if(kist->current==kist->bottom) return True;
	else return False;
}

void InsertAfterCurrentKist(KistHndl kist,Data *data){
	/* leaves current unchanged */
	assert(kist);
	assert(data);

    Unit *unit = (Unit *)malloc(sizeof(Unit));

    unit->data = data;

    if(kist->Nunits > 0){
      assert(kist->current);

      unit->prev = kist->current;
      unit->next = kist->current->next;

      if(kist->current == kist->bottom) kist->bottom = unit;
      else kist->current->next->prev=unit;
      kist->current->next=unit;
    }else{  // empty kist case
      kist->current=unit;
      kist->top=unit;
      kist->bottom=unit;

      unit->prev = NULL;
      unit->next = NULL;
    }

    kist->Nunits++;
    return;
}

void InsertBeforeCurrentKist(KistHndl kist,Data *data){
	// leaves current unchanged
	assert(kist);
	assert(data);

	Unit *unit;

	unit=(Unit *)malloc(sizeof(Unit));
	unit->data=data;

    if(kist->Nunits > 0){
      assert(kist->current);

      unit->prev=kist->current->prev;
      unit->next=kist->current;
      if(kist->current == kist->top) kist->top=unit;
      else kist->current->prev->next=unit;
      kist->current->prev=unit;
    }else{  /* empty kist case */
      kist->current=unit;
      kist->top=unit;
      kist->bottom=unit;

      unit->prev = NULL;
      unit->next = NULL;
    }

    kist->Nunits++;
    return;
}

void MoveCurrentToBottomKist(KistHndl kist){
	assert(kist);

	Data *data;
	// leaves current one above former current

	data=kist->current->data;
	kist->current->data=kist->bottom->data;
	kist->bottom->data=data;
	MoveUpKist(kist);
}

Data *TakeOutCurrentKist(KistHndl kist){
	assert(kist);

    Data *data;
    Unit *unit;

    /* takes out current data and set current to data previous */
    /* Except at top where current is set to new top */
    /* returns dataer to removed data */

    if(kist == NULL || kist->Nunits <= 0) return NULL;

    assert(kist->current);
    assert(kist->top);
    assert(kist->bottom);

    data = kist->current->data;
    unit = kist->current;

    if(kist->top == kist->bottom){  /* last data */
      kist->current=NULL;
      kist->top=NULL;
      kist->bottom=NULL;
    }else if(kist->current==kist->top){
      kist->top=kist->top->next;
      kist->current=kist->top;
      kist->top->prev=NULL;
    } else if(kist->current==kist->bottom){
      kist->bottom=kist->bottom->prev;
      kist->current=kist->bottom;
      kist->bottom->next=NULL;
    }else{
      kist->current->prev->next=kist->current->next;
      kist->current->next->prev=kist->current->prev;
      kist->current=kist->current->prev;
    }
    
    free(unit);
    kist->Nunits--;

    return data;
}

void EmptyKist(KistHndl kist){
	// reduce kist to no elements.
	// Note: Does not destroy data.
	assert(kist);
	while(kist->Nunits > 0) TakeOutCurrentKist(kist);
	return;
}

void freeKist(KistHndl kist){
	// deallocate memory for kist
	// Note: Does not destroy data.  A handle must still point
	//    to the data
	EmptyKist(kist);
	free(kist);
	return;
}

void FreeAllKist(KistHndl kist){
	// deallocate memory for kist and all points in list
	// Note: Does not destroy data.  A handle must still point
	//    to the data
	assert(kist);
	Data *data[kist->Nunits];
	unsigned long Nheads=0,Ndata=0,i,Nunits;

	Nunits = kist->Nunits;
	while(kist->Nunits > 0){
		data[i] = TakeOutCurrentKist(kist);
		Ndata += data[i]->head;
		if(data[i]->head) ++Nheads;
	}

	for(i = 0; i < Nheads ; ++i){
		free(data[i]->x);
		free(data[i]);
	}

	if(Nunits != Ndata){
		ERROR_MESSAGE();
		printf("FreeAllKist freeed no all of data in kist\n");
		exit(0);
	}
	return;
}

/*
 *  Moving through kist
 */

void JumpDownKist(KistHndl kist,int jump){
	assert(kist);
	int i;

	if(jump > 0) for(i=0;i<jump;++i) MoveDownKist(kist);
	if(jump < 0) for(i=0;i<abs(jump);++i) MoveUpKist(kist);
	return ;
}

Boolean MoveDownKist(KistHndl kist){
	assert(kist);

	if(kist->Nunits == 0) return False;
	if(kist->current==kist->bottom) return False;
	kist->current=kist->current->next;

	return True;
}

Boolean MoveUpKist(KistHndl kist){
	assert(kist);

	if(kist->Nunits == 0) return False;
	if(kist->current==kist->top) return False;
	kist->current=kist->current->prev;

	return True;
}


void MoveToTopKist(KistHndl kist){
	assert(kist);
	kist->current=kist->top;
}
void MoveToBottomKist(KistHndl kist){
	assert(kist);
	kist->current=kist->bottom;
}

void FillKist(KistHndl kist,Data *data_array,unsigned long N){
  unsigned long i;
  /* add N data to to end of kist */
  /* id numbers are given in order from idmin */
  /* this is used to initialize kist */

   for(i=1;i<N;++i){
    InsertAfterCurrentKist(kist,&(data_array[i]));
    MoveDownKist(kist);
  }
}

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

/*********************************/
/*  data extraction routines */
/*********************************/
Data *getCurrentKist(KistHndl kist){
	assert(kist);
	assert(kist->current);
	assert(kist->current->data);

	return kist->current->data;
}

/*********************************
 * specific to image points
 * ******************************/

void TranformPlanesKist(KistHndl kist){
	assert(kist);

	if(kist->Nunits == 0) return;
	MoveToTopKist(kist);
	do{
		assert(kist->current->data);
		assert(kist->current->data->image);
		kist->current->data = kist->current->data->image;
	}while(MoveDownKist(kist));

	return;
}

Boolean AreDataUniqueKist(KistHndl kist){
	assert(kist);
	if(kist->Nunits < 2) return True;

	unsigned long i,j;
	KistHndl tmpkist = NewKist();

	// clone kist
	//  the two lists use the same units as well
	//  as data.
	tmpkist->Nunits = kist->Nunits;
	tmpkist->top = kist->top;
	tmpkist->bottom = kist->bottom;

	MoveToTopKist(kist);
	for(i=0;i<kist->Nunits;++i,MoveDownKist(kist)){
		MoveToTopKist(tmpkist);
		for(j=0;j<i;++j,MoveDownKist(tmpkist)){
			if( getCurrentKist(kist) == getCurrentKist(tmpkist) ){
				free(tmpkist);
				return False;
			}
		}
	}

	free(tmpkist);
	return True;
}
