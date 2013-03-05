/*
 * KistDriver.h
 *
 *  Created on: Nov 16, 2010
 *      Author: bmetcalf
 */

#ifndef KIST_DRIVER_H_
#define KIST_DRIVER_H_

#include "Tree.h"

void SetInImage(KistHndl kist,Boo value);
void DirtyFoF(ImageInfo *imageinfo ,int *Nimages ,double linkinglength ,int MaxNimages );
void _DirtyFoF(KistHndl neighbors,KistHndl wholekist,double linkinglength);
void DirtyDivider(ImageInfo *imageinfo,int *Nimages ,int MaxNimages ,int Ngroup);
void _DirtyDivider(KistHndl neighbors,KistHndl wholekist,int Ngroup);

#endif
