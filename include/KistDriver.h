/*
 * KistDriver.h
 *
 *  Created on: Nov 16, 2010
 *      Author: bmetcalf
 */

#ifndef KIST_DRIVER_H_
#define KIST_DRIVER_H_

#include "Tree.h"

void SetInImage(Kist<Point> * kist,Boo value);
void DirtyFoF(std::vector<ImageInfo> &imageinfo ,int *Nimages ,double linkinglength );
void _DirtyFoF(Kist<Point> * neighbors,Kist<Point> * wholekist,double linkinglength);
void DirtyDivider(ImageInfo *imageinfo,int *Nimages ,int MaxNimages ,int Ngroup);
void _DirtyDivider(Kist<Point> * neighbors,Kist<Point> * wholekist,int Ngroup);

#endif
