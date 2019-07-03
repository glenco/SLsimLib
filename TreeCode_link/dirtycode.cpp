/*
 * dirtycode.c
 *
 *  Created on: Jan 31, 2012
 *      Author: bmetcalf
 */

#include "slsimlib.h"

/** 
 *
 *  \brief Divides kist of points into friends-of-friends groups with.  Does not need
 *  a tree to be built, but takes N^2 time so it is not good for a large number of points.
 *
 *  On entry all points must be in imageinfo[0].imagekist.  On exit they are divided into
 *  imageinfo[0...*Nimages=1].imagekist
 *
 *  If linkinglength > 0.0 a fixed linkinglength is used.  If <= 0.0 the SQRT(2)*(gridsize + gridsize)
 *  is used for a linking length.
 */

void DirtyFoF(
    std::vector<ImageInfo> &imageinfo  /// Contains the kists of points in each image
		,int *Nimages        /// Number of images
		,PosType linkinglength /// linking length, If it is <= 0 the gridsize's for the points is used.
		){

	Kist<Point> * wholekist = new Kist<Point>;
	unsigned long i;

	*Nimages = 0;

	while(imageinfo[0].imagekist->Nunits() > 0) wholekist->InsertAfterCurrent(imageinfo[0].imagekist->TakeOutCurrent());

	i=0;
	while(wholekist->Nunits() > 0){
    if(i > imageinfo.size()-1) imageinfo.resize(i+2);
		imageinfo[i].imagekist->Empty();
		imageinfo[i].imagekist->InsertAfterCurrent(wholekist->TakeOutCurrent());
		if(wholekist->Nunits() > 0) _DirtyFoF(imageinfo[i].imagekist,wholekist,linkinglength);
		++(*Nimages);
		++i;
	}

	delete wholekist;
	return;
}
void _DirtyFoF(Kist<Point> * neighbors,Kist<Point> * wholekist,PosType linkinglength){
	PosType ll2 = linkinglength*linkinglength;
	bool check = false;

	assert(neighbors->Nunits() == 1);
	do{
		wholekist->MoveToTop();
		do{
			if(check){
				wholekist->Up();
				check=false;
			}
			if(linkinglength <= 0) ll2 = 2.01*pow(neighbors->getCurrent()->gridsize + wholekist->getCurrent()->gridsize,2);
			if( ll2 > pow(neighbors->getCurrent()->x[0] - wholekist->getCurrent()->x[0],2)
					+ pow(neighbors->getCurrent()->x[1] - wholekist->getCurrent()->x[1],2)   ){
				if(wholekist->AtTop()) check=true;
				neighbors->InsertAfterCurrent(wholekist->TakeOutCurrent());
			}
		}while(wholekist->Down());
	}while(neighbors->Down() && wholekist->Nunits() > 0);

	return ;
}

/** 
 *
 *  \brief Divides kist of points into groups.
 *
 *  Does not need a tree to be built, but takes N^2 time so it is not good for a large number of points.
 *
 *  On entry all points must be in imageinfo[0].imagekist.  On exit they are divided into
 *  imageinfo[0...*Nimages=1].imagekist with each containing approximately <= Ngroup points.  The left most point
 *  is used as a seed to start the first group and the closest Ngroup points are found.  If there are not Ngroup
 *  points left, they are all put into one group.
 *
 */

void DirtyDivider(
		ImageInfo *imageinfo  /// Contains the kists of points in each image
		,int *Nimages        /// Number of images
		,int MaxNimages       /// Maximum size of imageinfo array
		,int Ngroup
		){

	if(imageinfo->imagekist->Nunits() <= Ngroup){
		if(imageinfo->imagekist->Nunits() == 0) *Nimages = 0;
		else *Nimages = 1;
		return;
	}

	Kist<Point> * wholekist = new Kist<Point>;
	unsigned long i;
	PosType xmin;

	*Nimages = 0;

	xmin = imageinfo->imagekist->getCurrent()->x[0];
	while(imageinfo->imagekist->Nunits() > 0){
		if(xmin < imageinfo->imagekist->getCurrent()->x[0] ){  // keeps most left point first
			xmin = imageinfo->imagekist->getCurrent()->x[0];
			wholekist->MoveToTop();
		}
		wholekist->InsertBeforeCurrent(imageinfo->imagekist->TakeOutCurrent());
	}

	i=0;
	while(wholekist->Nunits() > 0 && *Nimages < MaxNimages){
		imageinfo[i].imagekist->Empty();
		imageinfo[i].imagekist->InsertAfterCurrent(wholekist->TakeOutCurrent());
		if(wholekist->Nunits() > 0) _DirtyDivider(imageinfo[i].imagekist,wholekist,Ngroup);
		++(*Nimages);
		++i;
	}

	assert(*Nimages < MaxNimages);

	delete wholekist;
	return;
}
void _DirtyDivider(Kist<Point> * neighbors,Kist<Point> * wholekist,int Ngroup){
	bool check = false;
	PosType xinit[2],rr[Ngroup+1],ll;
	long i,j;

	if(wholekist->Nunits() <= Ngroup){
		while(wholekist->Nunits() > 0) neighbors->InsertBeforeCurrent(wholekist->TakeOutCurrent());
		return;
	}

	xinit[0] = neighbors->getCurrent()->x[0];
	xinit[1] = neighbors->getCurrent()->x[1];

	for(i=0;i<Ngroup+1;++i) rr[i] = 1.0e200;

	wholekist->MoveToTop();
	do{
		if(check){
			wholekist->Up();
			check=false;
		}

		ll = pow(xinit[0]-wholekist->getCurrent()->x[0],2)
				+ pow(xinit[1]-wholekist->getCurrent()->x[1],2);

		if(ll < rr[Ngroup]){
			neighbors->MoveToTop();
			for(i=0;i<Ngroup+1;++i,neighbors->Down()){
				if(ll < rr[i]){
					neighbors->InsertBeforeCurrent(wholekist->TakeOutCurrent());
					for(j=Ngroup;j>i;--j) rr[j] = rr[j-1];
					rr[i] = ll;
					break;
				}
			}
			if(wholekist->AtTop()) check=true;
		}
	}while(wholekist->Down());

	// put nearby points at the beginning of wholekist
	neighbors->MoveToBottom();
	wholekist->MoveToTop();
	while(neighbors->Nunits() > Ngroup){
		wholekist->InsertBeforeCurrent(neighbors->TakeOutCurrent());
		wholekist->MoveToTop();
	}

	return ;
}
