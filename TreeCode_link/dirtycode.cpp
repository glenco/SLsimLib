/*
 * dirtycode.c
 *
 *  Created on: Jan 31, 2012
 *      Author: bmetcalf
 */

#include <slsimlib.h>

/** \ingroup ImageFindingL2
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
		ImageInfo *imageinfo  /// Contains the kists of points in each image
		,int *Nimages        /// Number of images
		,double linkinglength /// linking length, If it is <= 0 the gridsize's for the points is used.
		,int MaxNimages       /// Maximum size of imageinfo array
		){

	KistHndl wholekist = new Kist;
	unsigned long i;

	*Nimages = 0;

	while(imageinfo->imagekist->Nunits() > 0) InsertAfterCurrentKist(wholekist,TakeOutCurrentKist(imageinfo->imagekist));

	i=0;
	while(wholekist->Nunits() > 0 && *Nimages < MaxNimages){
		EmptyKist(imageinfo[i].imagekist);
		InsertAfterCurrentKist(imageinfo[i].imagekist,TakeOutCurrentKist(wholekist));
		if(wholekist->Nunits() > 0) _DirtyFoF(imageinfo[i].imagekist,wholekist,linkinglength);
		++(*Nimages);
		++i;
	}

	assert(*Nimages < MaxNimages);

	delete wholekist;
	return;
}
void _DirtyFoF(KistHndl neighbors,KistHndl wholekist,double linkinglength){
	double ll2 = linkinglength*linkinglength;
	bool check = false;

	assert(neighbors->Nunits() == 1);
	do{
		MoveToTopKist(wholekist);
		do{
			if(check){
				MoveUpKist(wholekist);
				check=false;
			}
			if(linkinglength <= 0) ll2 = 2.01*pow(getCurrentKist(neighbors)->gridsize + getCurrentKist(wholekist)->gridsize,2);
			if( ll2 > pow(getCurrentKist(neighbors)->x[0] - getCurrentKist(wholekist)->x[0],2)
					+ pow(getCurrentKist(neighbors)->x[1] - getCurrentKist(wholekist)->x[1],2)   ){
				if(AtTopKist(wholekist)) check=true;
				InsertAfterCurrentKist(neighbors,TakeOutCurrentKist(wholekist));
			}
		}while(MoveDownKist(wholekist));
	}while(MoveDownKist(neighbors) && wholekist->Nunits() > 0);

	return ;
}

/** \ingroup ImageFindingL2
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

	KistHndl wholekist = new Kist;
	unsigned long i;
	double xmin;

	*Nimages = 0;

	xmin = getCurrentKist(imageinfo->imagekist)->x[0];
	while(imageinfo->imagekist->Nunits() > 0){
		if(xmin < getCurrentKist(imageinfo->imagekist)->x[0] ){  // keeps most left point first
			xmin = getCurrentKist(imageinfo->imagekist)->x[0];
			MoveToTopKist(wholekist);
		}
		InsertBeforeCurrentKist(wholekist,TakeOutCurrentKist(imageinfo->imagekist));
	}

	i=0;
	while(wholekist->Nunits() > 0 && *Nimages < MaxNimages){
		EmptyKist(imageinfo[i].imagekist);
		InsertAfterCurrentKist(imageinfo[i].imagekist,TakeOutCurrentKist(wholekist));
		if(wholekist->Nunits() > 0) _DirtyDivider(imageinfo[i].imagekist,wholekist,Ngroup);
		++(*Nimages);
		++i;
	}

	assert(*Nimages < MaxNimages);

	delete wholekist;
	return;
}
void _DirtyDivider(KistHndl neighbors,KistHndl wholekist,int Ngroup){
	bool check = false;
	double xinit[2],rr[Ngroup+1],ll;
	long i,j;

	if(wholekist->Nunits() <= Ngroup){
		while(wholekist->Nunits() > 0) InsertBeforeCurrentKist(neighbors,TakeOutCurrentKist(wholekist));
		return;
	}

	xinit[0] = getCurrentKist(neighbors)->x[0];
	xinit[1] = getCurrentKist(neighbors)->x[1];

	for(i=0;i<Ngroup+1;++i) rr[i] = 1.0e200;

	MoveToTopKist(wholekist);
	do{
		if(check){
			MoveUpKist(wholekist);
			check=false;
		}

		ll = pow(xinit[0]-getCurrentKist(wholekist)->x[0],2)
				+ pow(xinit[1]-getCurrentKist(wholekist)->x[1],2);

		if(ll < rr[Ngroup]){
			MoveToTopKist(neighbors);
			for(i=0;i<Ngroup+1;++i,MoveDownKist(neighbors)){
				if(ll < rr[i]){
					InsertBeforeCurrentKist(neighbors,TakeOutCurrentKist(wholekist));
					for(j=Ngroup;j>i;--j) rr[j] = rr[j-1];
					rr[i] = ll;
					break;
				}
			}
			if(AtTopKist(wholekist)) check=true;
		}
	}while(MoveDownKist(wholekist));

	// put nearby points at the beginning of wholekist
	MoveToBottomKist(neighbors);
	MoveToTopKist(wholekist);
	while(neighbors->Nunits() > Ngroup){
		InsertBeforeCurrentKist(wholekist,TakeOutCurrentKist(neighbors));
		MoveToTopKist(wholekist);
	}

	return ;
}
