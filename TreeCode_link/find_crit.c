/*
 * find_crit.c
 *
 *  Created on: Sep 8, 2009
 *      Author: R.B. Metcalf
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../../Library/Recipes/nrutil.h"
#include "Tree.h"
#define NMAXCRITS 100

ImageInfo *find_crit(TreeHndl s_tree,TreeHndl i_tree,int *Ncrits,double resolution
		,Boolean *orderingsuccess,Boolean ordercurve,Boolean verbose){
  /* resolution is the resolution on the image plane */
  /* the inner out outer boundaries of the result are the estimated critical curves */
  /* OUTPUT: each critical curve is in a array of IamgeInfo's    */
  /*         result.parity = 1 tangential caustic, 2 radial, 0 not enough points to determine */

  Point *minpoint;
  ImageInfo *critcurve,*critexport;
  //unsigned long j,k,m,jold;
  unsigned long Npoints,i=0;
  short refinements;
  //short spur,closed;
  double maxgridsize,mingridsize,x[2];
  ListHndl negpointlist;

  negpointlist=NewList();
  critcurve=NewImageInfo(NMAXCRITS);
  minpoint=NewPoint(x,0);
  minpoint->invmag=1.0e99;
  /*point=NmewPoint(x,0);*/

  for(;;){

	  // find list of points with negative magnification
	  EmptyList(negpointlist);
	  MoveToTopList(i_tree->pointlist);
	  for(i=0,minpoint->kappa=0;i<i_tree->pointlist->Npoints;++i){
		  if(i_tree->pointlist->current->invmag < 0){

			  InsertAfterCurrent(negpointlist,i_tree->pointlist->current->x,i_tree->pointlist->current->id
					  ,i_tree->pointlist->current->image);
			  MoveDownList(negpointlist);
			  PointCopyData(negpointlist->current,i_tree->pointlist->current);
		  }

		  // record point of maximum kappa
		  if(i_tree->pointlist->current->kappa > minpoint->kappa) PointCopyData(minpoint,i_tree->pointlist->current);
		  MoveDownList(i_tree->pointlist);
	  }

	  Npoints=negpointlist->Npoints;
	  critcurve[0].Npoints=Npoints;

	  if(Npoints == 0){
		  if(minpoint->gridsize <= resolution){  // no caustic found at this resolution
			  *Ncrits=0;
			  return critcurve;
		  }

		  /* if there is no negative magnification points use maximum mag point */
		  ++Npoints;
		  //critcurve[0].points=(Point *) malloc(sizeof(Point));
		  //critcurve[0].points->head=1;
		  critcurve[0].points=NewPointArray(1,False);
		  critcurve[0].points->in_image=False;
		  critcurve[0].Npoints=1;
		  PointCopyData(critcurve[0].points,minpoint);
	  }else{
		  /*critcurve[0].points=(Point *) malloc(Npoints*sizeof(Point));
		  critcurve[0].Npoints=Npoints;
		  critcurve[0].points->head=Npoints;
		   */

		  critcurve[0].points=NewPointArray(Npoints,False);

		  MoveToTopList(negpointlist);
		  for(i=0;i<negpointlist->Npoints;++i){
			  PointCopyData(&(critcurve[0].points[i]),negpointlist->current);
    	  //printf("     negpointlist = %e %e \n       critcurve[0].point[%i]=%e %e\n",negpointlist->current->x[0]
    	 // 	       ,negpointlist->current->x[1],i,critcurve[0].points[i].x[0],critcurve[0].points[i].x[1]);
			  MoveDownList(negpointlist);
		  }
	  }

    if(verbose) printf("find_crit, going into findborders 1\n");
    findborders2(i_tree,&critcurve[0]);
    if(verbose) printf("find_crit, came out of findborders 1\n");

    /* make inner border the image */
    MoveToTopKist(critcurve[0].innerborder);
    for(i=0,maxgridsize=0.0,mingridsize=1.0e99;i<critcurve[0].innerborder->Nunits;++i){
    	PointCopyData(&(critcurve[0].points[i]),getCurrentKist(critcurve[0].innerborder));
    	if(critcurve[0].points[i].gridsize > maxgridsize) maxgridsize=critcurve[0].points[i].gridsize;
    	if(critcurve[0].points[i].gridsize < mingridsize) mingridsize=critcurve[0].points[i].gridsize;
    	MoveDownKist(critcurve[0].innerborder);
    }
    critcurve[0].Npoints=critcurve[0].innerborder->Nunits;

    // find borders again to properly define outer border
    //printf("going into findborders 2\n");
	findborders2(i_tree,critcurve);
	//printf("came out of findborders 2\n");

	if(verbose) printf("find_crit, going into refine_grid\n");
     //printf("  Npoints=%i\n",critcurve->Npoints);
	refinements=refine_grid(i_tree,s_tree,critcurve,1,resolution,2,False);
    if(verbose) printf("find_crit, came out of refine_grid\n");

     if(refinements==0){
      break;
    }else free(critcurve[0].points);
  }

  if(verbose) printf("find_crit, number of caustic points: %li\n",critcurve[0].Npoints);

// order points in curve
  if(ordercurve) split_order_curve4(critcurve,NMAXCRITS,Ncrits);
  else if(critcurve->Npoints > 0) *Ncrits=1;
  if(critcurve->Npoints == 0) *Ncrits=0;

/*
//   print out the critical curves and caustics
  printf("Ncrits=%i\n",*Ncrits);
  for(j=0;j<*Ncrits;++j){
	printf("%li\n",critcurve[j].Npoints);
	for(i=0;i<critcurve[j].Npoints;++i)
		printf("%e %e\n",critcurve[j].points[i].x[0]
	                    ,critcurve[j].points[i].x[1]);
  }
  printf("Ncrits=%i\n",*Ncrits);
  for(j=0;j<*Ncrits;++j){
	printf("%li\n",critcurve[j].Npoints);
	for(i=0;i<critcurve[j].Npoints;++i)
		printf("%e %e\n",critcurve[j].points[i].image->x[0]
	                    ,critcurve[j].points[i].image->x[1]);
  }
  exit(0);
	*/

  if(*Ncrits==0 && critcurve->Npoints > 0 ){
	  *Ncrits=1;
	  *orderingsuccess=False;
  }else{ *orderingsuccess=True;}

  if(*Ncrits > NMAXCRITS){ERROR_MESSAGE(); printf("ERROR: in find_crit, too many critical curves Ncrits=%i > NMAXCRITS gridsize=%e\n"
			       ,*Ncrits,critcurve[0].points[0].gridsize); exit(1);}

  /* find area of critical curves */
  x[0]=x[1]=0.0;
  for(i=0;i<*Ncrits;++i){
	  if(critcurve[i].Npoints > 5){
			   windings(x,critcurve[i].points,critcurve[i].Npoints,&(critcurve[i].area),0);
			   //printf("critarea = %e\n",critcurve[i].area);
	  }else critcurve[i].area=0;
  }

  EmptyList(negpointlist);
  free(negpointlist);
  free(minpoint);

  // resize crit array so it doesn't use more mem than necessary
  critexport=NewImageInfo(*Ncrits);
  for(i=0;i<*Ncrits;++i){
	  critexport[i].Nencircled=critcurve[i].Nencircled;
	  critexport[i].Npoints=critcurve[i].Npoints;
	  critexport[i].area=critcurve[i].area;
	  critexport[i].area_error=critcurve[i].area_error;
	  critexport[i].points=critcurve[i].points;
  }
  freeImageInfo(critcurve,NMAXCRITS);
  return critexport;
}

