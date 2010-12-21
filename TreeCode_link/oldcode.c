/*
 * oldcode.c
 *
 *  Created on: Jan 15, 2010
 *      Author: R.B. Metcalf
 *
 *      this code is no longer used because it has been superseded
 *      by replacement code.  It was fully functional before moving here.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../../Library/Recipes/nrutil.h"
#include "../../Library/RecipesD/nrD.h"
#include "Tree.h"
#define pi  3.1415926


void findborders(TreeHndl i_tree,ImageInfo *imageinfo){
/*   finds borders of images by using nearest neighbor routine
 *      it has been replaced by findborders2 which uses cell neighbors.
 *      this routine still works when there is more than 1 point in each
 *      leaf in the tree
 */

  double rmin,rmax,r[8],x[2];
  short mark;
  int j,k,m,larger,smaller;
  ListHndl neighborlist;
  Point *point;
  short dir[8],sum;
  time_t to,t1;

  time(&to);
  /*printf("In findborders\n");*/

  neighborlist=NewList();
  EmptyList(imageinfo->innerborder);
  EmptyList(imageinfo->outerborder);

  imageinfo->gridrange[2] = 1.0e99; /* minimum grid size in image */
  imageinfo->gridrange[0] = 0.0; /* maximum grid size in outerborder */
  imageinfo->gridrange[1] = 0.0;      /* maximum grid size in image */

  for(j=0;j<imageinfo->Npoints;++j){
    mark=0;

    if(imageinfo->gridrange[1] < imageinfo->points[j].gridsize)
      imageinfo->gridrange[1] = imageinfo->points[j].gridsize;
    if(imageinfo->gridrange[2] > imageinfo->points[j].gridsize)
      imageinfo->gridrange[2] = imageinfo->points[j].gridsize;

    point=NearestNeighbor(i_tree,imageinfo->points[j].x,9,neighborlist,0);

    /* determine if point j is on the edge of a grid */
    MoveToTopList(neighborlist);
    point=TakeOutCurrent(neighborlist); /* take out jth point itself */
    free(point);

    for(k=0,rmin=1.0e99,rmax=0,larger=0,smaller=0;k<8;++k){ /* loop through neighbors to point j */

      r[k]=sqrt(pow(imageinfo->points[j].x[0]-neighborlist->current->x[0],2)
	     + pow(imageinfo->points[j].x[1]-neighborlist->current->x[1],2));
      /*printf("in findborders\n");*/
      if(r[k]<rmin) rmin=r[k];
      if(r[k]>rmax) rmax=r[k];

      if(neighborlist->current->gridsize > 1.01*imageinfo->points[j].gridsize) ++larger;
      if(neighborlist->current->gridsize < 0.99*imageinfo->points[j].gridsize) ++smaller;

      MoveDownList(neighborlist);
    }
    for(k=0;k<8;++k) r[k]/=rmin;

    if(smaller){ /* on border of larger grid */

      EmptyList(neighborlist);
      NearestNeighbor(i_tree,imageinfo->points[j].x,1,neighborlist,1);
      NearestNeighbor(i_tree,imageinfo->points[j].x,1,neighborlist,2);
      NearestNeighbor(i_tree,imageinfo->points[j].x,1,neighborlist,3);
      NearestNeighbor(i_tree,imageinfo->points[j].x,1,neighborlist,4);

    }else if( fabs( (r[4]+r[5]+r[6]+r[7])/4 - 1.41421 ) < 1.0e-2){ /* in grid */


    }else{ /* edge of smaller grid */

/*       printf("Hello on edge %e rsecond %e j=%i\n",rmin,rsecond,j); */
/*       for(m=0;m<8;++m) printf("r[%i]=%e\n",m,r[m]);  */
      /* find which grid neighbors are missing */

      for(m=0;m<8;++m) dir[m]=0;
      MoveToTopList(neighborlist);
      for(m=0;m<neighborlist->Npoints;++m){
    	  dir[0]+=(sqrt(pow(imageinfo->points[j].x[0]+rmin-neighborlist->current->x[0],2)
		      +pow(imageinfo->points[j].x[1]-neighborlist->current->x[1],2)) < 1.0e-3*rmin);
    	  dir[1]+=(sqrt(pow(imageinfo->points[j].x[0]+rmin-neighborlist->current->x[0],2)
		      +pow(imageinfo->points[j].x[1]+rmin-neighborlist->current->x[1],2)) < 1.0e-3*rmin);
    	  dir[2]+=(sqrt(pow(imageinfo->points[j].x[0]-neighborlist->current->x[0],2)
		      +pow(imageinfo->points[j].x[1]+rmin-neighborlist->current->x[1],2)) < 1.0e-3*rmin);
    	  dir[3]+=(sqrt(pow(imageinfo->points[j].x[0]-rmin-neighborlist->current->x[0],2)
		      +pow(imageinfo->points[j].x[1]+rmin-neighborlist->current->x[1],2)) < 1.0e-3*rmin);
    	  dir[4]+=(sqrt(pow(imageinfo->points[j].x[0]-rmin-neighborlist->current->x[0],2)
		      +pow(imageinfo->points[j].x[1]-neighborlist->current->x[1],2)) < 1.0e-3*rmin);
    	  dir[5]+=(sqrt(pow(imageinfo->points[j].x[0]-rmin-neighborlist->current->x[0],2)
		      +pow(imageinfo->points[j].x[1]-rmin-neighborlist->current->x[1],2)) < 1.0e-3*rmin);
    	  dir[6]+=(sqrt(pow(imageinfo->points[j].x[0]-neighborlist->current->x[0],2)
		      +pow(imageinfo->points[j].x[1]-rmin-neighborlist->current->x[1],2)) < 1.0e-3*rmin);
    	  dir[7]+=(sqrt(pow(imageinfo->points[j].x[0]+rmin-neighborlist->current->x[0],2)
		      +pow(imageinfo->points[j].x[1]-rmin-neighborlist->current->x[1],2)) < 1.0e-3*rmin);
    	  MoveDownList(neighborlist);
      }

      for(m=0,sum=0;m<8;++m) sum+=dir[m];
      if(sum == 8){
    	  printf("ERROR in findborders, false grid edge \n");
    	  printf("center %i %e %e\n",imageinfo->points[j].id,imageinfo->points[j].x[0],imageinfo->points[j].x[1]);
    	  PrintList(neighborlist);
    	  for(m=0;m<8;++m) printf("r[%i]=%e\n",m,r[m]);
    	  exit(1);
      }

      if(sum == 7){ /* inner corner */
    	  MoveToBottomList(neighborlist);
    	  point=TakeOutCurrent(neighborlist);
    	  free(point);

    	  x[0]=( (1-dir[1]) - (1-dir[3]) - (1-dir[5]) + (1-dir[7]) )*2*rmin;
    	  x[1]=( (1-dir[1]) + (1-dir[3]) - (1-dir[5]) - (1-dir[7]) )*2*rmin;

    	  x[0]+=imageinfo->points[j].x[0];
    	  x[1]+=imageinfo->points[j].x[1];

    	  NearestNeighbor(i_tree,&x[0],1,neighborlist,-1);
      }

      if(sum == 6){ /* next to inner corner */

	MoveToBottomList(neighborlist);
	point=TakeOutCurrent(neighborlist);
	free(point);
	point=TakeOutCurrent(neighborlist);
	free(point);

	x[0]=( (1-dir[1]) - (1-dir[3]) - (1-dir[5]) + (1-dir[7]) )*rmin
	  + ( (1-dir[0]) - (1-dir[4])  )*rmin;
	x[1]=( (1-dir[1]) + (1-dir[3]) - (1-dir[5]) - (1-dir[7]) )*rmin
	  + ( (1-dir[2]) - (1-dir[6])  )*rmin;

	x[0]+=imageinfo->points[j].x[0];
	x[1]+=imageinfo->points[j].x[1];

	NearestNeighbor(i_tree,x,1,neighborlist,-1);
      }
      if(sum == 5){ /* straight edge */

	EmptyList(neighborlist);
	NearestNeighbor(i_tree,imageinfo->points[j].x,1,neighborlist,1);
	NearestNeighbor(i_tree,imageinfo->points[j].x,1,neighborlist,2);
	NearestNeighbor(i_tree,imageinfo->points[j].x,1,neighborlist,3);
	NearestNeighbor(i_tree,imageinfo->points[j].x,1,neighborlist,4);
      }
      if(sum == 3){ /* outer corner */

	EmptyList(neighborlist);
	NearestNeighbor(i_tree,imageinfo->points[j].x,1,neighborlist,1);
	NearestNeighbor(i_tree,imageinfo->points[j].x,1,neighborlist,2);
	NearestNeighbor(i_tree,imageinfo->points[j].x,1,neighborlist,3);
	NearestNeighbor(i_tree,imageinfo->points[j].x,1,neighborlist,4);

	x[0]=( (1-dir[1]) - (1-dir[3]) - (1-dir[5]) + (1-dir[7]) )*2*rmin;
	x[1]=( (1-dir[1]) + (1-dir[3]) - (1-dir[5]) - (1-dir[7]) )*2*rmin;

	x[0]+=imageinfo->points[j].x[0];
	x[1]+=imageinfo->points[j].x[1];

	NearestNeighbor(i_tree,x,1,neighborlist,-1);/**/
      }

    }

    /*PrintList(neighborlist);*/
    /* see if it is a border point */
    MoveToTopList(neighborlist);
    for(k=0;k<neighborlist->Npoints;++k){ /* loop through neighbors to point j */
      for(m=0;m<imageinfo->Npoints;++m) if(imageinfo->points[m].id == neighborlist->current->id) break;

	/* current is a neighbor to the image */
    if(m == imageinfo->Npoints){
    	mark=1;
	  /* check to make sure this is not already in the borderlist */
	MoveToTopList(imageinfo->outerborder);
	for(m=0;m<imageinfo->outerborder->Npoints;++m){
	  if(imageinfo->outerborder->current->id == neighborlist->current->id) break;
	  MoveDownList(imageinfo->outerborder);
	}

	if(m == imageinfo->outerborder->Npoints){
	  InsertAfterCurrent(imageinfo->outerborder,neighborlist->current->x
			     ,neighborlist->current->id
			     ,neighborlist->current->image);
	  MoveDownList(imageinfo->outerborder);
	  PointCopyData(imageinfo->outerborder->current,neighborlist->current);
	  if(imageinfo->gridrange[0] < neighborlist->current->gridsize)
	    imageinfo->gridrange[0] = neighborlist->current->gridsize;
 	  if(imageinfo->gridrange[2] > neighborlist->current->gridsize)
 	    imageinfo->gridrange[2] = neighborlist->current->gridsize;
	}

      }
      MoveDownList(neighborlist);
    }
    /* enter point as inner border */
    if(mark){
      /* make sure it isn't already in list */
      MoveToTopList(imageinfo->innerborder);
      for(m=0;m<imageinfo->innerborder->Npoints;++m){
	if(imageinfo->innerborder->current->id == imageinfo->points[j].id) break;
	MoveDownList(imageinfo->innerborder);
      }

      if(m == imageinfo->innerborder->Npoints){
	InsertAfterCurrent(imageinfo->innerborder,imageinfo->points[j].x
			   ,imageinfo->points[j].id,imageinfo->points[j].image);
	MoveDownList(imageinfo->innerborder);
	PointCopyData(imageinfo->innerborder->current,&(imageinfo->points[j]));
      }
    }

  }/* loop through image points */


/*   printf("Nouterborder=%i Ninnerborder=%i\n",imageinfo->outerborder->Npoints,imageinfo->innerborder->Npoints); */
/*   exit(0); */
/* */
  EmptyList(neighborlist);  // un-commented not sure why it was commented
  free(neighborlist);

  /*printf(" grid range on leaving find_borders %e %e\n",imageinfo->gridrange[0],imageinfo->gridrange[1]);*/
  time(&t1);
  // printf("                    time in findborders:  %f sec\n",difftime(t1,to));
}
