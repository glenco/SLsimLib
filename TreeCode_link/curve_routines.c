/*
 * curve_routines.c
 *
 *  Created on: Jan 15, 2010
 *      Author: R.B. Metcalf
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <omp.h>
#include <nrutil.h>
#include <nrD.h>
#include "Tree.h"
#ifndef pi
#define pi  3.1415926
#endif

void split_order_curve4(ImageInfo *curves,int Maxcurves,int *Ncurves){
/*  orders points in a curve, separates disconnected curves
 *   curves[0...Maxcurves] must be allocated before
 *
 *	uses neighbors-of-neighbors to split into curves and then
 *	uses sorts by angle and then walks the list sorting
 *
 *  can break down for crescent curves
 */

	long i,m,j,end;
	//short spur,closed,attach;
	unsigned long NpointsTotal;
	double center[2],*theta;
	//Boolean delta,tmp,step;
	//ListHndl reservoir,orderedlist;
	Point *newpointarray;

	//printf("entering split_order_curve\n");
	if(curves[0].Npoints==0){
		*Ncurves=0;
		return;
	}

	NpointsTotal=curves[0].Npoints;

	splitter(curves,Maxcurves,Ncurves);

	/*
	reservoir=NewList();
	orderedlist=NewList();

	// copy points into a list
	for(i=0;i<NpointsTotal;++i) InsertPointAfterCurrent(reservoir,&(curves[0].points[i]));
	m=0;
	i=0;

	// divide points into disconnected curves using neighbors-of-neighbors
	while(reservoir->Npoints > 0 && i < Maxcurves){
		curves[i].points=TakeOutCurrent(reservoir);
		MoveToBottomList(orderedlist);
		InsertPointAfterCurrent(orderedlist,curves[i].points);
		MoveDownList(orderedlist);
		NeighborsOfNeighbors(orderedlist,reservoir);
		curves[i].Npoints=orderedlist->Npoints - m;
		//printf("curves[%i].Npoints=%i %i %i\n",i,curves[i].Npoints,orderedlist->Npoints
		//		,reservoir->Npoints);
		m+=curves[i].Npoints;
		++i;
	}
	*Ncurves=i;

	// copy list back into array
	point=curves[0].points;
	newpointarray=NewPointArray(NpointsTotal,False);
	MoveToTopList(orderedlist);
	m=0;
	i=0;
	do{
		if(i < *Ncurves && curves[i].points==orderedlist->current){
			curves[i].points=&(newpointarray[m]);
			++i;
		}
		PointCopyData(&(newpointarray[m]),orderedlist->current);
		++m;
	}while(MoveDownList(orderedlist));

	assert(m == NpointsTotal);
*/
	// order curve points
	theta=(double *)malloc(NpointsTotal*sizeof(double));
	assert(theta);
	for(i=0;i<*Ncurves;++i){

		if(curves[i].Npoints > 3){
			// sort points by angle around center point
			center[0]=center[1]=0;
			for(m=0;m<curves[i].Npoints;++m){
				center[0]+=curves[i].points[m].x[0];
				center[1]+=curves[i].points[m].x[1];
			}
			center[0]/=curves[i].Npoints;
			center[1]/=curves[i].Npoints;

			for(m=0;m<curves[i].Npoints;++m){
				theta[m]=atan2(curves[i].points[m].x[1]-center[1]
			              ,curves[i].points[m].x[0]-center[0]);
			}
			quicksortPoints(curves[i].points,theta,curves[i].Npoints);

			// check to make sure the center is inside the curves
			//assert(abs(windings(center,curves[i].points,curves[i].Npoints,&tmp1,0)));

			//printf("N=%i\n",curves[i].Npoints);
			//assert(abs(windings(center,curves[i].points,curves[i].Npoints,&tmp1,0)) > 0 );

			// find the last point that is a neighbor of the first point
			m=curves[i].Npoints-1;
			while(!AreBoxNeighbors(&(curves[i].points[0]),&(curves[i].points[m]))) --m;
			end=m;
			//assert(m > 0);

		//printf("windings = %i\n",windings(center,curves[i].points,curves[i].Npoints,&tmp1,0));
		// walk curve to remove shadowing effect
			j=0;
		//end=0;
			m=0;
			while(j < curves[i].Npoints-1){
				walkcurve(curves[i].points,curves[i].Npoints,&j,&end);
				//printf("i = %i Npoints = %i end+1 = %i j = %i\n",i,curves[i].Npoints,end+1,j);
				if(j < curves[i].Npoints-1)
					backtrack(curves[i].points,curves[i].Npoints,&j,-1,&end);
				++m;
				assert(m < curves[i].Npoints);
			}
			//printf("i = %i Npoints = %i end+1 = %i j=%i\n",i,curves[i].Npoints,end+1,j);
			curves[i].Npoints=end+1;
		}
	}

	free(theta);
	//free(reservoir);
	//free(orderedlist);
	//free(point);

	return ;
}

void split_order_curve3(ImageInfo *curves,int Maxcurves,int *Ncurves){
/*  orders points in a curve, separates disconnected curves
 *   curves[0...Maxcurves] must be allocated before
 *
 *	 uses list instead of arrays
 *   does not attempt to remove spurs
 */

	long i,k=0,m;
	short spur,closed,attach;
	unsigned long Npoints;
	Boolean delta,tmp,step;
	ListHndl lists[Maxcurves+1];
	Point *point,*newpointarray;

	lists[Maxcurves]=NewList();
	for(i=0;i<curves[0].Npoints;++i) InsertPointAfterCurrent(lists[Maxcurves],&(curves[0].points[i]));

	//printf("entering split_order_curve\n");
	if(curves[0].Npoints==0){
		*Ncurves=0;
		return;
	}

	// separate critical curves
	Npoints=curves[0].Npoints;
	//printf("number of curve points: %i\n",curves[0].Npoints);

	for(k=0;lists[Maxcurves]->Npoints > 0;++k){
		printf("lists[Maxcurves]->Npoints=%li  \n",lists[Maxcurves]->Npoints);
		lists[k]=NewList();
	    if(k==Maxcurves){ ERROR_MESSAGE(); printf("ERROR: in split_order_curve3, not large enough curves array k=%li\n",k); break;}

		// start new curve
		point=TakeOutCurrent(lists[Maxcurves]);
 		InsertPointAfterCurrent(lists[k],point);
		MoveDownList(lists[k]);

	    // go through remaining points ordering them as we go down the list
	    do{
	    	spur=0;
	    	MoveToTopList(lists[Maxcurves]);
	    	tmp=True;
	    	step=True;

	    	if(lists[Maxcurves]->Npoints > 0){
	    		do{
	    			if(step){
	    				delta = sqrt( pow(lists[k]->current->x[0]-lists[Maxcurves]->current->x[0],2)
		    			        + pow(lists[k]->current->x[1]-lists[Maxcurves]->current->x[1],2) )
		    			        < 1.01*(lists[k]->current->gridsize + lists[Maxcurves]->current->gridsize)/2;
	    			}else{
	    				delta = AreBoxNeighbors(lists[k]->current,lists[Maxcurves]->current);
	    			}

	    			//printf("move down\n");
	    			if(delta){
	    				point=TakeOutCurrent(lists[Maxcurves]);
	    				MoveToTopList(lists[Maxcurves]);
	    				InsertPointAfterCurrent(lists[k],point);
	    				MoveDownList(lists[k]);
	    				step=True;
	    			}else{
	    				tmp=MoveDownList(lists[Maxcurves]);
	    			}

	    			if( tmp==False && step){
	    				MoveToTopList(lists[Maxcurves]);
	    				tmp=True;
	    				step=False;
	    			}

	    		}while(tmp && lists[Maxcurves]->Npoints > 0);
	    	}

	    	// check if curve is closed
	    	if( lists[k]->Npoints < 3 || AreBoxNeighbors(lists[k]->top,lists[k]->bottom ) == False ){
			  // curve is not closed

	    		closed=0;

	    		// work backward along curve to find point with another neighbor
	    		if(lists[Maxcurves]->Npoints > 0){
	    			for(m=-1;(m>-lists[k]->Npoints && spur==0);--m){

	    				MoveUpList(lists[k]);
	    				// reverse the order of the spur
	    				// this insures that if the code returns to the start
	    				// the points will still be in a order
	    				/*JumpDownList(lists[k],m);
	    				MoveCurrentToBottom(lists[k]);
	    				MoveToBottomList(lists[k]);
	*/
	    				MoveToTopList(lists[Maxcurves]);
	    				do{
	    					if( AreBoxNeighbors(lists[k]->current,lists[Maxcurves]->current) ){
	    						point=TakeOutCurrent(lists[Maxcurves]);
	    						InsertPointAfterCurrent(lists[k],point);
	    						MoveDownList(lists[k]);
	    						spur=1;
	    					}
	    				}while(spur==0 && MoveDownList(lists[Maxcurves]));
	    			}
	    		}

	    		// if spur is not found see if curve is attached to a previous curves
	    		if(spur==0){
	    			attach = 0;
	    			for(m=0;m<k;++m){
	    				if(lists[m]->Npoints > 0){
	    					MoveToTopList(lists[m]);
	    					do{
	    						if( AreBoxNeighbors(lists[m]->current,lists[k]->bottom) ){
	    							ShiftList(lists[m]);
	    							MoveToBottomList(lists[k]);
	    							InsertListAfterCurrent(lists[k],lists[m]);
	    							MoveToBottomList(lists[k]);
	    							//MoveToTopList(lists[Maxcurves]);
	    							//InsertListBeforeCurrent(lists[Maxcurves],lists[m]);
	    							lists[m]->Npoints=0;
	    							lists[m]->top = lists[m]->bottom = lists[m]->current = NULL;
	    							attach=1;
	    							spur=1;
	    						}
	    					}while(attach == 0 && MoveDownList(lists[m]) );
	    				}
	    			}
	    		}

	      }else{closed=1;}

	    }while(spur && lists[Maxcurves]->Npoints > 0);

	  }

	//printf("  end of loop j=%i k=%i i=%i Npoints=%i\n",j,k,i,Npoints);

	point=curves[0].points;
	newpointarray=NewPointArray(curves[0].Npoints,False);
	for(i=0,*Ncurves=0,m=0;i<k;++i){
		printf("lists[%li]->Npoints = %li\n",i,lists[i]->Npoints);
		if(lists[i]->Npoints > 0){
			MoveToTopList(lists[i]);
			PointCopyData(&(newpointarray[m]),lists[i]->current);
			curves[*Ncurves].points=&(newpointarray[m]);
			curves[*Ncurves].Npoints=lists[i]->Npoints;
			++m;

			while(MoveDownList(lists[i])){
				PointCopyData(&(newpointarray[m]),lists[i]->current);
				++m;
			}

			*Ncurves=*Ncurves+1;
		}
	}
	for(i=0;i<k;++i) free(lists[i]);
	free(lists[Maxcurves]);
	free(point);
	return ;
}

void split_order_curve(ImageInfo *curves,int Maxcurves,int *Ncurves){
/*  orders points in a curve, separates disconnected curves
 *   curves[0...Maxcurves] must be allocated before
 *
 *   does not attempt to remove spurs
 */

	long j,k=0,jold,end=0;
	short spur,closed;
	unsigned long Npoints,Maxpoint;

	//printf("entering split_order_curve\n");
	if(curves[0].Npoints==0){
		*Ncurves=0;
		return;
	}

	// separate critical curves
	Npoints=curves[0].Npoints;
	//printf("number of curve points: %i\n",curves[0].Npoints);

	for(j=0,k=0,jold=0;j<Npoints;++k){
	    if(k==Maxcurves){ ERROR_MESSAGE();printf("ERROR: in split_order_curve, not large enough curves array k=%li j=%li\n",k,j); break;}
	    Maxpoint=Npoints;

	    // go through remaining points ordering them as we go down the list
	    do{
	    	spur=0;
	    	walkcurve(curves[0].points,Maxpoint,&j,&end);

	      // check if curve is closed

		  if( AreBoxNeighbors(&(curves[k].points[0]),&(curves[0].points[j]) ) == False &&
				  (j != jold && j < Maxpoint) ){
			  // curve is not closed

			  closed=0;
			  // work backward along curve to find point with another neighbor
			  spur=backtrack(curves[0].points,Maxpoint,&j,jold,&end);

		  }else{closed=1;}

	    }while(spur);

		++j;
		curves[k].Npoints=j-jold;
		jold=j;
		if(j < Npoints && k+1<Maxcurves) curves[k+1].points=&(curves[0].points[j]);
	}


	//printf("  end of loop j=%i k=%i i=%i Npoints=%i\n",j,k,i,Npoints);

	*Ncurves=k;
	//printf("exiting split_order_curve\n");

}

void split_order_curve2(ImageInfo *curves,int Maxcurves,int *Ncurves){
/*  orders points in a curve, separates disconnected curves, cuts off spurs
 *   curves[0...Maxcurves] must be allocated before
 *
 *   Attempts to remove spurs,  does not always work
 */

	long i,j,k=0,jold,m,end;
	short spur,closed;
	unsigned long Npoints,Maxpoint;

	//printf("entering split_order_curve\n");
	if(curves[0].Npoints==0){
		*Ncurves=0;
		return;
	}

	// separate critical curves
	Npoints=curves[0].Npoints;
	//printf("number of curve points: %i\n",curves[0].Npoints);

	for(j=0,k=0,jold=0;j<Npoints;++k){
	    if(k==Maxcurves){ ERROR_MESSAGE(); printf("ERROR: in split_order_curve2, not large enough curves array k=%li j=%li\n",k,j); break;}
	    Maxpoint=Npoints;

	    // go through remaining points ordering them as we go down the list
	    do{

	    	spur=0;
	    	walkcurve(curves[0].points,Maxpoint,&j,&end);

	    	/*
	    	for(i=j+1,step=1,rmin=1.0e99;i<Maxpoint;++i){
	    		if( j+1 > Npoints-1 || i==j){ ERROR_MESSAGE(); printf("ERROR: i=%l j=%i\n",i,j); exit(1);}

	    		if(step){
	    			delta = sqrt( pow(curves[0].points[j].x[0]-curves[0].points[i].x[0],2)
		    			  + pow(curves[0].points[j].x[1]-curves[0].points[i].x[1],2) )
		    			  < 1.01*(curves[0].points[j].gridsize + curves[0].points[i].gridsize)/2;
	    		}else{
	    			delta = AreBoxNeighbors(&(curves[0].points[i]),&(curves[0].points[j]) );
	    		}

	    		if(delta){
	    			SwapPointsInArray( &(curves[0].points[j+1]) , &(curves[0].points[i]) );
	    			++j;
	    			i=j;
	    			step=1;
	    		}

	    		// try larger linking length
	    		if( (i == Maxpoint-1) && step){
	    			i=j;
	    			step=0;
	    		}
	    	}*/

	      // check if curve is closed

		  if( AreBoxNeighbors(&(curves[k].points[0]),&(curves[0].points[j]) ) == False &&
				  (j != jold && j < Maxpoint) ){
			  // curve is not closed

			  closed=0;

	    	  // work backward along curve to find point with another neighbor
			  //spur=backtrack(curves[0].points,Maxpoint,&j,jold,&end);

	    	  for(m=j-1;(m>jold && spur==0);--m){
	    		  for(i=j+1;(i<Maxpoint && spur==0);++i){

	    			  if(AreBoxNeighbors(&(curves[0].points[i]),&(curves[0].points[m])) ){
	    				  // add neighbor to end of curve and restart
	    				  SwapPointsInArray(&(curves[0].points[m+1]),&(curves[0].points[i]) );
	    				  spur=1;
	    	    	  }
	    		  }
 	    	  }

	    	  // move spur out of list for this k-th image
	    	  if(spur){
		    	  // move last point out of possible curve list
	    		  --Maxpoint;
    			  SwapPointsInArray( &(curves[0].points[i-1]),&(curves[0].points[Maxpoint]) );
	    		  for(i=m+3;i<=j;++i){
	    			  --Maxpoint;
	    			  SwapPointsInArray( &(curves[0].points[i]),&(curves[0].points[Maxpoint]) );
	    		  }
	    		  j=m+2;
	    	  }

	      }else{closed=1;}
	    }while(spur);

	    // prune spurs in curve k
	    // Note: this may make separating close images inefficient
		Npoints=Maxpoint;

	    // prune open curves by moving them out of list for all curves
	    if(closed==0 || j-jold < 2 ){
	    	// move spur to end of array, reset to end of last curve and reset total number of points
	    	for(i=jold,m=Npoints-1;i<=j;++i,--m){
		    	//printf("%i %e %e\n",i,curves[0].points[i].x[2],curves[0].points[i].x[1]);
	    		SwapPointsInArray(&(curves[0].points[i]),&(curves[0].points[m]) );
	    	}
	    	--k;
	    	Npoints-=j-jold+1;
	    	j=jold;
		}else{
	    	++j;
	    	curves[k].Npoints=j-jold;
	    	jold=j;
	    	if(j < Npoints && k+1<Maxcurves) curves[k+1].points=&(curves[0].points[j]);
	    }
	    // printf("2 Npoints=%i\n",curves[k].Npoints);
	  }
	//printf("  end of loop j=%i k=%i i=%i Npoints=%i\n",j,k,i,Npoints);

	*Ncurves=k;
	//printf("exiting split_order_curve\n");

}

void order_curve(ImageInfo *curve){
/*  orders points in a curve,
 *
 *   must already be a connected set of points
 *
 *   does not attempt to remove spurs
 */

	long j,jold,end=0;
	short spur,closed;

	//printf("entering split_order_curve\n");
	if(curve->Npoints==0) return;

	// separate critical curves
	//printf("number of curve points: %i\n",curves[0].Npoints);

	j=0; jold=0;

	// go through remaining points ordering them as we go down the list
	do{
		spur=0;
		walkcurve(curve->points,curve->Npoints,&j,&end);

		// check if curve is closed
		if( AreBoxNeighbors(&(curve->points[0]),&(curve->points[j]) ) == False &&
				(j != 0 && j < curve->Npoints) ){
				  // curve is not closed

			closed=0;
			// work backward along curve to find point with another neighbor
			spur=backtrack(curve->points,curve->Npoints,&j,0,&end);

		}else{closed=1;}

	}while(spur);

	//printf("  end of loop j=%i k=%i i=%i Npoints=%i\n",j,k,i,Npoints);
	//printf("exiting split_order_curve\n");

}

void walkcurve(Point *points,long Npoints,long *j,long *end){
	/*
	 * orders curve by finding closest point ahead in list
	 *   tries x/y jump before diagonal jump
	 *   exits when it cannot find a cell neighbor ahead in array
	 *   leaves j as the last point in curves
	 */

	long i,k;
	short step;
	Boolean delta;

	//if((*j)==0) *end=Npoints-1;
	//printf("end = %i\n",*end);

	for(i=(*j)+1,step=1;i<Npoints ;++i){
		if( (*j)+1 > Npoints-1 || i==(*j)){
			printf("ERROR: i=%li (*j)=%li\n",i,(*j));
			exit(1);
		}

		if(step){
			delta = sqrt( pow(points[(*j)].x[0]-points[i].x[0],2)
					+ pow(points[(*j)].x[1]-points[i].x[1],2) )
					< 1.01*(points[(*j)].gridsize + points[i].gridsize)/2;
		}else{
			delta = AreBoxNeighbors(&(points[i]),&(points[(*j)]) );
		}

		if(delta){
			if(i==(*end)) *end=(*j)+1;
			for(k=i;k>(*j)+1;--k) SwapPointsInArray( &(points[k]) , &(points[k-1]) );
			//SwapPointsInArray( &(points[(*j)+1]) , &(points[i]) );
			++(*j);
			i=(*j);
			step=1;
		}

		// try larger linking length
		if( (i == Npoints-1) && step){
			i=(*j);
			step=0;
		}
	}
}

short backtrack(Point *points,long Npoints,long *j,long jold,long *end){
	  /* work backward along curve to find point with another neighbor
	   * further along in the array
	   * end - the index of the point where end has been moved to
	   *       if the end point of the path is already known points > end
	   *       are spurs
	   *  Does not change order of any point < j
	   *  will not go back past jold
	   */

	long m,i,k;

	  for(m=(*j)-1;m>jold;--m){
		  for(i=(*j)+1;i<Npoints;++i){

			  if(AreBoxNeighbors(&(points[i]),&(points[m])) ){
				  // add neighbor to end of curve and restart
				  if(i==*end) *end=*j+1;
				  // move point to position after last one
				  for(k=i;k>(*j)+1;--k) SwapPointsInArray( &(points[k]),&(points[k-1]) );
				//SwapPointsInArray(&(points[(*j)+1]),&(points[i]) );
				  ++(*j);
				  return 1;
	    	  }
		  }
	  }

	  return 0;
}

void nesting_curve(ImageInfo *curves,int Ncurves){
	/* determine which curves are enclosed in another curve
	 */
	int i,k;

	if(Ncurves==1){
		curves->Nencircled=0;
		return;
	}

	for(k=0;k<Ncurves;++k){
		// calculate how many curves enclose k-th curve
		curves[k].Nencircled=0;
		for(i=0;i<Ncurves;++i){
			if( (k != i) &&
				(abs(windings(curves[k].points[0].x,curves[i].points,curves[i].Npoints,&(curves[i].area),0)) > 0) ){
				++curves[k].Nencircled;
			}
		}
	}

}

void split_images(TreeHndl i_tree,ImageInfo *images,int Maximages
		,int *Nimages,Boolean sortallpoints){
	/*
	 * splits images by edge method to avoid FOF
	 *    start: all image points should be under image[0].points
	 *
	 *    Note: points in tree must be marked as inside or outside image previously
	 *    this can be changed by changing findborder3 to findborder2
	 */

	if(images->Npoints == 0){
		*Nimages=0;
		return ;
	}
	unsigned long Npoints,Npoints_tot,i,j,k,jold;
	ImageInfo *borders;
	int TmpNimages=0,maxN,m;
	double *image_number_array,tmp;

	Npoints_tot=images->Npoints;

	//printf("entering split_images\n");
	//checkTree(i_tree);

	// find borders of total image
	//findborders2(i_tree,images);
	findborders3(i_tree,images);

	// copy outer borders to borders array
	borders=(ImageInfo *) malloc(Maximages*sizeof(ImageInfo));
	assert(borders);
	Npoints=images->outerborder->Nunits;

	// make a copy of the border points in a point array
	borders->points=NewPointArray(Npoints,False);
	// copy outer borders of images to border
	MoveToTopKist(images->outerborder);
	for(i=0;i<images->outerborder->Nunits;++i){
		PointCopyData(&(borders[0].points[i]),getCurrentKist(images->outerborder));
		MoveDownKist(images->outerborder);
	}
	borders->Npoints = images->outerborder->Nunits;

	printf("borders  %i  images %i  outerborder %i\n",borders->Npoints,images->Npoints
			,images->outerborder->Nunits);

	split_order_curve(borders,Maximages,&TmpNimages);

//	splitter(borders,Maximages,&TmpNimages);

	printf("*Nimages=%i\n",TmpNimages);
	// split up into separable curves
	// order points in each curve
//	for(i=0;i<TmpNimages;++i){ order_curve(&(borders[i])); printf("%i\n",borders[i].Npoints);}

	//printf("number of borders = %i\n",TmpNimages);

	// classify curves into those that are inside another curve or not
	nesting_curve(borders,TmpNimages);
	// if border is encircled by an even number of other borders count it as an
	//      outer border of an image
	for(i=0,*Nimages=0;i<TmpNimages;++i) if(borders[i].Nencircled % 2 == 0){
		images[*Nimages].area=borders[i].area;  // only an estimate of the area
		++*Nimages;
	}


	// divide points into separate images
	if(sortallpoints){
		if(*Nimages > 1){
			image_number_array=(double *)malloc(Npoints_tot*sizeof(double));
			assert(image_number_array);
			for(j=0,maxN=0;j<TmpNimages;++j) if(borders[j].Nencircled % 2 == 0 && borders[j].Nencircled > maxN) maxN=borders[j].Nencircled;
			// sort points into images
			//#pragma omp parallel for private(j)

			// mark each point with which image it is in
			for(k=0;k<Npoints_tot;++k){
				image_number_array[k]=-1;
				for(j=0;j<TmpNimages;++j){
					if(borders[j].Nencircled % 2 == 0){  // j is a outer image border

						if( abs(windings(images->points[k].x,borders[j].points,borders[j].Npoints
								,&(images[i].area),0) ) == 1){
							// point is in j
							if(image_number_array[k] == -1 ||
									borders[(int)(image_number_array[k]+0.5)].Nencircled < borders[j].Nencircled){
								image_number_array[k]=(double) j;
							// if point can't be in another move on to next image point
								if(borders[j].Nencircled == maxN) j=TmpNimages;
							}
						}
					}
				}
				// check that an image is found for every image point
				assert(image_number_array[k] > -1);
			}

		//printf("Npoints=%i\n",Npoints);
			if(*Nimages > 2){
				// sort points into images
				double_sort_points(Npoints_tot,image_number_array-1,images->points);

				for(k=0,jold=0,i=0;k<TmpNimages;++k){
					//for(i=1,jold=0;i<*Nimages;++i){
					if(borders[k].Nencircled % 2 == 0){
						if(i>0){
							locateD(image_number_array-1,Npoints_tot,k-0.5,&j);
							images[i].points=&(images[0].points[j]);
							//printf("%e %e %i %i\n",images[i].points[0].x[0],images[i].points[0].x[1],j,k);
							images[i-1].Npoints=j-jold;
							jold=j;
						}
						++i;
					}
				}
				images[*Nimages-1].Npoints=Npoints_tot-jold;

			}else{
				// when there are two images there is a faster way of sorting
				for(i=0;i<TmpNimages;++i) if(borders[i].Nencircled % 2 == 0) break;
				//printf("i=%i\n",i);
				for(k=0,j=images[0].Npoints-1,images[1].Npoints=0;k<j;++k){
					if(image_number_array[k] > i+0.5){
					//printf("image_number_array[%i]=%f\n",k,image_number_array[k]);
						SwapPointsInArray(&(images[0].points[j]),&(images[0].points[k]));
						tmp=image_number_array[k];
						image_number_array[k]=image_number_array[j];
						image_number_array[j]=tmp;
						images[1].points=&(images[0].points[j]);
						++images[1].Npoints;
						--images[0].Npoints;
						--j;
						--k;
					}
				}
			}


			free(image_number_array);

			// find borders of individual images, and take out images
			//     that of no points within them
			for(i=0;i<*Nimages;++i){
				if(images[i].Npoints==0){
					for(j=i;j<*Nimages-1;++j){
						images[j].Nencircled=images[j+1].Nencircled;
						images[j].Npoints=images[j+1].Npoints;
						images[j].points=images[j+1].points;
						images[j].area=images[j+1].area;
						images[j].area_error=images[j+1].area_error;
						for(m=0;m<3;++m) images[j].gridrange[m]=images[j+1].gridrange[m];
					}
					--(*Nimages);
				}
				//findborders2(i_tree,&(images[i]));
				printf("images[%i].Npoints=%i\n",i,images[i].Npoints);
				findborders3(i_tree,&(images[i]));
			}

		}
	}


	free(borders->points);
	free(borders);

	//checkTree(i_tree);
	//printf("exiting split_images\n");
}

void split_images3(TreeHndl i_tree,ImageInfo *images,int Maximages
		,int *Nimages,Boolean sortallpoints){
	/*
	 * splits images by edge method to avoid FOF
	 *    start: all image points should be under image[0].points
	 *
	 *    Note: points in tree must be marked as inside or outside image previously
	 *    this can be changed by changing findborder3 to findborder2
	 *
	 *    Warning: works only with simply connect images
	 */

	if(images->Npoints == 0){
		*Nimages=0;
		return ;
	}
	unsigned long Npoints,Npoints_tot,i,j,k,jold;
	ImageInfo *borders;
	int TmpNimages=0,maxN,m;
	double *image_number_array,tmp,r,rmin;

	Npoints_tot=images->Npoints;

	// find borders of total image
	findborders3(i_tree,images);

	// copy inner borders to borders array
	borders=(ImageInfo *) malloc(Maximages*sizeof(ImageInfo));
	assert(borders);
	Npoints=images->innerborder->Nunits;

	// make a copy of the border points in a point array
	borders->points=NewPointArray(Npoints,False);
	// copy inner borders of images to border
	MoveToTopKist(images->innerborder);
	for(i=0;i<images->innerborder->Nunits;++i){
		PointCopyData(&(borders[0].points[i]),getCurrentKist(images->innerborder));
		MoveDownKist(images->innerborder);
	}
	borders->Npoints = images->innerborder->Nunits;

	//printf("borders  %li  images %li  innerborder %li\n",borders->Npoints,images->Npoints
	//		,images->innerborderkist->Nunits);
	// split up into separable curves
	splitter(borders,Maximages,Nimages);

	TmpNimages=*Nimages;
	//printf("*Nimages=%i\n",TmpNimages);

	// order points in each curve
	/*for(i=0;i<TmpNimages;++i){
		order_curve(&(borders[i]));
		printf("%i\n",borders[i].Npoints);
	}*/
	//printf("number of borders = %i\n",TmpNimages);
	// classify curves into those that are inside another curve or not
	//nesting_curve(borders,TmpNimages);
	// if border is encircled by an even number of other borders count it as an
	//      outer border of an image
	//for(i=0,*Nimages=0;i<TmpNimages;++i) if(borders[i].Nencircled % 2 == 0){
	//	images[*Nimages].area=borders[i].area;  // only an estimate of the area
	//	++*Nimages;
	//}

	// divide points into separate images
	if(sortallpoints){
		if(*Nimages > 1){
			image_number_array=(double *)malloc(Npoints_tot*sizeof(double));
			assert(image_number_array);
			// sort points into images
			//#pragma omp parallel for private(j)

			// mark each point with which image it is in
			for(k=0;k<Npoints_tot;++k){
				image_number_array[k]=-1;
				for(j=0,rmin=1.0e99;j<TmpNimages;++j){
					for(i=0;i<borders[j].Npoints;++i){
						r = pow(borders[j].points[i].x[0] - images->points[k].x[0],2)
						  + pow(borders[j].points[i].x[1] - images->points[k].x[1],2);
						if(r<rmin){ rmin=r; image_number_array[k]=j;}
					}
				}
			}

		//printf("Npoints=%i\n",Npoints);
			if(*Nimages > 2){
				// sort points into images
				//double_sort_points(Npoints_tot,image_number_array-1,images->points);
				quicksortPoints(images->points,image_number_array,Npoints_tot);

				for(k=1,jold=0;k<TmpNimages;++k){
					//for(i=1,jold=0;i<*Nimages;++i){
					//if(borders[k].Nencircled % 2 == 0){
					locateD(image_number_array-1,Npoints_tot,k-0.5,&j);
					images[k].points=&(images[0].points[j]);
					//printf("%e %e %i %i\n",images[i].points[0].x[0],images[i].points[0].x[1],j,k);
					images[k-1].Npoints=j-jold;
					jold=j;
				}

				images[*Nimages-1].Npoints=Npoints_tot-jold;

			}else{
				// when there are two images there is a faster way of sorting
				//for(i=0;i<TmpNimages;++i) if(borders[i].Nencircled % 2 == 0) break;
				//printf("i=%i\n",i);
				for(k=0,j=images[0].Npoints-1,images[1].Npoints=0;k<j;++k){
					if(image_number_array[k] > i+0.5){
					//printf("image_number_array[%i]=%f\n",k,image_number_array[k]);
						SwapPointsInArray(&(images[0].points[j]),&(images[0].points[k]));
						tmp=image_number_array[k];
						image_number_array[k]=image_number_array[j];
						image_number_array[j]=tmp;
						images[1].points=&(images[0].points[j]);
						++images[1].Npoints;
						--images[0].Npoints;
						--j;
						--k;
					}
				}
			}


			free(image_number_array);

			// find borders of individual images, and take out images
			//     that of no points within them
			for(i=0;i<*Nimages;++i){
				if(images[i].Npoints==0){
					for(j=i;j<*Nimages-1;++j){
						images[j].Nencircled=images[j+1].Nencircled;
						images[j].Npoints=images[j+1].Npoints;
						images[j].points=images[j+1].points;
						images[j].area=images[j+1].area;
						images[j].area_error=images[j+1].area_error;
						for(m=0;m<3;++m) images[j].gridrange[m]=images[j+1].gridrange[m];
					}
					--(*Nimages);
				}
				//findborders2(i_tree,&(images[i]));
				//printf("images[%i].Npoints=%i\n",i,images[i].Npoints);
				findborders3(i_tree,&(images[i]));
			}

		}
	}

	free(borders->points);
	free(borders);

	//checkTree(i_tree);
	//printf("exiting split_images\n");
}


void split_images2(TreeHndl i_tree,ImageInfo *images,int Maximages,int *Nimages){
	int i;
	/* brute force method of splitting images by doing neighbors-of-neighbors
	 * on all the points.  Dependable, but slow.
	 *
	 * in_image markers must be set
	 */

	if(images[0].Npoints==0){
		*Nimages=0;
		return;
	}

	//printf("In split_image2\n");
	splitter(images,Maximages,Nimages);
	//printf("Images split\n");
	 //find borders of each image
	for(i=0;i<*Nimages;++i) findborders3(i_tree,&(images[i]));
	//printf("Out split_image2\n");

	return ;
}

void splitter(ImageInfo *images,int Maximages,int *Nimages){
	/* meant to be a sure fire way to split all points into separate images or
	 *   curves into separate curves
	 */
	long i,m,j;
	ListHndl imagelist=NewList();
	unsigned long NpointsTotal=0;
	Point *point,*newpointarray;

	//assert(images[0].Npoints);

	if(images[0].Npoints==0){
		*Nimages=0;
		return;
	}

	NpointsTotal = images[0].Npoints;

	// copy image points into a list
	for(i=0;i<images[0].Npoints;++i) InsertPointAfterCurrent(imagelist,&(images[0].points[i]));

	assert(imagelist->Npoints == images[0].Npoints);
	//printf("imagelist = %il\n",imagelist->Npoints);

	splitlist(imagelist,images,Nimages,Maximages);
	//printf("imagelist = %il NpointsTotal = %il\n",imagelist->Npoints,NpointsTotal);
	assert(imagelist->Npoints == NpointsTotal);

	// copy list back into array
	point = images[0].points;
	newpointarray = NewPointArray(NpointsTotal,False);
	MoveToTopList(imagelist);
	m=0;
	i=0;
	do{
		if(i < *Nimages && images[i].points==imagelist->current){
			images[i].points=&(newpointarray[m]);
			++i;
		}
		PointCopyData(&(newpointarray[m]),imagelist->current);
		++m;
	}while(MoveDownList(imagelist));

	assert(m == NpointsTotal);

	free(imagelist);
	free(point);

	// take out images that have no points in them
	for(i=0;i<*Nimages;++i){
		if(images[i].Npoints==0){
			for(j=i;j<*Nimages-1;++j){
				images[j].Nencircled=images[j+1].Nencircled;
				images[j].Npoints=images[j+1].Npoints;
				images[j].points=images[j+1].points;
				images[j].area=images[j+1].area;
				images[j].area_error=images[j+1].area_error;
				for(m=0;m<3;++m) images[j].gridrange[m]=images[j+1].gridrange[m];
			}
			--(*Nimages);
			--i;
		}
	}

	return ;
}

void splitlist(ListHndl imagelist,ImageInfo *images,int *Nimages,int Maximages){
	/* reorders imagelist into separate images using reliable
	 *      neighbor-of-neighbor method
	 *
	 * in images the number of points in each image is updated
	 *      and the pointer to the first point in each image
	 *
	 */
	unsigned long i=0,m=0;
	ListHndl orderedlist = NewList();

	//printf("imagelist = %li\n",imagelist->Npoints);
	// divide images into disconnected curves using neighbors-of-neighbors

	while(imagelist->Npoints > 0 && i < Maximages){
		images[i].points = TakeOutCurrent(imagelist);
		MoveToBottomList(orderedlist);
		InsertPointAfterCurrent(orderedlist,images[i].points);
		MoveDownList(orderedlist);
		NeighborsOfNeighbors(orderedlist,imagelist);
		images[i].Npoints = orderedlist->Npoints - m;
		m += images[i].Npoints;
		++i;
	}
	*Nimages = i;

	// if there are too many images put the extra points in one last image
	if(i == Maximages && imagelist->Npoints > 0){
		Point *point;
		do{
			point = TakeOutCurrent(imagelist);
			MoveToBottomList(orderedlist);
			InsertPointAfterCurrent(orderedlist,point);
			MoveDownList(orderedlist);
			++(images[Maximages-1].Npoints);
		}while(imagelist->Npoints > 0);
	}

	if(i >= Maximages){
		ERROR_MESSAGE();
		printf("Too many images for Note enough images\n");
	}

	assert(imagelist->Npoints == 0);

	imagelist->Npoints = orderedlist->Npoints;
	imagelist->top = orderedlist->top;
	imagelist->bottom = orderedlist->bottom;
	imagelist->current = orderedlist->current;

	//EmptyList(imagelist);
	//free(imagelist);
	//imagelist = orderedlist;

	free(orderedlist);
	//printf("imagelist = %il\n",imagelist->Npoints);
	return;
}
