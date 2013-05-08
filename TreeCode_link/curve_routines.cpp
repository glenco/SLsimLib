/*
 * curve_routines.c
 *
 *  Created on: Jan 15, 2010
 *      Author: R.B. Metcalf
 */

#include "slsimlib.h"

/**  orders points in a curve, separates disconnected curves
 *   curves[0...Maxcurves] must be allocated before
 *
 *	uses neighbors-of-neighbors to split into curves and then
 *	uses sorts by angle and then walks the ist sorting
 *
 *  can break down for crescent curves
 */

void split_order_curve4(OldImageInfo *curves,int Maxcurves,int *Ncurves){

	//OldImageInfo* curves = new OldImageInfo(curves_in);

	long i,m,j,end;
	//short spur,closed,attach;
	unsigned long NpointsTotal;
	double center[2],*theta;
	//bool delta,tmp,step;
	//ListHndl reservoir,orderedlist;
	//Point *newpointarray;

	//std::printf("entering split_order_curve\n");
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
		//std::printf("curves[%i].Npoints=%i %i %i\n",i,curves[i].Npoints,orderedlist->Npoints
		//		,reservoir->Npoints);
		m+=curves[i].Npoints;
		++i;
	}
	*Ncurves=i;

	// copy list back into array
	point=curves[0].points;
	newpointarray=NewPointArray(NpointsTotal,false);
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
			Utilities::quicksortPoints(curves[i].points,theta,curves[i].Npoints);

			// check to make sure the center is inside the curves
			//assert(abs(windings(center,curves[i].points,curves[i].Npoints,&tmp1,0)));

			//std::printf("N=%i\n",curves[i].Npoints);
			//assert(abs(windings(center,curves[i].points,curves[i].Npoints,&tmp1,0)) > 0 );

			// find the last point that is a neighbor of the first point
			m=curves[i].Npoints-1;
			while(!AreBoxNeighbors(&(curves[i].points[0]),&(curves[i].points[m]))) --m;
			end=m;
			//assert(m > 0);

		//std::printf("windings = %i\n",windings(center,curves[i].points,curves[i].Npoints,&tmp1,0));
		// walk curve to remove shadowing effect
			j=0;
		//end=0;
			m=0;
			while(j < curves[i].Npoints-1){
				walkcurve(curves[i].points,curves[i].Npoints,&j,&end);
				//std::printf("i = %i Npoints = %i end+1 = %i j = %i\n",i,curves[i].Npoints,end+1,j);
				if(j < curves[i].Npoints-1)
					backtrack(curves[i].points,curves[i].Npoints,&j,-1,&end);
				++m;
				assert(m < curves[i].Npoints);
			}
			//std::printf("i = %i Npoints = %i end+1 = %i j=%i\n",i,curves[i].Npoints,end+1,j);
			curves[i].Npoints=end+1;
		}
	}

	free(theta);
	//free(reservoir);
	//free(orderedlist);
	//free(point);

	return ;
}
namespace Utilities{
/** \ingroup Utill
 *
 * \brief Orders points on a closed curve.
 *
 * The algorithm first finds the "center" of the curve. It then does a rough ordering according to the angle
 * around this center. It then walks the curve jumping to a neighbor cell each step choosing a neighbor along
 * one of the x or y-axis before taking a diagonal step.  If it comes to the point where there is no more neighbors
 * (as may occur after going through a self-intersection and then returning to it) the algorithm backtracks until
 * it finds a point in the ordered list that is also a neighbor to a point in the not yet ordered list and attaches
 * this to the end of the ordered list and continuous to walk.  This algorithm works well at finding a closed loop.  It
 * can cut off points from the curve that are either in loops or if four cells intersect and are all on the curve as
 * can happen when there is a lot of structure in the curve that is not resolved at the gridsize used.  The points that are
 * cut of are at the end of the array in no guaranteed order.
 *
 * Returns true if a closed loop is found which incorporates all of the points.
 *
 * This algorithm could be improve by inserting the remaining points, if any, into the existing curve and recursively calling itself.
 */
unsigned long order_curve4(Point *curve,long Npoints){

	long m,j,end;
	double center[2],*theta;

	//std::printf("entering split_order_curve\n");
	if(Npoints < 3) return false;


	// order curve points
	theta=(double *)malloc(Npoints*sizeof(double));
	assert(theta);


	// sort points by angle around center point
	center[0]=center[1]=0;
	for(m=0;m<Npoints;++m){
		center[0]+=curve[m].x[0];
		center[1]+=curve[m].x[1];
	}
	center[0]/=Npoints;
	center[1]/=Npoints;

	for(m=0;m<Npoints;++m){
		theta[m]=atan2(curve[m].x[1]-center[1]
		              ,curve[m].x[0]-center[0]);
	}
	Utilities::quicksortPoints(curve,theta,Npoints);

	// check to make sure the center is inside the curves

	// find the last point that is a neighbor of the first point
	m=Npoints-1;
	while(!AreBoxNeighbors(&(curve[0]),&(curve[m]))) --m;
	end=m;
	//assert(m > 0);

	//std::printf("windings = %i\n",windings(center,curve,Npoints,&tmp1,0));
	// walk curve to remove shadowing effect
	j=0;
	m=0;
	while(j < Npoints-1){
		walkcurve(curve,Npoints,&j,&end);
		//std::printf("i = %i Npoints = %i end+1 = %i j = %i\n",i,Npoints,end+1,j);
		if(j < Npoints-1)
			backtrack(curve,Npoints,&j,-1,&end);
		++m;
        // TODO This is a kluge.  There should be a better exit stratagy.
        if(m >= Npoints) return false;
		assert(m < Npoints);
	}

	free(theta);

	return end+1;
}
/// Overloads and is dependent on version that takes a point array
bool order_curve4(Kist<Point> * curve){
	Point *tmpcurve = NewPointArray(curve->Nunits(),false);
	unsigned long i=0,Npoints = curve->Nunits();
	bool bo;

	curve->MoveToTop();
	do{
		PointCopyData(&tmpcurve[i++],curve->getCurrent());
	}while(curve->Down());

	bo = order_curve4(tmpcurve,curve->Nunits());

	// resort points in imagekist
	for(i=0;i<Npoints;++i){
		curve->MoveToTop();
		do{
			if(tmpcurve[i].id == curve->getCurrent()->id)
				  curve->MoveCurrentToBottom();
		}while(curve->Down());
	}

	FreePointArray(tmpcurve,false);

	return bo;
}
}
/**
 *
 *
 * I can think of one pathological case in which this routine would fail
 */
/*bool order_ExteriorBoundary(
		Point *curve           /// Array of points representing the curve
		,long Npoints          /// Number of points in curve
		,long *NewNpoints      /// Number of points in the exterior boundary, *NewNpoints <= Npoints
		,double *area          /// Area within exterior boundary
		){

	long m,j,end,k;
	double center[2],*theta,tmp;

	cout << AreBoxNeighbors(&(curve[0]),&(curve[Npoints-1])) << endl;

	//std::printf("entering split_order_curve\n");
	if(Npoints < 3) return false;


	// order curve points
	theta=(double *)malloc(Npoints*sizeof(double));
	assert(theta);


	// sort points by angle around center point
	center[0]=center[1]=0;
	for(m=0;m<Npoints;++m){
		center[0]+=curve[m].x[0];
		center[1]+=curve[m].x[1];
	}
	center[0]/=Npoints;
	center[1]/=Npoints;

	// check to make sure the center is inside the curves
	assert(abs(windings(center,curve,Npoints,&tmp,0)));

//	for(m=0;m<Npoints;++m){
//		theta[m]=atan2(curve[m].x[1]-center[1]
//		              ,curve[m].x[0]-center[0]);
//	}
//	quicksortPoints(curve,theta,Npoints);
	free(theta);

	cout << AreBoxNeighbors(&(curve[0]),&(curve[Npoints-1])) << endl;

	// make bottom most point the first point
	double ymin = curve[0].x[1];
	long imin = 0;
	for(m=1;m<Npoints;++m){
		if(ymin > curve[m].x[1]){
			ymin = curve[m].x[1];
			imin = m;
		}
	}
	// Cyclic permutation of points until bottom point is first
	for(m=0;m<imin;++m){
		for(k=0;k<Npoints-1;++k) SwapPointsInArrayData( &(curve[k]) , &(curve[k+1]) );
	}

	cout << AreBoxNeighbors(&(curve[0]),&(curve[Npoints-1])) << endl;
	cout << AreBoxNeighbors(&(curve[Npoints-imin]),&(curve[Npoints-imin-1])) << endl;

	// find the last point that is a neighbor of the first point
	j=0;
	do{
		m = Npoints-1;
		while(!AreBoxNeighbors(&(curve[j]),&(curve[m]))) --m;
		end = m;
		++j;
	}while(end == 0);

	// walk curve to remove shadowing effect
	j=0;
	m=0;
	while(j < end){
		walkcurveRight(curve,Npoints,&j,&end);
		//std::printf("i = %i Npoints = %i end+1 = %i j = %i\n",i,Npoints,end+1,j);
		if(j < Npoints-1)
			backtrack(curve,Npoints,&j,0,&end);
		++m;
		if(m > Npoints-1){
			*NewNpoints = 0;
			return false;  // did not succeed
		}
	}
	//std::printf("i = %i Npoints = %i end+1 = %i j=%i\n",i,Npoints,end+1,j);
	*NewNpoints = end + 1;

	windings(curve[0].x,curve,*NewNpoints,area,0);

	return true;
}*/
/* \ingroup Utill
 *
 * \brief Finds area within a curve by summing every cell.
 *
 * Should work every time provided the curve ordering is correct.
 * Testing if each cell is inside the curve can be slow.
 */
/*double findAreaOfCurve(TreeHndl tree,ImageInfo *curve,int NimageMax){

	if(curve->imagekist->Nunits() < 3) return 0.0;

	int Nimages,i,imax=0;
	double xcm[2],area,tmp;
	ImageInfo *borders = new ImageInfo[NimageMax];

	// find borders of curve
	findborders4(tree,curve);
	// copy outer borders into borders->imagekist
	curve->outerborder->MoveToTop();
	do{
		borders->imagekist->InsertAfterCurrent(curve->outerborder->getCurrent());
	}while(curve->outerborder->Down());

	borders->imagekist->Print();
exit(0);
	// split borders up
	divide_images_kist(tree,borders,&Nimages,NimageMax);



	area = 0.0;
	xcm[0]=xcm[1]=0.0;
	for(i=0;i<Nimages;++i){
		// order the curve
		order_curve4(borders[i].imagekist);
		windings(xcm,borders[i].imagekist,&tmp,0);
		if(tmp < area){
			tmp = area;
			imax=i;
		}
	}

	borders[imax].imagekist->Print();

	exit(0);

	delete[] borders;

	return area;
}
	/*******
	// find "center" of curve
	xcm[0]=xcm[1]=0.0;
	curve->MoveToTop();
	do{
		xcm[0] += curve->getCurrent()->x[0];
		xcm[1] += curve->getCurrent()->x[1];
	}while(curve->Down());
	xcm[0] /= curve->Nunits();
	xcm[1] /= curve->Nunits();

	//cout << "windings " << windings(xcm,curve,&tmp,0) << endl;

	// find farthest point on curve from center
	rmax = 0.0;
	curve->MoveToTop();
	do{
		r = pow(xcm[0] - curve->getCurrent()->x[0],2)
				+ pow(xcm[1] - curve->getCurrent()->x[1],2);
		if(r>rmax){
			rmax=r;
		}
	}while(curve->Down());

	PointsWithinKist(tree,xcm,sqrt(rmax),subkist,0);

	area = 0.0;
	subkist->MoveToTop();
	do{
		if(windings(subkist->getCurrent()->x,curve,&tmp,0)){
			area += (subkist->getCurrent()->leaf->boundary_p2[0] - subkist->getCurrent()->leaf->boundary_p1[0])
				*(subkist->getCurrent()->leaf->boundary_p2[1] - subkist->getCurrent()->leaf->boundary_p1[1]);

		}else{
			subkist->TakeOutCurrent();
		}
	}while(subkist->Down());

	subkist->MoveToTop();
	cout << subkist->Nunits() << endl;
	do{
		cout << subkist->getCurrent()->leaf->boundary_p1[0]
		        << "  " << subkist->getCurrent()->leaf->boundary_p1[1]
				<< "   " << subkist->getCurrent()->leaf->boundary_p2[0]
				<< "  " << subkist->getCurrent()->leaf->boundary_p2[1]
				                          << endl;
	}while(subkist->Down());
	assert(area >= 0.0);
	delete subkist;

	return area;
}
*/

/*  orders points in a curve, separates disconnected curves
 *   curves[0...Maxcurves] must be allocated before
 *
 *	 uses list instead of arrays
 *   does not attempt to remove spurs
 */
/*void split_order_curve3(OldImageInfo *curves,int Maxcurves,int *Ncurves){

	long i,k=0,m;
	short spur,closed,attach;
	unsigned long Npoints;
	bool delta,tmp,step;
	ListHndl lists[Maxcurves+1];
	Point *point,*newpointarray;

	lists[Maxcurves]=NewList();
	for(i=0;i<curves[0].Npoints;++i) InsertPointAfterCurrent(lists[Maxcurves],&(curves[0].points[i]));

	//std::printf("entering split_order_curve\n");
	if(curves[0].Npoints==0){
		*Ncurves=0;
		return;
	}

	// separate critical curves
	Npoints=curves[0].Npoints;
	//std::printf("number of curve points: %i\n",curves[0].Npoints);

	for(k=0;lists[Maxcurves]->Npoints > 0;++k){
		std::printf("lists[Maxcurves]->Npoints=%li  \n",lists[Maxcurves]->Npoints);
		lists[k]=NewList();
	    if(k==Maxcurves){ ERROR_MESSAGE(); std::printf("ERROR: in split_order_curve3, not large enough curves array k=%li\n",k); break;}

		// start new curve
		point=TakeOutCurrent(lists[Maxcurves]);
 		InsertPointAfterCurrent(lists[k],point);
		MoveDownList(lists[k]);

	    // go through remaining points ordering them as we go down the list
	    do{
	    	spur=0;
	    	MoveToTopList(lists[Maxcurves]);
	    	tmp=true;
	    	step=true;

	    	if(lists[Maxcurves]->Npoints > 0){
	    		do{
	    			if(step){
	    				delta = sqrt( pow(lists[k]->current->x[0]-lists[Maxcurves]->current->x[0],2)
		    			        + pow(lists[k]->current->x[1]-lists[Maxcurves]->current->x[1],2) )
		    			        < 1.01*(lists[k]->current->gridsize + lists[Maxcurves]->current->gridsize)/2;
	    			}else{
	    				delta = AreBoxNeighbors(lists[k]->current,lists[Maxcurves]->current);
	    			}

	    			//std::printf("move down\n");
	    			if(delta){
	    				point=TakeOutCurrent(lists[Maxcurves]);
	    				MoveToTopList(lists[Maxcurves]);
	    				InsertPointAfterCurrent(lists[k],point);
	    				MoveDownList(lists[k]);
	    				step=true;
	    			}else{
	    				tmp=MoveDownList(lists[Maxcurves]);
	    			}

	    			if( tmp==false && step){
	    				MoveToTopList(lists[Maxcurves]);
	    				tmp=true;
	    				step=false;
	    			}

	    		}while(tmp && lists[Maxcurves]->Npoints > 0);
	    	}

	    	// check if curve is closed
	    	if( lists[k]->Npoints < 3 || AreBoxNeighbors(lists[k]->top,lists[k]->bottom ) == false ){
			  // curve is not closed

	    		closed=0;

	    		// work backward along curve to find point with another neighbor
	    		if(lists[Maxcurves]->Npoints > 0){
	    			for(m=-1;(m>-lists[k]->Npoints && spur==0);--m){

	    				MoveUpList(lists[k]);
	    				// reverse the order of the spur
	    				// this insures that if the code returns to the start
	    				// the points will still be in a order
	    				//JumpDownList(lists[k],m);
	    				//MoveCurrentToBottom(lists[k]);
	    				//MoveToBottomList(lists[k]);

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

	//std::printf("  end of loop j=%i k=%i i=%i Npoints=%i\n",j,k,i,Npoints);

	point=curves[0].points;
	newpointarray=NewPointArray(curves[0].Npoints,false);
	for(i=0,*Ncurves=0,m=0;i<k;++i){
		std::printf("lists[%li]->Npoints = %li\n",i,lists[i]->Npoints);
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
}/*

/*  orders points in a curve, separates disconnected curves
 *   curves[0...Maxcurves] must be allocated before
 *
 *   does not attempt to remove spurs
 */
/*void split_order_curve(OldImageInfo *curves,int Maxcurves,int *Ncurves){

	long j,k=0,jold,end=0;
	short spur,closed;
	unsigned long Npoints,Maxpoint;

	//std::printf("entering split_order_curve\n");
	if(curves[0].Npoints==0){
		*Ncurves=0;
		return;
	}

	// separate critical curves
	Npoints=curves[0].Npoints;
	//std::printf("number of curve points: %i\n",curves[0].Npoints);

	for(j=0,k=0,jold=0;j<Npoints;++k){
	    if(k==Maxcurves){ ERROR_MESSAGE();std::printf("ERROR: in split_order_curve, not large enough curves array k=%li j=%li\n",k,j); break;}
	    Maxpoint=Npoints;

	    // go through remaining points ordering them as we go down the list
	    do{
	    	spur=0;
	    	walkcurve(curves[0].points,Maxpoint,&j,&end);

	      // check if curve is closed

		  if( AreBoxNeighbors(&(curves[k].points[0]),&(curves[0].points[j]) ) == false &&
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


	//std::printf("  end of loop j=%i k=%i i=%i Npoints=%i\n",j,k,i,Npoints);

	*Ncurves=k;
	//std::printf("exiting split_order_curve\n");

}*/

/*  orders points in a curve, separates disconnected curves, cuts off spurs
 *   curves[0...Maxcurves] must be allocated before
 *
 *   Attempts to remove spurs,  does not always work
 */
/*void split_order_curve2(OldImageInfo *curves,int Maxcurves,int *Ncurves){

	long i,j,k=0,jold,m,end;
	short spur,closed;
	unsigned long Npoints,Maxpoint;

	//std::printf("entering split_order_curve\n");
	if(curves[0].Npoints==0){
		*Ncurves=0;
		return;
	}

	// separate critical curves
	Npoints=curves[0].Npoints;
	//std::printf("number of curve points: %i\n",curves[0].Npoints);

	for(j=0,k=0,jold=0;j<Npoints;++k){
	    if(k==Maxcurves){ ERROR_MESSAGE(); std::printf("ERROR: in split_order_curve2, not large enough curves array k=%li j=%li\n",k,j); break;}
	    Maxpoint=Npoints;

	    // go through remaining points ordering them as we go down the list
	    do{

	    	spur=0;
	    	walkcurve(curves[0].points,Maxpoint,&j,&end);


//	    	for(i=j+1,step=1,rmin=1.0e99;i<Maxpoint;++i){
//	    		if( j+1 > Npoints-1 || i==j){ ERROR_MESSAGE(); std::printf("ERROR: i=%l j=%i\n",i,j); exit(1);}
//
//	    		if(step){
//	    			delta = sqrt( pow(curves[0].points[j].x[0]-curves[0].points[i].x[0],2)
//		    			  + pow(curves[0].points[j].x[1]-curves[0].points[i].x[1],2) )
//		    			  < 1.01*(curves[0].points[j].gridsize + curves[0].points[i].gridsize)/2;
//	    		}else{
//	    			delta = AreBoxNeighbors(&(curves[0].points[i]),&(curves[0].points[j]) );
//	    		}
//
//	    		if(delta){
//	    			SwapPointsInArray( &(curves[0].points[j+1]) , &(curves[0].points[i]) );
//	    			++j;
//	    			i=j;
//	    			step=1;
//	    		}
//
//	    		// try larger linking length
//	    		if( (i == Maxpoint-1) && step){
//	    			i=j;
//	    			step=0;
//	    		}
//	    	}


	      // check if curve is closed

		  if( AreBoxNeighbors(&(curves[k].points[0]),&(curves[0].points[j]) ) == false &&
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
		    	//std::printf("%i %e %e\n",i,curves[0].points[i].x[2],curves[0].points[i].x[1]);
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
	    // std::printf("2 Npoints=%i\n",curves[k].Npoints);
	  }
	//std::printf("  end of loop j=%i k=%i i=%i Npoints=%i\n",j,k,i,Npoints);

	*Ncurves=k;
	//std::printf("exiting split_order_curve\n");

}*/

/*  orders points in a curve,
 *
 *   must already be a connected set of points
 *
 *   does not attempt to remove spurs
 */
/*void order_curve(OldImageInfo *curve){

	long j,jold,end=0;
	short spur,closed;

	//std::printf("entering split_order_curve\n");
	if(curve->Npoints==0) return;

	// separate critical curves
	//std::printf("number of curve points: %i\n",curves[0].Npoints);

	j=0; jold=0;

	// go through remaining points ordering them as we go down the list
	do{
		spur=0;
		walkcurve(curve->points,curve->Npoints,&j,&end);

		// check if curve is closed
		if( AreBoxNeighbors(&(curve->points[0]),&(curve->points[j]) ) == false &&
				(j != 0 && j < curve->Npoints) ){
				  // curve is not closed

			closed=0;
			// work backward along curve to find point with another neighbor
			spur=backtrack(curve->points,curve->Npoints,&j,0,&end);

		}else{closed=1;}

	}while(spur);

	//std::printf("  end of loop j=%i k=%i i=%i Npoints=%i\n",j,k,i,Npoints);
	//std::printf("exiting split_order_curve\n");

}*/

/**
 * orders curve by finding closest point ahead in list
 *   tries x/y jump before diagonal jump
 *   exits when it cannot find a cell neighbor ahead in array
 *   leaves j as the last point in curves
 *
 *   If the *end point is found to be in the curve it is moved with its point.  The
 *   walking continues beyond end in which case *j > *end.
 *
 */
void walkcurve(Point *points,long Npoints,long *j,long *end){

	long i,k;
	short step;
	bool delta;

	//if((*j)==0) *end=Npoints-1;
	//std::printf("end = %i\n",*end);

	for(i=(*j)+1,step=1;i<Npoints ;++i){
		if( (*j)+1 > Npoints-1 || i==(*j)){
			std::printf("ERROR: i=%li (*j)=%li\n",i,(*j));
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
/*
 * orders curve by finding closest point ahead in array.
 *   Tries to make the most right hand turn as possible at
 *   each step in the sense that the angle between the rays connecting
 *   the end point to the one behind it and connecting the end point to
 *   the next point is the minimum of the endpoint's neighbors.
 *
 *   If at the beginning of the array, it will try to walk down and to
 *   the right first.
 *
 *   If the *end point is found to be in the curve it is moved with its point.  The
 *   walking continues beyond end in which case *j > *end.
 */

/*void walkcurveRight(Point *points,long Npoints,long *j,long *end){

	long i,k,i_next;
	double mintheta,x,y,phi;

	if(*j == 0) phi = pi/2;
	else phi = atan2(points[*j].x[1] - points[*j-1].x[1] , points[*j].x[0] - points[*j-1].x[0]);
	i_next = -1;

	while(*j < Npoints){
		for(i=(*j)+1;i<Npoints ;++i){

			// find neighbor with minimum forward angle with rerspect to the last two points
			if( AreBoxNeighbors(&(points[i]),&(points[(*j)])) ){
				x = (points[i].x[0] - points[*j].x[0])*cos(phi) + (points[i].x[1] - points[*j].x[1])*sin(phi);
				y = (points[i].x[0] - points[*j].x[0])*sin(phi) - (points[i].x[1] - points[*j].x[1])*cos(phi);
				if(mintheta > atan2(y,x)){
					i_next = i;
					mintheta = atan2(y,x);
				}
			}
		}
		if(i_next != -1){  // did find the next point
			if(i_next == (*end)) *end=(*j)+1;
			for(k=i_next;k>(*j)+1;--k) SwapPointsInArray( &(points[k]) , &(points[k-1]) );
			++(*j);
			mintheta = 2*pi;
			phi = atan2(points[*j].x[1] - points[*j-1].x[1] , points[*j].x[0] - points[*j-1].x[0]);
		}else{
			break;  // did not find next point
		}
		i_next = -1;
	}
}*/

/* work backward along curve to find point with another neighbor
 * further along in the array
 * end - the index of the point where end has been moved to
 *       if the end point of the path is already known points > end
 *       are spurs
 *  Does not change order of any point < j
 *  will not go back past jold
 *
 *   If the *end point is found to be the next point in the curve it is moved
 *   with its point and *j = *end.
 */
short backtrack(Point *points,long Npoints,long *j,long jold,long *end){

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

/*void nesting_curve(OldImageInfo *curves,int Ncurves){
	// determine which curves are enclosed in another curve

	int i,k;

	if(Ncurves==1){
		curves->ShouldNotRefine=0;
		return;
	}

	for(k=0;k<Ncurves;++k){
		// calculate how many curves enclose k-th curve
		curves[k].ShouldNotRefine=0;
		for(i=0;i<Ncurves;++i){
			if( (k != i) &&
				(abs(windings(curves[k].points[0].x,curves[i].points,curves[i].Npoints,&(curves[i].area),0)) > 0) ){
				++curves[k].ShouldNotRefine;
			}
		}
	}

}*/

/*
 * splits images by edge method to avoid FOF
 *    start: all image points should be under image[0].points
 *
 *    Note: points in tree must be marked as inside or outside image previously
 *    this can be changed by changing findborder3 to findborder2
 */
/*void split_images(TreeHndl i_tree,OldImageInfo *images,int Maximages
		,int *Nimages,bool sortallpoints){

	if(images->Npoints == 0){
		*Nimages=0;
		return ;
	}
	unsigned long Npoints,Npoints_tot,i,j,k,jold;
	OldImageInfo *borders;
	int TmpNimages=0,maxN,m;
	double *image_number_array,tmp;

	Npoints_tot=images->Npoints;

	//std::printf("entering split_images\n");
	//checkTree(i_tree);

	// find borders of total image
	//findborders2(i_tree,images);
	findborders3(i_tree,images);

	// copy outer borders to borders array
	borders=(OldImageInfo *) malloc(Maximages*sizeof(OldImageInfo));
	assert(borders);
	Npoints=images->outerborder->Nunits();

	// make a copy of the border points in a point array
	borders->points=NewPointArray(Npoints,false);
	// copy outer borders of images to border
	MoveToTopKist(images->outerborder);
	for(i=0;i<images->outerborder->Nunits();++i){
		PointCopyData(&(borders[0].points[i]),getCurrentKist(images->outerborder));
		MoveDownKist(images->outerborder);
	}
	borders->Npoints = images->outerborder->Nunits();

	std::printf("borders  %i  images %i  outerborder %i\n",borders->Npoints,images->Npoints
			,images->outerborder->Nunits());

	split_order_curve(borders,Maximages,&TmpNimages);

//	splitter(borders,Maximages,&TmpNimages);

	std::printf("*Nimages=%i\n",TmpNimages);
	// split up into separable curves
	// order points in each curve
//	for(i=0;i<TmpNimages;++i){ order_curve(&(borders[i])); std::printf("%i\n",borders[i].Npoints);}

	//std::printf("number of borders = %i\n",TmpNimages);

	// classify curves into those that are inside another curve or not
	nesting_curve(borders,TmpNimages);
	// if border is encircled by an even number of other borders count it as an
	//      outer border of an image
	for(i=0,*Nimages=0;i<TmpNimages;++i) if(borders[i].ShouldNotRefine % 2 == 0){
		images[*Nimages].area=borders[i].area;  // only an estimate of the area
		++*Nimages;
	}


	// divide points into separate images
	if(sortallpoints){
		if(*Nimages > 1){
			image_number_array=(double *)malloc(Npoints_tot*sizeof(double));
			assert(image_number_array);
			for(j=0,maxN=0;j<TmpNimages;++j) if(borders[j].ShouldNotRefine % 2 == 0 && borders[j].ShouldNotRefine > maxN) maxN=borders[j].ShouldNotRefine;
			// sort points into images

			// mark each point with which image it is in
			for(k=0;k<Npoints_tot;++k){
				image_number_array[k]=-1;
				for(j=0;j<TmpNimages;++j){
					if(borders[j].ShouldNotRefine % 2 == 0){  // j is a outer image border

						if( abs(windings(images->points[k].x,borders[j].points,borders[j].Npoints
								,&(images[i].area),0) ) == 1){
							// point is in j
							if(image_number_array[k] == -1 ||
									borders[(int)(image_number_array[k]+0.5)].ShouldNotRefine < borders[j].ShouldNotRefine){
								image_number_array[k]=(double) j;
							// if point can't be in another move on to next image point
								if(borders[j].ShouldNotRefine == maxN) j=TmpNimages;
							}
						}
					}
				}
				// check that an image is found for every image point
				assert(image_number_array[k] > -1);
			}

		//std::printf("Npoints=%i\n",Npoints);
			if(*Nimages > 2){
				// sort points into images
				double_sort_points(Npoints_tot,image_number_array-1,images->points);

				for(k=0,jold=0,i=0;k<TmpNimages;++k){
					//for(i=1,jold=0;i<*Nimages;++i){
					if(borders[k].ShouldNotRefine % 2 == 0){
						if(i>0){
							locateD(image_number_array-1,Npoints_tot,k-0.5,&j);
							images[i].points=&(images[0].points[j]);
							//std::printf("%e %e %i %i\n",images[i].points[0].x[0],images[i].points[0].x[1],j,k);
							images[i-1].Npoints=j-jold;
							jold=j;
						}
						++i;
					}
				}
				images[*Nimages-1].Npoints=Npoints_tot-jold;

			}else{
				// when there are two images there is a faster way of sorting
				for(i=0;i<TmpNimages;++i) if(borders[i].ShouldNotRefine % 2 == 0) break;
				//std::printf("i=%i\n",i);
				for(k=0,j=images[0].Npoints-1,images[1].Npoints=0;k<j;++k){
					if(image_number_array[k] > i+0.5){
					//std::printf("image_number_array[%i]=%f\n",k,image_number_array[k]);
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
						images[j].ShouldNotRefine=images[j+1].ShouldNotRefine;
						images[j].Npoints=images[j+1].Npoints;
						images[j].points=images[j+1].points;
						images[j].area=images[j+1].area;
						images[j].area_error=images[j+1].area_error;
						for(m=0;m<3;++m) images[j].gridrange[m]=images[j+1].gridrange[m];
					}
					--(*Nimages);
				}
				//findborders2(i_tree,&(images[i]));
				std::printf("images[%i].Npoints=%i\n",i,images[i].Npoints);
				findborders3(i_tree,&(images[i]));
			}

		}
	}


	free(borders->points);
	free(borders);

	//checkTree(i_tree);
	//std::printf("exiting split_images\n");
}*/

/*
 * splits images by edge method to avoid FOF
 *    start: all image points should be under image[0].points
 *
 *    Note: points in tree must be marked as inside or outside image previously
 *    this can be changed by changing findborder3 to findborder2
 *
 *    Warning: works only with simply connect images
 */
/*void split_images3(TreeHndl i_tree,OldImageInfo *images,int Maximages
		,int *Nimages,bool sortallpoints){

	if(images->Npoints == 0){
		*Nimages=0;
		return ;
	}
	unsigned long Npoints,Npoints_tot,i,j,k,jold;
	OldImageInfo *borders;
	int TmpNimages=0,maxN,m;
	double *image_number_array,tmp,r,rmin;

	Npoints_tot=images->Npoints;

	// find borders of total image
	findborders3(i_tree,images);

	// copy inner borders to borders array
	borders=(OldImageInfo *) malloc(Maximages*sizeof(OldImageInfo));
	assert(borders);
	Npoints=images->innerborder->Nunits();

	// make a copy of the border points in a point array
	borders->points=NewPointArray(Npoints,false);
	// copy inner borders of images to border
	MoveToTopKist(images->innerborder);
	for(i=0;i<images->innerborder->Nunits();++i){
		PointCopyData(&(borders[0].points[i]),getCurrentKist(images->innerborder));
		MoveDownKist(images->innerborder);
	}
	borders->Npoints = images->innerborder->Nunits();

	//std::printf("borders  %li  images %li  innerborder %li\n",borders->Npoints,images->Npoints
	//		,images->innerborderkist->Nunits());
	// split up into separable curves
	splitter(borders,Maximages,Nimages);

	TmpNimages=*Nimages;
	//std::printf("*Nimages=%i\n",TmpNimages);

	// order points in each curve
//	for(i=0;i<TmpNimages;++i){
//		order_curve(&(borders[i]));
//		std::printf("%i\n",borders[i].Npoints);
//	}
	//std::printf("number of borders = %i\n",TmpNimages);
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

		//std::printf("Npoints=%i\n",Npoints);
			if(*Nimages > 2){
				// sort points into images
				//double_sort_points(Npoints_tot,image_number_array-1,images->points);
				quicksortPoints(images->points,image_number_array,Npoints_tot);

				for(k=1,jold=0;k<TmpNimages;++k){
					//for(i=1,jold=0;i<*Nimages;++i){
					//if(borders[k].Nencircled % 2 == 0){
					locateD(image_number_array-1,Npoints_tot,k-0.5,&j);
					images[k].points=&(images[0].points[j]);
					//std::printf("%e %e %i %i\n",images[i].points[0].x[0],images[i].points[0].x[1],j,k);
					images[k-1].Npoints=j-jold;
					jold=j;
				}

				images[*Nimages-1].Npoints=Npoints_tot-jold;

			}else{
				// when there are two images there is a faster way of sorting
				//for(i=0;i<TmpNimages;++i) if(borders[i].Nencircled % 2 == 0) break;
				//std::printf("i=%i\n",i);
				for(k=0,j=images[0].Npoints-1,images[1].Npoints=0;k<j;++k){
					if(image_number_array[k] > i+0.5){
					//std::printf("image_number_array[%i]=%f\n",k,image_number_array[k]);
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
						images[j].ShouldNotRefine=images[j+1].ShouldNotRefine;
						images[j].Npoints=images[j+1].Npoints;
						images[j].points=images[j+1].points;
						images[j].area=images[j+1].area;
						images[j].area_error=images[j+1].area_error;
						for(m=0;m<3;++m) images[j].gridrange[m]=images[j+1].gridrange[m];
					}
					--(*Nimages);
				}
				//findborders2(i_tree,&(images[i]));
				//std::printf("images[%i].Npoints=%i\n",i,images[i].Npoints);
				findborders3(i_tree,&(images[i]));
			}

		}
	}

	free(borders->points);
	free(borders);

	//checkTree(i_tree);
	//std::printf("exiting split_images\n");
}*/


/* brute force method of splitting images by doing neighbors-of-neighbors
 * on all the points.  Dependable, but slow.
 *
 * in_image markers must be set
 */
/*void split_images2(TreeHndl i_tree,OldImageInfo *images,int Maximages,int *Nimages){
	int i;

	if(images[0].Npoints==0){
		*Nimages=0;
		return;
	}

	//std::printf("In split_image2\n");
	splitter(images,Maximages,Nimages);
	//std::printf("Images split\n");
	 //find borders of each image
	for(i=0;i<*Nimages;++i) findborders3(i_tree,&(images[i]));
	//std::printf("Out split_image2\n");

	return ;
}*/

/** meant to be a sure fire way to split all points into separate images or
 *   curves into separate curves
 */
void splitter(OldImageInfo *images,int Maximages,int *Nimages){
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
	//std::printf("imagelist = %il\n",imagelist->Npoints);

	splitlist(imagelist,images,Nimages,Maximages);
	//std::printf("imagelist = %il NpointsTotal = %il\n",imagelist->Npoints,NpointsTotal);
	assert(imagelist->Npoints == NpointsTotal);

	// copy list back into array
	point = images[0].points;
	newpointarray = NewPointArray(NpointsTotal,false);
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
				images[j].ShouldNotRefine=images[j+1].ShouldNotRefine;
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

/** reorders imagelist into separate images using reliable
 *      neighbor-of-neighbor method
 *
 * in images the number of points in each image is updated
 *      and the pointer to the first point in each image
 *
 */
void splitlist(ListHndl imagelist,OldImageInfo *images,int *Nimages,int Maximages){
	unsigned long i=0,m=0;
	ListHndl orderedlist = NewList();

	//std::printf("imagelist = %li\n",imagelist->Npoints);
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
		std::printf("Too many images for Note enough images\n");
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
	//std::printf("imagelist = %il\n",imagelist->Npoints);
	return;
}
/** \brief The Utilities namespace contains functions for wide use in many classes that
 *  perform generic tasks.
 */
namespace Utilities{
/** \ingroup Utill
 * \brief windings(): winding number test for a point in a polygon
 * Returns: Number of times a curves winds around the point x.
 *
 * The number of times the curve loops around a point is calculated.
 *
 * The area of a self-intersecting curve will be the area of the regions encircled in
 * a clockwise direction minus the regions encircled in a counterclockwise direction - an
 * infinity symbol has zero area.
 */
int windings(
		double *x              /// Point for which the winding number is calculated
		,Point *points         /// The points on the border.  These must be ordered.
		,unsigned long Npoints /// number of points in curve
		,double *area          /// returns absolute the area within the curve with oriented border
		,short image           /// if == 1 the image of the curve is uses as the curve
		){
	int wn=0;
	unsigned long k,i;

	*area=0.0;
	if(Npoints < 3) return 0;

	if(image){
		for(i=0;i<Npoints;++i){
			k= i < Npoints-1 ? i+1 : 0;
			*area+=(points[i].image->x[0] + points[k].image->x[0])
					*(points[i].image->x[1] - points[k].image->x[1]);

			if(points[i].image->x[1] <= x[1]){
				if(points[k].image->x[1] > x[1])
					if( isLeft(points[i].image,points[k].image,x) > 0) ++wn;
			}else{
				if(points[k].image->x[1] <= x[1])
					if( isLeft(points[i].image,points[k].image,x) < 0) --wn;
			}
		}
	}else{

		for(i=0;i<Npoints;++i){
			k= i < Npoints-1 ? i+1 : 0;
			*area+=(points[i].x[0] + points[k].x[0])*(points[i].x[1] - points[k].x[1]);

			if(points[i].x[1] <= x[1]){
				if(points[k].x[1] > x[1])
					if( isLeft(&points[i],&points[k],x) > 0) ++wn;
			}else{
				if(points[k].x[1] <= x[1])
					if( isLeft(&points[i],&points[k],x) < 0) --wn;
			}
		}

	}

	*area = fabs(*area)*0.5;
	//std::printf("wn = %i\n",wn);
	//if(abs(wn) > 0) exit(0);
	return wn;
}
int windings(
		double *x              /// Point for which the winding number is calculated
		,Kist<Point> * kist         /// Kist of points on the border.  These must be ordered.
		,double *area          /// returns absolute the area within the curve with oriented border
		,short image           /// if == 1 the image of the curve is uses as the curve
		){
	int wn=0;
	unsigned long k,i;
	unsigned long Npoints = kist->Nunits();

	*area=0.0;
	if(Npoints < 3) return 0;

	Point **points = new Point*[Npoints];

	kist->MoveToTop();
	for(i=0;i<Npoints;++i,kist->Down()){
		if(image) points[i] = kist->getCurrent()->image;
		else points[i] = kist->getCurrent();
	}

	for(i=0;i<Npoints;++i){
		k = i < Npoints-1 ? i+1 : 0;
		*area += (points[i]->x[0] + points[k]->x[0])*(points[i]->x[1] - points[k]->x[1]);

		if(points[i]->x[1] <= x[1]){
			if(points[k]->x[1] > x[1])
				if( isLeft(points[i],points[k],x) > 0) ++wn;
		}else{
			if(points[k]->x[1] <= x[1])
				if( isLeft(points[i],points[k],x) < 0) --wn;
		}
	}

	*area = fabs(*area)*0.5;
	delete points;

	return wn;
}

/** \ingroup Utill
 *
 *  \brief determines the whether a point is inside a curve, that has been stretched 1.2 times
 *  returns the area of the stretched curve
 */

int windings2(
		double *x              /// Point for which the winding number is calculated
		,Point *points_original         /// The points on the border.  These must be ordered.
		,unsigned long Npoints /// number of points in curve
		,double *area          /// returns absolute the area within the curve with oriented border
		,short image           /// if == 0 the image of the curve is uses as the curve
		){
	int wn=0;
	unsigned long k,i;
	double center[2];

	center[0] = center[1] = 0.0;

	*area=0.0;
	if(Npoints < 3) return 0;

	Point **points = new Point*[Npoints];

	for(i=0;i<Npoints;++i){
		if(image) points[i] = points_original[i].image;
		else points[i] = &points_original[i];

		center[0] += points[i]->x[0];
		center[1] += points[i]->x[1];
	}

	center[0] /= Npoints;
	center[1] /= Npoints;

	for(i=0;i<Npoints;++i){
		points[i]->x[0] -= center[0];
		points[i]->x[1] -= center[1];

		points[i]->x[0] *= 1.2;
		points[i]->x[1] *= 1.2;

		points[i]->x[0] += center[0];
		points[i]->x[1] += center[1];

		k = i < Npoints-1 ? i+1 : 0;
		*area+=(points[i]->x[0] + points[k]->x[0])*(points[i]->x[1] - points[k]->x[1]);

		if(points[i]->x[1] <= x[1]){
			if(points[k]->x[1] > x[1])
				if( isLeft(points[i],points[k],x) > 0) ++wn;
		}else{
			if(points[k]->x[1] <= x[1])
				if( isLeft(points[i],points[k],x) < 0) --wn;
		}
	}

	*area = fabs(*area)*0.5;
	delete points;

	return wn;
}

/*
 * writes in four files the critical curves and the caustics for all the curves found and also for a
 * specified one (ind_causic)
 */
void writeCurves(int m			/// part of te filename, could be the number/index of the main lens
		,ImageInfo *critical	/// the crit curve
		,int Ncrit				/// the number of crit curves
		,int ind_caustic		/// the index of the cuvre of interest
		){


		  // caustic points of the fov

		  std::stringstream f;
		  f << "caustic_points_fov_" << m << ".data";
		  std::string filecrit = f.str();
		  std::ofstream file_crit(filecrit.c_str());
		  if(!file_crit){
			  std::cout << "unable to open file " << filecrit << std::endl;
			exit(1);
		  }

		  int ii,l;

		  for(ii=0;ii<Ncrit;++ii){
			  for(critical[ii].imagekist->MoveToTop(),l=0;l<critical[ii].imagekist->Nunits();l++,critical[ii].imagekist->Down()){
				  file_crit << 180/pi*3600.*critical[ii].imagekist->getCurrent()->image->x[0] << " ";
				  file_crit << 180/pi*3600.*critical[ii].imagekist->getCurrent()->image->x[1] << std::endl;
			  }
		  }

		  file_crit.close();

		  // critical points of the fov

		  f.clear();
		  f.str("");
		  f << "critical_points_fov_" << m << ".data";
		  filecrit = f.str();
		  file_crit.open(filecrit.c_str());
		  if(!file_crit){
			  std::cout << "unable to open file " << filecrit << std::endl;
			  exit(1);
		  }

	 	  for(ii=0;ii<Ncrit;++ii){
			  for(critical[ii].imagekist->MoveToTop(),l=0;l<critical[ii].imagekist->Nunits();l++,critical[ii].imagekist->Down()){
				  file_crit << 180/pi*3600.*critical[ii].imagekist->getCurrent()->x[0] << " ";
				  file_crit << 180/pi*3600.*critical[ii].imagekist->getCurrent()->x[1] << std::endl;
			  }
		  }
		  file_crit.close();

		  // caustic points of the analens

		  f.clear();
		  f.str("");
		  f << "caustic_points_" << m << ".data";
		  filecrit = f.str();
		  file_crit.open(filecrit.c_str());
		  if(!file_crit){
			  std::cout << "unable to open file " << filecrit << std::endl;
			  exit(1);
		  }

		  for(critical[ind_caustic].imagekist->MoveToTop(),l=0;l<critical[ind_caustic].imagekist->Nunits();l++,critical[ind_caustic].imagekist->Down()){
			  file_crit << 180/pi*3600.*critical[ind_caustic].imagekist->getCurrent()->image->x[0] << " ";
			  file_crit << 180/pi*3600.*critical[ind_caustic].imagekist->getCurrent()->image->x[1] << std::endl;
		  }

		  file_crit.close();

		  // critical points of the analens

		  f.clear();
		  f.str("");
		  f << "critical_points_" << m << ".data";
		  filecrit = f.str();
		  file_crit.open(filecrit.c_str());
		  if(!file_crit){
			  std::cout << "unable to open file " << filecrit << std::endl;
			  exit(1);
		  }

		  for(critical[ind_caustic].imagekist->MoveToTop(),l=0;l<critical[ind_caustic].imagekist->Nunits();l++,critical[ind_caustic].imagekist->Down()){
			  file_crit << 180/pi*3600.*critical[ind_caustic].imagekist->getCurrent()->x[0] << " ";
			  file_crit << 180/pi*3600.*critical[ind_caustic].imagekist->getCurrent()->x[1] << std::endl;
		  }

		  file_crit.close();

}
}
