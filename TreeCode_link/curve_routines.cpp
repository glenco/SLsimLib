/*
 * curve_routines.c
 *
 *  Created on: Jan 15, 2010
 *      Author: R.B. Metcalf
 */

#include "slsimlib.h"
#include "concave_hull.h"

/**  orders points in a curve, separates disconnected curves
 *   curves[0...Maxcurves] must be allocated before
 *
 *	uses neighbors-of-neighbors to split into curves and then
 *	uses sorts by angle and then walks the ist sorting
 *
 *  can break down for crescent curves
 */

//void split_order_curve4(OldImageInfo *curves,int Maxcurves,int *Ncurves){
//
//  //OldImageInfo* curves = new OldImageInfo(curves_in);
//
//  long i,m,j,end;
//  //short spur,closed,attach;
//  unsigned long NpointsTotal;
//  PosType center[2],*theta;
//  //bool delta,tmp,step;
//  //ListHndl reservoir,orderedlist;
//  //Point *newpointarray;
//
//  //std::printf("entering split_order_curve\n");
//  if(curves[0].Npoints==0){
//    *Ncurves=0;
//    return;
//  }
//
//  NpointsTotal=curves[0].Npoints;
//
//  splitter(curves,Maxcurves,Ncurves);
//
//  /*
//   reservoir=NewList();
//   orderedlist=NewList();
//
//   // copy points into a list
//   for(i=0;i<NpointsTotal;++i) InsertPointAfterCurrent(reservoir,&(curves[0].points[i]));
//   m=0;
//   i=0;
//
//   // divide points into disconnected curves using neighbors-of-neighbors
//   while(reservoir->Npoints > 0 && i < Maxcurves){
//   curves[i].points=TakeOutCurrent(reservoir);
//   MoveToBottomList(orderedlist);
//   InsertPointAfterCurrent(orderedlist,curves[i].points);
//   MoveDownList(orderedlist);
//   NeighborsOfNeighbors(orderedlist,reservoir);
//   curves[i].Npoints=orderedlist->Npoints - m;
//   //std::printf("curves[%i].Npoints=%i %i %i\n",i,curves[i].Npoints,orderedlist->Npoints
//   //		,reservoir->Npoints);
//   m+=curves[i].Npoints;
//   ++i;
//   }
//   *Ncurves=i;
//
//   // copy list back into array
//   point=curves[0].points;
//   newpointarray=NewPointArray(NpointsTotal);
//   MoveToTopList(orderedlist);
//   m=0;
//   i=0;
//   do{
//   if(i < *Ncurves && curves[i].points==orderedlist->current){
//			curves[i].points=&(newpointarray[m]);
//			++i;
//   }
//   PointCopyData(&(newpointarray[m]),orderedlist->current);
//   ++m;
//   }while(MoveDownList(orderedlist));
//
//   assert(m == NpointsTotal);
//   */
//  // order curve points
//  theta=(PosType *)malloc(NpointsTotal*sizeof(PosType));
//  assert(theta);
//  for(i=0;i<*Ncurves;++i){
//
//    if(curves[i].Npoints > 3){
//      // sort points by angle around center point
//      center[0]=center[1]=0;
//      for(m=0;m<curves[i].Npoints;++m){
//        center[0]+=curves[i].points[m].x[0];
//        center[1]+=curves[i].points[m].x[1];
//      }
//      center[0]/=curves[i].Npoints;
//      center[1]/=curves[i].Npoints;
//
//      for(m=0;m<curves[i].Npoints;++m){
//        theta[m]=atan2(curves[i].points[m].x[1]-center[1]
//                       ,curves[i].points[m].x[0]-center[0]);
//      }
//      Utilities::quicksortPoints_multithread<4>(curves[i].points,theta,curves[i].Npoints);
//
//      // check to make sure the center is inside the curves
//      //assert(abs(windings(center,curves[i].points,curves[i].Npoints,&tmp1,0)));
//
//      //std::printf("N=%i\n",curves[i].Npoints);
//      //assert(abs(windings(center,curves[i].points,curves[i].Npoints,&tmp1,0)) > 0 );
//
//      // find the last point that is a neighbor of the first point
//      m=curves[i].Npoints-1;
//      while(!AreBoxNeighbors(&(curves[i].points[0]),&(curves[i].points[m]))) --m;
//      end=m;
//      //assert(m > 0);
//
//      //std::printf("windings = %i\n",windings(center,curves[i].points,curves[i].Npoints,&tmp1,0));
//      // walk curve to remove shadowing effect
//      j=0;
//      //end=0;
//      m=0;
//      while(j < curves[i].Npoints-1){
//        walkcurve(curves[i].points,curves[i].Npoints,&j,&end);
//        //std::printf("i = %i Npoints = %i end+1 = %i j = %i\n",i,curves[i].Npoints,end+1,j);
//        if(j < curves[i].Npoints-1)
//          backtrack(curves[i].points,curves[i].Npoints,&j,-1,&end);
//        ++m;
//        assert(m < curves[i].Npoints);
//      }
//      //std::printf("i = %i Npoints = %i end+1 = %i j=%i\n",i,curves[i].Npoints,end+1,j);
//      curves[i].Npoints=end+1;
//    }
//  }
//
//  free(theta);
//  //free(reservoir);
//  //free(orderedlist);
//  //free(point);
//
//  return ;
//}
namespace Utilities{

  /** 
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
   * cut off are at the end of the array in no guaranteed order.
   *
   * Returns the number of point that have been ordered - total number minus the cuttout points.
   *
   * This algorithm could be improve by inserting the remaining points, if any, into the existing curve and recursively calling itself.
   */
  unsigned long order_curve4(Point *curve,long Npoints){
    
    long m,j,end;
    PosType center[2],*theta;
    
    //std::printf("entering split_order_curve\n");
    if(Npoints < 3) return Npoints;
    
    
    // order curve points
    theta=(PosType *)malloc(Npoints*sizeof(PosType));
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
    Utilities::quicksortPoints_multithread<4>(curve,theta,Npoints);
    
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
      // TODO: This is a kluge.  There should be a better exit stratagy.
      if(m >= Npoints) return false;
      assert(m < Npoints);
    }
    
    free(theta);
    
    return end+1;
  }
  /** \brief Overloads and is dependent on version that takes a point array.  Returns number of points that have been ordered.
   * Remaining, unordered points are left at the end of the kist.
   */
  unsigned long order_curve4(Kist<Point> * curve){
    unsigned long i=0,Npoints = curve->Nunits(),newnumber;
    
    if(Npoints < 3) return Npoints;
    
    MemmoryBank<Point> pointfactory;
    Point *tmpcurve = pointfactory(curve->Nunits());
    
    curve->MoveToTop();
    do{
      PointCopyData(&tmpcurve[i++],curve->getCurrent());
    }while(curve->Down());
    
    newnumber = order_curve4(tmpcurve,curve->Nunits());
    
    // resort points in imagekist to match tmpcurve
    for(int i=0;i<Npoints;++i){
      curve->JumpDown(i);
      do{
        if(tmpcurve[Npoints-1-i].id == curve->getCurrent()->id){
          curve->MoveCurrentToTop();
          break;
        }
      }while(curve->Down());
    }
    
    //FreePointArray(tmpcurve,false);
    
    return newnumber;
  }
  
  /**
   *  \brief For odering the curve by the convex hull method.  Warning: Does not work very well.
   *
   *   The convex hull is found for the points in the kist.  Then each additional point
   *   is inserted into the curve where it will increase the length of the curve the least.
   *   This method leaves loops where they shouldn't be and probably doesn't handle self-intersections
   *   well.
   */
  unsigned long order_curve5(Kist<Point> * curve){
    std::vector<Point *> copy;
    for(curve->MoveToTop();!(curve->OffBottom());curve->Down()){
      copy.push_back(curve->getCurrent());
    }
    
    std::vector<Point *> hull = Utilities::convex_hull(copy);
    std::vector<Point *>::iterator hit,it_min;
    
    PosType ro,d1,d2,rmin,sq2=0.99999*sqrt(2.);
    size_t i;
    bool tag;
    while(hull.size() < copy.size()){
      for(std::vector<Point *>::iterator cit = copy.begin() ; cit != copy.end() ; ++cit){
        tag = true;
        rmin = std::numeric_limits<PosType>::max();
        
        for(hit = hull.begin(),i=0 ; i < hull.size()-1 ; ++hit,++i){
          
          if(*cit == hull[i] || *cit == hull[i+1]){
            tag = false;
            break;
          }
          
          
          d1 = sqrt( ((*cit)->x[0] - hull[i]->x[0])*((*cit)->x[0] - hull[i]->x[0])
                    + ((*cit)->x[1] - hull[i]->x[1])*((*cit)->x[1] - hull[i]->x[1]) );
          
          d2 = sqrt( (hull[i+1]->x[0] - (*cit)->x[0])*(hull[i+1]->x[0] - (*cit)->x[0])
                    + (hull[i+1]->x[1] - (*cit)->x[1])*(hull[i+1]->x[1] - (*cit)->x[1]) );
          
          if(d1 > ((*cit)->gridsize + hull[i]->gridsize)/sq2 && d2 > ((*cit)->gridsize + hull[i+1]->gridsize)/sq2){
            tag = false;
            break;
          }
          
          ro = sqrt( (hull[i+1]->x[0] - hull[i]->x[0])*(hull[i+1]->x[0] - hull[i]->x[0])
                    + (hull[i+1]->x[1] - hull[i]->x[1])*(hull[i+1]->x[1] - hull[i]->x[1]) );
          
          if(rmin > (d1 + d2 - ro)){
            rmin = (d1 + d2 - ro);
            it_min = hit;
          }
        }
        
        if(tag){
          
          d1 = sqrt( ((*cit)->x[0] - hull[i]->x[0])*((*cit)->x[0] - hull[i]->x[0])
                    + ((*cit)->x[1] - hull[i]->x[1])*((*cit)->x[1] - hull[i]->x[1]) );
          d2 = sqrt( (hull[0]->x[0] - (*cit)->x[0])*(hull[0]->x[0] - (*cit)->x[0])
                    + (hull[0]->x[1] - (*cit)->x[1])*(hull[0]->x[1] - (*cit)->x[1]) );
          
          if(d1 > ((*cit)->gridsize + hull[i]->gridsize)/sq2 && d2 > ((*cit)->gridsize + hull[0]->gridsize)/sq2){
            tag = false;
            break;
          }
          
          ro = sqrt( (hull[0]->x[0] - hull[i]->x[0])*(hull[0]->x[0] - hull[i]->x[0])
                    + (hull[0]->x[1] - hull[i]->x[1])*(hull[0]->x[1] - hull[i]->x[1]) );
          
          if(rmin > (d1 + d2 - ro)){
            rmin = (d1 + d2 - ro);
            it_min = hit;
          }
          
          hull.insert((it_min+1),*cit);
        }
      }
    }
    assert(hull.size() == copy.size());
    
    curve->copy(hull);
    
    return hull.size();
  }
  
  /// Replaces curve->imagekist with its convex hull.  The number of points will change.
  void ordered_convexhull(Kist<Point> * curve){
    //int i;
    std::vector<Point *> copy = curve->copytovector();
    //copy.resize(curve->Nunits());
    //for(i=0,curve->MoveToTop();!(curve->OffBottom());curve->Down(),++i){
    //  copy[i] = curve->getCurrent();
    //}
    
    std::vector<Point *> hull = Utilities::convex_hull(copy);
    curve->copy(hull);
    
    return;
  }
  
  /// Replaces curve with its convex hull.  The number of points will change.
/*  void ordered_concavehull(Kist<Point> * curve){
    std::vector<Point *> copy = curve->copytovector();
    std::vector<Point *> hull = Utilities::concave_hull(copy,10);
    curve->copy(hull);
    
    return;
  }*/
  
  /// gives the area within the convex hull of the curve
  PosType ConvexHullArea(Kist<Point> * curve){
    Kist<Point> tmpkist;
    PosType area = 0;
    
    tmpkist.copy(curve);
    
    ordered_convexhull(curve);
    
    windings((*tmpkist)->x, curve, &area);
    
    return fabs(area);
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
 ,PosType *area          /// Area within exterior boundary
 ){
 
	long m,j,end,k;
	PosType center[2],*theta,tmp;
 
	cout << AreBoxNeighbors(&(curve[0]),&(curve[Npoints-1])) << endl;
 
	//std::printf("entering split_order_curve\n");
	if(Npoints < 3) return false;
 
 
	// order curve points
	theta=(PosType *)malloc(Npoints*sizeof(PosType));
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
	PosType ymin = curve[0].x[1];
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
/* 
 *
 * \brief Finds area within a curve by summing every cell.
 *
 * Should work every time provided the curve ordering is correct.
 * Testing if each cell is inside the curve can be slow.
 */
/*PosType findAreaOfCurve(TreeHndl tree,ImageInfo *curve,int NimageMax){
 
	if(curve->imagekist->Nunits() < 3) return 0.0;
 
	int Nimages,i,imax=0;
	PosType xcm[2],area,tmp;
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
	/ *******
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
	newpointarray=NewPointArray(curves[0].Npoints);
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
 }*/

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
      for(k=i;k>(*j)+1;--k){
        assert(k < Npoints);
        SwapPointsInArray( &(points[k]) , &(points[k-1]) );
      }
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
	PosType mintheta,x,y,phi;
 
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
          for(k=i;k>(*j)+1;--k){
            assert(k < Npoints);
            SwapPointsInArray( &(points[k]),&(points[k-1]) );
          }
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
	PosType *image_number_array,tmp;
 
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
	borders->points=NewPointArray(Npoints);
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
 image_number_array=(PosType *)malloc(Npoints_tot*sizeof(PosType));
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
 image_number_array[k]=(PosType) j;
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
	PosType *image_number_array,tmp,r,rmin;
 
	Npoints_tot=images->Npoints;
 
	// find borders of total image
	findborders3(i_tree,images);
 
	// copy inner borders to borders array
	borders=(OldImageInfo *) malloc(Maximages*sizeof(OldImageInfo));
	assert(borders);
	Npoints=images->innerborder->Nunits();
 
	// make a copy of the border points in a point array
	borders->points=NewPointArray(Npoints);
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
 image_number_array=(PosType *)malloc(Npoints_tot*sizeof(PosType));
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
//void splitter(OldImageInfo *images,int Maximages,int *Nimages){
//
//  PointList imagelist;
//  unsigned long NpointsTotal=0;
//  //Point *newpointarray;
//
//  //assert(images[0].Npoints);
//
//  if(images[0].Npoints==0){
//    *Nimages=0;
//    return;
//  }
//
//  NpointsTotal = images[0].Npoints;
//
//  PointList::iterator imagelist_it;
//  imagelist_it.current = imagelist.Bottom();
//  // copy image points into a list
//  for(long i=0;i<images[0].Npoints;++i) imagelist.InsertPointAfterCurrent(imagelist_it,&(images[0].points[i]));
//
//  assert(imagelist.size() == images[0].Npoints);
//  //std::printf("imagelist = %il\n",imagelist->Npoints);
//
//  splitlist(&imagelist,images,Nimages,Maximages);
//  //std::printf("imagelist = %il NpointsTotal = %il\n",imagelist->Npoints,NpointsTotal);
//  assert(imagelist.size() == NpointsTotal);
//
//  // copy list back into array
//  //MemmoryBank<Point> factory;
//  //newpointarray = NewPointArray(NpointsTotal);
//  //newpointarray = factory(NpointsTotal);
//  std::vector<Point> newpointarray(NpointsTotal);
//
//  imagelist_it.current = imagelist.Top();
//  {
//    long i=0,m=0;
//    do{
//      if(i < *Nimages && images[i].points == *imagelist_it){
//        images[i].points=&(newpointarray[m]);
//        ++i;
//      }
//      PointCopyData(&(newpointarray[m]),*imagelist_it);
//      ++m;
//    }while(--imagelist_it);
//
//    assert(m == NpointsTotal);
//  }
//  // take out images that have no points in them
//  for(long i=0;i<*Nimages;++i){
//    if(images[i].Npoints==0){
//      for(long j=i;j<*Nimages-1;++j){
//        images[j].ShouldNotRefine=images[j+1].ShouldNotRefine;
//        images[j].Npoints=images[j+1].Npoints;
//        images[j].points=images[j+1].points;
//        images[j].area=images[j+1].area;
//        images[j].area_error=images[j+1].area_error;
//        for(int m=0;m<3;++m) images[j].gridrange[m]=images[j+1].gridrange[m];
//      }
//      --(*Nimages);
//      --i;
//    }
//  }
//
//  return ;
//}

/** reorders imagelist into separate images using reliable
 *      neighbor-of-neighbor method
 *
 * in images the number of points in each image is updated
 *      and the pointer to the first point in each image
 *
 */
void splitlist(ListHndl imagelist,OldImageInfo *images,int *Nimages,int Maximages){
  unsigned long i=0,m=0;
  ListHndl orderedlist = new PointList;
  
  //std::printf("imagelist = %li\n",imagelist->Npoints);
  // divide images into disconnected curves using neighbors-of-neighbors
  
  PointList::iterator imagelist_it;
  imagelist_it.current = imagelist->Top();
  
  PointList::iterator orderedlist_it;
  
  while(imagelist->size() > 0 && i < Maximages){
    images[i].points = imagelist->TakeOutCurrent(imagelist_it);
    orderedlist_it.current = orderedlist->Bottom();
    orderedlist->InsertPointAfterCurrent(orderedlist_it,images[i].points);
    --orderedlist_it;
    NeighborsOfNeighbors(orderedlist,imagelist);
    images[i].Npoints = orderedlist->size() - m;
    m += images[i].Npoints;
    ++i;
  }
  *Nimages = i;
  
  // if there are too many images put the extra points in one last image
  if(i == Maximages && imagelist->size() > 0){
    Point *point;
    do{
      point = imagelist->TakeOutCurrent(imagelist_it);
      //MoveToBottomList(orderedlist);
      orderedlist_it.current = orderedlist->Bottom();
      orderedlist->InsertPointAfterCurrent(orderedlist_it,point);
      --orderedlist_it;
      ++(images[Maximages-1].Npoints);
    }while(imagelist->size() > 0);
  }
  
  if(i >= Maximages){
    ERROR_MESSAGE();
    std::printf("Too many images for Note enough images\n");
  }
  
  assert(imagelist->size() == 0);
  
  imagelist->setN(orderedlist->size());
  imagelist->setTop(orderedlist->Top());
  imagelist->setBottom(orderedlist->Bottom());
  
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
  /** 
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
               PosType *x              /// Point for which the winding number is calculated
               ,Point *points         /// The points on the border.  These must be ordered.
               ,unsigned long Npoints /// number of points in curve
               ,PosType *area          /// returns absolute the area within the curve with oriented border
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
               PosType *x              /// Point for which the winding number is calculated
               ,Point **points         /// The points on the border.  These must be ordered.
               ,unsigned long Npoints /// number of points in curve
               ,PosType *area          /// returns absolute the area within the curve with oriented border
               ,short image           /// if == 1 the image of the curve is uses as the curve
  ){
    int wn=0;
    unsigned long k,i;
    
    *area=0.0;
    if(Npoints < 3) return 0;
    
    if(image){
      for(i=0;i<Npoints;++i){
        k= i < Npoints-1 ? i+1 : 0;
        *area+=(points[i]->image->x[0] + points[k]->image->x[0])
        *(points[i]->image->x[1] - points[k]->image->x[1]);
        
        if(points[i]->image->x[1] <= x[1]){
          if(points[k]->image->x[1] > x[1])
            if( isLeft(points[i]->image,points[k]->image,x) > 0) ++wn;
        }else{
          if(points[k]->image->x[1] <= x[1])
            if( isLeft(points[i]->image,points[k]->image,x) < 0) --wn;
        }
      }
    }else{
      
      for(i=0;i<Npoints;++i){
        k= i < Npoints-1 ? i+1 : 0;
        *area+=(points[i]->x[0] + points[k]->x[0])*(points[i]->x[1] - points[k]->x[1]);
        
        if(points[i]->x[1] <= x[1]){
          if(points[k]->x[1] > x[1])
            if( isLeft(points[i],points[k],x) > 0) ++wn;
        }else{
          if(points[k]->x[1] <= x[1])
            if( isLeft(points[i],points[k],x) < 0) --wn;
        }
      }
      
    }
    
    *area = fabs(*area)*0.5;
    //std::printf("wn = %i\n",wn);
    //if(abs(wn) > 0) exit(0);
    return wn;
  }
  int windings(
               const Point_2d &x              /// Point for which the winding number is calculated
               ,const std::vector<Point_2d> &point         /// The points on the border.  These must be ordered.
               ,PosType *area          /// returns absolute the area within the curve with oriented border
  ){
    int wn=0;
    unsigned long k,i;
    size_t Npoints = point.size();
    
    *area=0.0;
    if(Npoints < 3) return 0;

    for(i=0;i<Npoints;++i){
      k= i < Npoints-1 ? i+1 : 0;
      *area+=(point[i][0] + point[k][0])*(point[i][1] - point[k][1]);
      
      if(point[i][1] <= x[1]){
        if(point[k][1] > x[1])
          if( isLeft(point[i],point[k],x) > 0) ++wn;
      }else{
        if(point[k][1] <= x[1])
          if( isLeft(point[i],point[k],x) < 0) --wn;
      }
    }
    *area = fabs(*area)*0.5;
    //std::printf("wn = %i\n",wn);
    //if(abs(wn) > 0) exit(0);
    return wn;
  }
int windings(
             Point_2d &x              /// Point for which the winding number is calculated
             ,std::vector<RAY> &point         /// The points on the border.  Uses image plane points
             ,PosType *area          /// returns absolute the area within the curve with oriented border
){
  int wn=0;
  unsigned long k,i;
  size_t Npoints = point.size();
  
  *area=0.0;
  if(Npoints < 3) return 0;

  for(i=0;i<Npoints;++i){
    k= i < Npoints-1 ? i+1 : 0;
    *area+=(point[i].x[0] + point[k].x[0])*(point[i].x[1] - point[k].x[1]);
    
    if(point[i].x[1] <= x[1]){
      if(point[k].x[1] > x[1])
        if( isLeft(point[i].x,point[k].x,x) > 0) ++wn;
    }else{
      if(point[k].x[1] <= x[1])
        if( isLeft(point[i].x,point[k].x,x) < 0) --wn;
    }
  }
  *area = fabs(*area)*0.5;
  //std::printf("wn = %i\n",wn);
  //if(abs(wn) > 0) exit(0);
  return wn;
}
  int windings(
               PosType *x              /// Point for which the winding number is calculated
               ,Kist<Point> * kist         /// Kist of points on the border.  These must be ordered.
               ,PosType *area          /// returns absolute the area within the curve with oriented border
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
      delete[] points;
      
      return wn;
    }
  
  /** 
   *
   *  \brief determines whether a point is inside a curve, that has been stretched 1.2 times
   *  returns the area of the stretched curve
   */
  int windings2(
                PosType *x              /// Point for which the winding number is calculated
                ,Point *points_original         /// The points on the border.  These must be ordered.
                ,unsigned long Npoints /// number of points in curve
                ,PosType *area          /// returns absolute the area within the curve with oriented border
                ,short image           /// if == 0 the image of the curve is used as the curve
		){
      int wn=0;
      unsigned long k,i;
      PosType center[2];
      
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
      delete[] points;
      
      return wn;
    }
  

  /**
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
          file_crit << 180/PI*3600.*critical[ii].imagekist->getCurrent()->image->x[0] << " ";
          file_crit << 180/PI*3600.*critical[ii].imagekist->getCurrent()->image->x[1] << std::endl;
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
          file_crit << 180/PI*3600.*critical[ii].imagekist->getCurrent()->x[0] << " ";
          file_crit << 180/PI*3600.*critical[ii].imagekist->getCurrent()->x[1] << std::endl;
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
        file_crit << 180/PI*3600.*critical[ind_caustic].imagekist->getCurrent()->image->x[0] << " ";
        file_crit << 180/PI*3600.*critical[ind_caustic].imagekist->getCurrent()->image->x[1] << std::endl;
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
        file_crit << 180/PI*3600.*critical[ind_caustic].imagekist->getCurrent()->x[0] << " ";
        file_crit << 180/PI*3600.*critical[ind_caustic].imagekist->getCurrent()->x[1] << std::endl;
      }
      
      file_crit.close();
      
    }
  
  // 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
  // Returns a positive value, if OAB makes a counter-clockwise turn,
  // negative for clockwise turn, and zero if the points are collinear.
  PosType cross(const Point *O, const Point *A, const Point *B)
  {
    return (A->x[0] - O->x[0]) * (B->x[1] - O->x[1]) - (A->x[1] - O->x[1]) * (B->x[0] - O->x[0]);
  }
  
  bool xorder(Point *p1,Point *p2){
    return p1->x[0] < p2->x[0];
  }
  bool yorder(Point *p1,Point *p2){
    return p1->x[1] < p2->x[1];
  }
  
  PosType crossD(const double *O, const double *A, const double *B)
  {
    return (A[0] - O[0]) * (B[1] - O[1]) - (A[1] - O[1]) * (B[0] - O[0]);
  }
  PosType crossD(const Point_2d &O,const Point_2d &A,const Point_2d &B)
  {
    return (A[0] - O[0]) * (B[1] - O[1]) - (A[1] - O[1]) * (B[0] - O[0]);
  }
  bool xorderD(double *p1,double *p2){
    return p1[0] < p2[0];
  }
  bool yorderD(double *p1,double *p2){
    return p1[1] < p2[1];
  }
}
/// Returns a vector of points on the convex hull in counter-clockwise order.
std::vector<Point *> Utilities::convex_hull(std::vector<Point *> &P)
{
  
  if(P.size() <= 3){
    std::vector<Point *> H = P;
    return H;
  }
  
  size_t n = P.size();
  size_t k = 0;
  std::vector<Point *> H(2*n);
  
  // Sort points lexicographically
  std::sort(P.begin(), P.end(), [](Point *p1,Point *p2){
    return p1->x[0] < p2->x[0];});

  
  // Build lower hull
  for (size_t i = 0; i < n; i++) {
    while (k >= 2 && cross(H[k-2], H[k-1], P[i]) <= 0){
      k--;
    }
    H[k++] = P[i];
  }
  
  // Build upper hull
  for (long i = n-2, t = k+1; i >= 0; i--) {
    while (k >= t && cross(H[k-2], H[k-1], P[i]) <= 0){
      k--;
    }
    H[k++] = P[i];
  }
  
  
  H.resize(k);
  H.pop_back();
  
  return H;
}



/// Returns a vector of points on the convex hull in counter-clockwise order.
std::vector<double *> Utilities::convex_hull(std::vector<double *> &P)
{
  
  if(P.size() <= 3){
    std::vector<double *> H = P;
    return H;
  }
  
  size_t n = P.size();
  size_t k = 0;
  std::vector<double *> H(2*n);
  
  // Sort points lexicographically
  std::sort(P.begin(), P.end(),
            [](double *p1,double *p2){return p1[0] < p2[0];});
  
  
  // Build lower hull
  for (size_t i = 0; i < n; i++) {
    while (k >= 2 && crossD(H[k-2], H[k-1], P[i]) <= 0){
      k--;
    }
    H[k++] = P[i];
  }
  
  // Build upper hull
  for (long i = n-2, t = k+1; i >= 0; i--) {
    while (k >= t && crossD(H[k-2], H[k-1], P[i]) <= 0){
      k--;
      assert(k > 1);
    }
    H[k++] = P[i];
  }
  
  
  H.resize(k);
  H.pop_back();
  
  return H;
}

double Utilities::interior_mass(const std::vector<Point_2d> &alpha
                     ,const std::vector<Point_2d> &x){
  
  long n=x.size();
  double ans = 0;
  for(int i=0 ; i<n-1 ; ++i){
    ans += (alpha[i+1]+alpha[i])^(x[i+1]-x[i]);
  }
  ans += (alpha[0]+alpha[n-1])^(x[0]-x[n-1]);

  return ans / 4;
}

/// this returns area within the curve x average kappa iwithin the curve
double Utilities::interior_mass(const std::vector<RAY> &rays){
  
  long n=rays.size();
  double ans = 0;
  for(int i=0 ; i<n-1 ; ++i){
    ans += ( rays[i+1].alpha()+rays[i].alpha() )^( rays[i+1].x-rays[i].x );
  }
  ans += ( rays[0].alpha()+rays[n-1].alpha() )^( rays[0].x-rays[n-1].x );

  return ans / 4;
}

/*
 double * convex_hull(double P1[],double P2[],std::vector<double *> P)
 {
 
 if(P.size() == 0) return NULL;
 
 if(P.size() == 1) return P[0];
 
 size_t n = P.size();
 size_t k = 0;
 std::vector<double *> H(2*n);
 
 // Sort points lexicographically
 std::sort(P.begin(), P.end(), xorderD);
 
 
 
 
 // Build lower hull
 for (size_t i = 0; i < n; i++) {
 while (k >= 2 && crossD(H[k-2], H[k-1], P[i]) <= 0){
 k--;
 }
 H[k++] = P[i];
 }
 
 // Build upper hull
 for (long i = n-2, t = k+1; i >= 0; i--) {
 while (k >= t && crossD(H[k-2], H[k-1], P[i]) <= 0){
 k--;
 }
 H[k++] = P[i];
 }
 
 
 H.resize(k);
 H.pop_back();
 
 return H[0];
 }
 */
/** \brief Returns a vector of points on the convcave hull in counter-clockwise order.
 
 This uses a K-nearest neighbour adapted from Moreira & Santos (GRAPP 2007 conference
 proceedings).  This is a modified gift wrap algorithm using k neighbours.  The value
 of k will automatically increase when certain special cases are encountered.
 
 This is an overloaded version of the other concave_hull()
 
 */
std::vector<double *> Utilities::concave_hull(std::vector<double *> &P,int k )
{
  
  if(P.size() <= 3){
    std::vector<double *> hull(P);
    return hull;
  }
  
  if(k  < 3) k =3;
  if( k  > P.size() ) k  = P.size();
  
  size_t j = 0;
  std::vector<double *> hull;
  
  std::vector<double> slist(k);
  std::vector<size_t> indexlist(k);
  double s,v1[2],v2[2];
  int tmp_k ;
  bool intersects = true,failed=false;
  
  
  
  do{
    std::vector<double *> Res(P);
    hull.resize(0);
    std::sort(Res.rbegin(), Res.rend(),
              [](double *p1,double *p2){return p1[1] < p2[1];});

    hull.push_back(Res.back());
    
    //std::cout << Res[0][1] << "  " << Res.back()->x[1] << std::endl;
    
    tmp_k = k;
    failed = false;
    j=0;
    while(Res.size() > 0 && (hull[0] != hull.back() || hull.size() == 1)){
      
      if(tmp_k  > Res.size()) tmp_k  = Res.size();
      size_t imin;
      
      slist.resize(tmp_k);
      indexlist.resize(tmp_k);
      // find list of nearest neighbors
      for(size_t kk=0;kk<tmp_k;++kk){
        slist[kk] = (Res[kk][0]- hull[j][0])*(Res[kk][0]- hull[j][0])
        + (Res[kk][1]- hull[j][1])*(Res[kk][1]- hull[j][1]);
      }
      Utilities::sort_indexes(slist,indexlist);
      std::sort(slist.begin(),slist.end());
      
      for(size_t kk=tmp_k ; kk<Res.size() ; ++kk){
        s = (Res[kk][0]- hull[j][0])*(Res[kk][0]- hull[j][0])
        + (Res[kk][1]- hull[j][1])*(Res[kk][1]- hull[j][1]);
        if(s < slist.back() && Res[kk] != hull[j]){
          slist.back() = s;
          indexlist.back() = kk;
          int ii = slist.size() - 1 ;
          while(ii > 0 && slist[ii] < slist[ii-1]){
            std::swap(slist[ii],slist[ii-1]);
            std::swap(indexlist[ii],indexlist[ii-1]);
            --ii;
          }
        }
      }
      
      //std::cout << "closest points form hull point " << j << std::endl;
      //for(int ii=0;ii<tmp_k;++ii) std::cout << indexlist[ii] << "  " << slist[ii] << " x: "
      //  << Res[indexlist[ii]][0] << "  " << Res[indexlist[ii]][1] << std::endl;
      
      // find which one has the furthest right hand turn
      if(j > 0){
        v1[0] = hull[j][0] - hull[j-1][0];
        v1[1] = hull[j][1] - hull[j-1][1];
      }else{
        v1[0] = 0;
        v1[1] = 1.0;
      }
      
      intersects = true;
      
      while(intersects && tmp_k > 0){
        
        v2[0] = Res[indexlist[0]][0] - hull[j][0];
        v2[1] = Res[indexlist[0]][1] - hull[j][1];
        
        double smin = Utilities::Geometry::AngleBetween2d( v1, v2 );
        imin = indexlist[0];
        int i_min = 0;
        
        for(int i=1;i<tmp_k;++i){
          v2[0] = Res[indexlist[i]][0] - hull[j][0];
          v2[1] = Res[indexlist[i]][1] - hull[j][1];
          
          s = Utilities::Geometry::AngleBetween2d( v1, v2 );
          if(s <= smin ){
            imin = indexlist[i];
            smin = s;
            i_min = i;
            
          }
        }
        // check for self intersection
        intersects = false;
        if(hull.size() > 2){
          for(int ii=0;ii<hull.size()-2;++ii){
            if(Utilities::Geometry::intersect(hull.back(),Res[imin], hull[ii], hull[ii+1])){
              intersects = true;
              
              //move this point to back and try again with lower tmp_k
              std::swap(indexlist[i_min],indexlist[tmp_k-1]);
              --tmp_k;
              break;
            }
          }
        }
      }
      
      
      if(!intersects){
        // add point to Hull
        std::swap(Res[imin],Res.back());
        hull.push_back(Res.back());
        ++j;
        Res.pop_back();
        tmp_k = k;
      }else{
        failed = true;
        ++k;
        
        
        break;
      }
    }
    if(!failed){
      // check to make sure all remaining points are within the hull
      
      for(size_t ii=0;ii<Res.size();++ii){
        //std::cout << Res[ii][0] << " " << Res[ii][1] << std::endl;
        if(Utilities::inhull(Res[ii],hull) == 0){
          failed = true;
          ++k;
          
          break;
        }
      }
    }
    // try again with larger k if not all the points are within the hull or no non-self-intersecting path was found
  }while(failed && k < P.size()-1);
  
  // if all else fails use the convex hull
  if(failed) hull = Utilities::convex_hull(P);
  else hull.pop_back();
  
  return hull;
}
/** \brief Returns a vector of points on the convcave hull in counter-clockwise order.
 
 This uses a K-nearest neighbour adapted from Moreira & Santos (GRAPP 2007 conference
 proceedings).  This is a modified gift wrap algorithm using k neighbours.  The value
 of k will automatically increase when certain special cases are encountered.
 
 */
std::vector<Point *> Utilities::concave_hull(std::vector<Point *> &P,int k,bool test )
{
  PixelMap<float> *testmap;
  int nt=0;
  
  
  if(P.size() <= 3){
    std::vector<Point *> hull(P);
    return hull;
  }
  
  if(k  < 3) k =3;
  if( k  > P.size() ) k  = P.size();
  
  size_t j = 0;
  std::vector<Point *> hull;
  
  std::vector<double> slist(k);
  std::vector<size_t> indexlist(k);
  double s,v1[2],v2[2];
  int tmp_k ;
  bool intersects = true,failed=false;
  
  if(test){
    PosType bound_p1[2],bound_p2[2],center[2];
    bound_p1[0] = bound_p2[0] = P[0]->x[0];
    bound_p1[1] = bound_p2[1] = P[0]->x[1];
    for(size_t ii = 1; ii<P.size();++ii){
      bound_p1[0] = MIN(P[ii]->x[0],bound_p1[0]);
      bound_p2[0] = MAX(P[ii]->x[0],bound_p2[0]);
      bound_p1[1] = MIN(P[ii]->x[1],bound_p1[1]);
      bound_p2[1] = MAX(P[ii]->x[1],bound_p2[1]);
    }
    center[0] = (bound_p2[0] + bound_p1[0])/2;
    center[1] = (bound_p2[1] + bound_p1[1])/2;
    testmap = new PixelMap<float>(center,512
                           ,1.05*MAX((bound_p2[0] - bound_p1[0]),(bound_p2[1] - bound_p1[1]))/512);
    
    testmap->drawPoints(P,0,1.0);
    testmap->printFITS("!test_"+std::to_string(nt++)+".fits");
  }
  
  
  do{
    std::vector<Point *> Res(P);
    hull.resize(0);
    std::sort(Res.rbegin(), Res.rend(),
              [](Point *p1,Point *p2){return p1->x[1] < p2->x[1];});

    hull.push_back(Res.back());
    
    //std::cout << Res[0]->x[1] << "  " << Res.back()->x[1] << std::endl;
    
    tmp_k = k;
    failed = false;
    j=0;
    while(Res.size() > 0 && (hull[0] != hull.back() || hull.size() == 1)){
      
      if(tmp_k  > Res.size()) tmp_k  = Res.size();
      size_t imin;
      
      slist.resize(tmp_k);
      indexlist.resize(tmp_k);
      // find list of nearest neighbors
      for(size_t kk=0;kk<tmp_k;++kk){
        slist[kk] = (Res[kk]->x[0]- hull[j]->x[0])*(Res[kk]->x[0]- hull[j]->x[0])
        + (Res[kk]->x[1]- hull[j]->x[1])*(Res[kk]->x[1]- hull[j]->x[1]);
      }
      Utilities::sort_indexes(slist,indexlist);
      std::sort(slist.begin(),slist.end());
      
      for(size_t kk=tmp_k ; kk<Res.size() ; ++kk){
        s = (Res[kk]->x[0]- hull[j]->x[0])*(Res[kk]->x[0]- hull[j]->x[0])
        + (Res[kk]->x[1]- hull[j]->x[1])*(Res[kk]->x[1]- hull[j]->x[1]);
        if(s < slist.back() && Res[kk] != hull[j]){
          slist.back() = s;
          indexlist.back() = kk;
          int ii = slist.size() - 1 ;
          while(ii > 0 && slist[ii] < slist[ii-1]){
            std::swap(slist[ii],slist[ii-1]);
            std::swap(indexlist[ii],indexlist[ii-1]);
            --ii;
          }
        }
      }
      
      //std::cout << "closest points form hull point " << j << std::endl;
      //for(int ii=0;ii<tmp_k;++ii) std::cout << indexlist[ii] << "  " << slist[ii] << " x: "
      //  << Res[indexlist[ii]]->x[0] << "  " << Res[indexlist[ii]]->x[1] << std::endl;
      
      // find which one has the furthest right hand turn
      if(j > 0){
        v1[0] = hull[j]->x[0] - hull[j-1]->x[0];
        v1[1] = hull[j]->x[1] - hull[j-1]->x[1];
      }else{
        v1[0] = 0;
        v1[1] = 1.0;
      }
      
      intersects = true;
      
      while(intersects && tmp_k > 0){
        
        v2[0] = Res[indexlist[0]]->x[0] - hull[j]->x[0];
        v2[1] = Res[indexlist[0]]->x[1] - hull[j]->x[1];
        
        double smin = Utilities::Geometry::AngleBetween2d( v1, v2 );
        imin = indexlist[0];
        int i_min = 0;
        
        for(int i=1;i<tmp_k;++i){
          v2[0] = Res[indexlist[i]]->x[0] - hull[j]->x[0];
          v2[1] = Res[indexlist[i]]->x[1] - hull[j]->x[1];
          
          s = Utilities::Geometry::AngleBetween2d( v1, v2 );
          if(s <= smin ){
            imin = indexlist[i];
            smin = s;
            i_min = i;
            
          }
        }
        // check for self intersection
        intersects = false;
        if(hull.size() > 2){
          for(int ii=0;ii<hull.size()-1;++ii){
            if(Utilities::Geometry::intersect(hull.back()->x,Res[imin]->x, hull[ii]->x, hull[ii+1]->x)){
              intersects = true;
              
              if(test){
                std::cout << "smin/pi = " << smin/PI << std::endl;
                std::cout << "intersecting line segments ii = " << ii << std::endl;
                std::cout << hull.back()->x[0] << " " << hull.back()->x[1] << " -- "
                << Res[imin]->x[0] << " " << Res[imin]->x[1] << std::endl;
                
                std::cout << hull[ii]->x[0] << " " << hull[ii]->x[1] << " -- "
                << hull[ii+1]->x[0] << " " << hull[ii+1]->x[1] << std::endl;
              }
              
              //move this point to back and try again with lower tmp_k
              std::swap(indexlist[i_min],indexlist[tmp_k-1]);
              --tmp_k;
              break;
            }
          }
        }
      }
      
      
      if(!intersects){
 
        assert(tmp_k >0);
        // add point to Hull
        std::swap(Res[imin],Res.back());
        hull.push_back(Res.back());
        ++j;
        Res.pop_back();
        tmp_k = k;
        
        if(test){
          for(int ii=0;ii<hull.size()-1;++ii){
            for(int jj=ii+1;jj<hull.size()-1;++jj){
              if(Utilities::Geometry::intersect(hull[jj]->x,hull[jj+1]->x, hull[ii]->x, hull[ii+1]->x)){
                std::cout << "ii: " << ii << "  jj: " << jj << std::endl;
                throw std::runtime_error("Hull should not self intersect");
              }
            }
          }
        }

      }else{
        assert(tmp_k == 0);

        failed = true;
        ++k;
        tmp_k = k;
        
        if(test){
          testmap->Clean();
          testmap->drawPoints(P,0,1.0);
          testmap->drawCurve(hull,2.0);
          testmap->printFITS("!test_"+std::to_string(nt)+".fits");
          
          assert(hull[j] != hull[j-1]);
          PixelMap<float> map(hull.back()->x,512,
                       MAX(fabs(hull[j]->x[0] - hull[j-1]->x[0]),fabs(hull[j]->x[1] - hull[j-1]->x[1]))/5 );
          
          map.drawPoints(P,0,1.0);
          map.drawCurve(hull,2.0);
          map.printFITS("!test_"+std::to_string(nt++)+".1.fits");
          for(int ii=0;ii<hull.size()-1;++ii) std::cout << hull[ii]->x[0] << " " << hull[ii]->x[1] << std::endl;
          
        }
        
        break;
      }
    }
    if(!failed){
      // check to make sure all remaining points are within the hull
      
      for(size_t ii=0;ii<Res.size();++ii){
        //std::cout << Res[ii]->x[0] << " " << Res[ii]->x[1] << std::endl;
        if(inhull(Res[ii]->x,hull) == 0){
          failed = true;
          ++k;
          
          if(test){
            
            testmap->Clean();
            testmap->drawPoints(P,0,1.0);
            testmap->drawCurve(hull,2.0);
            testmap->printFITS("!test_"+std::to_string(nt)+".fits");
            
            PixelMap<float> map(Res[ii]->x,512,
                         MAX(fabs(hull[j]->x[0] - hull[j-1]->x[0]),fabs(hull[j]->x[1] - hull[j-1]->x[1]))/5 );
            
            map.drawPoints(P,0,1.0);
            //map.drawCurve(hull,2.0);
            map.drawPoints(hull,0,2.0);
            map.printFITS("!test_"+std::to_string(nt++)+".1.fits");
            for(int ii=0;ii<hull.size()-1;++ii) std::cout << hull[ii]->x[0] << " " << hull[ii]->x[1] << std::endl;
            
          }
          
          break;
        }
      }
    }
    // try again with larger k if not all the points are within the hull or no non-self-intersecting path was found
  }while(failed && k < P.size()-1);
  
  
  // if all else fails use the convex hull
  if(failed) hull = Utilities::convex_hull(P);
  else hull.pop_back();
  
  if(test){
    for(int ii=0;ii<hull.size()-1;++ii){
      for(int jj=ii+1;jj<hull.size()-1;++jj){
        if(Utilities::Geometry::intersect(hull[jj]->x,hull[jj+1]->x, hull[ii]->x, hull[ii+1]->x)){
          std::cout << "ii: " << ii << "  jj: " << jj << std::endl;
          throw std::runtime_error("Hull should not self intersect");
        }
      }
    }
  }
  return hull;
}


/** \brief Returns axis ratio, area and points of an ellipse engulfed by some contour (e.g. a contour of same convergence calculated with find_contour).
 
    The axis ratio of the ellipse b/a is equal to the ratio of the distances between center and the nearest contour point (i.e. b) and between center and the farthest contour point (i.e. a).
    NOTE that the center used to calculate a and b is an input parameter. The definition of the center is crucial to the meaning of above output parameters. The center of the convex_hull
    produces for even slightly distorted hulls significant offsets resulting in overestimated major axis values (a).
    The function Utilities::contour_center() calculates the center as the midpoint between the two points in the contour that are farthest apart from one another, which gives already more reliable results.
    The output vector describing the ellipse is resized to match the size of the contour vector.
 
 */

void Utilities::contour_ellipse(std::vector<Point_2d> &P, Point_2d center, unsigned long Npoints ,std::vector<Point_2d> &C, double *ellipticity, double *ellipse_area){
  double nearest_contour_pnt=1E12;
  double mostdistant_contour_pnt=0.0;
  double tmp_dist=1E12;
  double t,pa,farpoint[2],nearpoint[2];
  
  
  for(size_t jj=0;jj<Npoints;++jj){
    tmp_dist=sqrt((P[jj].x[0]-center.x[0])*(P[jj].x[0]-center.x[0])+(P[jj].x[1]-center.x[1])*(P[jj].x[1]-center.x[1]));
    if (nearest_contour_pnt>tmp_dist && tmp_dist!=0){nearest_contour_pnt=tmp_dist;nearpoint[0]=P[jj].x[0];nearpoint[1]=P[jj].x[1];}
    if (mostdistant_contour_pnt<tmp_dist){mostdistant_contour_pnt=tmp_dist;farpoint[0]=P[jj].x[0];farpoint[1]=P[jj].x[1];}
  }
  *ellipticity=nearest_contour_pnt/mostdistant_contour_pnt;
  *ellipse_area=nearest_contour_pnt*mostdistant_contour_pnt*PI;
  pa=atan2f(farpoint[1]-center.x[1],farpoint[0]-center.x[0]);
  
  C.resize(Npoints);
  
  for(size_t jj=0;jj<Npoints;++jj){
      t=jj*2.*PI/Npoints;
      C[jj].x[0]=cos(pa)*mostdistant_contour_pnt*cos(t)-sin(pa)*nearest_contour_pnt*sin(t);
      C[jj].x[1]=sin(pa)*mostdistant_contour_pnt*cos(t)+cos(pa)*nearest_contour_pnt*sin(t);
  }
}
void Utilities::contour_ellipse(std::vector<RAY> &P, Point_2d center, unsigned long Npoints ,std::vector<Point_2d> &C, double *ellipticity, double *ellipse_area){
  double nearest_contour_pnt=1E12;
  double mostdistant_contour_pnt=0.0;
  double tmp_dist=1E12;
  double t,pa,farpoint[2],nearpoint[2];
  
  
  for(size_t jj=0;jj<Npoints;++jj){
    tmp_dist=sqrt((P[jj].x[0]-center.x[0])*(P[jj].x[0]-center.x[0])+(P[jj].x[1]-center.x[1])*(P[jj].x[1]-center.x[1]));
    if (nearest_contour_pnt>tmp_dist && tmp_dist!=0){nearest_contour_pnt=tmp_dist;nearpoint[0]=P[jj].x[0];nearpoint[1]=P[jj].x[1];}
    if (mostdistant_contour_pnt<tmp_dist){mostdistant_contour_pnt=tmp_dist;farpoint[0]=P[jj].x[0];farpoint[1]=P[jj].x[1];}
  }
  *ellipticity=nearest_contour_pnt/mostdistant_contour_pnt;
  *ellipse_area=nearest_contour_pnt*mostdistant_contour_pnt*PI;
  pa=atan2f(farpoint[1]-center.x[1],farpoint[0]-center.x[0]);
  
  C.resize(Npoints);
  
  for(size_t jj=0;jj<Npoints;++jj){
      t=jj*2.*PI/Npoints;
      C[jj].x[0]=cos(pa)*mostdistant_contour_pnt*cos(t)-sin(pa)*nearest_contour_pnt*sin(t);
      C[jj].x[1]=sin(pa)*mostdistant_contour_pnt*cos(t)+cos(pa)*nearest_contour_pnt*sin(t);
  }
}

/** \brief Returns the center of a contour defined as the midpoint between the two points in the contour that are farthest apart from one another.
 
    The performance of the algorithm is ~O(N^2). Less naive methods go like O(N) at best. Most commonly a combined convex hull plus rotating calipers algorithm is used.
    Since we have the convex_hull already, we only need to implement the latter algorithm.
 */


Point_2d Utilities::contour_center(std::vector<Point_2d> &P, unsigned long Npoints){
  double dPx,dPy,dist,jx[2],ix[2];
  double maxdist=0.0;
  Point_2d center;
  
  for(size_t jj=0;jj<Npoints;++jj){
    for(size_t ii=0;ii<Npoints;++ii){
      dPx=P[jj].x[0]-P[ii].x[0];
      dPy=P[jj].x[1]-P[ii].x[1];
      dist=sqrt(dPx*dPx+dPy*dPy);
      if (dist>maxdist){maxdist=dist;jx[0]=P[jj].x[0];jx[1]=P[jj].x[1];ix[0]=P[ii].x[0];ix[1]=P[ii].x[1];}
    }
  }
  
  center.x[0]=0.5*(ix[0]+jx[0]);
  center.x[1]=0.5*(ix[1]+jx[1]);
  return center;
}

Point_2d Utilities::contour_center(std::vector<RAY> &P, unsigned long Npoints){
  double dPx,dPy,dist,jx[2],ix[2];
  double maxdist=0.0;
  Point_2d center;
  
  for(size_t jj=0;jj<Npoints;++jj){
    for(size_t ii=0;ii<Npoints;++ii){
      dPx=P[jj].x[0]-P[ii].x[0];
      dPy=P[jj].x[1]-P[ii].x[1];
      dist=sqrt(dPx*dPx+dPy*dPy);
      if (dist>maxdist){maxdist=dist;jx[0]=P[jj].x[0];jx[1]=P[jj].x[1];ix[0]=P[ii].x[0];ix[1]=P[ii].x[1];}
    }
  }
  
  center.x[0]=0.5*(ix[0]+jx[0]);
  center.x[1]=0.5*(ix[1]+jx[1]);
  return center;
}

void Utilities::find_islands(std::vector<bool> &bitmap  // = true inside
                  ,long nx  // number of pixels in x direction
                  ,std::vector<std::vector<long> > &indexes
                  ,std::vector<bool> &hits_edge
                  ,bool add_to_vector
                  ){
  
  size_t n = bitmap.size();
  long ny = n/nx;
  
  if(n != nx*ny){
    std::cerr << "Wrong sizes in Utilities::find_boundaries." << std::endl;
    throw std::invalid_argument("invalid size");
  }
  
  std::vector<bool> not_used(n,true);
  
  // pad edge of field with bitmap=false
  for(size_t i=0 ; i<nx ; ++i) bitmap[i]=false;
  size_t j = nx*(ny-1);
  for(size_t i=0 ; i<nx ; ++i) bitmap[i + j]=false;
  for(size_t i=0 ; i<ny ; ++i) bitmap[i*nx]=false;
  j = nx-1;
  for(size_t i=0 ; i<ny ; ++i) bitmap[j + i*nx]=false;
  
  std::list< std::set<long> > contours;  // a set ensures that it is in assending order and no element is repeated
  
  if(!add_to_vector){
    hits_edge.resize(0);
  }
  
  bool done = false;
  long kfirst_in_bound = -1;
  while(!done){
    // find first cell in edge
    size_t k=0,k1,k2,k3;
    int type;
    for( k = kfirst_in_bound + 1 ; k < n - nx ; ++k){
      if(k % nx != nx-1){ // one less cells than points
        type = 0;
        if(bitmap[k] ) type +=1;
        if(bitmap[k+1]) type += 10;
        if(bitmap[k + nx]) type += 100;
        if(bitmap[k + nx + 1]) type += 1000;
        
        if(type > 0
           && type != 1111
           && not_used[k]
           ) break;
      }
    }
    
    kfirst_in_bound = k;
    
    if(k == n-nx){
      done=true;
    }else{ // found an edge
      
      contours.resize(contours.size() + 1);
      std::set<long> &contour = contours.back();
      hits_edge.push_back(false);
      
      
      int type;
      int face_in=0;
      size_t n_edge = 0;
      
      // follow edge until we return to the first point
      while(k != kfirst_in_bound || n_edge==0){
        
        if(k%nx == 0 || k%nx == nx-2) hits_edge.back() = true;
        if(k/nx == 0 || k/nx == ny-2) hits_edge.back() = true;
        
        not_used[k] = false;
        
        ++n_edge;
        type = 0;
        
        k1 = k+1;
        k2 = k+nx;
        k3 = k+nx+1;
        
        // find type of cell
        if(bitmap[k]) type +=1;
        if(bitmap[k1]) type += 10;
        if(bitmap[k2]) type += 100;
        if(bitmap[k3]) type += 1000;
        
        if(type == 0 || type == 1111){  // all in or all out
          throw std::runtime_error("off edge!!");
        }else if(type == 1 || type == 1110){ // lower left only
          
          if(type == 1) contour.insert(k);
          if(type == 1110){
            contour.insert(k1);
            contour.insert(k2);
            contour.insert(k3);
          }
          
          if(face_in==0){
            face_in=1;
            k -= nx;
          }else{
            face_in=2;
            k -= 1;
          }
          
        }else if(type == 10 || type == 1101){ // lower right only
          
          if(type == 10) contour.insert(k1);
          if(type == 1101){
            contour.insert(k);
            contour.insert(k2);
            contour.insert(k3);
          }
          
          if(face_in==2){
            face_in=1;
            k -= nx;
          }else{
            face_in=0;
            k += 1;
          }
          
        }else if(type == 100 || type == 1011){ // upper left only
          
          if(type == 100) contour.insert(k2);
          if(type == 1011){
            contour.insert(k);
            contour.insert(k1);
            contour.insert(k3);
          }
          
          if(face_in==0){
            face_in=3;
            k += nx;
          }else{
            face_in=2;
            k -= 1;
          }
          
        }else if(type == 1000 || type == 111){ // upper right only
          
          if(type == 1000) contour.insert(k3);
          if(type == 1011){
            contour.insert(k);
            contour.insert(k1);
            contour.insert(k2);
          }
          
          if(face_in==1){
            face_in=0;
            k += 1;
          }else{
            face_in=3;
            k += nx;
          }
          
        }else if(type == 11 || type == 1100){ // lower two
          
          if(type == 11){
            contour.insert(k);
            contour.insert(k1);
          }
          if(type == 1100){
            contour.insert(k2);
            contour.insert(k3);
          }
          
          if(face_in==0){
            k += 1;
          }else{
            face_in = 2;
            k -= 1;
          }
          
        }else if(type == 1010 || type == 101){ // right two
          
          if(type == 1010){
            contour.insert(k1);
            contour.insert(k3);
          }
          if(type == 101){
            contour.insert(k);
            contour.insert(k2);
          }
          
          if(face_in==1){
            k -= nx;
          }else{
            face_in = 3;
            k += nx;
          }
          
        }else if(type == 1001){ // lower left upper right
          
          contour.insert(k);
          contour.insert(k3);
          
          if(face_in==0){
            face_in=3;
            k += nx;
          }else if(face_in==1){
            face_in=2;
            k -= 1;
          }else if(face_in==2){
            face_in=1;
            k -= nx;
          }else{
            face_in=0;
            k += 1;
          }
          
        }else if(type == 110){ // upper left lower right
          
          contour.insert(k1);
          contour.insert(k2);
          
          if(face_in==0){
            face_in=1;
            k -= nx;
          }else if(face_in==1){
            face_in=0;
            k += 1;
          }else if(face_in==2){
            face_in=3;
            k += nx;
          }else{
            face_in=2;
            k -= 1;
          }
        }
      }
      
   
      // fill in island
      indexes.emplace_back();
      std::vector<long> &index = indexes.back();
      
      long k=0;
      for(auto it = contour.begin()
            ; it != contour.end()
            ; ++it
            ){
        if(k<*it){
          k=*it;
          while(bitmap[k]){
            index.push_back(k);
            ++k;
          }
        }
      }
      
      
    }
  }
}

Point_2d Utilities::line_intersection(const Point_2d &v1,const Point_2d &v2,
                            const Point_2d &w1,const Point_2d &w2){
  Point_2d dv = v2-v1;
  Point_2d dw = w2-w1;
  double s = ( dw^(w1-v1) ) / (dw^dv);
  return v1 + dv*s;
}

bool Utilities::circleIntersetsCurve(const Point_2d &x,double r,const std::vector<Point_2d> &v){
  
  int n=v.size();
  double rr,B;
  Point_2d dx,dl;
  for(int i=0 ; i<n ; ++i){
    dl = v[i]-x;
    rr = dl.length();
    if(rr <= r) return true;
    dx = v[(i+1)%n]-v[i];
    B = -(dx*dl)/dx.length_sqr();
    if(B>0 && B<1){
      if( dx*dx*B*B + dx*dl*2*B + rr*rr <= r*r) return true;
    }
  }
  
  return false;
}

bool Utilities::circleOverlapsCurve(const Point_2d &x,double r,const std::vector<Point_2d> &v){

  double area;
  // is center of circle inside curve
  if(Utilities::windings(x,v,&area)) return true;
  // if curve entirely within circle this will be true
  if( (x-v[0]).length() < r) return true;
  
  return Utilities::circleIntersetsCurve(x,r,v);
}

std::vector<std::vector<Point_2d> > Utilities::thicken_poly(
                                const std::vector<Point_2d> &v
                                ,double R){
  
  std::vector<std::vector<Point_2d> > output;
  if(v.size() == 0) return output;
  if(R <= 0){
    output.push_back(v);
    return output;
  }
  
  // find bounding box
  Point_2d ll,ur;
  ll = ur = v[0];
  for(const Point_2d &p : v){
    ll[0] = MIN(p[0],ll[0]);
    ll[1] = MIN(p[1],ll[1]);
    
    ur[0] = MAX(p[0],ur[0]);
    ur[1] = MAX(p[1],ur[1]);
  }
  
  double resolution=R/10,f=1.2;
  double caustic_size = MIN(ur[0]-ll[0],ur[1]-ll[1]);
  
  //resolution = MAX(resolution,caustic_size/100);
  
  long nx = (ur[0]-ll[0])/resolution;
  long ny = (ur[1]-ll[1])/resolution;
  long nv = v.size();
  
  if(nx > 500 || ny > 500){
    output.resize(2);
    output[0].resize(2*nv);
    output[1].resize(2*nv);
    for(long i=0 ; i<nv ; ++i){
      Point_2d d = v[ (i+1)%nv ] - v[i];
      d *= R/d.length();
      std::swap(d[0],d[1]);
      d[0] *= -1;
      
      Point_2d x = (v[ (i+1)%nv ] - v[i])*0.1  + v[i];
      output[0][2*i] = x + d;
      output[1][2*i] = x - d;
      
      x = (v[ (i+1)%nv ] - v[i])*0.9  + v[i];
      output[0][2*i + 1] = x + d;
      output[1][2*i + 1] = x - d;
 
    }
    
  }else{
    
    bool too_small = true;
    while(too_small ){
      too_small= false;
      ll[0] -= f*R;
      ll[1] -= f*R;
      ur[0] += f*R;
      ur[1] += f*R;
      
      nx = (ur[0]-ll[0])/resolution;
      ny = (ur[1]-ll[1])/resolution;
      
      long count=0;
      Point_2d p = ll;
      std::vector<bool> bitmap(nx*ny,false);
      for(long j=0 ; j<ny ; ++j,p[1] += resolution){
        p[0] = ll[0];
        for(long i=0 ; i<nx ; ++i,p[0] += resolution){
          long m = i + nx*j;
          for(long k=0 ; k < nv ; ++k){
            if(R > Utilities::distance_to_segment(p, v[k],v[ (k+1)%nv ]) ){
              bitmap[m] = true;
              ++count;
              break;
            }
          }
        }
      }
      
      std::vector<bool> hits_edge;
      Utilities::find_boundaries(bitmap,nx,output,hits_edge);
      // check that it doesn't touch an edge
      for(bool a : hits_edge){
        if(a){
          too_small = true;
          f *= 1.5;
        }
      }
      
    }
    
    // rescale
    for(int i = 0 ; i<output.size() ; ++i){
      for(Point_2d &p : output[i]) p = ll + p * resolution;
    }
  }
  return output;
}

std::vector<Point_2d> Utilities::envelope2(const std::vector<Point_2d> &v
                               ,const std::vector<Point_2d> &w){
  
  size_t nv = v.size();
  size_t nw = w.size();
  //std::cout << nv << " " << nw << std::endl;
  if(nv<3) return w;
  if(nw<3) return v;
  
  Utilities::Geometry::CYCLIC cycv(nv);
  Utilities::Geometry::CYCLIC cycw(nw);
  Point_2d inter_p;
  
  bool intersecting = false;
  // check if they intersect
  long i,j;
  for(i=0; i<nv && !intersecting ; ++i){
    for(j=0; j<nw && !intersecting ; ++j){
      if(Utilities::Geometry::intersect(v[i].x,v[ cycv[i+1] ].x
                                        ,w[j].x,w[ cycw[j+1] ].x))
        intersecting = true;
        inter_p = line_intersection(v[i].x,v[ cycv[i+1] ].x
                        ,w[j].x,w[ cycw[j+1] ].x);
    }
  }
  
  if(intersecting == false){
    if( Utilities::inCurve<Point_2d>(w[0],v)){
      return v;
    }
    
    if( Utilities::inCurve<Point_2d>(v[0],w)){
      return w;
    }

    // case where they do not overlap
    return std::vector<Point_2d>(0);
  }
  
  // concatenate curves
  std::vector<Point_2d> curve(w.size() + v.size() + 2);
  long ii=0;
  while(ii < nw ){
    curve[ii++] = w[ cycw[j++] ];
  }
  curve[ii++] = inter_p;  // insert extra point where curves intersect
  while(ii < curve.size() - 1 ){
    curve[ii++] = v[ cycv[i++] ];
  }
  
  curve[ii] = inter_p + ( curve[ii-1] - inter_p) * 0.01;  // insert extra point close to where curves intersect

  try{
    return Utilities::TighterHull(curve);
  }catch(...){
    return Utilities::convex_hull(curve);
  }
}

std::vector<Point_2d> Utilities::envelope(const std::vector<Point_2d> &v
                               ,const std::vector<Point_2d> &w){
  
  size_t nv = v.size();
  size_t nw = w.size();
  //std::cout << nv << " " << nw << std::endl;
  if(nv<3) return w;
  if(nw<3) return v;
  
  Utilities::Geometry::CYCLIC cycv(nv);
  Utilities::Geometry::CYCLIC cycw(nw);
  Point_2d inter_p;
  
  bool intersecting = false;
  // check if they intersect
  long i,j;
  for(i=0; i<nv && !intersecting ; ++i){
    for(j=0; j<nw && !intersecting ; ++j){
      if(Utilities::Geometry::intersect(v[i].x,v[ cycv[i+1] ].x
                                        ,w[j].x,w[ cycw[j+1] ].x))
        intersecting = true;
        inter_p = line_intersection(v[i].x,v[ cycv[i+1] ].x
                        ,w[j].x,w[ cycw[j+1] ].x);
    }
  }
  
  if(intersecting == false){
    if( Utilities::inCurve<Point_2d>(w[0],v)){
      return v;
    }
    
    if( Utilities::inCurve<Point_2d>(v[0],w)){
      return w;
    }

    // case where they do not overlap
    return std::vector<Point_2d>(0);
  }
  
  // this is basicly TightestHull but using the segments of both curves
  
  // find bounding box
  Point_2d ll,ur;
  ll = ur = v[0];
  for(const Point_2d &p : v){
    ll[0] = MIN(p[0],ll[0]);
    ll[1] = MIN(p[1],ll[1]);
    
    ur[0] = MAX(p[0],ur[0]);
    ur[1] = MAX(p[1],ur[1]);
  }
  for(const Point_2d &p : w){
    ll[0] = MIN(p[0],ll[0]);
    ll[1] = MIN(p[1],ll[1]);
    
    ur[0] = MAX(p[0],ur[0]);
    ur[1] = MAX(p[1],ur[1]);
  }

  std::vector<Point_2d> output;
  double resolution =  MIN(ur[0]-ll[0],ur[1]-ll[1]) /100;
  if(resolution==0){
    output.push_back(ll);
    output.push_back(ur);
    return output;
  }
    
  double R = resolution*sqrt(2);
  
  ll[0] -= 2*R;
  ll[1] -= 2*R;
  ur[0] += 2*R;
  ur[1] += 2*R;
  
  long nx = (ur[0]-ll[0])/resolution;
  long ny = (ur[1]-ll[1])/resolution;
  
  long count=0;
  Point_2d p = ll;
  std::vector<bool> bitmap(nx*ny,false);
  for(long j=0 ; j<ny ; ++j,p[1] += resolution){
    p[0] = ll[0];
    for(long i=0 ; i<nx ; ++i,p[0] += resolution){
      long m = i + nx*j;
      
      for(long k=0 ; k < nv ; ++k){
        if(R > Utilities::distance_to_segment(p, v[k],v[ (k+1)%nv ] ) ){
          bitmap[m] = true;
          ++count;
          break;
        }
      }
      if( bitmap[m] == false){
        for(long k=0 ; k < nw ; ++k){
          if(R > Utilities::distance_to_segment(p, w[k],w[ (k+1)%nw ] ) ){
            bitmap[m] = true;
            ++count;
            break;
          }
        }
      }
    }
  }
  
  if(count > 0){
    std::vector<std::vector<Point_2d> > envelopes;
    std::vector<bool> hits_edge;
    
    Utilities::find_boundaries(bitmap,nx,envelopes,hits_edge,false,true);
    
    for(Point_2d &p : envelopes[0]) p = ll + p * resolution;
    
    size_t n = envelopes[0].size();
    output.resize(n);
    
    double min_d,tmp_d;
    for(long i=0 ; i<n ; ++i){
      
      Point_2d p,closest_point;
      double min_d = Utilities::distance_to_segment(envelopes[0][i],v[0],v[1],closest_point);
      //assert(!isnan(closest_point[0]));
      //assert(!isnan(closest_point[1]));

      for(long k=1 ; k < nv ; ++k){
        tmp_d = Utilities::distance_to_segment(envelopes[0][i],v[k],v[(k+1)%nv],p);
        if(tmp_d<min_d){
          min_d=tmp_d;
          closest_point = p;
        }
      }
      
      for(long k=0 ; k < nw ; ++k){
        tmp_d = Utilities::distance_to_segment(envelopes[0][i],w[k],w[(k+1)%nw],p);
        if(tmp_d<min_d){
          min_d=tmp_d;
          closest_point = p;
        }
      }

      
      output[i] = closest_point;
      //assert(!isnan(output[i][0]));
      //assert(!isnan(output[i][1]));
    }
  }
  
  return output;
}

std::vector<Point_2d> Utilities::TighterHull(const std::vector<Point_2d> &vv){

  if(vv.size() < 4 ) return vv;
 
  std::vector<Point_2d> v = vv;
  long nv = v.size();
  
  // find length of largest edge
  double length_scale=0;
  for(long i=0 ; i<nv ; ++i){
    double tmp = (v[i] - v[(i+1)%nv]).length();
    if(tmp > length_scale) length_scale = tmp;
  }
  
  length_scale *= 0.01;
  
  // remove repeated consecutive points
  for(long i=0 ; i<nv-1 ; ++i){
    if( (v[i] - v[i+1]).length() < length_scale ){
      for(long k=i+1 ; k<nv-1 ; ++k){
        std::swap(v[k],v[k+1]);
      }
      --nv;
    }
  }
  while( ( (v[nv-1] - v[0]).length() < length_scale ) && nv>1) --nv;
  assert(nv>1);
  
  v.resize(nv);
  
  Utilities::Geometry::CYCLIC cycv(nv);
  
  // displace points that are the same
  {
    Utilities::RandomNumbers_NR ran;
  
    for(long i=0 ; i<nv-1 ; ++i){
      for(long j=i+1 ; j<nv ; ++j){
        if( (v[i] - v[j]).length() < length_scale ){
          Point_2d dw = v[cycv[j+1]]-v[j];
          //Point_2d dv = v[cycv[i+1]]-v[i];
          v[j] += dw * length_scale *(2*ran()-1)/dw.length();
        }
      }
    }
  }
 
//  // check for colinear segments
//  for(long i=0 ; i<nv-1 ; ++i){
//    Point_2d dv = v[i+1]-v[i];
//    for(long j=i+1 ; j<nv ; ++j){
//
//      long jp = cycv[j+1];
//      if( ( dv^( v[j]-v[i]) ) == 0 && ( dv^( v[jp]-v[i]) ) == 0 ){  // colinear
//        Point_2d dw = v[jp]-v[j];
//
//        double s =( dv*(v[j]-v[i]) )/dv.length_sqr();
//        double s2 =( dv*(v[jp]-v[i]) )/dv.length_sqr();
//        double s3 =( dw*(v[i]-v[j]) )/dw.length_sqr();
//
//        if(  (s >= 0 && s <= 1)
//           || (s2 >= 0 && s2 <= 1)
//           || (s3 >= 0 && s3 <= 1)
//           ){  // overlapping
//          // add small perturbation perpendiculare to segmant
//          std::swap(dw[0],dw[1]);
//          dw[0] *= -1;
//          v[i] += dw*0.001;
//          v[i+1] -= dw*0.001;
//        }
//      }
//    }
//  }

  // find most left
  size_t imax = 0;
  {
    double tmp = v[0][0];
    for(int ii=1 ; ii<v.size() ; ++ii){
      if(tmp > v[ii][0]){
        tmp=v[ii][0];
        imax = ii;
      }
    }
  }
  size_t i = imax;
 
  //write_csv("v_pruned.csv",v);
  
  std::vector<Point_2d> env;

  env.push_back(v[i]);

  // orientation
  int o = sign( (v[cycv[imax-1]]-v[imax])^(v[cycv[imax+1]]-v[imax])  );

  int n_intersect=0,step=0;
  long last_i = cycv[imax-o];
  if(o==-1) last_i = imax;

  while( env.size() <= 20*nv ){
    ++step;

    long j_int=-1;
    double s=-1;//,p;
    Point_2d dio = v[ cycv[i+o] ] - env.back();
    n_intersect=0;
    
    double s_max = 1;
    int jj_parallel = -1;
    int o_par =0;
    
    for(int jj=0; jj<nv ; ++jj){
      if(  jj != i && jj != last_i
         && jj != cycv[i+o]
         ){
        
        long jp = cycv[jj+1];
        Point_2d dj = v[ jp ] - v[jj];
        double pji = dj^dio;
        
        if(pji == 0  // parallel
           && ( ( v[jj]-env.back() )^dio ) == 0   // colinear
           ){ // colinear
          
          double sj = (v[jj]-env.back())*dio /dio.length_sqr();
          double sjp = (v[jp]-env.back())*dio /dio.length_sqr();
          
          if(MAX(sj,sjp) > s_max  // new edge extends beyond current edge
          && MIN(sj,sjp) < 1      // overlaps with current edge
             ){
            
            if(sjp > sj){
              jj_parallel = jj;
              o_par = 1;
              s_max = sjp;
            }else{
              jj_parallel = jp;
              o_par = -1;
              s_max = sj;
            }
          }
          
        }else if(
            Utilities::Geometry::intersect(env.back().x,v[ cycv[i+o] ].x
                                           ,v[jj].x,v[ jp ].x)  // does not include sheared endpoints
           
          ){
           ++n_intersect;
             
            Point_2d inter_p = line_intersection(env.back().x,v[ cycv[i+o] ].x
                                ,v[jj].x,v[ jp ].x);
             
            double s2 = ( (inter_p-env.back())*dio )/dio.length_sqr();
            if(s2>1.0) s2=1.0; // roundoff error occurs sometimes.
             
             if( s2 > 0){
               if(j_int==-1 ){
                 j_int=jj;
                 
                 //p = pji;
                 if( v[ cycv[i+o] ] == v[jj]){
                   s=1;  // avoid numerical error
                 }else{
                   s = s2;
                 }
               }else{ // case of two crossings
                 
                 //double p2 = pji;
                 //if( v[ cycv[i+o] ] == v[jj]){
                 //  s2=1;  // avoid numerical error
                 //}else{
                 //  s2 = ( (inter_p-env.back())*di )/di.length_sqr();
                 // }
                 
                 //assert(s <= 1);
                 if(s > s2){
                   j_int=jj;
                   s=s2;
                   //p=p2;
                 }
               }
               assert(s>0);
             }
           
         }
      }
    }

    if(o==1) last_i = i;
    else last_i = cycv[i-1];

    if(j_int != -1){
      
      //Point_2d  v1 = v[ cycv[i+o] ] - env.back();
      Point_2d  w1 = v[ cycv[j_int] ] - env.back();
      Point_2d  w2 = v[ cycv[j_int+1] ] - env.back();
 
      env.push_back(env.back() + dio*s );
      assert(s >= 0);
      assert(s <= 1);

      //double tmp = w1^dio;
      //Point_2d p = v[i];
      
      if( (w1^dio) < 0){
        i = cycv[j_int+1];
        o = -1;
      }else if( (w2^dio) < 0){
        i = j_int;
        o = 1;
      }else{ // case where neither point is above the line
        i = cycv[i+o];
        env.push_back(v[i]);
        if(i==imax) break;
      }
      
    }else if(jj_parallel != -1){
      i=jj_parallel;
      o = o_par;
    }else{
      // case of moving along in same direction
      i = cycv[i+o];
      env.push_back(v[i]);
      if(i==imax) break;
    }

    //write_csv("hull.csv",env);
  }

  if(env.size() > 20*nv){
    write_csv("test_curve.csv",v);
    write_csv("test_hull.csv",env);
    std::cerr << "Failure of TighterHull. See test_hull.csv and test_curve.csv." << std::endl;
    throw std::runtime_error("Failure in TighterHull");
  }
 
  //assert(env.size() > 2);

  env.pop_back();
  return env;
}


std::vector<Point_2d> Utilities::TightestHull(const std::vector<Point_2d> &v){
  if(v.size() <= 3) return v;
  
  // find bounding box
  Point_2d ll,ur;
  ll = ur = v[0];
  for(const Point_2d &p : v){
    ll[0] = MIN(p[0],ll[0]);
    ll[1] = MIN(p[1],ll[1]);
    
    ur[0] = MAX(p[0],ur[0]);
    ur[1] = MAX(p[1],ur[1]);
  }
  
  std::vector<Point_2d> output;
  double resolution =  MIN(ur[0]-ll[0],ur[1]-ll[1]) /100;
  if(resolution==0){
    output.push_back(ll);
    output.push_back(ur);
    return output;
  }
    
  double R = resolution*sqrt(2.5);
  
  ll[0] -= 2*R;
  ll[1] -= 2*R;
  ur[0] += 2*R;
  ur[1] += 2*R;
  
  long nx = (ur[0]-ll[0])/resolution;
  long ny = (ur[1]-ll[1])/resolution;
  long nv = v.size();
  
  long count=0;
  Point_2d p = ll;
  std::vector<bool> bitmap(nx*ny,false);
  for(long j=0 ; j<ny ; ++j,p[1] += resolution){
    p[0] = ll[0];
    for(long i=0 ; i<nx ; ++i,p[0] += resolution){
      long m = i + nx*j;
      
      for(long k=0 ; k < nv ; ++k){
        if(R > Utilities::distance_to_segment(p, v[k],v[ (k+1)%nv ] ) ){
          bitmap[m] = true;
          ++count;
          break;
        }
      }
    }
  }
  
  if(count > 0){
    std::vector<std::vector<Point_2d> > envelopes;
    std::vector<bool> hits_edge;
    
    Utilities::find_boundaries(bitmap,nx,envelopes,hits_edge,false,true);
    
    for(Point_2d &p : envelopes[0]) p = ll + p * resolution;
    
    size_t n = envelopes[0].size();
    output.resize(n);
    Point_2d p,closest_point;
    double min_d,tmp_d;
    for(long i=0 ; i<n ; ++i){
      
      min_d = Utilities::distance_to_segment(envelopes[0][i],v[0],v[1],closest_point);
      
      for(long k=1 ; k < nv ; ++k){
        tmp_d = Utilities::distance_to_segment(envelopes[0][i],v[k],v[(k+1)%nv],p);
        if(tmp_d<min_d){
          min_d=tmp_d;
          closest_point = p;
        }
      }
      
      output[i] = closest_point;
    }
  }
  return output;
}

Point_2d Utilities::RandomInTriangle(const Point_2d &x1,
                                     const Point_2d &x2,
                                     const Point_2d &x3,
                                     Utilities::RandomNumbers_NR &ran
                                     ){
  double u1=ran(),u2=ran();
  if(u1+u2 > 1){
    u1=1-u1;
    u2=1-u2;
  }
   
  return x1*u1 + x2*u2 + x3*(1-u1-u2);
}

Point_2d Utilities::RandomInConvexPoly(const std::vector<Point_2d> &pp,
                                     Utilities::RandomNumbers_NR &ran
               ){
  assert(pp.size() >0);
  int n=pp.size();
  if(n==2){
    return pp[0] + (pp[1]-pp[0])*ran();
  }else if(n==3){
    return RandomInTriangle(pp[0],pp[1],pp[2],ran);
  }
  std::vector<double> areas(n-1,0);
  for(int i=1 ; i<n-1 ; ++i){
    areas[i] =  areas[i-1] + abs( (pp[i]-pp[0])^(pp[i+1]-pp[0]) );
  }
  double tmp = ran()*areas.back();
  
  int i=1;
  while(tmp > areas[i]) ++i;
  
  assert(i<n-1);
  return RandomInTriangle(pp[i],pp[i+1],pp[0],ran);
}

std::vector<Point_2d> Utilities::RandomInConvexPoly(const std::vector<Point_2d> &pp,
                                                    int N,
                                                    Utilities::RandomNumbers_NR &ran
               ){
  assert(pp.size() >0);
  std::vector<Point_2d> output(N);
  int n=pp.size();
  if(n==2){
    for(Point_2d &p : output) p = pp[0] + (pp[1]-pp[0])*ran();
    return output;
  }else if(n==3){
    for(Point_2d &p : output) p =  RandomInTriangle(pp[0],pp[1],pp[2],ran);
    return output;
  }
  std::vector<double> areas(n-1,0);
  for(int i=1 ; i<n-1 ; ++i){
    areas[i] =  areas[i-1] + abs( (pp[i]-pp[0])^(pp[i+1]-pp[0]) );
  }
  for(Point_2d &p : output){
    double tmp = ran()*areas.back();
  
    int i=1;
    while(tmp > areas[i]) ++i;
  
    assert(i<n-1);
    p = RandomInTriangle(pp[i],pp[i+1],pp[0],ran);
  }
  return output;
}

Point_2d Utilities::RandomInPoly(std::vector<Point_2d> &pp,
                                Utilities::RandomNumbers_NR &ran
               ){
  assert(pp.size() >0);
  int n=pp.size();
  if(n==2){
    return pp[0] + (pp[1]-pp[0])*ran();
  }else if(n==3){
    return RandomInTriangle(pp[0],pp[1],pp[2],ran);
  }
  
  std::vector<Point_2d> hull = convex_hull(pp);
  
  n=hull.size();
  std::vector<double> areas(n-1,0);
  for(int i=1 ; i<n-1 ; ++i){
    areas[i] =  areas[i-1] + abs( (hull[i]-hull[0])^(hull[i+1]-hull[0]) );
  }
  
  Point_2d p;
  do{
    double tmp = ran()*areas.back();
  
    int i=1;
    while(tmp > areas[i]) ++i;
  
    assert(i<n-1);
    p = RandomInTriangle(hull[i],hull[i+1],hull[0],ran);
  }while( !inCurve(p,pp) );
         
  return p;
}

std::vector<Point_2d> Utilities::RandomInPoly(std::vector<Point_2d> &pp,
                                              int N,
                                              Utilities::RandomNumbers_NR &ran
               ){
  assert(pp.size() >0);
  std::vector<Point_2d> output(N);
  int n=pp.size();
  if(n==2){
    for(Point_2d &p : output) p = pp[0] + (pp[1]-pp[0])*ran();
    return output;
  }else if(n==3){
    for(Point_2d &p : output) p =  RandomInTriangle(pp[0],pp[1],pp[2],ran);
    return output;
  }
  
  std::vector<Point_2d> hull = convex_hull(pp);
  
  n=hull.size();
  std::vector<double> areas(n-1,0);
  for(int i=1 ; i<n-1 ; ++i){
    areas[i] =  areas[i-1] + abs( (hull[i]-hull[0])^(hull[i+1]-hull[0]) );
  }
  
  for(Point_2d &p : output){
    do{
      double tmp = ran()*areas.back();
      
      int i=1;
      while(tmp > areas[i]) ++i;
      
      assert(i<n-1);
      p = RandomInTriangle(hull[i],hull[i+1],hull[0],ran);
    }while( !inCurve(p,pp) );
  }
  return output;
}

Point_2d Utilities::RandomNearPoly(std::vector<Point_2d> &pp,
                                   double R,
                                   Utilities::RandomNumbers_NR &ran
               ){
  Point_2d p = RandomInPoly(pp,ran);
  
  double r = R*sqrt(ran());
  double theta = 2*PI*ran();
  
  p[0] = p[0] + r*cos(theta);
  p[1] = p[1] + r*sin(theta);
  
  return p;
}

std::vector<Point_2d> Utilities::RandomNearPoly(std::vector<Point_2d> &pp,
                                               int N,
                                               double R,
                                              Utilities::RandomNumbers_NR &ran
                                             ){
  std::vector<Point_2d> ps = RandomInPoly(pp,N,ran);
  
  for(Point_2d &p : ps){
    double r = R*sqrt(ran());
    double theta = 2*PI*ran();
    
    p[0] = p[0] + r*cos(theta);
    p[1] = p[1] + r*sin(theta);
  }
  
  return ps;
}
