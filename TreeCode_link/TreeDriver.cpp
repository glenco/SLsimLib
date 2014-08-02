
#include "slsimlib.h"

#include <nrutil.h>

/* median_cut determines how the cells are subdivided */
/*    if ==0  equal volume cuts, Warning this option causes an error*/
/*    if ==1  median point cuts */
//static int incell;

//static PosType realray[2];
//Point *point_global;

/** \ingroup ImageFundingL2
 *
 * \brief Finds nearest neighbor points to ray.
 *
 *    This is a kludge that relies on NearestNeighbor which uses a List and translates
 *    the list to a kist.  Could be rewritten.
 *
 *    Warning: The number of neighbor points in neighborkist will be less than Nneighbors when
 *             the number of points in the tree is less than Nneighbors
 */
Point *TreeStruct::NearestNeighborKist(const PosType* center,int Nneighbors,Kist<Point>* neighborkist){
  /* nearest neighbor points in a given direction, direction != 0 should not */
  /*    be used except on a grid */
  /* direction = 0 distance */
  /*           = 1 to right */
  /*           = 2 to left */
  /*           = 3 to up */
  /*           = 4 to down */
  unsigned long i;
  //static int oldNneighbors=-1;
  //static PosType *rneighbors;
  //static Point **neighborpoints;
  short direction = 0;

  /*std::printf("entering NN\n");*/

  if(top->npoints <= Nneighbors){
	  //ERROR_MESSAGE();
	  //std::printf("ERROR: in NearestNeighbor, number of neighbors > total number of points\n");
	  //throw std::out_of_range(std::string() + "ERROR: in NearestNeighbor, number of neighbors > total number of points");
    Nneighbors = top->npoints-1;
  }
//  EmptyList(neighborlist);
  neighborkist->Empty();
  if(Nneighbors <= 0) return NULL;
/*
  if(count==0){
    rneighbors= (PosType *)malloc((Nneighbors+Nbucket)*sizeof(PosType));
    assert(rneighbors);
    neighborpoints=(Point **)malloc((Nneighbors+Nbucket)*sizeof(Point *));
    assert(neighborpoints);
    //temp_points = (Point **)malloc((Nneighbors+Nbucket)*sizeof(Point *));
    assert(temp_points);
    
    ++count;
    oldNneighbors=Nneighbors;

  }else if(oldNneighbors < Nneighbors){ // if the number of nearest neighbors goes up get more mem
    rneighbors = (PosType *)realloc(rneighbors,(Nneighbors+Nbucket)*sizeof(PosType));
    neighborpoints=(Point **)realloc(neighborpoints,(Nneighbors+Nbucket)*sizeof(Point *));
    //temp_points = (Point **)realloc(temp_points,(Nneighbors+Nbucket)*sizeof(Point *));
    oldNneighbors=Nneighbors;
  }
  */
  
  PosType rneighbors[Nneighbors+Nbucket];
  Point *neighborpoints[Nneighbors+Nbucket];
  tmp_point.resize(Nneighbors+Nbucket);

  /*   std::printf("Nneighbors=%i\n",Nneighbors); */
  /*   std::printf("array sizes=%i\n",Nneighbors+Nbucket); */

  /* initalize distance to neighbors to a large number */
  for(i=0;i<(Nbucket+Nneighbors);++i){
    rneighbors[i]=10*(top->boundary_p2[0]-top->boundary_p1[0]);
  }

  moveTop();
  //   std::printf("p1= [%f,%f]\n", current->boundary_p1[0],current->boundary_p1[1]);
  //   std::printf("p2= [%f,%f]\n", current->boundary_p2[0],current->boundary_p2[1]);

  PosType ray[2] = { center[0], center[1] };
	
  realray[0]=ray[0];
  realray[1]=ray[1];

  if( inbox(ray,current->boundary_p1,current->boundary_p2) == 0 ){
    std::printf("Warning: in NearestNeighbor, ray is not inside the simulation box\n    should work in any case\n      ray= %e %e\n",ray[0],ray[1]);

    ray[0]=DMAX(ray[0],current->boundary_p1[0]);
    ray[0]=DMIN(ray[0],current->boundary_p2[0]);

    ray[1]=DMAX(ray[1],current->boundary_p1[1]);
    ray[1]=DMIN(ray[1],current->boundary_p2[1]);
  }
  incell=1;

  //if(direction==0) EmptyList(neighborlist);
  if(direction==-1) direction=0;

  _NearestNeighbor(ray,Nneighbors,neighborpoints,rneighbors,&direction);

  /* convert from point array to exported point list */

  for(i=0;i<Nneighbors;++i){
    neighborkist->InsertAfterCurrent(neighborpoints[i]);
    neighborkist->Down();
  }
/*
  for(i=0;i<Nneighbors;++i){
    InsertAfterCurrent(neighborlist,neighborpoints[i]->x,neighborpoints[i]->id,neighborpoints[i]->image);
    MoveDownList(neighborlist);
    PointCopyData(neighborlist->current,neighborpoints[i]);
  }
*/
    
  neighborkist->MoveToTop();
  return neighborkist->getCurrent();
}

void TreeStruct::_NearestNeighbor(PosType* ray,int Nneighbors,Point **neighborpoints,PosType *rneighbors,short *direction){

  int i,incell2=1;
  unsigned long index[Nneighbors+Nbucket];
  PosType dx,dy;

  //std::printf("**************************************\nlevel %i\n",current->level);
  //for(i=0;i<current->npoints;++i) std::printf("   %i\n",current->points[i]);

  if(incell){  /* not found cell yet */

    if( inbox(ray,current->boundary_p1,current->boundary_p2) ){

    	// found the box small enough */
    	if( current->npoints <= (Nneighbors+Nbucket) ){
    		incell=0;
    		/*std::printf("found box with %i points\n",current->npoints);*/

    		//if( (ray[0]!=realray[0])*(ray[1]!=realray[1]) ){ std::printf("ray != realray _NearestNeighbor\n"); ERROR_MESSAGE();}
    		//assert( (ray[0] == realray[0]) && (ray[1] == realray[1]) );

    		// this sets ray back to real value once closest leaf box is found
    		ray[0]=realray[0];
    		ray[1]=realray[1];

    		/* calculate the distance to all the points in cell */
    		if(current->points != NULL ) pointlist->current=current->points;
    		for(i=0;i<current->npoints;++i){

    			dx=pointlist->current->x[0]-ray[0];
    			dy=pointlist->current->x[1]-ray[1];

    			switch(*direction){
    			case 0: /* distance */
    				rneighbors[i] = sqrt(dx*dx+dy*dy);
    				break;
    			case 1: /* right */
    				if(dx>0 && fabs(dy/dx)<0.999) rneighbors[i]=sqrt(dx*dx+dy*dy);
    				else rneighbors[i]=1.0e99;
    				break;
    			case 2: /* left */
    				if(dx<0 && fabs(dy/dx)<0.999) rneighbors[i]=sqrt(dx*dx+dy*dy);
    				else rneighbors[i]=1.0e99;
    				break;
    			case 3: /* up */
    				if(dy>0 && fabs(dx/dy)<0.999) rneighbors[i]=sqrt(dx*dx+dy*dy);
    				else rneighbors[i]=1.0e99;
    				break;
    			case 4: /* down */
    				if(dy<0 && fabs(dx/dy)<0.999) rneighbors[i]=sqrt(dx*dx+dy*dy);
    				else rneighbors[i]=1.0e99;
    				break;
    			}

     			index[i]=i;
    			//temp_points[i]=pointlist->current;
          tmp_point[i] = pointlist->current;
    			MoveDownList(pointlist);
    		}

    		if(current->npoints > 0){
    			Utilities::double_sort(current->npoints,rneighbors-1,index-1);
    			//for(i=0;i<current->npoints;++i) neighborpoints[i] = temp_points[index[i]];
    			for(i=0;i<current->npoints;++i) neighborpoints[i] = tmp_point[index[i]];
    		}

      }else{  // keep going down the tree

    	  if(current->child1 !=NULL){
    		  moveToChild(1);
    		  _NearestNeighbor(ray,Nneighbors,neighborpoints,rneighbors,direction);
    		  moveUp();

    		  incell2=incell;
    	  }

    	  if(current->child2 !=NULL){
    		  moveToChild(2);
    		  _NearestNeighbor(ray,Nneighbors,neighborpoints,rneighbors,direction);
    		  moveUp();
    	  }

    	  // if ray found in second child go back to first to search for neighbors
    	  if( (incell2==1) && (incell==0) ){
    		  if(current->child1 !=NULL){
    			  moveToChild(1);
    			  _NearestNeighbor(ray,Nneighbors,neighborpoints,rneighbors,direction);
    			  moveUp();
    		  }
    	  }

      }
    } // not in the box

  }else{  /* already found cell */

	  // does radius cut into the box
	  if( Utilities::cutbox(ray,current->boundary_p1,current->boundary_p2,rneighbors[Nneighbors-1]) ){

		  if( atLeaf() ){  /* leaf case */

			  /* combine found neighbors with points in box and resort */
			  if(current->points != NULL) pointlist->current=current->points;

			  for(i=0;i<Nneighbors;++i){
				  index[i]=i;
				  //temp_points[i]=neighborpoints[i];
				  tmp_point[i]=neighborpoints[i];
			  }

			  for(i=Nneighbors;i<(current->npoints+Nneighbors);++i){

				  dx=pointlist->current->x[0]-ray[0];
				  dy=pointlist->current->x[1]-ray[1];

				  switch(*direction){
				  case 0: /* distance */
					  rneighbors[i] = sqrt(dx*dx+dy*dy);
					  break;
				  case 1: /* right */
					  if(dx>0 && fabs(dy/dx)<0.999) rneighbors[i]=sqrt(dx*dx+dy*dy);
					  else rneighbors[i]=1.0e99;
					  break;
				  case 2: /* left */
					  if(dx<0 && fabs(dy/dx)<0.999) rneighbors[i]=sqrt(dx*dx+dy*dy);
					  else rneighbors[i]=1.0e99;
					  break;
				  case 3: /* up */
					  if(dy>0 && fabs(dx/dy)<0.999) rneighbors[i]=sqrt(dx*dx+dy*dy);
					  else rneighbors[i]=1.0e99;
					  break;
				  case 4: /* down */
					  if(dy<0 && fabs(dx/dy)<0.999) rneighbors[i]=sqrt(dx*dx+dy*dy);
					  else rneighbors[i]=1.0e99;
					  break;
				  }

				  index[i]=i;
				  //temp_points[i]=pointlist->current;
				  tmp_point[i]=pointlist->current;
				  MoveDownList(pointlist);
			  }

			  // sort the points found so far
			  if(current->npoints > 0){
				  Utilities::double_sort(current->npoints+Nneighbors,rneighbors-1,index-1);
				  //for(i=0;i<(current->npoints+Nneighbors);++i) neighborpoints[i] = temp_points[index[i]];
				  for(i=0;i<(current->npoints+Nneighbors);++i) neighborpoints[i] = tmp_point[index[i]];
			  }

		  }else{

			  if(current->child1 !=NULL){
				  moveToChild(1);
				  _NearestNeighbor(ray,Nneighbors,neighborpoints,rneighbors,direction);
				  moveUp();
			  }

			  if(current->child2 !=NULL){
				  moveToChild(2);
				  _NearestNeighbor(ray,Nneighbors,neighborpoints,rneighbors,direction);
				  /*std::printf("moving up from level %i\n",current->level);*/
				  moveUp();
			  }
		  }

	  }

  }

  return;
}

// returns true if branch1 is fully inside barnch2
bool boxinbox(Branch *branch1,Branch *branch2){

	if( branch1 == branch2 ) return true;

	if(inbox(branch1->boundary_p1,branch2->boundary_p1,branch2->boundary_p2) == 0) return false;
	if(inbox(branch1->boundary_p2,branch2->boundary_p1,branch2->boundary_p2) == 0) return false;

	return true;
}
PosType BoxIntersection(Branch *branch1,Branch *branch2){
	// returns area of intersection between two branches
	PosType area=0;

	area = MIN(branch1->boundary_p2[0],branch2->boundary_p2[0])
	     - MAX(branch1->boundary_p1[0],branch2->boundary_p1[0]);
	if(area < 0) return 0.0;

	area *= MIN(branch1->boundary_p2[1],branch2->boundary_p2[1])
	      - MAX(branch1->boundary_p1[1],branch2->boundary_p1[1]);
	if(area < 0) return 0.0;

	return area;
}

// return 1 (0) if box is (not) within rmax of ray
int cutbox(PosType *ray,PosType *p1,PosType *p2,PosType rmax){
	/*  returns:  0 if whole box is outside rmax from ray[]
	 *            1 if whole box is inside circle but ray is not in the box
	 *            2 if ray[] is inside box
	 *            3 if box intersects circle but ray[] is not inside box
	 */
  short i,tick=0;
  PosType close[2],rtmp;
  
  // find closest point on box borders to ray[]
  for(i=0;i<2;++i){
    if( ray[i] < p1[i] ){
      close[i]=p1[i];
    }else if(ray[i] > p2[i]){
      close[i]=p2[i];
    }else{
      close[i]=ray[i];
      ++tick;
    }
  }
  
  if(tick==2) return 2;  // ray is inside box

  for(i=0,rtmp=0;i<2;++i) rtmp += pow(ray[i] - close[i],2);
 
  if(rtmp>rmax*rmax) return 0;  // box is all outside circle

  // find farthest point on box border from ray[]
  for(i=0,rtmp=0;i<2;++i) rtmp+=DMAX(pow(ray[i]-p1[i],2),pow(ray[i]-p2[i],2));

  if(rtmp<rmax*rmax) return 1;  // box is all inside circle

  return 3;  // box intersects circle
}

/**
 * \brief If the circle centered at ray with radius is entirely within the box
 * returns true.
 *
 */
bool CircleInBox(const PosType* center,PosType radius,PosType *p1,PosType *p2){

	if(!inbox(center,p1,p2)) return false;

	if((center[0] + radius) > p2[0] ) return false;
	if((center[0] - radius) < p1[0] ) return false;
	if((center[1] + radius) > p2[1] ) return false;
	if((center[1] - radius) < p1[1] ) return false;

	return true;
}

bool BoxInCircle(const PosType* center,PosType radius,PosType *p1,PosType *p2){

	PosType rad2 = radius*radius;

	if((pow(p1[0] - center[0],2) + pow(p1[1] - center[1],2)) > rad2) return false;
	if((pow(p2[0] - center[0],2) + pow(p2[1] - center[1],2)) > rad2) return false;
	if((pow(p1[0] - center[0],2) + pow(p2[1] - center[1],2)) > rad2) return false;
	if((pow(p2[0] - center[0],2) + pow(p1[1] - center[1],2)) > rad2) return false;

	return true;
}

/**
 *   Finds the leaf the ray is in and adds Nadd to all of is parent leaves
*/

void TreeStruct::_FindLeaf(const PosType* ray,unsigned long Nadd){

	bool contin;

	do{
		current->npoints += Nadd;

		/* change center of mass */
		//current->center[0]=(current->center[0]*(current->npoints-1) + ray[0])/current->npoints;
		//current->center[1]=(current->center[1]*(current->npoints-1) + ray[1])/current->npoints;

		contin=false;
		if(current->child1 != NULL &&
				inbox(ray,current->child1->boundary_p1,current->child1->boundary_p2) ){
			moveToChild(1);
			contin=true;
		}else if(current->child2 != NULL &&
				inbox(ray,current->child2->boundary_p1,current->child2->boundary_p2) ){
			moveToChild(2);
			contin=true;
		}
	}while(contin);

  return;
}

bool AreBoxNeighbors(Point *point1,Point *point2){

	if(    point1->leaf->boundary_p1[0] <= point2->leaf->boundary_p2[0]
	    && point1->leaf->boundary_p2[0] >= point2->leaf->boundary_p1[0]
	    && point1->leaf->boundary_p1[1] <= point2->leaf->boundary_p2[1]
	    && point1->leaf->boundary_p2[1] >= point2->leaf->boundary_p1[1] ) return true;
	return false;
}

bool AreBoxNeighbors(Branch *branch1,Branch *branch2){

	if(    branch1->boundary_p1[0] <= branch2->boundary_p2[0]
	    && branch1->boundary_p2[0] >= branch2->boundary_p1[0]
	    && branch1->boundary_p1[1] <= branch2->boundary_p2[1]
	    && branch1->boundary_p2[1] >= branch2->boundary_p1[1] ) return true;
	return false;
}

/** return the point that is in the same box as ray[2]
 * if Nbuck > 1 the head of the point array is returned
 */
void TreeStruct::FindBoxPoint(const PosType* ray,Point *point){
	Branch *branch;

	branch=current;
	moveTop();
	   // check if ray is outside initial box
	if( inbox(ray,current->boundary_p1,current->boundary_p2) == 0 ){
		std::printf("FindBox: ray outside of grided range\n");
		return;
	}

	_FindBox(ray);
	PointCopyData(point,current->points);

	//if(foundpoint == 0){ std::printf("FindBoxPoint failed to find point\n"); exit(1);}
	//PrintPoint(point);

	// error check
	if(fabs(ray[0]-point->x[0]) > point->gridsize/2
			|| fabs(ray[1]-point->x[1]) > point->gridsize/2){
		ERROR_MESSAGE();
		std::printf("ERROR: FindBox did not find box\n  ray = %e %e\n  Delta/gridsize = %e %e\n"
				,ray[0],ray[1]
		        ,2*(ray[0]-point->x[0])/point->gridsize
				,2*(ray[1]-point->x[1])/point->gridsize);
		exit(0);
	}

	current=branch;
	//return point;
}

void TreeStruct::_FindBox(const PosType* ray){
	bool contin;

	// replaced recursion with iteration
	do{
		contin=false;
		if(current->child1 !=NULL &&
				inbox(ray,current->child1->boundary_p1,current->child1->boundary_p2) ){
			moveToChild(1);
			contin=true;
		}
		if(current->child2 !=NULL &&
				inbox(ray,current->child2->boundary_p1,current->child2->boundary_p2) ){
			moveToChild(2);
			contin=true;
		}
	}while(contin);

	return;
}


Point *sortList(long n, PosType *arr,ListHndl list,Point *firstpointin){
	/** simple sort for points in linked list **/
	/** slow for > 20 points **/
  long i,j;
  PosType a;
  Point *point,*firstpoint;

  if(n <= 1) return firstpointin;

  firstpoint=firstpointin;

  for(j=1;j<n;j++){
    a=arr[j];

    list->current=firstpoint;
    JumpDownList(list,j);

    point=TakeOutCurrent(list);

    i=j-1;
    while(i>-1 && arr[i] > a){
      arr[i+1]=arr[i];
      i--;
      MoveUpList(list);
      /*std::printf("      current= %i %f %f\n",list->current->id,list->current->x[0],list->current->x[1]);*/
    }
    arr[i+1]=a;

    if( list->current==list->top && i==-1) InsertPointBeforeCurrent(list,point);
    else InsertPointAfterCurrent(list,point);

    if(i == -1) firstpoint=point;
  }
  return firstpoint;
}

/**
 *  finds all points in tree that lie within rmax of the point ray[]
 *   markpoints = 0  does not change in_image variable in any point, gives a list of neighbors
 *              = 1  makes in_image=true for all points in image, gives no list of neighbors
 *              = -1 makes in_image=false for all points in image to reset, gives no list of neighbors
 *
void TreeStruct::PointsWithin(PosType *ray,float rmax,ListHndl neighborlist,short markpoints){

  if(markpoints==0) EmptyList(neighborlist);
  
  realray[0]=ray[0];
  realray[1]=ray[1];

  moveTop();
  if( inbox(ray,current->boundary_p1,current->boundary_p2) == 0 ){
    std::printf("Warning: in PointsWithin, ray is not inside the simulation box\n    should work in any case\n      ray= %e %e\n     boundary p1 = %e %e p2 = %e %e\n",ray[0],ray[1]
	   ,current->boundary_p1[0],current->boundary_p1[1]
	   ,current->boundary_p2[0],current->boundary_p2[1]);

    ray[0]=DMAX(ray[0],current->boundary_p1[0]);
    ray[0]=DMIN(ray[0],current->boundary_p2[0]);

    ray[1]=DMAX(ray[1],current->boundary_p1[1]);
    ray[1]=DMIN(ray[1],current->boundary_p2[1]);
  }
  incell=1;

  _PointsWithin(ray,&rmax,neighborlist,markpoints);
}

void TreeStruct::PointsWithin_iter(PosType *ray,float rmax,ListHndl neighborlist,short markpoints){
	bool descend;
	short desition;
	unsigned long i,j;
	PosType radius,length;

	assert(neighborlist);

	EmptyList(neighborlist);

	moveTop();

	do{
		descend = true;
		desition = cutbox(ray,current->boundary_p1,current->boundary_p2,rmax);
		length = FurthestBorder(ray,current->boundary_p1,current->boundary_p2);
		if( FurthestBorder(ray,current->boundary_p1,current->boundary_p2) < rmax ){

			//ClosestBorder(ray,current->boundary_p1,current->boundary_p2) < rmax) ){
			// whole box is outside circle
			descend = false;

			// put all the points into neighborlist
		   	  pointlist->current=current->points;
    		  for(i=0;i<current->npoints;++i){
     			  if(markpoints == 1){
     				  pointlist->current->in_image=TRUE;
     				  pointlist->current->image->in_image=TRUE;
     			  }else if(markpoints == -1){
     				  pointlist->current->in_image=FALSE;
     				  pointlist->current->image->in_image=FALSE;
					  pointlist->current->surface_brightness = pointlist->current->image->surface_brightness = 0.0;
    			  }else if(markpoints == 0){
     				  InsertAfterCurrent(neighborlist,pointlist->current->x
							  ,pointlist->current->id,pointlist->current->image);

					  MoveDownList(neighborlist);
					  PointCopyData(neighborlist->current,pointlist->current);
					  MoveUpList(neighborlist);
				  }

     			  MoveDownList(pointlist);
     		  }

		}else if(cutbox(ray,current->boundary_p1,current->boundary_p2,rmax) == 0 ){  // whole box is outside circle
			descend = false;
		}else if(current->child1 == NULL){  // leaf case

		   	  pointlist->current=current->points;
	   		  for(i=0;i<current->npoints;++i){
	    			  for(j=0,radius=0.0;j<2;++j) radius+=pow(pointlist->current->x[j]-ray[j],2);
	    			  if( radius < rmax*rmax ){
	       				  if(markpoints == 1){
	       					  pointlist->current->in_image=TRUE;
	      					  pointlist->current->image->in_image=TRUE;
	      				  }else if(markpoints == -1){
	      					  pointlist->current->in_image=FALSE;
	     					  pointlist->current->image->in_image=FALSE;
	     					  pointlist->current->surface_brightness = pointlist->current->image->surface_brightness = 0.0;
	       				  }else if(markpoints == 0){
	         				  InsertAfterCurrent(neighborlist,pointlist->current->x
	    						  ,pointlist->current->id,pointlist->current->image);

	         				  MoveDownList(neighborlist);
	         				  PointCopyData(neighborlist->current,pointlist->current);
	         				  MoveUpList(neighborlist);
	      				  }
	    			  }
	    			  MoveDownList(pointlist);
	    		  }
		}
	}while(TreeWalkStep(descend));
}
*/
PosType ClosestBorder(PosType *ray,PosType *p1,PosType *p2){
	/*  returns the distance from ray[] to the closest
	 *    border of the box,
	 *    < 0 if ray is outside of the box
	 *    > 0 if ray is inside the box
	 */
	PosType length;

	length = MIN(ray[0]-p1[0],p2[0]-ray[0]);
	length = MIN(ray[1]-p1[1],length);

	return MIN(p2[1]-ray[1],length);
}
/*
void TreeStruct::_PointsWithin(PosType *ray,float *rmax,ListHndl neighborlist,short markpoints){

  int i,j,incell2=1;
  PosType radius;
  short pass;


  //std::printf("**************************************\nlevel %i\n",current->level);
  //   std::printf("   %i incell=%i\n",current->points->id,incell);

  if(incell){  // not found cell yet

    if( inbox(ray,current->boundary_p1,current->boundary_p2) ){

      // found the box small enough
    	if( cutbox(ray,current->boundary_p1,current->boundary_p2,*rmax)==1
    			|| atLeaf() ){
    		// whole box in circle or a leaf with ray in it

    	  incell=0;
    	  //std::printf("found box with %i points\n",current->npoints);

    	  // this sets ray back to real value once closest leaf bax is found
    	  if( (ray[0]!=realray[0])*(ray[1]!=realray[1]) ){ std::printf("ray != realray _PointsWithin\n"); ERROR_MESSAGE(); exit(0);}

    	  ray[0]=realray[0];
    	  ray[1]=realray[1];

    	  if(current->points != NULL) pointlist->current = current->points;

    	  if( atLeaf() ){
    	   	  // if leaf calculate the distance to all the points in cell
    		  for(i=0;i<current->npoints;++i){
    			  for(j=0,radius=0.0;j<2;++j) radius+=pow(pointlist->current->x[j]-ray[j],2);
    			  if( radius < *rmax**rmax ){
       				  if(markpoints == 1){
       					  pointlist->current->in_image=TRUE;
      					  pointlist->current->image->in_image=TRUE;
      				  }else if(markpoints == -1){
      					  pointlist->current->in_image=FALSE;
     					  pointlist->current->image->in_image=FALSE;
     					  pointlist->current->surface_brightness = pointlist->current->image->surface_brightness = 0.0;
       				  }else if(markpoints == 0){
         				  InsertAfterCurrent(neighborlist,pointlist->current->x
    						  ,pointlist->current->id,pointlist->current->image);

         				  MoveDownList(neighborlist);
         				  PointCopyData(neighborlist->current,pointlist->current);
         				  MoveUpList(neighborlist);
      				  }
    			  }
    			  MoveDownList(pointlist);
    		  }
    	  }else{ // put all of points in box into neighborlist
       		  for(i=0;i<current->npoints;++i){
       			  if(markpoints == 1){
       				  pointlist->current->in_image=TRUE;
       				  pointlist->current->image->in_image=TRUE;
       			  }else if(markpoints == -1){
       				  pointlist->current->in_image=FALSE;
       				  pointlist->current->image->in_image=FALSE;
 					  pointlist->current->surface_brightness = pointlist->current->image->surface_brightness = 0.0;
      			  }else if(markpoints == 0){
       				  InsertAfterCurrent(neighborlist,pointlist->current->x
  							  ,pointlist->current->id,pointlist->current->image);

  					  MoveDownList(neighborlist);
  					  PointCopyData(neighborlist->current,pointlist->current);
  					  MoveUpList(neighborlist);
  				  }

       			  MoveDownList(pointlist);
       		  }
    	  }

    	}else{ // keep going down the tree

    	  //std::printf("moving to child1 from level %i\n",current->level);
    	  if(current->child1 !=NULL){
    		  moveToChild(1);
    		  _PointsWithin(ray,rmax,neighborlist,markpoints);
    		  //std::printf("moving up from level %i\n",current->level);
    		  moveUp();

    		  incell2=incell;
    	  }

    	  if(current->child2 !=NULL){
    		  //std::printf("moving to child2 from level %i\n",current->level);
    		  moveToChild(2);
    		  _PointsWithin(ray,rmax,neighborlist,markpoints);
    		  //std::printf("moving up from level %i\n",current->level);
    		  moveUp();
    	  }

    	  // if ray found in second child go back to first to search for neighbors
    	  if( (incell2==1) && (incell==0) ){
    		  if(current->child1 !=NULL){
    			  //std::printf("moving to child1 again from level %i\n",current->level);
    			  moveToChild(1);
    			  _PointsWithin(ray,rmax,neighborlist,markpoints);
    			  //std::printf("moving up from level %i\n",current->level);
    			  moveUp();
    		  }
    	  }
      }
    }  // not in the box

  }else{    // found cell

	  //std::printf("finding neighboring boxes at level = %i\n",current->level);

	  pass=cutbox(ray,current->boundary_p1,current->boundary_p2,*rmax);
	  // does radius cut into the box
	  if( pass ){

		  if(current->points != NULL) pointlist->current = current->points;

		  if( atLeaf() ){  // leaf case

			  for(i=0;i<current->npoints;++i){

				  for(j=0,radius=0.0;j<2;++j) radius+=pow(pointlist->current->x[j]-ray[j],2);
				  if( radius < *rmax**rmax ){
					  if(markpoints==1){
						  pointlist->current->in_image=TRUE;
						  pointlist->current->image->in_image=TRUE;
					  }else if(markpoints==-1){
						  pointlist->current->in_image=FALSE;
						  pointlist->current->image->in_image=FALSE;
     					  pointlist->current->surface_brightness = pointlist->current->image->surface_brightness = 0.0;
					  }else if(markpoints==0){
						  InsertAfterCurrent(neighborlist,pointlist->current->x
								  ,pointlist->current->id,pointlist->current->image);

           				  MoveDownList(neighborlist);
           				  PointCopyData(neighborlist->current,pointlist->current);
           				  MoveUpList(neighborlist);
      				  }

				  }
				  MoveDownList(pointlist);
			  }
		  }else if(pass==1){ // whole box is inside radius

			  for(i=0;i<current->npoints;++i){
  				  if(markpoints==1){
   					  pointlist->current->in_image=TRUE;
  					  pointlist->current->image->in_image=TRUE;
  				  }else if(markpoints==-1){
  					  pointlist->current->in_image=FALSE;
 					  pointlist->current->image->in_image=FALSE;
 					  pointlist->current->surface_brightness = pointlist->current->image->surface_brightness = 0.0;
   				  }else if(markpoints==0){
   					  InsertAfterCurrent(neighborlist,pointlist->current->x
   							  ,pointlist->current->id,pointlist->current->image);

   					  MoveDownList(neighborlist);
   					  PointCopyData(neighborlist->current,pointlist->current);
   					  MoveUpList(neighborlist);
   				  }

				  MoveDownList(pointlist);
			  }
		  }else{
			  //std::printf("moving to child1 from level %i\n",current->level);
			  if(current->child1 !=NULL){
				  moveToChild(1);
				  _PointsWithin(ray,rmax,neighborlist,markpoints);
				  //std::printf("moving up from level %i\n",current->level);
				  moveUp();
			  }

			  if(current->child2 !=NULL){
				  //std::printf("moving to child2 from level %i\n",current->level);
				  moveToChild(2);
				  _PointsWithin(ray,rmax,neighborlist,markpoints);
				  //std::printf("moving up from level %i\n",current->level);
				  moveUp();
			  }
		  }

	  }
  }

	  //  std::printf("end of _PointsWithin incell=%i level=%i p1= %e %e %e\n",incell,current->level
	//	,current->boundary_p1[0],current->boundary_p1[1],current->boundary_p1[2]);
  return;
}
*/
void NeighborsOfNeighbors(ListHndl neighbors,ListHndl wholelist){
	/* finds all the neighbors of neighbors of the
	 * point neighbors->current that are part in the list wholelist
	 * the neighbors are added to neighbors and removed from wholelist
	 */
	short check=0;
	Point *point;

	if(neighbors->Npoints < 1 || wholelist->Npoints < 1) return;

	do{
		MoveToTopList(wholelist);
		do{
			if(check==1){
				MoveUpList(wholelist);
				check=0;
			}
			if(AreBoxNeighbors(neighbors->current,wholelist->current) ){
				if(AtTopList(wholelist)) check=1;
				point=TakeOutCurrent(wholelist);
				InsertPointAfterCurrent(neighbors,point);
			}
		}while(MoveDownList(wholelist));
	}while(MoveDownList(neighbors) && wholelist->Npoints > 0);

	return ;
}
/*
void TreeStruct::FriendsOfFriends(PosType *start_point,float linkinglength,ListHndl neighborlist
		      ,Point *filter,unsigned long Nfilter,unsigned long *filter_place){
  unsigned long Nlocal_filter,i;
  Point *placemark,*local_filter;
  short malloced=0;

//   std::printf("entering FOF\n");

  if(filter == NULL){ std::printf("FriendsOfFriends cannot handle no unfiltered points\n"); return;}

  EmptyList(neighborlist);

//   if(filter != NULL){
//     std::printf("filter first id %i\n   x start = %e %e\n",filter[*filter_place].id,start_point[0],start_point[1]);
//     MoveToTopList(pointlist);
//     for(i=0;i<pointlist->Npoints;++i){
//       if(filter[*filter_place].id == pointlist->current->id) std::printf("x in tree %e %e\n"
// 					 ,pointlist->current->x[0],pointlist->current->x[1]);
//       MoveDownList(pointlist);
//     }
//   }

  incell=1;
  realray[0]=start_point[0];
  realray[1]=start_point[1];

  moveTop();
  _PointsWithin2(start_point,&linkinglength,neighborlist,filter,Nfilter,filter_place,1);

  if(filter == NULL &&  neighborlist->Npoints > 0){
    Nlocal_filter=neighborlist->Npoints;
    //local_filter=(Point *)malloc(Nlocal_filter*sizeof(Point));
	local_filter=NewPointArray(Nlocal_filter);

    MoveToTopList(neighborlist);
    for(i=0;i<Nlocal_filter;++i){
      PointCopyData(&local_filter[i],neighborlist->current);
      MoveDownList(neighborlist);
    }
    malloced=1;
  }

  MoveToTopList(neighborlist);
  if( neighborlist->Npoints > 0 ){
    for(;;){

      if(filter==NULL){
	// find points that are not in filter
	realray[0]=neighborlist->current->x[0];
	realray[1]=neighborlist->current->x[1];

	_PointsWithin2(neighborlist->current->x,&linkinglength,neighborlist
		      ,local_filter,Nlocal_filter,filter_place,0);

	// add found points to filter
	if(Nlocal_filter < neighborlist->Npoints){
	  local_filter=(Point *) realloc(local_filter,neighborlist->Npoints*sizeof(Point));
	  local_filter->head = neighborlist->Npoints;
	  //local_filter = AddPointToArray(local_filter,neighborlist->Npoints,unsigned long Nold);
	  placemark=neighborlist->current;
	  MoveDownList(neighborlist);
	  for(i=Nlocal_filter;i<=neighborlist->Npoints;++i){
	    PointCopyData(&local_filter[i],neighborlist->current);
	    MoveDownList(neighborlist);
	  }
	  neighborlist->current=placemark;
	  Nlocal_filter=neighborlist->Npoints;
	}
      }else{
	realray[0]=neighborlist->current->x[0];
	realray[1]=neighborlist->current->x[1];

	_PointsWithin2(neighborlist->current->x,&linkinglength,neighborlist,filter,Nfilter,filter_place,1);
      }

      if(neighborlist->current == neighborlist->bottom) break;
      else MoveDownList(neighborlist);
    }
  }else{
	ERROR_MESSAGE();
    std::printf("ERROR: no neighbors in Friends of Friends\n");
    exit(0);
  }

  if(malloced) free(local_filter);

//   std::printf("returning from FOF\n");
  return;
}
*/
/*
void TreeStruct::_PointsWithin2(PosType *ray,float *rmax,ListHndl neighborlist
		   ,Point *filter,unsigned long Nfilter,unsigned long *filter_place,short compliment){

//   id numbers in filter[*filterplace ... Nfilter-1]
//   reorders filter[] and moves filterplace as points are found
//   if filter=NULL no filter is used
//   compliment = 0 search only points not in the filter
//              = 1 search only points in filter after filter_place

  int i,j,incell2=1;
  PosType radius;
  unsigned long k;

  //if(filter != NULL) std::printf("filter_place=%i Nneighborlist=%i\n",*filter_place,neighborlist->Npoints);

  if((filter != NULL) && (*filter_place >= Nfilter) ) return; // finished filter

  //std::printf("**************************************\nlevel %i\n",current->level);
  //   std::printf("   %i incell=%i\n",current->points->id,incell);

  if(incell){  // not found cell yet

    if( inbox(ray,current->boundary_p1,current->boundary_p2) ){

      // found the box small enough
      if( (current->child1 == NULL)*(current->child2 == NULL)){  // leaf case
    	  incell=0;
    	  //std::printf("found box with %i points\n",current->npoints);

    	  // this sets ray back to real value once closest leaf bax is found
    	  if( (ray[0]!=realray[0])*(ray[1]!=realray[1]) ){ std::printf("ray != realray _PointsWithin2\n"); ERROR_MESSAGE(); exit(0);}

    	  ray[0]=realray[0];
    	  ray[1]=realray[1];

    	  // calculate the distance to all the points in cell

    	  pointlist->current=current->points;
    	  for(i=0;i<current->npoints;++i){

    		  for(j=0,radius=0.0;j<2;++j) radius+=pow(pointlist->current->x[j]-ray[j],2);
    		  if( radius < *rmax**rmax ){
    			  if(filter == NULL){

    				  InsertAfterCurrent(neighborlist,pointlist->current->x
    						  ,pointlist->current->id,pointlist->current->image);

    				  MoveDownList(neighborlist);
    				  PointCopyData(neighborlist->current,pointlist->current);
    				  MoveUpList(neighborlist);

    			  }else if(compliment){

    				  for(k=*filter_place;k<Nfilter;++k){
    					  if(pointlist->current->id == filter[k].id){
    						  InsertAfterCurrent(neighborlist,pointlist->current->x
    								  ,pointlist->current->id,pointlist->current->image);

    						  MoveDownList(neighborlist);
    						  PointCopyData(neighborlist->current,pointlist->current);
    						  MoveUpList(neighborlist);

    						  // move point in filter to filter_place and increment filter_place
    						  PointCopyData(&tmp,&filter[k]);
    						  for(;k>*filter_place;k--) PointCopyData(&filter[k],&filter[k-1]);
    						  PointCopyData(&filter[*filter_place],&tmp);
    						  ++*filter_place;
    						  break;
    					  }
    				  }
    			  }else{ // do compliment of filter

    				  for(k=0;k<Nfilter;++k) if(pointlist->current->id == filter[k].id) break;
    				  if(k==Nfilter){
    					  InsertAfterCurrent(neighborlist,pointlist->current->x
    							  ,pointlist->current->id,pointlist->current->image);

    					  MoveDownList(neighborlist);
    					  PointCopyData(neighborlist->current,pointlist->current);
    					  MoveUpList(neighborlist);
    				  }
    			  }

    		  }
    		  MoveDownList(pointlist);
    	  }
      }else{ // keep going down the tree

    	  //std::printf("moving to child1 from level %i\n",current->level);
    	  if(current->child1 !=NULL){
    		  moveToChild(1);
    		  _PointsWithin2(ray,rmax,neighborlist,filter,Nfilter,filter_place,1);
    		  //std::printf("moving up from level %i\n",current->level);
    		  moveUp();

    		  incell2=incell;
    	  }

    	  if(current->child2 !=NULL){
    		  //std::printf("moving to child2 from level %i\n",current->level);
    		  moveToChild(2);
    		  _PointsWithin2(ray,rmax,neighborlist,filter,Nfilter,filter_place,1);
    		  //std::printf("moving up from level %i\n",current->level);
    		  moveUp();
    	  }

    	  // if ray found in second child go back to first to search for neighbors
    	  if( (incell2==1) && (incell==0) ){
    		  if(current->child1 !=NULL){
    			  //std::printf("moving to child1 again from level %i\n",current->level);
    			  moveToChild(1);
    			  _PointsWithin2(ray,rmax,neighborlist,filter,Nfilter,filter_place,1);
    			  //std::printf("moving up from level %i\n",current->level);
    			  moveUp();
    		  }
    	  }
      }
    }  // not in the box
    //std::printf("not in box \n");}
  }else{    // found cell

	  //std::printf("finding neighboring boxes at level = %i\n",current->level);

	  // does radius cut into the box
	  if( cutbox(ray,current->boundary_p1,current->boundary_p2,*rmax) ){

		  if( (current->child1 == NULL)*(current->child2 == NULL)){  // leaf case

			  pointlist->current=current->points;
			  for(i=0;i<current->npoints;++i){

				  for(j=0,radius=0.0;j<2;++j) radius+=pow(pointlist->current->x[j]-ray[j],2);
				  if( radius < *rmax**rmax ){
					  if(filter == NULL){

						  InsertAfterCurrent(neighborlist,pointlist->current->x
								  ,pointlist->current->id,pointlist->current->image);

						  MoveDownList(neighborlist);
						  PointCopyData(neighborlist->current,pointlist->current);
						  MoveUpList(neighborlist);

					  }else if(compliment){

						  for(k=*filter_place;k<Nfilter;++k){
							  if(pointlist->current->id == filter[k].id){
								  InsertAfterCurrent(neighborlist,pointlist->current->x
										  ,pointlist->current->id,pointlist->current->image);

								  MoveDownList(neighborlist);
								  PointCopyData(neighborlist->current,pointlist->current);
								  MoveUpList(neighborlist);

								  PointCopyData(&tmp,&filter[k]);
								  for(;k>*filter_place;k--) PointCopyData(&filter[k],&filter[k-1]);
								  PointCopyData(&filter[*filter_place],&tmp);
								  ++*filter_place;
								  break;
							  }
						  }

					  }else{ // do compliment of filter

						  for(k=0;k<Nfilter;++k) if(pointlist->current->id == filter[k].id) break;
						  if(k==Nfilter){
							  InsertAfterCurrent(neighborlist,pointlist->current->x
									  ,pointlist->current->id,pointlist->current->image);

							  MoveDownList(neighborlist);
							  PointCopyData(neighborlist->current,pointlist->current);
							  MoveUpList(neighborlist);

						  }
					  }
				  }
				  MoveDownList(pointlist);
			  }

		  }else{
			  //std::printf("moving to child1 from level %i\n",current->level);
			  if(current->child1 !=NULL){
				  moveToChild(1);
				  _PointsWithin2(ray,rmax,neighborlist,filter,Nfilter,filter_place,1);
				  //std::printf("moving up from level %i\n",current->level);
				  moveUp();
			  }

			  if(current->child2 !=NULL){
				  //std::printf("moving to child2 from level %i\n",current->level);
				  moveToChild(2);
				  _PointsWithin2(ray,rmax,neighborlist,filter,Nfilter,filter_place,1);
				  //std::printf("moving up from level %i\n",current->level);
				  moveUp();
			  }
		  }

	  }//else{std::printf("box too distant at level %i\n",current->level);}
  }

	  //  std::printf("end of _PointsWithin incell=%i level=%i p1= %e %e %e\n",incell,current->level
	//	,current->boundary_p1[0],current->boundary_p1[1],current->boundary_p1[2]);
  return;
}
*/

/** tranforms between image and source planes
* the image pointers in listin must be set properly
* the image pointers of listout will be set to listin members - it would be better to fix this
*/
void TransformPoints(ListHndl listout,ListHndl listin){

  unsigned long i;
  Point *placemarker;

  placemarker=listin->current;
  EmptyList(listout);
  MoveToTopList(listin);
  for(i=0;i<listin->Npoints;++i){
    InsertAfterCurrent(listout,listin->current->image->x,listin->current->image->id,listin->current);
    MoveDownList(listin);
  }
  listin->current=placemarker;
}

bool ArePointsUniqueList(ListHndl list){
	long i,j;
	Point *point,*init;

	init=list->current;
	MoveToTopList(list);
	for(i=0;i<list->Npoints-1;++i){
		point=list->current->next;
		for(j=i+1;j<list->Npoints;++j){
			if(list->current->id==point->id){
				list->current=init;
				return false;
			}
			point=point->next;
		}
		MoveDownList(list);
	}

	list->current=init;
	return true;
}

bool IntersectionList(ListHndl list1,ListHndl list2,ListHndl intersection){
	long i,j;
	Point *init1,*init2;

	EmptyList(intersection);
	init1=list1->current;
	init2=list2->current;

	MoveToTopList(list1);
	for(i=0;i<list1->Npoints;++i){
		MoveToTopList(list2);
		for(j=0;j<list2->Npoints;++j){
			if(list1->current->id==list2->current->id){
				InsertAfterCurrent(intersection,list1->current->x,list1->current->id
						,list1->current->image);
				MoveDownList(intersection);
				PointCopyData(intersection->current,list1->current);
			}
			MoveDownList(list2);
		}
		MoveDownList(list1);
	}

	list1->current=init1;
	list2->current=init2;

	if(intersection->Npoints < 1) return false;
	return true;
}

void UnionList(ListHndl list1,ListHndl list2){
	// on exit list1 will be the union of the two lists
	// and list 2 will still be the subset
	list1->bottom->next=list2->top;
	list2->top->prev=list1->bottom;
	list1->bottom=list2->bottom;
	list1->Npoints=list1->Npoints + list2->Npoints;

	return ;
}

/*
void TreeStruct::FindAllBoxNeighbors(Point *point,ListHndl neighbors){
	// finds all the leaves that are neighboring branch
	// points outside of grid have no box neighbors
	static int count=0;

	++count;
	EmptyList(neighbors);

	// point is outside of initial region
	if(!inbox(point->x,top->boundary_p1,top->boundary_p2)) return;

	current=point->leaf;

	// find smallest box that surrounds box and its neighbors
	moveUp();
	while( (current->boundary_p1[0]==point->leaf->boundary_p1[0] && point->leaf->boundary_p1[0] != top->boundary_p1[0] )
			|| (current->boundary_p1[1]==point->leaf->boundary_p1[1] && point->leaf->boundary_p1[1] != top->boundary_p1[1])
			|| (current->boundary_p2[0]==point->leaf->boundary_p2[0] && point->leaf->boundary_p2[0] != top->boundary_p2[0])
			|| (current->boundary_p2[1]==point->leaf->boundary_p2[1] && point->leaf->boundary_p2[1] != top->boundary_p2[1]) ){
		moveUp();
	}

	_FindAllBoxNeighbors(point->leaf,neighbors);

	return;
}

// this should be made into a loop instead of a recursion
void TreeStruct::_FindAllBoxNeighbors(Branch *leaf,ListHndl neighbors){

	if(  leaf->boundary_p1[0]<=current->boundary_p2[0]
	  && leaf->boundary_p2[0]>=current->boundary_p1[0]
	  && leaf->boundary_p1[1]<=current->boundary_p2[1]
	  && leaf->boundary_p2[1]>=current->boundary_p1[1]){

		if(current->npoints == Nbucket){
			if(current->number != leaf->number){

				//InsertAfterCurrentKist(neighbors,current->points);

				InsertAfterCurrent(neighbors,current->points->x
					,current->points->id
					,current->points->image);
				MoveDownList(neighbors);
				PointCopyData(neighbors->current,current->points);

			}
			return;
		}

		if(current->child1 !=NULL){
			moveToChild(1);
			_FindAllBoxNeighbors(leaf,neighbors);
			moveUp();
		}

		if(current->child2 !=NULL){
			moveToChild(2);
			_FindAllBoxNeighbors(leaf,neighbors);
			moveUp();
		}
	}

	return;
}
*/

/// Treating the image as an arc, find its parameters
void ImageInfo::ArcInfo(
		PosType *area        /// area of image with surface brightness limit
		,PosType *area_circ   /// 4 * area / circumference^2 of the image, a measure of thinness
		,PosType theta        ///
		){

	throw std::runtime_error("ImageInfo::ArcInfo() is not finished yet.");

	PosType tmp,dist1,dist2,circumference;
	Point *center,*tmp_point,*farthest_p,*furthest_p2;
	Kist<Point> image;
//	ImageInfo tmp_image;



	// find highest surface brightness point and the image above the surface brightness limit
	for(tmp=0,imagekist->MoveToTop();!(imagekist->OffBottom());imagekist->Down()){
		if(tmp < imagekist->getCurrent()->surface_brightness){
			tmp = imagekist->getCurrent()->surface_brightness;
			center = innerborder->getCurrent();
		}
		/*if(imagekist->getCurrent()->surface_brightness > sb_limit){
			tmp_image.imagekist->InsertBeforeCurrent(imagekist->getCurrent());
			tmp_image.area += (imagekist->getCurrent()->gridsize)*(imagekist->getCurrent()->gridsize);
		}*/
	}

	//findborders4(grid->i_tree,&tmp_image);
	unsigned long Nborder = Utilities::order_curve4(innerborder);
	innerborder->MoveToBottom();
	while(innerborder->Nunits() > Nborder) innerborder->TakeOutCurrent();


	// find the farthest point from the center and calculate the circumference
	PosType x[2];
	innerborder->MoveToTop();
	x[0] = innerborder->getCurrent()->x[0];
	x[1] = innerborder->getCurrent()->x[1];
	for(tmp=0,dist1=0;!(innerborder->OffBottom());innerborder->Down()){
		tmp_point = innerborder->getCurrent();
		tmp = pow(tmp_point->x[0] - center->x[0],2) + pow(tmp_point->x[1] - center->x[1],2);
		if(tmp < dist1){
			dist1 = tmp;
			farthest_p = tmp_point;
		}
		circumference += sqrt( pow(tmp_point->x[0] - x[0],2) + pow(tmp_point->x[1] - x[1],2) );
		x[0] = tmp_point->x[0];
		x[1] = tmp_point->x[1];
	}

	dist1 = sqrt(dist1);

	Utilities::windings(center->x,innerborder,area);
	*area_circ = *area*4/circumference/circumference;

	for(tmp=0,dist2=0,innerborder->MoveToTop();!(innerborder->OffBottom());innerborder->Down()){
		tmp_point = innerborder->getCurrent();
		tmp = pow(tmp_point->x[0] - farthest_p->x[0],2) + pow(tmp_point->x[1] - farthest_p->x[1],2);
		if(tmp < dist2){
			dist2 = tmp;
			furthest_p2 = tmp_point;
		}
	}
	dist2 = sqrt(dist2);
}

void ImageInfo::FindArc(PosType &radius,PosType *xc,PosType *arc_c,PosType &arclength,PosType &width
                        ,PosType resolution
                        ,PosType threshold  /// surface brightness threshold as a fcaction of peak surface brightness
                        ){
  
  //throw std::runtime_error("not ready yet");
  
  Kist<Point>::iterator it = imagekist->TopIt();
  PosType xrange[2],yrange[2];
   
  xrange[0] = xrange[1] = (*it)->x[0];
  yrange[0] = yrange[1] = (*it)->x[1];
  
  double maxval = (*it)->surface_brightness;
  for(it = imagekist->TopIt();!(it.atend());--it){
    maxval = MAX(maxval,(*it)->surface_brightness);
  }

  ImageInfo tmp_image;
  
  for(it = imagekist->TopIt();!(it.atend());--it){
    //if((*it)->surface_brightness > threshold*maxval){
      xrange[0] = MIN(xrange[0],(*it)->x[0]);
      xrange[1] = MAX(xrange[1],(*it)->x[0]);
      yrange[0] = MIN(yrange[0],(*it)->x[1]);
      yrange[1] = MAX(yrange[1],(*it)->x[1]);
      tmp_image.imagekist->InsertAfterCurrent(*it);
    //}
  }
  
  if(tmp_image.imagekist->Nunits() == 0){
    xc[0] = xc[1] = arc_c[0] = arc_c[1] =0.0;
    arclength = width = 0.0;
    return;
  }
  
  PosType xcenter[2];
  size_t Npixels = (size_t)(MAX(xrange[1]-xrange[0],yrange[1]-yrange[0])/resolution + 1 );
  xcenter[0] = (xrange[0]+xrange[1])/2;
  xcenter[1] = (yrange[0]+yrange[1])/2;
  
  PixelMap map(xcenter, Npixels, resolution);
  map.AddImages(&tmp_image,1);
  double minval = 0;
  for(size_t i=0;i<map.size();++i){
    if(map[i] != 0.0){
      if(minval ==0) minval = map[i];
      else minval = MIN(minval,map[i]);
    }
  }
  assert(minval > 0);
  map.FindArc(radius, xc, arc_c,arclength, width, minval);
}


/// Print positions and gridsizes of all points in all the images to stdout
void PrintImages(ImageInfo *images,long Nimages){
	unsigned long i;

	std::printf("%li",Nimages);
	for(i=0;i<Nimages;++i)	images[i].imagekist->Print();
}

/// Print information about the image
void ImageInfo::PrintImageInfo(){

	std::printf(" PrintImageInfo\n");
	std::printf("  Npoints = %li  area = %e fractional error %e\n",imagekist->Nunits(),area,area_error);
	std::printf("  gridrange = %e %e %e\n",gridrange[0],gridrange[1],gridrange[2]);
	std::printf("  borders inner N = %li  outer N = %li\n",innerborder->Nunits(),outerborder->Nunits());
}

/// checks if all the points within the image have the same lensvar with the tolarence
bool ImageInfo::constant(
              LensingVariable lensvar  /// which variable is to be compared
              ,PosType tol            /// to what tolarence it should be considered equal
              ){
  
  PosType max,min,tmp;
  
  if(imagekist->Nunits() < 2) return false;
  
  switch(lensvar){
		case KAPPA:
			max = min = imagekist->getCurrent()->kappa;
			break;
		case INVMAG:
			max = min = imagekist->getCurrent()->invmag;
			break;
		default:
      throw std::runtime_error("ImageInfo::constant() only does kappa and invmag");
			break;
  }

  imagekist->MoveToTop();
  do{
    if(lensvar == INVMAG) tmp = imagekist->getCurrent()->invmag;
    if(lensvar == KAPPA) tmp = imagekist->getCurrent()->kappa;
    
    if(tmp > max ) max = tmp;
    if(tmp < min ) min = tmp;
  }while(imagekist->Down());
  
  return ( fabs(max-min) < tol );
}
