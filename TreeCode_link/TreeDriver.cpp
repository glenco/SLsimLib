
#include <slsimlib.h>

/* median_cut determines how the cells are subdivided */
/*    if ==0  equal volume cuts, Warning this option causes an error*/
/*    if ==1  median point cuts */
static int incell;
static Point **temp_points,tmp;
static double realray[2];
//Point *point_global;

Point *NearestNeighbor(TreeHndl tree,double *ray,int Nneighbors,ListHndl neighborlist,short direction){
  /* nearest neighbor points in a given direction, direction != 0 should not */
  /*    be used except on a grid */
  /* direction = 0 distance */
  /*           = 1 to right */
  /*           = 2 to left */
  /*           = 3 to up */
  /*           = 4 to down */
  unsigned long i;
  void _NearestNeighbor(TreeHndl tree,double *ray,int Nneighbors
			,Point **neighborpoints,double *rneighbor,short *direction);
  static int count=0,oldNneighbors=-1;
  static double *rneighbors;
  static Point **neighborpoints;

  /*std::printf("entering NN\n");*/

  if(tree->top->npoints <= Nneighbors){
	  ERROR_MESSAGE();
	  std::printf("ERROR: in NearestNeighbor, number of neighbors > total number of points\n");
	  exit(1);
  }

  if(Nneighbors <= 0){
	  EmptyList(neighborlist);
	  return NULL;
  }

  if(count==0){
    /*std::printf("allocating memory\n");*/
    rneighbors=(double *)malloc((Nneighbors+tree->Nbucket)*sizeof(double));
    assert(rneighbors);
    neighborpoints=(Point **)malloc((Nneighbors+tree->Nbucket)*sizeof(Point *));
    assert(neighborpoints);
    temp_points=(Point **)malloc((Nneighbors+tree->Nbucket)*sizeof(Point *));
    assert(temp_points);
    ++count;
    oldNneighbors=Nneighbors;

  }else if(oldNneighbors < Nneighbors){ /* if the number of nearest neighbors goes up get more mem */
    /*std::printf("re-allocating memory\n");*/
    rneighbors=(double *)realloc(rneighbors,(Nneighbors+tree->Nbucket)*sizeof(double));
    neighborpoints=(Point **)realloc(neighborpoints,(Nneighbors+tree->Nbucket)*sizeof(Point *));
    temp_points=(Point **)realloc(temp_points,(Nneighbors+tree->Nbucket)*sizeof(Point *));
    oldNneighbors=Nneighbors;
  }


  /*   std::printf("Nneighbors=%i\n",Nneighbors); */
  /*   std::printf("array sizes=%i\n",Nneighbors+tree->Nbucket); */

  /* initalize distance to neighbors to a large number */
  for(i=0;i<(tree->Nbucket+Nneighbors);++i){
    rneighbors[i]=10*(tree->top->boundary_p2[0]-tree->top->boundary_p1[0]);
  }

  moveTop(tree);
  //   std::printf("p1= [%f,%f]\n", tree->current->boundary_p1[0],tree->current->boundary_p1[1]);
  //   std::printf("p2= [%f,%f]\n", tree->current->boundary_p2[0],tree->current->boundary_p2[1]);

  realray[0]=ray[0];
  realray[1]=ray[1];

  if( inbox(ray,tree->current->boundary_p1,tree->current->boundary_p2) == 0 ){
    std::printf("Warning: in NearestNeighbor, ray is not inside the simulation box\n    should work in any case\n      ray= %e %e\n",ray[0],ray[1]);

    ray[0]=DMAX(ray[0],tree->current->boundary_p1[0]);
    ray[0]=DMIN(ray[0],tree->current->boundary_p2[0]);

    ray[1]=DMAX(ray[1],tree->current->boundary_p1[1]);
    ray[1]=DMIN(ray[1],tree->current->boundary_p2[1]);
  }
  incell=1;

  if(direction==0) EmptyList(neighborlist);
  if(direction==-1) direction=0;

  _NearestNeighbor(tree,ray,Nneighbors,neighborpoints,rneighbors,&direction);

  /* convert from point array to exported point list */

  for(i=0;i<Nneighbors;++i){
    InsertAfterCurrent(neighborlist,neighborpoints[i]->x,neighborpoints[i]->id,neighborpoints[i]->image);
    MoveDownList(neighborlist);
    PointCopyData(neighborlist->current,neighborpoints[i]);/**/
  }

  return neighborpoints[0];
}

void _NearestNeighbor(TreeHndl tree,double *ray,int Nneighbors,Point **neighborpoints,double *rneighbors,short *direction){

  int i,incell2=1;
  unsigned long index[Nneighbors+tree->Nbucket];
  double dx,dy;

  //std::printf("**************************************\nlevel %i\n",tree->current->level);
  //for(i=0;i<tree->current->npoints;++i) std::printf("   %i\n",tree->current->points[i]);

  if(incell){  /* not found cell yet */

    if( inbox(ray,tree->current->boundary_p1,tree->current->boundary_p2) ){

    	// found the box small enough */
    	if( tree->current->npoints <= (Nneighbors+tree->Nbucket) ){
    		incell=0;
    		/*std::printf("found box with %i points\n",tree->current->npoints);*/

    		/* this sets ray back to real value once closest leaf bax is found */
    		if( (ray[0]!=realray[0])*(ray[1]!=realray[1]) ){ std::printf("ray != realray _NearestNeighbor\n"); exit(0);}

    		ray[0]=realray[0];
    		ray[1]=realray[1];

    		/* calculate the distance to all the points in cell */
    		if(tree->current->points != NULL ) tree->pointlist->current=tree->current->points;
    		for(i=0;i<tree->current->npoints;++i){

    			dx=tree->pointlist->current->x[0]-ray[0];
    			dy=tree->pointlist->current->x[1]-ray[1];

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
    			temp_points[i]=tree->pointlist->current;
    			MoveDownList(tree->pointlist);
    		}

    		if(tree->current->npoints > 0){
    			double_sort(tree->current->npoints,rneighbors-1,index-1);
    			for(i=0;i<tree->current->npoints;++i) neighborpoints[i] = temp_points[index[i]];
    		}

      }else{  // keep going down the tree

    	  if(tree->current->child1 !=NULL){
    		  moveToChild(tree,1);
    		  _NearestNeighbor(tree,ray,Nneighbors,neighborpoints,rneighbors,direction);
    		  moveUp(tree);

    		  incell2=incell;
    	  }

    	  if(tree->current->child2 !=NULL){
    		  moveToChild(tree,2);
    		  _NearestNeighbor(tree,ray,Nneighbors,neighborpoints,rneighbors,direction);
    		  moveUp(tree);
    	  }

    	  // if ray found in second child go back to first to search for neighbors
    	  if( (incell2==1) && (incell==0) ){
    		  if(tree->current->child1 !=NULL){
    			  moveToChild(tree,1);
    			  _NearestNeighbor(tree,ray,Nneighbors,neighborpoints,rneighbors,direction);
    			  moveUp(tree);
    		  }
    	  }

      }
    } // not in the box

  }else{  /* already found cell */

	  // does radius cut into the box
	  if( cutbox(ray,tree->current->boundary_p1,tree->current->boundary_p2,rneighbors[Nneighbors-1]) ){

		  if( atLeaf(tree) ){  /* leaf case */

			  /* combine found neighbors with points in box and resort */
			  if(tree->current->points != NULL) tree->pointlist->current=tree->current->points;

			  for(i=0;i<Nneighbors;++i){
				  index[i]=i;
				  temp_points[i]=neighborpoints[i];
			  }

			  for(i=Nneighbors;i<(tree->current->npoints+Nneighbors);++i){

				  dx=tree->pointlist->current->x[0]-ray[0];
				  dy=tree->pointlist->current->x[1]-ray[1];

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
				  temp_points[i]=tree->pointlist->current;
				  MoveDownList(tree->pointlist);
			  }

			  // sort the points found so far
			  if(tree->current->npoints > 0){
				  double_sort(tree->current->npoints+Nneighbors,rneighbors-1,index-1);
				  for(i=0;i<(tree->current->npoints+Nneighbors);++i) neighborpoints[i] = temp_points[index[i]];
			  }

		  }else{

			  if(tree->current->child1 !=NULL){
				  moveToChild(tree,1);
				  _NearestNeighbor(tree,ray,Nneighbors,neighborpoints,rneighbors,direction);
				  moveUp(tree);
			  }

			  if(tree->current->child2 !=NULL){
				  moveToChild(tree,2);
				  _NearestNeighbor(tree,ray,Nneighbors,neighborpoints,rneighbors,direction);
				  /*std::printf("moving up from level %i\n",tree->current->level);*/
				  moveUp(tree);
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
double BoxIntersection(Branch *branch1,Branch *branch2){
	// returns area of intersection between two branches
	double area=0;

	area = MIN(branch1->boundary_p2[0],branch2->boundary_p2[0])
	     - MAX(branch1->boundary_p1[0],branch2->boundary_p1[0]);
	if(area < 0) return 0.0;

	area *= MIN(branch1->boundary_p2[1],branch2->boundary_p2[1])
	      - MAX(branch1->boundary_p1[1],branch2->boundary_p1[1]);
	if(area < 0) return 0.0;

	return area;
}

// return 1 (0) if box is (not) within rmax of ray
int cutbox(double *ray,double *p1,double *p2,double rmax){
	/*  returns:  0 if whole box is outside rmax from ray[]
	 *            1 if whole box is inside circle but ray is not in the box
	 *            2 if ray[] is inside box
	 *            3 if box intersects circle but ray[] is not inside box
	 */
  short i,tick=0;
  double close[2],rtmp;
  
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
bool CircleInBox(double *ray,double radius,double *p1,double *p2){

	if(!inbox(ray,p1,p2)) return false;

	double close[2];
	int i;

	if((ray[0] + radius) > p2[0] ) return false;
	if((ray[0] - radius) < p1[0] ) return false;
	if((ray[1] + radius) > p2[1] ) return false;
	if((ray[1] - radius) < p1[1] ) return false;

	return true;
}

bool BoxInCircle(double *ray,double radius,double *p1,double *p2){

	double rad2 = radius*radius;

	if((pow(p1[0] - ray[0],2) + pow(p1[1] - ray[1],2)) > rad2) return false;
	if((pow(p2[0] - ray[0],2) + pow(p2[1] - ray[1],2)) > rad2) return false;
	if((pow(p1[0] - ray[0],2) + pow(p2[1] - ray[1],2)) > rad2) return false;
	if((pow(p2[0] - ray[0],2) + pow(p1[1] - ray[1],2)) > rad2) return false;

	return true;
}

/**
 *   Finds the leaf the ray is in and adds Nadd to all of is parent leaves
*/

void _FindLeaf(TreeHndl tree,double *ray,unsigned long Nadd){

	bool contin;

	assert(inbox(ray,tree->current->boundary_p1,tree->current->boundary_p2) );
	do{
		tree->current->npoints += Nadd;

		/* change center of mass */
		tree->current->center[0]=(tree->current->center[0]*(tree->current->npoints-1) + ray[0])/tree->current->npoints;
		tree->current->center[1]=(tree->current->center[1]*(tree->current->npoints-1) + ray[1])/tree->current->npoints;

		contin=false;
		if(tree->current->child1 != NULL &&
				inbox(ray,tree->current->child1->boundary_p1,tree->current->child1->boundary_p2) ){
			moveToChild(tree,1);
			contin=true;
		}else if(tree->current->child2 != NULL &&
				inbox(ray,tree->current->child2->boundary_p1,tree->current->child2->boundary_p2) ){
			moveToChild(tree,2);
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

void FindBoxPoint(TreeHndl tree,double *ray,Point *point){
	// return the point that is in the same box as ray[2]
	// if Nbuck > 1 the head of the point array is returned
	void _FindBox(TreeHndl tree,double *ray);
	Branch *branch;

	//point=(Point *)malloc(sizeof(Point));
	//moveTop(tree);
	//printTree(tree);

	branch=tree->current;
	moveTop(tree);
	   // check if ray is outside initial box
	if( inbox(ray,tree->current->boundary_p1,tree->current->boundary_p2) == 0 ){
		std::printf("FindBox: ray outside of grided range\n");
		return;
	}

	_FindBox(tree,ray);
	PointCopyData(point,tree->current->points);

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

	tree->current=branch;
	//return point;
}

void _FindBox(TreeHndl tree,double *ray){
	bool contin;

	// replaced recursion with iteration
	do{
		contin=false;
		if(tree->current->child1 !=NULL &&
				inbox(ray,tree->current->child1->boundary_p1,tree->current->child1->boundary_p2) ){
			moveToChild(tree,1);
			contin=true;
		}
		if(tree->current->child2 !=NULL &&
				inbox(ray,tree->current->child2->boundary_p1,tree->current->child2->boundary_p2) ){
			moveToChild(tree,2);
			contin=true;
		}
	}while(contin);

	return;
}


Point *sortList(long n, double *arr,ListHndl list,Point *firstpointin){
	/** simple sort for points in linked list **/
	/** slow for > 20 points **/
  long i,j;
  double a;
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


void PointsWithin(TreeHndl tree,double *ray,float rmax,ListHndl neighborlist,short markpoints){
	void _PointsWithin(TreeHndl tree,double *ray,float *rmax,ListHndl neighborlist,short markpoints);
/*
 *  finds all points in tree that lie within rmax of the point ray[]
 *   markpoints = 0  does not change in_image variable in any point, gives a list of neighbors
 *              = 1  makes in_image=true for all points in image, gives no list of neighbors
 *              = -1 makes in_image=false for all points in image to reset, gives no list of neighbors
 */

  if(markpoints==0) EmptyList(neighborlist);
  
  realray[0]=ray[0];
  realray[1]=ray[1];

  moveTop(tree);
  if( inbox(ray,tree->current->boundary_p1,tree->current->boundary_p2) == 0 ){
    std::printf("Warning: in PointsWithin, ray is not inside the simulation box\n    should work in any case\n      ray= %e %e\n     boundary p1 = %e %e p2 = %e %e\n",ray[0],ray[1]
	   ,tree->current->boundary_p1[0],tree->current->boundary_p1[1]
	   ,tree->current->boundary_p2[0],tree->current->boundary_p2[1]);

    ray[0]=DMAX(ray[0],tree->current->boundary_p1[0]);
    ray[0]=DMIN(ray[0],tree->current->boundary_p2[0]);

    ray[1]=DMAX(ray[1],tree->current->boundary_p1[1]);
    ray[1]=DMIN(ray[1],tree->current->boundary_p2[1]);
  }
  incell=1;

  _PointsWithin(tree,ray,&rmax,neighborlist,markpoints);
}

void PointsWithin_iter(TreeHndl tree,double *ray,float rmax,ListHndl neighborlist,short markpoints){
	bool descend;
	short desition;
	unsigned long i,j;
	double radius,length;

	assert(neighborlist);
	assert(tree);

	EmptyList(neighborlist);

	moveTop(tree);

	do{
		descend = true;
		desition = cutbox(ray,tree->current->boundary_p1,tree->current->boundary_p2,rmax);
		length = FurthestBorder(ray,tree->current->boundary_p1,tree->current->boundary_p2);
		if( FurthestBorder(ray,tree->current->boundary_p1,tree->current->boundary_p2) < rmax ){

			//ClosestBorder(ray,tree->current->boundary_p1,tree->current->boundary_p2) < rmax) ){
			// whole box is outside circle
			descend = false;

			// put all the points into neighborlist
		   	  tree->pointlist->current=tree->current->points;
    		  for(i=0;i<tree->current->npoints;++i){
     			  if(markpoints == 1){
     				  tree->pointlist->current->in_image=TRUE;
     				  tree->pointlist->current->image->in_image=TRUE;
     			  }else if(markpoints == -1){
     				  tree->pointlist->current->in_image=FALSE;
     				  tree->pointlist->current->image->in_image=FALSE;
					  tree->pointlist->current->surface_brightness = tree->pointlist->current->image->surface_brightness = 0.0;
    			  }else if(markpoints == 0){
     				  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
							  ,tree->pointlist->current->id,tree->pointlist->current->image);

					  MoveDownList(neighborlist);
					  PointCopyData(neighborlist->current,tree->pointlist->current);
					  MoveUpList(neighborlist);
				  }

     			  MoveDownList(tree->pointlist);
     		  }

		}else if(cutbox(ray,tree->current->boundary_p1,tree->current->boundary_p2,rmax) == 0 ){  // whole box is outside circle
			descend = false;
		}else if(tree->current->child1 == NULL){  // leaf case

		   	  tree->pointlist->current=tree->current->points;
	   		  for(i=0;i<tree->current->npoints;++i){
	    			  for(j=0,radius=0.0;j<2;++j) radius+=pow(tree->pointlist->current->x[j]-ray[j],2);
	    			  if( radius < rmax*rmax ){
	       				  if(markpoints == 1){
	       					  tree->pointlist->current->in_image=TRUE;
	      					  tree->pointlist->current->image->in_image=TRUE;
	      				  }else if(markpoints == -1){
	      					  tree->pointlist->current->in_image=FALSE;
	     					  tree->pointlist->current->image->in_image=FALSE;
	     					  tree->pointlist->current->surface_brightness = tree->pointlist->current->image->surface_brightness = 0.0;
	       				  }else if(markpoints == 0){
	         				  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
	    						  ,tree->pointlist->current->id,tree->pointlist->current->image);

	         				  MoveDownList(neighborlist);
	         				  PointCopyData(neighborlist->current,tree->pointlist->current);
	         				  MoveUpList(neighborlist);
	      				  }
	    			  }
	    			  MoveDownList(tree->pointlist);
	    		  }
		}
	}while(TreeWalkStep(tree,descend));
}

double ClosestBorder(double *ray,double *p1,double *p2){
	/*  returns the distance from ray[] to the closest
	 *    border of the box,
	 *    < 0 if ray is outside of the box
	 *    > 0 if ray is inside the box
	 */
	double length;

	length = MIN(ray[0]-p1[0],p2[0]-ray[0]);
	length = MIN(ray[1]-p1[1],length);

	return MIN(p2[1]-ray[1],length);
}

void _PointsWithin(TreeHndl tree,double *ray,float *rmax,ListHndl neighborlist,short markpoints){

  int i,j,incell2=1;
  double radius;
  short pass;


  //std::printf("**************************************\nlevel %i\n",tree->current->level);
  //   std::printf("   %i incell=%i\n",tree->current->points->id,incell);

  if(incell){  // not found cell yet

    if( inbox(ray,tree->current->boundary_p1,tree->current->boundary_p2) ){

      // found the box small enough
    	if( cutbox(ray,tree->current->boundary_p1,tree->current->boundary_p2,*rmax)==1
    			|| atLeaf(tree) ){
    		// whole box in circle or a leaf with ray in it

    	  incell=0;
    	  //std::printf("found box with %i points\n",tree->current->npoints);

    	  // this sets ray back to real value once closest leaf bax is found
    	  if( (ray[0]!=realray[0])*(ray[1]!=realray[1]) ){ std::printf("ray != realray _PointsWithin\n"); exit(0);}

    	  ray[0]=realray[0];
    	  ray[1]=realray[1];

    	  if(tree->current->points != NULL) tree->pointlist->current = tree->current->points;

    	  if( atLeaf(tree) ){
    	   	  // if leaf calculate the distance to all the points in cell
    		  for(i=0;i<tree->current->npoints;++i){
    			  for(j=0,radius=0.0;j<2;++j) radius+=pow(tree->pointlist->current->x[j]-ray[j],2);
    			  if( radius < *rmax**rmax ){
       				  if(markpoints == 1){
       					  tree->pointlist->current->in_image=TRUE;
      					  tree->pointlist->current->image->in_image=TRUE;
      				  }else if(markpoints == -1){
      					  tree->pointlist->current->in_image=FALSE;
     					  tree->pointlist->current->image->in_image=FALSE;
     					  tree->pointlist->current->surface_brightness = tree->pointlist->current->image->surface_brightness = 0.0;
       				  }else if(markpoints == 0){
         				  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
    						  ,tree->pointlist->current->id,tree->pointlist->current->image);

         				  MoveDownList(neighborlist);
         				  PointCopyData(neighborlist->current,tree->pointlist->current);
         				  MoveUpList(neighborlist);
      				  }
    			  }
    			  MoveDownList(tree->pointlist);
    		  }
    	  }else{ // put all of points in box into neighborlist
       		  for(i=0;i<tree->current->npoints;++i){
       			  if(markpoints == 1){
       				  tree->pointlist->current->in_image=TRUE;
       				  tree->pointlist->current->image->in_image=TRUE;
       			  }else if(markpoints == -1){
       				  tree->pointlist->current->in_image=FALSE;
       				  tree->pointlist->current->image->in_image=FALSE;
 					  tree->pointlist->current->surface_brightness = tree->pointlist->current->image->surface_brightness = 0.0;
      			  }else if(markpoints == 0){
       				  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
  							  ,tree->pointlist->current->id,tree->pointlist->current->image);

  					  MoveDownList(neighborlist);
  					  PointCopyData(neighborlist->current,tree->pointlist->current);
  					  MoveUpList(neighborlist);
  				  }

       			  MoveDownList(tree->pointlist);
       		  }
    	  }

    	}else{ // keep going down the tree

    	  //std::printf("moving to child1 from level %i\n",tree->current->level);
    	  if(tree->current->child1 !=NULL){
    		  moveToChild(tree,1);
    		  _PointsWithin(tree,ray,rmax,neighborlist,markpoints);
    		  //std::printf("moving up from level %i\n",tree->current->level);
    		  moveUp(tree);

    		  incell2=incell;
    	  }

    	  if(tree->current->child2 !=NULL){
    		  //std::printf("moving to child2 from level %i\n",tree->current->level);
    		  moveToChild(tree,2);
    		  _PointsWithin(tree,ray,rmax,neighborlist,markpoints);
    		  //std::printf("moving up from level %i\n",tree->current->level);
    		  moveUp(tree);
    	  }

    	  // if ray found in second child go back to first to search for neighbors
    	  if( (incell2==1) && (incell==0) ){
    		  if(tree->current->child1 !=NULL){
    			  //std::printf("moving to child1 again from level %i\n",tree->current->level);
    			  moveToChild(tree,1);
    			  _PointsWithin(tree,ray,rmax,neighborlist,markpoints);
    			  //std::printf("moving up from level %i\n",tree->current->level);
    			  moveUp(tree);
    		  }
    	  }
      }
    }  // not in the box

  }else{    // found cell

	  //std::printf("finding neighboring boxes at level = %i\n",tree->current->level);

	  pass=cutbox(ray,tree->current->boundary_p1,tree->current->boundary_p2,*rmax);
	  // does radius cut into the box
	  if( pass ){

		  if(tree->current->points != NULL) tree->pointlist->current = tree->current->points;

		  if( atLeaf(tree) ){  /* leaf case */

			  for(i=0;i<tree->current->npoints;++i){

				  for(j=0,radius=0.0;j<2;++j) radius+=pow(tree->pointlist->current->x[j]-ray[j],2);
				  if( radius < *rmax**rmax ){
					  if(markpoints==1){
						  tree->pointlist->current->in_image=TRUE;
						  tree->pointlist->current->image->in_image=TRUE;
					  }else if(markpoints==-1){
						  tree->pointlist->current->in_image=FALSE;
						  tree->pointlist->current->image->in_image=FALSE;
     					  tree->pointlist->current->surface_brightness = tree->pointlist->current->image->surface_brightness = 0.0;
					  }else if(markpoints==0){
						  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
								  ,tree->pointlist->current->id,tree->pointlist->current->image);

           				  MoveDownList(neighborlist);
           				  PointCopyData(neighborlist->current,tree->pointlist->current);
           				  MoveUpList(neighborlist);
      				  }

				  }
				  MoveDownList(tree->pointlist);
			  }
		  }else if(pass==1){ // whole box is inside radius

			  for(i=0;i<tree->current->npoints;++i){
  				  if(markpoints==1){
   					  tree->pointlist->current->in_image=TRUE;
  					  tree->pointlist->current->image->in_image=TRUE;
  				  }else if(markpoints==-1){
  					  tree->pointlist->current->in_image=FALSE;
 					  tree->pointlist->current->image->in_image=FALSE;
 					  tree->pointlist->current->surface_brightness = tree->pointlist->current->image->surface_brightness = 0.0;
   				  }else if(markpoints==0){
   					  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
   							  ,tree->pointlist->current->id,tree->pointlist->current->image);

   					  MoveDownList(neighborlist);
   					  PointCopyData(neighborlist->current,tree->pointlist->current);
   					  MoveUpList(neighborlist);
   				  }

				  MoveDownList(tree->pointlist);
			  }
		  }else{
			  //std::printf("moving to child1 from level %i\n",tree->current->level);
			  if(tree->current->child1 !=NULL){
				  moveToChild(tree,1);
				  _PointsWithin(tree,ray,rmax,neighborlist,markpoints);
				  //std::printf("moving up from level %i\n",tree->current->level);
				  moveUp(tree);
			  }

			  if(tree->current->child2 !=NULL){
				  //std::printf("moving to child2 from level %i\n",tree->current->level);
				  moveToChild(tree,2);
				  _PointsWithin(tree,ray,rmax,neighborlist,markpoints);
				  //std::printf("moving up from level %i\n",tree->current->level);
				  moveUp(tree);
			  }
		  }

	  }
  }

	  //  std::printf("end of _PointsWithin incell=%i level=%i p1= %e %e %e\n",incell,tree->current->level
	//	,tree->current->boundary_p1[0],tree->current->boundary_p1[1],tree->current->boundary_p1[2]);/**/
  return;
}

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

void FriendsOfFriends(TreeHndl tree,double *start_point,float linkinglength,ListHndl neighborlist
		      ,Point *filter,unsigned long Nfilter,unsigned long *filter_place){
  unsigned long Nlocal_filter,i;
  Point *placemark,*local_filter;
  short malloced=0;

/*   std::printf("entering FOF\n"); */

  if(filter == NULL){ std::printf("FriendsOfFriends cannot handle no unfiltered points\n"); return;}

  EmptyList(neighborlist);

/*   if(filter != NULL){   */
/*     std::printf("filter first id %i\n   x start = %e %e\n",filter[*filter_place].id,start_point[0],start_point[1]); */
/*     MoveToTopList(tree->pointlist); */
/*     for(i=0;i<tree->pointlist->Npoints;++i){ */
/*       if(filter[*filter_place].id == tree->pointlist->current->id) std::printf("x in tree %e %e\n" */
/* 					 ,tree->pointlist->current->x[0],tree->pointlist->current->x[1]); */
/*       MoveDownList(tree->pointlist); */
/*     } */
/*   } */

  incell=1;
  realray[0]=start_point[0];
  realray[1]=start_point[1];

  moveTop(tree);
  _PointsWithin2(tree,start_point,&linkinglength,neighborlist,filter,Nfilter,filter_place,1);

  if(filter == NULL &&  neighborlist->Npoints > 0){
    Nlocal_filter=neighborlist->Npoints;
    //local_filter=(Point *)malloc(Nlocal_filter*sizeof(Point));
	local_filter=NewPointArray(Nlocal_filter,false);

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
	/* find points that are not in filter */
	realray[0]=neighborlist->current->x[0];
	realray[1]=neighborlist->current->x[1];

	_PointsWithin2(tree,neighborlist->current->x,&linkinglength,neighborlist
		      ,local_filter,Nlocal_filter,filter_place,0);

	/* add found points to filter */
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

	_PointsWithin2(tree,neighborlist->current->x,&linkinglength,neighborlist,filter,Nfilter,filter_place,1);
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

/*   std::printf("returning from FOF\n"); */
  return;
}

void _PointsWithin2(TreeHndl tree,double *ray,float *rmax,ListHndl neighborlist
		   ,Point *filter,unsigned long Nfilter,unsigned long *filter_place,short compliment){

/*   id numbers in filter[*filterplace ... Nfilter-1] */
/*   reorders filter[] and moves filterplace as points are found */
/*   if filter=NULL no filter is used */
/*   compliment = 0 search only points not in the filter */
/*              = 1 search only points in filter after filter_place */

  int i,j,incell2=1;
  double radius;
  unsigned long k;

  //if(filter != NULL) std::printf("filter_place=%i Nneighborlist=%i\n",*filter_place,neighborlist->Npoints);

  if((filter != NULL) && (*filter_place >= Nfilter) ) return; // finished filter

  //std::printf("**************************************\nlevel %i\n",tree->current->level);
  //   std::printf("   %i incell=%i\n",tree->current->points->id,incell);

  if(incell){  // not found cell yet

    if( inbox(ray,tree->current->boundary_p1,tree->current->boundary_p2) ){

      // found the box small enough
      if( (tree->current->child1 == NULL)*(tree->current->child2 == NULL)){  // leaf case
    	  incell=0;
    	  //std::printf("found box with %i points\n",tree->current->npoints);

    	  // this sets ray back to real value once closest leaf bax is found
    	  if( (ray[0]!=realray[0])*(ray[1]!=realray[1]) ){ std::printf("ray != realray _PointsWithin\n"); exit(0);}

    	  ray[0]=realray[0];
    	  ray[1]=realray[1];

    	  // calculate the distance to all the points in cell

    	  tree->pointlist->current=tree->current->points;
    	  for(i=0;i<tree->current->npoints;++i){

    		  for(j=0,radius=0.0;j<2;++j) radius+=pow(tree->pointlist->current->x[j]-ray[j],2);
    		  if( radius < *rmax**rmax ){
    			  if(filter == NULL){

    				  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
    						  ,tree->pointlist->current->id,tree->pointlist->current->image);

    				  MoveDownList(neighborlist);
    				  PointCopyData(neighborlist->current,tree->pointlist->current);
    				  MoveUpList(neighborlist);

    			  }else if(compliment){

    				  for(k=*filter_place;k<Nfilter;++k){
    					  if(tree->pointlist->current->id == filter[k].id){
    						  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
    								  ,tree->pointlist->current->id,tree->pointlist->current->image);

    						  MoveDownList(neighborlist);
    						  PointCopyData(neighborlist->current,tree->pointlist->current);
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

    				  for(k=0;k<Nfilter;++k) if(tree->pointlist->current->id == filter[k].id) break;
    				  if(k==Nfilter){
    					  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
    							  ,tree->pointlist->current->id,tree->pointlist->current->image);

    					  MoveDownList(neighborlist);
    					  PointCopyData(neighborlist->current,tree->pointlist->current);
    					  MoveUpList(neighborlist);
    				  }
    			  }

    		  }
    		  MoveDownList(tree->pointlist);
    	  }
      }else{ // keep going down the tree

    	  //std::printf("moving to child1 from level %i\n",tree->current->level);
    	  if(tree->current->child1 !=NULL){
    		  moveToChild(tree,1);
    		  _PointsWithin2(tree,ray,rmax,neighborlist,filter,Nfilter,filter_place,1);
    		  //std::printf("moving up from level %i\n",tree->current->level);
    		  moveUp(tree);

    		  incell2=incell;
    	  }

    	  if(tree->current->child2 !=NULL){
    		  //std::printf("moving to child2 from level %i\n",tree->current->level);
    		  moveToChild(tree,2);
    		  _PointsWithin2(tree,ray,rmax,neighborlist,filter,Nfilter,filter_place,1);
    		  //std::printf("moving up from level %i\n",tree->current->level);
    		  moveUp(tree);
    	  }

    	  // if ray found in second child go back to first to search for neighbors
    	  if( (incell2==1) && (incell==0) ){
    		  if(tree->current->child1 !=NULL){
    			  //std::printf("moving to child1 again from level %i\n",tree->current->level);
    			  moveToChild(tree,1);
    			  _PointsWithin2(tree,ray,rmax,neighborlist,filter,Nfilter,filter_place,1);
    			  //std::printf("moving up from level %i\n",tree->current->level);
    			  moveUp(tree);
    		  }
    	  }
      }
    }  // not in the box
    //std::printf("not in box \n");}
  }else{    // found cell

	  //std::printf("finding neighboring boxes at level = %i\n",tree->current->level);

	  // does radius cut into the box
	  if( cutbox(ray,tree->current->boundary_p1,tree->current->boundary_p2,*rmax) ){

		  if( (tree->current->child1 == NULL)*(tree->current->child2 == NULL)){  /* leaf case */

			  tree->pointlist->current=tree->current->points;
			  for(i=0;i<tree->current->npoints;++i){

				  for(j=0,radius=0.0;j<2;++j) radius+=pow(tree->pointlist->current->x[j]-ray[j],2);
				  if( radius < *rmax**rmax ){
					  if(filter == NULL){

						  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
								  ,tree->pointlist->current->id,tree->pointlist->current->image);

						  MoveDownList(neighborlist);
						  PointCopyData(neighborlist->current,tree->pointlist->current);
						  MoveUpList(neighborlist);

					  }else if(compliment){

						  for(k=*filter_place;k<Nfilter;++k){
							  if(tree->pointlist->current->id == filter[k].id){
								  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
										  ,tree->pointlist->current->id,tree->pointlist->current->image);

								  MoveDownList(neighborlist);
								  PointCopyData(neighborlist->current,tree->pointlist->current);
								  MoveUpList(neighborlist);

								  PointCopyData(&tmp,&filter[k]);
								  for(;k>*filter_place;k--) PointCopyData(&filter[k],&filter[k-1]);
								  PointCopyData(&filter[*filter_place],&tmp);
								  ++*filter_place;
								  break;
							  }
						  }

					  }else{ // do compliment of filter

						  for(k=0;k<Nfilter;++k) if(tree->pointlist->current->id == filter[k].id) break;
						  if(k==Nfilter){
							  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
									  ,tree->pointlist->current->id,tree->pointlist->current->image);

							  MoveDownList(neighborlist);
							  PointCopyData(neighborlist->current,tree->pointlist->current);
							  MoveUpList(neighborlist);

						  }
					  }
				  }
				  MoveDownList(tree->pointlist);
			  }

		  }else{
			  //std::printf("moving to child1 from level %i\n",tree->current->level);
			  if(tree->current->child1 !=NULL){
				  moveToChild(tree,1);
				  _PointsWithin2(tree,ray,rmax,neighborlist,filter,Nfilter,filter_place,1);
				  //std::printf("moving up from level %i\n",tree->current->level);
				  moveUp(tree);
			  }

			  if(tree->current->child2 !=NULL){
				  //std::printf("moving to child2 from level %i\n",tree->current->level);
				  moveToChild(tree,2);
				  _PointsWithin2(tree,ray,rmax,neighborlist,filter,Nfilter,filter_place,1);
				  //std::printf("moving up from level %i\n",tree->current->level);
				  moveUp(tree);
			  }
		  }

	  }//else{std::printf("box too distant at level %i\n",tree->current->level);}
  }

	  //  std::printf("end of _PointsWithin incell=%i level=%i p1= %e %e %e\n",incell,tree->current->level
	//	,tree->current->boundary_p1[0],tree->current->boundary_p1[1],tree->current->boundary_p1[2]);/**/
  return;
}

void TransformPoints(ListHndl listout,ListHndl listin){
  /* tranforms between image and source planes */
  /* the image pointers in listin must be set properly */
  /* the image pointers of listout will be set to listin members - it would be better to fix this */

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

void ClearAllMarks(TreeHndl tree){
	unsigned long i;

	MoveToTopList(tree->pointlist);
	for(i=0;i<tree->pointlist->Npoints;++i){
		tree->pointlist->current->in_image=FALSE;
		tree->pointlist->current->image->in_image=FALSE;
		MoveDownList(tree->pointlist);
	}
}


void FindAllBoxNeighbors(TreeHndl tree,Point *point,ListHndl neighbors){
	// finds all the leaves that are neighboring branch
	// points outside of grid have no box neighbors
	void _FindAllBoxNeighbors(TreeHndl tree,Branch *leaf,ListHndl neighbors);
	static int count=0;

	++count;
	EmptyList(neighbors);

	// point is outside of initial region
	if(!inbox(point->x,tree->top->boundary_p1,tree->top->boundary_p2)) return;

	tree->current=point->leaf;

	// find smallest box that surrounds box and its neighbors
	moveUp(tree);
	while( (tree->current->boundary_p1[0]==point->leaf->boundary_p1[0] && point->leaf->boundary_p1[0] != tree->top->boundary_p1[0] )
			|| (tree->current->boundary_p1[1]==point->leaf->boundary_p1[1] && point->leaf->boundary_p1[1] != tree->top->boundary_p1[1])
			|| (tree->current->boundary_p2[0]==point->leaf->boundary_p2[0] && point->leaf->boundary_p2[0] != tree->top->boundary_p2[0])
			|| (tree->current->boundary_p2[1]==point->leaf->boundary_p2[1] && point->leaf->boundary_p2[1] != tree->top->boundary_p2[1]) ){
		moveUp(tree);
	}

	_FindAllBoxNeighbors(tree,point->leaf,neighbors);

	return;
}

// this should be made into a loop instead of a recursion
void _FindAllBoxNeighbors(TreeHndl tree,Branch *leaf,ListHndl neighbors){

	if(  leaf->boundary_p1[0]<=tree->current->boundary_p2[0]
	  && leaf->boundary_p2[0]>=tree->current->boundary_p1[0]
	  && leaf->boundary_p1[1]<=tree->current->boundary_p2[1]
	  && leaf->boundary_p2[1]>=tree->current->boundary_p1[1]){

		if(tree->current->npoints == tree->Nbucket){
			if(tree->current->number != leaf->number){

				//InsertAfterCurrentKist(neighbors,tree->current->points);

				InsertAfterCurrent(neighbors,tree->current->points->x
					,tree->current->points->id
					,tree->current->points->image);
				MoveDownList(neighbors);
				PointCopyData(neighbors->current,tree->current->points);

			}
			return;
		}

		if(tree->current->child1 !=NULL){
			moveToChild(tree,1);
			_FindAllBoxNeighbors(tree,leaf,neighbors);
			moveUp(tree);
		}

		if(tree->current->child2 !=NULL){
			moveToChild(tree,2);
			_FindAllBoxNeighbors(tree,leaf,neighbors);
			moveUp(tree);
		}
	}

	return;
}

void PrintImages(ImageInfo *images,long Nimages){
	unsigned long i,j;

	std::printf("%li",Nimages);
	for(i=0;i<Nimages;++i){
		std::printf("%li\n",images[i].imagekist->Nunits());
		for(j=0,MoveToTopKist(images[i].imagekist);j<images[i].imagekist->Nunits();++j,MoveDownKist(images[i].imagekist))
			std::printf("%e %e  %e\n",getCurrentKist(images[i].imagekist)->x[0],getCurrentKist(images[i].imagekist)->x[1]
			                         ,getCurrentKist(images[i].imagekist)->gridsize);
	}
}

void PrintImageInfo(ImageInfo *image){

	std::printf(" PrintImageInfo\n");
	std::printf("  Npoints = %li  area = %e +/- %e\n",image->imagekist->Nunits(),image->area,image->area_error);
	std::printf("  gridrange = %e %e %e\n",image->gridrange[0],image->gridrange[1],image->gridrange[2]);
	std::printf("  borders inner N = %li  outer N = %li\n",image->innerborder->Nunits(),image->outerborder->Nunits());
}
