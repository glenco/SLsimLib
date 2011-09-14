/*
 * KistDriver.c
 *
 *  Created on: Nov 16, 2010
 *      Author: bmetcalf
 *
 *      routines imported from TreeDriver.c made to use Kist instead of List
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <Tree.h>
#include <Kist.h>
#include <KistDriver.h>

static int Nbucket = 1;   // must be =1 if each leaf is to coincide with each cell
static int incell;
static double realray[2];

void FindAllBoxNeighborsKist(TreeHndl tree,Point *point,KistHndl neighbors){
	// finds all the leaves that are neighboring branch
	// points outside of grid have no box neighbors
	void _FindAllBoxNeighborsKist(TreeHndl tree,Branch *leaf,KistHndl neighbors);
	static int count=0;

	++count;
	EmptyKist(neighbors);

	// point is outside of initial region
	if(!inbox(point->x,tree->top->boundery_p1,tree->top->boundery_p2)) return;

	tree->current = point->leaf;

	// find smallest box that surrounds box and its neighbors
	moveUp(tree);
	while( (tree->current->boundery_p1[0]==point->leaf->boundery_p1[0] && point->leaf->boundery_p1[0] != tree->top->boundery_p1[0] )
			|| (tree->current->boundery_p1[1]==point->leaf->boundery_p1[1] && point->leaf->boundery_p1[1] != tree->top->boundery_p1[1])
			|| (tree->current->boundery_p2[0]==point->leaf->boundery_p2[0] && point->leaf->boundery_p2[0] != tree->top->boundery_p2[0])
			|| (tree->current->boundery_p2[1]==point->leaf->boundery_p2[1] && point->leaf->boundery_p2[1] != tree->top->boundery_p2[1]) ){
		moveUp(tree);
	}

	assert(inbox(point->x,tree->current->boundery_p1,tree->current->boundery_p2));
	_FindAllBoxNeighborsKist_iter(tree,point->leaf,neighbors);

	return;
}

// There is an iterative version of this recursive function below
//    This one has been known to cause stack overflow.
void _FindAllBoxNeighborsKist(TreeHndl tree,Branch *leaf,KistHndl neighbors){

	if(  leaf->boundery_p1[0] <= tree->current->boundery_p2[0]
	  && leaf->boundery_p2[0] >= tree->current->boundery_p1[0]
	  && leaf->boundery_p1[1] <= tree->current->boundery_p2[1]
	  && leaf->boundery_p2[1] >= tree->current->boundery_p1[1]){

		if(tree->current->npoints == Nbucket){
			if(tree->current->number != leaf->number){

				InsertAfterCurrentKist(neighbors,tree->current->points);
				MoveDownKist(neighbors);
			}
			return;
		}

		if(tree->current->child1 !=NULL){
			moveToChild(tree,1);
			_FindAllBoxNeighborsKist(tree,leaf,neighbors);
			moveUp(tree);
		}

		if(tree->current->child2 !=NULL){
			moveToChild(tree,2);
			_FindAllBoxNeighborsKist(tree,leaf,neighbors);
			moveUp(tree);
		}
	}

	return;
}

/* Iterative instead of recursive method for finding neighbors
 *   the tree->current must be preset so that leaf is within it.
 *
 *   The recursive routine _FindAllBoxNeighborsKist has caused stack overflow.
 */

void _FindAllBoxNeighborsKist_iter(TreeHndl tree,Branch *leaf,KistHndl neighbors){

	Boolean allowDescent = True;
	long level = tree->current->level;
	//Branch *brother = tree->current->brother;
	//unsigned long count = 0;

	EmptyKist(neighbors);

	while(TreeWalkStep(tree,allowDescent) && tree->current->level > level ){

		if(  leaf->boundery_p1[0] <= tree->current->boundery_p2[0]
		  && leaf->boundery_p2[0] >= tree->current->boundery_p1[0]
	      && leaf->boundery_p1[1] <= tree->current->boundery_p2[1]
	      && leaf->boundery_p2[1] >= tree->current->boundery_p1[1]){

			if(tree->current->npoints == Nbucket && tree->current != leaf){
				InsertAfterCurrentKist(neighbors,tree->current->points);
				MoveDownKist(neighbors);
			}

			allowDescent = True;
		}else{
			allowDescent = False;
		}

		//++count;
		//if(tree->Nbranches > 0) assert(count <= tree->Nbranches );

	}

	return;
}

void PointsWithinKist(TreeHndl tree,double *ray,float rmax,KistHndl neighborkist,short markpoints){
/*
 *  finds all points in tree that lie within rmax of the point ray[]
 *   markpoints = 0  does not change in_image variable in any point, gives a list of neighbors
 *              = 1  makes in_image=True for all points in image, gives no list of neighbors
 *              = -1 makes in_image=False for all points in image to reset, gives no list of neighbors
 */

  if(markpoints==0) EmptyKist(neighborkist);

  realray[0]=ray[0];
  realray[1]=ray[1];

  moveTop(tree);
  if( inbox(ray,tree->current->boundery_p1,tree->current->boundery_p2) == 0 ){
    printf("Warning: in PointsWithinKist, ray is not inside the simulation box\n    should work in any case\n      ray= %e %e\n     boundery p1 = %e %e p2 = %e %e\n",ray[0],ray[1]
	   ,tree->current->boundery_p1[0],tree->current->boundery_p1[1]
	   ,tree->current->boundery_p2[0],tree->current->boundery_p2[1]);

    ray[0]=MAX(ray[0],tree->current->boundery_p1[0]);
    ray[0]=MIN(ray[0],tree->current->boundery_p2[0]);

    ray[1]=MAX(ray[1],tree->current->boundery_p1[1]);
    ray[1]=MIN(ray[1],tree->current->boundery_p2[1]);
  }
  incell=1;

  _PointsWithinKist(tree,ray,&rmax,neighborkist,markpoints);
}

void _PointsWithinKist(TreeHndl tree,double *ray,float *rmax,KistHndl neighborkist
		,short markpoints){

  int i,j,incell2=1;
  double radius;
  short pass;


  //printf("**************************************\nlevel %i\n",tree->current->level);
  //   printf("   %i incell=%i\n",tree->current->points->id,incell);

  if(incell){  // not found cell yet

    if( inbox(ray,tree->current->boundery_p1,tree->current->boundery_p2) ){

      // found the box small enough
    	if( cutbox(ray,tree->current->boundery_p1,tree->current->boundery_p2,*rmax)==1
    			|| (tree->current->child1 == NULL)*(tree->current->child2 == NULL) ){
    		// whole box in circle or a leaf with ray in it

    	  incell=0;
    	  //printf("found box with %i points\n",tree->current->npoints);

    	  // this sets ray back to real value once closest leaf bax is found
    	  if( (ray[0]!=realray[0])*(ray[1]!=realray[1]) ){ printf("ray != realray _PointsWithinKist\n"); exit(0);}

    	  ray[0]=realray[0];
    	  ray[1]=realray[1];

    	  tree->pointlist->current=tree->current->points;

    	  if((tree->current->child1 == NULL)*(tree->current->child2 == NULL)){
    	   	  // if leaf calculate the distance to all the points in cell
    		  for(i=0;i<tree->current->npoints;++i){
    			  for(j=0,radius=0.0;j<2;++j) radius+=pow(tree->pointlist->current->x[j]-ray[j],2);
    			  if( radius < *rmax**rmax ){
       				  if(markpoints == 1){
       					  tree->pointlist->current->in_image=True;
      					  tree->pointlist->current->image->in_image=True;
      				  }else if(markpoints == -1){
      					  tree->pointlist->current->in_image=False;
     					  tree->pointlist->current->image->in_image=False;
     					  tree->pointlist->current->surface_brightness = tree->pointlist->current->image->surface_brightness = 0.0;
       				  }else if(markpoints == 0){
         				  InsertAfterCurrentKist(neighborkist,tree->pointlist->current);
      				  }
    			  }
    			  MoveDownList(tree->pointlist);
    		  }
    	  }else{ // put all of points in box into getCurrentKist(imagekist)
       		  for(i=0;i<tree->current->npoints;++i){
       			  if(markpoints == 1){
       				  tree->pointlist->current->in_image=True;
       				  tree->pointlist->current->image->in_image=True;
       			  }else if(markpoints == -1){
       				  tree->pointlist->current->in_image=False;
       				  tree->pointlist->current->image->in_image=False;
 					  tree->pointlist->current->surface_brightness = tree->pointlist->current->image->surface_brightness = 0.0;
      			  }else if(markpoints == 0){
       				  InsertAfterCurrentKist(neighborkist,tree->pointlist->current);
 				  }

       			  MoveDownList(tree->pointlist);
       		  }
    	  }

    	}else{ // keep going down the tree

    	  //printf("moving to child1 from level %i\n",tree->current->level);
    	  if(tree->current->child1 !=NULL){
    		  moveToChild(tree,1);
    		  _PointsWithinKist(tree,ray,rmax,neighborkist,markpoints);
    		  //printf("moving up from level %i\n",tree->current->level);
    		  moveUp(tree);

    		  incell2=incell;
    	  }

    	  if(tree->current->child2 !=NULL){
    		  //printf("moving to child2 from level %i\n",tree->current->level);
    		  moveToChild(tree,2);
    		  _PointsWithinKist(tree,ray,rmax,neighborkist,markpoints);
    		  //printf("moving up from level %i\n",tree->current->level);
    		  moveUp(tree);
    	  }

    	  // if ray found in second child go back to first to search for neighbors
    	  if( (incell2==1) && (incell==0) ){
    		  if(tree->current->child1 !=NULL){
    			  //printf("moving to child1 again from level %i\n",tree->current->level);
    			  moveToChild(tree,1);
    			  _PointsWithinKist(tree,ray,rmax,neighborkist,markpoints);
    			  //printf("moving up from level %i\n",tree->current->level);
    			  moveUp(tree);
    		  }
    	  }
      }
    }  // not in the box

  }else{    // found cell

	  //printf("finding neighboring boxes at level = %i\n",tree->current->level);

	  pass=cutbox(ray,tree->current->boundery_p1,tree->current->boundery_p2,*rmax);
	  // does radius cut into the box
	  if( pass ){

		  if( (tree->current->child1 == NULL)*(tree->current->child2 == NULL)  ){  /* leaf case */

			  tree->pointlist->current=tree->current->points;
			  for(i=0;i<tree->current->npoints;++i){

				  for(j=0,radius=0.0;j<2;++j) radius+=pow(tree->pointlist->current->x[j]-ray[j],2);
				  if( radius < *rmax**rmax ){
					  if(markpoints==1){
						  tree->pointlist->current->in_image=True;
						  tree->pointlist->current->image->in_image=True;
					  }else if(markpoints==-1){
						  tree->pointlist->current->in_image=False;
						  tree->pointlist->current->image->in_image=False;
     					  tree->pointlist->current->surface_brightness = tree->pointlist->current->image->surface_brightness = 0.0;
					  }else if(markpoints==0){
						  InsertAfterCurrentKist(neighborkist,tree->pointlist->current);
     				  }

				  }
				  MoveDownList(tree->pointlist);
			  }
		  }else if(pass==1){ // whole box is inside radius
			  tree->pointlist->current=tree->current->points;
			  for(i=0;i<tree->current->npoints;++i){
  				  if(markpoints==1){
   					  tree->pointlist->current->in_image=True;
  					  tree->pointlist->current->image->in_image=True;
  				  }else if(markpoints==-1){
  					  tree->pointlist->current->in_image=False;
 					  tree->pointlist->current->image->in_image=False;
 					  tree->pointlist->current->surface_brightness = tree->pointlist->current->image->surface_brightness = 0.0;
   				  }else if(markpoints==0){
   					  InsertAfterCurrentKist(neighborkist,tree->pointlist->current);
  				  }

				  MoveDownList(tree->pointlist);
			  }
		  }else{
			  //printf("moving to child1 from level %i\n",tree->current->level);
			  if(tree->current->child1 !=NULL){
				  moveToChild(tree,1);
				  _PointsWithinKist(tree,ray,rmax,neighborkist,markpoints);
				  //printf("moving up from level %i\n",tree->current->level);
				  moveUp(tree);
			  }

			  if(tree->current->child2 !=NULL){
				  //printf("moving to child2 from level %i\n",tree->current->level);
				  moveToChild(tree,2);
				  _PointsWithinKist(tree,ray,rmax,neighborkist,markpoints);
				  //printf("moving up from level %i\n",tree->current->level);
				  moveUp(tree);
			  }
		  }

	  }
  }

	  //  printf("end of _PointsWithinKist incell=%i level=%i p1= %e %e %e\n",incell,tree->current->level
	//	,tree->current->boundery_p1[0],tree->current->boundery_p1[1],tree->current->boundery_p1[2]);/**/
  return;
}

Point *NearestNeighborKist(TreeHndl tree,double *ray,int Nneighbors,KistHndl neighborkist){
	/* nearest neighbor points
	 *    this is a kludge that relies on NearestNeighbor and translates the list to a kist
	 */
	ListHndl neighborlist = NewList();
	Point *point = 0;
	unsigned long i;

	point = NearestNeighbor(tree,ray,Nneighbors,neighborlist,0);

	// convert from point array to exported point kist
	EmptyKist(neighborkist);
	MoveToTopList(neighborlist);
	for(i = 0; i < neighborlist->Npoints ;++i){
		InsertAfterCurrentKist(neighborkist,neighborlist->current->image->image);
		MoveDownKist(neighborkist);
		MoveDownList(neighborlist);
	}

	assert(neighborlist->Npoints == neighborkist->Nunits);

	EmptyList(neighborlist);
	free(neighborlist);

	return point;
}
