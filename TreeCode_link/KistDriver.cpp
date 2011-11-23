/*! \file
 * KistDriver.c
 *
 *  Created on: Nov 16, 2010
 *      Author: bmetcalf
 *
 *      routines imported from TreeDriver.c made to use Kist instead of List
 */
/*#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <Tree.h>
#include <Kist.h>
#include <KistDriver.h>*/

#include <slsimlib.h>

//static int Nbucket = 1;   // must be =1 if each leaf is to coincide with each cell
static int incell;
static double realray[2];

/** \ingroup ImageFundingL2
 *
 * \brief Finds all the leaves that are neighboring a point.
*
* Points outside of grid have no box neighbors
* Warning: Does not take empty leaves into account.
*/

void FindAllBoxNeighborsKist(TreeHndl tree,Point *point,KistHndl neighbors){
	void _FindAllBoxNeighborsKist(TreeHndl tree,Branch *leaf,KistHndl neighbors);
	static int count=0;

	++count;
	EmptyKist(neighbors);

	// point is outside of initial region
	if(!inbox(point->x,tree->top->boundary_p1,tree->top->boundary_p2)) return;

	tree->current = point->leaf;

	// find smallest box that surrounds box and its neighbors
	moveUp(tree);
	while( (tree->current->boundary_p1[0]==point->leaf->boundary_p1[0] && point->leaf->boundary_p1[0] != tree->top->boundary_p1[0] )
			|| (tree->current->boundary_p1[1]==point->leaf->boundary_p1[1] && point->leaf->boundary_p1[1] != tree->top->boundary_p1[1])
			|| (tree->current->boundary_p2[0]==point->leaf->boundary_p2[0] && point->leaf->boundary_p2[0] != tree->top->boundary_p2[0])
			|| (tree->current->boundary_p2[1]==point->leaf->boundary_p2[1] && point->leaf->boundary_p2[1] != tree->top->boundary_p2[1]) ){
		moveUp(tree);
	}

	assert(inbox(point->x,tree->current->boundary_p1,tree->current->boundary_p2));
	_FindAllBoxNeighborsKist_iter(tree,point->leaf,neighbors);

	return;
}
/**  \ingroup LowLevel
 * A recessive function that was used in FindAllBoxNeighborsKist().
*    It has been known to cause stack overflow. Use _FindAllBoxNeighborsKist_iter instead.
*/
void _FindAllBoxNeighborsKist(TreeHndl tree,Branch *leaf,KistHndl neighbors){

	if(  leaf->boundary_p1[0] <= tree->current->boundary_p2[0]
	  && leaf->boundary_p2[0] >= tree->current->boundary_p1[0]
	  && leaf->boundary_p1[1] <= tree->current->boundary_p2[1]
	  && leaf->boundary_p2[1] >= tree->current->boundary_p1[1]){

		if( atLeaf(tree) ){
			assert(tree->current->npoints <= tree->Nbucket);
			//if(tree->current->npoints == Nbucket){
			// What if number is > than Nbucket and it is not in a leaf
			if(tree->current->number != leaf->number && tree->current->npoints > 0){
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
/**  \ingroup LowLevel
 * Used in FindAllBoxNeighborsKist to walk tree for neighbors.
*/
void _FindAllBoxNeighborsKist_iter(TreeHndl tree,Branch *leaf,KistHndl neighbors){

	/** Iterative instead of recursive method for finding neighbors
	 *   the tree->current must be preset so that leaf is within it.
	 *
	 *   The recursive routine _FindAllBoxNeighborsKist has caused stack overflow.
	 */

	bool allowDescent = true;
	long level = tree->current->level;

	EmptyKist(neighbors);

	while(TreeWalkStep(tree,allowDescent) && tree->current->level > level ){

		if(  leaf->boundary_p1[0] <= tree->current->boundary_p2[0]
		  && leaf->boundary_p2[0] >= tree->current->boundary_p1[0]
	      && leaf->boundary_p1[1] <= tree->current->boundary_p2[1]
	      && leaf->boundary_p2[1] >= tree->current->boundary_p1[1]){

			if( atLeaf(tree) ){
				if(tree->current != leaf && tree->current->npoints > 0){
					InsertAfterCurrentKist(neighbors,tree->current->points);
					MoveDownKist(neighbors);
				}
			}

			allowDescent = true;
		}else{
			allowDescent = false;
		}
	}

	return;
}
/** \ingroup ImageFundingL2
 * \brief Finds points within an ellipse
 *
 * This becomes less efficient when the ellipse is very elongated.  Could
 * be improved by incorporating the test of it being in the ellipse into the
 * tree walk.
 *
 * The
 */
void PointsWithinEllipKist(
	TreeHndl tree    /// tree of points
	,double *ray     /// center of ellipse
	,float rmax      /// major axis
	,float rmin     /// minor axis
	,float posangle  /// position angle of major axis, smallest angle between the x-axis and the long axis
	,KistHndl neighborkist  /// output neighbor kist, will be emptied if it contains anything on entry
	){
	unsigned long i,Ntmp;
	double *xtmp,x,y,cs,sn;

	EmptyKist(neighborkist);

	if(rmax <=0.0 || rmin <= 0.0) return;

	// find point within a circle circumscribes the ellipse
	PointsWithinKist(tree,ray,rmax,neighborkist,false);

	cs = cos(posangle);
	sn = sin(posangle);
	Ntmp = neighborkist->Nunits;
	for(i=0,MoveToTopKist(neighborkist);i<Ntmp;++i,MoveDownKist(neighborkist) ){
		xtmp = getCurrentKist(neighborkist)->x;
		x = xtmp[0]*cs - xtmp[1]*sn;
		y = xtmp[0]*sn + xtmp[1]*cs;
		if( pow(x/rmax,2) + pow(y/rmin,2) > 1)	TakeOutCurrentKist(neighborkist);
	}
	return;
}
/** \ingroup ImageFundingL2
 * \brief Finds all points in tree that lie within rmax of the point ray[]
 *
 *   markpoints = 0  does not change in_image variable in any point, gives a list of neighbors
 *              = 1  makes in_image=true for all points in image, gives no list of neighbors
 *              = -1 makes in_image=false for all points in image to reset, gives no list of neighbors
 */
void PointsWithinKist(
		TreeHndl tree    /// tree of points
		,double *ray     /// center of circle
		,float rmax      /// radius of circle
		,KistHndl neighborkist  /// output neighbor kist, will be emptied if it contains anything on entry
		,short markpoints       /// see comment
		){

  if(markpoints==0) EmptyKist(neighborkist);

  realray[0]=ray[0];
  realray[1]=ray[1];

  moveTop(tree);
  if( inbox(ray,tree->current->boundary_p1,tree->current->boundary_p2) == 0 ){
    printf("Warning: in PointsWithinKist, ray is not inside the simulation box\n    should work in any case\n      ray= %e %e\n     boundary p1 = %e %e p2 = %e %e\n",ray[0],ray[1]
	   ,tree->current->boundary_p1[0],tree->current->boundary_p1[1]
	   ,tree->current->boundary_p2[0],tree->current->boundary_p2[1]);

    ray[0]=MAX(ray[0],tree->current->boundary_p1[0]);
    ray[0]=MIN(ray[0],tree->current->boundary_p2[0]);

    ray[1]=MAX(ray[1],tree->current->boundary_p1[1]);
    ray[1]=MIN(ray[1],tree->current->boundary_p2[1]);
  }
  incell=1;

  _PointsWithinKist(tree,ray,&rmax,neighborkist,markpoints);

  return;
}
/** \ingroup LowLevel
 * Used in PointsWithinKist() to walk tree.*/
void _PointsWithinKist(TreeHndl tree,double *ray,float *rmax,KistHndl neighborkist
		,short markpoints){

  int i,j,incell2=1;
  double radius;
  short pass;


  //printf("**************************************\nlevel %i\n",tree->current->level);
  //   printf("   %i incell=%i\n",tree->current->points->id,incell);

  if(incell){  // not found cell yet

    if( inbox(ray,tree->current->boundary_p1,tree->current->boundary_p2) ){

      // found the box small enough
    	if( cutbox(ray,tree->current->boundary_p1,tree->current->boundary_p2,*rmax)==1
    			|| atLeaf(tree) ){
    		// whole box in circle or a leaf with ray in it

    	  incell=0;

    	  // this sets ray back to real value once closest leaf bax is found
    	  if( (ray[0]!=realray[0])*(ray[1]!=realray[1]) ){ printf("ray != realray _PointsWithinKist\n"); exit(0);}

    	  ray[0]=realray[0];
    	  ray[1]=realray[1];

    	  if(tree->current->points != NULL) tree->pointlist->current=tree->current->points;

    	  if( atLeaf(tree) ){
    	   	  // if leaf calculate the distance to all the points in cell
    		  for(i=0;i<tree->current->npoints;++i){
    			  for(j=0,radius=0.0;j<2;++j) radius+=pow(tree->pointlist->current->x[j]-ray[j],2);
    			  if( radius < *rmax**rmax ){
       				  if(markpoints == 1){
       					  tree->pointlist->current->in_image=true;
      					  tree->pointlist->current->image->in_image=true;
      				  }else if(markpoints == -1){
      					  tree->pointlist->current->in_image=false;
     					  tree->pointlist->current->image->in_image=false;
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
       				  tree->pointlist->current->in_image=true;
       				  tree->pointlist->current->image->in_image=true;
       			  }else if(markpoints == -1){
       				  tree->pointlist->current->in_image=false;
       				  tree->pointlist->current->image->in_image=false;
 					  tree->pointlist->current->surface_brightness = tree->pointlist->current->image->surface_brightness = 0.0;
      			  }else if(markpoints == 0){
       				  InsertAfterCurrentKist(neighborkist,tree->pointlist->current);
 				  }

       			  MoveDownList(tree->pointlist);
       		  }
    	  }

    	}else{ // keep going down the tree

     	  if(tree->current->child1 !=NULL){
    		  moveToChild(tree,1);
    		  _PointsWithinKist(tree,ray,rmax,neighborkist,markpoints);
    		  moveUp(tree);

    		  incell2=incell;
    	  }

    	  if(tree->current->child2 !=NULL){
    		  moveToChild(tree,2);
    		  _PointsWithinKist(tree,ray,rmax,neighborkist,markpoints);
    		  moveUp(tree);
    	  }

    	  // if ray found in second child go back to first to search for neighbors
    	  if( (incell2==1) && (incell==0) ){
    		  if(tree->current->child1 !=NULL){
    			  moveToChild(tree,1);
    			  _PointsWithinKist(tree,ray,rmax,neighborkist,markpoints);
    			  moveUp(tree);
    		  }
    	  }
      }
    }  // not in the box

  }else{    // found cell

	  pass=cutbox(ray,tree->current->boundary_p1,tree->current->boundary_p2,*rmax);
	  // does radius cut into the box
	  if( pass ){

    	  if(tree->current->points != NULL) tree->pointlist->current=tree->current->points;

		  if( atLeaf(tree)  ){  /* leaf case */

			  for(i=0;i<tree->current->npoints;++i){

				  for(j=0,radius=0.0;j<2;++j) radius+=pow(tree->pointlist->current->x[j]-ray[j],2);
				  if( radius < *rmax**rmax ){
					  if(markpoints==1){
						  tree->pointlist->current->in_image=true;
						  tree->pointlist->current->image->in_image=true;
					  }else if(markpoints==-1){
						  tree->pointlist->current->in_image=false;
						  tree->pointlist->current->image->in_image=false;
     					  tree->pointlist->current->surface_brightness = tree->pointlist->current->image->surface_brightness = 0.0;
					  }else if(markpoints==0){
						  InsertAfterCurrentKist(neighborkist,tree->pointlist->current);
     				  }

				  }
				  MoveDownList(tree->pointlist);
			  }
		  }else if(pass==1){ // whole box is inside radius

			  for(i=0;i<tree->current->npoints;++i){
  				  if(markpoints==1){
   					  tree->pointlist->current->in_image=true;
  					  tree->pointlist->current->image->in_image=true;
  				  }else if(markpoints==-1){
  					  tree->pointlist->current->in_image=false;
 					  tree->pointlist->current->image->in_image=false;
 					  tree->pointlist->current->surface_brightness = tree->pointlist->current->image->surface_brightness = 0.0;
   				  }else if(markpoints==0){
   					  InsertAfterCurrentKist(neighborkist,tree->pointlist->current);
  				  }

				  MoveDownList(tree->pointlist);
			  }
		  }else{
			  if(tree->current->child1 !=NULL){
				  moveToChild(tree,1);
				  _PointsWithinKist(tree,ray,rmax,neighborkist,markpoints);
				  moveUp(tree);
			  }

			  if(tree->current->child2 !=NULL){
				  moveToChild(tree,2);
				  _PointsWithinKist(tree,ray,rmax,neighborkist,markpoints);
				  moveUp(tree);
			  }
		  }

	  }
  }

  return;
}

/** \ingroup ImageFundingL2
 *
 * \brief Finds nearest neighbor points to ray.
 *
 *    This is a kludge that relies on NearestNeighbor which uses a List and translates
 *    the list to a kist.  Could be rewritten.
 */
Point *NearestNeighborKist(TreeHndl tree,double *ray,int Nneighbors,KistHndl neighborkist){
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