/*! \file
 * KistDriver.c
 *
 *  Created on: Nov 16, 2010
 *      Author: bmetcalf
 *
 *      routines imported from TreeDriver.c made to use Kist instead of List
 */

#include "slsimlib.h"

bool TreeStruct::Test(){
  
  
  //Branch *treecur = &branch;
  Point *listcur = pointlist->current;
  TreeStruct::iterator treeit(top);
//  moveTop();
  do{
    //pointlist->current = current->points;
    //for(int k = 0; k < current->npoints ; ++k,MoveDownList(pointlist)){
      pointlist->current = (*treeit)->points;
      for(int k = 0; k < (*treeit)->npoints ; ++k,MoveDownList(pointlist)){
      assert( inbox(pointlist->current->x,(*treeit)->boundary_p1,(*treeit)->boundary_p2) );
    }
  }while(treeit.TreeWalkStep(true));

  //current = treecur;
  pointlist->current = listcur;

  return true;
}

/** \ingroup ImageFundingL2
 *
 * \brief Finds all the leaves that are neighboring a point.
*
* Points outside of grid have no box neighbors
* Warning: Does not take empty leaves into account.
*/

void TreeStruct::FindAllBoxNeighborsKist(Point *point,Kist<Point> * neighbors) const{
	neighbors->Empty();

	//Kist<Point>* testkist = new Kist<Point>;

	for(std::list<Branch *>::iterator it = point->leaf->neighbors.begin();
			it != point->leaf->neighbors.end() ; ++it){
		assert((*it)->npoints <= 1);
		//std::cout << point->leaf->neighbors.size() << std::endl;
		if((*it)->npoints == 1) neighbors->InsertAfterCurrent((*it)->points);
	}

/*
	// point is outside of initial region
	if(!inbox(point->x,top->boundary_p1,top->boundary_p2)) return;

	assert(point->leaf);
	current = point->leaf;

	//std::cout << "tree current " << current << std::endl;

	// find smallest box that surrounds box and its neighbors
	//printTree(tree);
	moveUp();
	while( (current->boundary_p1[0]==point->leaf->boundary_p1[0] && point->leaf->boundary_p1[0] != top->boundary_p1[0] )
			|| (current->boundary_p1[1]==point->leaf->boundary_p1[1] && point->leaf->boundary_p1[1] != top->boundary_p1[1])
			|| (current->boundary_p2[0]==point->leaf->boundary_p2[0] && point->leaf->boundary_p2[0] != top->boundary_p2[0])
			|| (current->boundary_p2[1]==point->leaf->boundary_p2[1] && point->leaf->boundary_p2[1] != top->boundary_p2[1]) ){
		moveUp();
	}

	assert(inbox(point->x,current->boundary_p1,current->boundary_p2));
	_FindAllBoxNeighborsKist_iter(point->leaf,neighbors);
*/
/* TODO Test lines
	bool found;
	neighbors->MoveToTop();
	do{
		found = false;
		testkist->MoveToTop();
		do{
			if(testkist->getCurrent() == neighbors->getCurrent()){
				found=true;
				break;
			}
		}while(testkist->Down());
		assert(found);
	}while(neighbors->Down());
	delete testkist;

// *************************************************/
	return;
}
/**  \ingroup LowLevel
 * A recessive function that was used in FindAllBoxNeighborsKist().
*    It has been known to cause stack overflow. Use _FindAllBoxNeighborsKist_iter instead.
*/
void TreeStruct::_FindAllBoxNeighborsKist(Branch *leaf,TreeStruct::iterator &current,Kist<Point> * neighbors) const{

	if(  leaf->boundary_p1[0] <= (*current)->boundary_p2[0]
	  && leaf->boundary_p2[0] >= (*current)->boundary_p1[0]
	  && leaf->boundary_p1[1] <= (*current)->boundary_p2[1]
	  && leaf->boundary_p2[1] >= (*current)->boundary_p1[1]){

		if(current.atLeaf() ){
			assert((*current)->npoints <= Nbucket);
			//if(current->npoints == Nbucket){
			// What if number is > than Nbucket and it is not in a leaf
			if((*current)->number != leaf->number && (*current)->npoints > 0){
				neighbors->InsertAfterCurrent((*current)->points);
				neighbors->Down();
			}
			return;
		}

		if((*current)->child1 !=NULL){
			current.down(1);
			_FindAllBoxNeighborsKist(leaf,current,neighbors);
			current.up();
		}

		if((*current)->child2 !=NULL){
			current.down(2);
			_FindAllBoxNeighborsKist(leaf,current,neighbors);
			current.up();
		}
	}

	return;
}
/**  \ingroup LowLevel
 * Used in FindAllBoxNeighborsKist to walk tree for neighbors.
*/
void TreeStruct::_FindAllBoxNeighborsKist_iter(Branch *leaf,TreeStruct::iterator &current,Kist<Point> * neighbors) const{

	/** Iterative instead of recursive method for finding neighbors
	 *   the current must be preset so that leaf is within it.
	 *
	 *   The recursive routine _FindAllBoxNeighborsKist has caused stack overflow.
	 */

	bool allowDescent = true;
	long level = (*current)->level;

	neighbors->Empty();

	while(current.TreeWalkStep(allowDescent) && (*current)->level > level ){

		if(  leaf->boundary_p1[0] <= (*current)->boundary_p2[0]
		  && leaf->boundary_p2[0] >= (*current)->boundary_p1[0]
	      && leaf->boundary_p1[1] <= (*current)->boundary_p2[1]
	      && leaf->boundary_p2[1] >= (*current)->boundary_p1[1]){

			if( current.atLeaf() ){
				if((*current) != leaf && (*current)->npoints > 0){
					neighbors->InsertAfterCurrent((*current)->points);
					neighbors->Down();
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
void TreeStruct::PointsWithinEllipKist(
	const PosType* center         /// center of ellipse
	,float rmax                  /// major axis
	,float rmin                  /// minor axis
	,float posangle              /// position angle of major axis, smallest angle between the x-axis and the long axis
	,Kist<Point> * neighborkist  /// output neighbor kist, will be emptied if it contains anything on entry
	) const {
	unsigned long i,Ntmp;
	PosType *xtmp,x,y,cs,sn;

	neighborkist->Empty();

	if(rmax <=0.0 || rmin <= 0.0) return;

	// find point within a circle circumscribes the ellipse
	PointsWithinKist(center,rmax,neighborkist,false);

	cs = cos(posangle);
	sn = sin(posangle);
	Ntmp = neighborkist->Nunits();
	for(i=0,neighborkist->MoveToTop();i<Ntmp;++i,neighborkist->Down() ){
		xtmp = neighborkist->getCurrent()->x;
		x = xtmp[0]*cs - xtmp[1]*sn;
		y = xtmp[0]*sn + xtmp[1]*cs;
		if( (x*x/rmax/rmax) + (y*y/rmin/rmin) > 1) neighborkist->TakeOutCurrent();
	}
	return;
}
/** \ingroup ImageFindingL2
 * \brief Finds all points in tree that lie within rmax of the point ray[]
 *
 *   markpoints = 0  does not change in_image variable in any point, gives a list of neighbors
 *              = 1  makes in_image=YES for all points and their images in image, gives no list of neighbors
 *              = -1 makes in_image=NO for all points in image to reset, gives no list of neighbors
 *
 *   Returns the largest gridsize of the points within the circle.  Note that this is the gridsize stored
 *   in the point.  If finding points on the source plane the i_point->gridsize must be set to the same as the
 *   image point to get the largest gridsize on the image plane.
 */
PosType TreeStruct::PointsWithinKist(
		const PosType* center         /// center of circle
		,PosType rmax                  /// radius of circle
		,Kist<Point> * neighborkist  /// output neighbor kist, will be emptied if it contains anything on entry
		,short markpoints            /// see comment
		) const
{

	PosType maxgridsize;

	if(markpoints==0) neighborkist->Empty();

  TreeStruct::Globals globs;
  globs.ray[0] = globs.realray[0] = center[0];
  globs.ray[1] = globs.realray[1] = center[1];
  
	//moveTop();
  TreeStruct::iterator current(top);
	if( inbox(globs.ray,(*current)->boundary_p1,(*current)->boundary_p2) == 0 ){
		std::printf("Warning: in PointsWithinKist, ray is not inside the simulation box\n    should work in any case\n      ray= %e %e\n     boundary p1 = %e %e p2 = %e %e\n",globs.ray[0],globs.ray[1]
	   ,(*current)->boundary_p1[0],(*current)->boundary_p1[1]
	   ,(*current)->boundary_p2[0],(*current)->boundary_p2[1]);

		globs.ray[0]=MAX(globs.ray[0],(*current)->boundary_p1[0]);
		globs.ray[0]=MIN(globs.ray[0],(*current)->boundary_p2[0]);

		globs.ray[1]=MAX(globs.ray[1],(*current)->boundary_p1[1]);
		globs.ray[1]=MIN(globs.ray[1],(*current)->boundary_p2[1]);
	}
	globs.incell=1;

	maxgridsize = 0;
  _PointsWithinKist(current,&rmax,neighborkist,markpoints,&maxgridsize,globs);

	return maxgridsize;
}
/** \ingroup LowLevel
 * Used in PointsWithinKist() to walk tree.*/
void TreeStruct::_PointsWithinKist(TreeStruct::iterator &current,PosType *rmax,Kist<Point> * neighborkist
                                   ,short markpoints,PosType *maxgridsize,TreeStruct::Globals &globs) const{

  int i,j,incell2=1;
  PosType radius;
  short pass;

  if((*current)->npoints == 0) return;

  //std::printf("**************************************\nlevel %i\n",(*current)->level);
  //   std::printf("   %i incell=%i\n",(*current)->points->id,incell);

  if(globs.incell){  // not found cell yet

    if( inbox(globs.ray,(*current)->boundary_p1,(*current)->boundary_p2) ){

      // found the box small enough
    	if( Utilities::cutbox(globs.ray,(*current)->boundary_p1,(*current)->boundary_p2,*rmax)==1
    			|| current.atLeaf() ){
    		// whole box in circle or a leaf with ray in it

    	  globs.incell=0;

    	  // this sets ray back to real value once closest leaf box is found
    	  //assert((ray[0] == realray[0])*(ray[1] == realray[1]));
    	  //if( (ray[0]!=realray[0])*(ray[1]!=realray[1]) ){ std::printf("ray != realray _PointsWithinKist\n"); exit(0);}

    	  globs.ray[0] = globs.realray[0];
    	  globs.ray[1] = globs.realray[1];

    	  if((*current)->points != NULL) pointlist->current=(*current)->points;

    	  if( current.atLeaf() ){
    	   	  // if leaf calculate the distance to all the points in cell
    		  for(i=0;i<(*current)->npoints;++i){
    			  for(j=0,radius=0.0;j<2;++j) radius+=pow(pointlist->current->x[j] - globs.ray[j],2);
    			  if( radius < *rmax**rmax ){
       				  if(markpoints == 1){
       					  pointlist->current->in_image = YES;
      					  pointlist->current->image->in_image = YES;
      				  }else if(markpoints == -1){
      					  pointlist->current->in_image=NO;
     					  pointlist->current->image->in_image=NO;
     					  pointlist->current->surface_brightness = pointlist->current->image->surface_brightness = 0.0;
       				  }else if(markpoints == 0){
         				  neighborkist->InsertAfterCurrent(pointlist->current);
      				  }
       				  if(*maxgridsize < pointlist->current->gridsize) *maxgridsize = pointlist->current->gridsize;
    			  }
    			  MoveDownList(pointlist);
    		  }
    	  }else{ // put all of points in box into getCurrentKist(imagekist)
       		  for(i=0;i<(*current)->npoints;++i){
       			  if(markpoints == 1){
       				  pointlist->current->in_image=YES;
       				  pointlist->current->image->in_image=YES;
       			  }else if(markpoints == -1){
       				  pointlist->current->in_image=NO;
       				  pointlist->current->image->in_image=NO;
 					  pointlist->current->surface_brightness = pointlist->current->image->surface_brightness = 0.0;
      			  }else if(markpoints == 0){
       				  neighborkist->InsertAfterCurrent(pointlist->current);
 				  }
  				  if(*maxgridsize < pointlist->current->gridsize) *maxgridsize = pointlist->current->gridsize;

       			  MoveDownList(pointlist);
       		  }
    	  }

    	}else{ // keep going down the tree

     	  if((*current)->child1 !=NULL){
    		  current.down(1);
    		  _PointsWithinKist(current,rmax,neighborkist,markpoints,maxgridsize,globs);
          current.up();

    		  incell2 = globs.incell;
    	  }

    	  if((*current)->child2 !=NULL){
    		  current.down(2);
    		  _PointsWithinKist(current,rmax,neighborkist,markpoints,maxgridsize,globs);
    		  current.up();
    	  }

    	  // if ray found in second child go back to first to search for neighbors
    	  if( (incell2==1) && (globs.incell==0) ){
    		  if((*current)->child1 !=NULL){
            current.down(1);
    			  _PointsWithinKist(current,rmax,neighborkist,markpoints,maxgridsize,globs);
            current.up();
    		  }
    	  }
      }
    }  // not in the box

  }else{    // found cell

	  pass=Utilities::cutbox(globs.ray,(*current)->boundary_p1,(*current)->boundary_p2,*rmax);
	  // does radius cut into the box
	  if( pass ){

    	  if((*current)->points != NULL) pointlist->current=(*current)->points;

		  if( current.atLeaf()  ){  /* leaf case */

			  for(i=0;i<(*current)->npoints;++i){

				  for(j=0,radius=0.0;j<2;++j) radius+=pow(pointlist->current->x[j] - globs.ray[j],2);
				  if( radius < *rmax**rmax ){
					  if(markpoints==1){
						  pointlist->current->in_image=YES;
						  pointlist->current->image->in_image=YES;
					  }else if(markpoints==-1){
						  pointlist->current->in_image=NO;
						  pointlist->current->image->in_image=NO;
     					  pointlist->current->surface_brightness = pointlist->current->image->surface_brightness = 0.0;
					  }else if(markpoints==0){
						  neighborkist->InsertAfterCurrent(pointlist->current);
     				  }
      				  if(*maxgridsize < pointlist->current->gridsize) *maxgridsize = pointlist->current->gridsize;

				  }
				  MoveDownList(pointlist);
			  }
		  }else if(pass==1){ // whole box is inside radius

			  pointlist->current = (*current)->points;
			  for(i=0;i<(*current)->npoints;++i){

				  //assert( inbox(pointlist->current->x,(*current)->boundary_p1,(*current)->boundary_p2) );
				  //assert( *rmax**rmax >= (pow(pointlist->current->x[0] - ray[0],2) + pow(pointlist->current->x[1] - ray[1],2) ));

				  if(markpoints==1){
   					  pointlist->current->in_image=YES;
  					  pointlist->current->image->in_image=YES;
  				  }else if(markpoints==-1){
  					  pointlist->current->in_image=NO;
 					  pointlist->current->image->in_image=NO;
 					  pointlist->current->surface_brightness = pointlist->current->image->surface_brightness = 0.0;
   				  }else if(markpoints==0){
   					  neighborkist->InsertAfterCurrent(pointlist->current);
  				  }
  				  if(*maxgridsize < pointlist->current->gridsize) *maxgridsize = pointlist->current->gridsize;

				  MoveDownList(pointlist);
			  }
		  }else{
			  if((*current)->child1 !=NULL){
				  current.down(1);
				  _PointsWithinKist(current,rmax,neighborkist,markpoints,maxgridsize,globs);
          current.up();
			  }

			  if((*current)->child2 !=NULL){
				  current.down(2);
				  _PointsWithinKist(current,rmax,neighborkist,markpoints,maxgridsize,globs);
          current.up();
			  }
		  }

	  }
  }

  return;
}

/** \ingroup ImageFindingL2
 *
 * \brief Finds all points within a circle.  Much simpler, iterative algorithm.
 *
 */

void TreeStruct::PointsWithinKist_iter(const PosType* center,float rmin,float rmax,Kist<Point> * neighborkist) const {
	bool decend;
	unsigned long i;

  TreeStruct::iterator current(top);
	neighborkist->Empty();

	if(rmax <= 0.0) return;
	assert(rmax >= rmin);
	if(rmax <= rmin) return;

	if( CircleInBox(center,rmax,top->boundary_p1,top->boundary_p2) ){
		_FindLeaf(current,center,0);
		// Move up the tree till the whole circle is inside the box
		while(!CircleInBox(center,rmax,(*current)->boundary_p1,(*current)->boundary_p2) && current.up());
	}

	Branch *top = (*current);

	if(rmin <= 0.0){

		while( (*current) != top->brother){

			decend = true;

			if(BoxInCircle(center,rmax,(*current)->boundary_p1,(*current)->boundary_p2)  // box is all inside outer circle
			){

				decend = false;
				if((*current)->points != NULL) pointlist->current = (*current)->points;
				for(i=0;i<(*current)->npoints;++i){
					neighborkist->InsertAfterCurrent(pointlist->current);
					MoveDownList(pointlist);
				}

			}else if(Utilities::cutbox(center,(*current)->boundary_p1,(*current)->boundary_p2,rmax) == 0  // box is all outside outer circle
			){

				decend = false;

			}else if(current.atLeaf()){      // box is a leaf that intersects the circle

				if((*current)->points != NULL) pointlist->current = (*current)->points;
				for(i=0;i<(*current)->npoints;++i){
					if(rmax*rmax >= pow(pointlist->current->x[0] - center[0],2) + pow(pointlist->current->x[1] - center[1],2) )
						neighborkist->InsertAfterCurrent(pointlist->current);
					MoveDownList(pointlist);
				}
			}

			if(!(current.TreeWalkStep(decend))) break;
		}

	}else{  // rmin > 0
		PosType r2;

		while(*current != top->brother){

			decend = true;

			if(current.atLeaf()){      // box is a leaf that intersects the circle

				if((*current)->points != NULL){
					pointlist->current = (*current)->points;
				}

				for(i=0;i<(*current)->npoints;++i){

					r2 = pow(pointlist->current->x[0] - center[0],2) + pow(pointlist->current->x[1] - center[1],2);
					if(rmax*rmax >= r2 && rmin*rmin <= r2){
						neighborkist->InsertAfterCurrent(pointlist->current);
					}
					MoveDownList(pointlist);
				}

			}else if(BoxInCircle(center,rmax,(*current)->boundary_p1,(*current)->boundary_p2)  // box is all inside outer circle
					&& Utilities::cutbox(center,(*current)->boundary_p1,(*current)->boundary_p2,rmin) == 0  // box is all outside inner circle
			){

				decend = false;
				if((*current)->points != NULL){
					pointlist->current = (*current)->points;
				}

				for(i=0;i<(*current)->npoints;++i){
					neighborkist->InsertAfterCurrent(pointlist->current);
					MoveDownList(pointlist);
				}

			}else if(Utilities::cutbox(center,(*current)->boundary_p1,(*current)->boundary_p2,rmax) == 0  // box is all outside outer circle
					|| BoxInCircle(center,rmin,(*current)->boundary_p1,(*current)->boundary_p2) // box is all inside inner circle
			){

				decend = false;

			}
			if(!(current.TreeWalkStep(decend))) break;
		}

	}
}

/** \ingroup ImageFundingL2
 *
 * \brief Finds nearest neighbor points to ray.
 *
 *    This is a kludge that relies on NearestNeighbor which uses a List and translates
 *    the list to a kist.  Could be rewritten.
 * 
 *    Warning: The number of neighbor points in neighborkist will be less than Nneighbors when
 *             the number of points in the tree is less than Nneighbors
 *
Point * TreeStruct::NearestNeighborKist(const PosType* center,int Nneighbors,Kist<Point> * neighborkist){
	ListHndl neighborlist = NewList();
	Point *point = 0;
	unsigned long i;

	//TODO: BEN Make this better!  NearestNeighbor() should be replaced.
	point = NearestNeighbor(center,Nneighbors,neighborlist,0);

	// convert from point array to exported point kist
	neighborkist->Empty();
	MoveToTopList(neighborlist);
	for(i = 0; i < neighborlist->Npoints ;++i){
		neighborkist->InsertAfterCurrent(neighborlist->current->image->image);
		neighborkist->Down();
		MoveDownList(neighborlist);
	}

	assert(neighborlist->Npoints == neighborkist->Nunits());

	EmptyList(neighborlist);
	free(neighborlist);

	return point;
}
*/