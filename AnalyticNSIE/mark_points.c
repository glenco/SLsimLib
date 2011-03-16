/*
 * mark_points.c
 *
 *  Created on: Oct 24, 2010
 *      Author: bmetcalf
 */

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <Tree.h>
#include"../AnalyticNSIE/analytic_lens.h"

/* MarkPoints sets the point.in_image to True
 * if the source point is within lens->source_r of lens->source_x
 * AND (if sbcut == True) the surface brightness at that point is > 0.
 * Both the source point and its image point are marked.
 *
 * Setting sb_cut = False speeds the code up for sources without sharp
 * edges.
 *
 * if invert == -1 the same points are unmarked
 *
 * Note: The whole source must be within lens->r_source.
 * Note: The source parameters must not be changed before inverting
 *       to unmark the points.
 */
static short incell,invertg=0;

void MarkPoints(TreeHndl s_tree,AnaLens *lens,Boolean sb_cut,short invert){

	assert(s_tree); assert(lens);

	incell = 1;
	invertg = invert;

	moveTop( s_tree);
	_MarkPoints( s_tree,lens,&sb_cut);

  return;
}

void _MarkPoints(TreeHndl s_tree,AnaLens *lens,Boolean *sbcut){

  int i,incell2=1;
  short pass;

  if(incell){  // not found cell yet

    if( inbox(lens->source_x, s_tree->current->boundery_p1, s_tree->current->boundery_p2) ){

      // found the box small enough
    	if( cutbox(lens->source_x, s_tree->current->boundery_p1, s_tree->current->boundery_p2,lens->source_r) == 1
    			|| ( s_tree->current->child1 == NULL)*( s_tree->current->child2 == NULL) ){

    		// whole box in circle or a leaf with ray in it
    		incell=0;

    		// if leaf calculate the distance to all the points in cell
       		 s_tree->pointlist->current =  s_tree->current->points;
       		for(i=0;i< s_tree->current->npoints;++i){
       			if(invertg == -1)
       				s_tree->pointlist->current->in_image = s_tree->pointlist->current->image->in_image = 0;
       			else
       				s_tree->pointlist->current->in_image = s_tree->pointlist->current->image->in_image
    					= InSource(s_tree->pointlist->current->x,lens,*sbcut);

    			MoveDownList( s_tree->pointlist);
    		}

    	}else{ // keep going down the  s_tree

    		if( s_tree->current->child1 !=NULL){
    			moveToChild( s_tree,1);
    			_MarkPoints( s_tree,lens,sbcut);
    			moveUp( s_tree);

    			incell2=incell;
    		}

    		if( s_tree->current->child2 !=NULL){
    			moveToChild( s_tree,2);
    			_MarkPoints( s_tree,lens,sbcut);
    			moveUp( s_tree);
    		}

    		// if ray found in second child go back to first to search for neighbors
    		if( (incell2==1) && (incell==0) ){
    			if( s_tree->current->child1 !=NULL){
    				moveToChild( s_tree,1);
    				_MarkPoints( s_tree,lens,sbcut);
    				moveUp( s_tree);
    			}
    		}

    	}

    }  // not in the box

  }else{    // found cell

	  pass=cutbox(lens->source_x, s_tree->current->boundery_p1, s_tree->current->boundery_p2,lens->source_r);
	  // does radius cut into the box
	  if( pass ){

		  if( ( s_tree->current->child1 == NULL)*( s_tree->current->child2 == NULL) || pass == 1 ){  /* leaf case */
			   s_tree->pointlist->current= s_tree->current->points;
			  for(i=0;i< s_tree->current->npoints;++i){
	       			if(invertg == -1)
	       				s_tree->pointlist->current->in_image = s_tree->pointlist->current->image->in_image = 0;
	       			else
	       				s_tree->pointlist->current->in_image = s_tree->pointlist->current->image->in_image
	    					= InSource(s_tree->pointlist->current->x,lens,*sbcut);
		   			MoveDownList( s_tree->pointlist);
			  }
		  }else{
			  if( s_tree->current->child1 !=NULL){
				  moveToChild( s_tree,1);
				  _MarkPoints( s_tree,lens,sbcut);
				  moveUp( s_tree);
			  }

			  if( s_tree->current->child2 !=NULL){
				  moveToChild( s_tree,2);
				  _MarkPoints( s_tree,lens,sbcut);
				  moveUp( s_tree);
			  }
		  }

	  }
  }
  return;
}

Boolean InSource(double *ray,AnaLens *lens,Boolean surfacebright){
	double r[2];

	r[0] = lens->source_x[0] - ray[0];
	r[1] = lens->source_x[1] - ray[1];

	if(r[0]*r[0]+r[1]*r[1] < lens->source_r*lens->source_r){
		if(surfacebright && lens->source_sb_func(r) <= 0.0) return False;
		else return True;
	}

	return False;
}
