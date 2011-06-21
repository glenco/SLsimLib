/*
 * grid_initialization.c
 *
 *  Created on: Apr 12, 2011
 *      Author: bmetcalf
 */
#include <stdlib.h>
#include <Tree.h>

void init_grid(GridHndl grid,unsigned long Ngrid,double *center,double range){
	Point *i_points,*s_points;

	i_points = NewPointArray(Ngrid*Ngrid,True);
	xygridpoints(i_points,range,center,Ngrid,0);
	s_points=LinkToSourcePoints(i_points,Ngrid*Ngrid);
	grid->i_tree=NULL;
	rayshooterInternal(Ngrid*Ngrid,i_points,grid->i_tree,False);

	// build trees
	grid->i_tree=BuildTree(i_points,Ngrid*Ngrid);
	grid->s_tree=BuildTree(s_points,Ngrid*Ngrid);

	grid->initialized = True;

	return;
}

void reinit_grid(GridHndl grid,unsigned long Ngrid,double *center,double range){
	Point *i_points,*s_points;

	if(grid->initialized){
	  // free old tree to speed up image finding
	  emptyTree(grid->i_tree);
	  emptyTree(grid->s_tree);

	  // build new initale grid
	  i_points = NewPointArray(Ngrid*Ngrid,True);
	  xygridpoints(i_points,range,center,Ngrid,0);
	  s_points = LinkToSourcePoints(i_points,Ngrid*Ngrid);
	  rayshooterInternal(Ngrid*Ngrid,i_points,grid->i_tree,False);
	  // fill trees
	  FillTree(grid->i_tree,i_points,Ngrid*Ngrid);
	  FillTree(grid->s_tree,s_points,Ngrid*Ngrid);
	}else{
		init_grid(grid,Ngrid,center,range);
	}

	return;
}
