/*
 * grid_maintenance.c
 *
 *  Created on: Oct 12, 2011
 *      Author: bmetcalf
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <point.h>
#include <Tree.h>
#include <tree_maintenance.h>
#include <grid_maintenance.h>

/** \ingroup Constructor
 * \brief Constructor for initializing grid.
 *
 * Note: Defelection solver must be specified before creating a Grid.
 */
GridHndl NewGrid(int Ngrid,double center[2],double range){
	GridHndl grid = (Grid *)malloc(sizeof(Grid));
	Point *i_points,*s_points;

	grid->Ngrid = Ngrid;

	i_points = NewPointArray(Ngrid*Ngrid,true);
	xygridpoints(i_points,range,center,Ngrid,0);
	s_points=LinkToSourcePoints(i_points,Ngrid*Ngrid);
	rayshooterInternal(Ngrid*Ngrid,i_points,true);
	// Build trees
	grid->i_tree = BuildTree(i_points,Ngrid*Ngrid);
	grid->s_tree = BuildTree(s_points,Ngrid*Ngrid);

	return grid;
}

/** \ingroup Constructor
 * \brief Destructor for a Grid.  Frees all memory.
 */
void FreeGrid(GridHndl grid){
	freeTree(grid->i_tree);
	freeTree(grid->s_tree);
	free(grid);

	return;
}
/** \ingroup ImageFinding
 *  \brief Reinitializes the grid so that it is back to the original coarse grid, but if
 *  the lens has changed the source positions will be updated.
 */
void ReInitalizeGrid(GridHndl grid){

	Point *i_points,*s_points;
	double range,center[2];
	int Ngrid;

	range = grid->i_tree->top->boundary_p2[0] - grid->i_tree->top->boundary_p1[0];
	center[0] = (grid->i_tree->top->boundary_p2[0] + grid->i_tree->top->boundary_p1[0])/2;
	center[1] = (grid->i_tree->top->boundary_p2[1] + grid->i_tree->top->boundary_p1[1])/2;
	Ngrid = grid->Ngrid;

	//////////////////////////////
	  // redo grid with stars in it
	  // free old tree to speed up image finding
	emptyTree(grid->i_tree);
	emptyTree(grid->s_tree);


	// build new initial grid
	i_points = NewPointArray(Ngrid*Ngrid,true);
	xygridpoints(i_points,range,center,Ngrid,0);
	s_points=LinkToSourcePoints(i_points,Ngrid*Ngrid);
	rayshooterInternal(Ngrid*Ngrid,i_points,true);
	// fill trees
	FillTree(grid->i_tree,i_points,Ngrid*Ngrid);
	FillTree(grid->s_tree,s_points,Ngrid*Ngrid);

	return;
}
/** \ingroup ImageFinding
 *
 */
void TrimGrid(GridHndl grid,double highestres,bool useSB){

	return;
}
