/*
 * grid_maintenance.c
 *
 *  Created on: Oct 12, 2011
 *      Author: bmetcalf
 */
/*#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <point.h>
#include <Tree.h>
#include <tree_maintenance.h>
#include <grid_maintenance.h>
*/

#include <slsimlib.h>

/** \ingroup Constructor
 * \brief Constructor for initializing grid.
 *
 * Note: Defelection solver must be specified before creating a Grid.
 */
GridHndl NewGrid(LensHndl lens, int Ngrid,double center[2],double range){
	GridHndl grid = (Grid *)malloc(sizeof(Grid));
	Point *i_points,*s_points;

	grid->Ngrid = Ngrid;

	i_points = NewPointArray(Ngrid*Ngrid,true);
	xygridpoints(i_points,range,center,Ngrid,0);
	s_points=LinkToSourcePoints(i_points,Ngrid*Ngrid);
	lens->rayshooterInternal(Ngrid*Ngrid,i_points,true);
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
void ReInitializeGrid(LensHndl lens, GridHndl grid){

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
	lens->rayshooterInternal(Ngrid*Ngrid,i_points,true);
	// fill trees
	FillTree(grid->i_tree,i_points,Ngrid*Ngrid);
	FillTree(grid->s_tree,s_points,Ngrid*Ngrid);

	return;
}
/** \ingroup ImageFinding
 *
 * \brief DOES NOT WORK YET !!!!
 */
void TrimGrid(GridHndl grid,double highestres,bool useSB){

	PruneTrees(grid->i_tree,grid->s_tree,highestres,useSB);
	return;
}

/** \ingroup ImageFinding
 * \brief Recalculate surface brightness at every point without changing the positions of the grid or any lens properties.
 *
 *  Recalculate the surface brightness at all points on the grid.
 * This is useful when changing the source model while preserving
 * changes in the grid.
 * Both i_tree and s_tree are both changed although only s_tree shows up here.
 */
void RefreshSurfaceBrightnesses(GridHndl grid,ModelHndl model){
	double y[2];

	CosmoHndl cosmo;
	SourceHndl source;

	cosmo = model->cosmo;
	source = model->source;

	MoveToTopList(grid->s_tree->pointlist);
	do{
		y[0] = grid->s_tree->pointlist->current->x[0] - source->source_x[0];
		y[1] = grid->s_tree->pointlist->current->x[1] - source->source_x[1];
		grid->s_tree->pointlist->current->surface_brightness = grid->s_tree->pointlist->current->image->surface_brightness
				= (model->*source_sb_func)(y);
		assert(grid->s_tree->pointlist->current->surface_brightness >= 0.0);
		grid->s_tree->pointlist->current->in_image = grid->s_tree->pointlist->current->image->in_image
				= false;
	}while( MoveDownList(grid->s_tree->pointlist) );

	return;
}

/** \ingroup ImageFinding
 * \brief Returns number of points on image plane.
 */
unsigned long NumberOfPoints(GridHndl grid){
	assert(grid->i_tree->top->npoints == grid->s_tree->top->npoints);
	return grid->i_tree->top->npoints;
}
