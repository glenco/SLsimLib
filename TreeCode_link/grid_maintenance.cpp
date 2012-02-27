/*
 * grid_maintenance.c
 *
 *  Created on: Oct 12, 2011
 *      Author: bmetcalf
 */

#include <slsimlib.h>

/** \ingroup Constructor
 * \brief Constructor for initializing grid.
 *
 * Note: Defelection solver must be specified before creating a Grid.
 */
Grid::Grid(
		LensHndl lens      /// lens model for initializing grid
		,int N1d           /// Initial number of grid points in each dimension.
		,double center[2]  /// Center of grid.
		,double range      /// Full width of grid in whatever units will be used.
		 ){

	Point *i_points,*s_points;

	Ngrid = N1d;

	i_points = NewPointArray(Ngrid*Ngrid,true);
	xygridpoints(i_points,range,center,Ngrid,0);
	s_points=LinkToSourcePoints(i_points,Ngrid*Ngrid);
	lens->rayshooterInternal(Ngrid*Ngrid,i_points,false);
	// Build trees
	i_tree = BuildTree(i_points,Ngrid*Ngrid);
	s_tree = BuildTree(s_points,Ngrid*Ngrid);
}
/*
GridHndl NewGrid(LensHndl lens, int Ngrid,double center[2],double range){
	GridHndl grid = (Grid *)malloc(sizeof(Grid));
	Point *i_points,*s_points;

	grid->Ngrid = Ngrid;

	i_points = NewPointArray(Ngrid*Ngrid,true);
	xygridpoints(i_points,range,center,Ngrid,0);
	s_points=LinkToSourcePoints(i_points,Ngrid*Ngrid);
	lens->rayshooterInternal(Ngrid*Ngrid,i_points,false);
	// Build trees
	grid->i_tree = BuildTree(i_points,Ngrid*Ngrid);
	grid->s_tree = BuildTree(s_points,Ngrid*Ngrid);

	return grid;
}
*/
/** \ingroup Constructor
 * \brief Destructor for a Grid.  Frees all memory.
 */
Grid::~Grid(){
	freeTree(i_tree);
	freeTree(s_tree);
	return;
}
/*
void FreeGrid(GridHndl grid){
	freeTree(grid->i_tree);
	freeTree(grid->s_tree);
	free(grid);

	return;
}
*/
/** \ingroup ImageFinding
 *  \brief Reinitializes the grid so that it is back to the original coarse grid, but if
 *  the lens has changed the source positions will be updated.
 */
void Grid::ReInitializeGrid(LensHndl lens){

	Point *i_points,*s_points;
	double range,center[2];
	int Ngrid;

	range = i_tree->top->boundary_p2[0] - i_tree->top->boundary_p1[0];
	center[0] = (i_tree->top->boundary_p2[0] + i_tree->top->boundary_p1[0])/2;
	center[1] = (i_tree->top->boundary_p2[1] + i_tree->top->boundary_p1[1])/2;
	Ngrid = Ngrid;

	//////////////////////////////
	  // redo grid with stars in it
	  // free old tree to speed up image finding
	emptyTree(i_tree);
	emptyTree(s_tree);


	// build new initial grid
	i_points = NewPointArray(Ngrid*Ngrid,true);
	xygridpoints(i_points,range,center,Ngrid,0);
	s_points=LinkToSourcePoints(i_points,Ngrid*Ngrid);
	lens->rayshooterInternal(Ngrid*Ngrid,i_points,false);
	// fill trees
	FillTree(i_tree,i_points,Ngrid*Ngrid);
	FillTree(s_tree,s_points,Ngrid*Ngrid);
}

/** \ingroup ImageFinding
 *
 * \brief DOES NOT WORK YET !!!!
 */
void Grid::TrimGrid(double highestres,bool useSB){

	PruneTrees(i_tree,s_tree,highestres,useSB);
}

/** \ingroup ImageFinding
 * \brief Recalculate surface brightness at every point without changing the positions of the grid or any lens properties.
 *
 *  Recalculate the surface brightness at all points on the grid.
 * This is useful when changing the source model while preserving
 * changes in the grid.
 * Both i_tree and s_tree are both changed although only s_tree shows up here.
 */
void Grid::RefreshSurfaceBrightnesses(SourceHndl source){
	double y[2];

	MoveToTopList(s_tree->pointlist);
	do{
		y[0] = s_tree->pointlist->current->x[0] - source->source_x[0];
		y[1] = s_tree->pointlist->current->x[1] - source->source_x[1];
		s_tree->pointlist->current->surface_brightness = s_tree->pointlist->current->image->surface_brightness
				= (source->source_sb_func)(y);
		assert(s_tree->pointlist->current->surface_brightness >= 0.0);
		s_tree->pointlist->current->in_image = s_tree->pointlist->current->image->in_image
				= false;
	}while( MoveDownList(s_tree->pointlist) );

}

/** \ingroup ImageFinding
 * \brief Returns number of points on image plane.
 */
unsigned long Grid::NumberOfPoints(){
	assert(i_tree->top->npoints == s_tree->top->npoints);
	return i_tree->top->npoints;
}
