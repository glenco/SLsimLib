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
 * Note: Deflection solver must be specified before creating a Grid.
 */
Grid::Grid(
		LensHndl lens      /// lens model for initializing grid
		,int N1d           /// Initial number of grid points in each dimension.
		,double center[2]  /// Center of grid.
		,double range      /// Full width of grid in whatever units will be used.
		 ){

	Point *i_points,*s_points;

	assert(N1d > 0);
	assert(range > 0);

	if(N1d <= 0){ERROR_MESSAGE(); std::cout << "cannot make Grid with no points" << std::endl; exit(1);}
	if(range <= 0){ERROR_MESSAGE(); std::cout << "cannot make Grid with no range" << std::endl; exit(1);}

	Ngrid = N1d;
	Ngrid_block = 3;  // never been tested with anything other than 3
	

	i_points = NewPointArray(Ngrid*Ngrid,true);
	xygridpoints(i_points,range,center,Ngrid,0);
	s_points=LinkToSourcePoints(i_points,Ngrid*Ngrid);
	lens->rayshooterInternal(Ngrid*Ngrid,i_points,false);
	// Build trees
	i_tree = BuildTree(i_points,Ngrid*Ngrid);
	s_tree = BuildTree(s_points,Ngrid*Ngrid);  // make tree on source plane a area splitting tree

	trashkist = new Kist;
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
	
	delete trashkist;
	
	return;
}

/** \ingroup ImageFinding
 *  \brief Reinitializes the grid so that it is back to the original coarse grid, but if
 *  the lens has changed the source positions will be updated.
 */
void Grid::ReInitializeGrid(LensHndl lens){

	Point *i_points,*s_points;
	double range,center[2];
	unsigned long i;

	range = i_tree->top->boundary_p2[0] - i_tree->top->boundary_p1[0];
	center[0] = (i_tree->top->boundary_p2[0] + i_tree->top->boundary_p1[0])/2;
	center[1] = (i_tree->top->boundary_p2[1] + i_tree->top->boundary_p1[1])/2;

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

	// need to resize root of source tree.  It can change in size
	s_tree->top->boundary_p1[0]=s_points[0].x[0]; s_tree->top->boundary_p1[1]=s_points[0].x[1];
	s_tree->top->boundary_p2[0]=s_points[0].x[0]; s_tree->top->boundary_p2[1]=s_points[0].x[1];

	for(i=0;i<Ngrid*Ngrid;++i){

	    /* find X boundary */
		if(s_points[i].x[0] < s_tree->top->boundary_p1[0] ) s_tree->top->boundary_p1[0]=s_points[i].x[0];
	    if(s_points[i].x[0] > s_tree->top->boundary_p2[0] ) s_tree->top->boundary_p2[0]=s_points[i].x[0];

	    /* find Y boundary */
	    if(s_points[i].x[1] < s_tree->top->boundary_p1[1] ) s_tree->top->boundary_p1[1]=s_points[i].x[1];
	    if(s_points[i].x[1] > s_tree->top->boundary_p2[1] ) s_tree->top->boundary_p2[1]=s_points[i].x[1];
	  }

	  // a little extra room for future points
	  s_tree->top->boundary_p1[0] -=  range/Ngrid;
	  s_tree->top->boundary_p1[1] -=  range/Ngrid;
	  s_tree->top->boundary_p2[0] +=  range/Ngrid;
	  s_tree->top->boundary_p2[1] +=  range/Ngrid;

	  s_tree->top->center[0] = (s_tree->top->boundary_p1[0]+s_tree->top->boundary_p2[0])/2;
	  s_tree->top->center[1] = (s_tree->top->boundary_p1[1]+s_tree->top->boundary_p2[1])/2;

	// fill trees
	FillTree(i_tree,i_points,Ngrid*Ngrid);
	FillTree(s_tree,s_points,Ngrid*Ngrid);

	return;
}

/** \ingroup ImageFinding
 *
 * \brief DOES NOT WORK YET !!!!
 */

/** \ingroup ImageFinding
 * \brief Recalculate surface brightness at every point without changing the positions of the grid or any lens properties.
 *
 *  Recalculate the surface brightness at all points on the grid.
 * This is useful when changing the source model while preserving
 * changes in the grid.
 * Both i_tree and s_tree are both changed although only s_tree shows up here.
 *
 * returns the sum of the surface brightnesses
 */
double Grid::RefreshSurfaceBrightnesses(SourceHndl source){
	double y[2],total=0,tmp;

	MoveToTopList(s_tree->pointlist);
	do{
		y[0] = s_tree->pointlist->current->x[0] - source->source_x[0];
		y[1] = s_tree->pointlist->current->x[1] - source->source_x[1];
		tmp = source->source_sb_func(y);
		s_tree->pointlist->current->surface_brightness = s_tree->pointlist->current->image->surface_brightness
				= tmp;
		total += tmp;
		assert(s_tree->pointlist->current->surface_brightness >= 0.0);
		s_tree->pointlist->current->in_image = s_tree->pointlist->current->image->in_image
				= FALSE;
	}while( MoveDownList(s_tree->pointlist) );

	return total;
}

/** \ingroup ImageFinding
 * \brief Returns number of points on image plane.
 */
unsigned long Grid::NumberOfPoints(){
	assert(i_tree->top->npoints == s_tree->top->npoints);
	assert(i_tree->top->npoints == i_tree->pointlist->Npoints);
	assert(s_tree->top->npoints == s_tree->pointlist->Npoints);

	return i_tree->top->npoints;
}
