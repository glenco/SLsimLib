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

	Ngrid_init = N1d;
	Ngrid_block = 3;  // never been tested with anything other than 3
	

	i_points = NewPointArray(Ngrid_init*Ngrid_init,true);
	xygridpoints(i_points,range,center,Ngrid_init,0);
	s_points=LinkToSourcePoints(i_points,Ngrid_init*Ngrid_init);
	lens->rayshooterInternal(Ngrid_init*Ngrid_init,i_points,false);
	// Build trees
	i_tree = BuildTree(i_points,Ngrid_init*Ngrid_init);
	s_tree = BuildTree(s_points,Ngrid_init*Ngrid_init);  // make tree on source plane a area splitting tree

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
	i_points = NewPointArray(Ngrid_init*Ngrid_init,true);
	xygridpoints(i_points,range,center,Ngrid_init,0);
	s_points=LinkToSourcePoints(i_points,Ngrid_init*Ngrid_init);
	lens->rayshooterInternal(Ngrid_init*Ngrid_init,i_points,false);

	// need to resize root of source tree.  It can change in size
	s_tree->top->boundary_p1[0]=s_points[0].x[0]; s_tree->top->boundary_p1[1]=s_points[0].x[1];
	s_tree->top->boundary_p2[0]=s_points[0].x[0]; s_tree->top->boundary_p2[1]=s_points[0].x[1];

	for(i=0;i<Ngrid_init*Ngrid_init;++i){

	    /* find X boundary */
		if(s_points[i].x[0] < s_tree->top->boundary_p1[0] ) s_tree->top->boundary_p1[0]=s_points[i].x[0];
	    if(s_points[i].x[0] > s_tree->top->boundary_p2[0] ) s_tree->top->boundary_p2[0]=s_points[i].x[0];

	    /* find Y boundary */
	    if(s_points[i].x[1] < s_tree->top->boundary_p1[1] ) s_tree->top->boundary_p1[1]=s_points[i].x[1];
	    if(s_points[i].x[1] > s_tree->top->boundary_p2[1] ) s_tree->top->boundary_p2[1]=s_points[i].x[1];
	  }

	  // a little extra room for future points
	  s_tree->top->boundary_p1[0] -=  range/Ngrid_init;
	  s_tree->top->boundary_p1[1] -=  range/Ngrid_init;
	  s_tree->top->boundary_p2[0] +=  range/Ngrid_init;
	  s_tree->top->boundary_p2[1] +=  range/Ngrid_init;

	  s_tree->top->center[0] = (s_tree->top->boundary_p1[0]+s_tree->top->boundary_p2[0])/2;
	  s_tree->top->center[1] = (s_tree->top->boundary_p1[1]+s_tree->top->boundary_p2[1])/2;

	// fill trees
	FillTree(i_tree,i_points,Ngrid_init*Ngrid_init);
	FillTree(s_tree,s_points,Ngrid_init*Ngrid_init);

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
	for(unsigned long i=0;i<s_tree->pointlist->Npoints;++i,MoveDownList(s_tree->pointlist)){
		//y[0] = s_tree->pointlist->current->x[0]; - source->getX()[0];
		//y[1] = s_tree->pointlist->current->x[1]; - source->getX()[1];
		tmp = source->SurfaceBrightness(s_tree->pointlist->current->x);
		s_tree->pointlist->current->surface_brightness = s_tree->pointlist->current->image->surface_brightness
				= tmp;
		total += tmp;//*pow( s_tree->pointlist->current->gridsize,2);
		assert(s_tree->pointlist->current->surface_brightness >= 0.0);
		s_tree->pointlist->current->in_image = s_tree->pointlist->current->image->in_image
				= FALSE;
	}

	return total;
}

/** \ingroup ImageFinding
 * \brief Returns number of points on image plane.
 */
unsigned long Grid::getNumberOfPoints(){
	assert(i_tree->top->npoints == s_tree->top->npoints);
	assert(i_tree->top->npoints == i_tree->pointlist->Npoints);
	assert(s_tree->top->npoints == s_tree->pointlist->Npoints);

	return i_tree->top->npoints;
}

/**  \ingroup ImageFindingL2
 *
 * \brief Fundamental function used to divide a leaf in the tree into nine subcells.
 *
 * Source and image points are created, linked, shot and added to the trees.  The leaf
 * pointers of the points including the input are assigned.
 *
 * If some of the of the points are outside the original grid they will not be added in
 * which case THERE WILL BE LESS THEN Ngrid*Ngrid-1 points added.  The true number will
 * always be result->head.
 *
 * Returns a pointer to the list of image points that have been added.  This array can then be
 * used for calculating the surface brightness or marking them as in the image.
 */

Point * Grid::RefineLeaf(LensHndl lens,Point *point,bool kappa_off){
	//Point * RefineLeaf(LensHndl lens,TreeHndl i_tree,TreeHndl s_tree,Point *point,int Ngrid,bool kappa_off){

	Point *i_points = NewPointArray(Ngrid_block*Ngrid_block-1,true);
	Point *s_points;
	int Nout,kk;

	assert(point->leaf->child1 == NULL && point->leaf->child2 == NULL);
	assert(point->gridsize > pow(10.,-DBL_DIG) ); // If cells are too small they will cause problems.

	point->leaf->refined = true;
	xygridpoints(i_points,point->gridsize*(Ngrid_block-1)/Ngrid_block
	      ,point->x,Ngrid_block,1);
	point->gridsize /= Ngrid_block;
	point->image->gridsize /= Ngrid_block;

	// take out points that are outside of original grid
	Nout = 0;
	if( (point->x[0] == i_tree->top->boundary_p1[0]) || (point->x[0] == i_tree->top->boundary_p2[0])
			|| (point->x[1] == i_tree->top->boundary_p1[1]) || (point->x[1] == i_tree->top->boundary_p2[1]) ){

		  // remove the points that are outside initial image grid
		  for(kk=0,Nout=0;kk < (Ngrid_block*Ngrid_block-1);++kk){
			  if( !inbox(i_points[kk - Nout].x,i_tree->top->boundary_p1,i_tree->top->boundary_p2) ){
				  SwapPointsInArray(&i_points[kk - Nout],&i_points[Ngrid_block*Ngrid_block - 2 - Nout]);
				  ++Nout;
			  }
		  }
		  assert(Nout > 0);
	}

	//if(Nout > 0) i_points = AddPointToArray(i_points,Ngrid_block*Ngrid_block-1-Nout,Ngrid_block*Ngrid_block-1);

	int  Ntemp = Ngrid_block*Ngrid_block-1-Nout;

	s_points = LinkToSourcePoints(i_points,Ntemp);
	lens->rayshooterInternal(Ntemp,i_points,kappa_off);

	// remove the points that are outside initial source grid
	for(kk=0,Nout=0;kk < Ntemp;++kk){
		assert(s_points[kk - Nout].x[0] == s_points[kk - Nout].x[0]);
		if( !inbox(s_points[kk - Nout].x,s_tree->top->boundary_p1,s_tree->top->boundary_p2) ){
			SwapPointsInArray(&i_points[kk - Nout],&i_points[Ntemp - 1 - Nout]);
			SwapPointsInArray(&s_points[kk - Nout],&s_points[Ntemp - 1 - Nout]);
			++Nout;
		}
	}

	assert(i_points->head == Ngrid_block*Ngrid_block-1);
	assert(s_points->head == Ntemp);

	// free memory of points that where outside image and source regions
	Nout = Ngrid_block*Ngrid_block - 1 - Ntemp + Nout;
	if(Nout > 0){
		//i_points = AddPointToArray(i_points,Ngrid_block*Ngrid_block-1-Nout,Ngrid_block*Ngrid_block-1);
		//s_points = AddPointToArray(s_points,Ngrid_block*Ngrid_block-1-Nout,Ntemp);

		i_points = AddPointToArray(i_points,Ngrid_block*Ngrid_block-1-Nout,i_points->head);
		s_points = AddPointToArray(s_points,Ngrid_block*Ngrid_block-1-Nout,s_points->head);
	}
	assert(i_points->head == s_points->head);

	//*** these could be mode more efficient by starting at the current in tree
	//TODO test line
	AddPointsToTree(i_tree,i_points,i_points->head);
	AddPointsToTree(s_tree,s_points,s_points->head);

	assert(s_points->head > 0);
	//AddPointsToTree(i_tree,i_points,Ngrid_block*Ngrid_block-1-Nout);
	//AddPointsToTree(s_tree,s_points,Ngrid_block*Ngrid_block-1-Nout);

	// re-assign leaf of point that was to be refined
	assert(inbox(point->x,i_tree->top->boundary_p1,i_tree->top->boundary_p2));
	i_tree->current = point->leaf;
	//TODO This line should not be necessary!! It is repairing the leaf that has been assigned incorrectly somewhere
	if(!inbox(point->x,i_tree->current->boundary_p1,i_tree->current->boundary_p2) ) moveTop(i_tree);
	_FindLeaf(i_tree,point->x,0);
	point->leaf = i_tree->current;

	assert(inbox(point->image->x,s_tree->top->boundary_p1,s_tree->top->boundary_p2));
	s_tree->current = point->image->leaf;
	//TODO This line should not be necessary!! It is repairing the leaf that has been assigned incorrectly somewhere
	if(!inbox(point->image->x,s_tree->current->boundary_p1,s_tree->current->boundary_p2) ) moveTop(s_tree);
	_FindLeaf(s_tree,point->image->x,0);
	point->image->leaf = s_tree->current;

	assert(point->leaf->child1 == NULL && point->leaf->child2 == NULL);

	return i_points;
}

