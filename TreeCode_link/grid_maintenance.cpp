/*
 * grid_maintenance.c
 *
 *  Created on: Oct 12, 2011
 *      Author: bmetcalf
 */

#include "slsimlib.h"

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
	neighbors = new Kist;
	maglimit = 1.0e-4;
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
	delete neighbors;
	
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

	/*for(i=0;i<Ngrid_init*Ngrid_init;++i){
		assert(i_points[i].leaf->child1 == NULL && i_points[i].leaf->child2 == NULL);
		assert(s_points[i].leaf->child1 == NULL && s_points[i].leaf->child2 == NULL);
	}*/
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
 * be result->head or, if no points are added, result = NULL.
 *
 * Returns a pointer to the list of image points that have been added.  This array can then be
 * used for calculating the surface brightness or marking them as in the image.
 */

Point * Grid::RefineLeaf(LensHndl lens,Point *point,bool kappa_off){
	//Point * RefineLeaf(LensHndl lens,TreeHndl i_tree,TreeHndl s_tree,Point *point,int Ngrid,bool kappa_off){

	Point *i_points = NewPointArray(Ngrid_block*Ngrid_block-1,true);
	Point *s_points;
	int Nout,kk;

	/* Test lines
	if(!testLeafs(i_tree)){ERROR_MESSAGE(); exit(1);}
	if(!testLeafs(s_tree)){ERROR_MESSAGE(); exit(1);}
	*/

	assert(point->leaf->child1 == NULL && point->leaf->child2 == NULL);
	assert(point->image->leaf->child1 == NULL && point->image->leaf->child2 == NULL);

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
	if(Ngrid_block*Ngrid_block-1-Nout <=0){
		FreePointArray(i_points);
		FreePointArray(s_points);
		point->leaf->refined = false;
		point->gridsize *= Ngrid_block;
		point->image->gridsize *= Ngrid_block;

		return NULL;
	}

	if(Nout > 0){
		//i_points = AddPointToArray(i_points,Ngrid_block*Ngrid_block-1-Nout,Ngrid_block*Ngrid_block-1);
		//s_points = AddPointToArray(s_points,Ngrid_block*Ngrid_block-1-Nout,Ntemp);

		i_points = AddPointToArray(i_points,Ngrid_block*Ngrid_block-1-Nout,i_points->head);
		s_points = AddPointToArray(s_points,Ngrid_block*Ngrid_block-1-Nout,s_points->head);
	}
	assert(i_points->head == s_points->head);

	//*** these could be mode more efficient by starting at the current in tree
	AddPointsToTree(i_tree,i_points,i_points->head);
	/* test lines ////////////////////////////////////////////////
	for(int i=0;i<i_points->head;++i){
		assert(i_points[i].leaf->child1 == NULL && i_points[i].leaf->child2 == NULL);
		assert(inbox(i_points[i].x,i_points[i].leaf->boundary_p1,i_points[i].leaf->boundary_p2));
	}
	if(!testLeafs(i_tree)){ERROR_MESSAGE(); std::cout << "point id "<< point->id << std::endl; exit(1);}
	///////////////////////////////////////////////*/
	AddPointsToTree(s_tree,s_points,s_points->head);
	/* test lines ////////////////////////////////////////////////
	for(int i=0;i<s_points->head;++i){
		assert(s_points[i].leaf->child1 == NULL && s_points[i].leaf->child2 == NULL);
		assert(inbox(s_points[i].x,s_points[i].leaf->boundary_p1,s_points[i].leaf->boundary_p2));
	}
	for(int i=0;i<i_points->head;++i)
		assert(i_points[i].image->leaf->child1 == NULL && i_points[i].image->leaf->child2 == NULL);
	if(!testLeafs(s_tree)){ERROR_MESSAGE(); std::cout << "point id "<< point->image->id << std::endl; exit(1);}
	///////////////////////////////////////////////*/

	assert(s_points->head > 0);
	//AddPointsToTree(i_tree,i_points,Ngrid_block*Ngrid_block-1-Nout);
	//AddPointsToTree(s_tree,s_points,Ngrid_block*Ngrid_block-1-Nout);

	// re-assign leaf of point that was to be refined
	assert(inbox(point->x,i_tree->top->boundary_p1,i_tree->top->boundary_p2));
	i_tree->current = point->leaf;
	assert(inbox(point->x,i_tree->current->boundary_p1,i_tree->current->boundary_p2));
	// This line should not be necessary!! It is repairing the leaf that has been assigned incorrectly somewhere
	//if(!inbox(point->x,i_tree->current->boundary_p1,i_tree->current->boundary_p2) ) moveTop(i_tree);
	_FindLeaf(i_tree,point->x,0);
	point->leaf = i_tree->current;

	assert(inbox(point->image->x,s_tree->top->boundary_p1,s_tree->top->boundary_p2));
	assert(inbox(point->image->x,s_tree->top->boundary_p1,s_tree->top->boundary_p2));
	s_tree->current = point->image->leaf;
	assert(inbox(point->image->x,s_tree->current->boundary_p1,s_tree->current->boundary_p2));
	// This line should not be necessary!! It is repairing the leaf that has been assigned incorrectly somewhere
	//if(!inbox(point->image->x,s_tree->current->boundary_p1,s_tree->current->boundary_p2) ) moveTop(s_tree);
	_FindLeaf(s_tree,point->image->x,0);
	point->image->leaf = s_tree->current;

	//Test lines
	assert(point->leaf->child1 == NULL && point->leaf->child2 == NULL);
	assert(point->image->leaf->child1 == NULL && point->image->leaf->child2 == NULL);
	/* Test lines
	if(!testLeafs(i_tree)){ERROR_MESSAGE(); std::cout << "point id "<< point->id; exit(1);}
	if(!testLeafs(s_tree)){ERROR_MESSAGE(); std::cout << "point id "<< point->image->id; exit(1);}
	*/
	return i_points;
}
/**
 * Same as RefineLeaf() but multiple points can be passed.  The rays are shot all together so that more
 * parallelization can be achieved in the rayshooting.
 */
Point * Grid::RefineLeaves(LensHndl lens,std::vector<Point *>& points,bool kappa_off){

	if(points.size() == 0) return NULL;

	long Nleaves = points.size();
	Point *i_points = NewPointArray((Ngrid_block*Ngrid_block-1)*Nleaves,true);
	Point *s_points;
	int Nout,kk,ii,addedtocell[Nleaves];
	long Nadded,Nout_tot;

	Nout_tot=0;
	for(ii=0,Nadded=0;ii<Nleaves;++ii){
		assert(points[ii]->leaf->child1 == NULL && points[ii]->leaf->child2 == NULL);
		assert(points[ii]->image->leaf->child1 == NULL && points[ii]->image->leaf->child2 == NULL);

		assert(points[ii]->gridsize > pow(10.,-DBL_DIG) ); // If cells are too small they will cause problems.

		points[ii]->leaf->refined = true;
		xygridpoints(&i_points[Nadded],points[ii]->gridsize*(Ngrid_block-1)/Ngrid_block
				,points[ii]->x,Ngrid_block,1);
		points[ii]->gridsize /= Ngrid_block;
		points[ii]->image->gridsize /= Ngrid_block;

		// take out points that are outside of original grid
		Nout = 0;
		if( (points[ii]->x[0] == i_tree->top->boundary_p1[0]) || (points[ii]->x[0] == i_tree->top->boundary_p2[0])
			|| (points[ii]->x[1] == i_tree->top->boundary_p1[1]) || (points[ii]->x[1] == i_tree->top->boundary_p2[1]) ){

			// remove the points that are outside initial image grid
			for(kk=0,Nout=0;kk < (Ngrid_block*Ngrid_block-1);++kk){
				if( !inbox(i_points[Nadded + kk - Nout].x,i_tree->top->boundary_p1,i_tree->top->boundary_p2) ){
					//SwapPointsInArray(&i_points[Nadded + kk - Nout],&i_points[(Ngrid_block*Ngrid_block-1)*Nleaves - 1 - Nout]);

					// This maintains the ordering in parent cells, but is rather inefficient
					for(unsigned long nn=Nadded + kk - Nout ; nn < (Ngrid_block*Ngrid_block-1)*Nleaves - 1 - Nout - Nout_tot ; ++nn)
						SwapPointsInArray(&i_points[nn],&i_points[nn + 1]);
					++Nout;
					//std::cout << "Nout_tot = " << Nout_tot << std::endl;
				}
			}
			assert(Nout > 0);
		}

		Nout_tot += Nout;
		Nadded += Ngrid_block*Ngrid_block-1 - Nout;
		addedtocell[ii] = Ngrid_block*Ngrid_block-1 - Nout;
		//if(Nout > 0) i_points = AddPointToArray(i_points,Ngrid_block*Ngrid_block-1-Nout,Ngrid_block*Ngrid_block-1);
	}
	assert(Nadded == (Ngrid_block*Ngrid_block-1)*Nleaves-Nout_tot);

	s_points = LinkToSourcePoints(i_points,Nadded);

	// Here the points that are in uniform magnification regions are calculated by interpolation and marked with in_image = MAYBE
	{
		double aa[4],dx[2];
		/*double a1[4],ss;
		bool tmp;
		Point *i_point = NewPointArray(1,true);
		Point *s_point = LinkToSourcePoints(i_point,1);
		 */

		for(ii=0,kk=0;ii<Nleaves;++ii){

			if(uniform_mag_from_deflect(aa,points[ii])){
				/*
				 tmp = uniform_mag_from_shooter(a1,points[ii]);

				std::cout << "shooter    a = " << a1[0] << " " << a1[1] << " " << a1[2] << " " << a1[3] << std::endl;
				std::cout << "deflection a = " << a2[0] << " " << a2[1] << " " << a2[2] << " " << a2[3] << std::endl;
				std::cout << "          da = " << (a2[0]-a1[0]) << " " << (a2[1]-a1[1]) << " " << (a2[2]-a1[2])
						<< " " << (a2[3]-a1[3]) << std::endl;
				*/

				for(unsigned long jj = kk; jj < addedtocell[ii] + kk; ++jj){
					dx[0] = i_points[jj].x[0] - points[ii]->x[0];
					dx[1] = i_points[jj].x[1] - points[ii]->x[1];

					i_points[jj].image->x[0] = points[ii]->image->x[0] + aa[0]*dx[0] + aa[2]*dx[1];
					i_points[jj].image->x[1] = points[ii]->image->x[1] + aa[3]*dx[0] + aa[1]*dx[1];
					/*
					i_point->x[0] = i_points[jj].x[0];
					i_point->x[1] = i_points[jj].x[1];
					lens->rayshooterInternal(1,i_point,false);

					ss = sqrt( pow(i_point->image->x[0] - points[ii]->image->x[0],2)
							+ pow(i_point->image->x[1] - points[ii]->image->x[1],2) );

					std::cout << (i_point->image->x[0] - i_points[jj].image->x[0])/ss << "  "
							  << (i_point->image->x[1] - i_points[jj].image->x[1])/ss  << std::endl;
					 */

					if(!kappa_off){
						i_points[jj].kappa = i_points[jj].image->kappa = 1-0.5*(aa[0]+aa[1]);
						i_points[jj].gamma[0] = i_points[jj].image->gamma[0] = -0.5*(aa[0]-aa[1]);
						i_points[jj].gamma[1] = i_points[jj].image->gamma[1] = -0.5*(aa[2]+aa[3]);
						i_points[jj].gamma[2] = i_points[jj].image->gamma[2] = -0.5*(aa[2]-aa[3]);
						i_points[jj].invmag = i_points[jj].image->invmag = (1-i_points[jj].kappa)*(1-i_points[jj].kappa)
								     - i_points[jj].gamma[0]*i_points[jj].gamma[0]
								     - i_points[jj].gamma[1]*i_points[jj].gamma[1]
								     + i_points[jj].gamma[2]*i_points[jj].gamma[2];
					}
					i_points[jj].in_image = MAYBE;

				}

				//std::cout << " " << std::endl;
			}
			kk += addedtocell[ii];
		}

		//FreePointArray(i_point);
		//FreePointArray(s_point);

	}

	lens->rayshooterInternal(Nadded,i_points,kappa_off);

	// remove the points that are outside initial source grid
	int j,Noutcell;
	for(ii=0,kk=0,Nout=0;ii<Nleaves;++ii){
		for(j = 0,Noutcell=0; j < addedtocell[ii]; ++j){
			if( !inbox(s_points[kk - Nout].x,s_tree->top->boundary_p1,s_tree->top->boundary_p2) ){
				//SwapPointsInArray(&i_points[kk - Nout],&i_points[Nadded - 1 - Nout]);
				//SwapPointsInArray(&s_points[kk - Nout],&s_points[Nadded - 1 - Nout]);
				for(unsigned long nn = kk - Nout; nn < Nadded - 2 - Nout;++nn){
					SwapPointsInArray(&i_points[nn],&i_points[nn+1]);
					SwapPointsInArray(&s_points[nn],&s_points[nn+1]);
				}
				++Nout;
				++Noutcell;
			}
			++kk;
		}
		addedtocell[ii] -= Noutcell;
		if(addedtocell[ii] == 0){  // case where all of the parent cell is out of source plane region
			points[ii]->leaf->refined = false;
			points[ii]->gridsize *= Ngrid_block;
			points[ii]->image->gridsize *= Ngrid_block;
		}
	}

	assert(i_points->head == (Ngrid_block*Ngrid_block-1)*Nleaves);
	assert(s_points->head == Nadded);

	Nadded -= Nout;
/*
	for(kk=0,Nout=0;kk < Nadded;++kk){
		assert(s_points[kk - Nout].x[0] == s_points[kk - Nout].x[0]);
		if( !inbox(s_points[kk - Nout].x,s_tree->top->boundary_p1,s_tree->top->boundary_p2) ){
			SwapPointsInArray(&i_points[kk - Nout],&i_points[Nadded - 1 - Nout]);
			SwapPointsInArray(&s_points[kk - Nout],&s_points[Nadded - 1 - Nout]);
			++Nout;
		}
	}
*/

	// free memory of points that where outside image and source regions
	if(Nadded == 0){
		FreePointArray(i_points);
		FreePointArray(s_points);
		return NULL;
	}

	if(Nadded < (Ngrid_block*Ngrid_block-1)*Nleaves){
		i_points = AddPointToArray(i_points,Nadded,i_points->head);
		s_points = AddPointToArray(s_points,Nadded,s_points->head);
	}
	assert(i_points->head == s_points->head);

	//*** these could be mode more efficient by starting at the current in tree
	AddPointsToTree(i_tree,i_points,i_points->head);
	AddPointsToTree(s_tree,s_points,s_points->head);

	assert(s_points->head > 0);

	//********* This repairs the leaf pointer which for some reason does not point to a leaf for s_points on some occasions ***********
	for(ii=0;ii < s_points->head;++ii){

		if(i_points[ii].leaf->child1 != NULL || i_points[ii].leaf->child2 != NULL){
			//std::cout << ii << std::endl;
			//i_points[ii].print();
			//i_points[ii].leaf->print();
			i_tree->current = i_points[ii].leaf;
			_FindLeaf(i_tree,i_points[ii].x,0);
			assert(i_tree->current->npoints == 1);
			i_points[ii].leaf = i_tree->current;
		}

		if(s_points[ii].leaf->child1 != NULL || s_points[ii].leaf->child2 != NULL){
			//std::cout << ii << std::endl;
			//s_points[ii].print();
			//s_points[ii].leaf->print();
			s_tree->current = s_points[ii].leaf;
			_FindLeaf(s_tree,s_points[ii].x,0);
			assert(s_tree->current->npoints == 1);
			s_points[ii].leaf = s_tree->current;
		}
	}
	//*********************************************************************

	// This loop should not be necessary!! It is repairing the leaf that has been assigned incorrectly somewhere
	//if(!inbox(point[ii].x,i_tree->current->boundary_p1,i_tree->current->boundary_p2) ) moveTop(i_tree);
	/*for(ii=0;ii < Nleaves;++ii){

		// re-assign leaf of point that was to be refined
		assert(inbox(points[ii]->x,i_tree->top->boundary_p1,i_tree->top->boundary_p2));
		i_tree->current = points[ii]->leaf;
		assert(inbox(points[ii]->x,i_tree->current->boundary_p1,i_tree->current->boundary_p2));
		_FindLeaf(i_tree,points[ii]->x,0);
		points[ii]->leaf = i_tree->current;

		assert(inbox(points[ii]->image->x,s_tree->top->boundary_p1,s_tree->top->boundary_p2));
		assert(inbox(points[ii]->image->x,s_tree->top->boundary_p1,s_tree->top->boundary_p2));
		s_tree->current = points[ii]->image->leaf;
		assert(inbox(points[ii]->image->x,s_tree->current->boundary_p1,s_tree->current->boundary_p2));
		_FindLeaf(s_tree,points[ii]->image->x,0);
		points[ii]->image->leaf = s_tree->current;

		//Test lines
		assert(points[ii]->leaf->child1 == NULL && points[ii]->leaf->child2 == NULL);
		assert(points[ii]->image->leaf->child1 == NULL && points[ii]->image->leaf->child2 == NULL);
	}*/

	return i_points;
}

/// Rest all in_image markers to False.
void Grid::ClearAllMarks(){
	unsigned long i;

	MoveToTopList(i_tree->pointlist);
	for(i=0;i<i_tree->pointlist->Npoints;++i){
		i_tree->pointlist->current->in_image=FALSE;
		i_tree->pointlist->current->image->in_image=FALSE;
		MoveDownList(i_tree->pointlist);
	}
}
/**
 *\brief Find the magnification given three points.
 *
 * The points must have attached source/image points and
 * they must not be colinear. Returns false if they are colinear
 *
 * 	 a[0] = a11, a[1] = a22, a[2] = a12, a[3] = a21
 *
 */
bool Grid::find_mag_matrix(double *a,Point *p0,Point *p1,Point *p2){
	double y1[2],y2[2],x1[2],x2[2];

	x1[0] = p1->x[0] - p0->x[0];
	x1[1] = p1->x[1] - p0->x[1];
	x2[0] = p2->x[0] - p0->x[0];
	x2[1] = p2->x[1] - p0->x[1];

	y1[0] = p1->image->x[0] - p0->image->x[0];
	y1[1] = p1->image->x[1] - p0->image->x[1];
	y2[0] = p2->image->x[0] - p0->image->x[0];
	y2[1] = p2->image->x[1] - p0->image->x[1];

	double det = x1[0]*x2[1] - x1[1]*x2[0];
	if(det == 0) return false;

	// a[0] = a11, a[1] = a22, a[2] = a12, a[3] = a21
	a[0] = (  x2[1]*y1[0] - x1[1]*y2[0] )/det;
	a[3] = ( -x2[0]*y1[0] + x1[0]*y2[0] )/det;
	a[2] = (  x2[1]*y1[1] - x1[1]*y2[1] )/det;
	a[1] = ( -x2[0]*y1[1] + x1[0]*y2[1] )/det;

	return true;
}
/** \brief Test if point is in a region of uniform magnification using only the deflections of the point and
 * its neighbors.
 *
 * An estimate of the magnification matrix is returned if it returns true.  Otherwise the magnification
 * matrix is unspecified.
 *
 * Magnification matrix elements are considered equal if their difference is smaller than maglimit which is
 * set in the Grid constructor.
 */
bool Grid::uniform_mag_from_deflect(
		double *a           /// Returned magnification matrix, not specified if returning false
		,Point *point       /// point to be tested
		){
	Point *point2;
	double ao[4];
	int count=0;

    FindAllBoxNeighborsKist(i_tree,point,neighbors);
    if(neighbors->Nunits() <= 3) return false;  // This is the case where the point does not have enough neighbors.
    neighbors->MoveToTop();
    point2 = neighbors->getCurrent();
    neighbors->Down();
	while(!find_mag_matrix(a,point,point2,neighbors->getCurrent())) neighbors->Down();
	//std::cout << "deflection neighbors a" << std::endl;
	while(neighbors->Down()){
		if(find_mag_matrix(ao,point,point2,neighbors->getCurrent())){
			if( !( (fabs(a[0]-ao[0]) < maglimit)*(fabs(a[1]-ao[1]) < maglimit)
    				*(fabs(a[2]-ao[2]) < maglimit)*(fabs(a[3]-ao[3]) < maglimit) )) return false;
			//std::cout  << ao[0] << " " << ao[1] << " " << ao[2] << " " << ao[3] << std::endl;

			a[0] = (count*a[0] + ao[0])/(count+1);
			a[1] = (count*a[1] + ao[1])/(count+1);
			a[2] = (count*a[2] + ao[2])/(count+1);
			a[3] = (count*a[3] + ao[3])/(count+1);

    		++count;
    	}
    }

	if(count == 0) return false;  // This is the case where the point does not have enough non-colinear neighbors.
    return true;
}
/** \brief Test if point is in a region of uniform magnification using the kappa and gamma calculated from the
 * rayshooter.
 *
 * An estimate of the magnification matrix is returned if it returns true.  Otherwise the magnification
 * matrix is unspecified.
 *
 * Magnification matrix elements are considered equal if their difference is smaller than maglimit which is
 * set in the Grid constructor.
 *
 */
bool Grid::uniform_mag_from_shooter(
		double *a           /// Returned magnification matrix, not specified if returning false
		,Point *point       /// point to be tested
		){
	Point *point2;

    FindAllBoxNeighborsKist(i_tree,point,neighbors);
    assert(neighbors->Nunits() > 1);
    neighbors->MoveToTop();
 	do{
 		point2 = neighbors->getCurrent();
 		if( !( (fabs(point->kappa-point2->kappa) < maglimit)*(fabs(point->gamma[0]-point2->gamma[0]) < maglimit)
 				*(fabs(point->gamma[1]-point2->gamma[1]) < maglimit)*(fabs(point->gamma[2]-point2->gamma[2]) < maglimit) ) ) return false;
    }while(neighbors->Down());

    neighbors->MoveToTop();
    std::cout << "shooter neighbors " << std::endl;
  	do{
  		point2 = neighbors->getCurrent();
 		std::cout << point->kappa-point2->kappa << "  " << point->gamma[0]-point2->gamma[0] << "   " <<
  				point->gamma[1]-point2->gamma[1] << "  " << point->gamma[2]-point2->gamma[2] << std::endl;
     }while(neighbors->Down());


	a[0] = 1 - point->kappa - point->gamma[0];
	a[1] = 1 - point->kappa + point->gamma[0];
	a[2] = -point->gamma[1] - point->gamma[2];
	a[3] = -point->gamma[1] + point->gamma[2];

    return true;
}
void Grid::test_mag_matrix(){
	double aa[4];
	Point *point2;
	double ao[4];
	int count;
	double kappa, gamma[3], invmag;

	MoveToTopList(i_tree->pointlist);
	do{
		count=0;
		FindAllBoxNeighborsKist(i_tree,i_tree->pointlist->current,neighbors);
		assert(neighbors->Nunits() >= 3);
		neighbors->MoveToTop();
		point2 = neighbors->getCurrent();
		neighbors->Down();
		while(!find_mag_matrix(aa,i_tree->pointlist->current,point2,neighbors->getCurrent())) neighbors->Down();
		//std::cout << "deflection neighbors a" << std::endl;
		while(neighbors->Down()){
			if(find_mag_matrix(ao,i_tree->pointlist->current,point2,neighbors->getCurrent())){
				aa[0] = (count*aa[0] + ao[0])/(count+1);
				aa[1] = (count*aa[1] + ao[1])/(count+1);
				aa[2] = (count*aa[2] + ao[2])/(count+1);
				aa[3] = (count*aa[3] + ao[3])/(count+1);

				++count;
			}
		}

		kappa = 1-0.5*(aa[0]+aa[1]);
		gamma[0] = -0.5*(aa[0]-aa[1]);
		gamma[1] = -0.5*(aa[2]+aa[3]);
		gamma[2] = -0.5*(aa[2]-aa[3]);
		invmag = (1-kappa)*(1-kappa)-gamma[0]*gamma[0]-gamma[1]*gamma[1]+gamma[2]*gamma[2];

		i_tree->pointlist->current->kappa =  kappa/i_tree->pointlist->current->kappa - 1.0;
		i_tree->pointlist->current->gamma[0] = gamma[0]/i_tree->pointlist->current->gamma[0] - 1.0;
		i_tree->pointlist->current->gamma[1] = gamma[1]/i_tree->pointlist->current->gamma[1] - 1.0;
		i_tree->pointlist->current->gamma[2] = gamma[2]/i_tree->pointlist->current->gamma[2] - 1.0;
		i_tree->pointlist->current->invmag = invmag/i_tree->pointlist->current->invmag - 1.0;

	}while(MoveDownList(i_tree->pointlist));
}
