/*
 * grid_maintenance.c
 *
 *  Created on: Oct 12, 2011
 *      Author: bmetcalf
 */

#include "slsimlib.h"
#include "Tree.h"
#include <mutex>
#include <thread>
#include "grid_maintenance.h"
#include "source.h"

std::mutex Grid::grid_mutex;

/** 
 * \brief Constructor for initializing square grid.
 *
 * Note: Deflection solver must be specified before creating a Grid.
 */
Grid::Grid(
		LensHndl lens               /// lens model for initializing grid
		,unsigned long N1d          /// Initial number of grid points in each dimension.
		,const PosType center[2]    /// Center of grid (usually in radian units)
		,PosType range              /// Full width of grid in whatever units will be used.
           ): Ngrid_init(N1d),Ngrid_init2(N1d),Ngrid_block(3),axisratio(1.0){

	Point *i_points,*s_points;
    pointID = 0;

	assert(N1d > 0);
	assert(range > 0);

	if(N1d <= 0){ERROR_MESSAGE(); std::cout << "cannot make Grid with no points" << std::endl; exit(1);}
	if(range <= 0){ERROR_MESSAGE(); std::cout << "cannot make Grid with no range" << std::endl; exit(1);}
  //if( (N1d & (N1d-1)) != 0 ){ERROR_MESSAGE(); std::printf("ERROR: Grid cannot be initialized with pixels less than a power of 2\n"); exit(1);}

	i_points = NewPointArray(Ngrid_init*Ngrid_init);
	xygridpoints(i_points,range,center,Ngrid_init,0);
	s_points=LinkToSourcePoints(i_points,Ngrid_init*Ngrid_init);

  {
    std::lock_guard<std::mutex> hold(grid_mutex);
  	lens->rayshooterInternal(Ngrid_init*Ngrid_init,i_points);
  }
  
	// Build trees
	i_tree = new TreeStruct(i_points,Ngrid_init*Ngrid_init);
	s_tree = new TreeStruct(s_points,Ngrid_init*Ngrid_init,1,range);  // make tree on source plane with a buffer
             
	trashkist = new Kist<Point>;
	neighbors = new Kist<Point>;
	maglimit = 1.0e-4;
  
}

/** 
 * \brief Constructor for initializing rectangular grid.
 *
 * Cells of grid will always be square with initial resolution rangeX/(Nx-1).
 * The Y range may not be exactly rangeY, but will be the nearest value that 
 * is a whole number of cells.
 *
 * Note: Deflection solver must be specified before creating a Grid.
 */
Grid::Grid(
           LensHndl lens       /// lens model for initializing grid
           ,unsigned long Nx   /// Initial number of grid points in X dimension.
           ,const PosType center[2]  /// Center of grid.
           ,PosType rangeX   /// Full width of grid in x direction in whatever units will be used.
           ,PosType rangeY  /// Full width of grid in y direction in whatever units will be used.
           ): Ngrid_init(Nx),Ngrid_block(3),axisratio(rangeY/rangeX){
    
	Point *i_points,*s_points;
    pointID = 0;
    
	assert(Nx > 0);
	assert(rangeX > 0 && rangeY >0);
    
	if(Nx <= 0){ERROR_MESSAGE();
        std::cout << "cannot make Grid with no points" << std::endl;
        throw std::runtime_error("");
    }
	if(rangeX <= 0 || rangeY <= 0 ){ERROR_MESSAGE();
        throw std::runtime_error("");
        std::cout << "cannot make Grid with no range" << std::endl;
    }

    //if(Ngrid_init % 2 == 1 ) ++Ngrid_init;
	
    Ngrid_init2 = (int)(Ngrid_init*axisratio);
    //if(Ngrid_init2 % 2 == 1) ++Ngrid_init2;
    
    i_points = NewPointArray(Ngrid_init*Ngrid_init2);
    
    int i;
    // set grid of positions
    for(int ii=0;ii<Ngrid_init;++ii){
        for(int jj=0;jj<Ngrid_init2;++jj){

            i = ii + jj*Ngrid_init;
            
            i_points[i].id=pointID;
            ++pointID;
            
            i_points[i].x[0] = center[0] + rangeX*ii/(Ngrid_init-1) - 0.5*rangeX;
            i_points[i].x[1] = center[1] + rangeX*jj/(Ngrid_init-1) - 0.5*rangeY;
            i_points[i].gridsize=rangeX/(Ngrid_init-1);
        }
    }

    assert(i == Ngrid_init*Ngrid_init2-1);
 
    s_points=LinkToSourcePoints(i_points,Ngrid_init*Ngrid_init2);
    
    {
        std::lock_guard<std::mutex> hold(grid_mutex);
        lens->rayshooterInternal(Ngrid_init*Ngrid_init2,i_points);
    }
              
	// Build trees
	i_tree = new TreeStruct(i_points,Ngrid_init*Ngrid_init2);
	s_tree = new TreeStruct(s_points,Ngrid_init*Ngrid_init2,1,MAX(rangeX,rangeY));  // make tree on source plane with a buffer
    
	trashkist = new Kist<Point>;
	neighbors = new Kist<Point>;
	maglimit = 1.0e-4;
    
}

/** 
 * \brief Destructor for a Grid.  Frees all memory.
 */
Grid::~Grid(){
	delete i_tree;
	delete s_tree;
	
	delete trashkist;
	delete neighbors;
  
  return;
}

/** Finding
 *  \brief Reinitializes the grid so that it is back to the original coarse grid, but if
 *  the lens has changed the source positions will be updated.
 */
void Grid::ReInitializeGrid(LensHndl lens){
  
	Point *i_points,*s_points;
	PosType rangeX,rangeY,center[2];
	unsigned long i;
  
	rangeX = (i_tree->getTop()->boundary_p2[0] - i_tree->getTop()->boundary_p1[0]);
  rangeY = (i_tree->getTop()->boundary_p2[1] - i_tree->getTop()->boundary_p1[1]);
	center[0] = (i_tree->getTop()->boundary_p2[0] + i_tree->getTop()->boundary_p1[0])/2;
	center[1] = (i_tree->getTop()->boundary_p2[1] + i_tree->getTop()->boundary_p1[1])/2;
  
	//////////////////////////////
  // redo grid with stars in it
  // free old tree to speed up image finding
  
	i_tree->emptyTree();
	s_tree->emptyTree();
  

	// build new initial grid
	i_points = NewPointArray(Ngrid_init*Ngrid_init2);
	if(Ngrid_init == Ngrid_init2){
    xygridpoints(i_points,rangeX,center,Ngrid_init,0);
  }else{
    int i;
    // set grid of positions
    for(int ii=0;ii<Ngrid_init;++ii){
      for(int jj=0;jj<Ngrid_init2;++jj){
        
        i = ii + jj*Ngrid_init;
        
        i_points[i].id=pointID;
        ++pointID;
        
        i_points[i].x[0] = center[0] + rangeX*ii/(Ngrid_init-1) - 0.5*rangeX;
        i_points[i].x[1] = center[1] + rangeX*jj/(Ngrid_init-1) - 0.5*rangeY;
        i_points[i].gridsize=rangeX/(Ngrid_init-1);
      }
    }

  }
	s_points=LinkToSourcePoints(i_points,Ngrid_init*Ngrid_init2);
  
  {
    std::lock_guard<std::mutex> hold(grid_mutex);
	  lens->rayshooterInternal(Ngrid_init*Ngrid_init2,i_points);
  }
	// need to resize root of source tree.  It can change in size
	s_tree->getTop()->boundary_p1[0]=s_points[0].x[0]; s_tree->getTop()->boundary_p1[1]=s_points[0].x[1];
	s_tree->getTop()->boundary_p2[0]=s_points[0].x[0]; s_tree->getTop()->boundary_p2[1]=s_points[0].x[1];
  
	for(i=0;i<Ngrid_init*Ngrid_init2;++i){
    
    /* find X boundary */
		if(s_points[i].x[0] < s_tree->getTop()->boundary_p1[0] ) s_tree->getTop()->boundary_p1[0]=s_points[i].x[0];
    if(s_points[i].x[0] > s_tree->getTop()->boundary_p2[0] ) s_tree->getTop()->boundary_p2[0]=s_points[i].x[0];
    
    /* find Y boundary */
    if(s_points[i].x[1] < s_tree->getTop()->boundary_p1[1] ) s_tree->getTop()->boundary_p1[1]=s_points[i].x[1];
    if(s_points[i].x[1] > s_tree->getTop()->boundary_p2[1] ) s_tree->getTop()->boundary_p2[1]=s_points[i].x[1];
  }
  
  // a little extra room for future points
  s_tree->getTop()->boundary_p1[0] -=  MAX(rangeX,rangeY)/Ngrid_init;
  s_tree->getTop()->boundary_p1[1] -=  MAX(rangeX,rangeY)/Ngrid_init;
  s_tree->getTop()->boundary_p2[0] +=  MAX(rangeX,rangeY)/Ngrid_init;
  s_tree->getTop()->boundary_p2[1] +=  MAX(rangeX,rangeY)/Ngrid_init;
  
  s_tree->getTop()->center[0] = (s_tree->getTop()->boundary_p1[0]+s_tree->getTop()->boundary_p2[0])/2;
  s_tree->getTop()->center[1] = (s_tree->getTop()->boundary_p1[1]+s_tree->getTop()->boundary_p2[1])/2;
  
	// fill trees
	i_tree->FillTree(i_points,Ngrid_init*Ngrid_init2);
	s_tree->FillTree(s_points,Ngrid_init*Ngrid_init2);
  
	/*for(i=0;i<Ngrid_init*Ngrid_init;++i){
   assert(i_points[i].leaf->child1 == NULL && i_points[i].leaf->child2 == NULL);
   assert(s_points[i].leaf->child1 == NULL && s_points[i].leaf->child2 == NULL);
   }*/
	return;
}


/** Finding
 *  \brief Reshoot the rays with the same image postions.
 *
 *  The source positions and source tree are updated to the current lens model.
 *  The advantage over Grid::ReInitializeGrid() is that the image plane refinements
 *  are preserved.
 */
void Grid::ReShoot(LensHndl lens){
  
	Point *i_points,*s_points;
	PosType range,center[2];
	unsigned long i;
  
	range = i_tree->getTop()->boundary_p2[0] - i_tree->getTop()->boundary_p1[0];
	center[0] = (i_tree->getTop()->boundary_p2[0] + i_tree->getTop()->boundary_p1[0])/2;
	center[1] = (i_tree->getTop()->boundary_p2[1] + i_tree->getTop()->boundary_p1[1])/2;
  
  // clear source tree
  delete s_tree;
  s_points = NewPointArray(i_tree->pointlist->size());
  
	// build new initial grid
  PointList::iterator i_tree_pointlist_it;
  i_tree_pointlist_it.current = i_tree->pointlist->Top();
  size_t k;
  for(i=0,k=0;i<i_tree->pointlist->size();++i){
    i_points = *i_tree_pointlist_it;
    if(i_points->head > 0){

      // link source and image points
      for(size_t j=0;j<i_points->head;++j,++k){
        i_points[j].image = &s_points[k];
        s_points[k].image = &i_points[j];
        s_points[k].id = i_points[j].id;
        s_points[k].gridsize = i_points[j].gridsize;
      };

      {
      // reshoot the rays
        std::lock_guard<std::mutex> hold(grid_mutex);
         lens->rayshooterInternal(i_points->head,i_points);
      }
    }
    
    --i_tree_pointlist_it;
  }
  
  s_tree = new TreeStruct(s_points,s_points->head,1,(i_tree->getTop()->boundary_p2[0] - i_tree->getTop()->boundary_p1[0])/10 );
	return;
}


/** Finding
 * \brief Recalculate surface brightness at every point without changing the positions of the grid or any lens properties.
 *
 *  Recalculate the surface brightness at all points on the grid.
 * This is useful when changing the source model while preserving
 * changes in the grid.
 * Both i_tree and s_tree are both changed although only s_tree shows up here.
 *
 * returns total flux
 */
PosType Grid::RefreshSurfaceBrightnesses(SourceHndl source){
  PosType total=0,tmp;
  
  SourcePoint *sp = dynamic_cast<SourcePoint *>(source);
  
  if(sp){
    ClearSurfaceBrightnesses();
    total = mark_point_source_images(sp->getTheta(),sp->getRadius(),sp->getTotalFlux());
  }else{
    PointList::iterator s_tree_pointlist_it;
    s_tree_pointlist_it.current = (s_tree->pointlist->Top());
    for(unsigned long i=0;i<s_tree->pointlist->size();++i,--s_tree_pointlist_it){
      tmp = source->SurfaceBrightness((*s_tree_pointlist_it)->x);
      (*s_tree_pointlist_it)->surface_brightness = (*s_tree_pointlist_it)->image->surface_brightness
      = tmp;
      total += tmp*pow( (*s_tree_pointlist_it)->image->gridsize ,2);
      assert((*s_tree_pointlist_it)->surface_brightness >= 0.0);
      (*s_tree_pointlist_it)->in_image = (*s_tree_pointlist_it)->image->in_image
      = NO;
    }
  }
  return total;
}
/**
 * \brief Recalculate surface brightness just like Grid::RefreshSurfaceBrightness but
 * the new source is added to any sources that were already there.  
 *
 * returns the sum of the surface brightnesses from the new source
 */
PosType Grid::AddSurfaceBrightnesses(SourceHndl source){
  PosType total=0,tmp;
  
  SourcePoint *sp = dynamic_cast<SourcePoint *>(source);
 
  if(sp){
    total = mark_point_source_images(sp->getTheta(),sp->getRadius(),sp->getTotalFlux());
  }else{
    PointList::iterator s_tree_pointlist_it;
    s_tree_pointlist_it.current = (s_tree->pointlist->Top());
    for(unsigned long i=0;i<s_tree->pointlist->size();++i,--s_tree_pointlist_it){
      tmp = source->SurfaceBrightness((*s_tree_pointlist_it)->x);
      (*s_tree_pointlist_it)->surface_brightness += tmp;
      (*s_tree_pointlist_it)->image->surface_brightness += tmp;
      total += tmp;//*pow( s_tree->pointlist->current->gridsize,2);
      assert((*s_tree_pointlist_it)->surface_brightness >= 0.0);
      (*s_tree_pointlist_it)->in_image = (*s_tree_pointlist_it)->image->in_image
      = NO;
    }
  }
  
  return total;
}

PosType Grid::EinsteinArea() const {
  PosType total=0;
  PointList::iterator it;
  it = (i_tree->pointlist->Top());
  size_t N = i_tree->pointlist->size();
  for(unsigned long i=0 ; i < N ; ++i,--it){
    if( (*it)->invmag() < 0) total += (*it)->gridsize * (*it)->gridsize;
  }
  
  return total;
}

//PosType Grid::magnification2() const{
//  double mag = 0,flux = 0;
//
//  PointList::iterator it;
//  it = (i_tree->pointlist->Top());
//  size_t N = i_tree->pointlist->size();
//  for(unsigned long i=0 ; i < N ; ++i,--it){
//    double f = (*it)->surface_brightness * (*it)->gridsize * (*it)->gridsize;
//    assert(f >= 0);
//    if(f > 0){
//      mag += f;
//      flux += f/fabs((*it)->invmag());
//    }
//  }
//
//  return flux/mag;
//}

PosType Grid::magnification() const {
  
  return LensedFlux()/UnlensedFlux();
  
  double mag = 0,flux = 0;
  
  PointList::iterator it;
  it = (i_tree->pointlist->Top());
  size_t N = i_tree->pointlist->size();
  for(unsigned long i=0 ; i < N ; ++i,--it){
    double f = (*it)->surface_brightness * (*it)->gridsize * (*it)->gridsize;
    assert(f >= 0);
    if(f > 0){
      mag += f;
      flux += f / mag_from_deflect((*it));
    }
  }

  return mag / flux;
}


//PosType Grid::magnification() const{
//  double mag = 0,flux = 0;
//
//  PointList::iterator it;
//  it = (i_tree->pointlist->Top());
//  size_t N = i_tree->pointlist->size();
//  for(unsigned long i=0 ; i < N ; ++i,--it){
//    double f = (*it)->surface_brightness * (*it)->gridsize * (*it)->gridsize;
//    assert(f >= 0);
//    if(f > 0){
//      mag += f*fabs((*it)->invmag());
//      flux += f;
//    }
//  }
//
//  return flux/mag;
//}

Point_2d Grid::centroid() const{
  double flux = 0;
  Point_2d centroid(0,0);
  
  PointList::iterator it;
  it = (i_tree->pointlist->Top());
  size_t N = i_tree->pointlist->size();
  for(unsigned long i=0 ; i < N ; ++i,--it){
    double area = (*it)->gridsize * (*it)->gridsize;
    centroid += *(*it) * (*it)->surface_brightness*area;
    flux += (*it)->surface_brightness*area;
  }

  return centroid/flux;
}

PosType Grid::UnlensedFlux() const{
  
  double flux = 0;
  PointList::iterator s_tree_pointlist_it;
  s_tree_pointlist_it.current = (s_tree->pointlist->Top());
  for(unsigned long i=0;i<s_tree->pointlist->size();++i,--s_tree_pointlist_it){
    flux += (*s_tree_pointlist_it)->surface_brightness*(*s_tree_pointlist_it)->leaf->area();
  }
  
  return flux;
}

PosType Grid::LensedFlux() const{
  
  double flux = 0;
  PointList::iterator it;
  it = (i_tree->pointlist->Top());
  size_t N = i_tree->pointlist->size();
  for(unsigned long i=0 ; i < N ; ++i,--it){
    flux += (*it)->surface_brightness * (*it)->gridsize * (*it)->gridsize;
  }
  
  return flux;
}

/**
 *  \brief Reset the surface brightness and in_image flag in every point on image and source planes to zero (false)
 */
PosType Grid::ClearSurfaceBrightnesses(){
	PosType total=0;
  
  PointList::iterator s_tree_pointlist_it;
  s_tree_pointlist_it.current = (s_tree->pointlist->Top());
	for(unsigned long i=0;i<s_tree->pointlist->size();++i,--s_tree_pointlist_it){
		(*s_tree_pointlist_it)->surface_brightness = (*s_tree_pointlist_it)->image->surface_brightness
    = 0.0;
		(*s_tree_pointlist_it)->in_image = (*s_tree_pointlist_it)->image->in_image
    = NO;
	}
  
	return total;
}

/** Finding
 * \brief Returns number of points on image plane.
 */
unsigned long Grid::getNumberOfPoints() const{
	assert(i_tree->getTop()->npoints == s_tree->getTop()->npoints);
	assert(i_tree->getTop()->npoints == i_tree->pointlist->size());
	assert(s_tree->getTop()->npoints == s_tree->pointlist->size());

	return i_tree->getTop()->npoints;
}

/**  
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
 *
 * i_tree current is left in one of the new subcells.
 */

Point * Grid::RefineLeaf(LensHndl lens,Point *point){

	Point *i_points = NewPointArray(Ngrid_block*Ngrid_block-1);
	Point *s_points;
	int Nout,kk;

	/*************** TODO Test lines *********************************************
	  ERROR_MESSAGE();
	if(!testLeafs(i_tree)){ERROR_MESSAGE(); std::cout << "point id "<< point->id; exit(1);}
	  ERROR_MESSAGE();
	if(!testLeafs(s_tree)){ERROR_MESSAGE(); std::cout << "point id "<< point->image->id; exit(1);}
	/ *****************************************************************************/

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
	if( (point->x[0] == i_tree->getTop()->boundary_p1[0]) || (point->x[0] == i_tree->getTop()->boundary_p2[0])
			|| (point->x[1] == i_tree->getTop()->boundary_p1[1]) || (point->x[1] == i_tree->getTop()->boundary_p2[1]) ){

		  // remove the points that are outside initial image grid
		  for(kk=0,Nout=0;kk < (Ngrid_block*Ngrid_block-1);++kk){
			  if( !inbox(i_points[kk - Nout].x,i_tree->getTop()->boundary_p1,i_tree->getTop()->boundary_p2) ){
				  SwapPointsInArray(&i_points[kk - Nout],&i_points[Ngrid_block*Ngrid_block - 2 - Nout]);
				  ++Nout;
			  }
		  }
		  assert(Nout > 0);
	}

	//if(Nout > 0) i_points = AddPointToArray(i_points,Ngrid_block*Ngrid_block-1-Nout,Ngrid_block*Ngrid_block-1);

	int  Ntemp = Ngrid_block*Ngrid_block-1-Nout;

	s_points = LinkToSourcePoints(i_points,Ntemp);
  
  {
    std::lock_guard<std::mutex> hold(grid_mutex);
	  lens->rayshooterInternal(Ntemp,i_points);
  }
  
	// remove the points that are outside initial source grid
	for(kk=0,Nout=0;kk < Ntemp;++kk){
		assert(s_points[kk - Nout].x[0] == s_points[kk - Nout].x[0]);
		if( !inbox(s_points[kk - Nout].x,s_tree->getTop()->boundary_p1,s_tree->getTop()->boundary_p2) ){
			assert(Ntemp - 1 - Nout < i_points[0].head);
			assert(kk - Nout < i_points[0].head);
			SwapPointsInArray(&i_points[kk - Nout],&i_points[Ntemp - 1 - Nout]);
			assert(Ntemp - 1 - Nout < s_points[0].head);
			assert(kk - Nout < s_points[0].head);
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
	i_tree->AddPointsToTree(i_points,i_points->head);
	s_tree->AddPointsToTree(s_points,s_points->head);

	assert(s_points->head > 0);

	// re-assign leaf of point that was to be refined
	assert(inbox(point->x,i_tree->getTop()->boundary_p1,i_tree->getTop()->boundary_p2));
  
  TreeStruct::iterator i_tree_it(i_tree);
  i_tree_it = point->leaf;
	assert(inbox(point->x,(*i_tree_it)->boundary_p1,(*i_tree_it)->boundary_p2));
	// This line should not be necessary!! It is repairing the leaf that has been assigned incorrectly somewhere
	//if(!inbox(point->x,(*i_tree_it)->boundary_p1,(*i_tree_it)->boundary_p2) ) moveTop(i_tree);
	i_tree->_FindLeaf(i_tree_it,point->x,0);
	point->leaf = *i_tree_it;

	assert(inbox(point->image->x,s_tree->getTop()->boundary_p1,s_tree->getTop()->boundary_p2));
	assert(inbox(point->image->x,s_tree->getTop()->boundary_p1,s_tree->getTop()->boundary_p2));
  
  TreeStruct::iterator s_tree_it(s_tree);
	s_tree_it = point->image->leaf;
	assert(inbox(point->image->x,(*s_tree_it)->boundary_p1,(*s_tree_it)->boundary_p2));
	// This line should not be necessary!! It is repairing the leaf that has been assigned incorrectly somewhere
	//if(!inbox(point->image->x,s_tree->current->boundary_p1,s_tree->current->boundary_p2) ) moveTop(s_tree);
	s_tree->_FindLeaf(s_tree_it,point->image->x,0);
	point->image->leaf = *s_tree_it;

	return i_points;
}
/**
 * \brief Same as RefineLeaf() but multiple points can be passed.  The rays are shot all together so that more
 * parallelization can be achieved in the rayshooting.
 */
Point * Grid::RefineLeaves(LensHndl lens,std::vector<Point *>& points){

	if(points.size() == 0) return NULL;

	size_t Nleaves = points.size();
	Point *i_points = NewPointArray((Ngrid_block*Ngrid_block-1)*Nleaves);
	Point *s_points;
	size_t Nout,kk,ii;
	size_t Nadded,Nout_tot;
  std::vector<size_t> addedtocell(points.size());

	Nout_tot=0;
	for(ii=0,Nadded=0;ii<Nleaves;++ii){
		assert(points[ii]->leaf->child1 == NULL && points[ii]->leaf->child2 == NULL);
		assert(points[ii]->image->leaf->child1 == NULL && points[ii]->image->leaf->child2 == NULL);
    
    // If cells are too small they will cause problems
		if(points[ii]->gridsize < pow(10.,-DBL_DIG)*MAX(fabs(points[ii]->x[0]),fabs(points[ii]->x[1])) ){
      ERROR_MESSAGE();
      std::cout << "Cell size is too small for double precision " << std::endl;
      throw std::runtime_error("Cell size is too small for double precision");
    };
    //assert(points[ii]->gridsize > 1.0e-10*MAX(fabs(points[ii]->x[0]),fabs(points[ii]->x[1])) ); // If cells are too small they will cause problems.

		points[ii]->leaf->refined = true;
		xygridpoints(&i_points[Nadded],points[ii]->gridsize*(Ngrid_block-1)/Ngrid_block
				,points[ii]->x,Ngrid_block,1);
		points[ii]->gridsize /= Ngrid_block;
		points[ii]->image->gridsize /= Ngrid_block;

    assert(points[ii]->gridsize > 0.0);
    
		// take out points that are outside of original grid
		Nout = 0;
		//if( (points[ii]->x[0] <= i_tree->getTop()->boundary_p1[0]) || (points[ii]->x[0] >= i_tree->getTop()->boundary_p2[0])
		//	|| (points[ii]->x[1] <= i_tree->getTop()->boundary_p1[1]) || (points[ii]->x[1] >= i_tree->getTop()->boundary_p2[1]) ){
      
    Point *point;
    // remove the points that are outside initial image grid
    for(kk=0,Nout=0;kk < (Ngrid_block*Ngrid_block-1);++kk){
      point = &i_points[Nadded + kk - Nout];
      if( !inbox(point->x,i_tree->getTop()->boundary_p1,i_tree->getTop()->boundary_p2)
           || point->gridsize < 1.0e-10*MAX(fabs(point->x[0]),fabs(point->x[1])) ){
        
					// This maintains the ordering in parent cells, but is rather inefficient
        for(unsigned long nn = Nadded + kk - Nout
							; nn < (Ngrid_block*Ngrid_block-1)*Nleaves - 1 - Nout - Nout_tot ; ++nn){
          assert(nn+1 < i_points[0].head);
          SwapPointsInArray(&i_points[nn],&i_points[nn + 1]);
        }
        ++Nout;
					//std::cout << "Nout_tot = " << Nout_tot << std::endl;
      }
    }

		Nout_tot += Nout;
		Nadded += Ngrid_block*Ngrid_block-1 - Nout;
		addedtocell[ii] = Ngrid_block*Ngrid_block-1 - Nout;
    
    if(addedtocell[ii] == 0){
      points[ii]->leaf->refined = false;
      points[ii]->gridsize *= Ngrid_block;
      points[ii]->image->gridsize *= Ngrid_block;
    }
		//if(Nout > 0) i_points = AddPointToArray(i_points,Ngrid_block*Ngrid_block-1-Nout,Ngrid_block*Ngrid_block-1);
	} 

	assert(Nadded == (Ngrid_block*Ngrid_block-1)*Nleaves-Nout_tot);
  
	if(Nadded == 0){
		FreePointArray(i_points,true);
		return NULL;
	}

	s_points = LinkToSourcePoints(i_points,Nadded);

  // This interpolation does not work for some reason and needs to be fixed some time.
	/*/ Here the points that are in uniform magnification regions are calculated by interpolation and marked with in_image = MAYBE
	{
		PosType aa[4],dx[2];
		//PosType a1[4],ss;
		//bool tmp;
		//Point *i_point = NewPointArray(1);
		//Point *s_point = LinkToSourcePoints(i_point,1);

		for(ii=0,kk=0;ii<Nleaves;++ii){

			if(uniform_mag_from_deflect(aa,points[ii])){

				 //tmp = uniform_mag_from_shooter(a1,points[ii]);

				//std::cout << "shooter    a = " << a1[0] << " " << a1[1] << " " << a1[2] << " " << a1[3] << std::endl;
				//std::cout << "deflection a = " << a2[0] << " " << a2[1] << " " << a2[2] << " " << a2[3] << std::endl;
				//std::cout << "          da = " << (a2[0]-a1[0]) << " " << (a2[1]-a1[1]) << " " << (a2[2]-a1[2])
				//		<< " " << (a2[3]-a1[3]) << std::endl;

				for(unsigned long jj = kk; jj < addedtocell[ii] + kk; ++jj){
					dx[0] = i_points[jj].x[0] - points[ii]->x[0];
					dx[1] = i_points[jj].x[1] - points[ii]->x[1];

					i_points[jj].image->x[0] = points[ii]->image->x[0] + aa[0]*dx[0] + aa[2]*dx[1];
					i_points[jj].image->x[1] = points[ii]->image->x[1] + aa[3]*dx[0] + aa[1]*dx[1];

					//i_point->x[0] = i_points[jj].x[0];
					//i_point->x[1] = i_points[jj].x[1];
					//lens->rayshooterInternal(1,i_point,false);

					//ss = sqrt( pow(i_point->image->x[0] - points[ii]->image->x[0],2)
					//		+ pow(i_point->image->x[1] - points[ii]->image->x[1],2) );

					//std::cout << (i_point->image->x[0] - i_points[jj].image->x[0])/ss << "  "
					//		  << (i_point->image->x[1] - i_points[jj].image->x[1])/ss  << std::endl;


					
						i_points[jj].kappa = i_points[jj].image->kappa = 1-0.5*(aa[0]+aa[1]);
						i_points[jj].gamma[0] = i_points[jj].image->gamma1() = -0.5*(aa[0]-aa[1]);
						i_points[jj].gamma[1] = i_points[jj].image->gamma2() = -0.5*(aa[2]+aa[3]);
						i_points[jj].gamma[2] = i_points[jj].image->gamma[2] = -0.5*(aa[2]-aa[3]);
						i_points[jj].invmag = i_points[jj].image->invmag = (1-i_points[jj].kappa)*(1-i_points[jj].kappa)
								     - i_points[jj].gamma[0]*i_points[jj].gamma[0]
								     - i_points[jj].gamma[1]*i_points[jj].gamma[1]
								     + i_points[jj].gamma[2]*i_points[jj].gamma[2];
					
					i_points[jj].in_image = MAYBE;

				}

				//std::cout << " " << std::endl;
			}
			kk += addedtocell[ii];
		}

		//FreePointArray(i_point);
		//FreePointArray(s_point);

	}*/

  {
    std::lock_guard<std::mutex> hold(grid_mutex);
    lens->rayshooterInternal(Nadded,i_points);
  }

	/*********************** TODO test line *******************************
	for(ii=0;ii<Nadded;++ii){
		assert(i_points[ii].image->image == &i_points[ii]);
		assert(s_points[ii].image->image == &s_points[ii]);
	}
	/ ******************************************************/
	// remove the points that are outside initial source grid
	int j,Noutcell;
	for(ii=0,kk=0,Nout=0;ii<Nleaves;++ii){
    if(addedtocell[ii] == 0) continue;
		for(j = 0,Noutcell=0; j < addedtocell[ii]; ++j){
			if( !inbox(s_points[kk - Nout].x,s_tree->getTop()->boundary_p1,s_tree->getTop()->boundary_p2) ){
				//SwapPointsInArray(&i_points[kk - Nout],&i_points[Nadded - 1 - Nout]);
				//SwapPointsInArray(&s_points[kk - Nout],&s_points[Nadded - 1 - Nout]);
				for(long nn = kk - Nout; nn < Nadded - 1 - Nout;++nn){
					assert(nn+1 < s_points[0].head);
					assert(nn+1 < i_points[0].head);
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
		if( !inbox(s_points[kk - Nout].x,s_tree->getTop()->boundary_p1,s_tree->getTop()->boundary_p2) ){
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
  assert(i_points->head == Nadded);
  //***** TODO: test line
  //for(long jj=0;jj<Nadded;++jj){ assert(inbox(i_points[jj].x,i_tree->getTop()->boundary_p1,i_tree->getTop()->boundary_p2));}
  //for(long jj=0;jj<Nadded;++jj){ assert(i_points[jj].gridsize > 0 );}
  //**************************/
 
    //*****************************************************************************/
    //*** these could be mode more efficient by starting at the current in tree
    i_tree->AddPointsToTree(i_points,i_points->head);
    s_tree->AddPointsToTree(s_points,s_points->head);
  
  TreeStruct::iterator i_tree_it(i_tree);
  TreeStruct::iterator s_tree_it(s_tree);

    assert(s_points->head > 0);

    //********* This repairs the leaf pointer which for some reason does not point to a leaf for s_points on some occasions ***********
    for(ii=0;ii < s_points->head;++ii){

      if(i_points[ii].leaf->child1 != NULL || i_points[ii].leaf->child2 != NULL){
        //std::cout << ii << std::endl;
        //i_points[ii].print();
        //i_points[ii].leaf->print();
        i_tree_it = i_points[ii].leaf;
        i_tree->_FindLeaf(i_tree_it,i_points[ii].x,0);
        assert((*i_tree_it)->npoints == 1);
        i_points[ii].leaf = *i_tree_it;
        assert(i_points[ii].prev != NULL || i_points[ii].next != NULL || s_points->head == 1);
        assert(i_points[ii].image->prev != NULL || i_points[ii].image->next != NULL || s_points->head == 1);
      }

      if(s_points[ii].leaf->child1 != NULL || s_points[ii].leaf->child2 != NULL){
        //std::cout << ii << std::endl;
        //s_points[ii].print();
        //s_points[ii].leaf->print();
        s_tree_it = s_points[ii].leaf;
        s_tree->_FindLeaf(s_tree_it,s_points[ii].x,0);
        assert((*s_tree_it)->npoints == 1);
        s_points[ii].leaf = *s_tree_it;
        assert(s_points[ii].prev != NULL || s_points[ii].next != NULL || s_points->head == 1);
        assert(s_points[ii].image->prev != NULL || s_points[ii].image->next != NULL || s_points->head == 1);
      }
    }
  
  for(ii=0;ii<Nleaves;++ii){
		assert(points[ii]->leaf->child1 == NULL && points[ii]->leaf->child2 == NULL);
		if(points[ii]->image->leaf->child1 != NULL || points[ii]->image->leaf->child2 != NULL){
      s_tree_it = points[ii]->image->leaf;
      s_tree_it.movetop();
      //assert(inbox(points[ii]->x,s_tree->getTop()->boundary_p1,s_tree->getTop()->boundary_p2));
      assert(inbox(points[ii]->x,(*s_tree_it)->boundary_p1,(*s_tree_it)->boundary_p2));
      s_tree->_FindLeaf(s_tree_it,points[ii]->x,0);
      assert((*s_tree_it)->npoints == 1);
      points[ii]->image->leaf = *s_tree_it;
    }
  }
	//*********************************************************************

	// This loop should not be necessary!! It is repairing the leaf that has been assigned incorrectly somewhere
	//if(!inbox(point[ii].x,(*i_tree_it)->boundary_p1,(*i_tree_it)->boundary_p2) ) moveTop(i_tree);
	/*for(ii=0;ii < Nleaves;++ii){

		// re-assign leaf of point that was to be refined
		assert(inbox(points[ii]->x,i_tree->getTop()->boundary_p1,i_tree->getTop()->boundary_p2));
		(*i_tree_it) = points[ii]->leaf;
		assert(inbox(points[ii]->x,(*i_tree_it)->boundary_p1,(*i_tree_it)->boundary_p2));
		_FindLeaf(i_tree,points[ii]->x,0);
		points[ii]->leaf = (*i_tree_it);

		assert(inbox(points[ii]->image->x,s_tree->getTop()->boundary_p1,s_tree->getTop()->boundary_p2));
		assert(inbox(points[ii]->image->x,s_tree->getTop()->boundary_p1,s_tree->getTop()->boundary_p2));
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

  if(i_tree->pointlist->size() > 0){
    PointList::iterator i_tree_pointlist_it;
    i_tree_pointlist_it.current = (i_tree->pointlist->Top());
    
    for(i=0;i<i_tree->pointlist->size();++i){
      (*i_tree_pointlist_it)->in_image=NO;
      (*i_tree_pointlist_it)->image->in_image=NO;
      --i_tree_pointlist_it;
    }
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
bool Grid::find_mag_matrix(PosType *a,Point *p0,Point *p1,Point *p2){
	PosType y1[2],y2[2],x1[2],x2[2];

	x1[0] = p1->x[0] - p0->x[0];
	x1[1] = p1->x[1] - p0->x[1];
	x2[0] = p2->x[0] - p0->x[0];
	x2[1] = p2->x[1] - p0->x[1];

	y1[0] = p1->image->x[0] - p0->image->x[0];
	y1[1] = p1->image->x[1] - p0->image->x[1];
	y2[0] = p2->image->x[0] - p0->image->x[0];
	y2[1] = p2->image->x[1] - p0->image->x[1];

	PosType det = x1[0]*x2[1] - x1[1]*x2[0];
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
                                    PosType *a           /// Returned magnification matrix, not specified if returning false
                                    ,Point *point       /// point to be tested
){
  Point *point2;
  PosType ao[4];
  int count=0;
  
  i_tree->FindAllBoxNeighborsKist(point,neighbors);
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

double Grid::mag_from_deflect(Point *point       /// point to be tested
) const {
 
  std::vector<Point *> neighbors;
  i_tree->FindAllBoxNeighborsKist(point,neighbors);
 
  Point_2d i_center = *point;
  std::sort(neighbors.begin(),neighbors.end()
            ,[&i_center](Point *p1,Point *p2){return atan2((*p1)[1]-i_center[1],(*p1)[0]-i_center[0]) < atan2((*p2)[1]-i_center[1],(*p2)[0]-i_center[0]);});
 
  Point_2d s_center = *(point->image);
 
  double image_area = 0;
  int n=neighbors.size();
  for(int i=1 ; i < n ; ++i){
    image_area += ( ( *(neighbors[i-1]) - i_center)^( *(neighbors[i]) - i_center) ) / 2;
  }
  image_area += ( ( *(neighbors.back()) - i_center)^( *(neighbors[0]) - i_center) ) / 2;
 
  assert(image_area > 0);
  
  // switch to source plane
  for(int i=0 ; i< n ; ++i) neighbors[i] = neighbors[i]->image;
  
  int npos=0,nneg=0;
  for(int i=0 ; i< n ; ++i){
    if(neighbors[i]->invmag() >= 0){
      ++npos;
    }else{
      ++nneg;
    }
  }
 
  double source_area = 0;
  if(npos==0 || nneg==0){
 
    for(int i=1 ; i < n ; ++i){
        source_area += ( ( *(neighbors[i-1]) - s_center)^( *(neighbors[i]) - s_center) ) / 2;
    }
    source_area += ( ( *(neighbors.back()) - s_center)^( *(neighbors[0]) - s_center) ) / 2;
  }else{
 
    double neg_area = 0,pos_area = 0;
    
    for(int i=1 ; i < n ; ++i){
      if(neighbors[i]->invmag() * neighbors[i-1]->invmag() > 0 ){
    
        if(neighbors[i]->invmag() > 0){
          pos_area += ( ( *(neighbors[i-1]) - s_center)^( *(neighbors[i]) - s_center) ) / 2;
        }else{
          neg_area += ( ( *(neighbors[i-1]) - s_center)^( *(neighbors[i]) - s_center) ) / 2;
        }
      }
    }
    
    source_area = fabs(neg_area) + fabs(pos_area);
  }
 
  assert(source_area != 0);
  return fabs(image_area / source_area);
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
//bool Grid::uniform_mag_from_shooter(
//                                    PosType *a           /// Returned magnification matrix, not specified if returning false
//                                    ,Point *point       /// point to be tested
//){
//  Point *point2;
//  
//  i_tree->FindAllBoxNeighborsKist(point,neighbors);
//  assert(neighbors->Nunits() > 1);
//  neighbors->MoveToTop();
//  do{
//    point2 = neighbors->getCurrent();
//    if( !( (fabs(point->kappa() - point2->kappa()) < maglimit)*(fabs(point->gamma1()-point2->gamma1()) < maglimit)
//          *(fabs(point->gamma2()-point2->gamma2()) < maglimit)*(fabs(point->gamma3()-point2->gamma3()) < maglimit) ) ) return false;
//  }while(neighbors->Down());
//  
//  neighbors->MoveToTop();
//  std::cout << "shooter neighbors " << std::endl;
//  do{
//    point2 = neighbors->getCurrent();
//    std::cout << point->kappa() - point2->kappa() << "  " << point->gamma1()-point2->gamma1() << "   " <<
//    point->gamma2()-point2->gamma2() << "  " << point->gamma3()-point2->gamma3() << std::endl;
//  }while(neighbors->Down());
//  
//  
//  // this gamma1 = -gamma1 as defined in Matrix2x2 class
//  a[0] = 1 - point->kappa() - point->gamma1();
//  a[1] = 1 - point->kappa() + point->gamma1();
//  a[2] = -point->gamma2() - point->gamma3();
//  a[3] = -point->gamma2() + point->gamma3();
//  
//  return true;
//}
//void Grid::test_mag_matrix(){
//	PosType aa[4];
//	Point *point2;
//	PosType ao[4];
//	int count;
//	PosType kappa, gamma[3], invmag;
//
//  PointList::iterator i_tree_pointlist_it;
//  i_tree_pointlist_it.current = (i_tree->pointlist->Top());
//	do{
//		count=0;
//		i_tree->FindAllBoxNeighborsKist(*i_tree_pointlist_it,neighbors);
//		assert(neighbors->Nunits() >= 3);
//		neighbors->MoveToTop();
//		point2 = neighbors->getCurrent();
//		neighbors->Down();
//		while(!find_mag_matrix(aa,*i_tree_pointlist_it,point2,neighbors->getCurrent())) neighbors->Down();
//		//std::cout << "deflection neighbors a" << std::endl;
//		while(neighbors->Down()){
//			if(find_mag_matrix(ao,*i_tree_pointlist_it,point2,neighbors->getCurrent())){
//				aa[0] = (count*aa[0] + ao[0])/(count+1);
//				aa[1] = (count*aa[1] + ao[1])/(count+1);
//				aa[2] = (count*aa[2] + ao[2])/(count+1);
//				aa[3] = (count*aa[3] + ao[3])/(count+1);
//
//				++count;
//			}
//		}
//
//		kappa = 1-0.5*(aa[0]+aa[1]);
//		gamma[0] = -0.5*(aa[0]-aa[1]);
//		gamma[1] = -0.5*(aa[2]+aa[3]);
//		gamma[2] = -0.5*(aa[2]-aa[3]);
//		invmag = (1-kappa)*(1-kappa) - gamma[0]*gamma[0] - gamma[1]*gamma[1] + gamma[2]*gamma[2];
//
//
//    //(*i_tree_pointlist_it)->kappa =  kappa/(*i_tree_pointlist_it)->kappa - 1.0;
//    //(*i_tree_pointlist_it)->gamma[0] = gamma[0]/(*i_tree_pointlist_it)->gamma[0] - 1.0;
//    //(*i_tree_pointlist_it)->gamma[1] = gamma[1]/(*i_tree_pointlist_it)->gamma[1] - 1.0;
//    //(*i_tree_pointlist_it)->gamma[2] = gamma[2]/(*i_tree_pointlist_it)->gamma[2] - 1.0;
//    //(*i_tree_pointlist_it)->invmag = invmag/(*i_tree_pointlist_it)->invmag - 1.0;
//
//    //(*i_tree_pointlist_it)->invmag = invmag/(*i_tree_pointlist_it)->invmag - 1.0;
//
//	}while(--i_tree_pointlist_it);
//}

/**
 * \brief quickly refines the grid down to a specific scale at a given point
 *
 *   top is an optional argument that allows for the zooming to start part way
 *   down the tree.  Default is to start at the root.  If the point is not within
 *   top or the root nothing is done.  The point will not necessarily be in the center
 *   of the smallest branch.
 */
void Grid::zoom(
		LensHndl lens
		,PosType *center      /// center of point where grid is refined
		,PosType min_scale    /// the smallest grid size to which the grid is refined
		,Branch *top         /// where on the tree to start, if NULL it will start at the root
		){

  TreeStruct::iterator i_tree_it(i_tree);
      
	if(top==NULL){
		if(!inbox(center,i_tree->getTop()->boundary_p1,i_tree->getTop()->boundary_p2)) return;
		i_tree_it.movetop();
	}else{
		if(!inbox(center,top->boundary_p1,top->boundary_p2)) return;
		i_tree_it = top;
	}
	Branch *tmp=NULL;
	Point *newpoints;

	i_tree->_FindLeaf(i_tree_it,center,0);
	while((*i_tree_it)->points->gridsize > min_scale ){
		tmp = *i_tree_it;
		newpoints = RefineLeaf(lens,(*i_tree_it)->points);
		if(newpoints == NULL) break;  // case where all the new points where outside the region
		i_tree_it = tmp;
		assert(inbox(center,(*i_tree_it)->boundary_p1,(*i_tree_it)->boundary_p2));
		i_tree->_FindLeaf(i_tree_it,center,0);
		assert((*i_tree_it)->npoints == 1);
	};

	return;
}

/// Outputs a fits image of a lensing variable of choice
void Grid::writeFits(
                     const PosType center[]     /// center of image
                     ,size_t Npixels           /// number of pixels in image in on dimension
                     ,PosType resolution        /// resolution of image in radians
                     ,LensingVariable lensvar  /// which quantity is to be displayed
                     ,std::string filename     /// file name for image -- .kappa.fits, .gamma1.fits, etc will be appended
                     ){
  writeFits(center,Npixels,Npixels,resolution,lensvar,filename);
}
  /// Outputs a fits image of a lensing variable of choice
void Grid::writeFits(
                     const PosType center[]     /// center of image
                     ,size_t Nx           /// number of pixels in image in x dimension
                     ,size_t Ny           /// number of pixels in image in y dimension
                     ,PosType resolution        /// resolution of image in radians
                     ,LensingVariable lensvar  /// which quantity is to be displayed
                     ,std::string filename     /// file name for image -- .kappa.fits, .gamma1.fits, etc will be appended
                     ){
  PixelMap map(center, Nx,Ny, resolution);
  std::string tag;
  
  switch (lensvar) {
    case DELAYT:
      tag = ".dt.fits";
      break;
    case ALPHA1:
      tag = ".alpha1.fits";
      break;
    case ALPHA2:
      tag = ".alpha2.fits";
      break;
    case ALPHA:
      tag = ".alpha.fits";
      break;
    case KAPPA:
      tag = ".kappa.fits";
      break;
    case GAMMA1:
      tag = ".gamma1.fits";
      break;
    case GAMMA2:
      tag = ".gamma2.fits";
      break;
    case GAMMA3:
      tag = ".gamma3.fits";
      break;
    case GAMMA:
      tag = ".gamma.fits";
      break;
    case INVMAG:
      tag = ".invmag.fits";
      break;
    default:
      break;
  }
  
  map.AddGrid(*this,lensvar);
  map.printFITS(filename + tag);

  return;
}

/// Outputs a PixelMap of the lensing quantities of a fixed grid
PixelMap Grid::writePixelMap(
                             const PosType center[]     /// center of image
                             ,size_t Npixels           /// number of pixels in image in on dimension
                             ,PosType resolution        /// resolution of image in radians
                             ,LensingVariable lensvar  /// which quantity is to be displayed
                             ){
  
  return writePixelMap(center,Npixels,Npixels,resolution,lensvar);
}
/// Outputs a PixelMap of the lensing quantities of a fixed grid
PixelMap Grid::writePixelMap(
                             const PosType center[]     /// center of image
                             ,size_t Nx           /// number of pixels in image in on dimension
                             ,size_t Ny           /// number of pixels in image in on dimension
                             ,PosType resolution        /// resolution of image in radians
                             ,LensingVariable lensvar  /// which quantity is to be displayed
){
  PixelMap map(center, Nx, Ny, resolution);
  map.AddGrid(*this,lensvar);
  
  return map;
}
/// Outputs a PixelMap of the lensing quantities of a fixed grid
PixelMap Grid::writePixelMap(
                             LensingVariable lensvar  /// which quantity is to be displayed
){
  
  Branch *branch = i_tree->getTop();
  double resolution = (branch->boundary_p2[0] - branch->boundary_p1[0])/Ngrid_init;
  PixelMap map(branch->center, Ngrid_init, Ngrid_init2, resolution);
  map.AddGrid(*this,lensvar);
  
  return map;
}

PixelMap  Grid::MapSurfaceBrightness(double resolution){
  Branch *branch = i_tree->getTop();
  int Nx = (int)( (branch->boundary_p2[0] - branch->boundary_p1[0])/resolution );
  int Ny = (int)( (branch->boundary_p2[1] - branch->boundary_p1[1])/resolution );
  
  PixelMap map(branch->center,Nx,Ny,resolution);
  map.AddGridBrightness(*this);

  return map;
}

/** \brief Make a fits map that is automatically centered on the grid and has approximately the same range as the grid.  Nx can be used to change the resolution.  Nx = grid.getInitNgrid() will give the initial grid resolution
 */
void Grid::writePixelFits(
                         size_t Nx           /// number of pixels in image in x dimension
                         ,LensingVariable lensvar  /// which quantity is to be displayed
                         ,std::string filename     /// file name for image -- .kappa.fits, .gamma1.fits, etc will be appended
){
  
  Point_2d center = getInitCenter();
  PosType resolution =  getInitRange()/Nx;
  size_t Ny = (size_t)(Nx*axisratio);
  writeFits(center.x,Nx,Ny,resolution, lensvar, filename);
  
  return;
}

/** \brief Output a fits map of the without distribution the pixels.
 *
 *  This will be faster than Grid::writePixelMap() and Grid::writeFits().
 *  But it puts each grid pixel in one pixelmap pixel and if there are two
 *  grid pixels in one pixelmap pixel it uses one at random.  This is meant
 *  for uniform maps to make equal sized PixelMaps.
 */
void Grid::writeFitsUniform(
                                const PosType center[]  /// center of image
                                ,size_t Nx       /// number of pixels in image in on dimension
                                ,size_t Ny       /// number of pixels in image in on dimension
                                ,LensingVariable lensvar  /// which quantity is to be displayed
                                ,std::string filename     /// file name for image -- .kappa.fits, .gamma1.fits, etc will be appended
                                ){
  std::string tag;
  
  switch (lensvar) {
    case DELAYT:
      tag = ".dt.fits";
      break;
    case ALPHA1:
      tag = ".alpha1.fits";
      break;
    case ALPHA2:
      tag = ".alpha2.fits";
      break;
    case ALPHA:
      tag = ".alpha.fits";
      break;
    case KAPPA:
      tag = ".kappa.fits";
      break;
    case GAMMA1:
      tag = ".gamma1.fits";
      break;
    case GAMMA2:
      tag = ".gamma2.fits";
      break;
    case GAMMA3:
      tag = ".gamma3.fits";
      break;
    case GAMMA:
      tag = ".gamma.fits";
      break;
    case INVMAG:
      tag = ".invmag.fits";
      break;
    default:
      break;
  }

  PixelMap map = writePixelMapUniform(center,Nx,Ny,lensvar);
  map.printFITS(filename + tag);
}

/** \brief Make a Pixel map of the without distribution the pixels.
 *
 *  This will be faster than Grid::writePixelMap() and Grid::writeFits().
 *  But it puts each grid pixel in one pixelmap pixel and if there are two 
 *  grid pixels in one pixelmap pixel it uses one at random.  This is meant 
 *  for uniform maps to make equal sized PixelMaps.
 */
PixelMap Grid::writePixelMapUniform(
                                    const PosType center[]  /// center of image
                                    ,size_t Nx       /// number of pixels in image in on dimension
                                    ,size_t Ny       /// number of pixels in image in on dimension
                                    ,LensingVariable lensvar  /// which quantity is to be displayed
                                    ){
  
  if(getNumberOfPoints() ==0 ) return PixelMap();
  PixelMap map(center, Nx, Ny,i_tree->pointlist->Top()->gridsize);
  map.Clean();
  
  //int Nblocks = Utilities::GetNThreads();
  int Nblocks = 16;
                                      
  //std::vector<PointList> lists(Nblocks);
                                      
  std::vector<Point *> heads(Nblocks);
  std::vector<size_t> sizes(Nblocks,0);
  
  bool allowDecent;
  TreeStruct::iterator i_tree_it(i_tree);
  int i = 0;
                                      
  do{
    if((*i_tree_it)->level == 4){
      
      heads[i] = (*i_tree_it)->points;
      sizes[i] = (*i_tree_it)->npoints;
      
      //lists[i].setTop( (*i_tree_it)->points );
      //lists[i].setN( (*i_tree_it)->npoints );
      ++i;
      allowDecent = false;
    }else{
      allowDecent = true;
    }
  }while(i_tree_it.TreeWalkStep(allowDecent) && i < Nblocks);
  
  std::vector<std::thread> thrs;
  for(int ii = 0; ii < i ;++ii){
  //writePixelMapUniform_(heads[ii],sizes[ii],&map,lensvar);
  //thrs.push_back(std::thread(&Grid::writePixelMapUniform_,this,lists[ii],&map,lensvar));

    thrs.push_back(std::thread(&Grid::writePixelMapUniform_,this,heads[ii],sizes[ii],&map,lensvar));
  }
  for(int ii = 0; ii < i ;++ii) thrs[ii].join();
  
  return map;
}

void Grid::writePixelMapUniform(
                                    PixelMap &map
                                    ,LensingVariable lensvar  /// which quantity is to be displayed
                                    ){
  
  if(getNumberOfPoints() ==0 ) return;
  
  map.Clean();
  int Nblocks = 16;
  //std::vector<PointList> lists(Nblocks);
  TreeStruct::iterator i_tree_it(i_tree);

  std::vector<Point *> heads(Nblocks);
  std::vector<size_t> sizes(Nblocks,0);

  
  bool allowDecent;
  i_tree_it.movetop();
  int i = 0;
  do{
    if((*i_tree_it)->level == 4){
      assert(i < 16);
      //lists[i].setTop( (*i_tree_it)->points );
      //lists[i].setN( (*i_tree_it)->npoints );
      heads[i] = (*i_tree_it)->points;
      sizes[i] = (*i_tree_it)->npoints;
      
      ++i;
      allowDecent = false;
    }else{
      allowDecent = true;
    }
  }while(i_tree_it.TreeWalkStep(allowDecent) && i < Nblocks);
  
  std::thread thr[16];
  for(int ii = 0; ii < i ;++ii){
    //thr[ii] = std::thread(&Grid::writePixelMapUniform_,this,lists[ii],&map,lensvar);
    thr[ii] = std::thread(&Grid::writePixelMapUniform_,this,heads[ii],sizes[ii],&map,lensvar);
  }
  for(int ii = 0; ii < i ;++ii) thr[ii].join();

}

//void Grid::writePixelMapUniform_(const PointList &list,PixelMap *map,LensingVariable val){
//  double tmp;
//  PosType tmp2[2];
//  long index;
//
//  PointList::iterator list_it;
//  list_it.current = (list.Top());
//  for(size_t i = 0; i< list.size(); ++i){
//    switch (val) {
//      case ALPHA:
//        tmp2[0] = (*list_it)->x[0] - (*list_it)->image->x[0];
//        tmp2[1] = (*list_it)->x[1] - (*list_it)->image->x[1];
//        tmp = sqrt(tmp2[0]*tmp2[0] + tmp2[1]*tmp2[1]);
//        break;
//      case ALPHA1:
//        tmp = ((*list_it)->x[0] - (*list_it)->image->x[0]);
//        break;
//      case ALPHA2:
//        tmp = ((*list_it)->x[1] - (*list_it)->image->x[1]);
//        break;
//      case KAPPA:
//        tmp = (*list_it)->kappa();
//        break;
//      case GAMMA:
//        tmp2[0] = (*list_it)->gamma1();
//        tmp2[1] = (*list_it)->gamma2();
//        tmp = sqrt(tmp2[0]*tmp2[0] + tmp2[1]*tmp2[1]);
//        break;
//      case GAMMA1:
//        tmp = (*list_it)->gamma1();
//        break;
//      case GAMMA2:
//        tmp = (*list_it)->gamma2();
//        break;
//      case GAMMA3:
//        tmp = (*list_it)->gamma[2];
//        break;
//      case INVMAG:
//        tmp = (*list_it)->invmag;
//        break;
//      case DELAYT:
//        tmp = (*list_it)->dt;
//        break;
//      default:
//        std::cerr << "PixelMap::AddGrid() does not work for the input LensingVariable" << std::endl;
//        throw std::runtime_error("PixelMap::AddGrid() does not work for the input LensingVariable");
//        break;
//        // If this list is to be expanded to include ALPHA or GAMMA take care to add them as vectors
//    }
//
//    index = map->find_index((*list_it)->x);
//    if(index != -1)(*map)[index] = tmp;
//
//    --list_it;
//  }
//}

void Grid::writePixelMapUniform_(Point *head,size_t N,PixelMap *map,LensingVariable val){
  double tmp;
  PosType tmp2[2];
  long index;
  
  Point *ppoint = head;
  
  for(size_t i = 0; i< N; ++i){
    
    switch (val) {
      case ALPHA:
        tmp2[0] = ppoint->x[0] - ppoint->image->x[0];
        tmp2[1] = ppoint->x[1] - ppoint->image->x[1];
        tmp = sqrt(tmp2[0]*tmp2[0] + tmp2[1]*tmp2[1]);
        break;
      case ALPHA1:
        tmp = (ppoint->x[0] - ppoint->image->x[0]);
        break;
      case ALPHA2:
        tmp = (ppoint->x[1] - ppoint->image->x[1]);
        break;
      case KAPPA:
        tmp = ppoint->kappa();
        break;
      case GAMMA:
        tmp2[0] = ppoint->gamma1();
        tmp2[1] = ppoint->gamma2();
        tmp = sqrt(tmp2[0]*tmp2[0] + tmp2[1]*tmp2[1]);
        break;
      case GAMMA1:
        tmp = ppoint->gamma1();
        break;
      case GAMMA2:
        tmp = ppoint->gamma2();
        break;
      case GAMMA3:
        tmp = ppoint->gamma3();
        break;
      case INVMAG:
        tmp = ppoint->invmag();
        break;
      case DELAYT:
        tmp = ppoint->dt;
        break;
      default:
        std::cerr << "PixelMap::AddGrid() does not work for the input LensingVariable" << std::endl;
        throw std::runtime_error("PixelMap::AddGrid() does not work for the input LensingVariable");
        break;
        // If this list is to be expanded to include ALPHA or GAMMA take care to add them as vectors
    }
    
    index = map->find_index(ppoint->x);
    if(index != -1)(*map)[index] = tmp;
    
    ppoint = ppoint->next;
  }
}

/// Outputs a fits file for making plots of vector fields
void Grid::writeFitsVector(
                     const PosType center[]     /// center of image
                     ,size_t Npixels           /// number of pixels in image in on dimension
                     ,PosType resolution        /// resolution of image in radians
                     ,LensingVariable lensvar  /// which quantity is to be displayed
                     ,std::string filename     /// file name for image -- .kappa.fits, .gamma1.fits, etc will be appended
                     ){
  //throw std::runtime_error("Not done yet!");
  
  PosType range = Npixels*resolution,tmp_x[2];
  ImageInfo tmp_image,tmp_image_theta;
  size_t i;
  std::string tag;
  
  i_tree->PointsWithinKist(center,range/sqrt(2.),tmp_image.imagekist,0);
  i_tree->PointsWithinKist(center,range/sqrt(2.),tmp_image_theta.imagekist,0);
  
  std::vector<PosType> tmp_sb_vec(tmp_image.imagekist->Nunits());
  
  for(tmp_image.imagekist->MoveToTop(),i=0;i<tmp_sb_vec.size();++i,tmp_image.imagekist->Down()){
    tmp_sb_vec[i] = tmp_image.imagekist->getCurrent()->surface_brightness;
    switch (lensvar) {
      case ALPHA1:
        tmp_x[0] = tmp_image.imagekist->getCurrent()->x[0]
            - tmp_image.imagekist->getCurrent()->image->x[0];

        tmp_x[1] = tmp_image.imagekist->getCurrent()->x[1]
            - tmp_image.imagekist->getCurrent()->image->x[1];
      
        tmp_image.imagekist->getCurrent()->surface_brightness = sqrt( tmp_x[0]*tmp_x[0] + tmp_x[1]*tmp_x[1]);
        tmp_image_theta.imagekist->getCurrent()->surface_brightness = atan2(tmp_x[1],tmp_x[0]);
            
        tag = ".alphaV.fits";
        break;
      case GAMMA:
        
        tmp_x[0] = tmp_image.imagekist->getCurrent()->gamma1();
        tmp_x[1] = tmp_image.imagekist->getCurrent()->gamma2();

        tmp_image.imagekist->getCurrent()->surface_brightness = sqrt( tmp_x[0]*tmp_x[0] + tmp_x[1]*tmp_x[1]);
        tmp_image_theta.imagekist->getCurrent()->surface_brightness = atan2(tmp_x[1],tmp_x[0])/2;
            
        tag = ".gammaV.fits";
        break;
      default:
        std::cout << "Grid::writeFitsVector() does not support the LensVariable you are using." << std::endl;
        return;
    }
  }
  
  PixelMap map_m(center, Npixels, resolution),map_t(center,Npixels,resolution);
  
  map_m.Clean();
  map_m.AddImages(&tmp_image,1,-1);
  map_m = PixelMap(map_m,4);
  map_m = PixelMap(map_m,1/4.);
    
  map_t.Clean();
  map_t.AddImages(&tmp_image_theta,1,-1);

  map_m.printFITS(filename + tag);
  
  for(tmp_image.imagekist->MoveToTop(),i=0;i<tmp_sb_vec.size();++i,tmp_image.imagekist->Down())
    tmp_image.imagekist->getCurrent()->surface_brightness = tmp_sb_vec[i];
}

void Grid::writeFits(double strech,LensingVariable lensvar ,std::string filename){
  Point_2d center = getInitCenter();
  size_t N1 = (size_t)(strech*Ngrid_init);
  size_t N2 = (size_t)(strech*Ngrid_init2);
  
  writeFits(center.x,N1,N2,getInitRange()/N1,lensvar,filename);
}


