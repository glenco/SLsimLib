/*
 * Code Name:     tree.c                                       
 * Programmer:    R Ben Metcalf
 * Last Revised:  Nov, 2005                                   
 * Discription:  2 way tree data structure with 2 branches
 * Comments:                           
 */

#include "slsimlib.h"

/***** Structs *****/

/* Branch: Private struct, not exported */ 


/************************************************************************
 * NewBranch
 * Returns pointer to new Branch struct.  Initializes children pointers to NULL,
 * and sets data field to input.  Private.
 ************************************************************************/
/** \ingroup ConstructorL2
 *
 * Gives each branch a unique number even if branches are destroed.
 */
/*Branch *NewBranch(Point *points,unsigned long npoints
		  ,PosType boundary_p1[2],PosType boundary_p2[2]
		  ,PosType center[2],int level){

    //Branch *branch;
    int i;
    static unsigned long number = 0;

    //branch = (Branch *)malloc(sizeof(Branch));
    if (!branch){
      ERROR_MESSAGE(); std::cout << "allocation failure in NewBranch()" << std::endl;
      assert(branch);
      exit(1);
    }
    branch->points = points;
    branch->npoints = npoints;

    branch->center[0]=center[0];
    branch->center[1]=center[1];
    branch->level=level;

    for(i=0;i<2;++i){
      branch->boundary_p1[i]= boundary_p1[i];
      branch->boundary_p2[i]= boundary_p2[i];
    }

    branch->number = number;
    ++number;

    branch->child1 = NULL;
    branch->child2 = NULL;
    branch->brother = NULL;
    branch->prev = NULL;
    branch->refined = false;

    //return branch;
}*/

Point::Point():Point_2d(0,0){
  head = 0;
  in_image = NO;
  surface_brightness = 0;
  leaf = nullptr;
  image = nullptr;
  next = prev = nullptr;
  kappa = dt = gridsize = 0;
  gamma[0] = gamma[1] = gamma[2] = 0;
  invmag = 1;
  flag = false;
}


/// print out all member data for testing purposes
void Point::Print(){
	 std::cout << "Point Data : " << std::endl;
	 std::cout << "  next " << next << std::endl;
	 std::cout << "  prev " << prev << std::endl;
	 std::cout << "  image " << image << std::endl;
	 std::cout << "  id " << id << std::endl;
	 std::cout << "  x " << x[0] << "    " << x[1] << std::endl;
	 std::cout << "  head " << head << std::endl;
	 std::cout << "  in_image " << in_image << std::endl;
	 std::cout << "  kappa " << kappa << std::endl;
	 std::cout << "  gamma " << gamma[0] << " " << gamma[1] << " " << gamma[2] << std::endl;
	 std::cout << "  dt " << dt << std::endl;
	 std::cout << "  invmag " << invmag << std::endl;
	 std::cout << "  inverted " << inverted() << std::endl;
	 std::cout << "  gridsize " << gridsize << std::endl;
	 std::cout << "  surface_brightness " << surface_brightness << std::endl;
	 std::cout << "  leaf " << leaf << std::endl;
}

/// print just position and gridsize
std::ostream &operator<<(std::ostream &os, Point const &p) {
  return os << p.x[0] << " " << p.x[1];
}


unsigned long Branch::countID = 0;

Branch::Branch(Point *my_points,unsigned long my_npoints
		  ,PosType my_boundary_p1[2],PosType my_boundary_p2[2]
		  ,PosType my_center[2],int my_level){

    points = my_points;
    npoints = my_npoints;

    center[0]=my_center[0];
    center[1]=my_center[1];
    level=my_level;

    for(int i=0;i<2;++i){
      boundary_p1[i]= my_boundary_p1[i];
      boundary_p2[i]= my_boundary_p2[i];
    }

    number = countID++;

    child1 = NULL;
    child2 = NULL;
    brother = NULL;
    prev = NULL;
    refined = false;
}
Branch::~Branch(){
	neighbors.clear();
}
/// print out all member data for testing purposes
void Branch::print(){
	 std::cout << "Branch Data : " << std::endl;
	 std::cout << "  points " << points << std::endl;
	 std::cout << "  npoints " << npoints << std::endl;
	 std::cout << "  level " << level << std::endl;
	 std::cout << "  center " << center[0] << "    " << center[1] << std::endl;
	 std::cout << "  number " << number << std::endl;
	 std::cout << "  boundary_p1 " << boundary_p1[0] << "    " << boundary_p1[1] << std::endl;
	 std::cout << "  boundary_p1 " << boundary_p2[0] << "    " << boundary_p2[1] << std::endl;
	 std::cout << "  child1 " << child1 << std::endl;
	 std::cout << "  child2 " << child2 << std::endl;
	 std::cout << "  brother " << brother << std::endl;
	 std::cout << "  prev " << prev << std::endl;
	 std::cout << "  refined " << refined << std::endl;
}

/** \ingroup ConstructorL2
*
void FreeBranch(Branch *branch){

    assert( branch != NULL);
    free(branch);

    return;
}*/

  /** \ingroup ConstructorL2
 **/
Point *NewPointArray(
		unsigned long N  /// number of points in array
		){

  if(N <= 0) return NULL;
  Point *points = new Point[N];
  points[0].head = N;

  for(unsigned long i = 1; i < N; i++){
      points[i].head = 0;
  }

  return points;
}

/** \ingroup ConstructorL2
 *
 */
void FreePointArray(Point *array,bool NewXs){
  /* Note: this deallocates positions!! */

  if(array[0].head){
	  //if(NewXs) for(i=0;i<array[0].head;++i) free(array[i].x);
	  //free(array);
	  //if(NewXs) for(i=0;i<array[0].head;++i) delete[] array[i].x;
	  delete[] array;
  }else{
	  ERROR_MESSAGE();
	  std::cout << "ERROR: FreePointArray, miss aligned attempt to free point array" << std::endl;
	  exit(1);
  }
}

std::mutex TreeStruct::mutex;
/**
 *  \brief  Make a new tree and the linked list of points in it.  Does
 *  not build the tree structure.  The other constructor should be used
 *  to build the whole tree.
 */
TreeStruct::TreeStruct(
		Point *xp   /// array of points to be added to the tree
		,unsigned long npoints   /// number of points
		,PosType boundary_p1[2]   /// bottom left hand corner of root
		,PosType boundary_p2[2]   /// upper right hand corner of root
		,PosType center[2]        /// center of root (this could be the center of mass)
		,int my_Nbucket             /// maximum number of points allowed in a leaf
		){
	construct_root(xp,npoints,boundary_p1,boundary_p2,center,my_Nbucket);
}
/// Basic construction of root with all particles in it but no children
void TreeStruct::construct_root(
		Point *xp   /// array of points to be added to the tree
		,unsigned long npoints   /// number of points
		,PosType boundary_p1[2]   /// bottom left hand corner of root
		,PosType boundary_p2[2]   /// upper right hand corner of root
		,PosType center[2]        /// center of root (this could be the center of mass)
		,int my_Nbucket             /// maximum number of points allowed in a leaf
		){
  unsigned long i;

  /*TreeStruct *tree;
  tree = (TreeStruct *)malloc(sizeof(TreeStruct));
  if (!tree){
    ERROR_MESSAGE();
    std::cout << "ERROR: allocation failure in NewTree()" << std::endl;
    exit(1);
  }*/

    /* make linked list of points */
      //pointlist=NewList();
      pointlist= new PointList;
   //EmptyList(pointlist);
      PointList::iterator pointlist_current;
  for(i=0;i<npoints;++i){
    pointlist->InsertPointAfterCurrent(pointlist_current,&xp[i]);
    --pointlist_current;
  }

  top = new Branch(pointlist->Top(),npoints,boundary_p1,boundary_p2
		      ,center,0);

  Nbranches = 1;
  //current = top;

  Nbucket = my_Nbucket;
  //return(tree);
}

/** \ingroup ConstructorL2
 * \brief Free tree and the linked list of points in it.
 */
TreeStruct::~TreeStruct(){

	emptyTree();
	//free(current);
  free(top);
	--Nbranches;

	free(pointlist);
}


/***** Access functions *****/

/************************************************************************
 * isEmpty
 * Returns "true" if the Tree is empty and "false" otherwise.  Exported.
 ************************************************************************/
bool TreeStruct::isEmpty(){

    return(Nbranches == 0);
}

/************************************************************************
 * noChild
 * Returns "true" if the child of the current branch does not exist and "false" otherwise.
 * Exported.
 * Pre: !isEmpty(tree)
 ************************************************************************/
bool TreeStruct::iterator::noChild(){

    if( (current->child1 == NULL) || (current->child2 == NULL) ) return true;
    return false;
}

/************************************************************************
 * offEnd
 * Returns "true" if current is off end and "false" otherwise.  Exported.
 ************************************************************************/
bool TreeStruct::iterator::offEnd(){
    return(current == NULL);
}

bool TreeStruct::iterator::IsSquareBranch(){
	if( fabs(1 - (current->boundary_p2[0] - current->boundary_p1[0])
			    /(current->boundary_p2[1] - current->boundary_p1[1]) )
			< 0.05){
		return true;
	}

	return false;
}

/************************************************************************
 * getCurrent
 * Returns the points of current.  Exported.
 * Pre: !offEnd(tree)
 ************************************************************************
void TreeStruct::getCurrent(Point *points,unsigned long *npoints){

    if( offEnd() ){
    	ERROR_MESSAGE();
    	std::cout << "Tree Error: calling getCurrent() when current is off end" << std::endl;
    	exit(1);
    }

    *npoints=current->npoints;
    points=current->points;

    return;
}*/

/************************************************************************
 * getNbranches
 * Returns the Nbranches of tree.  Exported.
 ************************************************************************/
unsigned long TreeStruct::getNbranches(){

    return(Nbranches);
}

/***** Manipulation procedures *****/

/************************************************************************
 * movePrev
 * Moves current to the branch before it in tree.  This can move current
 * off end.  Exported.
 * Pre: !offEnd(tree)
 ************************************************************************
bool TreeStruct::moveToChild(int child){
    
    assert(current != NULL);

    if(child==1){
      if( current->child1 == NULL ) return false;
      current = current->child1;
      return true;
    }
    if(child==2){
      if( current->child2 == NULL ) return false;
      current = current->child2;
      return true;
    }
    return false;
}*/

/*bool TreeStruct::moveUp(){

    assert(!offEnd());

    if( current == top ) return false;
    assert(current->prev);
    current = current->prev;  // can move off end
    return true;
}*/

/************************************************************************
 * moveToChild
 * Moves current to child branch after it in tree.  This can move current off
 * end.  Exported.
 * Pre: !offEnd(tree)
 ************************************************************************/
/************************************************************************
 * insertAfterCurrent
 * Inserts a new Branch after the current branch in the tree and sets the
 * data field of the new Branch to input.  Exported.
 * Pre: !offEnd(tree)
 ************************************************************************/
/*void insertChildToCurrent(TreeHndl tree,Point *points,unsigned long npoints
			  ,PosType boundary_p1[2],PosType boundary_p2[2]
			  ,PosType center[2],int child){
    
    Branch *branch;

    //branch = NewBranch(points,npoints,boundary_p1,boundary_p2,center
 	//	       ,tree->current->level+1);
    branch = new Branch(points,npoints,boundary_p1,boundary_p2,center
 		       ,tree->current->level+1);

    assert(tree != NULL);
    
    if( offEnd(tree) ){
    	ERROR_MESSAGE();
        std::cout << "Tree Error: calling insertChildToCurrent() when current is off end" << std::endl;
    	exit(1);
    }

    branch->prev = tree->current;

    if(child==1){
      if(tree->current->child1 != NULL){
    	  ERROR_MESSAGE();
          std::cout << "Tree Error: calling insertChildToCurrent() when child1 already exists" << std::endl;
    	  exit(1);
      }
      tree->current->child1 = branch;
      tree->current->child1->brother = tree->current->child2;
    }
    if(child==2){
      if(tree->current->child2 != NULL){
    	  ERROR_MESSAGE();
          std::cout << "Tree Error: calling insertChildToCurrent() when child2 already exists" << std::endl;
          exit(1);
      }
      tree->current->child2 = branch;      
      tree->current->child2->brother = tree->current->brother;
    }

    tree->Nbranches++;

    return;
}*/
void TreeStruct::insertChildToCurrent(Branch *current,Branch *branch,int child){

    assert(branch->boundary_p1[0] >= current->boundary_p1[0]);
    assert(branch->boundary_p1[1] >= current->boundary_p1[1]);
    assert(branch->boundary_p2[0] <= current->boundary_p2[0]);
    assert(branch->boundary_p2[1] <= current->boundary_p2[1]);


    branch->prev = current;

    if(child==1){
      if(current->child1 != NULL){
    	  ERROR_MESSAGE();
          std::cout << "Tree Error: calling insertChildToCurrent() when child1 already exists" << std::endl;
    	  exit(1);
      }
      current->child1 = branch;
      current->child1->brother = current->child2;
    }
    if(child==2){
      if(current->child2 != NULL){
    	  ERROR_MESSAGE();
          std::cout << "Tree Error: calling insertChildToCurrent() when child2 already exists" << std::endl;
          exit(1);
      }
      current->child2 = branch;
      current->child2->brother = current->brother;
    }

    if(branch->npoints > 0){
      PointList::iterator pointlist_current(branch->points);
    	for(unsigned long i=0;i<branch->npoints;++i){
    		(*pointlist_current)->leaf = branch;
        --pointlist_current;
    	}
    }

    Nbranches++;

    return;
}

  /* same as above but takes a branch structure */
/*
void attachChildToCurrent(TreeHndl tree,Branch data,int child){

	  //insertChildToCurrent(tree,data.points,data.npoints,data.boundary_p1,data.boundary_p2,data.center,child);
	  insertChildToCurrent(tree,data,child);
  return;
}*/

void TreeStruct::attachChildrenToCurrent(Branch *current,Branch* child1,Branch* child2){
	// this is an addition that keeps assigns the brother pointers

	assert(current->child1 == NULL);
	insertChildToCurrent(current,child1,1);
	assert(current->child2 == NULL);
	insertChildToCurrent(current,child2,2);

	current->child1->brother = current->child2;
	current->child2->brother = current->brother;

	// update lists of branch neighbors
	Branch *neighbor;
	current->child1->neighbors.push_back(current->child2);
	current->child2->neighbors.push_back(current->child1);
	std::list<Branch *>::iterator jt,it;
	for( it=current->neighbors.begin() ; it != current->neighbors.end() ; ++it){
		neighbor = *it;//current->neighbors[i];
		for(jt=neighbor->neighbors.begin();jt != neighbor->neighbors.end();++jt){
			if(*jt == current){
				if(AreBoxNeighbors(neighbor,current->child1)){
					current->child1->neighbors.push_back(neighbor);
					neighbor->neighbors.insert(jt,current->child1);
				}
				if(AreBoxNeighbors(neighbor,current->child2)){
					current->child2->neighbors.push_back(neighbor);
					neighbor->neighbors.insert(jt,current->child2);
				}
				assert(*jt == current);
				jt = neighbor->neighbors.erase(jt);
				--jt;
			}
		}
	}

	current->neighbors.clear();
	return;
}

/***** Other operations *****/


/************************************************************************
 * printTree
 * Prints the contents of tree to stdout.  The current element, if there
 * is one, will be enclosed in []'s.  Currently, only to be used when the
 * TreeElements is are integers.  Exported.
 ************************************************************************/  
void TreeStruct::printTree(TreeStruct::iterator &current){
  int i;

    printBranch(*current);
  PointList::iterator pointlist_current((*current)->points);
    for(i=0;i<(*current)->npoints;++i){
      std::cout << (*pointlist_current)->id << " " << (*pointlist_current)->x[0] <<
      " " << (*pointlist_current)->x[1] << std::endl;
      --pointlist_current;
    }
    if((*current)->child1 == NULL) return;

    if( ((*current)->boundary_p1[0]==(*current)->boundary_p2[0]) ||
    		((*current)->boundary_p1[0]==(*current)->boundary_p2[0])	){
    	ERROR_MESSAGE();
    	std::cout << "ERROR: zero area branch" << std::endl;
    	exit(0);
    }
  current.down(1);
  printTree(current);

  current.up();

  current.down(2);
    printTree(current);

  current.up();
  
  return;
}

void printBranch(Branch *data){

  std::cout << "******* branch *******" << std::endl;
  std::cout << "level=" << data->level << " number=" << data->number << std::endl;
  std::cout << "center = [" << data->center[0] << "," << data->center[1] << "]" << std::endl;
  std::cout << "p1 = [" << data->boundary_p1[0] << "," << data->boundary_p1[1] << "]" << std::endl;
  std::cout << "p2 = [" << data->boundary_p2[0] << "," << data->boundary_p2[1] << "]" << std::endl;
  std::cout<< "number of points = " << data->npoints << std::endl;
}
/****************************************************************
 *  code for testing tree
 *****************************************************************/

void TreeStruct::checkTree(){

	unsigned long count=0;

  TreeStruct::iterator current(top);

	_checkTree(current,&count);

	if(count != Nbranches){ std::cout << "checkTree did not reach all branches" << std::endl; exit(0);}
  
	return;
}

void TreeStruct::_checkTree(TreeStruct::iterator &current,unsigned long *count){
	int checkBranch(Branch *branch);

	//std::printf("     hello\n");
	++*count;
	if(checkBranch(*current)) exit(1);

	if((*current)->child1 != NULL){
    current.down(1);
		_checkTree(current,count);
    current.up();
	}

	if((*current)->child2 != NULL){
    current.down(2);
		_checkTree(current,count);
    current.up();
	}

    return;
}

// put in test for each branch
int checkBranch(Branch *branch){
   	if(branch->boundary_p1[0] >= branch->boundary_p2[0]
   	  || branch->boundary_p1[1] >= branch->boundary_p2[1] ){
    		std::cout << std::endl << "this branch is screwed" << std::endl;
    		printBranch(branch);
    		PrintPoint(branch->points);
    		return 1;
   	}
   	return 0;
}

Point *AddPointToArray(Point *points,unsigned long N,unsigned long Nold){

  if(N==Nold) return points;

  unsigned long i;
  Point *newpoints;

  if(Nold==0){
	  newpoints = NewPointArray(N);
  }else{
	  if(points[0].head != Nold){ ERROR_MESSAGE(); std::cout << "ERROR: AddPointToArray head not set correctly" << std::endl; exit(0);}
	  newpoints = NewPointArray(N);
	  //for(i=N;i<Nold;++i) free(points[i].x);

	  assert(points);
	  for(i=0;(i<N && i<Nold);++i){
		  PointCopy(&newpoints[i],&points[i]);
	  }
	  //for(i=Nold;i<N;++i) newpoints[i].x = (PosType *) malloc(2*sizeof(PosType));
  }

  FreePointArray(points,false);
  return newpoints;
}

Point *NewPoint(PosType *x,unsigned long id){
  Point *point;

  point=(Point *) malloc(sizeof(Point));
  assert(point);

  point->head = 1;
  point->id=id;
  //point->x=x;
    point->x[0] = x[0];
    point->x[1] = x[1];
  point->in_image=NO;

  if (!point){
    ERROR_MESSAGE(); std::cout << "allocation failure in NewPoint()" << std::endl;
    exit(1);
  }
  point->next = point->prev = NULL;
  point->leaf=NULL;

  return(point);
}

/** SWAPS information in points without changing
* pointers to prev and next
* changes links of image points to follow */
void SwapPointsInArray(Point *p1,Point *p2){
  Point pt;

  if(p1==p2) return;

  PointCopy(&pt,p1);
  PointCopy(p1,p2);
  PointCopy(p2,&pt);
}
/** SWAPS information in points without changing
* pointers to prev and next
* Does not change links of image points to follow */
void SwapPointsInArrayData(Point *p1,Point *p2){
  Point pt;

  if(p1==p2) return;

  PointCopyData(&pt,p1);
  PointCopyData(p1,p2);
  PointCopyData(p2,&pt);
}

/** copies information in point without copying
 * pointers to prev and next, but moving the link
 * to the image point, does not touch head */
void PointCopy(Point *pcopy,Point *pin){
  pcopy->id = pin->id;
  pcopy->image = pin->image;
  pcopy->invmag = pin->invmag;
  pcopy->kappa = pin->kappa;
  pcopy->gamma[0] = pin->gamma[0];
  pcopy->gamma[1] = pin->gamma[1];
  pcopy->gamma[2] = pin->gamma[2];
  pcopy->dt = pin->dt;
  //pcopy->x = pin->x;
  pcopy->x[0] = pin->x[0];
  pcopy->x[1] = pin->x[1];
  pcopy->gridsize = pin->gridsize;
  pcopy->in_image = pin->in_image;
  pcopy->surface_brightness = pin->surface_brightness;
  pcopy->leaf = pin->leaf;

  if((pin->image != NULL) && (pin->image->image == pin)) pin->image->image = pcopy;
}
/** copies information in point without copying
* pointers to prev and next
* does copy image pointer, does not touch head does not 
* copy the ->image->image pointer */
void PointCopyData(Point *pcopy,Point *pin){
  pcopy->id = pin->id;
  pcopy->image = pin->image;
  pcopy->invmag = pin->invmag;
  pcopy->kappa = pin->kappa;
  pcopy->gamma[0] = pin->gamma[0];
  pcopy->gamma[1] = pin->gamma[1];
  pcopy->gamma[2] = pin->gamma[2];
  pcopy->dt = pin->dt;
  //pcopy->x = pin->x;
  pcopy->x[0] = pin->x[0];
  pcopy->x[1] = pin->x[1];
  pcopy->gridsize = pin->gridsize;
  pcopy->in_image = pin->in_image;
  pcopy->surface_brightness = pin->surface_brightness;
  pcopy->leaf = pin->leaf;
}

/*********************************/
/*  point extraction routines */
/*********************************
void TreeStruct::PointsInCurrent(unsigned long *ids,PosType **x){
  Point *point;
  unsigned long i;

  pointlist->current=current->points;

  point=current->points;
  for(i=0;i<current->npoints;++i){
    x[i]=pointlist->current->x;
    ids[i]=pointlist->current->id;
    MoveDownList(pointlist);
  }

  return;
}*/

void PrintPoint(Point *point){
  std::cout << "Point id = " << point->id << std::endl;
  std::cout << "   x= " <<  point->x[0] << " " << point->x[1] << " gridsize = " << point->gridsize << std::endl;
  std::cout << "   kappa= " << point->kappa;
  std::cout << " gamma = " << point->gamma[0] << " " << point->gamma[1];
  std::cout << " invmag " << point->in_image << std::endl;
}
/** \ingroup Constructor
 * Make an array of imageinfo types.
 */
ImageInfo::ImageInfo(){

  imagekist = new Kist<Point>;
  innerborder = new Kist<Point>;
  outerborder = new Kist<Point>;
  
  ShouldNotRefine = 0;
  uniform_mag = unchecked;

  area = area_error = 0.0;
  centroid[0] = centroid[1] = 0.0;
  gridrange[0] = gridrange[1] = gridrange[2] = 0.0;
}

/*ImageInfo *NewImageInfo(int Nimages){
  ImageInfo *imageinfo;
  int i;

  imageinfo=(ImageInfo *) malloc(Nimages*sizeof(ImageInfo));
  assert(imageinfo);

  for(i = 0;i < Nimages; i++)
  {
	imageinfo[i].imagekist = new Kist<Point>;
	imageinfo[i].Npoints=0;
    imageinfo[i].innerborder = new Kist<Point>;
    imageinfo[i].outerborder = new Kist<Point>;
  }

  return imageinfo;
}*/

/** \ingroup Constructor
 * Destructor of imageinfo types.
 */
ImageInfo::~ImageInfo(){
    delete imagekist;
    delete innerborder;
    delete outerborder;
}

/// return to original state after construction
void ImageInfo::Empty(){
  imagekist->Empty();
  innerborder->Empty();
  outerborder->Empty();
  centroid[0] = centroid[1] =0.0;
  area = area_error = 0.0;
  gridrange[0] = gridrange[1] = gridrange[2] = 0.0;
  ShouldNotRefine = 0;
  uniform_mag = unchecked;
}

/**
 * \brief Copy all information about the image including making copies of the imagekist,
 * innerborder and outerborder.  Previous information in the image will be
 */
void ImageInfo::copy(
                     const ImageInfo &image   /// image to be copied
                     ,bool copykists    /// false for turning off copying imagekist, innerborder and outerborder.
                                        /// This is useful when the grid has already been distroyed.
                     ){

  if(copykists){
    imagekist->copy(image.imagekist);
    innerborder->copy(image.innerborder);
    outerborder->copy(image.outerborder);
  }
  
	gridrange[0] = image.gridrange[0];
	gridrange[1] = image.gridrange[1];
	gridrange[2] = image.gridrange[2];

	centroid[0] = image.centroid[0];
	centroid[1] = image.centroid[1];

	area = image.area;
	area_error = image.area_error;
	ShouldNotRefine = image.ShouldNotRefine;
	uniform_mag = image.uniform_mag;
}
/*
void freeImageInfo(ImageInfo *imageinfo,int Nimages){
	int i;

	//freeKist(imageinfo->imagekist);
	for(i=0;i<Nimages;++i){
		delete imageinfo[i].imagekist;
	    delete imageinfo[i].innerborder;
	    delete imageinfo[i].outerborder;
	}
	free(imageinfo);
}
*/

OldImageInfo::OldImageInfo(){

  innerborder = new Kist<Point>;
  outerborder = new Kist<Point>;
}
OldImageInfo::~OldImageInfo(){
    delete innerborder;
    delete outerborder;
}

/** \ingroup LowLevel
 *  step for walking tree by iteration instead of recursion
 *
bool TreeStruct::TreeWalkStep(bool allowDescent){

	if(allowDescent && current->child1 != NULL){
		moveToChild(1);
		return true;
	}
	if(allowDescent && current->child2 != NULL){
		moveToChild(2);
		return true;
	}

	if(current->brother != NULL){
		current = current->brother;
		return true;
	}

	return false;
}*/

void ImageFinding::CriticalCurve::RandomSourceWithinCaustic(int N,std::vector<Point_2d> &y,Utilities::RandomNumbers_NR &rng)
{
  CausticRange(p1,p2);
  y.resize(N);
  for(int ii=0;ii<N;++ii){
    do{
      y[ii][0] = (p2[0] - p1[0])*rng() + p1[0];
      y[ii][1] = (p2[1] - p1[1])*rng() + p1[1];
    }while(!inCausticCurve(y[ii]));
  }
}

/**
 *  Returns N source positions which are not overlapping with the caustic lines.
 *  i.e. the distance between the caustic line and the center of the source is larger than radius of the source.
 */
void ImageFinding::CriticalCurve::RandomSourceStrictlyWithinCaustic(int N,std::vector<Point_2d> &y,Utilities::RandomNumbers_NR &rng, PosType sourceRadius, PosType distSourceToCaustic)
{
  assert(distSourceToCaustic > sourceRadius);
  CausticRange(p1,p2);
  y.resize(N);
  for(int ii=0;ii<N;++ii){
    do{
      y[ii][0] = (p2[0] - p1[0])*rng() + p1[0];
      y[ii][1] = (p2[1] - p1[1])*rng() + p1[1];
    }while(!EntirelyinCausticCurve(y[ii],distSourceToCaustic));
  }
}


void ImageFinding::CriticalCurve::CritRange(Point_2d &my_p1,Point_2d &my_p2){
  if(critical_curve.size() == 0){
    my_p1 *= 0;
    my_p2 *= 0;
    return;
  }
  my_p1 = my_p2 = critical_curve[0];
  
  for(size_t ii=1;ii < critical_curve.size();++ii){
    if(critical_curve[ii][0] < my_p1[0] ) my_p1[0] = critical_curve[ii][0];
    if(critical_curve[ii][0] > my_p2[0] ) my_p2[0] = critical_curve[ii][0];
    if(critical_curve[ii][1] < my_p1[1] ) my_p1[1] = critical_curve[ii][1];
    if(critical_curve[ii][1] > my_p2[1] ) my_p2[1] = critical_curve[ii][1];
  }
}
void ImageFinding::CriticalCurve::CausticRange(Point_2d &my_p1,Point_2d &my_p2){
  if(caustic_curve_outline.size() == 0){
    my_p1 *= 0;
    my_p2 *= 0;
    return;
  }
  my_p1 = my_p2 = caustic_curve_outline[0];
  
  for(size_t ii=1;ii < caustic_curve_outline.size();++ii){
    if(caustic_curve_outline[ii][0] < my_p1[0] ) my_p1[0] = caustic_curve_outline[ii][0];
    if(caustic_curve_outline[ii][0] > my_p2[0] ) my_p2[0] = caustic_curve_outline[ii][0];
    if(caustic_curve_outline[ii][1] < my_p1[1] ) my_p1[1] = caustic_curve_outline[ii][1];
    if(caustic_curve_outline[ii][1] > my_p2[1] ) my_p2[1] = caustic_curve_outline[ii][1];
  }
}




