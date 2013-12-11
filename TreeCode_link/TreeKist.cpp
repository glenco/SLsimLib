/*
 * Code Name:     tree.c                                       
 * Programmer:    R Ben Metcalf
 * Last Revised:  Nov, 2005                                   
 * Discription:  2 way tree data structure with 2 branches
 * Comments:                           
 */

#include "slsimlib.h"
#include "Tree.h"
#include "TreeKist.h"

/**
 *  \brief  Make a new tree and the linked list of points in it.  Does
 *  not build the tree structure.  The other constructor should be used
 *  to build the whole tree.
 */
TreeKist::TreeKist(
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
void TreeKist::construct_root(
		Point *xp   /// array of points to be added to the tree
		,unsigned long npoints   /// number of points
		,PosType boundary_p1[2]   /// bottom left hand corner of root
		,PosType boundary_p2[2]   /// upper right hand corner of root
		,PosType center[2]        /// center of root (this could be the center of mass)
		,int my_Nbucket             /// maximum number of points allowed in a leaf
		){
  unsigned long i;

    /* make linked list of points */
  for(i=0;i<npoints;++i){
      pointkist.InsertAfterCurrent(&xp[i]);
      ++pointkist;
  }

  pointkist.MoveToTop();
  top = new Branch(*pointkist,npoints,boundary_p1,boundary_p2
		      ,center,0);

  Nbranches = 1;
  current = top;

  Nbucket = my_Nbucket;
}

/** \ingroup ConstructorL2
 * \brief Free tree and the linked list of points in it.
 */
TreeKist::~TreeKist(){

	emptyTree();
	free(current);
	--Nbranches;

  pointkist.Clear();
}


/***** Access functions *****/

/************************************************************************
 * isEmpty
 * Returns "true" if the Tree is empty and "false" otherwise.  Exported.
 ************************************************************************/
bool TreeKist::isEmpty(){

    return(Nbranches == 0);
}

/************************************************************************
 * atTop
 * Returns "true" if current is the same as top and "false" otherwise.
 * Exported.
 * Pre: !isEmpty(tree)
 ************************************************************************/
bool TreeKist::atTop(){

    if( isEmpty() ){
	
	ERROR_MESSAGE();
    std::cout << "Tree Error: calling atTop() on empty tree" << std::endl;
	exit(1);
    }
    return(current == top);
}

/************************************************************************
 * noChild
 * Returns "true" if the child of the current branch does not exist and "false" otherwise.
 * Exported.
 * Pre: !isEmpty(tree)
 ************************************************************************/
bool TreeKist::noChild(){

	if( isEmpty() ){
		ERROR_MESSAGE();
		std::cout << "Tree Error: calling atTop() on empty tree" << std::endl;
		exit(1);
    }

    if( (current->child1 == NULL) || (current->child2 == NULL) ) return true;
    return false;
}

/************************************************************************
 * offEnd
 * Returns "true" if current is off end and "false" otherwise.  Exported.
 ************************************************************************/
bool TreeKist::offEnd(){
    return(current == NULL);
}

bool TreeKist::CurrentIsSquareBranch(){
	if( fabs(1 - (current->boundary_p2[0] - current->boundary_p1[0])
			    /(current->boundary_p2[1] - current->boundary_p1[1]) )
			< 0.05){
		return true;
	}

	return false;
}

/************************************************************************
 * getCurrent
 * Returns the first point in branch
 * Pre: !offEnd(tree)
 ************************************************************************/
Kist<Point>::iterator TreeKist::getCurrent(unsigned long *npoints){

    if( offEnd() ){
    	ERROR_MESSAGE();
    	std::cout << "Tree Error: calling getCurrent() when current is off end" << std::endl;
    	exit(1);
    }

    *npoints=current->npoints;

    return current->pointit;
}

/************************************************************************
 * getNbranches
 * Returns the Nbranches of tree.  Exported.
 ************************************************************************/
unsigned long TreeKist::getNbranches(){

    return(Nbranches);
}

/***** Manipulation procedures *****/

/************************************************************************
 * movePrev
 * Moves current to the branch before it in tree.  This can move current
 * off end.  Exported.
 * Pre: !offEnd(tree)
 ************************************************************************/
bool TreeKist::moveToChild(int child){
    
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
}

bool TreeKist::moveUp(){

    assert(!offEnd());
    /*if( offEnd(tree) ){
      ERROR_MESSAGE();
      std::cout << "Tree Error: calling moveUp() when current is off end" << std::endl;
      exit(1);
    }*/

    if( current == top ) return false;
    assert(current->prev);
    current = current->prev;  /* can move off end */
    return true;
}

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
void TreeKist::insertChildToCurrent(Branch *branch,int child){

    assert(branch->boundary_p1[0] >= current->boundary_p1[0]);
    assert(branch->boundary_p1[1] >= current->boundary_p1[1]);
    assert(branch->boundary_p2[0] <= current->boundary_p2[0]);
    assert(branch->boundary_p2[1] <= current->boundary_p2[1]);

    if( offEnd() ){
    	ERROR_MESSAGE();
        std::cout << "Tree Error: calling insertChildToCurrent() when current is off end" << std::endl;
    	exit(1);
    }

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
    	/*pointlist->current = branch->points;
    	for(unsigned long i=0;i<branch->npoints;++i){
    		pointlist->current->leaf = branch;
    		MoveDownList(pointlist);
    	}*/

      //KistIt<Point> it = branch->pointit;
        Kist<Point>::iterator it = branch->pointit;
    	for(unsigned long i=0;i<branch->npoints;++i,--it){
    		(*it)->leaf = branch;
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

void TreeKist::attachChildrenToCurrent(Branch* child1,Branch* child2){
	// this is an addition that keeps assigns the brother pointers

	/*
	insertChildToCurrent(tree,child1.points,child1.npoints
			,child1.boundary_p1,child1.boundary_p2,child1.center,1);
	insertChildToCurrent(tree,child2.points,child2.npoints
			,child2.boundary_p1,child2.boundary_p2,child2.center,2);
*/

	assert(current->child1 == NULL);
	insertChildToCurrent(child1,1);
	assert(current->child2 == NULL);
	insertChildToCurrent(child2,2);

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
void TreeKist::printTree(){
  int i;

    printBranch(current);
/*    pointlist->current = current->points;
    for(i=0;i<current->npoints;++i){
      std::cout << pointlist->current->id << " " << pointlist->current->x[0] << " " << pointlist->current->x[1] << std::endl;
      MoveDownList(pointlist);
    }
*/
    //KistIt<Point> it = current->pointit;
    Kist<Point>::iterator it = current->pointit;
  for(i=0;i<current->npoints;++i,--it){
    std::cout << (*it)->id << " " << (*it)->x[0] << " " << (*it)->x[1] << std::endl;
  }
  
  if(current->child1 == NULL) return;

    if( (current->boundary_p1[0]==current->boundary_p2[0]) ||
    		(current->boundary_p1[0]==current->boundary_p2[0])	){
    	ERROR_MESSAGE();
    	std::cout << "ERROR: zero area branch" << std::endl;
    	exit(0);
    }
    moveToChild(1);
    printTree();

    moveUp();

    moveToChild(2);
    printTree();

    moveUp();

    return;
}/****************************************************************
 *  code for testing tree
 *****************************************************************/

void TreeKist::checkTree(){
	Branch *initial;
	unsigned long count=0;

	initial=current;

	moveTop();
	_checkTree(&count);

	if(count != Nbranches){ std::cout << "checkTree did not reach all branches" << std::endl; exit(0);}
	current=initial;
	return;
}

void TreeKist::_checkTree(unsigned long *count){
	int checkBranch(Branch *branch);

	//std::printf("     hello\n");
	++*count;
	if(checkBranch(current)) exit(1);

	if(current->child1 != NULL){
		moveToChild(1);
		_checkTree(count);
		moveUp();
	}

	if(current->child2 != NULL){
		moveToChild(2);
		_checkTree(count);
		moveUp();
	}

    return;
}

/*********************************/
/*  point extraction routines */
/*********************************/
void TreeKist::PointsInCurrent(unsigned long *ids,PosType **x){
  unsigned long i;

  Kist<Point>::iterator it = current->pointit;

  
  for(i=0;i<current->npoints;--it){
    x[i]=(*it)->x;
    ids[i]=(*it)->id;
  }

  return;
}

/** \ingroup LowLevel
 *  step for walking tree by iteration instead of recursion
 */
bool TreeKist::TreeWalkStep(bool allowDescent){

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
}

/**
 * \brief Build a complete tree from a list of points.
 */

TreeKist::TreeKist(Point *xp,unsigned long Npoints,short my_median_cut,PosType buffer){
  //TreeHndl tree;
  unsigned long i;
  PosType p1[2],p2[2],center[2];
  const int my_Nbucket = 1;   // must be =1 if each leaf is to coincide with each cell
  
  if( (Npoints & (Npoints-1)) != 0){
	  ERROR_MESSAGE();
	  std::printf("ERROR: BuildTree, Npoints is not a power of 2\n");
	  exit(1);
  }
  
  p1[0]=xp[0].x[0]; p1[1]=xp[0].x[1];
  p2[0]=xp[0].x[0]; p2[1]=xp[0].x[1];
  
  for(i=0;i<Npoints;++i){
    
    /* find X boundary */
    if(xp[i].x[0] < p1[0] ) p1[0]=xp[i].x[0];
    if(xp[i].x[0] > p2[0] ) p2[0]=xp[i].x[0];
    
    /* find Y boundary */
    if(xp[i].x[1] < p1[1] ) p1[1]=xp[i].x[1];
    if(xp[i].x[1] > p2[1] ) p2[1]=xp[i].x[1];
  }
  
  center[0]=(p1[0]+p2[0])/2;
  center[1]=(p1[1]+p2[1])/2;
  
  p1[0] -= buffer;
  p1[1] -= buffer;
  p2[0] += buffer;
  p2[1] += buffer;
  
  /* Initialize tree root */
  construct_root(xp,Npoints,p1,p2,center,my_Nbucket);
  
  median_cut = my_median_cut;
  /* build the tree */
  _BuildTree();
  
  moveTop();
}
/** \ingroup  ImageFindingL2
 * \brief Fill a tree with points.  The previous tree structure will be destroyed.  Used for refilling.
 */
void TreeKist::FillTree(Point *xp,unsigned long Npoints){
  unsigned long i;
  
  current->npoints=Npoints;
  // link point array into point list
  pointkist.Clear();

  for(i=0;i<Npoints;++i){
    pointkist.InsertAfterCurrent(&xp[i]);
    pointkist.Down();
  }
  pointkist.MoveToTop();
  current->pointit = pointkist.getTopIt();
  
  // build the tree 
  // make sure there are no branches in the tree
  _freeBranches_iter();
  _BuildTree();
  
  return ;
}

/** \ingroup  ImageFindingL2
 * \brief Rebuilds the tree from the points that are already in the tree->pointlist
 *
 * This is not the best function because it copies all the points
 */
void TreeKist::RebuildTreeFromList(){
  /* Builds or rebuilds a tree that already has a tree->pointlist.
   */
  
	if(pointkist.Nunits() < 1) return;
  
	unsigned long Npoints = pointkist.Nunits(),i;
	PosType *tmp;
  
	// Make new array of points
	Point *points = NewPointArray(Npoints,true);
  
	pointkist.MoveToTop();
    Kist<Point>::iterator it = pointkist.getTopIt();
	for(i=0;i<Npoints;++i,--it){
		tmp = points[i].x;
		// PointCopy() copies the x pointer which is later freed in
		// emptyTree() so the information needs to be copied
    
		PointCopy(&(points[i]),*it);
		points[i].x = tmp;
		points[i].x[0] = (*it)->x[0];
		points[i].x[1] = (*it)->x[1];
    
		/* Below is for trash collection in PruneTree().
		 * When PruneTree has been used not all of the
		 * heads in the point arrays are guaranteed to be
		 * in the pointlist.
		 */
		(*it)->leaf = NULL;
	}
	//std::printf(" %e %e\n",points[0].x[0],points[0].x[1]);
	// emptry the tree and free all former points
	emptyTree();
  
	assert(Nbranches == 1);
	assert(top->npoints == 0);
	assert(pointkist.Nunits() == 0);
  
	FillTree(points,Npoints);
  
	return;
}

/** \brief Spawn a subtree with current as its top
 *
 *  The new tree contains all of the tree below the current.
 *  Warning:: Adding points to the new tree will not update the
 *  parent tree so it can become dangerously out of sync.
 */
TreeKist * TreeKist::spawn(){
  throw std::runtime_error("This is untested and could cause significant problems");
  
  TreeKist *newTree = new TreeKist;
  
  newTree->Nbucket = Nbucket;
  newTree->top = current;
  
  // !!!! TODO:need to cut out the point kist and without making a copy
  Kist<Point>::iterator it;
  size_t i;
  for(it = current->pointit,i=0;i<current->npoints;--it,++i){
      newTree->pointkist.InsertAfterCurrent(*it);
  }
 
  // count the number of branches below
  do{
    newTree->Nbranches++;
  }while(TreeWalkStep(true) && current != newTree->top->brother);
  current = newTree->top;
  
  return newTree;
}

/** \ingroup  ImageFindingL2
 * \brief Empty tree of all point leaving a tree with an empty root.
 *
 * The points are freed, but the list structure is not destroyed.
 *
 *  FillTree can then be used to regenerate tree.
 */
short TreeKist::emptyTree(){
    
  _freeBranches_iter();
  
  assert(Nbranches == 1);
  
  pointkist.Clear();
  
  top->npoints = 0;
  
  return 1;
}

/** \ingroup LowLevel
 * \brief Recursively free branches
 */
void TreeKist::_freeBranches(short child){
	Branch *branch;
  
	if(current->child1 != NULL){
		moveToChild(1);
		_freeBranches(1);
	}
  
  if(current->child2 != NULL){
    moveToChild(2);
    _freeBranches(2);
  }
  
  if( (current->child1 == NULL)*(current->child2 == NULL) ){
    
    if(atTop()){
      return;
    }
    
    branch=current;
    moveUp();
    delete branch;
        
    if(child==1) current->child1=NULL;
    if(child==2) current->child2=NULL;
    --Nbranches;
    
    return;
  }
  
  return;
}
/** \ingroup LowLevel
 * \brief Iteratively free branches
 *
 * Frees all the branches of the tree so there is only the stump.
 */
void TreeKist::_freeBranches_iter(){
	Branch *branch;
  
	moveTop();
  
	/*printBranch(current);*/
  
	while(Nbranches > 1){
    
		if(current->child1 != NULL){
			moveToChild(1);
		}else if(current->child2 != NULL){
			moveToChild(2);
		}else{
			branch = current;
			if(current->brother == current->prev->brother){
				moveUp();
				assert(current->child1 == NULL);
				current->child2 = NULL;
			}else{
				current = current->brother;
				current->prev->child1 = NULL;
			}
			delete branch;
			--Nbranches;
		}
	}
  
	// re-assign the leaf pointers in the particles to the root
  pointkist.MoveToTop();
  do{ pointkist.getCurrent()->leaf = top; }while(pointkist.Down());
  
  return;
}


/** \ingroup LowLevel
 * \brief Recursively build tree from points in its linked list.
 */
void TreeKist::_BuildTree(){
  /* pointlist must be both a linked list and an array of points in the */
  /* same order as the linked list */
  unsigned long i,cut,dimension;
  Branch *cbranch;//,branch1,branch2;
  PosType xcut;
  
  cbranch=current; /* pointer to current branch */
  
  /* leaf case */
  if(cbranch->npoints <= Nbucket){
	  //current->points->leaf = current;
    (*(current->pointit))->leaf = current;
	  if(cbranch->npoints == 0) cbranch->points = NULL;
	  return;
  }
  
	Branch* branch1 = new Branch(NULL,0,cbranch->boundary_p1,cbranch->boundary_p2
                               ,cbranch->center,cbranch->level+1);
	Branch* branch2 = new Branch(NULL,0,cbranch->boundary_p1,cbranch->boundary_p2
                               ,cbranch->center,cbranch->level+1);
  
  
  /* initialize boundaries to old boundaries
   for(i=0;i<2;++i){
   branch1->boundary_p1[i]=cbranch->boundary_p1[i];
   branch1->boundary_p2[i]=cbranch->boundary_p2[i];
   
   branch2->boundary_p1[i]=cbranch->boundary_p1[i];
   branch2->boundary_p2[i]=cbranch->boundary_p2[i];
   }*/
  
  /* set dimension to cut box */
  dimension=(cbranch->level % 2);
  
  PosType *x = new PosType[cbranch->npoints];
  assert(x);
  
  /* reorder points */
  
  Kist<Point>::iterator it = current->pointit;
  for(i=0;i<cbranch->npoints;++i,--it){
    x[i]=(*it)->x[dimension];
    /*points[i]=pointlist->current->id;*/
  }
  
  if(median_cut){
    //double_sort_points(cbranch->npoints,x-1,current->points);
    
    // TODO: needs to be fixed *********************************************************
	  Utilities::quicksortPoints(current->points,x,cbranch->npoints);
    
	  cut=cbranch->npoints/2;
    branch1->boundary_p2[dimension]=(x[cut]+x[cut-1])/2;
    branch2->boundary_p1[dimension]=(x[cut]+x[cut-1])/2;
    
  }else{
    
	  xcut=(cbranch->boundary_p1[dimension]+cbranch->boundary_p2[dimension])/2;
    branch1->boundary_p2[dimension]=xcut;
    branch2->boundary_p1[dimension]=xcut;
    
    // TODO: needs to be fixed *********************************************************
	  Utilities::quickPartitionPoints(xcut,&cut,current->points,x,cbranch->npoints);
    
    //locateD(x-1,cbranch->npoints,xcut,&cut);
  }
  
  
  /* set point numbers and pointers to points */
  branch1->npoints=cut;
  assert((*(current->pointit))->next || (*(current->pointit))->prev);
  branch1->points=current->points;
  branch1->pointit = current->pointit;
  
  branch2->npoints=cbranch->npoints - cut;
  //pointlist->current=current->points;
  pointkist.SetCurrentIt(current->pointit);
  //JumpDownList(pointlist,cut);
  pointkist.JumpDown(cut);
  branch2->points=pointkist.getCurrent();
  branch2->pointit = pointkist.getCurrentIt();
  
  delete[] x;
  
  
	/* use geometric center */
	branch1->center[0] = (branch1->boundary_p1[0] + branch1->boundary_p2[0])/2;
	branch1->center[1] = (branch1->boundary_p1[1] + branch1->boundary_p2[1])/2;
  
	branch2->center[0] = (branch2->boundary_p1[0] + branch2->boundary_p2[0])/2;
	branch2->center[1] = (branch2->boundary_p1[1] + branch2->boundary_p2[1])/2;
  
  /* centers of mass *
   
   for(i=0;i<2;++i) branch1->center[i]=0;
   pointlist->current=branch1->points;
   for(i=0;i<cut; ++i){
   branch1->center[0]+=pointlist->current->x[0]/branch1->npoints;
   branch1->center[1]+=pointlist->current->x[1]/branch1->npoints;
   MoveDownList(pointlist);
   }
   
   for(i=0;i<2;++i) branch2->center[i]=0;
   pointlist->current=branch2->points;
   for(i=cut;i<cbranch->npoints; ++i){
   branch2->center[0]+=pointlist->current->x[0]/branch2->npoints;
   branch2->center[1]+=pointlist->current->x[1]/branch2->npoints;
   MoveDownList(pointlist);
   }*/
  
  attachChildrenToCurrent(branch1,branch2);
  
  if( branch1->npoints > 0 ){
    //attachChildToCurrent(branch1,1);
    moveToChild(1);
    _BuildTree();
    moveUp();
  }
  
  if(branch2->npoints > 0 ){
    //attachChildToCurrent(branch2,2);
	  moveToChild(2);
	  _BuildTree();
	  moveUp();
  }
  
  /*std::printf("reached end of _BuildTree level=%i\n",current->level);*/
  return;
}

/** \ingroup  ImageFindingL2
 *  \brief Expands tree by adding points
 */
int TreeKist::AddPointsToTree(Point *xpoint,unsigned long Nadd){
  unsigned long j,Ntest;
  //Branch *parent_branch;
  
  if(Nadd==0) return 1;
  
  Ntest = pointkist.Nunits();
  
  for(j=0;j<Nadd;++j){
    
    // add only that are inside original grid
    if( !inbox(xpoint[j].x,top->boundary_p1,top->boundary_p2) ){
      ERROR_MESSAGE();
      std::printf("ERROR: in AddPointToTree, ray is not inside the simulation box x = %e %e Nadd=%li\n  not adding it to tree\n",
                  xpoint[j].x[0],xpoint[j].x[1],Nadd);
      //std::printf("root of tree\n");
      //printBranch(top);
      //exit(0);
      return 0;
    }else{
      
      moveTop();
      _FindLeaf(xpoint[j].x,1);
      assert(atLeaf());
      //parent_branch = current->prev->prev;
      
      //assert(inbox(xpoint[j].x,parent_branch->boundary_p1,parent_branch->boundary_p2));
      //assert(inbox(xpoint[j].x,current->boundary_p1,current->boundary_p2));
      
      // insert point into point list
      //if(current->points == NULL){
      if(current->npoints == 1){
        // case of no previous points in leaf
        
        current->points = &xpoint[j];
        current->points->leaf = current;
        // put point into right place in list
        //pointlist->current = current->prev->points;
        pointkist.SetCurrentIt(current->prev->pointit);
        if(current == current->prev->child1){
          Point *oldpoint = pointkist.getCurrent();
          Branch *tmp = current;
          
          // insert the point
          if(pointkist.Up()){
            pointkist.InsertAfterCurrent(current->points);
          }else{
            pointkist.InsertBeforeCurrent(current->points);
          }
          //InsertPointBeforeCurrent(pointlist,current->points);
          // reassign parent first points
          moveUp();
          while(current->points == oldpoint){
            current->points = &xpoint[j];
            moveUp();
          }
          current = tmp;
          
        }else{
          pointkist.JumpDown(current->prev->npoints-2);  // adds point to end of parent branches list
          pointkist.InsertAfterCurrent(current->points);
        }
        
        //pointlist->current = current->points;
        //current->points->leaf = current;
        
        assert(current->points->next || current->points->prev);
      }else{ // case where the leaf needs to be divided
        assert(current->npoints > 1);
        
        // case of points already in leaf
        assert(current->points->next || current->points->prev);
        // Test lines
        //if(!testLeafs()){ERROR_MESSAGE(); std::cout << "before _addPoint of AddPointsToTree "<< std::endl; exit(1);}
        
        // adds point to end of branches list, note that current->npoints has already been increased
        //pointlist->current = current->points;
        pointkist.SetCurrentIt(current->pointit);
        pointkist.JumpDown(current->npoints-2);
        pointkist.InsertAfterCurrent(&xpoint[j]);
        
        //pointlist->current = current->points;
        pointkist.SetCurrentIt(current->pointit);
        xpoint[j].leaf = current;
        
        /*/ Test lines
         pointlist->current = current->points;
         for(int k = 0; k< current->npoints ; ++k,MoveDownList(pointlist)){
         assert(inbox(pointlist->current->x
         ,current->boundary_p1,current->boundary_p2));
         }
         ///////////////////////////////////*/
        
        Branch *parent_branch = current;
        //unsigned long n=pointlist->Npoints;
        _AddPoint();
        //assert(n==pointlist->Npoints);
        
        assert(parent_branch->child1->points->leaf == parent_branch->child1);
        assert(parent_branch->child2->points->leaf == parent_branch->child2);
        assert(inbox(xpoint[j].x,xpoint[j].leaf->boundary_p1,xpoint[j].leaf->boundary_p2));
        /*/  Test lines
         if(!testLeafs()){
         ERROR_MESSAGE();
         pointlist->current = parent_branch->points;
         for(int i=0;i<parent_branch->npoints;++i,MoveDownList(pointlist)){
         std::cout << "points in split branch " << pointlist->current->id << std::endl;
         }
         std::cout << "Adding point " << xpoint[j].id << " branch " << parent_branch << std::endl;
         std::cout << "End of AddPointsToTree "<< std::endl; exit(1);
         }
         //////////////////////////////////////////*/
        
      }
      
    }
    
    while( !inbox(xpoint[j].x,current->boundary_p1,current->boundary_p2) ) moveUp();
    _FindLeaf(xpoint[j].x,0);
    assert(atLeaf());
    xpoint[j].leaf = current;
    assert(inbox(xpoint[j].x,xpoint[j].leaf->boundary_p1,xpoint[j].leaf->boundary_p2));
    /*/ Test lines
     if(!testLeafs()){ERROR_MESSAGE(); std::cout << "End of AddPointsToTree "<< std::endl; exit(1);}
     ///////////////////////*/
  }
  
  assert(pointkist.Nunits() == Ntest + Nadd);
  return 1;
}

void TreeKist::_AddPoint(){
  
  /*/ Test lines
   pointlist->current = current->points;
   for(int k = 0; k< current->npoints ; ++k,MoveDownList(pointlist)){
   assert(inbox(pointlist->current->x
   ,current->boundary_p1,current->boundary_p2));
   }
   ///////////////////////////////////*/
  
	unsigned long i;
	if(current->npoints == 0){
		return;
	}else if(current->npoints <= Nbucket){
		/**** test lines *******************
     pointlist->current = current->points;
     for(i=0;i<current->npoints;i++){
     pointlist->current->leaf = current;
     assert(inbox(pointlist->current->x
     ,pointlist->current->leaf->boundary_p1
     ,pointlist->current->leaf->boundary_p2));
     MoveDownList(pointlist);
     }
     / ****************************************/
		return;
	}
  
	if(atLeaf()){
		//Branch branch1,branch2;
		Branch *current_branch;
		unsigned long dimension,cut;
		PosType xcut;
		Point *oldfirstpoint,*newfirstpoint;
		PosType *x = new PosType[current->npoints];
    
		Branch* branch1 = new Branch(NULL,0,current->boundary_p1,current->boundary_p2
                                 ,current->center,current->level+1);
		Branch* branch2 = new Branch(NULL,0,current->boundary_p1,current->boundary_p2
                                 ,current->center,current->level+1);
    
    
		// initialize boundaries to old boundaries
		/*for(i=0;i<2;++i){
     branch1->boundary_p1[i] = current->boundary_p1[i];
     branch1->boundary_p2[i] = current->boundary_p2[i];
     
     branch2->boundary_p1[i] = current->boundary_p1[i];
     branch2->boundary_p2[i] = current->boundary_p2[i];
     }*/
    
		/*/ Test lines
     pointlist->current = current->points;
     for(int k = 0; k< current->npoints ; ++k,MoveDownList(pointlist)){
     assert(inbox(pointlist->current->x
     ,current->boundary_p1,current->boundary_p2));
     }
     pointlist->current = current->points;
     ///////////////////////////////////*/
    
		// set dimension to cut box
		dimension=(current->level % 2);
    
		// reorder points
		//pointlist->current = current->points;
    pointkist.SetCurrentIt(current->pointit);
		for(i=0;i<current->npoints;++i){
			x[i] = pointkist.getCurrent()->x[dimension];
			//MoveDownList(pointlist);
      --pointkist;
		}
		/*if(x[0] == x[1]){
     dimension = !dimension;
     pointlist->current = current->points;
     for(i=0;i<current->npoints;++i){
     x[i] = pointlist->current->x[dimension];
     MoveDownList(pointlist);
     }
     assert(current->npoints == 2);
     }*/
    
		/*/ Test lines
     pointlist->current = current->points;
     for(int k = 0; k< current->npoints ; ++k,MoveDownList(pointlist)){
     assert(inbox(pointlist->current->x
     ,current->boundary_p1,current->boundary_p2));
     }
     pointlist->current = current->points;
     ///////////////////////////////////*/
    
		oldfirstpoint = current->points;
		if(current->npoints == 2){
			//pointlist->current = current->points;
      pointkist.SetCurrentIt(current->pointit);
			if(pointkist.getCurrent()->x[dimension] > pointkist.getCurrent()->next->x[dimension]){
				//pointlist->current = current->points;
        pointkist.SetCurrentIt(current->pointit);
				bool attop = pointkist.AtTop();
				Point *point = pointkist.TakeOutCurrent();
				if(!attop) pointkist.Down();
				pointkist.InsertAfterCurrent(point);
				current->points = point->prev;
				PosType tmp = x[1];
				x[1] = x[0];
				x[0] = tmp;
			}
		}else{
			ERROR_MESSAGE();
			std::cout << "This is prone to errors and should never happen! npoints in this branch = " << current->npoints << std::endl;
			throw std::runtime_error("Not the right number of points in a leaf");
			//current->points = sortKist(current->npoints,x,pointlist,current->points);
		}
		newfirstpoint = current->points;
		/*/ Test lines
     pointlist->current = current->points;
     for(int k = 0; k< current->npoints ; ++k,MoveDownList(pointlist)){
     assert(inbox(pointlist->current->x
     ,current->boundary_p1,current->boundary_p2));
     }
     pointlist->current = current->points;
     ///////////////////////////////////*/
    
		cut = current->npoints/2;
    
		// check that median split in this dimension will split particles
		if(median_cut && x[cut] == x[cut-1]){
			// change dimension
      
			dimension=!dimension;
      
			//pointlist->current = current->points;
      pointkist.SetCurrentIt(current->pointit);
			for(i=0;i<current->npoints;i++){
				x[i]=pointkist.getCurrent()->x[dimension];
				pointkist.Down();
			}
      
			if(current->npoints == 2){
        pointkist.SetCurrentIt(current->pointit);
				if(pointkist.getCurrent()->x[dimension] > pointkist.getCurrent()->next->x[dimension]){
          pointkist.SetCurrentIt(current->pointit);
					bool attop = pointkist.AtTop();
					Point *point = pointkist.TakeOutCurrent();
					if(!attop) pointkist.Down();
					pointkist.InsertAfterCurrent(point);
					current->points = point->prev;
					PosType tmp = x[1];
					x[1] = x[0];
					x[0] = tmp;
				}
        if(x[0] == x[1]){
          pointkist.SetCurrentIt(current->pointit);
          for(i=0;i<current->npoints;i++){
            std::cout << std::scientific << pointkist.getCurrent()->x[0] << "  "
            << pointkist.getCurrent()->x[1] << " " << pointkist.getCurrent()->id << "      "
            << pointkist.getCurrent()->image->x[0] << "  " << pointkist.getCurrent()->image->x[1] << " "
            << pointkist.getCurrent()->image->id << std::endl;
            std::cout << "invmag " << pointkist.getCurrent()->invmag << " gridsize " << pointkist.getCurrent()->gridsize << std::endl;
            pointkist.Down();
          }
          std::cout << "top bounderies  " << top->boundary_p1[0] << " " << top->boundary_p1[1] << "      " << top->boundary_p2[0] << " " << top->boundary_p2[1] << std::endl;
          throw std::runtime_error("Points in grid are the same");
        }
			}else{
				ERROR_MESSAGE();
				std::cout << "This is prone to errors and this should never happen!" << std::endl;
				exit(1);
				//current->points = sortList(current->npoints,x,pointlist,current->points);
			}
			newfirstpoint=current->points;
		}
    
		assert(x[0] != x[1]);
		if(median_cut){
			cut=current->npoints/2;
      
			branch1->boundary_p2[dimension]=(x[cut]+x[cut-1])/2;
			branch2->boundary_p1[dimension]=(x[cut]+x[cut-1])/2;
      
		}else{
			xcut=(current->boundary_p1[dimension]+current->boundary_p2[dimension])/2;
			branch1->boundary_p2[dimension]=xcut;
			branch2->boundary_p1[dimension]=xcut;
      
			locateD(x-1,current->npoints,xcut,&cut);
		}
		delete[] x;
		/*/  Test lines
     pointlist->current = current->points;
     for(int k = 0; k< current->npoints ; ++k,MoveDownList(pointlist)){
     assert(inbox(pointlist->current->x
     ,current->boundary_p1,current->boundary_p2));
     }
     pointlist->current = current->points;
     ///////////////////////////////////*/
    
		assert(cut <= current->npoints);
    
		// set point numbers and pointers to points
		branch1->npoints = cut;
		if(branch1->npoints > 0) branch1->points = current->points;
		else branch1->points = NULL;
    
		branch2->npoints = current->npoints - cut;
    pointkist.SetCurrentIt(current->pointit);

		if(branch2->npoints > 0){
      pointkist.JumpDown(cut);
			//branch2->points = pointlist->current;
      branch2->pointit = pointkist.getCurrentIt();
		}else{
			branch2->pointit = NULL;
		}
		/*/ Test lines
     pointlist->current = current->points;
     for(int k = 0; k< current->npoints ; ++k,MoveDownList(pointlist)){
     assert(inbox(pointlist->current->x
     ,current->boundary_p1,current->boundary_p2));
     }
     pointlist->current = current->points;
     ///////////////////////////////////*/
    
		// use geometric center
		branch1->center[0] = (branch1->boundary_p1[0] + branch1->boundary_p2[0])/2;
		branch1->center[1] = (branch1->boundary_p1[1] + branch1->boundary_p2[1])/2;
    
		branch2->center[0] = (branch2->boundary_p1[0] + branch2->boundary_p2[0])/2;
		branch2->center[1] = (branch2->boundary_p1[1] + branch2->boundary_p2[1])/2;
    
		/*/ Test lines
     if(branch1->npoints > 0 && !inbox(branch1->points->x,branch1->boundary_p1,branch1->boundary_p2)){
     printBranch(current);
     printBranch(&branch1);
     printBranch(&branch2);
     pointlist->current = current->points;
     for(i=0;i<current->npoints;i++){
     std::cout << pointlist->current->x[0] << "  " << pointlist->current->x[1] << std::endl;
     MoveDownList(pointlist);
     }
     ERROR_MESSAGE();
     exit(0);
     }
     if(branch2->npoints > 0 && !inbox(branch2->points->x,branch2->boundary_p1,branch2->boundary_p2)){
     printBranch(current);
     printBranch(&branch1);
     printBranch(&branch2);
     pointlist->current = current->points;
     for(i=0;i<current->npoints;i++){
     std::cout << pointlist->current->x[0] << "  " << pointlist->current->x[1] << std::endl;
     MoveDownList(pointlist);
     }
     ERROR_MESSAGE();
     exit(0);
     }
     ////////////////////////////////////////////////////*/
		/* centers of mass *
     
     for(i=0;i<2;++i) branch1->center[i]=0;
     pointlist->current=branch1->points;
     for(i=0;i<cut; ++i){
     branch1->center[0]+=pointlist->current->x[0]/branch1->npoints;
     branch1->center[1]+=pointlist->current->x[1]/branch1->npoints;
     MoveDownList(pointlist);
     }
     
     for(i=0;i<2;++i) branch2->center[i]=0;
     pointlist->current=branch2->points;
     for(i=cut;i<current->npoints; ++i){
     branch2->center[0]+=pointlist->current->x[0]/branch2->npoints;
     branch2->center[1]+=pointlist->current->x[1]/branch2->npoints;
     MoveDownList(pointlist);
     }
     */
    
    assert(atLeaf());
		attachChildrenToCurrent(branch1,branch2);
    
		// TODO: take these out when problem is fixed
		assert(current->child1->child1 == NULL);
		assert(current->child1->child2 == NULL);
		assert(current->child2->child1 == NULL);
		assert(current->child2->child2 == NULL);
    
		assert( (current->child1->npoints + current->child2->npoints)
           == current->npoints );
		assert(current->child1->npoints == 1 && current->child2->npoints == 1);
    
		/******* test lines ***********
     pointlist->current = current->points;
     for(int n=0;n<current->npoints;++n){
     assert(pointlist->current->leaf != current);
     assert(pointlist->current->leaf == current->child1 || pointlist->current->leaf == current->child2);
     MoveDownList(pointlist);
     }
     / *******************************************/
    
    
		// If needed, change the particle pointer in all parent cells
		current_branch = current;
		moveUp();
		while(current->points == oldfirstpoint){
			current->points = newfirstpoint;
			if(!moveUp()) break;
		}
		current = current_branch;
	}
  
	// recursively descent to children
	if( current->child1 != NULL ){
		moveToChild(1);
		_AddPoint();
		moveUp();
	}
  
	if( current->child2 != NULL ){
		moveToChild(2);
		_AddPoint();
		moveUp();
	}
}



/** \ingroup ImageFindingL2
 *
 *  Frees all branches of the tree below the current branch in i_tree
 * if that branch is square and i_tree->current->refined == true.
 * If either of these are not true nothing happens.
 *
 * On exit: The i_tree->current is back to the original current.  If it is
 *          square it will have no children and contain one point.  The source
 *          points and branches are also removed.
 */

unsigned long FreeBranchesBelow(TreeHndl i_tree,TreeHndl s_tree,Kist<Point> * trashkist){
  
	if(!i_tree->CurrentIsSquareBranch()) return 0;
	if(i_tree->atLeaf()) return 0;
	if(i_tree->current->refined == false) return 0;
  
	assert( i_tree !=NULL);
	assert( s_tree !=NULL);
  
	Branch *branch,*headbranch;
	Point *point;
	unsigned long Ntmp,NtoRemove,i,count = 0,count2 = 0,count1;
	PosType center[2];
  
	//_freeBranches_iter(s_tree);  // s_tree will no longer be valid on exit.  This is to make sure it isn't used later without a rebuild.
  
	headbranch = i_tree->current;
	//i_tree->TreeWalkStep(true);
  
	while( (headbranch->child1 != NULL) || (headbranch->child2 != NULL) ){
    
		assert(boxinbox(i_tree->current,headbranch));
		if(i_tree->atLeaf()){
			//assert(i_tree->current->points->image->leaf);
			//s_tree->current = i_tree->current->points->image->leaf;  // set s_tree to source of current image cell
      
			/***************** test line  **************************/
			assert(i_tree->current->points->next || i_tree->current->points->prev);
      
			/***************** test line  **************************/
			branch = i_tree->current->prev;
			i = branch->npoints;
      
			if(i_tree->current != headbranch) i_tree->RemoveLeafFromTree(&Ntmp);
      
			/***************** test line  **************************/
			assert(i_tree->current == branch);
			assert(i == i_tree->current->npoints);
      
			/***************** test line  **************************/
			assert(i_tree->current->points->next || i_tree->current->points->prev);
      
			// in a square leaf cell take out extra points that have come up from below
      
			if(i_tree->atLeaf() && i_tree->current->refined){
        /***************** test line  **************************/
				assert(i_tree->current->points->next || i_tree->current->points->prev);
        
				//std::printf("  collecting points from removed leaves\n");
				assert(i_tree->current->points);
				i_tree->pointlist->current = i_tree->current->points;
				NtoRemove = i_tree->current->npoints;
				assert(NtoRemove == 9);
				center[0] = (i_tree->current->boundary_p1[0] + i_tree->current->boundary_p2[0])/2;
				center[1] = (i_tree->current->boundary_p1[1] + i_tree->current->boundary_p2[1])/2;
        
				for(i=0,count1=0,count2=0;i<NtoRemove;++i,MoveDownList(i_tree->pointlist)){
					// find central point and remove others
          
					if( (pow(center[0]-i_tree->pointlist->current->x[0],2)
               + pow(center[1]-i_tree->pointlist->current->x[1],2) )
             < pow(i_tree->pointlist->current->gridsize/2,2) ){
            
						++count1;
						// keep this central point
						i_tree->pointlist->current->gridsize *= 3;
						i_tree->pointlist->current->image->gridsize = i_tree->pointlist->current->gridsize;
						i_tree->current->points = i_tree->pointlist->current;
            
						/***************** test line  **************************/
						assert(i_tree->current->points->next || i_tree->current->points->prev);
						assert(i_tree->pointlist->current->leaf == i_tree->current);
					}else{
            
						++count;  // count of total number of points removed
						++count2;
            
						// reduce the number of particles in all parent cells
            
						/* First take points out of source plane
						 *   This is tricky because they are not ordered into
						 *   square blocks with 9 points in each.
						 */
            
						// Take point out of the source plane
            
						assert(i_tree->pointlist->current->image);
						point = i_tree->pointlist->current->image;
						assert(point->leaf);
						s_tree->current = point->leaf;
						//if(s_tree->current->npoints < 2) RemoveLeafFromTree(s_tree,&Ntmp);
						while(!s_tree->atTop()){
							--(s_tree->current->npoints);
							if(s_tree->current->npoints > 0 && s_tree->current->points == point)
								s_tree->current->points = point->next;
              
							if(s_tree->current->npoints == 0) s_tree->current->points = NULL;
              
							if(s_tree->current->npoints == 0 && s_tree->current->prev->npoints == 1){
								// only remove empty leaves if it will make its parent a leaf
								assert(s_tree->atLeaf());
								s_tree->RemoveLeafFromTree(&Ntmp);
								s_tree->TreeWalkStep(true);  // Go to other child.
								s_tree->RemoveLeafFromTree(&Ntmp);
							}else{
								s_tree->moveUp();
							}
						}
						assert(boxinbox(i_tree->current,headbranch));
            
						// Do it for top
						--(s_tree->top->npoints);
						if(s_tree->top->npoints > 0 && s_tree->top->points == point) s_tree->top->points = point->next;
            
						s_tree->pointlist->current = point;
						TakeOutCurrent(s_tree->pointlist);
						point->leaf = NULL;  // set leaf to NULL to indicate that point is no longer in tree
						if(point->head) trashkist->InsertAfterCurrent(point);  // collect heads for later trash collection
            
						assert(boxinbox(i_tree->current,headbranch));
            
						// take points out of image plane
						branch = i_tree->current;
						do{
							assert(i_tree->current->npoints);
							--(i_tree->current->npoints);
							if(i_tree->current->points == i_tree->pointlist->current) i_tree->current->points = NULL;
              
						}while(i_tree->moveUp());
						i_tree->current = branch;
						assert(boxinbox(i_tree->current,headbranch));
            
						if(i_tree->pointlist->current == i_tree->current->points) i_tree->current->points = i_tree->pointlist->current->next;
						point = TakeOutCurrent(i_tree->pointlist);
						point->leaf = NULL;
						// If point is a head of a memory block add it to trashlist for eventual trash collection
						if(point->head){
							assert(point->head == 8);
							trashkist->InsertAfterCurrent(point);
						}
					}
					assert(boxinbox(i_tree->current,headbranch));
          
				}  // loop through points in leaf
        
				// reassign first point in branches above the current
				branch = i_tree->current;
				do{
					assert(i_tree->current->npoints);
					if(i_tree->current->points == NULL) i_tree->current->points = branch->points;
				}while(i_tree->moveUp());
				i_tree->current = branch;
        
				assert(boxinbox(i_tree->current,headbranch));
        
				assert(count1 == 1);
				assert(count2 == 8);
				assert(i_tree->current->npoints == 1);
			}  // if current was leaf that was refined
      
		} // at tree leaf
		assert(i_tree->current->points);
		assert(boxinbox(i_tree->current,headbranch));
    
		if( !(i_tree->atLeaf()) ) i_tree->TreeWalkStep(true);
	}  // while entry current is not a leaf
  
	//assert(CurrentIsSquareTree(i_tree));
	assert(i_tree->current->npoints == 1);
	assert(i_tree->atLeaf());
	assert(i_tree->current == headbranch);
  
	i_tree->current->refined = false;
  
	// Free the memory for the points that have been removed.
	//CollectTrash(trashkist,false);
  
	//assert(trashkist->Nunits() == 0);
	//if(count) std::printf("FreeBranchesBelow() freed %li points and moved up %li points\n",count,count2);
  return count;
}

/** \ingroup LowLevel
 * Removes current from a tree if it is a leaf.
 *   Will not remove root of tree.
 *
 *  on output: Current is left at the father of the leaf that was removed.
 *             All the points in the leaf that was removed are in its father
 *             so the father might be a leaf without Nbucket points.
 *             The ->leaf pointer of these points are reassigned to the father.
 *
 *  returns: Pointer to first in list of points that were reassigned.
 *           *Npoints = number of points reassigned.
 */
Point * TreeKist::RemoveLeafFromTree(unsigned long *Npoints){
  
	Branch *branch;
	Point *point = NULL;
	unsigned long i;
  
	if(atTop() || !(atLeaf()) ) return NULL;
  
	branch = current;
	moveUp();
  
	//if(current->number == 2125226) branchaddress = NULL;
  
	if(branch == branch->prev->child1){
		branch->prev->child1 = NULL;
	}
  
	if(branch == branch->prev->child2){
		if(branch->prev->child1 != NULL) branch->prev->child1->brother = branch->prev->brother;
		branch->prev->child2 = NULL;
	}
  
	// leaves of points in the father
	if(branch->npoints >0){
    pointkist.SetCurrentIt(branch->pointit);

		assert(inbox(pointkist.getCurrent()->x,branch->boundary_p1,branch->boundary_p2));
    
		for(i=0;i<branch->npoints;++i,pointkist.Down()){
			pointkist.getCurrent()->leaf = branch->prev;
      
			assert(inbox(pointkist.getCurrent()->x,branch->boundary_p1,branch->boundary_p2));
			assert(inbox(pointkist.getCurrent()->x,branch->prev->boundary_p1,branch->prev->boundary_p2));
		}
    
		assert(boxinbox(branch,branch->prev));
    
		point = branch->points;
	}
	*Npoints = branch->npoints;
	delete branch;
	--Nbranches;
  
	return point;
}



