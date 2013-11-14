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
		,double boundary_p1[2]   /// bottom left hand corner of root
		,double boundary_p2[2]   /// upper right hand corner of root
		,double center[2]        /// center of root (this could be the center of mass)
		,int my_Nbucket             /// maximum number of points allowed in a leaf
		){
	construct_root(xp,npoints,boundary_p1,boundary_p2,center,my_Nbucket);
}
/// Basic construction of root with all particles in it but no children
void TreeKist::construct_root(
		Point *xp   /// array of points to be added to the tree
		,unsigned long npoints   /// number of points
		,double boundary_p1[2]   /// bottom left hand corner of root
		,double boundary_p2[2]   /// upper right hand corner of root
		,double center[2]        /// center of root (this could be the center of mass)
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
 * Returns the points of current.  Exported.
 * Pre: !offEnd(tree)
 ************************************************************************/
void TreeKist::getCurrent(Point *points,unsigned long *npoints){

    if( offEnd() ){
    	ERROR_MESSAGE();
    	std::cout << "Tree Error: calling getCurrent() when current is off end" << std::endl;
    	exit(1);
    }

    *npoints=current->npoints;
    points=current->points;

    return;
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
			  ,double boundary_p1[2],double boundary_p2[2]
			  ,double center[2],int child){
    
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
    	pointlist->current = branch->points;
    	for(unsigned long i=0;i<branch->npoints;++i){
    		pointlist->current->leaf = branch;
    		MoveDownList(pointlist);
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
    pointlist->current = current->points;
    for(i=0;i<current->npoints;++i){
      std::cout << pointlist->current->id << " " << pointlist->current->x[0] << " " << pointlist->current->x[1] << std::endl;
      MoveDownList(pointlist);
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
void TreeKist::PointsInCurrent(unsigned long *ids,double **x){
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

