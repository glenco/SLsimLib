/*
 * tree_maintenance.c
 *
 *  Created on: Sep 29, 2011
 *      Author: bmetcalf
 *
 *      This file contains routines for building, adding to, and removing grid points from
 *      the tree.
 */

#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision

#include "slsimlib.h"

/************************** test routine *****************************/
bool tree_count_test(TreeHndl tree){
	/*static int init=0;
	static PosType p1[2],p2[2];

	if(init==0){
		p1[0] = tree->top->boundary_p1[0];
		p1[1] = tree->top->boundary_p1[1];
		p2[0] = tree->top->boundary_p2[0];
		p2[1] = tree->top->boundary_p2[1];
	}
	assert(p1[0] == tree->top->boundary_p1[0]);
	assert(p1[1] == tree->top->boundary_p1[1]);
	assert(p2[0] == tree->top->boundary_p2[0]);
	assert(p2[1] == tree->top->boundary_p2[1]);
	 */

/*	do{
		tree->pointlist->current = tree->current->points;
		for(i=0;i<tree->current->npoints;++i,MoveDownList(tree->pointlist))
			assert(inbox(tree->pointlist->current->x,tree->current->boundary_p1,tree->current->boundary_p2));

		++nbranches;
		assert(nbranches <= tree->Nbranches);

	}while(tree->TreeWalkStep(true) && tree->current != branch->brother);

	tree->current = branch;
*/
	return true;
}

bool testLeafs(TreeHndl tree){
	Branch *leaf;
	Point *point = tree->pointlist->current;
	MoveToTopList(tree->pointlist);
	for(unsigned long i=0;i<tree->pointlist->Npoints;MoveDownList(tree->pointlist),++i){
		leaf = tree->pointlist->current->leaf;
		if(leaf->child1 != NULL || leaf->child2 != NULL){
			std::cout << "a point " << tree->pointlist->current->id << "'s leaf is not a leaf!"
					<< " leaf " << leaf << std::endl;
			tree->pointlist->current = point;
			return false;
		}
		if(!inbox(tree->pointlist->current->x,leaf->boundary_p1,leaf->boundary_p2)){
			std::cout << "point " << tree->pointlist->current->id << " is not in its leaf!"
					<< " leaf " << leaf << std::endl;
			tree->pointlist->current = point;
			return false;
		}
		assert(tree->pointlist->current->prev != NULL || tree->pointlist->current->next != NULL );
		if(!(tree->pointlist->current->image->prev != NULL || tree->pointlist->current->image->next != NULL) ){
			std::cout << "point " << tree->pointlist->current->id << " image of point is not in image tree list!"
					<< " leaf " << leaf << std::endl;
			std::cout << "image " << tree->pointlist->current->image->id
					<< " leaf " << tree->pointlist->current->image->leaf << std::endl;
			tree->pointlist->current = point;
			return false;
		}
		/*if(leaf->npoints != 1){
			std::cout << "a point's leaf has more than one point in it!" << std::endl;
			tree->pointlist->current = point;
			return false;
		}*/
	}
	tree->pointlist->current = point;

	return true;
}
/******************************************************************/

//static int median_cut=1;

/**
 * \brief Build a complete tree from a list of points.

 <p>
* median_cut determines how the cells are subdivided 
*    if ==0  equal volume cuts, Warning this option causes an error
*    if ==1  pseudo-median point cuts, never cuts through a point, but near the median
<\p>
 
 */
//TreeHndl BuildTree(Point *xp,unsigned long Npoints,short my_median_cut){
TreeStruct::TreeStruct(Point *xp,unsigned long Npoints,short my_median_cut,PosType buffer){
  //TreeHndl tree;
  unsigned long i;
  PosType p1[2],p2[2],center[2];
  const int my_Nbucket = 1;   // must be =1 if each leaf is to coincide with each cell
/*
  if( (Npoints & (Npoints-1)) != 0){
	  ERROR_MESSAGE();
	  std::printf("ERROR: BuildTree, Npoints is not a power of 2\n");
	  exit(1);
  }
*/
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
  //tree=NewTree(xp,Npoints,p1,p2,center,Nbucket);
  construct_root(xp,Npoints,p1,p2,center,my_Nbucket);

  median_cut = my_median_cut;
 /* build the tree */
  TreeStruct::iterator current(top);
  _BuildTree(current);

  /************************** test routine *****************************/
  //moveTop();
  //tree_count_test(tree);

  //return tree;
}
/** \ingroup  ImageFindingL2
 * \brief Fill a tree with points.  The previous tree structure will be destroyed.  Used for refilling.
 */
void TreeStruct::FillTree(Point *xp,unsigned long Npoints){
  unsigned long i;

  top->points=xp;
  top->npoints=Npoints;
  // link point array into point list
  EmptyList(pointlist);
  for(i=0;i<Npoints;++i){
    InsertPointAfterCurrent(pointlist,&xp[i]);
    MoveDownList(pointlist);
  }

  MoveToTopList(pointlist);
 /* build the tree */
  // make sure there are no branches in the tree
  _freeBranches_iter();

  TreeStruct::iterator current(top);
  _BuildTree(current);

  return ;
}

/** \ingroup  ImageFindingL2
 * \brief Rebuilds the tree from the points that are already in the tree->pointlist
 *
 * This is not the best function because it copies all the points
 */
void TreeStruct::RebuildTreeFromList(){
 /* Builds or rebuilds a tree that already has a tree->pointlist.
  */

	if(pointlist->Npoints < 1) return;

	unsigned long Npoints = pointlist->Npoints,i;

	// Make new array of points
	Point *points = NewPointArray(Npoints);

	MoveToTopList(pointlist);
	for(i=0;i<Npoints;++i,MoveDownList(pointlist)){
		//tmp = points[i].x;
		// PointCopy() copies the x pointer which is later freed in
		// emptyTree() so the information needs to be copied

		PointCopy(&(points[i]),pointlist->current);
		//points[i].x = tmp;
		//points[i].x[0] = pointlist->current->x[0];
		//points[i].x[1] = pointlist->current->x[1];

		/* Below is for trash collection in PruneTree().
		 * When PruneTree has been used not all of the
		 * heads in the point arrays are guaranteed to be
		 * in the pointlist.
		 */
		pointlist->current->leaf = NULL;
	}
	//std::printf(" %e %e\n",points[0].x[0],points[0].x[1]);
	// emptry the tree and free all former points
	emptyTree();

	assert(Nbranches == 1);
	assert(top->npoints == 0);
	assert(pointlist->Npoints == 0);

	FillTree(points,Npoints);

	return;
}

/** \brief Spawn a subtree with current as its top
 *
 *  The new tree contains all of the tree below the current.
 *  Warning:: Adding points to the new tree will not update the 
 *  parent tree so it can become dangerously out of sync.
 */
TreeStruct * TreeStruct::spawn(TreeStruct::iterator &current){
  throw std::runtime_error("This is untested and could cause significant problems");
  
  TreeStruct *newTree = new TreeStruct;

  newTree->Nbucket = Nbucket;
  newTree->top = *current;
  newTree->pointlist->top = (*current)->points;
  newTree->pointlist->Npoints = (*current)->npoints;
  pointlist->current = (*current)->points;
  for(size_t i=0 ; i < (*current)->npoints-1 ; ++i) MoveDownList(pointlist);
  newTree->pointlist->bottom = pointlist->current;
  
  // count the number of branches below
  do{
    newTree->Nbranches++;
  }while(current.TreeWalkStep(true) && (*current) != newTree->top->brother);
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
short TreeStruct::emptyTree(){
  Point **heads;
  unsigned long i,j,count;

  heads = (Point **) malloc(pointlist->Npoints*sizeof(Point*));  // maximum number of pointers that could be needed

  _freeBranches_iter();

  assert(Nbranches == 1);

  MoveToTopList(pointlist);
  for(i=0,j=0,count=0;i<pointlist->Npoints;++i){
	  if(pointlist->current->head){
		  heads[j] = pointlist->current;
		  ++j;
		  count += pointlist->current->head;
	  }
	  pointlist->current->leaf = NULL;  // This is for future trash collection if PruneTree() has been used.
	  MoveDownList(pointlist);
  }

  //std::printf("freed %i arrays out of %i points in array, %i freed\n",j,i,count);
  for(i=0;i<j;++i) FreePointArray(heads[i]);
  top->npoints = 0;

  //printTree;

  pointlist->Npoints=0;
  pointlist->top = pointlist->bottom = pointlist->current=NULL;

  free(heads);
  return 1;
}

/** \ingroup LowLevel
* \brief Recursively free branches
*/
void TreeStruct::_freeBranches(TreeStruct::iterator &current,short child){
	Branch *branch;

	/*printBranch(current);*/

	if((*current)->child1 != NULL){
    current.down(1);
		_freeBranches(current,1);
	}

    if((*current)->child2 != NULL){
      current.down(2);
      _freeBranches(current,2);
    }

    if( ((*current)->child1 == NULL)*((*current)->child2 == NULL) ){

    	if(current.atTop()){
        return;
    	}

    	branch = *current;
      current.up();
      delete branch;

      if(child==1) (*current)->child1=NULL;
    	if(child==2) (*current)->child2=NULL;
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
void TreeStruct::_freeBranches_iter(){
	Branch *branch;

	//moveTop();
	/*printBranch(current);*/

  TreeStruct::iterator current(top);
  
	while(Nbranches > 1){
		if((*current)->child1 != NULL){
			current.down(1);
		}else if((*current)->child2 != NULL){
			current.down(2);
		}else{
			branch = *current;
			if((*current)->brother == (*current)->prev->brother){
        current.up();
				assert((*current)->child1 == NULL);
				(*current)->child2 = NULL;
			}else{
				current = (*current)->brother;
				(*current)->prev->child1 = NULL;
			}
			delete branch;
			--Nbranches;
		}
	}

	// re-assign the leaf pointers in the particles to the root
	 MoveToTopList(pointlist);
	 do{ pointlist->current->leaf = top; }while(MoveDownList(pointlist));

    return;
}


/** \ingroup LowLevel
* \brief Recursively build tree from points in its linked list.
*/
void TreeStruct::_BuildTree(TreeStruct::iterator &current){
  /* pointlist must be both a linked list and an array of points in the */
  /* same order as the linked list */
  unsigned long i,cut,dimension;
  Branch *cbranch;//,branch1,branch2;
  PosType xcut;

  cbranch = *current; /* pointer to current branch */

    /* leaf case */
  if(cbranch->npoints <= Nbucket){
	  (*current)->points->leaf = *current;
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
  pointlist->current = (*current)->points;
  for(i=0;i<cbranch->npoints;++i){
    x[i]=pointlist->current->x[dimension];
    /*points[i]=pointlist->current->id;*/
    MoveDownList(pointlist);
  }
  
  /*PrintList(pointlist);*/

  //double_sort(cbranch->npoints,x-1,points-1);
  //double_sort_points(cbranch->npoints,x-1,current->points);

  /* copy information back to points in new order */
/*   pointlist->current=current->points; */
/*   for(i=0;i<cbranch->npoints;++i){ */
    /*     pointlist->current->x[0]=xp[points[i]].x[0]; */
    /*     pointlist->current->x[1]=xp[points[i]].x[1]; */
/*     pointlist->current->id=points[i]; */
/*     MoveDownList(pointlist); */
/*   } */

  if(median_cut){
    //double_sort_points(cbranch->npoints,x-1,current->points);
	  Utilities::quicksortPoints((*current)->points,x,cbranch->npoints);
        
    if(x[0] == x[cbranch->npoints-1]){
      dimension = !dimension;
      /* reorder points */
      pointlist->current = (*current)->points;
      for(i=0;i<cbranch->npoints;++i){
        x[i]=pointlist->current->x[dimension];
        MoveDownList(pointlist);
      }

      Utilities::quicksortPoints((*current)->points,x,cbranch->npoints);
    }

	  cut=cbranch->npoints/2;
    size_t ii = cut-1;
    while(x[cut] == x[ii] && ii > 0 ) --ii;  // find closest unique value
    while(x[cut] == x[ii] && ii < cbranch->npoints - 1 ) ++ii;  // try at higher index
   
    if(ii < cut){
      cut = ii + 1;
    }else{
      cut = ii;
      --ii;
    }
    branch1->boundary_p2[dimension]=(x[cut]+x[ii])/2;
    branch2->boundary_p1[dimension]=(x[cut]+x[ii])/2;
    
  }else{

	  xcut=(cbranch->boundary_p1[dimension]+cbranch->boundary_p2[dimension])/2;
    branch1->boundary_p2[dimension]=xcut;
    branch2->boundary_p1[dimension]=xcut;

	  Utilities::quickPartitionPoints(xcut,&cut,(*current)->points,x,cbranch->npoints);

  }

  /*/ Test lines /////////////////////////////////////////////////
  pointlist->current = current->points;
  for(int k = 0; k< current->npoints-1 ; ++k,MoveDownList(pointlist)){
    assert(pointlist->current->x[dimension] <= pointlist->current->next->x[dimension]);
    assert(inbox(pointlist->current->x,current->boundary_p1,current->boundary_p2));
  }
  //////////////////////////////////////////////////////////////*/

  /* set point numbers and pointers to points */
  branch1->npoints=cut;
  assert((*current)->points->next || (*current)->points->prev);
  branch1->points = (*current)->points;

  /*/ Test lines /////////////////////////////////////////////////
  pointlist->current = branch1->points;
  for(int k = 0; k< branch1->npoints ; ++k,MoveDownList(pointlist)){
    assert(inbox(pointlist->current->x,branch1->boundary_p1,branch1->boundary_p2));
  }
   //////////////////////////////////////////////////////////////*/

  
  branch2->npoints=cbranch->npoints - cut;
  pointlist->current = (*current)->points;
  JumpDownList(pointlist,cut);
  branch2->points=pointlist->current;

  /*/ Test lines /////////////////////////////////////////////////
  pointlist->current = branch2->points;
  for(int k = 0; k< branch2->npoints ; ++k,MoveDownList(pointlist)){
    assert(inbox(pointlist->current->x,branch2->boundary_p1,branch2->boundary_p2));
  }
  //////////////////////////////////////////////////////////////*/
  
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

  attachChildrenToCurrent(*current,branch1,branch2);

  if( branch1->npoints > 0 ){
    current.down(1);
     _BuildTree(current);
    current.up();
  }

  if(branch2->npoints > 0 ){
    current.down(2);
	  _BuildTree(current);
    current.up();
  }

 /*std::printf("reached end of _BuildTree level=%i\n",current->level);*/
  return;
}

/** \ingroup  ImageFindingL2
 *  \brief Expands tree by adding points
*/
int TreeStruct::AddPointsToTree(Point *xpoint,unsigned long Nadd){
  unsigned long j,Ntest;
  //Branch *parent_branch;

  if(Nadd==0) return 1;

  Ntest = pointlist->Npoints;
  
  TreeStruct::iterator current(top);
  
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

        current.movetop();
     		_FindLeaf(current,xpoint[j].x,1);
     		assert(current.atLeaf());

    		assert(inbox(xpoint[j].x,(*current)->boundary_p1,(*current)->boundary_p2));
    		assert(inbox(xpoint[j].x,(*current)->prev->boundary_p1,(*current)->prev->boundary_p2));
        
        /*/ Test lines /////////////////////////////////////////////////////////
        pointlist->current = current->points;
        for(int k = 0; k < current->npoints - 1 ; ++k,MoveDownList(pointlist)){
          assert( inbox(pointlist->current->x,current->boundary_p1,current->boundary_p2) );
        }
        pointlist->current = current->points;
        /////////////////////////////////////////////////////////////////////*/

    		// insert point into point list
    		//if((*current)->points == NULL){
    		if((*current)->npoints == 1){
    			// case of no previous points in leaf

    			(*current)->points = &xpoint[j];
    			(*current)->points->leaf = *current;
    			// put point into right place in list
          pointlist->current = (*current)->prev->points;
    			if(*current == (*current)->prev->child1){
       				Point *oldpoint = pointlist->current;
       				Branch *tmp = *current;

       				// insert the point
       				if(MoveUpList(pointlist)){
       					InsertPointAfterCurrent(pointlist,(*current)->points);
       				}else{
       					InsertPointBeforeCurrent(pointlist,(*current)->points);
       				}
  					  // InsertPointBeforeCurrent(pointlist,(*current)->points);
         			// re-assign parent first points
            current.up();
              while((*current)->points == oldpoint){
                (*current)->points = &xpoint[j];
                current.up();
              }
              current = tmp;

    			}else{
              JumpDownList(pointlist,(*current)->prev->npoints-2);  // adds point to end of parent branches list
        			InsertPointAfterCurrent(pointlist,(*current)->points);
    			}

    			//pointlist->current = (*current)->points;
          //(*current)->points->leaf = current;

          assert((*current)->points->next || (*current)->points->prev);
          
          /*/ Test lines //////////////////////////////////////////////////////
    			pointlist->current = (*current)->points;
    			for(int k = 0; k< (*current)->npoints ; ++k,MoveDownList(pointlist)){
    				assert(inbox(pointlist->(*current)->x
                         ,(*current)->boundary_p1,(*current)->boundary_p2));
    			}
          pointlist->current = (*current)->points;
    			/////////////////////////////////////////////////////////////////////*/

    		}else{ // case where the leaf needs to be divided
    			assert((*current)->npoints > 1);

    			/*/ Test lines /////////////////////////////////////////////////////////
    			pointlist->current = (*current)->points;
    			for(int k = 0; k< (*current)->npoints - 1 ; ++k,MoveDownList(pointlist)){
    				assert(inbox(pointlist->current->x
                         ,(*current)->boundary_p1,(*current)->boundary_p2));
    			}
          pointlist->current = (*current)->points;
    			/////////////////////////////////////////////////////////////////////*/

    			// case of points already in leaf
          assert((*current)->points->next || (*current)->points->prev);
          // Test lines
          //if(!testLeafs()){ERROR_MESSAGE(); std::cout << "before _addPoint of AddPointsToTree "<< std::endl; exit(1);}

          // adds point to end of branches list, note that (*current)->npoints has already been increased
    			pointlist->current = (*current)->points;
    			JumpDownList(pointlist,(*current)->npoints-2);
    			InsertPointAfterCurrent(pointlist,&xpoint[j]);
          
          assert(inbox(xpoint[j].x,(*current)->boundary_p1,(*current)->boundary_p2));

    			pointlist->current = (*current)->points;
    			xpoint[j].leaf = *current;

          Branch *parent_branch = *current;
          //unsigned long n=pointlist->Npoints;
          _AddPoint(current);
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

    /*/ Test lines /////////////////////////////////////////////////////////
    pointlist->current = (*current)->points;
    for(int k = 0; k < (*current)->npoints ; ++k,MoveDownList(pointlist)){
      assert(inbox(pointlist->current->x
                   ,(*current)->boundary_p1,(*current)->boundary_p2));
    }
    pointlist->current = (*current)->points;
    /////////////////////////////////////////////////////////////////////*/

    while( !inbox(xpoint[j].x,(*current)->boundary_p1,(*current)->boundary_p2) ) current.up();
    _FindLeaf(current,xpoint[j].x,0);
    assert(current.atLeaf());
    xpoint[j].leaf = *current;
    assert(inbox(xpoint[j].x,xpoint[j].leaf->boundary_p1,xpoint[j].leaf->boundary_p2));
  }

  assert(pointlist->Npoints == Ntest + Nadd);
  return 1;
}

void TreeStruct::_AddPoint(TreeStruct::iterator &current){

		/*/ Test lines /////////////////////////////////////////////////
		pointlist->current = current->points;
		for(int k = 0; k< current->npoints ; ++k,MoveDownList(pointlist)){
			assert(inbox(pointlist->current->x
				,current->boundary_p1,current->boundary_p2));
		}
		//////////////////////////////////////////////////////////////*/

	unsigned long i;
	if((*current)->npoints == 0){
		return;
	}else if((*current)->npoints <= Nbucket){
		/**** test lines *******************
		pointlist->current = (*current)->points;
		for(i=0;i<(*current)->npoints;i++){
			pointlist->current->leaf = current;
			assert(inbox(pointlist->current->x
					,pointlist->current->leaf->boundary_p1
					,pointlist->current->leaf->boundary_p2));
			MoveDownList(pointlist);
		}
		/ ****************************************/
		return;
	}

	if(current.atLeaf()){
		//Branch branch1,branch2;
		Branch *current_branch;
		unsigned long dimension,cut;
		PosType xcut;
		Point *oldfirstpoint,*newfirstpoint;
		PosType *x = new PosType[(*current)->npoints];

		/*// Test lines
		pointlist->current = (*current)->points;
		for(int k = 0; k< (*current)->npoints ; ++k,MoveDownList(pointlist)){
			assert(inbox(pointlist->current->x
                   ,top->boundary_p1,top->boundary_p2));
			assert(inbox(pointlist->current->x
                   ,(*current)->boundary_p1,(*current)->boundary_p2));
		}
		pointlist->current = (*current)->points;
		///////////////////////////////////*/
  
    Branch* branch1 = new Branch(NULL,0,(*current)->boundary_p1,(*current)->boundary_p2
				,(*current)->center,(*current)->level+1);
		Branch* branch2 = new Branch(NULL,0,(*current)->boundary_p1,(*current)->boundary_p2
				,(*current)->center,(*current)->level+1);


		// initialize boundaries to old boundaries
		/*for(i=0;i<2;++i){
			branch1->boundary_p1[i] = (*current)->boundary_p1[i];
			branch1->boundary_p2[i] = (*current)->boundary_p2[i];

			branch2->boundary_p1[i] = (*current)->boundary_p1[i];
			branch2->boundary_p2[i] = (*current)->boundary_p2[i];
		}*/

		/*/ Test lines
		pointlist->current = (*current)->points;
		for(int k = 0; k< (*current)->npoints ; ++k,MoveDownList(pointlist)){
			assert(inbox(pointlist->current->x
				,(*current)->boundary_p1,(*current)->boundary_p2));
		}
		pointlist->current = (*current)->points;
		///////////////////////////////////*/

		// set dimension to cut box
		dimension=((*current)->level % 2);

		// reorder points
		pointlist->current = (*current)->points;
		for(i=0;i<(*current)->npoints;++i){
			x[i] = pointlist->current->x[dimension];
			MoveDownList(pointlist);
		}
		/*if(x[0] == x[1]){
			dimension = !dimension;
			pointlist->current = (*current)->points;
			for(i=0;i<(*current)->npoints;++i){
				x[i] = pointlist->current->x[dimension];
				MoveDownList(pointlist);
			}
			assert((*current)->npoints == 2);
		}*/

		/*/ Test lines
		pointlist->current = (*current)->points;
		for(int k = 0; k< (*current)->npoints ; ++k,MoveDownList(pointlist)){
			assert(inbox(pointlist->current->x
				,(*current)->boundary_p1,(*current)->boundary_p2));
		}
		pointlist->current = (*current)->points;
		///////////////////////////////////*/

		oldfirstpoint = (*current)->points;
		if((*current)->npoints == 2){
			pointlist->current = (*current)->points;
			if(pointlist->current->x[dimension] > pointlist->current->next->x[dimension]){
				pointlist->current = (*current)->points;
				bool attop=AtTopList(pointlist);
				Point *point = TakeOutCurrent(pointlist);
				if(!attop) MoveDownList(pointlist);
				InsertPointAfterCurrent(pointlist,point);
				(*current)->points = point->prev;
				PosType tmp = x[1];
				x[1] = x[0];
				x[0] = tmp;
			}
		}else{
			ERROR_MESSAGE();
			std::cout << "This is prone to errors and should never happen! npoints in this branch = " << (*current)->npoints << std::endl;
			throw std::runtime_error("Not the right number of points in a leaf");
			(*current)->points = sortList((*current)->npoints,x,pointlist,(*current)->points);
		}
		newfirstpoint = (*current)->points;
  
		cut = (*current)->npoints/2;

		// check that median split in this dimension will split particles
		if(median_cut && x[cut] == x[cut-1]){
			// change dimension

			dimension=!dimension;

			pointlist->current = (*current)->points;
			for(i=0;i<(*current)->npoints;i++){
				x[i]=pointlist->current->x[dimension];
				MoveDownList(pointlist);
			}

			if((*current)->npoints == 2){
				pointlist->current = (*current)->points;
				if(pointlist->current->x[dimension] > pointlist->current->next->x[dimension]){
					pointlist->current = (*current)->points;
					bool attop=AtTopList(pointlist);
					Point *point = TakeOutCurrent(pointlist);
					if(!attop) MoveDownList(pointlist);
					InsertPointAfterCurrent(pointlist,point);
					(*current)->points = point->prev;
					PosType tmp = x[1];
					x[1] = x[0];
					x[0] = tmp;
				}
        if(x[0] == x[1]){
          pointlist->current = (*current)->points;
          for(i=0;i<(*current)->npoints;i++){
            std::cout << std::scientific << std::setprecision(15) << pointlist->current->x[0] << "  " << pointlist->current->x[1] << " " << pointlist->current->id << "      "
                      << pointlist->current->image->x[0] << "  " << pointlist->current->image->x[1] << " " << pointlist->current->image->id << std::endl;
            std::cout << "invmag " << pointlist->current->invmag << " gridsize " << pointlist->current->gridsize << std::endl;
            MoveDownList(pointlist);
          }
          std::cout << std::setprecision(15) << "top bounderies  " << top->boundary_p1[0] << " " << top->boundary_p1[1] << "      " << top->boundary_p2[0] << " " << top->boundary_p2[1] << std::endl;
          throw std::runtime_error("Points in grid are the same");
        }
			}else{
				ERROR_MESSAGE();
				std::cout << "This is prone to errors and this should never happen!" << std::endl;
				exit(1);
				(*current)->points = sortList((*current)->npoints,x,pointlist,(*current)->points);
			}
			newfirstpoint=(*current)->points;
		}

		assert(x[0] != x[1]);
		if(median_cut){
			cut=(*current)->npoints/2;

			branch1->boundary_p2[dimension]=(x[cut]+x[cut-1])/2;
			branch2->boundary_p1[dimension]=(x[cut]+x[cut-1])/2;

		}else{
			xcut=((*current)->boundary_p1[dimension]+(*current)->boundary_p2[dimension])/2;
			branch1->boundary_p2[dimension]=xcut;
			branch2->boundary_p1[dimension]=xcut;

			locateD(x-1,(*current)->npoints,xcut,&cut);
		}
    assert(branch1->boundary_p1[0] < branch1->boundary_p2[0]);
    assert(branch1->boundary_p1[1] < branch1->boundary_p2[1]);
    assert(branch2->boundary_p1[0] < branch2->boundary_p2[0]);
    assert(branch2->boundary_p1[1] < branch2->boundary_p2[1]);
    assert((*current)->boundary_p1[0] < (*current)->boundary_p2[0]);
    assert((*current)->boundary_p1[1] < (*current)->boundary_p2[1]);
		delete[] x;
		/*/  Test lines
		pointlist->current = (*current)->points;
		for(int k = 0; k< (*current)->npoints ; ++k,MoveDownList(pointlist)){
			assert(inbox(pointlist->current->x
				,(*current)->boundary_p1,(*current)->boundary_p2));
		}
		pointlist->current = (*current)->points;
		///////////////////////////////////*/

		assert(cut <= (*current)->npoints);

		// set point numbers and pointers to points
		branch1->npoints = cut;
		if(branch1->npoints > 0) branch1->points = (*current)->points;
		else branch1->points = NULL;

		branch2->npoints = (*current)->npoints - cut;
		pointlist->current = (*current)->points;

		if(branch2->npoints > 0){
			JumpDownList(pointlist,cut);
			branch2->points = pointlist->current;
		}else{
			branch2->points = NULL;
		}
		/*/ Test lines
		pointlist->current = (*current)->points;
		for(int k = 0; k< (*current)->npoints ; ++k,MoveDownList(pointlist)){
			assert(inbox(pointlist->current->x
				,(*current)->boundary_p1,(*current)->boundary_p2));
		}
		pointlist->current = (*current)->points;
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
			pointlist->current = (*current)->points;
			for(i=0;i<(*current)->npoints;i++){
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
			pointlist->current = (*current)->points;
			for(i=0;i<(*current)->npoints;i++){
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
	for(i=cut;i<(*current)->npoints; ++i){
		branch2->center[0]+=pointlist->current->x[0]/branch2->npoints;
		branch2->center[1]+=pointlist->current->x[1]/branch2->npoints;
		MoveDownList(pointlist);
	}
	*/

    assert(current.atLeaf());
		attachChildrenToCurrent(*current,branch1,branch2);

		// TODO: take these out when problem is fixed
		assert((*current)->child1->child1 == NULL);
		assert((*current)->child1->child2 == NULL);
		assert((*current)->child2->child1 == NULL);
		assert((*current)->child2->child2 == NULL);
    
    /*/ Test lines /////////////////////////////////////////////////////////
    pointlist->current = (*current)->child1->points;
    for(int k = 0; k< (*current)->child1->npoints - 1 ; ++k,MoveDownList(pointlist)){
      assert(inbox(pointlist->current->x
                   ,(*current)->child1->boundary_p1,(*current)->child1->boundary_p2));
    }
    pointlist->current = (*current)->child2->points;
    for(int k = 0; k< (*current)->child2->npoints - 1 ; ++k,MoveDownList(pointlist)){
      assert(inbox(pointlist->current->x
                   ,(*current)->child2->boundary_p1,(*current)->child2->boundary_p2));
    }
    pointlist->current = (*current)->points;
    /////////////////////////////////////////////////////////////////////*/

		assert( ((*current)->child1->npoints + (*current)->child2->npoints)
				== (*current)->npoints );
		assert((*current)->child1->npoints == 1 && (*current)->child2->npoints == 1);

		/******* test lines ***********
		pointlist->current = (*current)->points;
		for(int n=0;n<(*current)->npoints;++n){
			assert(pointlist->current->leaf != current);
			assert(pointlist->current->leaf == (*current)->child1 || pointlist->(*current)->leaf == (*current)->child2);
			MoveDownList(pointlist);
		}
		/ *******************************************/


		// If needed, change the particle pointer in all parent cells
		current_branch = *current;
    current.up();
		while((*current)->points == oldfirstpoint){
			(*current)->points = newfirstpoint;
			if(!(current.up())) break;
		}
		current = current_branch;
	}

	// recursively descent to children
	if( (*current)->child1 != NULL ){
    current.down(1);
		_AddPoint(current);
    current.up();
	}

	if( (*current)->child2 != NULL ){
    current.down(2);
		_AddPoint(current);
    current.up();
	}
}

/** \ingroup  ImageFinding
 *  \brief THIS DOES NOT WORK YET!!!
 *
 *  Reduces the size of the tree by removing points and branches that are no longer needed.
 *
 *
*/
unsigned long Grid::PruneTrees(
		PosType resolution  /// Maximum size of a cell to be removed.
		,bool useSB   /// If true it will not remove any point that has a flux above fluxlimit.
		,PosType fluxlimit ///  flux limit threshold
		){
		
	long i,Ntmp,count = 0;
	PosType res,initres;
	bool go;

	assert(trashkist);

	if(i_tree == NULL) return 0;
	if(s_tree == NULL) return 0;

	Ntmp = i_tree->pointlist->Npoints;

	//i_tree->moveTop();
      
      TreeStruct::iterator i_tree_current(i_tree);
      
	initres = (i_tree->getTop()->boundary_p2[0]-i_tree->getTop()->boundary_p1[0]);
	if(resolution > initres/3 || resolution <= 0.0) return 0;  // do not allow pruning up to the initial grid size

	// walk tree
	i=0;
	go = true;
	do{
		assert((*i_tree_current)->points->next || (*i_tree_current)->points->prev);

		res = ((*i_tree_current)->boundary_p2[0]-(*i_tree_current)->boundary_p1[0]);
		if( (res <= resolution && i_tree_current.IsSquareBranch() )
				 &&  (*i_tree_current)->refined){

			if(useSB){
				go = true;
				// Check if surface brightness of all points in cell are zero.
				i_tree->pointlist->current = (*i_tree_current)->points;
				for(i=0; i < (*i_tree_current)->npoints;++i,MoveDownList(i_tree->pointlist) ){
					if(i_tree->pointlist->current->surface_brightness*pow(i_tree->pointlist->current->gridsize,2) > fluxlimit ){
						go = false;
						break;
					}
				}
			}

			 // remove all lower branches and make current a leaf
			if(go && (*i_tree_current)->npoints > 1) count += FreeBranchesBelow(i_tree_current,i_tree,s_tree,trashkist);

		 }
	 }while(i_tree_current.TreeWalkStep(true));

	 // rebuild source tree from list.
	 //if(count > 0) RebuildTreeFromList(s_tree);

	 assert(count == (Ntmp - i_tree->pointlist->Npoints) );

	 // Trash collection
	 //CollectTrash(trashkist,true);

/*	 if(count > 10 && trashlist->Npoints > 10){
		 MoveToTopList(trashlist);
		 do{
			 // check to see if all points in the block have been removed from the trees
			 for(i=0;i<trashlist->current->head;++i) if(trashlist->current[i].leaf != NULL) break;
			 if(i == trashlist->current->head){
				 if(AtTopList(trashlist)) go = false; else go = true;
				 points = TakeOutCurrent(trashlist);
				 //std::printf("freeing memory!\n");
				 FreePointArray(points);
			 }
		 }while(MoveDownList(trashlist) && go);
	 }
*/

	 return count;
}



/** \ingroup ImageFindingL2
 *
 *   \brief Prune off points that are below a resolution and in an annulus on the
 *  source plane.
 *
 *  Used to keep the number of grid points limited while telescoping.
 *
 *  The points that are removed have cells that do not overlap the inner circle and centers
 *  that are within the outer circle.  Thus some points will be outside of the inner circle
 *  and some cells that are not removed may intersect with the outer circle.
 *
*/

unsigned long Grid::PrunePointsOutside(
		PosType resolution  /// Maximum size of a cell to be removed.
		,PosType *y   /// Center on source plane
		,PosType r_in /// Inner radius of annulus on the source plane
		,PosType r_out /// Outer radius of annulus on the source plane
		){

	if(i_tree == NULL) return 0;
	if(s_tree == NULL) return 0;
	if(r_in > r_out || resolution <= 0) return 0;
	if(r_out <= 0.0) return 0.0;

	PosType res,dr2;

	long i,count = 0;
	Point *point;
	Branch *branch;
	Kist<Point> * subkist = new Kist<Point>;
	bool go = true;
	unsigned long Ntmp;
	//Unit *unit;

	assert(trashkist);

	// Make a kist of all points in the annulus.
	//PointsWithinKist(s_tree,y,r_out+resolution,subkist,0);
	s_tree->PointsWithinKist_iter(y,r_in,r_out+resolution,subkist);

	std::printf("number of points after PointsWithin %li\n",subkist->Nunits());
	/*/*** test lines *******
	if(subkist->Nunits() > 0){
	MoveToTopKist(subkist);
	i=0;
	do{

		point = getCurrentKist(subkist);
		unit = subkist->current;

		assert(pow(r_out+resolution,2) >= pow(point->x[0] - y[0],2) + pow(point->x[1] - y[1],2)
				&& r_in*r_in <= pow(point->x[0] - y[0],2) + pow(point->x[1] - y[1],2));
		j=i;
		while(MoveDownKist(subkist)){
			++j;
			assert( getCurrentKist(subkist) != point );
		}
		subkist->current=unit;
		++i;
	}while(MoveDownKist(subkist));
	}
	/ *****************************************/

	if(r_in > 0){
		// take out points that are within inner circle
		subkist->MoveToTop();
		Ntmp = subkist->Nunits();
		for(i=0;i<Ntmp;++i){
			go = true;
			point = subkist->getCurrent();
			if(point->gridsize*Ngrid_block > resolution){
				if(subkist->AtTop()) go = false;
				subkist->TakeOutCurrent();
			}else if( Utilities::cutbox(y,point->leaf->boundary_p1,point->leaf->boundary_p2,r_in) ){
				if(subkist->AtTop()) go = false;
				subkist->TakeOutCurrent();
			}

			if(go) subkist->Down();
		}
	}
	std::printf("number of points after culling %li\n",subkist->Nunits());

	if(subkist->Nunits() == 0){
		delete subkist;
		return 0;
	}

	// move from source plane to image plane
	subkist->TranformPlanes();

	/*/*** test lines *******
	MoveToTopKist(subkist);
	do{
		point = getCurrentKist(subkist);
		unit = subkist->current;

		while(MoveDownKist(subkist)){
			assert( getCurrentKist(subkist) != point );
		}
		subkist->current=unit;

	}while(MoveDownKist(subkist));
	/ *****************************************/

      TreeStruct::iterator i_tree_current(i_tree);
      
	// Take out all points that are not at the center of their parent refined cell
	subkist->MoveToTop();
	Ntmp = subkist->Nunits();
	for(i = 0; i < Ntmp ; ++i){
		go = true;
		point = subkist->getCurrent();
		i_tree_current = point->leaf;
		// Move up to nearest ancestor that was refined.
		while(!((*i_tree_current)->refined) && i_tree_current.up() );

		res = ((*i_tree_current)->boundary_p2[0]-(*i_tree_current)->boundary_p1[0]);

		if((*i_tree_current)->npoints != Ngrid_block*Ngrid_block  || res > resolution){
			// Take out the point if it is not in a parent block that has been refined than once or if the parent block is about the resolution limit
			if(subkist->AtTop()) go = false;
			subkist->TakeOutCurrent();
		}else{
			assert(inbox(point->x,(*i_tree_current)->boundary_p1,(*i_tree_current)->boundary_p2));

			//assert(fabs((*i_tree_current)->center[0] - ((*i_tree_current)->boundary_p2[0]+(*i_tree_current)->boundary_p1[0])/2) < point->gridsize/2);
			//assert(fabs((*i_tree_current)->center[1] - ((*i_tree_current)->boundary_p2[1]+(*i_tree_current)->boundary_p1[1])/2) < point->gridsize/2 );
			assert(point->gridsize < fabs((*i_tree_current)->boundary_p2[0]-(*i_tree_current)->boundary_p1[0]) );
			assert(point->gridsize < fabs((*i_tree_current)->boundary_p2[1]-(*i_tree_current)->boundary_p1[1]) );

			 if( (point->gridsize)/2 < fabs(point->x[0] - (*i_tree_current)->center[0])
					 || (point->gridsize)/2 < fabs(point->x[1] - (*i_tree_current)->center[1])
								){ // Take out point if it is not at the center of it's parent block
				 if(subkist->AtTop()) go = false;
				 subkist->TakeOutCurrent();
			 }
		}

		if(go) subkist->Down();
	}

	/*/ ******** Test *************
	MoveToTopKist(subkist);
	do{
		point = getCurrentKist(subkist);
		unit = subkist->current;

		while(MoveDownKist(subkist)){
			assert( getCurrentKist(subkist) != point );
		}
		subkist->current=unit;

	}while(MoveDownKist(subkist));

	MoveToTopKist(subkist);
	do{
		i_tree->current = getCurrentKist(subkist)->leaf;
		// Move up to nearest ancestor that was refined.
		while(!((*i_tree_current)->refined) && moveUp(i_tree) );
		assert((*i_tree_current)->npoints == Ngrid_block*Ngrid_block);
	    unit = subkist->current;
	    branch = i_tree->current;
		while(MoveDownKist(subkist)){
			assert(unit->data != getCurrentKist(subkist));
			assert( !inbox(getCurrentKist(subkist)->x,(*i_tree_current)->boundary_p1,(*i_tree_current)->boundary_p2) );
			assert((*i_tree_current)->npoints == Ngrid_block*Ngrid_block);
		}
		subkist->current = unit;

		branch = i_tree->current;
		point = getCurrentKist(subkist);
		unit = subkist->current;
		while(!((*i_tree_current)->refined) && moveUp(i_tree) );
		assert((*i_tree_current)->npoints == Ngrid_block*Ngrid_block);
		branch2 = i_tree->current;
		while(MoveDownKist(subkist)){
			i_tree->current = getCurrentKist(subkist)->leaf;
			while(!((*i_tree_current)->refined) && moveUp(i_tree) );

			assert(branch2 != i_tree->current);
			assert((*i_tree_current)->npoints == Ngrid_block*Ngrid_block);
		}
		subkist->current = unit;
		i_tree->current = branch;


	}while(MoveDownKist(subkist));

	/ ******************************/

	assert(subkist->AtBottom());

	subkist->MoveToTop();
	while(subkist->Nunits() > 0){
		i_tree_current = subkist->getCurrent()->leaf;
		// Move up to nearest ancestor that was refined.
		while(!((*i_tree_current)->refined) && i_tree_current.up() );

		//assert((*i_tree_current)->npoints == Ngrid_block*Ngrid_block);

		/*/************* Test lines ******************
		point = getCurrentKist(subkist);
	    unit = subkist->current;
		MoveDownKist(subkist);
		while(MoveDownKist(subkist)){
			assert( !inbox(getCurrentKist(subkist)->x,(*i_tree_current)->boundary_p1,(*i_tree_current)->boundary_p2) );
		}
		subkist->current = unit;
		/ *******************************************/

		// make sure that none of the child points are outside the annulus
		i_tree->pointlist->current = (*i_tree_current)->points;
		for(i=0;i < Ngrid_block*Ngrid_block ; ++i,MoveDownList(i_tree->pointlist)){
			dr2 = pow(i_tree->pointlist->current->image->x[0] - y[0],2) + pow(i_tree->pointlist->current->image->x[1] - y[1],2);
			if(dr2 > r_out*r_out || dr2 < r_in*r_in) break;
			//if(dr2 < r_in*r_in) break;
		}

		if(i == Ngrid_block*Ngrid_block){
			branch = *i_tree_current;

			/*/************* Test lines ******************
			assert(branch == i_tree->current);

			if(subkist->Nunits() > 1){
				branch = i_tree->current;
				point = getCurrentKist(subkist);
				unit = subkist->current;
				while(!((*i_tree_current)->refined) && moveUp(i_tree) );
				assert((*i_tree_current)->npoints == Ngrid_block*Ngrid_block);
				branch2 = i_tree->current;
				while(MoveDownKist(subkist)){
					i_tree->current = getCurrentKist(subkist)->leaf;
					while(!((*i_tree_current)->refined) && moveUp(i_tree) );

					assert(branch2 != i_tree->current);
					assert((*i_tree_current)->npoints == Ngrid_block*Ngrid_block);
				}
				subkist->current = unit;
				i_tree->current = branch;

				point = getCurrentKist(subkist);
				unit = subkist->current;
				while(MoveDownKist(subkist)){
					assert( !inbox(getCurrentKist(subkist)->x,(*i_tree_current)->boundary_p1,(*i_tree_current)->boundary_p2) );
				}
				subkist->current = unit;

			}
			/ *******************************************/

			count += FreeBranchesBelow(i_tree_current,i_tree,s_tree,trashkist);

			/*/************* Test lines ******************
			assert(branch == i_tree->current);

			if(subkist->Nunits() > 1){
				branch = i_tree->current;
				point = getCurrentKist(subkist);
				unit = subkist->current;
				while(!((*i_tree_current)->refined) && moveUp(i_tree) );
				branch2 = i_tree->current;
				while(MoveDownKist(subkist)){
					i_tree->current = getCurrentKist(subkist)->leaf;
					while(!((*i_tree_current)->refined) && moveUp(i_tree) );

					assert(branch2 != i_tree->current);
					assert((*i_tree_current)->npoints == Ngrid_block*Ngrid_block);
				}
				subkist->current = unit;
				i_tree->current = branch;

				point = getCurrentKist(subkist);
				unit = subkist->current;
				while(MoveDownKist(subkist)){
					assert( !inbox(getCurrentKist(subkist)->x,(*i_tree_current)->boundary_p1,(*i_tree_current)->boundary_p2) );
				}
				subkist->current = unit;

			}
			/ *******************************************/

			assert((*i_tree_current)->npoints == 1);
			assert( fabs( 1 - ((*i_tree_current)->boundary_p2[0] - (*i_tree_current)->boundary_p1[0])/(*i_tree_current)->points->gridsize) < 1.0e-4 );
			assert((*i_tree_current)->refined == false);

		}

		subkist->TakeOutCurrent();

	}

	delete subkist;

	return count;
}

/** \ingroup ImageFindingL2
 *  \brief Empty trash points.
 *
 *  Frees point arrays whose heads are stored in trashlist.
 *  If check=true if will only free arrays where all the points have NULL leafs.
 *     check=false all the point arrays are freed..
 *
 */
void CollectTrash(Kist<Point> * trashkist,bool check){
	bool go;
	unsigned long i,j,Ntmp;
	Point *points;

	if(trashkist->Nunits() == 0) return;

	Ntmp = trashkist->Nunits();
	for(j=0,trashkist->MoveToTop();j<Ntmp;++j){

		 if(!(trashkist->getCurrent()->head)){  // point is not the head of an array
			 if(trashkist->AtTop()) go = false; else go = true;
			 trashkist->TakeOutCurrent();
		 }else{

			 // check to see if all points in the block have been removed from the trees
			 if(check){
				 for(i=0;i<trashkist->getCurrent()->head;++i) if(trashkist->getCurrent()[i].leaf != NULL) break;
			 }else{
				 i = trashkist->getCurrent()->head;
			 }

			 if(i == trashkist->getCurrent()->head){
				 if(trashkist->AtTop()) go = false; else go = true;
				 points = trashkist->TakeOutCurrent();
				 FreePointArray(points);
			 }else{
				 go = true;
			 }
		 }

		if(go) trashkist->Down();
	 }

	 return;
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

unsigned long FreeBranchesBelow(TreeStruct::iterator &i_tree_current,TreeHndl i_tree,TreeHndl s_tree,Kist<Point> * trashkist){

	if(!i_tree_current.IsSquareBranch()) return 0;
	if(i_tree_current.atLeaf()) return 0;
	if((*i_tree_current)->refined == false) return 0;

  TreeStruct::iterator s_tree_current(s_tree);
	assert( s_tree !=NULL);

	Branch *branch,*headbranch;
	Point *point;
	unsigned long Ntmp,NtoRemove,i,count = 0,count2 = 0,count1;
	PosType center[2];

	//_freeBranches_iter(s_tree);  // s_tree will no longer be valid on exit.  This is to make sure it isn't used later without a rebuild.

	headbranch = *i_tree_current;
	//i_tree->TreeWalkStep(true);

	while( (headbranch->child1 != NULL) || (headbranch->child2 != NULL) ){

		assert(boxinbox(*i_tree_current,headbranch));
		if(i_tree_current.atLeaf()){
			//assert(i_tree->current->points->image->leaf);
			//s_tree->current = i_tree->current->points->image->leaf;  // set s_tree to source of current image cell

			branch = (*i_tree_current)->prev;
			i = branch->npoints;

			if((*i_tree_current) != headbranch) i_tree->RemoveLeafFromTree(i_tree_current,&Ntmp);

			/***************** test line  **************************/
			assert(*i_tree_current == branch);
			assert(i == (*i_tree_current)->npoints);

			/***************** test line  **************************/
			assert((*i_tree_current)->points->next || (*i_tree_current)->points->prev);

			// in a square leaf cell take out extra points that have come up from below

			if(i_tree_current.atLeaf() && (*i_tree_current)->refined){
					/***************** test line  **************************/
				assert((*i_tree_current)->points->next || (*i_tree_current)->points->prev);

				//std::printf("  collecting points from removed leaves\n");
				assert((*i_tree_current)->points);
				i_tree->pointlist->current = (*i_tree_current)->points;
				NtoRemove = (*i_tree_current)->npoints;
				assert(NtoRemove == 9);
				center[0] = ((*i_tree_current)->boundary_p1[0] + (*i_tree_current)->boundary_p2[0])/2;
				center[1] = ((*i_tree_current)->boundary_p1[1] + (*i_tree_current)->boundary_p2[1])/2;

				for(i=0,count1=0,count2=0;i<NtoRemove;++i,MoveDownList(i_tree->pointlist)){
					// find central point and remove others

					if( (pow(center[0]-i_tree->pointlist->current->x[0],2)
						+ pow(center[1]-i_tree->pointlist->current->x[1],2) )
						< pow(i_tree->pointlist->current->gridsize/2,2) ){

						++count1;
						// keep this central point
						i_tree->pointlist->current->gridsize *= 3;
						i_tree->pointlist->current->image->gridsize = i_tree->pointlist->current->gridsize;
						(*i_tree_current)->points = i_tree->pointlist->current;

						/***************** test line  **************************/
						assert((*i_tree_current)->points->next || (*i_tree_current)->points->prev);
						assert(i_tree->pointlist->current->leaf == *i_tree_current);
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
						s_tree_current = point->leaf;
						//if((*s_tree_current)->npoints < 2) RemoveLeafFromTree(s_tree,&Ntmp);
						while(!(s_tree_current.atTop())){
							--((*s_tree_current)->npoints);
							if((*s_tree_current)->npoints > 0 && (*s_tree_current)->points == point)
								(*s_tree_current)->points = point->next;

							if((*s_tree_current)->npoints == 0) (*s_tree_current)->points = NULL;

							if((*s_tree_current)->npoints == 0 && (*s_tree_current)->prev->npoints == 1){
								// only remove empty leaves if it will make its parent a leaf
								assert(s_tree_current.atLeaf());
								s_tree->RemoveLeafFromTree(s_tree_current,&Ntmp);
								s_tree_current.TreeWalkStep(true);  // Go to other child.
								s_tree->RemoveLeafFromTree(s_tree_current,&Ntmp);
							}else{
								s_tree_current.up();
							}
						}
						assert(boxinbox(*i_tree_current,headbranch));

						// Do it for top
						--(s_tree->getTop()->npoints);
						if(s_tree->getTop()->npoints > 0 && s_tree->getTop()->points == point) s_tree->getTop()->points = point->next;

						s_tree->pointlist->current = point;
						TakeOutCurrent(s_tree->pointlist);
						point->leaf = NULL;  // set leaf to NULL to indicate that point is no longer in tree
						if(point->head) trashkist->InsertAfterCurrent(point);  // collect heads for later trash collection

						assert(boxinbox(*i_tree_current,headbranch));

						// take points out of image plane
						branch = *i_tree_current;
						do{
							assert((*i_tree_current)->npoints);
							--((*i_tree_current)->npoints);
							if((*i_tree_current)->points == i_tree->pointlist->current) (*i_tree_current)->points = NULL;

						}while(i_tree_current.up());
						i_tree_current = branch;
						assert(boxinbox(*i_tree_current,headbranch));

						if(i_tree->pointlist->current == (*i_tree_current)->points) (*i_tree_current)->points = i_tree->pointlist->current->next;
						point = TakeOutCurrent(i_tree->pointlist);
						point->leaf = NULL;
						// If point is a head of a memory block add it to trashlist for eventual trash collection
						if(point->head){
							assert(point->head == 8);
							trashkist->InsertAfterCurrent(point);
						}
					}
					assert(boxinbox(*i_tree_current,headbranch));

				}  // loop through points in leaf

				// reassign first point in branches above the current
				branch = *i_tree_current;
				do{
					assert((*i_tree_current)->npoints);
					if((*i_tree_current)->points == NULL) (*i_tree_current)->points = branch->points;
				}while(i_tree_current.up());
				i_tree_current = branch;

				assert(boxinbox(*i_tree_current,headbranch));

				assert(count1 == 1);
				assert(count2 == 8);
				assert((*i_tree_current)->npoints == 1);
			}  // if current was leaf that was refined

		} // at tree leaf
		assert((*i_tree_current)->points);
		assert(boxinbox(*i_tree_current,headbranch));

		if( !(i_tree_current.atLeaf()) ) i_tree_current.TreeWalkStep(true);
	}  // while entry current is not a leaf

	//assert(CurrentIsSquareTree(i_tree));
	assert((*i_tree_current)->npoints == 1);
	assert(i_tree_current.atLeaf());
	assert(*i_tree_current == headbranch);

	(*i_tree_current)->refined = false;

	// Free the memory for the points that have been removed.
	CollectTrash(trashkist,false);

	assert(trashkist->Nunits() == 0);
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
Point * TreeStruct::RemoveLeafFromTree(TreeStruct::iterator &current,unsigned long *Npoints){

	Branch *branch;
	Point *point = NULL;
	unsigned long i;

	if(current.atTop() || !(current.atLeaf()) ) return NULL;

	branch = *current;
  current.up();

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
		pointlist->current = branch->points;
		assert(inbox(pointlist->current->x,branch->boundary_p1,branch->boundary_p2));

		for(i=0;i<branch->npoints;++i,MoveDownList(pointlist)){
			pointlist->current->leaf = branch->prev;

			assert(inbox(pointlist->current->x,branch->boundary_p1,branch->boundary_p2));
			assert(inbox(pointlist->current->x,branch->prev->boundary_p1,branch->prev->boundary_p2));
		}

		assert(boxinbox(branch,branch->prev));

		point = branch->points;
	}
	*Npoints = branch->npoints;
	delete branch;
	--Nbranches;

	return point;
}

/// Move up to the parent of current branch if current in not the root. Otherwise returns false.
bool TreeStruct::iterator::up(){
  if(current == NULL || current == top){
    return false;
  }else{
    current = current->prev;
    return true;
  }
}

/// Move to brother if it exists
bool TreeStruct::iterator::brother(){
  if(current->brother == NULL){
    return false;
  }else{
    current = current->brother;
    return true;
  }
}

/// Move to child
bool TreeStruct::iterator::down(short child){
  if(child == 1){
    if(current->child1 == NULL){
      return false;
    }else{
      current = current->child1;
      return true;
    }
  }
  if(child == 2){
    if(current->child2 == NULL){
      return false;
    }else{
      current = current->child2;
      return true;
    }
  }
  
  throw std::runtime_error("There are only two children!");
}

/**
 *  \brief step for walking tree by iteration instead of recursion.
 *
 *  This walk will not exit the tree defined by the descendants of
 *  the root that is set in the constructor of TreeIt.  If allowed 
 *  to, it will return to the root and return false.  If used again 
 *  after this it will repeat its walk.
 */
bool TreeStruct::iterator::TreeWalkStep(bool allowDescent){
  
	if(allowDescent && current->child1 != NULL){
		down(1);
		return true;
	}
	if(allowDescent && current->child2 != NULL){
		down(2);
		return true;
	}
  
  if(current->brother == top->brother){
    current = top;
    return false;
  }
     
	if(current->brother != NULL){
		current = current->brother;
		return true;
	}
  
	return false;
}

