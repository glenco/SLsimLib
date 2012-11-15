/*
 * tree_maintenance.c
 *
 *  Created on: Sep 29, 2011
 *      Author: bmetcalf
 *
 *      This file contains routines for building, adding to, and removing grid points from
 *      the tree.
 */

/* median_cut determines how the cells are subdivided */
/*    if ==0  equal volume cuts, Warning this option causes an error*/
/*    if ==1  median point cuts */

#include "slsimlib.h"
/************************** test routine *****************************/
bool tree_count_test(TreeHndl tree){
	int i,nbranches=0;
	Branch *branch = tree->current;
	/*static int init=0;
	static double p1[2],p2[2];

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

	}while(TreeWalkStep(tree,true) && tree->current != branch->brother);

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

/** \ingroup  ImageFindingL2
 * \brief Build a complete tree from a list of points.
 */
TreeHndl BuildTree(Point *xp,unsigned long Npoints,short my_median_cut){
  TreeHndl tree;
  unsigned long i;
  double p1[2],p2[2],center[2];
  const int Nbucket = 1;   // must be =1 if each leaf is to coincide with each cell

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

  /* Initialize tree root */
  tree=NewTree(xp,Npoints,p1,p2,center,Nbucket);

  tree->median_cut = my_median_cut;
 /* build the tree */
  _BuildTree(tree);

  /************************** test routine *****************************/
  moveTop(tree);
  //tree_count_test(tree);

  return tree;
}
/** \ingroup  ImageFindingL2
 * \brief Fill a tree with points.  The previous tree structure will be destroyed.  Used for refilling.
 */
void FillTree(TreeHndl tree,Point *xp,unsigned long Npoints){
  unsigned long i;
  void _BuildTree(TreeHndl tree);

  assert(tree != NULL);

  tree->current->points=xp;
  tree->current->npoints=Npoints;
  // link point array into point list
  EmptyList(tree->pointlist);
  for(i=0;i<Npoints;++i){
    InsertPointAfterCurrent(tree->pointlist,&xp[i]);
    MoveDownList(tree->pointlist);
  }

  MoveToTopList(tree->pointlist);
 /* build the tree */

  // make sure there are no branches in the tree
  _freeBranches_iter(tree);

  _BuildTree(tree);

  return ;
}

/** \ingroup  ImageFindingL2
 * \brief Rebuilds the tree from the points that are already in the tree->pointlist
 *
 * This is not the best function because it copies all the points
 */
void RebuildTreeFromList(TreeHndl tree){
 /* Builds or rebuilds a tree that already has a tree->pointlist.
  */
	assert(tree != NULL);
	assert(tree->pointlist);
	if(tree->pointlist->Npoints < 1) return;

	unsigned long Npoints = tree->pointlist->Npoints,i;
	double *tmp;

	// Make new array of points
	Point *points = NewPointArray(Npoints,true);

	MoveToTopList(tree->pointlist);
	for(i=0;i<Npoints;++i,MoveDownList(tree->pointlist)){
		tmp = points[i].x;
		// PointCopy() copies the x pointer which is later freed in emptyTree() so the information needs to be copied
		PointCopy(&(points[i]),tree->pointlist->current);
		points[i].x = tmp;
		points[i].x[0] = tree->pointlist->current->x[0];
		points[i].x[1] = tree->pointlist->current->x[1];

		/* Below is for trash collection in PruneTree().
		 * When PruneTree has been used not all of the
		 * heads in the point arrays are guaranteed to be
		 * in the tree->pointlist.
		 */
		tree->pointlist->current->leaf = NULL;
	}
	//std::printf(" %e %e\n",points[0].x[0],points[0].x[1]);
	// emptry the tree and free all former points
	emptyTree(tree);

	assert(tree->Nbranches == 1);
	assert(tree->top->npoints == 0);
	assert(tree->pointlist->Npoints == 0);

	FillTree(tree,points,Npoints);

	return;
}

/** \ingroup  ImageFindingL2
* \brief Empty tree of all point leaving a tree with an empty root.
*
* The points are freed, but the list structure is not destroyed.
*
*  FillTree can then be used to regenerate tree.
*/
short emptyTree(TreeHndl tree){
  Point **heads;
  unsigned long i,j,count;

  heads = (Point **) malloc(tree->pointlist->Npoints*sizeof(Point*));  // maximum number of pointers that could be needed

  if(tree == NULL) return 1;

  assert(tree);

  _freeBranches_iter(tree);
  //_freeTree(tree,0);
  //printTree(tree);

  assert(tree->Nbranches == 1);

  MoveToTopList(tree->pointlist);
  for(i=0,j=0,count=0;i<tree->pointlist->Npoints;++i){
	  if(tree->pointlist->current->head){
		  heads[j] = tree->pointlist->current;
		  ++j;
		  count += tree->pointlist->current->head;
	  }
	  tree->pointlist->current->leaf = NULL;  // This is for future trash collection if PruneTree() has been used.
	  MoveDownList(tree->pointlist);
  }
  //assert(count == tree->pointlist->Npoints);  After using PruneTree this will no longer be guaranteed.

  //std::printf("freed %i arrays out of %i points in array, %i freed\n",j,i,count);
  for(i=0;i<j;++i) FreePointArray(heads[i]);
  tree->top->npoints = 0;

  //printTree(tree);

  tree->pointlist->Npoints=0;
  tree->pointlist->top = tree->pointlist->bottom = tree->pointlist->current=NULL;

  free(heads);
  return 1;
}
/** \ingroup LowLevel
* \brief Recursively free branches
*/
void _freeBranches(TreeHndl tree,short child){
	Branch *branch;

	assert( tree !=NULL);

	/*printBranch(tree->current);*/

	if(tree->current->child1 != NULL){
		moveToChild(tree,1);
		_freeBranches(tree,1);
	}

    if(tree->current->child2 != NULL){
      moveToChild(tree,2);
      _freeBranches(tree,2);
    }

    if( (tree->current->child1 == NULL)*(tree->current->child2 == NULL) ){

    	if(atTop(tree)){
    		//free(tree->current);
    		//tree->current=NULL;
    		//tree->top=NULL;
    		//--tree->Nbranches;
    		return;
    	}

    	branch=tree->current;
    	moveUp(tree);
       	free(branch);

    	/*std::printf("*** removing branch %i number of branches %i\n",branch->number
			,tree->Nbranches-1);*/

       	if(child==1) tree->current->child1=NULL;
    	if(child==2) tree->current->child2=NULL;
    	--tree->Nbranches;

    	return;
    }

    return;
}
/** \ingroup LowLevel
* \brief Iteratively free branches
*
* Frees all the barches of the tree so there is only the stump.
 */
void _freeBranches_iter(TreeHndl tree){
	Branch *branch;

	assert( tree !=NULL);
	moveTop(tree);

	/*printBranch(tree->current);*/

	while(tree->Nbranches > 1){

		if(tree->current->child1 != NULL){
			moveToChild(tree,1);
		}else if(tree->current->child2 != NULL){
			moveToChild(tree,2);
		}else{
			branch = tree->current;
			if(tree->current->brother == tree->current->prev->brother){
				moveUp(tree);
				assert(tree->current->child1 == NULL);
				tree->current->child2 = NULL;
			}else{
				tree->current = tree->current->brother;
				tree->current->prev->child1 = NULL;
			}

			free(branch);
			--tree->Nbranches;
		}
	}

	// re-assign the leaf pointers in the particles to the root
	 MoveToTopList(tree->pointlist);
	 do{ tree->pointlist->current->leaf = tree->top; }while(MoveDownList(tree->pointlist));

    return;
}


/** \ingroup LowLevel
* \brief Recursively build tree from points in its linked list.
*/
void _BuildTree(TreeHndl tree){
  /* tree->pointlist must be both a linked list and an array of points in the */
  /* same order as the linked list */
  unsigned long i,cut,dimension;
  Branch *cbranch,branch1,branch2;
  double xcut;

  cbranch=tree->current; /* pointer to current branch */

    /* leaf case */
  if(cbranch->npoints <= tree->Nbucket){
	  tree->current->points->leaf = tree->current;
	  if(cbranch->npoints == 0) cbranch->points = NULL;
	  return;
  }


  /* initialize boundaries to old boundaries */
  for(i=0;i<2;++i){
      branch1.boundary_p1[i]=cbranch->boundary_p1[i];
      branch1.boundary_p2[i]=cbranch->boundary_p2[i];

      branch2.boundary_p1[i]=cbranch->boundary_p1[i];
      branch2.boundary_p2[i]=cbranch->boundary_p2[i];
  }

  /* set dimension to cut box */
  dimension=(cbranch->level % 2);

  double *x = new double[cbranch->npoints];
  assert(x);

   /* reorder points */
  tree->pointlist->current=tree->current->points;
  for(i=0;i<cbranch->npoints;++i){
    x[i]=tree->pointlist->current->x[dimension];
    /*points[i]=tree->pointlist->current->id;*/
    MoveDownList(tree->pointlist);
  }

  /*PrintList(tree->pointlist);*/

  //double_sort(cbranch->npoints,x-1,points-1);
  //double_sort_points(cbranch->npoints,x-1,tree->current->points);

  /* copy information back to points in new order */
/*   tree->pointlist->current=tree->current->points; */
/*   for(i=0;i<cbranch->npoints;++i){ */
/*     tree->pointlist->current->x=xp[points[i]].x; */
/*     tree->pointlist->current->id=points[i]; */
/*     MoveDownList(tree->pointlist); */
/*   } */

  if(tree->median_cut){
	  double_sort_points(cbranch->npoints,x-1,tree->current->points);
	  //quicksortPoints(tree->current->points,x,cbranch->npoints);

	  cut=cbranch->npoints/2;
      branch1.boundary_p2[dimension]=(x[cut]+x[cut-1])/2;
      branch2.boundary_p1[dimension]=(x[cut]+x[cut-1])/2;

  }else{

	  xcut=(cbranch->boundary_p1[dimension]+cbranch->boundary_p2[dimension])/2;
      branch1.boundary_p2[dimension]=xcut;
      branch2.boundary_p1[dimension]=xcut;

	  quickPartitionPoints(xcut,&cut,tree->current->points,x,cbranch->npoints);

      //locateD(x-1,cbranch->npoints,xcut,&cut);
  }


  /* set point numbers and pointers to points */
  branch1.npoints=cut;
  assert(tree->current->points->next || tree->current->points->prev);
  branch1.points=tree->current->points;

  branch2.npoints=cbranch->npoints - cut;
  tree->pointlist->current=tree->current->points;
  JumpDownList(tree->pointlist,cut);
  branch2.points=tree->pointlist->current;

  delete[] x;


	/* use geometric center */
	branch1.center[0] = (branch1.boundary_p1[0] + branch1.boundary_p2[0])/2;
	branch1.center[1] = (branch1.boundary_p1[1] + branch1.boundary_p2[1])/2;

	branch2.center[0] = (branch2.boundary_p1[0] + branch2.boundary_p2[0])/2;
	branch2.center[1] = (branch2.boundary_p1[1] + branch2.boundary_p2[1])/2;

  /* centers of mass *

  for(i=0;i<2;++i) branch1.center[i]=0;
  tree->pointlist->current=branch1.points;
  for(i=0;i<cut; ++i){
	  branch1.center[0]+=tree->pointlist->current->x[0]/branch1.npoints;
	  branch1.center[1]+=tree->pointlist->current->x[1]/branch1.npoints;
	  MoveDownList(tree->pointlist);
  }

  for(i=0;i<2;++i) branch2.center[i]=0;
  tree->pointlist->current=branch2.points;
  for(i=cut;i<cbranch->npoints; ++i){
	  branch2.center[0]+=tree->pointlist->current->x[0]/branch2.npoints;
	  branch2.center[1]+=tree->pointlist->current->x[1]/branch2.npoints;
	  MoveDownList(tree->pointlist);
  }/**/

  attachChildrenToCurrent(tree,branch1,branch2);

  if( branch1.npoints > 0 ){
     //attachChildToCurrent(tree,branch1,1);
     moveToChild(tree,1);
     _BuildTree(tree);
     moveUp(tree);
  }

  if(branch2.npoints > 0 ){
     //attachChildToCurrent(tree,branch2,2);
	  moveToChild(tree,2);
	  _BuildTree(tree);
	  moveUp(tree);
  }

 /*std::printf("reached end of _BuildTree level=%i\n",tree->current->level);*/
  return;
}

/** \ingroup  ImageFindingL2
 *  \brief Expands tree by adding points
*/
int AddPointsToTree(TreeHndl tree,Point *xpoint,unsigned long Nadd){
  unsigned long j;
  //Branch *parent_branch;

  if(Nadd==0) return 1;

  for(j=0;j<Nadd;++j){

	  // Test lines
	   //if(!testLeafs(tree)){ERROR_MESSAGE(); std::cout << "Beginning of AddPointsToTree "<< std::endl; exit(1);}

	   // add only that are inside original grid
    	if( !inbox(xpoint[j].x,tree->top->boundary_p1,tree->top->boundary_p2) ){
    		ERROR_MESSAGE();
    		std::printf("ERROR: in AddPointToTree, ray is not inside the simulation box x = %e %e Nadd=%li\n  not adding it to tree\n",
    				   xpoint[j].x[0],xpoint[j].x[1],Nadd);
    		//std::printf("root of tree\n");
       		//printBranch(tree->top);
        	//exit(0);
    		return 0;
    	}else{

    		moveTop(tree);
     		_FindLeaf(tree,xpoint[j].x,1);
     		assert(atLeaf(tree));
     		//parent_branch = tree->current->prev->prev;

    		//assert(inbox(xpoint[j].x,parent_branch->boundary_p1,parent_branch->boundary_p2));
    		//assert(inbox(xpoint[j].x,tree->current->boundary_p1,tree->current->boundary_p2));

    		// insert point into point list
    		//if(tree->current->points == NULL){
    		if(tree->current->npoints == 1){
    			// case of no previous points in leaf

    			tree->current->points = &xpoint[j];
    			tree->current->points->leaf = tree->current;
    			// put point into right place in list
       			tree->pointlist->current = tree->current->prev->points;
    			if(tree->current == tree->current->prev->child1){
       				Point *oldpoint = tree->pointlist->current;
       				Branch *tmp = tree->current;

       				// insert the point
       				MoveUpList(tree->pointlist);
          			InsertPointAfterCurrent(tree->pointlist,tree->current->points);

          			// reassign parent first points
          			moveUp(tree);
           			while(tree->current->points == oldpoint){
           				tree->current->points = &xpoint[j];
           				moveUp(tree);
           			}
           			tree->current = tmp;

    			}else{
           			JumpDownList(tree->pointlist,tree->current->prev->npoints-2);  // adds point to end of parent branches list
        			InsertPointAfterCurrent(tree->pointlist,tree->current->points);
    			}

    			//tree->pointlist->current = tree->current->points;
       			//tree->current->points->leaf = tree->current;

       			assert(tree->current->points->next || tree->current->points->prev);
    		}else{
    			assert(tree->current->npoints > 1);

    			// case of points already in leaf
       			assert(tree->current->points->next || tree->current->points->prev);
       	    	// Test lines
       	    	//if(!testLeafs(tree)){ERROR_MESSAGE(); std::cout << "before _addPoint of AddPointsToTree "<< std::endl; exit(1);}

       			// adds point to end of branches list, note that tree->current->npoints has already been increased
    			tree->pointlist->current = tree->current->points;
    			JumpDownList(tree->pointlist,tree->current->npoints-2);
    			InsertPointAfterCurrent(tree->pointlist,&xpoint[j]);

    			tree->pointlist->current = tree->current->points;
    			xpoint[j].leaf = tree->current;

    			/*/ Test lines
    			tree->pointlist->current = tree->current->points;
    			for(int k = 0; k< tree->current->npoints ; ++k,MoveDownList(tree->pointlist)){
    				assert(inbox(tree->pointlist->current->x
    					,tree->current->boundary_p1,tree->current->boundary_p2));
    			}
    			///////////////////////////////////*/

          		Branch *parent_branch = tree->current;
       			_AddPoint(tree);
      			assert(parent_branch->child1->points->leaf = parent_branch->child1);
      			assert(parent_branch->child2->points->leaf = parent_branch->child2);
      	    	assert(inbox(xpoint[j].x,xpoint[j].leaf->boundary_p1,xpoint[j].leaf->boundary_p2));
      	    	/*/  Test lines
      	    	if(!testLeafs(tree)){
      	    		ERROR_MESSAGE();
      	    		tree->pointlist->current = parent_branch->points;
      	    		for(int i=0;i<parent_branch->npoints;++i,MoveDownList(tree->pointlist)){
      	    			std::cout << "points in split branch " << tree->pointlist->current->id << std::endl;
      	    		}
      	    		std::cout << "Adding point " << xpoint[j].id << " branch " << parent_branch << std::endl;
      	    		std::cout << "End of AddPointsToTree "<< std::endl; exit(1);
      	    	}
      	    	//////////////////////////////////////////*/

    		}

    	}

    	while( !inbox(xpoint[j].x,tree->current->boundary_p1,tree->current->boundary_p2) ) moveUp(tree);
    	_FindLeaf(tree,xpoint[j].x,0);
    	xpoint[j].leaf = tree->current;
    	assert(inbox(xpoint[j].x,xpoint[j].leaf->boundary_p1,xpoint[j].leaf->boundary_p2));
    	/*/ Test lines
    	if(!testLeafs(tree)){ERROR_MESSAGE(); std::cout << "End of AddPointsToTree "<< std::endl; exit(1);}
    	///////////////////////*/
  }
    return 1;
}

void _AddPoint(TreeHndl tree){

		/*/ Test lines
		tree->pointlist->current = tree->current->points;
		for(int k = 0; k< tree->current->npoints ; ++k,MoveDownList(tree->pointlist)){
			assert(inbox(tree->pointlist->current->x
				,tree->current->boundary_p1,tree->current->boundary_p2));
		}
		///////////////////////////////////*/

	unsigned long i;
	if(tree->current->npoints == 0){
		return;
	}else if(tree->current->npoints <= tree->Nbucket){
		tree->pointlist->current = tree->current->points;
		for(i=0;i<tree->current->npoints;i++){
			tree->pointlist->current->leaf = tree->current;
			assert(inbox(tree->pointlist->current->x
					,tree->pointlist->current->leaf->boundary_p1
					,tree->pointlist->current->leaf->boundary_p2));
			MoveDownList(tree->pointlist);
		}
		return;
	}

	if(atLeaf(tree)){
		Branch branch1,branch2;
		Branch *current;
		unsigned long dimension,cut;
		double xcut;
		Point *oldfirstpoint,*newfirstpoint;
		double *x = new double[tree->current->npoints];

		// initialize boundaries to old boundaries
		for(i=0;i<2;++i){
			branch1.boundary_p1[i] = tree->current->boundary_p1[i];
			branch1.boundary_p2[i] = tree->current->boundary_p2[i];

			branch2.boundary_p1[i] = tree->current->boundary_p1[i];
			branch2.boundary_p2[i] = tree->current->boundary_p2[i];
		}

		/*/ Test lines
		tree->pointlist->current = tree->current->points;
		for(int k = 0; k< tree->current->npoints ; ++k,MoveDownList(tree->pointlist)){
			assert(inbox(tree->pointlist->current->x
				,tree->current->boundary_p1,tree->current->boundary_p2));
		}
		tree->pointlist->current = tree->current->points;
		///////////////////////////////////*/

		// set dimension to cut box
		dimension=(tree->current->level % 2);

		// reorder points
		tree->pointlist->current = tree->current->points;

		assert(tree->current->points->next || tree->current->points->prev);

		for(i=0;i<tree->current->npoints;++i){
			x[i] = tree->pointlist->current->x[dimension];
			MoveDownList(tree->pointlist);
		}

		/*/ Test lines
		tree->pointlist->current = tree->current->points;
		for(int k = 0; k< tree->current->npoints ; ++k,MoveDownList(tree->pointlist)){
			assert(inbox(tree->pointlist->current->x
				,tree->current->boundary_p1,tree->current->boundary_p2));
		}
		tree->pointlist->current = tree->current->points;
		///////////////////////////////////*/

		oldfirstpoint = tree->current->points;
		if(tree->current->npoints == 2){
			tree->pointlist->current = tree->current->points;
			if(tree->pointlist->current->x[dimension] > tree->pointlist->current->next->x[dimension]){
				tree->pointlist->current = tree->current->points;
				bool attop=AtTopList(tree->pointlist);
				Point *point = TakeOutCurrent(tree->pointlist);
				if(!attop) MoveDownList(tree->pointlist);
				InsertPointAfterCurrent(tree->pointlist,point);
				tree->current->points = point->prev;
				double tmp = x[1];
				x[1] = x[0];
				x[0] = tmp;
			}
		}else{
			ERROR_MESSAGE();
			std::cout << "This is prone to errors and this should never happen!" << std::endl;
			exit(1);
			tree->current->points = sortList(tree->current->npoints,x,tree->pointlist,tree->current->points);
		}
		newfirstpoint = tree->current->points;
		/*/ Test lines
		tree->pointlist->current = tree->current->points;
		for(int k = 0; k< tree->current->npoints ; ++k,MoveDownList(tree->pointlist)){
			assert(inbox(tree->pointlist->current->x
				,tree->current->boundary_p1,tree->current->boundary_p2));
		}
		tree->pointlist->current = tree->current->points;
		///////////////////////////////////*/

		cut = tree->current->npoints/2;

		// check that median split in this dimension will split particles
		if(tree->median_cut && x[cut] == x[cut-1]){
			// change dimension

			dimension=!dimension;

			tree->pointlist->current = tree->current->points;
			for(i=0;i<tree->current->npoints;i++){
				x[i]=tree->pointlist->current->x[dimension];
				MoveDownList(tree->pointlist);
			}

			if(tree->current->npoints == 2){
				tree->pointlist->current = tree->current->points;
				if(tree->pointlist->current->x[dimension] > tree->pointlist->current->next->x[dimension]){
					tree->pointlist->current = tree->current->points;
					bool attop=AtTopList(tree->pointlist);
					Point *point = TakeOutCurrent(tree->pointlist);
					if(!attop) MoveDownList(tree->pointlist);
					InsertPointAfterCurrent(tree->pointlist,point);
					tree->current->points = point->prev;
					double tmp = x[1];
					x[1] = x[0];
					x[0] = tmp;
				}
			}else{
				ERROR_MESSAGE();
				std::cout << "This is prone to errors and this should never happen!" << std::endl;
				exit(1);
				tree->current->points = sortList(tree->current->npoints,x,tree->pointlist,tree->current->points);
			}
			newfirstpoint=tree->current->points;
		}

		assert(x[0] != x[1]);
		if(tree->median_cut){
			cut=tree->current->npoints/2;

			branch1.boundary_p2[dimension]=(x[cut]+x[cut-1])/2;
			branch2.boundary_p1[dimension]=(x[cut]+x[cut-1])/2;

		}else{
			xcut=(tree->current->boundary_p1[dimension]+tree->current->boundary_p2[dimension])/2;
			branch1.boundary_p2[dimension]=xcut;
			branch2.boundary_p1[dimension]=xcut;

			locateD(x-1,tree->current->npoints,xcut,&cut);
		}
		delete[] x;
		/*/  Test lines
		tree->pointlist->current = tree->current->points;
		for(int k = 0; k< tree->current->npoints ; ++k,MoveDownList(tree->pointlist)){
			assert(inbox(tree->pointlist->current->x
				,tree->current->boundary_p1,tree->current->boundary_p2));
		}
		tree->pointlist->current = tree->current->points;
		///////////////////////////////////*/

		assert(cut <= tree->current->npoints);

		// set point numbers and pointers to points
		branch1.npoints = cut;
		if(branch1.npoints > 0) branch1.points = tree->current->points;
		else branch1.points = NULL;

		branch2.npoints = tree->current->npoints - cut;
		tree->pointlist->current = tree->current->points;

		if(branch2.npoints > 0){
			JumpDownList(tree->pointlist,cut);
			branch2.points = tree->pointlist->current;
		}else{
			branch2.points = NULL;
		}
		/*/ Test lines
		tree->pointlist->current = tree->current->points;
		for(int k = 0; k< tree->current->npoints ; ++k,MoveDownList(tree->pointlist)){
			assert(inbox(tree->pointlist->current->x
				,tree->current->boundary_p1,tree->current->boundary_p2));
		}
		tree->pointlist->current = tree->current->points;
		///////////////////////////////////*/

		// use geometric center
		branch1.center[0] = (branch1.boundary_p1[0] + branch1.boundary_p2[0])/2;
		branch1.center[1] = (branch1.boundary_p1[1] + branch1.boundary_p2[1])/2;

		branch2.center[0] = (branch2.boundary_p1[0] + branch2.boundary_p2[0])/2;
		branch2.center[1] = (branch2.boundary_p1[1] + branch2.boundary_p2[1])/2;

		/*/ Test lines
		if(branch1.npoints > 0 && !inbox(branch1.points->x,branch1.boundary_p1,branch1.boundary_p2)){
			printBranch(tree->current);
			printBranch(&branch1);
			printBranch(&branch2);
			tree->pointlist->current = tree->current->points;
			for(i=0;i<tree->current->npoints;i++){
				std::cout << tree->pointlist->current->x[0] << "  " << tree->pointlist->current->x[1] << std::endl;
				MoveDownList(tree->pointlist);
			}
			ERROR_MESSAGE();
			exit(0);
		}
		if(branch2.npoints > 0 && !inbox(branch2.points->x,branch2.boundary_p1,branch2.boundary_p2)){
			printBranch(tree->current);
			printBranch(&branch1);
			printBranch(&branch2);
			tree->pointlist->current = tree->current->points;
			for(i=0;i<tree->current->npoints;i++){
				std::cout << tree->pointlist->current->x[0] << "  " << tree->pointlist->current->x[1] << std::endl;
				MoveDownList(tree->pointlist);
			}
			ERROR_MESSAGE();
			exit(0);
		}
		/////////////////////////////////////////////////////
		/* centers of mass *

	for(i=0;i<2;++i) branch1.center[i]=0;
	tree->pointlist->current=branch1.points;
	for(i=0;i<cut; ++i){
		branch1.center[0]+=tree->pointlist->current->x[0]/branch1.npoints;
		branch1.center[1]+=tree->pointlist->current->x[1]/branch1.npoints;
		MoveDownList(tree->pointlist);
	}

	for(i=0;i<2;++i) branch2.center[i]=0;
	tree->pointlist->current=branch2.points;
	for(i=cut;i<tree->current->npoints; ++i){
		branch2.center[0]+=tree->pointlist->current->x[0]/branch2.npoints;
		branch2.center[1]+=tree->pointlist->current->x[1]/branch2.npoints;
		MoveDownList(tree->pointlist);
	}
	*/

		attachChildrenToCurrent(tree,branch1,branch2);

		assert( (tree->current->child1->npoints + tree->current->child2->npoints)
				== tree->current->npoints );
		assert(tree->current->child1->npoints == 1 && tree->current->child2->npoints == 1);

		// If needed, change the particle pointer in all parent cells
		current = tree->current;
		moveUp(tree);
		while(tree->current->points == oldfirstpoint){
			tree->current->points = newfirstpoint;
			if(!moveUp(tree)) break;
		}
		tree->current = current;
	}

	// recursively descent to children
	if( tree->current->child1 != NULL ){
		moveToChild(tree,1);
		_AddPoint(tree);
		moveUp(tree);
	}

	if( tree->current->child2 != NULL ){
		moveToChild(tree,2);
		_AddPoint(tree);
		moveUp(tree);
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
		double resolution  /// Maximum size of a cell to be removed.
		,bool useSB   /// If true it will not remove any point that has a flux above fluxlimit.
		,double fluxlimit ///  flux limit threshold
		){
		
	long i,Ntmp,count = 0;
	double res,initres;
	bool go;

	assert(trashkist);

	if(i_tree == NULL) return 0;
	if(s_tree == NULL) return 0;

	Ntmp = i_tree->pointlist->Npoints;

	moveTop(i_tree);
	initres = (i_tree->top->boundary_p2[0]-i_tree->top->boundary_p1[0]);
	if(resolution > initres/3 || resolution <= 0.0) return 0;  // do not allow pruning up to the initial grid size

	// walk tree
	i=0;
	go = true;
	do{
		assert(i_tree->current->points->next || i_tree->current->points->prev);

		res = (i_tree->current->boundary_p2[0]-i_tree->current->boundary_p1[0]);
		if( (res <= resolution && CurrentIsSquareTree(i_tree) )
				 &&  i_tree->current->refined){

			if(useSB){
				go = true;
				// Check if surface brightness of all points in cell are zero.
				i_tree->pointlist->current = i_tree->current->points;
				for(i=0; i < i_tree->current->npoints;++i,MoveDownList(i_tree->pointlist) ){
					if(i_tree->pointlist->current->surface_brightness*pow(i_tree->pointlist->current->gridsize,2) > fluxlimit ){
						go = false;
						break;
					}
				}
			}

			 // remove all lower branches and make current a leaf
			if(go && i_tree->current->npoints > 1) count += FreeBranchesBelow(i_tree,s_tree,trashkist);

		 }
	 }while(TreeWalkStep(i_tree,true));

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
	 moveTop(i_tree);
	 moveTop(s_tree);

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
		double resolution  /// Maximum size of a cell to be removed.
		,double *y   /// Center on source plane
		,double r_in /// Inner radius of annulus on the source plane
		,double r_out /// Outer radius of annulus on the source plane
		){

	if(i_tree == NULL) return 0;
	if(s_tree == NULL) return 0;
	if(r_in > r_out || resolution <= 0) return 0;
	if(r_out <= 0.0) return 0.0;

	double res,dr2;
	TreeHndl i_tree = i_tree,s_tree = s_tree;
	int Ngrid_block = getNgrid_block();

	long i,count = 0,j;
	Point *point;
	Branch *branch,*branch2;
	KistHndl subkist = new Kist;
	bool go = true;
	unsigned long Ntmp;
	Unit *unit;

	assert(trashkist);

	// Make a kist of all points in the annulus.
	//PointsWithinKist(s_tree,y,r_out+resolution,subkist,0);
	PointsWithinKist_iter(s_tree,y,r_in,r_out+resolution,subkist);

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
	/*****************************************/

	if(r_in > 0){
		// take out points that are within inner circle
		MoveToTopKist(subkist);
		Ntmp = subkist->Nunits();
		for(i=0;i<Ntmp;++i){
			go = true;
			point = getCurrentKist(subkist);
			if(point->gridsize*Ngrid_block > resolution){
				if(AtTopKist(subkist)) go = false;
				TakeOutCurrentKist(subkist);
			}else if( cutbox(y,point->leaf->boundary_p1,point->leaf->boundary_p2,r_in) ){
				if(AtTopKist(subkist)) go = false;
				TakeOutCurrentKist(subkist);
			}

			if(go) MoveDownKist(subkist);
		}
	}
	std::printf("number of points after culling %li\n",subkist->Nunits());

	if(subkist->Nunits() == 0){
		delete subkist;
		return 0;
	}

	// move from source plane to image plane
	TranformPlanesKist(subkist);

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
	/*****************************************/

	// Take out all points that are not at the center of their parent refined cell
	MoveToTopKist(subkist);
	Ntmp = subkist->Nunits();
	for(i = 0; i < Ntmp ; ++i){
		go = true;
		point = getCurrentKist(subkist);
		i_tree->current = point->leaf;
		// Move up to nearest ancestor that was refined.
		while(!(i_tree->current->refined) && moveUp(i_tree) );

		res = (i_tree->current->boundary_p2[0]-i_tree->current->boundary_p1[0]);

		if(i_tree->current->npoints != Ngrid_block*Ngrid_block  || res > resolution){
			// Take out the point if it is not in a parent block that has been refined than once or if the parent block is about the resolution limit
			if(AtTopKist(subkist)) go = false;
			TakeOutCurrentKist(subkist);
		}else{
			assert(inbox(point->x,i_tree->current->boundary_p1,i_tree->current->boundary_p2));

			//assert(fabs(i_tree->current->center[0] - (i_tree->current->boundary_p2[0]+i_tree->current->boundary_p1[0])/2) < point->gridsize/2);
			//assert(fabs(i_tree->current->center[1] - (i_tree->current->boundary_p2[1]+i_tree->current->boundary_p1[1])/2) < point->gridsize/2 );
			assert(point->gridsize < fabs(i_tree->current->boundary_p2[0]-i_tree->current->boundary_p1[0]) );
			assert(point->gridsize < fabs(i_tree->current->boundary_p2[1]-i_tree->current->boundary_p1[1]) );

			 if( (point->gridsize)/2 < fabs(point->x[0] - i_tree->current->center[0])
					 || (point->gridsize)/2 < fabs(point->x[1] - i_tree->current->center[1])
								){ // Take out point if it is not at the center of it's parent block
				 if(AtTopKist(subkist)) go = false;
				 TakeOutCurrentKist(subkist);
			 }
		}

		if(go) MoveDownKist(subkist);
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
		while(!(i_tree->current->refined) && moveUp(i_tree) );
		assert(i_tree->current->npoints == Ngrid_block*Ngrid_block);
	    unit = subkist->current;
	    branch = i_tree->current;
		while(MoveDownKist(subkist)){
			assert(unit->data != getCurrentKist(subkist));
			assert( !inbox(getCurrentKist(subkist)->x,i_tree->current->boundary_p1,i_tree->current->boundary_p2) );
			assert(i_tree->current->npoints == Ngrid_block*Ngrid_block);
		}
		subkist->current = unit;

		branch = i_tree->current;
		point = getCurrentKist(subkist);
		unit = subkist->current;
		while(!(i_tree->current->refined) && moveUp(i_tree) );
		assert(i_tree->current->npoints == Ngrid_block*Ngrid_block);
		branch2 = i_tree->current;
		while(MoveDownKist(subkist)){
			i_tree->current = getCurrentKist(subkist)->leaf;
			while(!(i_tree->current->refined) && moveUp(i_tree) );

			assert(branch2 != i_tree->current);
			assert(i_tree->current->npoints == Ngrid_block*Ngrid_block);
		}
		subkist->current = unit;
		i_tree->current = branch;


	}while(MoveDownKist(subkist));

	/******************************/

	assert(AtBottomKist(subkist));

	MoveToTopKist(subkist);
	while(subkist->Nunits() > 0){
		i_tree->current = getCurrentKist(subkist)->leaf;
		// Move up to nearest ancestor that was refined.
		while(!(i_tree->current->refined) && moveUp(i_tree) );

		//assert(i_tree->current->npoints == Ngrid_block*Ngrid_block);

		/*/************* Test lines ******************
		point = getCurrentKist(subkist);
	    unit = subkist->current;
		MoveDownKist(subkist);
		while(MoveDownKist(subkist)){
			assert( !inbox(getCurrentKist(subkist)->x,i_tree->current->boundary_p1,i_tree->current->boundary_p2) );
		}
		subkist->current = unit;
		//*******************************************/

		// make sure that none of the child points are outside the annulus
		i_tree->pointlist->current = i_tree->current->points;
		for(i=0;i < Ngrid_block*Ngrid_block ; ++i,MoveDownList(i_tree->pointlist)){
			dr2 = pow(i_tree->pointlist->current->image->x[0] - y[0],2) + pow(i_tree->pointlist->current->image->x[1] - y[1],2);
			if(dr2 > r_out*r_out || dr2 < r_in*r_in) break;
			//if(dr2 < r_in*r_in) break;
		}

		if(i == Ngrid_block*Ngrid_block){
			branch = i_tree->current;

			/*/************* Test lines ******************
			assert(branch == i_tree->current);

			if(subkist->Nunits() > 1){
				branch = i_tree->current;
				point = getCurrentKist(subkist);
				unit = subkist->current;
				while(!(i_tree->current->refined) && moveUp(i_tree) );
				assert(i_tree->current->npoints == Ngrid_block*Ngrid_block);
				branch2 = i_tree->current;
				while(MoveDownKist(subkist)){
					i_tree->current = getCurrentKist(subkist)->leaf;
					while(!(i_tree->current->refined) && moveUp(i_tree) );

					assert(branch2 != i_tree->current);
					assert(i_tree->current->npoints == Ngrid_block*Ngrid_block);
				}
				subkist->current = unit;
				i_tree->current = branch;

				point = getCurrentKist(subkist);
				unit = subkist->current;
				while(MoveDownKist(subkist)){
					assert( !inbox(getCurrentKist(subkist)->x,i_tree->current->boundary_p1,i_tree->current->boundary_p2) );
				}
				subkist->current = unit;

			}
			/*******************************************/

			count += FreeBranchesBelow(i_tree,s_tree,trashkist);

			/*/************* Test lines ******************
			assert(branch == i_tree->current);

			if(subkist->Nunits() > 1){
				branch = i_tree->current;
				point = getCurrentKist(subkist);
				unit = subkist->current;
				while(!(i_tree->current->refined) && moveUp(i_tree) );
				branch2 = i_tree->current;
				while(MoveDownKist(subkist)){
					i_tree->current = getCurrentKist(subkist)->leaf;
					while(!(i_tree->current->refined) && moveUp(i_tree) );

					assert(branch2 != i_tree->current);
					assert(i_tree->current->npoints == Ngrid_block*Ngrid_block);
				}
				subkist->current = unit;
				i_tree->current = branch;

				point = getCurrentKist(subkist);
				unit = subkist->current;
				while(MoveDownKist(subkist)){
					assert( !inbox(getCurrentKist(subkist)->x,i_tree->current->boundary_p1,i_tree->current->boundary_p2) );
				}
				subkist->current = unit;

			}
			/*******************************************/

			assert(i_tree->current->npoints == 1);
			assert( fabs( 1 - (i_tree->current->boundary_p2[0] - i_tree->current->boundary_p1[0])/i_tree->current->points->gridsize) < 1.0e-4 );
			assert(i_tree->current->refined == false);

		}

		TakeOutCurrentKist(subkist);

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
void CollectTrash(KistHndl trashkist,bool check){
	bool go;
	unsigned long i,j,Ntmp;
	Point *points;

	if(trashkist->Nunits() == 0) return;

	Ntmp = trashkist->Nunits();
	for(j=0,MoveToTopKist(trashkist);j<Ntmp;++j){

		 if(!(getCurrentKist(trashkist)->head)){  // point is not the head of an array
			 if(AtTopKist(trashkist)) go = false; else go = true;
			 TakeOutCurrentKist(trashkist);
		 }else{

			 // check to see if all points in the block have been removed from the trees
			 if(check){
				 for(i=0;i<getCurrentKist(trashkist)->head;++i) if(getCurrentKist(trashkist)[i].leaf != NULL) break;
			 }else{
				 i = getCurrentKist(trashkist)->head;
			 }

			 if(i == getCurrentKist(trashkist)->head){
				 if(AtTopKist(trashkist)) go = false; else go = true;
				 points = TakeOutCurrentKist(trashkist);
				 FreePointArray(points);
			 }else{
				 go = true;
			 }
		 }

		if(go) MoveDownKist(trashkist);
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

unsigned long FreeBranchesBelow(TreeHndl i_tree,TreeHndl s_tree,KistHndl trashkist){

	if(!CurrentIsSquareTree(i_tree)) return 0;
	if(atLeaf(i_tree)) return 0;
	if(i_tree->current->refined == false) return 0;

	assert( i_tree !=NULL);
	assert( s_tree !=NULL);

	Branch *branch,*headbranch;
	Point *point;
	unsigned long Ntmp,NtoRemove,i,count = 0,count2 = 0,count1,j;
	double center[2];

	//_freeBranches_iter(s_tree);  // s_tree will no longer be valid on exit.  This is to make sure it isn't used later without a rebuild.

	assert(i_tree->current->refined);
	assert(CurrentIsSquareTree(i_tree));
	headbranch = i_tree->current;
	//TreeWalkStep(i_tree,true);

	while( (headbranch->child1 != NULL) || (headbranch->child2 != NULL) ){

		assert(boxinbox(i_tree->current,headbranch));
		if(atLeaf(i_tree)){
			//assert(i_tree->current->points->image->leaf);
			//s_tree->current = i_tree->current->points->image->leaf;  // set s_tree to source of current image cell

			/***************** test line  **************************/
			assert(i_tree->current->points->next || i_tree->current->points->prev);

			/***************** test line  **************************/
			branch = i_tree->current->prev;
			i = branch->npoints;

			if(i_tree->current != headbranch) RemoveLeafFromTree(i_tree,&Ntmp);

			/***************** test line  **************************/
			assert(i_tree->current == branch);
			assert(i == i_tree->current->npoints);

			/***************** test line  **************************/
			assert(i_tree->current->points->next || i_tree->current->points->prev);

			// in a square leaf cell take out extra points that have come up from below

			if(atLeaf(i_tree) && i_tree->current->refined){
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
						while(!atTop(s_tree)){
							--(s_tree->current->npoints);
							if(s_tree->current->npoints > 0 && s_tree->current->points == point)
								s_tree->current->points = point->next;

							if(s_tree->current->npoints == 0) s_tree->current->points = NULL;

							if(s_tree->current->npoints == 0 && s_tree->current->prev->npoints == 1){
								// only remove empty leaves if it will make its parent a leaf
								assert(atLeaf(s_tree));
								RemoveLeafFromTree(s_tree,&Ntmp);
								TreeWalkStep(s_tree,true);  // Go to other child.
								RemoveLeafFromTree(s_tree,&Ntmp);
							}else{
								moveUp(s_tree);
							}
						}
						assert(boxinbox(i_tree->current,headbranch));

						// Do it for top
						--(s_tree->top->npoints);
						if(s_tree->top->npoints > 0 && s_tree->top->points == point) s_tree->top->points = point->next;

						s_tree->pointlist->current = point;
						TakeOutCurrent(s_tree->pointlist);
						point->leaf = NULL;  // set leaf to NULL to indicate that point is no longer in tree
						if(point->head) InsertAfterCurrentKist(trashkist,point);  // collect heads for later trash collection

						assert(boxinbox(i_tree->current,headbranch));

						// take points out of image plane
						branch = i_tree->current;
						do{
							assert(i_tree->current->npoints);
							--(i_tree->current->npoints);
							if(i_tree->current->points == i_tree->pointlist->current) i_tree->current->points = NULL;

						}while(moveUp(i_tree));
						i_tree->current = branch;
						assert(boxinbox(i_tree->current,headbranch));

						if(i_tree->pointlist->current == i_tree->current->points) i_tree->current->points = i_tree->pointlist->current->next;
						point = TakeOutCurrent(i_tree->pointlist);
						point->leaf = NULL;
						// If point is a head of a memory block add it to trashlist for eventual trash collection
						if(point->head){
							assert(point->head == 8);
							InsertAfterCurrentKist(trashkist,point);
						}
					}
					assert(boxinbox(i_tree->current,headbranch));

				}  // loop through points in leaf

				// reassign first point in branches above the current
				branch = i_tree->current;
				do{
					assert(i_tree->current->npoints);
					if(i_tree->current->points == NULL) i_tree->current->points = branch->points;
				}while(moveUp(i_tree));
				i_tree->current = branch;

				assert(boxinbox(i_tree->current,headbranch));

				assert(count1 == 1);
				assert(count2 == 8);
				assert(i_tree->current->npoints == 1);
			}  // if current was leaf that was refined

		} // at tree leaf
		assert(i_tree->current->points);
		assert(boxinbox(i_tree->current,headbranch));

		if( !(atLeaf(i_tree)) )TreeWalkStep(i_tree,true);
	}  // while entry current is not a leaf

	//assert(CurrentIsSquareTree(i_tree));
	assert(i_tree->current->npoints == 1);
	assert(atLeaf(i_tree));
	assert(i_tree->current == headbranch);

	i_tree->current->refined = false;

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
Point *RemoveLeafFromTree(TreeHndl tree,unsigned long *Npoints){

	Branch *branch;
	Point *point;
	unsigned long i;

	if(atTop(tree) || !(atLeaf(tree)) ) return NULL;

	branch = tree->current;
	moveUp(tree);

	//if(tree->current->number == 2125226) branchaddress = NULL;

	if(branch == branch->prev->child1){
		branch->prev->child1 = NULL;
	}

	if(branch == branch->prev->child2){
		if(branch->prev->child1 != NULL) branch->prev->child1->brother = branch->prev->brother;
		branch->prev->child2 = NULL;
	}

	// leaves of points in the father
	if(branch->npoints >0){
		tree->pointlist->current = branch->points;
		assert(inbox(tree->pointlist->current->x,branch->boundary_p1,branch->boundary_p2));

		for(i=0;i<branch->npoints;++i,MoveDownList(tree->pointlist)){
			tree->pointlist->current->leaf = branch->prev;

			assert(inbox(tree->pointlist->current->x,branch->boundary_p1,branch->boundary_p2));
			assert(inbox(tree->pointlist->current->x,branch->prev->boundary_p1,branch->prev->boundary_p2));
		}

		assert(boxinbox(branch,branch->prev));

		point = branch->points;
	}
	*Npoints = branch->npoints;
	free(branch);
	--(tree->Nbranches);

	return point;
}


