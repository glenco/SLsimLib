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
/*#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <nrD.h>
#include <Tree.h>
#include <nrutil.h>
#include <tree_maintenance.h>*/

#include <slsimlib.h>

static int median_cut=1;
//const int Ngrid_block = 3;

/** \ingroup  ImageFindingL2
 * \brief Build a complete tree from a list of points.
 */
TreeHndl BuildTree(Point *xp,unsigned long Npoints){
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

 /* build the tree */
  _BuildTree(tree);

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
 * This is not the best function because it copies all
 */
void RebuildTreeFromList(TreeHndl tree){
 /* Builds or rebuilds a tree that already has a tree->pointlist.
  */
	assert(tree != NULL);
	assert(tree->pointlist);
	if(tree->pointlist->Npoints < 1) return;

	unsigned long Npoints = tree->pointlist->Npoints,i,N=0;
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

	// fill point list with new (copied) points
	for(i=0;i<Npoints;++i){
		InsertPointAfterCurrent(tree->pointlist,&(points[i]));
	    MoveDownList(tree->pointlist);
	}
	tree->top->npoints = Npoints;
	tree->top->points = points;

	// build the tree
	 _BuildTree(tree);

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
  long i,j,count;

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
  double *x,xcut;

  cbranch=tree->current; /* pointer to current branch */

    /* leaf case */
  if(cbranch->npoints <= tree->Nbucket){
	  tree->current->points->leaf=tree->current;
	  return;
  }

  x=(double *)malloc(cbranch->npoints*sizeof(double));
  assert(x);

  /* initialize bounderies to old bounderies */
  for(i=0;i<2;++i){
      branch1.boundary_p1[i]=cbranch->boundary_p1[i];
      branch1.boundary_p2[i]=cbranch->boundary_p2[i];

      branch2.boundary_p1[i]=cbranch->boundary_p1[i];
      branch2.boundary_p2[i]=cbranch->boundary_p2[i];
  }

  /* set dimension to cut box */
  dimension=(cbranch->level % 2);

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

  if(median_cut){
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
  if(cut < cbranch->npoints) branch2.points=tree->pointlist->current;
  else branch2.points=NULL;

  free(x);

  /* centers of mass */

  for(i=0;i<2;++i) branch1.center[i]=0;
  tree->pointlist->current=branch1.points;
  for(i=0;i<cut; ++i){
/*    branch1.center[0]+=xp[points[i]][0]/branch1.npoints; */
/*    branch1.center[1]+=xp[points[i]][1]/branch1.npoints; */
	  branch1.center[0]+=tree->pointlist->current->x[0]/branch1.npoints;
	  branch1.center[1]+=tree->pointlist->current->x[1]/branch1.npoints;
	  MoveDownList(tree->pointlist);
  }

  for(i=0;i<2;++i) branch2.center[i]=0;
  tree->pointlist->current=branch2.points;
  for(i=cut;i<cbranch->npoints; ++i){
/*    branch2.center[0]+=xp[points[i]][0]/branch2.npoints; */
/*    branch2.center[1]+=xp[points[i]][1]/branch2.npoints; */
	  branch2.center[0]+=tree->pointlist->current->x[0]/branch2.npoints;
	  branch2.center[1]+=tree->pointlist->current->x[1]/branch2.npoints;
	  MoveDownList(tree->pointlist);
  }

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
  unsigned long i,j,cut,dimension;
  Branch branch1,branch2,*parent_branch;
  double *x,xcut;
  void _FindLeaf(TreeHndl tree,double *ray,unsigned long Nadd);
  Point *oldfirstpoint,*newfirstpoint;

  //checkTree(tree);

  if(Nadd==0) return 1;

   moveTop(tree);
   _FindLeaf(tree,xpoint->x,Nadd);
   parent_branch=tree->current;
   tree->current->npoints -= Nadd;  // subtract number that was added in _FindLeaf
   x=(double *) malloc(2*tree->Nbucket*sizeof(double));
   assert(x);

    for(j=0;j<Nadd;++j){

    	// add only that are inside original grid
    	if( inbox(xpoint[j].x,tree->top->boundary_p1,tree->top->boundary_p2) == 0 ){
    		ERROR_MESSAGE();
    		std::printf("ERROR: in AddPointToTree, ray is not inside the simulation box x = %e %e Nadd=%li\n  not adding it to tree\n",
    				   xpoint[j].x[0],xpoint[j].x[1],Nadd);
    		std::printf("root of tree\n");
       		printBranch(tree->top);
        		//exit(0);
    		//return 0;
    	}else{
    		tree->current = parent_branch;

    		if(inbox(xpoint[j].x,tree->current->boundary_p1,tree->current->boundary_p2)){
    			_FindLeaf(tree,xpoint[j].x,1);
    		}else{
    			//std::printf("going to other parent box\n");
    			while(inbox(xpoint[j].x,tree->current->boundary_p1,tree->current->boundary_p2)
    					== false){
    				if(atTop(tree)){ERROR_MESSAGE(); std::printf("ERROR: AddPointsToTree, point not in region\n   x=%e %e\n"
    						,xpoint[j].x[0],xpoint[j].x[1]); printBranch(tree->current); exit(1);}
    				moveUp(tree);
    				tree->current->npoints += j - Nadd;
    			}
    			_FindLeaf(tree,xpoint[j].x,Nadd-j);
    			tree->current->npoints += 1+j-Nadd;
    			parent_branch=tree->current;
    		}

			/***************** test lines  **************************

    		 if(tree->current->number == 4829221){
     			 std::printf("STOP!\n");
     		 }

			/***************** **************************/

    		if( !(atLeaf(tree)) ){
    			ERROR_MESSAGE();
    			std::printf("ERROR: _FindLeaf did not find a leaf for x = %e %e\n"
    					,xpoint[j].x[0],xpoint[j].x[1]);
    			printBranch(tree->current);
    			std::printf("\nchildren\n");
    			printBranch(tree->current->child1);
    			std::printf(" pointer = %p %e\n",&(tree->current->child1->boundary_p1[1])
    				,tree->current->child1->boundary_p1[1]);
    			printBranch(tree->current->child2);
    		}

    		// insert point into point list
    		if(tree->current->points == NULL){
    			// case of no previous points in leaf
    			assert(tree->current->npoints == 1);
    			tree->current->points = &xpoint[j];
     			tree->pointlist->current = tree->current->prev->points;
       			JumpDownList(tree->pointlist,tree->current->prev->npoints-2);  // adds point to end of parent branches list
       			InsertPointAfterCurrent(tree->pointlist,tree->current->points);
       			tree->pointlist->current = tree->current->points;
       			tree->current->points->leaf = tree->current;

       			assert(tree->current->points->next || tree->current->points->prev);
    		}else{
    			// case of point already in leaf
    			assert(atLeaf(tree));
    			tree->pointlist->current = tree->current->points;
    			JumpDownList(tree->pointlist,tree->current->npoints-2);  // adds point to end of branches list
    			InsertPointAfterCurrent(tree->pointlist,&xpoint[j]);
    			tree->pointlist->current = tree->current->points;
    			xpoint[j].leaf = tree->current;

    			assert(tree->current->points->next || tree->current->points->prev);
    		}

    		assert(tree->Nbucket == 1);
    		if( tree->current->npoints > tree->Nbucket ){ // create new leaves

    			// initialize boundaries to old boundaries
    			for(i=0;i<2;++i){
    				branch1.boundary_p1[i]=tree->current->boundary_p1[i];
    				branch1.boundary_p2[i]=tree->current->boundary_p2[i];

    				branch2.boundary_p1[i]=tree->current->boundary_p1[i];
    				branch2.boundary_p2[i]=tree->current->boundary_p2[i];
    			}

    			/* set dimension to cut box */
    			dimension=(tree->current->level % 2);

    			/* reorder points */
    			tree->pointlist->current = tree->current->points;
    			assert(tree->current->points->next || tree->current->points->prev);
    			for(i=0;i<tree->current->npoints;i++){
    				x[i]=tree->pointlist->current->x[dimension];
    				MoveDownList(tree->pointlist);
    			}

    			oldfirstpoint=tree->current->points;
    			/***************** test line  **************************/
    			assert(tree->current->points->next || tree->current->points->prev);

    			tree->current->points=sortList(tree->current->npoints,x,tree->pointlist,tree->current->points);

    			/***************** test line  **************************/
    			assert(tree->current->points->next || tree->current->points->prev);

    			newfirstpoint=tree->current->points;
    			/***************** test line  **************************/
    			assert(tree->current->points->next || tree->current->points->prev);

    			cut=tree->current->npoints/2;

    			// check that median split in this dimension will split particles
    			if(x[cut] == x[cut-1]){
    				// change dimension

    				dimension=!dimension;

    				tree->pointlist->current = tree->current->points;

    				/***************** test line  **************************/
        			assert(tree->current->points->next || tree->current->points->prev);

    				for(i=0;i<tree->current->npoints;i++){
    					x[i]=tree->pointlist->current->x[dimension];
    					MoveDownList(tree->pointlist);
    				}

    				tree->current->points=sortList(tree->current->npoints,x,tree->pointlist,tree->current->points);
    				newfirstpoint=tree->current->points;

    				/***************** test line  **************************/
    				assert(tree->current->points->next || tree->current->points->prev);

    			}

    			//     std::printf("top of branch list after sort id=%i\n",tree->current->points->id);
    			//     PrintList(tree->pointlist);

    			/*
         	std::printf("\n\nafter sortList n= %i\n",tree->current->npoints);
			tree->pointlist->current=tree->current->points;
			for(i=0;i<tree->current->npoints;i++){
				std::printf("%i  %f %f  x=%f\n",tree->pointlist->current->id
    				,tree->pointlist->current->x[0],tree->pointlist->current->x[1],x[i]);
				MoveDownList(tree->pointlist);
			}
    			 */

    			if(median_cut){
    				cut=tree->current->npoints/2;

    				branch1.boundary_p2[dimension]=(x[cut]+x[cut-1])/2;
    				branch2.boundary_p1[dimension]=(x[cut]+x[cut-1])/2;

    			}else{
    				xcut=(tree->current->boundary_p1[dimension]+tree->current->boundary_p2[dimension])/2;
    				branch1.boundary_p2[dimension]=xcut;
    				branch2.boundary_p1[dimension]=xcut;

    				locateD(x-1,tree->current->npoints,xcut,&cut);
    			}

    			/* set point numbers and pointers to points */
    			branch1.npoints=cut;
    			branch1.points=tree->current->points;
    			/***************** test line  **************************/
    			assert(tree->current->points->next || tree->current->points->prev);


    			branch2.npoints=tree->current->npoints - cut;
    			tree->pointlist->current=tree->current->points;
    			/***************** test line  **************************/
    			assert(tree->current->points->next || tree->current->points->prev);
    			JumpDownList(tree->pointlist,cut);
    			if(cut < tree->current->npoints) branch2.points = tree->pointlist->current;
    			else branch2.points=NULL;

    			/* centers of mass */

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

    			attachChildrenToCurrent(tree,branch1,branch2);
    			//attachChildToCurrent(tree,branch1,1);
    			//attachChildToCurrent(tree,branch2,2);

    			/***************** test lines  **************************

  	 		 if(tree->current->child1->number == 4829221 || tree->current->child2->number == 4829221){
    	  			 std::printf("STOP!\n");
    	  		 }

    			assert(tree->current->points->next || tree->current->points->prev);
       			assert(tree->current->child1->points->next || tree->current->child1->points->prev);
       			assert(tree->current->child2->points->next || tree->current->child2->points->prev);
			/*****************   **************************/

    			tree->current->child1->points->leaf = tree->current->child1;
    			tree->current->child2->points->leaf = tree->current->child2;

       			/***************** test lines  **************************
  			assert(newfirstpoint->next || newfirstpoint->prev);
    			assert(oldfirstpoint->next || oldfirstpoint->prev);

    			assert(tree->current->points == newfirstpoint);
			/*****************   **************************/

	    			/*** reset first particles in parent branches ***/
    			if( tree->current->points != oldfirstpoint){
 /*   				moveUp(tree);
    				while( tree->current->points == oldfirstpoint ){
    					tree->current->points = newfirstpoint;
    	    			assert(tree->current->points->next || tree->current->points->prev);
    					if(tree->current != tree->top) moveUp(tree);
    					else break;
    				}
  */
    				while(moveUp(tree)){
    					/***************** test lines  **************************
    					if(tree->current->number == 4829221){
    						std::printf("STOP!\n");
    					}
    					/*****************  **************************/

    					if(tree->current->points == oldfirstpoint) tree->current->points = newfirstpoint;
    					/***************** test line  **************************/
    					assert(tree->current->points->next || tree->current->points->prev);
    				}
    			}
    		}
    	}

    }

	/***************** test lines  **************************
 	 moveTop(tree);
 	 do{ assert(tree->current->points->next || tree->current->points->prev); }while(TreeWalkStep(tree,true));
	/*****************  **************************/


   	free(x);

   	//checkTree(tree);

    return 1;
}

/** \ingroup  ImageFindingL2
 *  \brief THIS DOES NOT WORK YET!!!
 *
 *  Reduces the size of the tree by removing points and branches that are no longer needed.
 *
 *
*/
unsigned long PruneTrees(
		TreeHndl i_tree     /// image plane tree
		,TreeHndl s_tree    /// source plane tree
		,double resolution  /// Maximum size of a cell to be removed.
		,bool useSB   /// If true it will not remove any point that has flux in it.
		){
	static ListHndl trashlist;
	static short init=1;
	Point *points;
	long i,Ntmp,count = 0;
	double res,initres;
	bool go;

	if(init){ trashlist = NewList(); init=0; }

	 if(i_tree == NULL) return 0;
	 if(s_tree == NULL) return 0;

	 Ntmp = i_tree->pointlist->Npoints;

	 moveTop(i_tree);
	 initres = (i_tree->top->boundary_p2[0]-i_tree->top->boundary_p1[0]);
	 if(resolution > initres/3 || resolution <= 0.0) return 0;  // do not allow pruning up to the initial grid size

	 // walk tree
	 go = true;
	 do{
		 assert(i_tree->current->points->next || i_tree->current->points->prev);

		 res = (i_tree->current->boundary_p2[0]-i_tree->current->boundary_p1[0]);
		 if( (res <= resolution && CurrentIsSquareTree(i_tree) )
				 &&  i_tree->current->npoints % 9 == 0){

			 if(useSB){
				 go = true;
				 // Check if surface brightness of all points in cell are zero.
				 i_tree->pointlist->current = i_tree->current->points;
				 for(i=0; i < i_tree->current->npoints;++i,MoveDownList(i_tree->pointlist) ){
					 if(i_tree->pointlist->current->surface_brightness > 0 ){
						 go = false;
						 break;
					 }
				 }
			 }

			 // remove all lower branches and make current a leaf
			 if(go && i_tree->current->npoints > 1) count += FreeBranchesBelow(i_tree,s_tree,trashlist);
			 assert(i_tree->current->points->next || i_tree->current->points->prev);

		 }
	 }while(TreeWalkStep(i_tree,true));

	 // rebuild source tree from list.
	 RebuildTreeFromList(s_tree);

		/***************** test lines  **************************
	 i=0;
	 moveTop(i_tree);
	 do{
		 if(atLeaf(i_tree)){
			 assert(i_tree->current->npoints == 1);
			 assert(i_tree->current->points->leaf == i_tree->current);
			 assert(i_tree->current->points->next || i_tree->current->points->prev);
			 ++i;
		 }else{
			 assert(i_tree->current->npoints > 1);
			 assert(i_tree->current->points->next || i_tree->current->points->prev);
		 }
	 }while(TreeWalkStep(i_tree,true));
	 assert(i == i_tree->top->npoints);


	 MoveToTopList(i_tree->pointlist);
	 i=0;
	 do{ assert(i_tree->pointlist->current->leaf); ++i;}while(MoveDownList(i_tree->pointlist));
	 MoveToTopList(s_tree->pointlist);
	 i=0;
	 do{ assert(s_tree->pointlist->current->leaf); ++i;}while(MoveDownList(s_tree->pointlist));
		/*****************  **************************/

	 assert(count == (Ntmp - i_tree->pointlist->Npoints) );

	 // Trash collection
	 if(count > 10 && trashlist->Npoints > 10){
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

	 moveTop(i_tree);
	 moveTop(s_tree);
	 //std::printf("count = %li\n",count);
	 return count;
}

/** \ingroup ImageFindingL2
 *
 *  Frees all branches of the tree below the current branch in i_tree
 * if that branch is square.  If current branch is not square nothing will happen.
 *
 * On exit: The i_tree->current is back to the original current.  If it is
 *          square it will have no children and contain one point.  The source
 *          points and branches are not removed.
 */

unsigned long FreeBranchesBelow(TreeHndl i_tree,TreeHndl s_tree,ListHndl trashlist){

	if(!CurrentIsSquareTree(i_tree)) return 0;
	if(atLeaf(i_tree)) return 0;

	assert( i_tree !=NULL);
	assert( s_tree !=NULL);

	Branch *branch,*top;
	Point *point;
	unsigned long Ntmp,NtoRemove,i,count = 0,count2 = 0,count1;
	double center[2];

	_freeBranches_iter(s_tree);  // s_tree will no longer be valid on exit.  This is to make sure it isn't used later without a rebuild.

	assert(i_tree->current->npoints % 9 == 0);
	top = i_tree->current;
	TreeWalkStep(i_tree,true);

	while( (top->child1 != NULL) || (top->child2 != NULL) ){

		if(atLeaf(i_tree)){
			//assert(i_tree->current->points->image->leaf);
			//s_tree->current = i_tree->current->points->image->leaf;  // set s_tree to source of current image cell

			/***************** test line  **************************/
			assert(i_tree->current->points->next || i_tree->current->points->prev);
			RemoveLeafFromTree(i_tree,&Ntmp);
			/***************** test line  **************************/
			assert(i_tree->current->points->next || i_tree->current->points->prev);

			// in a square leaf cell take out extra points that have come up from below
			if( ( CurrentIsSquareTree(i_tree) && atLeaf(i_tree) )
					&& i_tree->current->npoints == 9){
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
						< pow(i_tree->pointlist->current->gridsize/4,2) ){

						++count1;
						// keep this central point
						i_tree->pointlist->current->gridsize = i_tree->current->boundary_p2[0]-i_tree->current->boundary_p1[0];
						i_tree->pointlist->current->image->gridsize = i_tree->pointlist->current->gridsize;
						i_tree->current->points = i_tree->pointlist->current;
						/***************** test line  **************************/
						assert(i_tree->current->points->next || i_tree->current->points->prev);

					}else{

						++count;  // count of total number of points removed
						++count2;

						// reduce the number of particles in all parent cells

						/* First take points out of source plane
						 *   This is tricky because they are not ordered into
						 *   square blocks with 9 points in each.
						 */
/*
 * Below it reduces the s_tree branches, but later it was decided that s_tree should be rebuilt from scratch.
 *
 * 						assert(i_tree->pointlist->current->image->leaf);
						// go up the tree to subtract point from branch counts
						branch = s_tree->current = i_tree->pointlist->current->image->leaf;
						assert(s_tree->current->npoints == 1);
						assert(atLeaf(s_tree));
						--(s_tree->current->npoints);
						s_tree->current->points = NULL;

						while(moveUp(s_tree)){
							--(s_tree->current->npoints);
							if(s_tree->current->npoints == 0){
								// If the branch now has no points in it remove its children.
								moveToChild(s_tree,1);
								RemoveLeafFromTree(s_tree,&Ntmp);

								moveToChild(s_tree,2);
								RemoveLeafFromTree(s_tree,&Ntmp);

								s_tree->current->points = NULL;

							}else{
							// reassign first point in parent branches if necessary
								if(s_tree->current->points == i_tree->pointlist->current->image){
									s_tree->current->points = s_tree->current->points->next;
								}
							}
						}
*/
						// Take point out of the source point list

						assert(i_tree->pointlist->current->image);
						s_tree->pointlist->current = i_tree->pointlist->current->image;
						point = TakeOutCurrent(s_tree->pointlist);
						point->leaf = NULL;  // set leaf to NULL to indicate that point is no longer in tree
						if(point->head) InsertPointAfterCurrent(trashlist,point);  // collect heads for later trash collection


						// take points out of image plane
						branch = i_tree->current;
						do{
							assert(i_tree->current->npoints);
							--(i_tree->current->npoints);
							if(i_tree->current->points == i_tree->pointlist->current){
								assert(i_tree->current->points->next);
								i_tree->current->points = i_tree->current->points->next;
							}
							assert(i_tree->current->points->next || i_tree->current->points->prev);
						}while(moveUp(i_tree));
						i_tree->current = branch;

						if(i_tree->pointlist->current == i_tree->current->points) i_tree->current->points = i_tree->pointlist->current->next;
						point = TakeOutCurrent(i_tree->pointlist);
						point->leaf = NULL;
						// If point is a head of a memory block add it to trashlist for eventual trash collection
						if(point->head) InsertPointAfterCurrent(trashlist,point);
					}

				}
				assert(count1 == 1);
				assert(count2 == 8);
			}

		}
		assert(i_tree->current->points);
		if( !(atLeaf(i_tree)) )TreeWalkStep(i_tree,true);
	}

	i_tree->current = top;

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
	/***************** test line  **************************/
	assert(tree->current->points->next || tree->current->points->prev);

	if( atLeaf(tree) ){
		branch = tree->current;
		moveUp(tree);
		/***************** test line  **************************/
		assert(tree->current->points->next || tree->current->points->prev);

		if(branch == branch->prev->child1){
			branch->prev->child1 = NULL;
		}

		if(branch == branch->prev->child2){
			if(branch->prev->child1 != NULL) branch->prev->child1->brother = branch->prev->brother;
			branch->prev->child2 = NULL;
		}

		// leaves of points in the father
		tree->pointlist->current = branch->points;
		for(i=0;i<branch->npoints;++i,MoveDownList(tree->pointlist)) tree->pointlist->current->leaf = branch->prev;

		point = branch->points;
		*Npoints = branch->npoints;
		free(branch);
		--(tree->Nbranches);
		/***************** test line  **************************/
		 assert(tree->current->points->next || tree->current->points->prev);

		return point;
	}

	return NULL;
}


