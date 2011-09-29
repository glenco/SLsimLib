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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <nrD.h>
#include "Tree.h"
#include <nrutil.h>

static int median_cut=1;

TreeHndl BuildTree(Point *xp,unsigned long Npoints){
  TreeHndl tree;
  unsigned long i;
  double p1[2],p2[2],center[2];
  void _BuildTree(TreeHndl tree);
  const int Nbucket = 1;   // must be =1 if each leaf is to coincide with each cell

  if( (Npoints & (Npoints-1)) != 0){
	  ERROR_MESSAGE();
	  printf("ERROR: BuildTree, Npoints is not a power of 2\n");
	  exit(1);
  }

  p1[0]=xp[0].x[0]; p1[1]=xp[0].x[1];
  p2[0]=xp[0].x[0]; p2[1]=xp[0].x[1];

  for(i=0;i<Npoints;++i){

    /* find X boundery */
    if(xp[i].x[0] < p1[0] ) p1[0]=xp[i].x[0];
    if(xp[i].x[0] > p2[0] ) p2[0]=xp[i].x[0];

    /* find Y boundery */
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

  _BuildTree(tree);

  return ;
}

/* tree must be created and first branch must be set before */
/* start */

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
      branch1.boundery_p1[i]=cbranch->boundery_p1[i];
      branch1.boundery_p2[i]=cbranch->boundery_p2[i];

      branch2.boundery_p1[i]=cbranch->boundery_p1[i];
      branch2.boundery_p2[i]=cbranch->boundery_p2[i];
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

	  cut=cbranch->npoints/2;
      branch1.boundery_p2[dimension]=(x[cut]+x[cut-1])/2;
      branch2.boundery_p1[dimension]=(x[cut]+x[cut-1])/2;

  }else{

	  xcut=(cbranch->boundery_p1[dimension]+cbranch->boundery_p2[dimension])/2;
      branch1.boundery_p2[dimension]=xcut;
      branch2.boundery_p1[dimension]=xcut;

	  quickPartitionPoints(xcut,&cut
	  		,tree->current->points,x,cbranch->npoints);

      //locateD(x-1,cbranch->npoints,xcut,&cut);
  }

  /* set point numbers and pointers to points */
  branch1.npoints=cut;
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

 /*printf("reached end of _BuildTree level=%i\n",tree->current->level);*/
 return;
}

int AddPointsToTree(TreeHndl tree,Point *xpoint,unsigned long Nadd){
	/****************************************************************
	 * Expands tree by adding points
	 ****************************************************************/
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
   tree->current->npoints -= Nadd;
   x=(double *) malloc(2*tree->Nbucket*sizeof(double));
   assert(x);

    for(j=0;j<Nadd;++j){

    	// add only that are inside original grid
    	if( inbox(xpoint[j].x,tree->top->boundery_p1,tree->top->boundery_p2) == 0 ){
    		ERROR_MESSAGE();
    		printf("ERROR: in AddPointToTree, ray is not inside the simulation box x = %e %e Nadd=%li\n  not adding it to tree\n",
    				   xpoint[j].x[0],xpoint[j].x[1],Nadd);
    		printf("root of tree\n");
       		printBranch(tree->top);
        		//exit(0);
    		//return 0;
    	}else{
    		tree->current=parent_branch;

    		if(inbox(xpoint[j].x,tree->current->boundery_p1,tree->current->boundery_p2)){
    			_FindLeaf(tree,xpoint[j].x,1);
    		}else{
    			//printf("going to other parent box\n");
    			while(inbox(xpoint[j].x,tree->current->boundery_p1,tree->current->boundery_p2)
    					== False){
    				if(atTop(tree)){ERROR_MESSAGE(); printf("ERROR: AddPointsToTree, point not in region\n   x=%e %e\n"
    						,xpoint[j].x[0],xpoint[j].x[1]); printBranch(tree->current); exit(1);}
    				moveUp(tree);
    				tree->current->npoints += j - Nadd;
    			}
    			_FindLeaf(tree,xpoint[j].x,Nadd-j);
    			tree->current->npoints += 1+j-Nadd;
    			parent_branch=tree->current;
    		}

    		if( tree->current->child1 != NULL || tree->current->child2 != NULL){
    			ERROR_MESSAGE();
    			printf("ERROR: _FindLeaf did not find a leaf for x = %e %e\n"
    					,xpoint[j].x[0],xpoint[j].x[1]);
    			printBranch(tree->current);
    			printf("\nchildren\n");
    			printBranch(tree->current->child1);
    			printf(" pointer = %p %e\n",&(tree->current->child1->boundery_p1[1])
    				,tree->current->child1->boundery_p1[1]);
    			printBranch(tree->current->child2);
    		}

    		// insert point into point list
    		tree->pointlist->current=tree->current->points;
    		JumpDownList(tree->pointlist,tree->current->npoints-2);  // adds point to end of branches list
    		InsertPointAfterCurrent(tree->pointlist,&xpoint[j]);
    		tree->pointlist->current=tree->current->points;

    		if( tree->current->npoints > tree->Nbucket ){ // create new leaves

    			// initialize boundaries to old boundaries
    			for(i=0;i<2;++i){
    				branch1.boundery_p1[i]=tree->current->boundery_p1[i];
    				branch1.boundery_p2[i]=tree->current->boundery_p2[i];

    				branch2.boundery_p1[i]=tree->current->boundery_p1[i];
    				branch2.boundery_p2[i]=tree->current->boundery_p2[i];
    			}

    			/* set dimension to cut box */
    			dimension=(tree->current->level % 2);

    			/* reorder points */
    			tree->pointlist->current=tree->current->points;
    			for(i=0;i<tree->current->npoints;i++){
    				x[i]=tree->pointlist->current->x[dimension];
    				MoveDownList(tree->pointlist);
    			}

    			oldfirstpoint=tree->current->points;
    			tree->current->points=sortList(tree->current->npoints,x,tree->pointlist,tree->current->points);
    			newfirstpoint=tree->current->points;

    			cut=tree->current->npoints/2;

    			// check that median split in this dimension will split particles
    			if(x[cut] == x[cut-1]){
    				// change dimension

    				dimension=!dimension;

    				tree->pointlist->current=tree->current->points;
    				for(i=0;i<tree->current->npoints;i++){
    					x[i]=tree->pointlist->current->x[dimension];
    					MoveDownList(tree->pointlist);
    				}

    				tree->current->points=sortList(tree->current->npoints,x,tree->pointlist,tree->current->points);
    				newfirstpoint=tree->current->points;
    			}

    			//     printf("top of branch list after sort id=%i\n",tree->current->points->id);
    			//     PrintList(tree->pointlist);

    			/*
         	printf("\n\nafter sortList n= %i\n",tree->current->npoints);
			tree->pointlist->current=tree->current->points;
			for(i=0;i<tree->current->npoints;i++){
				printf("%i  %f %f  x=%f\n",tree->pointlist->current->id
    				,tree->pointlist->current->x[0],tree->pointlist->current->x[1],x[i]);
				MoveDownList(tree->pointlist);
			}
    			 */

    			if(median_cut){
    				cut=tree->current->npoints/2;

    				branch1.boundery_p2[dimension]=(x[cut]+x[cut-1])/2;
    				branch2.boundery_p1[dimension]=(x[cut]+x[cut-1])/2;

    			}else{
    				xcut=(tree->current->boundery_p1[dimension]+tree->current->boundery_p2[dimension])/2;
    				branch1.boundery_p2[dimension]=xcut;
    				branch2.boundery_p1[dimension]=xcut;

    				locateD(x-1,tree->current->npoints,xcut,&cut);
    			}

    			/* set point numbers and pointers to points */
    			branch1.npoints=cut;
    			branch1.points=tree->current->points;

    			branch2.npoints=tree->current->npoints - cut;
    			tree->pointlist->current=tree->current->points;
    			JumpDownList(tree->pointlist,cut);
    			if(cut < tree->current->npoints) branch2.points=tree->pointlist->current;
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

    			tree->current->child1->points->leaf = tree->current->child1;
    			tree->current->child2->points->leaf = tree->current->child2;

    			/*** reset first particles in parent branches ***/
    			if( tree->current->points != oldfirstpoint){
    				moveUp(tree);
    				while( tree->current->points == oldfirstpoint ){
    					tree->current->points = newfirstpoint;
    					if(tree->current != tree->top) moveUp(tree);
    					else break;
    				}
    			}
    		}
    	}
    }

   	free(x);

   	//checkTree(tree);

    return 1;
}

unsigned long PruneTree(TreeHndl i_tree,TreeHndl s_tree,double resolution){
	static ListHndl trashlist;
	static short init=1;
	Point *points;
	long i,Ntmp,count = 0;
	double res,initres;

	some ability to take surface brightness into account

	if(init){ trashlist = NewList(); init=0; }

	 if(i_tree == NULL) return 0;
	 if(s_tree == NULL) return 0;

	 Ntmp = i_tree->pointlist->Npoints;

	 moveTop(i_tree);
	 initres = (i_tree->top->boundery_p2[0]-i_tree->top->boundery_p1[0]);
	 if(resolution > initres/3 || resolution <= 0.0) return 0;  // do not allow pruning up to the initial grid size

	 // walk tree
	 do{
		 res = (i_tree->current->boundery_p2[0]-i_tree->current->boundery_p1[0]);
		 if(res <= resolution ){
			 // remove all lower branches and make current a leaf
			 count += FreeBranchesBelow(i_tree,s_tree,trashlist);
		 }
	 }while(TreeWalkStep(i_tree,True));

	 assert(count == (Ntmp - i_tree->pointlist->Npoints) );

	 // Trash collection
	 if(count > 10 && trashlist->Npoints > 10){
		 Boolean step;
		 MoveToTopList(trashlist);
		 do{
			 // check to see if all points in the block have been removed from the trees
			 for(i=0;i<trashlist->current->head;++i) if(trashlist->current[i].leaf != NULL) break;
			 if(i == trashlist->current->head){
				 if(AtTopList(trashlist)) step = False; else step = True;
				 points = TakeOutCurrent(trashlist);
				 FreePointArray(points);
			 }
		 }while(MoveDownList(trashlist) && step);
	 }

	 return count;
}


unsigned long FreeBranchesBelow(TreeHndl i_tree,TreeHndl s_tree,ListHndl trashlist){
	/* Frees all branches of the tree below the current branch in i_tree
	 * if that branch is square.  If current branch is not square nothing will happen.
	 *
	 * On exit: The i_tree->current is back to the original current.  If it is
	 *          square it will have no children and contain one point.  The source
	 *          points and branches are also removed.
	 */

	if(!CurrentIsSquareTree(i_tree)) return 0;
	if(atLeaf(i_tree)) return 0;

	assert( i_tree !=NULL);
	assert( s_tree !=NULL);

	Branch *branch,*top;
	Point *point;
	unsigned long Ntmp,i,count = 0;

	top = i_tree->current;
	TreeWalkStep(i_tree,True);
	while(i_tree->current != top->brother){

		if(atLeaf(i_tree)){
			s_tree->current = i_tree->current->points->image->leaf;  // set s_tree to source of current image cell
			RemoveLeafFromTree(s_tree,&Ntmp);
			RemoveLeafFromTree(i_tree,&Ntmp);

			// in a square leaf cell take out extra points that have come up from below
			if( CurrentIsSquareTree(i_tree) && atLeaf(i_tree) ){

				i_tree->pointlist->current = i_tree->current->points;
				Ntmp = i_tree->current->npoints;
				for(i=0;i<Ntmp;++i,MoveDownList(i_tree->pointlist)){
					// find central point and remove others
					if( (pow(i_tree->current->center[0]-i_tree->pointlist->current->x[0],2)
						+ pow(i_tree->current->center[1]-i_tree->pointlist->current->x[1],2) )
						< pow(i_tree->pointlist->current->gridsize,2) ){

						// keep this central point
						i_tree->pointlist->current->gridsize = i_tree->current->boundery_p2[0]-i_tree->current->boundery_p1[0];
						i_tree->pointlist->current->image->gridsize = i_tree->pointlist->current->gridsize;
						i_tree->current->points = i_tree->pointlist->current;

					}else{

						++count;

						// reduce the number of particles in all parent cells
						branch = s_tree->current = i_tree->pointlist->current->image->leaf;
						do{
							--(s_tree->current->npoints);
							if(s_tree->current->points == i_tree->pointlist->current->image){
								s_tree->current->points = s_tree->current->points->next;
							}
						}while(moveUp(s_tree));
						s_tree->current = branch;

						// Take point out of the source point list
						s_tree->pointlist->current = i_tree->pointlist->current->image;
						point = TakeOutCurrent(s_tree->pointlist);
						point->leaf = NULL;  // set leaf to NULL to indicate that point is no longer in tree
						if(point->head) InsertPointAfterCurrent(trashlist,point);  // save the head of memory blocks

						branch = i_tree->current;
						do{
							--(i_tree->current->npoints);
							if(i_tree->current->points == i_tree->pointlist->current){
								i_tree->current->points = i_tree->current->points->next;
							}
						}while(moveUp(i_tree));
						i_tree->current = branch;

						point = TakeOutCurrent(i_tree->pointlist);
						point->leaf = NULL;
						// If point is a head of a memory block add it to trashlist for eventual trash collection
						if(point->head) InsertPointAfterCurrent(trashlist,point);
					}

				}
			}

		}
		TreeWalkStep(i_tree,True);
	}

	i_tree->current = top;

    return count;
}

Point *RemoveLeafFromTree(TreeHndl tree,unsigned long *Npoints){
	/*
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

	Branch *branch;
	Point *point;
	unsigned long i;

	if(atTop(tree)) return NULL;

	if( atLeaf(tree) ){
		branch = tree->current;

		if(branch == branch->prev->child1){
			branch->prev->child1 = NULL;
		}

		if(branch == branch->prev->child2){
			if(branch->prev->child1 != NULL) branch->prev->child1->brother = branch->prev->brother;
			branch->prev->child2 = NULL;
		}

		// leaves of points to father
		tree->pointlist->current = branch->points;
		for(i=0;i<branch->npoints;++i,MoveDownList(tree->pointlist)) tree->pointlist->current->leaf = branch->prev;
		moveUp(tree);

		point = branch->points;
		*Npoints = branch->npoints;
		free(branch);
		--tree->Nbranches;

		return point;
	}

	return NULL;
}
