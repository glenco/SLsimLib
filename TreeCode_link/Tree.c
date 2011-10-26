/*
 * Code Name:     tree.c                                       
 * Programmer:    R Ben Metcalf
 * Last Revised:  Nov, 2005                                   
 * Discription:  2 way tree data structure with 2 branches
 * Comments:                           
 */

/*#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "Tree.h"
#include <Kist.h>
#include <List.h>
#include <tree_maintenance.h>*/

#include <slsimlib.h>

/***** Structs *****/

/* Branch: Private struct, not exported */ 


/************************************************************************
 * NewBranch
 * Returns pointer to new Branch struct.  Initializes children pointers to NULL,
 * and sets data field to input.  Private.
 ************************************************************************/
/** \ingroup ConstructorL2
 *
 * Gives each branch a unique number even if branches are destroyed.
 */
Branch *NewBranch(Point *points,unsigned long npoints
		  ,double boundary_p1[2],double boundary_p2[2]
		  ,double center[2],int level){

    Branch *branch;
    int i;
    static unsigned long number = 0;

    branch = (Branch *)malloc(sizeof(Branch));
    if (!branch){
      ERROR_MESSAGE(); fprintf(stderr,"allocation failure in NewBranch()\n");
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

    return branch;
}

/** \ingroup ConstructorL2
*/
void FreeBranch(Branch *branch){

    assert( branch != NULL);
    free(branch);

    return;
}

  /** \ingroup ConstructorL2
 **/
Point *NewPointArray(
		unsigned long N  /// number of points in array
		,bool NewXs   /// Allocate memory for point positions or assign pointer to existing position
		){
  Point *points;
  unsigned long i;

  if(N <= 0) return NULL;
  points = (Point *) calloc(N,sizeof(Point));
  if(NewXs) points[0].x = (double *) calloc(2,sizeof(double));
  points[0].head = N;
  points[0].in_image = false;
  points[0].leaf = NULL;
  for(i=1;i<N;++i){
	  if(NewXs) points[i].x = (double *) calloc(2,sizeof(double));
	  points[i].head = 0;
	  points[i].in_image = false;
	  points[i].surface_brightness = 0;
	  points[i].leaf = NULL;
  }

  return points;
}

/** \ingroup ConstructorL2
 *
 */
void FreePointArray(Point *array){
  /* Note: this deallocates positions!! */
  unsigned long i;

  if(array[0].head){
	  for(i=0;i<array[0].head;++i) free(array[i].x);
	  free(array);
  }else{
	  ERROR_MESSAGE();
	  printf("ERROR: FreePointArray, miss aligned attempt to free point array");
	  exit(1);
  }
}

/** \ingroup ConstructorL2
 *  \brief  Make a new tree and the linked list of points in it.  Does
 *  not build the tree structure.  BuildTree() calls this and should be used
 *  for building trees.
 */
TreeHndl NewTree(
		Point *xp   /// array of points to be added to the tree
		,unsigned long npoints   /// number of points
		,double boundary_p1[2]   /// bottom left hand corner of root
		,double boundary_p2[2]   /// upper right hand corner of root
		,double center[2]        /// center of root (this could be the center of mass)
		,int Nbucket             /// maximum number of points allowed in a leaf
		){
  unsigned long i;
  TreeStruct *tree;

  tree = (TreeStruct *)malloc(sizeof(TreeStruct));
  if (!tree){
    ERROR_MESSAGE(); fprintf(stderr,"allocation failure in NewTree()\n");
    exit(1);
  }

    /* make linked list of points */
  tree->pointlist=NewList();
   //EmptyList(tree->pointlist);
  for(i=0;i<npoints;++i){
    InsertPointAfterCurrent(tree->pointlist,&xp[i]);
    MoveDownList(tree->pointlist);
  }

  tree->top=NewBranch(tree->pointlist->top,npoints,boundary_p1,boundary_p2
		      ,center,0);

  tree->Nbranches = 1;
  tree->current = tree->top;

  tree->Nbucket = Nbucket;
  //printTree(tree);
  /*FillList(tree->pointlist,&xp[1],npoints-1,1);*/
  return(tree);
}

/** \ingroup ConstructorL2
 * \brief Free tree and the linked list of points in it.
 */

short freeTree(TreeHndl tree){

	emptyTree(tree);
	free(tree->current);
	--tree->Nbranches;

	free(tree->pointlist);
	free(tree);

	tree = NULL;
	return 1;
}


/***** Access functions *****/

/************************************************************************
 * isEmpty
 * Returns "true" if the Tree is empty and "false" otherwise.  Exported.
 ************************************************************************/
bool isEmpty(TreeHndl tree){

    assert(tree != NULL);
    return(tree->Nbranches == 0);
}

/************************************************************************
 * atTop
 * Returns "true" if current is the same as top and "false" otherwise.
 * Exported.
 * Pre: !isEmpty(tree)
 ************************************************************************/
bool atTop(TreeHndl tree){

    assert(tree != NULL);
    if( isEmpty(tree) ){
	
	ERROR_MESSAGE(); fprintf(stderr, "Tree Error: calling atTop() on empty tree\n");
	exit(1);
    }
    return(tree->current == tree->top);
}

inline bool atLeaf(TreeHndl tree){
    return( (tree->current->child1==NULL)*(tree->current->child2==NULL) );
}

/************************************************************************
 * noChild
 * Returns "true" if the child of the current branch does not exist and "false" otherwise.
 * Exported.
 * Pre: !isEmpty(tree)
 ************************************************************************/
bool noChild(TreeHndl tree){

    assert(tree != NULL);
    if( isEmpty(tree) ){
	
	ERROR_MESSAGE(); fprintf(stderr, "Tree Error: calling atTop() on empty tree\n");
	exit(1);
    }

    if( (tree->current->child1 == NULL) || (tree->current->child2 == NULL) ) return true;
    return false;
}

/************************************************************************
 * offEnd
 * Returns "true" if current is off end and "false" otherwise.  Exported.
 ************************************************************************/
bool offEnd(TreeHndl tree){

    assert(tree != NULL);
    return(tree->current == NULL);
}

bool CurrentIsSquareTree(TreeHndl tree){
	if( fabs(1 - (tree->current->boundary_p2[0] - tree->current->boundary_p1[0])
			    /(tree->current->boundary_p2[1] - tree->current->boundary_p1[1]) )
			< 0.001){
		return true;
	}

	return false;
}

/************************************************************************
 * getCurrent
 * Returns the points of current.  Exported.
 * Pre: !offEnd(tree)
 ************************************************************************/
void getCurrent(TreeHndl tree,Point *points,unsigned long *npoints){

    assert(tree != NULL);
    if( offEnd(tree) ){
	
	ERROR_MESSAGE(); fprintf(stderr, "Tree Error: calling getCurrent() when current is off end\n");
	exit(1);
    }

    *npoints=tree->current->npoints;
    points=tree->current->points;

    return;
}

/************************************************************************
 * getNbranches
 * Returns the Nbranches of tree.  Exported.
 ************************************************************************/
unsigned long getNbranches(TreeHndl tree){

    assert(tree != NULL);
    return(tree->Nbranches);
}

/***** Manipulation procedures *****/

/************************************************************************
 * moveTop
 * Moves current to the front of tree.  Exported.
 * Pre: !isEmpty(tree)
 ************************************************************************/
void moveTop(TreeHndl tree){
    
	//assert(tree != NULL);
    if( isEmpty(tree) ){
    	ERROR_MESSAGE(); fprintf(stderr, "Tree Error: calling moveTop() on empty tree\n");
    	exit(1);
    }

    tree->current = tree->top;
    return;
}

/************************************************************************
 * movePrev
 * Moves current to the branch before it in tree.  This can move current
 * off end.  Exported.
 * Pre: !offEnd(tree)
 ************************************************************************/
bool moveToChild(TreeHndl tree,int child){
    
    assert(tree != NULL);
    assert(tree->current != NULL);

    if(child==1){
      if( tree->current->child1 == NULL ) return false;
      tree->current = tree->current->child1;
      return true;
    }
    if(child==2){
      if( tree->current->child2 == NULL ) return false;
      tree->current = tree->current->child2;
      return true;
    }
    return false;
}

bool moveUp(TreeHndl tree){

    assert(tree != NULL);
    if( offEnd(tree) ){
      ERROR_MESSAGE(); fprintf(stderr, "Tree Error: call to moveUp() when current is off end\n");
      exit(1);
    }

    if( tree->current == tree->top ) return false;
    assert(tree->current->prev);
    tree->current = tree->current->prev;  /* can move off end */
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
void insertChildToCurrent(TreeHndl tree,Point *points,unsigned long npoints
			  ,double boundary_p1[2],double boundary_p2[2]
			  ,double center[2],int child){
    
    Branch *branch;

    /*printf("attaching child%i  current paricle number %i\n",child,tree->current->npoints);*/

    branch = NewBranch(points,npoints,boundary_p1,boundary_p2,center
		       ,tree->current->level+1);

    assert(tree != NULL);
    
    if( offEnd(tree) ){
    	ERROR_MESSAGE(); fprintf(stderr, "Tree Error: calling insertChildToCurrent() when current is off end\n");
    	exit(1);
    }

    branch->prev = tree->current;

    if(child==1){
      if(tree->current->child1 != NULL){
    	  ERROR_MESSAGE(); fprintf(stderr, "Tree Error: calling insertChildToCurrent() when child1 already exists\n");
    	  exit(1);
      }
      tree->current->child1 = branch;
      tree->current->child1->brother = tree->current->child2;
    }
    if(child==2){
      if(tree->current->child2 != NULL){
    	  ERROR_MESSAGE();
    	  fprintf(stderr, "Tree Error: calling insertChildToCurrent() when child2 already exists\n  current level=%i Nbranches=%li\n",tree->current->level,tree->Nbranches);
    	  exit(1);
      }
      tree->current->child2 = branch;      
      tree->current->child2->brother = tree->current->brother;
    }

    tree->Nbranches++;

    return;
}

  /* same as above but takes a branch structure */

void attachChildToCurrent(TreeHndl tree,Branch data,int child){

  insertChildToCurrent(tree,data.points,data.npoints,data.boundary_p1,data.boundary_p2,data.center,child);
  return;
}

void attachChildrenToCurrent(TreeHndl tree,Branch child1,Branch child2){
	// this is an addition that keeps assigns the brother pointers

	insertChildToCurrent(tree,child1.points,child1.npoints
			,child1.boundary_p1,child1.boundary_p2,child1.center,1);
	insertChildToCurrent(tree,child2.points,child2.npoints
			,child2.boundary_p1,child2.boundary_p2,child2.center,2);

	tree->current->child1->brother = tree->current->child2;
	tree->current->child2->brother = tree->current->brother;
	return;
}

/***** Other operations *****/


/************************************************************************
 * printTree
 * Prints the contents of tree to stdout.  The current element, if there
 * is one, will be enclosed in []'s.  Currently, only to be used when the
 * TreeElements is are integers.  Exported.
 ************************************************************************/  
void printTree(TreeHndl tree){
  int i;
    assert( tree !=NULL);

    printBranch(tree->current);
    tree->pointlist->current=tree->current->points;
    for(i=0;i<tree->current->npoints;++i){
      printf("    %li    %f %f\n",tree->pointlist->current->id,tree->pointlist->current->x[0],tree->pointlist->current->x[1]);
      MoveDownList(tree->pointlist);
    }
    if(tree->current->child1 == NULL) return;

    if( (tree->current->boundary_p1[0]==tree->current->boundary_p2[0]) ||
    		(tree->current->boundary_p1[0]==tree->current->boundary_p2[0])	){
    	ERROR_MESSAGE();
    	printf("ERROR: zero area branch");
    	exit(0);
    }
    moveToChild(tree,1);
    printTree(tree);

    moveUp(tree);

    moveToChild(tree,2);
    printTree(tree);

    moveUp(tree);

    return;
}

void printBranch(Branch *data){

  printf("******* branch *******\nlevel=%i number=%li\n",data->level,data->number);
  printf("center = [%e,%e]\n",data->center[0],data->center[1]);
  printf("p1 = [%e,%e] p2 = [%e,%e]\n"
	 ,data->boundary_p1[0],data->boundary_p1[1]
	 ,data->boundary_p2[0],data->boundary_p2[1]);
  printf("number of points = %li\n",data->npoints);
  /*for(i=0;i<data->npoints;++i) printf("%e %e %e\n",xp[data->points[i]][0] ,xp[data->points[i]][1],xp[data->points[i]][2]); */
}
/****************************************************************
 *  code for testing tree
 *****************************************************************/

void checkTree(TreeHndl tree){
	void _checkTree(TreeHndl tree,unsigned long *count);
	Branch *initial;
	unsigned long count=0;

	initial=tree->current;

	moveTop(tree);
	_checkTree(tree,&count);

	if(count != tree->Nbranches){ printf("checkTree did not reach all branches.\n"); exit(0);}
	tree->current=initial;
	return;
}

void _checkTree(TreeHndl tree,unsigned long *count){
	int checkBranch(Branch *branch);

	//printf("     hello\n");
	++*count;
	if(checkBranch(tree->current)) exit(1);

	if(tree->current->child1 != NULL){
		moveToChild(tree,1);
		_checkTree(tree,count);
		moveUp(tree);
	}

	if(tree->current->child2 != NULL){
		moveToChild(tree,2);
		_checkTree(tree,count);
		moveUp(tree);
	}

    return;
}

// put in test for each branch
int checkBranch(Branch *branch){
   	if(branch->boundary_p1[0] >= branch->boundary_p2[0]
   	  || branch->boundary_p1[1] >= branch->boundary_p2[1] ){
    		printf("\nthis branch is screwed\n");
    		printBranch(branch);
    		PrintPoint(branch->points);
    		return 1;
   	}
   	return 0;
}

Point *AddPointToArray(Point *points,unsigned long N,unsigned long Nold){
  unsigned long i;

  if(Nold==0){
	  points = NewPointArray(N,true);
  }else{
	  if(points[0].head != Nold){ ERROR_MESSAGE(); printf("ERROR: AddPointToArray head not set correctly\n"); exit(0);}
	  for(i=N;i<Nold;++i) free(points[i].x);
	  points=(Point *) realloc(points,N*sizeof(Point));
	  for(i=Nold;i<N;++i){
		  points[i].x=(double *) malloc(2*sizeof(double));
		  assert(points[i].x);
		  points[i].head=0;
		  points[i].in_image=false;
		  points[i].leaf=NULL;
	  }
	  if(N>0) points[0].head=N;
  }

  return points;
}

Point *NewPoint(double *x,unsigned long id){
  Point *point;

  point=(Point *) malloc(sizeof(Point));
  assert(point);

  point->head = 1;
  point->id=id;
  point->x=x;
  point->in_image=false;

  if (!point){
    ERROR_MESSAGE(); fprintf(stderr,"allocation failure in NewPoint()\n");
    exit(1);
  }
  point->next = point->prev = NULL;
  point->leaf=NULL;

  return(point);
}

void SwapPointsInArray(Point *p1,Point *p2){
  /* SWAPS information in points without changing */
  /* pointers to prev and next */
  /* chnages links of image points to fallow*/ 
  Point pt;

  if(p1==p2) return;

  PointCopy(&pt,p1);
  PointCopy(p1,p2);
  PointCopy(p2,&pt);
}

void PointCopy(Point *pcopy,Point *pin){
  /* copies information in point without copying
  * pointers to prev and next, but moving the link
  * to the image point */
  pcopy->id = pin->id;
  pcopy->image = pin->image;
  pcopy->invmag = pin->invmag;
  pcopy->kappa = pin->kappa;
  pcopy->gamma[0] = pin->gamma[0];
  pcopy->gamma[1] = pin->gamma[1];
  pcopy->dt = pin->dt;
  pcopy->x = pin->x;
  pcopy->gridsize = pin->gridsize;
  pcopy->in_image = pin->in_image;
  pcopy->surface_brightness = pin->surface_brightness;
  pcopy->leaf = pin->leaf;

  if((pin->image != NULL) && (pin->image->image == pin)) pin->image->image = pcopy;
}
void PointCopyData(Point *pcopy,Point *pin){
  /* copies information in point without copying */
  /* pointers to prev and next */
  /* does copy image pointer */
  pcopy->id = pin->id;
  pcopy->image = pin->image;
  pcopy->invmag = pin->invmag;
  pcopy->kappa = pin->kappa;
  pcopy->gamma[0] = pin->gamma[0];
  pcopy->gamma[1] = pin->gamma[1];
  pcopy->dt = pin->dt;
  pcopy->x = pin->x;
//  pcopy->x[0] = pin->x[0];
//  pcopy->x[1] = pin->x[1];
  pcopy->gridsize = pin->gridsize;
  pcopy->in_image = pin->in_image;
  pcopy->surface_brightness = pin->surface_brightness;
  pcopy->leaf = pin->leaf;
}

/*********************************/
/*  point extraction routines */
/*********************************/
void PointsInCurrent(TreeHndl tree,unsigned long *ids,double **x){
  Point *point;
  unsigned long i;

  tree->pointlist->current=tree->current->points;

  point=tree->current->points;
  for(i=0;i<tree->current->npoints;++i){
    x[i]=tree->pointlist->current->x;
    ids[i]=tree->pointlist->current->id;
    MoveDownList(tree->pointlist);
  }

  return;
}

void PrintPoint(Point *point){
  printf("Point id = %li\n",point->id);
  printf("   x= %e %e gridsize = %e \n",point->x[0],point->x[1]
	 ,point->gridsize);
  printf("   kappa= %e gamma = %e %e invmag =%e\n",point->kappa
	 ,point->gamma[0],point->gamma[1],point->invmag);
  printf("   in_image=%i\n",point->in_image);
}
/** \ingroup Constructor
 * Make an array of imageinfo types.
 */
ImageInfo *NewImageInfo(int Nimages){
  ImageInfo *imageinfo;
  int i;

  imageinfo=(ImageInfo *) malloc(Nimages*sizeof(ImageInfo));
  assert(imageinfo);

  for(i=0;i<Nimages;++i){
	imageinfo[i].imagekist = NewKist();
	imageinfo[i].Npoints=0;
    imageinfo[i].innerborder = NewKist();
    imageinfo[i].outerborder = NewKist();
  }

  return imageinfo;
}

/** \ingroup Constructor
 * Free an array of imageinfo types.
 */
void freeImageInfo(ImageInfo *imageinfo,int Nimages){
	int i;

	//freeKist(imageinfo->imagekist);
	for(i=0;i<Nimages;++i){
		freeKist(imageinfo[i].imagekist);
	    freeKist(imageinfo[i].innerborder);
	    freeKist(imageinfo[i].outerborder);
	}
	free(imageinfo);
}

/** \ingroup LowLevel
 *  step for walking tree by iteration instead of recursion
 */
bool TreeWalkStep(TreeHndl tree,bool allowDescent){

	if(allowDescent && tree->current->child1 != NULL){
		moveToChild(tree,1);
		return true;
	}
	if(allowDescent && tree->current->child2 != NULL){
		moveToChild(tree,2);
		return true;
	}

	if(tree->current->brother != NULL){
		tree->current = tree->current->brother;
		return true;
	}

	return false;
}


/**************************************************************************
   routines for saving and restoring tree structure
*************************************************************************/
/*
void saveTree(TreeHndl tree,char *filename){
  RelativeBranch *tree_arr;
  double *x,*y;
  unsigned long *point_id,i;
  void _saveTree(TreeHndl tree,RelativeBranch *tree_arr);
  FILE *file;

  tree_arr=(RelativeBranch *)malloc(tree->Nbranches*sizeof(RelativeBranch));
  x=(double *)malloc(tree->pointlist->Npoints*sizeof(double));
  y=(double *)malloc(tree->pointlist->Npoints*sizeof(double));
  point_id=(unsigned long *)malloc(tree->pointlist->Npoints*sizeof(unsigned long));

  // write linked list to arrays
  MoveToTopList(tree->pointlist);
  for(i=0;i<tree->pointlist->Npoints;++i){
    x[i]=tree->pointlist->current->x[0];
    y[i]=tree->pointlist->current->x[1];
    point_id[i]=tree->pointlist->current->id;
    MoveDownList(tree->pointlist);
  }

  file=fopen(filename,"w");
  fwrite(&(tree->top->npoints),sizeof(unsigned long),1,file);
  fwrite(&(tree->Nbranches),sizeof(unsigned long),1,file);
  fwrite(point_id,sizeof(unsigned long),tree->top->npoints,file);
  fwrite(x,sizeof(double),tree->top->npoints,file);
  fwrite(y,sizeof(double),tree->top->npoints,file);
   fclose(file);

  // free(tree_arr);
  free(x);
  free(y);
  free(point_id);

  return;
}*/

/* void _saveTree(TreeHndl tree,RelativeBranch *tree_arr){ */
/*   int i; */
/*   unsigned long current; */

/*   current=tree->current->number; */

/*   /\* count down list to branch's first point *\/ */
/*   MoveToTopList(tree->pointlist); */
/*   tree_arr[current].points=0; */
/*   while(tree->pointlist->current != tree->current->points){ */
/*     MoveDownList(tree->pointlist); */
/*     ++tree_arr[current].points; */
/*   } */

/*   tree_arr[current].npoints=tree->current->npoints; */
/*   for(i=0;i<2;++i){ */
/*     tree_arr[current].center[i]=tree->current->center[i]; */
/*     tree_arr[current].boundary_p1[i]=tree->current->boundary_p1[i]; */
/*     tree_arr[current].boundary_p2[i]=tree->current->boundary_p2[i]; */
/*   } */

/*   if(current>0) tree_arr[current].prev=tree->current->prev->number; */

/*   if(tree->current->child1 == NULL){ */
/*     tree_arr[current].child1=-1; */
/*     return; */
/*   }else{ */
/*     tree_arr[current].child1=tree->current->child1->number; */
/*     moveToChild(tree,1); */
/*     _saveTree(tree,tree_arr); */
/*     moveUp(tree); */
/*   } */

/*   if(tree->current->child2 == NULL){ */
/*     tree_arr[current].child2=-1; */
/*     return; */
/*   }else{ */
/*     tree_arr[current].child2=tree->current->child2->number; */
/*     moveToChild(tree,2); */
/*     _saveTree(tree,tree_arr); */
/*     moveUp(tree); */
/*   } */

/*   return; */
/* } */
/*
TreeHndl readTree(char *filename){
  RelativeBranch *tree_arr;
  TreeHndl tree;
  FILE *file;
  double *x,*y;
  Point *xp;
  unsigned long Npoints,Nbranches,*point_id,i;
  void _readTree(TreeHndl tree,RelativeBranch *tree_arr,unsigned long *points
		 ,unsigned long current);

  file=fopen(filename,"r");

  fread(&Npoints,sizeof(unsigned long),1,file);
  fread(&Nbranches,sizeof(unsigned long),1,file);
  printf("npoints %i   Nbranches %i\n",Npoints,Nbranches);

  x=(double *)malloc(Npoints*sizeof(double));
  y=(double *)malloc(Npoints*sizeof(double));
  point_id=(unsigned long *)malloc(Npoints*sizeof(unsigned long));

  fread(point_id,sizeof(unsigned long),Npoints,file);
  fread(x,sizeof(double),Npoints,file);
  fread(y,sizeof(double),Npoints,file);
  //fread(tree_arr,sizeof(RelativeBranch),Nbranches,file);
  fclose(file);

  xp=(Point *) malloc(Npoints*sizeof(Point));
  for(i=0;i<Npoints;++i){
    xp[i].x[0]=x[i];
    xp[i].x[1]=y[i];
  }
  free(x);
  free(y);


  // make linked list of points

  tree=BuildTree(xp,Npoints);

  /* re-assign particle ids */
/*   MoveToTopList(tree->pointlist); */
/*   for(i=0;i<tree->pointlist->Npoints;++i){ */
/*     printf("i=%i   %i\n",i,point_id[i]); */
/*      if( (tree->pointlist->current->x[0] == xp[i][0])*(tree->pointlist->current->x[1] == xp[i][1]) ){ */
/*       tree->pointlist->current->id=point_id[i]; */
/*      }else{ */
/*       ERROR_MESSAGE();
/*       printf("ERROR in readtree: reconstructed tree is not the same as saved one\n"); */
/*       exit(1); */
/*     } */
/*     MoveDownList(tree->pointlist); */
/*   } */

//  return tree;
//}
/*
void _readTree(TreeHndl tree,RelativeBranch *tree_arr,ListHndl points
	       ,unsigned long current){
  int i;

  tree->current->number=current;
  tree->current->points=&points[tree_arr[current].points];
  tree->current->npoints=tree_arr[current].npoints;
  tree->current->level=tree_arr[current].level;

  for(i=0;i<2;++i){
    tree->current->center[i]=tree_arr[current].center[i];
    tree->current->boundary_p1[i]=tree_arr[current].boundary_p1[i];
    tree->current->boundary_p2[i]=tree_arr[current].boundary_p2[i];
  }

  if(tree_arr[current].child1==-1){
    tree->current->child1=NULL;
    return;
  }else{
    insertChildToCurrent(tree,&current,0,&dummy,&dummy,&dummy,1);
    moveToChild(tree,1);
    _readTree(tree,tree_arr,points,tree_arr[current].child1);
    moveUp(tree);
  }

  if(tree_arr[current].child2==-1){
    tree->current->child2=NULL;
    return;
  }else{
    insertChildToCurrent(tree,&current,0,&dummy,&dummy,&dummy,2);
    moveToChild(tree,2);
    _readTree(tree,tree_arr,points,tree_arr[current].child2);
    moveUp(tree);
  }

  return;
}
*/
