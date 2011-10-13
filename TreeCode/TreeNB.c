/*
 * Code Name:     tree.c                                       
 * Programmer:    R Ben Metcalf
 * Last Revised:  Nov, 2005                                   
 * Discription:  2 way tree data structure with 2 branchNBes
 * Comments:                           
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "TreeNB.h"

double dummy;

/***** Structs *****/

/* BranchNB: Private struct, not exported */ 


/***** Constructors/Destructors *****/

/************************************************************************
 * NewBranchNB
 * Returns pointer to new BranchNB struct.  Initializes children pointers to NULL,
 * and sets data field to input.  Private.
 ************************************************************************/
BranchNB *NewBranchNB(IndexType *particles,IndexType nparticles
		  ,PosType boundary_p1[treeNBdim],PosType boundary_p2[treeNBdim]
		  ,PosType center[treeNBdim],int level,unsigned long branchNBnumber){

    BranchNB *branchNB;
    int i;

    branchNB = (BranchNB *)malloc(sizeof(BranchNB));
    if (!branchNB){
      ERROR_MESSAGE(); fprintf(stderr,"allocation failure in NewBranchNB()\n");
      exit(1);
    }
    branchNB->particles = particles;
    branchNB->nparticles = nparticles;

    for(i=0;i<treeNBdim;++i) branchNB->center[i]=center[i];
    branchNB->level=level;

    for(i=0;i<treeNBdim;++i){
      branchNB->boundary_p1[i]= boundary_p1[i];
      branchNB->boundary_p2[i]= boundary_p2[i];
    }

    branchNB->number=branchNBnumber;

    branchNB->child1 = NULL;
    branchNB->child2 = NULL;
    branchNB->prev = NULL;
    branchNB->brother = NULL;

    return(branchNB);
}

/************************************************************************
 * FreeBranchNB
 * Frees memory pointed to by branchNB.  Private.
 ************************************************************************/
void FreeBranchNB(BranchNB *branchNB){

    assert( branchNB != NULL);
    free(branchNB);

    return;
}

/************************************************************************
 * NewTreeNB
 * Returns pointer to new TreeNB struct.  Initializes top, last, and
 * current pointers to NULL.  Sets NbranchNBes field to 0.  Exported.
 ************************************************************************/
TreeNBHndl NewTreeNB(IndexType *particles,IndexType nparticles
		 ,PosType boundary_p1[treeNBdim],PosType boundary_p2[treeNBdim],
		     PosType center[treeNBdim],short Ndimensions){

    TreeNBHndl tree;
    
    tree = (TreeNBStruct *)malloc(sizeof(TreeNBStruct));
    if (!tree){
      ERROR_MESSAGE(); fprintf(stderr,"allocation failure in NewTreeNB()\n");
      exit(1);
    }

    tree->top= NewBranchNB(particles,nparticles,boundary_p1,boundary_p2,center,0,0);
    if (!(tree->top)){
      ERROR_MESSAGE(); fprintf(stderr,"allocation failure in NewTreeNB()\n");
      exit(1);
    }

    tree->Nbranches = 1;
    tree->current = tree->top;
    tree->Ndimensions=Ndimensions;

    return(tree);
}

void freeTreeNB(TreeNBHndl tree){
	/* free treeNB
	 *  does not free the particle positions, masses or rsph
	 */

	if(tree == NULL) return;

	emptyTreeNB(tree);
  	FreeBranchNB(tree->top);
	free(tree);

	return;
}

short emptyTreeNB(TreeNBHndl tree){

	moveTopNB(tree);
	_freeTreeNB(tree,0);

	assert(tree->Nbranches == 1);

	return 1;
}

void _freeTreeNB(TreeNBHndl tree,short child){
	BranchNB *branch;

	assert( tree );
	assert( tree->current);

	if(tree->current->child1 != NULL){
		moveToChildNB(tree,1);
		_freeTreeNB(tree,1);
	}

    if(tree->current->child2 != NULL){
      moveToChildNB(tree,2);
      _freeTreeNB(tree,2);
    }

    if( (tree->current->child1 == NULL)*(tree->current->child2 == NULL) ){

    	if(atTopNB(tree)) return;

    	branch = tree->current;
    	moveUpNB(tree);
       	FreeBranchNB(branch);

    	/*printf("*** removing branch %i number of branches %i\n",branch->number
			,tree->Nbranches-1);*/

       	if(child==1) tree->current->child1 = NULL;
    	if(child==2) tree->current->child2 = NULL;

    	--tree->Nbranches;

    	return;
    }

    return;
}
/***** Access functions *****/

/************************************************************************
 * isEmptyNB
 * Returns "True" if the TreeNB is empty and "False" otherwise.  Exported.
 ************************************************************************/
Boolean isEmptyNB(TreeNBHndl tree){

    assert(tree != NULL);
    return(tree->Nbranches == 0);
}

/************************************************************************
 * atTop
 * Returns "True" if current is the same as top and "False" otherwise.
 * Exported.
 * Pre: !isEmptyNB(tree)
 ************************************************************************/
Boolean atTopNB(TreeNBHndl tree){

    assert(tree != NULL);
    if( isEmptyNB(tree) ){
    	ERROR_MESSAGE();
    	fprintf(stderr, "TreeNB Error: calling atTop() on empty tree\n");
    	exit(1);
    }
    return(tree->current == tree->top);
}

/************************************************************************
 * noChild
 * Returns "True" if the child of the current branchNB does not exist and "False" otherwise.
 * Exported.
 * Pre: !isEmptyNB(tree)
 ************************************************************************/
Boolean noChildNB(TreeNBHndl tree){

    assert(tree != NULL);
    if( isEmptyNB(tree) ){
	
	ERROR_MESSAGE(); fprintf(stderr, "TreeNB Error: calling atTop() on empty tree\n");
	exit(1);
    }

    if( (tree->current->child1 == NULL) || (tree->current->child2 == NULL) ) return True;
    return False;
}

/************************************************************************
 * offEndNB
 * Returns "True" if current is off end and "False" otherwise.  Exported.
 ************************************************************************/
Boolean offEndNB(TreeNBHndl tree){

    assert(tree != NULL);
    return(tree->current == NULL);
}

/************************************************************************
 * getCurrentNB
 * Returns the particuls of current.  Exported.
 * Pre: !offEndNB(tree)
 ************************************************************************/
void getCurrentNB(TreeNBHndl tree,IndexType *particles,IndexType *nparticles){

    assert(tree != NULL);
    if( offEndNB(tree) ){
	
	ERROR_MESSAGE(); fprintf(stderr, "TreeNB Error: calling getCurrent() when current is off end\n");
	exit(1);
    }

    *nparticles=tree->current->nparticles;
    particles=tree->current->particles;

    return;
}

/************************************************************************
 * getNbranchesNB
 * Returns the NbranchNBes of tree.  Exported.
 ************************************************************************/
unsigned long getNbranchesNB(TreeNBHndl tree){

    assert(tree != NULL);
    return(tree->Nbranches);
}

/***** Manipulation procedures *****/

/************************************************************************
 * moveTopNB
 * Moves current to the front of tree.  Exported.
 * Pre: !isEmptyNB(tree)
 ************************************************************************/
void moveTopNB(TreeNBHndl tree){
    
    assert(tree != NULL);
    if( isEmptyNB(tree) ){
	
	ERROR_MESSAGE(); fprintf(stderr, "TreeNB Error: calling moveTopNB() on empty tree\n");
	exit(1);
    }

    tree->current = tree->top;
    return;
}

/************************************************************************
 * movePrev
 * Moves current to the branchNB before it in tree.  This can move current
 * off end.  Exported.
 * Pre: !offEndNB(tree)
 ************************************************************************/
void moveUpNB(TreeNBHndl tree){
    
    assert(tree != NULL);
    if( offEndNB(tree) ){
      ERROR_MESSAGE(); fprintf(stderr, "TreeNB Error: call to moveUpNB() when current is off end\n");
      exit(1);
    }
    if( tree->current == tree->top ){
      ERROR_MESSAGE(); fprintf(stderr, "TreeNB Error: call to moveUpNB() tried to move off the top\n");
      exit(1);
    }

    tree->current = tree->current->prev;  /* can move off end */
    return;
}

/************************************************************************
 * moveToChildNB
 * Moves current to child branchNB after it in tree.  This can move current off
 * end.  Exported.
 * Pre: !offEndNB(tree)
 ************************************************************************/
void moveToChildNB(TreeNBHndl tree,int child){
    
    assert(tree != NULL);
    if( offEndNB(tree) ){
	
	ERROR_MESSAGE(); fprintf(stderr, "TreeNB Error: calling moveChildren() when current is off end\n");
	exit(1);
    }
    if(child==1){
      if( tree->current->child1 == NULL ){
	ERROR_MESSAGE(); fprintf(stderr, "TreeNB Error: moveToChildNB() typing to move to child1 when it doesn't exist\n");
	exit(1);
      }
      tree->current = tree->current->child1;
    }
    if(child==2){
      if( tree->current->child2 == NULL ){
	ERROR_MESSAGE(); fprintf(stderr, "TreeNB Error: moveToChildNB() typing to move to child2 when it doesn't exist\n");
	exit(1);
      }
      tree->current = tree->current->child2;
    }
    return;
}

/************************************************************************
 * insertAfterCurrent
 * Inserts a new BranchNB after the current branchNB in the tree and sets the
 * data field of the new BranchNB to input.  Exported.
 * Pre: !offEndNB(tree)
 ************************************************************************/
void insertChildToCurrentNB(TreeNBHndl tree, IndexType *particles,IndexType nparticles
			  ,PosType boundary_p1[treeNBdim],PosType boundary_p2[treeNBdim]
			  ,PosType center[treeNBdim],int child){
    
    BranchNB *branchNB;

    /*printf("attaching child%i  current paricle number %i\n",child,tree->current->nparticles);*/

    branchNB = NewBranchNB(particles,nparticles,boundary_p1,boundary_p2,center
		       ,tree->current->level+1,tree->Nbranches);

    assert(tree != NULL);
    
    if( offEndNB(tree) ){
      
	ERROR_MESSAGE(); fprintf(stderr, "TreeNB Error: calling insertChildToCurrentNB() when current is off end\n");
	exit(1);
    }

    branchNB->prev = tree->current;

    if(child==1){
      if(tree->current->child1 != NULL){
	ERROR_MESSAGE(); fprintf(stderr, "TreeNB Error: calling insertChildToCurrentNB() when child1 alread exists\n");
	exit(1);
      }
      tree->current->child1 = branchNB;
    }
    if(child==2){
      if(tree->current->child2 != NULL){
    	  ERROR_MESSAGE();
    	  fprintf(stderr, "TreeNB Error: calling insertChildToCurrentNB() when child2 alread exists\n  current level=%i Nbranches=%li\n"
    			  ,tree->current->level,tree->Nbranches);
    	  exit(1);
      }
      tree->current->child2 = branchNB;      
    }

    tree->Nbranches++;

    return;
}

  /* same as above but takes a branchNB structure */

void attachChildToCurrentNB(TreeNBHndl tree,BranchNB data,int child){

  insertChildToCurrentNB(tree,data.particles,data.nparticles
		  ,data.boundary_p1,data.boundary_p2,data.center,child);
  return;
}

// step for walking tree by iteration instead of recursion
Boolean TreeNBWalkStep(TreeNBHndl tree,Boolean allowDescent){
	if(allowDescent && tree->current->child1 != NULL){
		moveToChildNB(tree,1);
		return True;
	}
	if(allowDescent && tree->current->child2 != NULL){
		moveToChildNB(tree,2);
		return True;
	}

	if(tree->current->brother != NULL){
		tree->current=tree->current->brother;
		return True;
	}
	return False;
}

/***** Other operations *****/

/************************************************************************
 * printTreeNB
 * Prints the contents of tree to stdout.  The current element, if there
 * is one, will be enclosed in []'s.  Currently, only to be used when the
 * TreeNBElements is are integers.  Exported.
 ************************************************************************/  
void printTreeNB(TreeNBHndl tree,PosType **xp){
    assert( tree !=NULL);

    printBranchNB(tree->current,xp,tree->Ndimensions);

    if(tree->current->child1 == NULL) return;

    moveToChildNB(tree,1);
    printTreeNB(tree,xp);

    moveUpNB(tree);

    if(tree->current->child2 == NULL) return;
    moveToChildNB(tree,2);
    printTreeNB(tree,xp);
    moveUpNB(tree);

    return;
}

void printBranchNB(BranchNB *data,PosType **xp,short Ndim){

  if(Ndim==3){
	  printf("******* branchNB *******\nlevel=%i\n",data->level);
	  printf("center = [%e,%e,%e]\n",data->center[0],data->center[1],data->center[2]);
	  printf("p1 = [%e,%e,%e] p2 = [%e,%e,%e]\n"
			  ,data->boundary_p1[0],data->boundary_p1[1],data->boundary_p1[2]
			  ,data->boundary_p2[0],data->boundary_p2[1],data->boundary_p2[2]);
  }
  if(Ndim==2){
	  printf("******* branchNB *******\nlevel=%i\n",data->level);
	  printf("center = [%e,%e]\n",data->center[0],data->center[1]);
	  printf("p1 = [%e,%e] p2 = [%e,%e]\n"
			  ,data->boundary_p1[0],data->boundary_p1[1]
			  ,data->boundary_p2[0],data->boundary_p2[1]);

  }
	  printf("number of particles = %li\n",data->nparticles);
  /*for(i=0;i<data->nparticles;++i) printf("%e %e %e\n",xp[data->particles[i]][0] ,xp[data->particles[i]][1],xp[data->particles[i]][2]); */
}

/**************************************************************************
   routines for saving and restoring tree structure
*************************************************************************

void saveSPHsmoothing(TreeNBHndl tree,IndexType *particles,float *rsph,char *filename){
	FILE *file;

  file=fopen(filename,"w");
  fwrite(&(tree->top->nparticles),sizeof(IndexType),1,file);
  fwrite(&(tree->Ndimensions),sizeof(short),1,file);
  fwrite(rsph,sizeof(float),tree->top->nparticles,file);
  fclose(file);

  return;
}

void saveTreeNB(TreeNBHndl tree,IndexType *particles,float *rsph,char *filename){
  RelativeBranchNB *tree_arr;
  void _saveTreeNB(TreeNBHndl tree,RelativeBranchNB *tree_arr,IndexType *particles);
  FILE *file;

  tree_arr=(RelativeBranchNB *)malloc(tree->NbranchNBes*sizeof(RelativeBranchNB));

  moveTopNB(tree);

  _saveTreeNB(tree,tree_arr,particles);
  
  file=fopen(filename,"w");
  fwrite(&(tree->top->nparticles),sizeof(IndexType),1,file);
  fwrite(&(tree->Ndimensions),sizeof(short),1,file);
  fwrite(&(tree->NbranchNBes),sizeof(IndexType),1,file);
  fwrite(tree_arr,sizeof(RelativeBranchNB),tree->NbranchNBes,file);
  fwrite(particles,sizeof(IndexType),tree->top->nparticles,file);
  fwrite(rsph,sizeof(float),tree->top->nparticles,file);
  fclose(file);

  free(tree_arr);

  return;
}

void _saveTreeNB(TreeNBHndl tree,RelativeBranchNB *tree_arr,IndexType *particles){
  int i;
  unsigned long current;

  current=tree->current->number;

  tree_arr[current].particles=tree->current->particles - tree->top->particles;
  tree_arr[current].nparticles=tree->current->nparticles;
  for(i=0;i<treeNBdim;++i){
    tree_arr[current].center[i]=tree->current->center[i];
    tree_arr[current].boundary_p1[i]=tree->current->boundary_p1[i];
    tree_arr[current].boundary_p2[i]=tree->current->boundary_p2[i];
  }

  if(current>0) tree_arr[current].prev=tree->current->prev->number;

  if(tree->current->child1 == NULL){
    tree_arr[current].child1=-1;
    return;
  }else{
    tree_arr[current].child1=tree->current->child1->number;//2*current;
    moveToChildNB(tree,1);
    _saveTreeNB(tree,tree_arr,particles);
    moveUpNB(tree);
  }

  if(tree->current->child2 == NULL){
    tree_arr[current].child2=-1;
    return;
  }else{
    tree_arr[current].child2=tree->current->child2->number;//=2*current+1;
    moveToChildNB(tree,2);
    _saveTreeNB(tree,tree_arr,particles);
    moveUpNB(tree);
  }

  return;
}

TreeNBHndl readTreeNB(IndexType *particles,float *rsph,IndexType Nparticles,char *filename){
  RelativeBranchNB *tree_arr;
  TreeNBHndl tree;
  FILE *file;
  IndexType nparticles;
  unsigned long NbranchNBes;
  short Ndimensions;
  void _readTreeNB(TreeNBHndl tree,RelativeBranchNB *tree_arr,IndexType *particles
		 ,unsigned long current);

  file=fopen(filename,"r");
  fread(&nparticles,sizeof(IndexType),1,file);
  fread(&Ndimensions,sizeof(short),1,file);
  //Ndimensions=3;
  fread(&NbranchNBes,sizeof(unsigned long),1,file);

  if(nparticles != Nparticles){
    ERROR_MESSAGE();
    printf("ERROR: input particle number does not match tree data file in readTreeNB\n");
    exit(1);
  }
  printf("nparticles %i   Nbranches %i  demensions=%i\n",nparticles,NbranchNBes,Ndimensions);

  tree_arr=(RelativeBranchNB *)malloc(NbranchNBes*sizeof(RelativeBranchNB));

  fread(tree_arr,sizeof(RelativeBranchNB),NbranchNBes,file);
  fread(particles,sizeof(IndexType),nparticles,file);
  fread(rsph,sizeof(float),nparticles,file);
  fclose(file);

  tree=NewTreeNB(particles,nparticles,tree_arr[0].boundary_p1,tree_arr[0].boundary_p2,tree_arr[0].center,Ndimensions);

  moveTopNB(tree);

  _readTreeNB(tree,tree_arr,particles,0);
  
  free(tree_arr);

  return tree;
}

void _readTreeNB(TreeNBHndl tree,RelativeBranchNB *tree_arr,IndexType *particles
	       ,unsigned long current){
  int i;

  tree->current->number=current;
  tree->current->particles=&particles[tree_arr[current].particles];
  tree->current->nparticles=tree_arr[current].nparticles;
  tree->current->level=tree_arr[current].level;

  for(i=0;i<treeNBdim;++i){
    tree->current->center[i]=tree_arr[current].center[i];
    tree->current->boundary_p1[i]=tree_arr[current].boundary_p1[i];
    tree->current->boundary_p2[i]=tree_arr[current].boundary_p2[i];
  }

  if(tree_arr[current].child1==-1){
    tree->current->child1=NULL;
    return;
  }else{
    insertChildToCurrentNB(tree,&current,0,&dummy,&dummy,&dummy,1);
    moveToChildNB(tree,1);
    _readTreeNB(tree,tree_arr,particles,tree_arr[current].child1);
    moveUpNB(tree);
  }

  if(tree_arr[current].child2==-1){
    tree->current->child2=NULL;
    return;
  }else{
    insertChildToCurrentNB(tree,&current,0,&dummy,&dummy,&dummy,2);
    moveToChildNB(tree,2);
    _readTreeNB(tree,tree_arr,particles,tree_arr[current].child2);
    moveUpNB(tree);
  }

  return;
}
*/
void readSmoothingNB(float *rsph,char *filename){
  /* just like readTreeNB only the full 3d tree is not returned  */
  /*   with the assumption that a 2d tree will be made later */

  //RelativeBranchNB *tree_arr;
  //TreeNBHndl tree;
  FILE *file;
  IndexType nparticles;
  short Ndimensions;


  file=fopen(filename,"r");
  fread(&nparticles,sizeof(IndexType),1,file);
  fread(&Ndimensions,sizeof(short),1,file);
  //fread(&NbranchNBes,sizeof(unsigned long),1,file);

  printf("nparticles %li  dimensions=%i\n",nparticles,Ndimensions);

  fread(rsph,sizeof(float),nparticles,file);
  fclose(file);
  

  //free(tree_arr);

  return;
}
