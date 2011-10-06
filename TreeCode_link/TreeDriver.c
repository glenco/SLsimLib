#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Tree.h"
#include <nrutil.h>
/*#include "double_sort.c"*/
#define Nbucket 1   // must be =1 if each leaf is to coincide with each cell

/* median_cut determines how the cells are subdivided */
/*    if ==0  equal volume cuts, Warning this option causes an error*/
/*    if ==1  median point cuts */
static int incell,median_cut=1;
static Point **temp_points,tmp;
static double realray[2];
Point *point_global;

TreeHndl BuildTree(Point *xp,unsigned long Npoints){
  TreeHndl tree;
  unsigned long i;
  double p1[2],p2[2],center[2];
  void _BuildTree(TreeHndl tree);

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
  tree=NewTree(xp,Npoints,p1,p2,center);

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
  if(cbranch->npoints <= Nbucket){
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

Point *NearestNeighbor(TreeHndl tree,double *ray,int Nneighbors,ListHndl neighborlist,short direction){
  /* nearest neighbor points in a given direction, direction != 0 should not */
  /*    be used except on a grid */
  /* direction = 0 distance */
  /*           = 1 to right */
  /*           = 2 to left */
  /*           = 3 to up */
  /*           = 4 to down */
  unsigned long i;
  void _NearestNeighbor(TreeHndl tree,double *ray,int Nneighbors
			,Point **neighborpoints,double *rneighbor,short *direction);
  static int count=0,oldNneighbors=-1;
  static double *rneighbors;
  static Point **neighborpoints;

  /*printf("entering NN\n");*/

  if(tree->top->npoints <= Nneighbors){
	  ERROR_MESSAGE();
    printf("ERROR: in NearestNeighbor, number of neighbors > total number of points\n");
    exit(1);
  }

  if(count==0){
    /*printf("allocating memory\n");*/
    rneighbors=(double *)malloc((Nneighbors+Nbucket)*sizeof(double));
    assert(rneighbors);
    neighborpoints=(Point **)malloc((Nneighbors+Nbucket)*sizeof(Point *));
    assert(neighborpoints);
    temp_points=(Point **)malloc((Nneighbors+Nbucket)*sizeof(Point *));
    assert(temp_points);
    ++count;
    oldNneighbors=Nneighbors;

  }else if(oldNneighbors < Nneighbors){ /* if the number of nearest neighbors goes up get more mem */
    /*printf("re-allocating memory\n");*/
    rneighbors=(double *)realloc(rneighbors,(Nneighbors+Nbucket)*sizeof(double));
    neighborpoints=(Point **)realloc(neighborpoints,(Nneighbors+Nbucket)*sizeof(Point *));
    temp_points=(Point **)realloc(temp_points,(Nneighbors+Nbucket)*sizeof(Point *));
    oldNneighbors=Nneighbors;
  }

  /*   printf("Nneighbors=%i\n",Nneighbors); */
  /*   printf("array sizes=%i\n",Nneighbors+Nbucket); */

  /* initalize distance to neighbors to a large number */
  for(i=0;i<Nbucket+Nneighbors;++i){
    rneighbors[i]=10*(tree->top->boundery_p2[0]-tree->top->boundery_p1[0]);
  }

  moveTop(tree);
//   printf("p1= [%f,%f]\n", tree->current->boundery_p1[0],tree->current->boundery_p1[1]);
//   printf("p2= [%f,%f]\n", tree->current->boundery_p2[0],tree->current->boundery_p2[1]);

  realray[0]=ray[0];
  realray[1]=ray[1];

  if( inbox(ray,tree->current->boundery_p1,tree->current->boundery_p2) == 0 ){
    printf("Warning: in NearestNeighbor, ray is not inside the simulation box\n    should work in any case\n      ray= %e %e\n",ray[0],ray[1]);

    ray[0]=DMAX(ray[0],tree->current->boundery_p1[0]);
    ray[0]=DMIN(ray[0],tree->current->boundery_p2[0]);

    ray[1]=DMAX(ray[1],tree->current->boundery_p1[1]);
    ray[1]=DMIN(ray[1],tree->current->boundery_p2[1]);
  }
  incell=1;

  if(direction==0) EmptyList(neighborlist);
  if(direction==-1) direction=0;

  _NearestNeighbor(tree,ray,Nneighbors,neighborpoints,rneighbors,&direction);

  /* convert from point array to exported point list */

  for(i=0;i<Nneighbors;++i){
/*     printf("hello i=%i Nneighbors=%i x= %e %e %i\n",i,Nneighbors */
/* 	   ,neighborpoints[i]->x[0],neighborpoints[i]->x[1],neighborpoints[i]->image->id); */
    InsertAfterCurrent(neighborlist,neighborpoints[i]->x,neighborpoints[i]->id,neighborpoints[i]->image);
/*     printf("did insert\n"); */
    MoveDownList(neighborlist);
    /*neighborlist->current->invmag=neighborpoints[i]->invmag;*/
    PointCopyData(neighborlist->current,neighborpoints[i]);/**/
  }

  /*  printf("returning from NN\n");*/

  return neighborpoints[0];
}

void _NearestNeighbor(TreeHndl tree,double *ray,int Nneighbors,Point **neighborpoints,double *rneighbors,short *direction){

  int i,incell2=1;
  unsigned long index[Nneighbors+Nbucket];
  double dx,dy;

  //printf("**************************************\nlevel %i\n",tree->current->level);
  //for(i=0;i<tree->current->npoints;++i) printf("   %i\n",tree->current->points[i]);

  if(incell){  /* not found cell yet */

    if( inbox(ray,tree->current->boundery_p1,tree->current->boundery_p2) ){

      /* found the box small enough */
      if( tree->current->npoints <= Nneighbors+Nbucket ){
	incell=0;
	/*printf("found box with %i points\n",tree->current->npoints);*/

	/* this sets ray back to real value once closest leaf bax is found */
	if( (ray[0]!=realray[0])*(ray[1]!=realray[1]) ){ printf("ray != realray _NearestNeighbor\n"); exit(0);}

	ray[0]=realray[0];
	ray[1]=realray[1];

	/* calculate the distance to all the points in cell */
	tree->pointlist->current=tree->current->points;
	for(i=0;i<tree->current->npoints;++i){

	  dx=tree->pointlist->current->x[0]-ray[0];
	  dy=tree->pointlist->current->x[1]-ray[1];

	  switch(*direction){
	  case 0: /* distance */
	    rneighbors[i] = sqrt(dx*dx+dy*dy);
	    break;
	  case 1: /* right */
	    if(dx>0 && fabs(dy/dx)<0.999) rneighbors[i]=sqrt(dx*dx+dy*dy);
	    else rneighbors[i]=1.0e99;
	    break;
	  case 2: /* left */
	    if(dx<0 && fabs(dy/dx)<0.999) rneighbors[i]=sqrt(dx*dx+dy*dy);
	    else rneighbors[i]=1.0e99;
	    break;
	  case 3: /* up */
	    if(dy>0 && fabs(dx/dy)<0.999) rneighbors[i]=sqrt(dx*dx+dy*dy);
	    else rneighbors[i]=1.0e99;
	    break;
	  case 4: /* down */
	    if(dy<0 && fabs(dx/dy)<0.999) rneighbors[i]=sqrt(dx*dx+dy*dy);
	    else rneighbors[i]=1.0e99;
	    break;
	  }

	  /*neighbors[i]=tree->pointlist->current->id;*/
	  index[i]=i;
	  temp_points[i]=tree->pointlist->current;
	  MoveDownList(tree->pointlist);
	}

	/*printf("first sort at level =%i\n",tree->current->level);*/
	/*printf("N=%i\n",tree->current->npoints);*/
	/*for(i=0;i<tree->current->npoints;++i) printf("  %i  %e\n",neighbors[i],rneighbors[i]);*/
	double_sort(tree->current->npoints,rneighbors-1,index-1);
	for(i=0;i<tree->current->npoints;++i) neighborpoints[i]=temp_points[index[i]];

      }else{ /* keep going down the tree */

	/*printf("moving to child1 from level %i\n",tree->current->level);*/
	  if(tree->current->child1 !=NULL){
	      moveToChild(tree,1);
	      _NearestNeighbor(tree,ray,Nneighbors,neighborpoints,rneighbors,direction);
	      /*printf("moving up from level %i\n",tree->current->level);*/
	      moveUp(tree);

	      incell2=incell;
	  }

	  if(tree->current->child2 !=NULL){
	      /*printf("moving to child2 from level %i\n",tree->current->level);*/
	      moveToChild(tree,2);
	      _NearestNeighbor(tree,ray,Nneighbors,neighborpoints,rneighbors,direction);
	      /*printf("moving up from level %i\n",tree->current->level);*/
	      moveUp(tree);
	  }

	/** if ray found in second child go back to first to search for neighbors **/
	if( (incell2==1) && (incell==0) ){
	  if(tree->current->child1 !=NULL){
	      /*printf("moving to child1 again from level %i\n",tree->current->level);*/
	      moveToChild(tree,1);
	      _NearestNeighbor(tree,ray,Nneighbors,neighborpoints,rneighbors,direction);
	      /*printf("moving up from level %i\n",tree->current->level);*/
	      moveUp(tree);
	  }
	}
      }
    } /* not in the box */

  }else{  /* already found cell */
    /*printf("finding neighboring boxes at level = %i\n",tree->current->level);*/

    /* does radius cut into the box */
    if( cutbox(ray,tree->current->boundery_p1,tree->current->boundery_p2,rneighbors[Nneighbors-1]) ){

      if( (tree->current->child1 == NULL)*(tree->current->child2 == NULL)){  /* leaf case */

	/*printf("entering leaf box level %i number of points %i Nbucket=%i\n"
	  ,tree->current->level,tree->current->npoints,Nbucket);*/

	/* combine found neighbors with points in box and resort */
	tree->pointlist->current=tree->current->points;
	for(i=0;i<Nneighbors;++i){
	  index[i]=i;
	  temp_points[i]=neighborpoints[i];
	}

	for(i=Nneighbors;i<(tree->current->npoints+Nneighbors);++i){
	  /*printf("i=%i\n",i);*/
/* 	  for(j=0,rneighbors[i]=0.0;j<2;++j){ */
/* 	    rneighbors[i]+=pow(tree->pointlist->current->x[j]-ray[j],2); */
/* 	  } */
/* 	  rneighbors[i]=sqrt( rneighbors[i] ); */

	  dx=tree->pointlist->current->x[0]-ray[0];
	  dy=tree->pointlist->current->x[1]-ray[1];

	  switch(*direction){
	  case 0: /* distance */
	    rneighbors[i] = sqrt(dx*dx+dy*dy);
	    break;
	  case 1: /* right */
	    if(dx>0 && fabs(dy/dx)<0.999) rneighbors[i]=sqrt(dx*dx+dy*dy);
	    else rneighbors[i]=1.0e99;
	    break;
	  case 2: /* left */
	    if(dx<0 && fabs(dy/dx)<0.999) rneighbors[i]=sqrt(dx*dx+dy*dy);
	    else rneighbors[i]=1.0e99;
	    break;
	  case 3: /* up */
	    if(dy>0 && fabs(dx/dy)<0.999) rneighbors[i]=sqrt(dx*dx+dy*dy);
	    else rneighbors[i]=1.0e99;
	    break;
	  case 4: /* down */
	    if(dy<0 && fabs(dx/dy)<0.999) rneighbors[i]=sqrt(dx*dx+dy*dy);
	    else rneighbors[i]=1.0e99;
	    break;
	  }

	  index[i]=i;
	  temp_points[i]=tree->pointlist->current;
	  MoveDownList(tree->pointlist);
	}
	/*for(i=(tree->current->npoints+Nneighbors);i<Nneighbors+Nbucket;++i) index[i]=i;*/

	/*for(i=0;i<(tree->current->npoints+Nneighbors);++i) 
	  printf("   i=%i rneighbor=%e neighbor=%i N=%i\n",i,rneighbors[i],neighbors[i],Nneighbors+Nbucket);*/

	double_sort(tree->current->npoints+Nneighbors,rneighbors-1,index-1);
	/*for(i=Nneighbors;i<(tree->current->npoints+Nneighbors);++i) neighborpoints[i]=temp_points[index[i]];*/
	/*printf("hello %i %i\n",Nneighbors+Nbucket,tree->current->npoints+Nneighbors);*/
	/*for(i=0;i<(tree->current->npoints+Nneighbors);++i) printf("index=%i\n",index[i]);*/
	for(i=0;i<(tree->current->npoints+Nneighbors);++i) neighborpoints[i]=temp_points[index[i]];
      }else{

	/*printf("moving to child1 from level %i\n",tree->current->level);*/
	    if(tree->current->child1 !=NULL){
		moveToChild(tree,1);
		_NearestNeighbor(tree,ray,Nneighbors,neighborpoints,rneighbors,direction);
		/*printf("moving up from level %i\n",tree->current->level);*/
		moveUp(tree);
	    }

	    if(tree->current->child2 !=NULL){
		/*printf("moving to child2 from level %i\n",tree->current->level);*/
		moveToChild(tree,2);
		_NearestNeighbor(tree,ray,Nneighbors,neighborpoints,rneighbors,direction);
		/*printf("moving up from level %i\n",tree->current->level);*/
		moveUp(tree);
	    }
      }

    }/*else{printf("box too distant at level %i\n",tree->current->level);}*/

  }

  /*printf("end of _NearestNeighbor incell=%i level=%i p1= %e %e %e\n",incell,tree->current->level
    ,tree->current->boundery_p1[0],tree->current->boundery_p1[1],tree->current->boundery_p1[2]);*/
  return;
}

/* return 1 (0) if ray is (not) in the cube */
inline int inbox(double *ray,double *p1,double *p2){

  return (ray[0]>=p1[0])*(ray[0]<=p2[0])*(ray[1]>=p1[1])*(ray[1]<=p2[1]);
}

// returns true if branch1 is fully inside barnch2
Boolean boxinbox(Branch *branch1,Branch *branch2){

	if(inbox(branch1->boundery_p1,branch2->boundery_p1,branch2->boundery_p2) == 0) return False;
	if(inbox(branch1->boundery_p2,branch2->boundery_p1,branch2->boundery_p2) == 0) return False;

	return True;
}
double BoxIntersection(Branch *branch1,Branch *branch2){
	// returns area of intersection between two branches
	double area=0;

	area = MIN(branch1->boundery_p2[0],branch2->boundery_p2[0])
	     - MAX(branch1->boundery_p1[0],branch2->boundery_p1[0]);
	if(area < 0) return 0.0;

	area *= MIN(branch1->boundery_p2[1],branch2->boundery_p2[1])
	      - MAX(branch1->boundery_p1[1],branch2->boundery_p1[1]);
	if(area < 0) return 0.0;

	return area;
}

// return 1 (0) if box is (not) within rmax of ray
int cutbox(double *ray,double *p1,double *p2,double rmax){
	/*  returns:  0 if whole box is outside rmax from ray[]
	 *            1 if whole box is inside circle but ray is not in the box
	 *            2 if ray[] is inside box
	 *            3 if box intersects circle but ray[] is not inside box
	 */
  short i,tick=0;
  double close[2],rtmp;
  
  // find closest point on box borders to ray[]
  for(i=0;i<2;++i){
    if( ray[i] < p1[i] ){
      close[i]=p1[i];
    }else if(ray[i] > p2[i]){
      close[i]=p2[i];
    }else{
      close[i]=ray[i];
      ++tick;
    }
  }
  
  if(tick==2) return 2;  // ray is inside box

  for(i=0,rtmp=0;i<2;++i) rtmp += pow(ray[i] - close[i],2);
 
  if(rtmp>rmax*rmax) return 0;  // box is all outside circle

  // find farthest point on box border from ray[]
  for(i=0,rtmp=0;i<2;++i) rtmp+=DMAX(pow(ray[i]-p1[i],2),pow(ray[i]-p2[i],2));

  if(rtmp<rmax*rmax) return 1;  // box is all inside circle

  return 3;  // box intersects circle
}

/****************************************************************
 *
 ****************************************************************/

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
    tree->current->npoints -= Nadd;
	x=(double *) malloc(2*Nbucket*sizeof(double));
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

    		if( tree->current->npoints > Nbucket ){ // create new leaves

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

void _FindLeaf(TreeHndl tree,double *ray,unsigned long Nadd){
	Boolean contin;
/*   printf("***********************\n"); */
/*   printf("level = %i  npoints = %i incell=%i\n", tree->current->level,tree->current->npoints,incell); */
/*   printf("p1= [%f,%f]\n", tree->current->boundery_p1[0],tree->current->boundery_p1[1]); */
/*   printf("p2= [%f,%f]\n", tree->current->boundery_p2[0],tree->current->boundery_p2[1]); */
/*   printf("first point = %i x= %f %f\n",tree->current->points->id,tree->current->points->x[0],tree->current->points->x[1]); */
/*   printf("ray = %f %f\n",ray[0],ray[1]); */

	assert(inbox(ray,tree->current->boundery_p1,tree->current->boundery_p2) );
	do{
		tree->current->npoints += Nadd;

		/* change center of mass */
		tree->current->center[0]=(tree->current->center[0]*(tree->current->npoints-1) + ray[0])/tree->current->npoints;
		tree->current->center[1]=(tree->current->center[1]*(tree->current->npoints-1) + ray[1])/tree->current->npoints;

		contin=False;
		if(tree->current->child1 != NULL &&
				inbox(ray,tree->current->child1->boundery_p1,tree->current->child1->boundery_p2) ){
			moveToChild(tree,1);
			contin=True;
		}else if(tree->current->child2 != NULL &&
				inbox(ray,tree->current->child2->boundery_p1,tree->current->child2->boundery_p2) ){
			moveToChild(tree,2);
			contin=True;
		}
	}while(contin);
/*
  if(tree->current->child1 !=NULL){
    if( inbox(ray,tree->current->child1->boundery_p1,tree->current->child1->boundery_p2) ){
      moveToChild(tree,1);
      _FindLeaf(tree,ray,Nadd);
    }
  }
  if(tree->current->child2 !=NULL){
    if( inbox(ray,tree->current->child2->boundery_p1,tree->current->child2->boundery_p2) ){
      moveToChild(tree,2);
      _FindLeaf(tree,ray,Nadd);
    }
  }
*/
  return;
}

Boolean AreBoxNeighbors(Point *point1,Point *point2){

	if(    point1->leaf->boundery_p1[0] <= point2->leaf->boundery_p2[0]
	    && point1->leaf->boundery_p2[0] >= point2->leaf->boundery_p1[0]
	    && point1->leaf->boundery_p1[1] <= point2->leaf->boundery_p2[1]
	    && point1->leaf->boundery_p2[1] >= point2->leaf->boundery_p1[1] ) return True;
	return False;
}

void FindBoxPoint(TreeHndl tree,double *ray,Point *point){
	// return the point that is in the same box as ray[2]
	// if Nbuck > 1 the head of the point array is returned
	void _FindBox(TreeHndl tree,double *ray);
	Branch *branch;

	//point=(Point *)malloc(sizeof(Point));
	//moveTop(tree);
	//printTree(tree);

	branch=tree->current;
	moveTop(tree);
	   // check if ray is outside initial box
	if( inbox(ray,tree->current->boundery_p1,tree->current->boundery_p2) == 0 ){
		printf("FindBox: ray outside of grided range\n");
		return;
	}

	_FindBox(tree,ray);
	PointCopyData(point,tree->current->points);

	//if(foundpoint == 0){ printf("FindBoxPoint failed to find point\n"); exit(1);}
	//PrintPoint(point);

	// error check
	if(fabs(ray[0]-point->x[0]) > point->gridsize/2
			|| fabs(ray[1]-point->x[1]) > point->gridsize/2){
		ERROR_MESSAGE();
		printf("ERROR: FindBox did not find box\n  ray = %e %e\n  Delta/gridsize = %e %e\n"
				,ray[0],ray[1]
		        ,2*(ray[0]-point->x[0])/point->gridsize
				,2*(ray[1]-point->x[1])/point->gridsize);
		exit(0);
	}

	tree->current=branch;
	//return point;
}

void _FindBox(TreeHndl tree,double *ray){
	Boolean contin;

	// replaced recursion with iteration
	do{
		contin=False;
		if(tree->current->child1 !=NULL &&
				inbox(ray,tree->current->child1->boundery_p1,tree->current->child1->boundery_p2) ){
			moveToChild(tree,1);
			contin=True;
		}
		if(tree->current->child2 !=NULL &&
				inbox(ray,tree->current->child2->boundery_p1,tree->current->child2->boundery_p2) ){
			moveToChild(tree,2);
			contin=True;
		}
	}while(contin);

	return;
}


/** simple sort for points in linked list **/
/** slow for > 20 points **/
Point *sortList(long n, double *arr,ListHndl list,Point *firstpointin){
  long i,j;
  double a;
  Point *point,*firstpoint;

  firstpoint=firstpointin;

  for(j=1;j<n;j++){
    a=arr[j];

    list->current=firstpoint;
    JumpDownList(list,j);

    point=TakeOutCurrent(list);

    i=j-1;
    while(i>-1 && arr[i] > a){
      arr[i+1]=arr[i];
      i--;
      MoveUpList(list);
      /*printf("      current= %i %f %f\n",list->current->id,list->current->x[0],list->current->x[1]);*/
    }
    arr[i+1]=a;

    if( list->current==list->top && i==-1) InsertPointBeforeCurrent(list,point);
    else InsertPointAfterCurrent(list,point);

    if(i == -1) firstpoint=point;
  }
  return firstpoint;
}


void PointsWithin(TreeHndl tree,double *ray,float rmax,ListHndl neighborlist,short markpoints){
	void _PointsWithin(TreeHndl tree,double *ray,float *rmax,ListHndl neighborlist,short markpoints);
/*
 *  finds all points in tree that lie within rmax of the point ray[]
 *   markpoints = 0  does not change in_image variable in any point, gives a list of neighbors
 *              = 1  makes in_image=True for all points in image, gives no list of neighbors
 *              = -1 makes in_image=False for all points in image to reset, gives no list of neighbors
 */

  if(markpoints==0) EmptyList(neighborlist);
  
  realray[0]=ray[0];
  realray[1]=ray[1];

  moveTop(tree);
  if( inbox(ray,tree->current->boundery_p1,tree->current->boundery_p2) == 0 ){
    printf("Warning: in PointsWithin, ray is not inside the simulation box\n    should work in any case\n      ray= %e %e\n     boundery p1 = %e %e p2 = %e %e\n",ray[0],ray[1]
	   ,tree->current->boundery_p1[0],tree->current->boundery_p1[1]
	   ,tree->current->boundery_p2[0],tree->current->boundery_p2[1]);

    ray[0]=DMAX(ray[0],tree->current->boundery_p1[0]);
    ray[0]=DMIN(ray[0],tree->current->boundery_p2[0]);

    ray[1]=DMAX(ray[1],tree->current->boundery_p1[1]);
    ray[1]=DMIN(ray[1],tree->current->boundery_p2[1]);
  }
  incell=1;

  _PointsWithin(tree,ray,&rmax,neighborlist,markpoints);
}

void PointsWithin_iter(TreeHndl tree,double *ray,float rmax,ListHndl neighborlist,short markpoints){
	Boolean descend;
	short desition;
	unsigned long i,j;
	double radius,length;

	assert(neighborlist);
	assert(tree);

	EmptyList(neighborlist);

	moveTop(tree);

//	if( inbox(ray,tree->current->boundery_p1,tree->current->boundery_p2) ){
//	_FindBox(tree,ray);
//	while(ClosestBorder(ray,tree->current->boundery_p1,tree->current->boundery_p2)
//			< rmax) moveUp(tree);

	do{
		descend = True;
		desition = cutbox(ray,tree->current->boundery_p1,tree->current->boundery_p2,rmax);
		length = FurthestBorder(ray,tree->current->boundery_p1,tree->current->boundery_p2);
		if( FurthestBorder(ray,tree->current->boundery_p1,tree->current->boundery_p2) < rmax ){

			//ClosestBorder(ray,tree->current->boundery_p1,tree->current->boundery_p2) < rmax) ){
			// whole box is outside circle
			descend = False;

			// put all the points into neighborlist
		   	  tree->pointlist->current=tree->current->points;
    		  for(i=0;i<tree->current->npoints;++i){
     			  if(markpoints == 1){
     				  tree->pointlist->current->in_image=True;
     				  tree->pointlist->current->image->in_image=True;
     			  }else if(markpoints == -1){
     				  tree->pointlist->current->in_image=False;
     				  tree->pointlist->current->image->in_image=False;
					  tree->pointlist->current->surface_brightness = tree->pointlist->current->image->surface_brightness = 0.0;
    			  }else if(markpoints == 0){
     				  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
							  ,tree->pointlist->current->id,tree->pointlist->current->image);

					  MoveDownList(neighborlist);
					  PointCopyData(neighborlist->current,tree->pointlist->current);
					  MoveUpList(neighborlist);
				  }

     			  MoveDownList(tree->pointlist);
     		  }

		}else if(cutbox(ray,tree->current->boundery_p1,tree->current->boundery_p2,rmax) == 0 ){  // whole box is outside circle
			descend = False;
		}else if(tree->current->child1 == NULL){  // leaf case

		   	  tree->pointlist->current=tree->current->points;
	   		  for(i=0;i<tree->current->npoints;++i){
	    			  for(j=0,radius=0.0;j<2;++j) radius+=pow(tree->pointlist->current->x[j]-ray[j],2);
	    			  if( radius < rmax*rmax ){
	       				  if(markpoints == 1){
	       					  tree->pointlist->current->in_image=True;
	      					  tree->pointlist->current->image->in_image=True;
	      				  }else if(markpoints == -1){
	      					  tree->pointlist->current->in_image=False;
	     					  tree->pointlist->current->image->in_image=False;
	     					  tree->pointlist->current->surface_brightness = tree->pointlist->current->image->surface_brightness = 0.0;
	       				  }else if(markpoints == 0){
	         				  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
	    						  ,tree->pointlist->current->id,tree->pointlist->current->image);

	         				  MoveDownList(neighborlist);
	         				  PointCopyData(neighborlist->current,tree->pointlist->current);
	         				  MoveUpList(neighborlist);
	      				  }
	    			  }
	    			  MoveDownList(tree->pointlist);
	    		  }
		}
	}while(TreeWalkStep(tree,descend));
}

double ClosestBorder(double *ray,double *p1,double *p2){
	/*  returns the distance from ray[] to the closest
	 *    border of the box,
	 *    < 0 if ray is outside of the box
	 *    > 0 if ray is inside the box
	 */
	double length;

	length = MIN(ray[0]-p1[0],p2[0]-ray[0]);
	length = MIN(ray[1]-p1[1],length);

	return MIN(p2[1]-ray[1],length);
}
inline double FurthestBorder(double *ray,double *p1,double *p2){
	/*  returns the distance from ray[] to the furthest point on the
	 *    border of the box,
	 */

	return sqrt( pow(MAX(ray[0]-p1[0],p2[0]-ray[0]),2) + pow(MAX(ray[1]-p1[1],p2[1]-ray[1]),2) );
}


void _PointsWithin(TreeHndl tree,double *ray,float *rmax,ListHndl neighborlist,short markpoints){

  int i,j,incell2=1;
  double radius;
  short pass;


  //printf("**************************************\nlevel %i\n",tree->current->level);
  //   printf("   %i incell=%i\n",tree->current->points->id,incell);

  if(incell){  // not found cell yet

    if( inbox(ray,tree->current->boundery_p1,tree->current->boundery_p2) ){

      // found the box small enough
    	if( cutbox(ray,tree->current->boundery_p1,tree->current->boundery_p2,*rmax)==1
    			|| (tree->current->child1 == NULL)*(tree->current->child2 == NULL) ){
    		// whole box in circle or a leaf with ray in it

    	  incell=0;
    	  //printf("found box with %i points\n",tree->current->npoints);

    	  // this sets ray back to real value once closest leaf bax is found
    	  if( (ray[0]!=realray[0])*(ray[1]!=realray[1]) ){ printf("ray != realray _PointsWithin\n"); exit(0);}

    	  ray[0]=realray[0];
    	  ray[1]=realray[1];

    	  tree->pointlist->current=tree->current->points;

    	  if((tree->current->child1 == NULL)*(tree->current->child2 == NULL)){
    	   	  // if leaf calculate the distance to all the points in cell
    		  for(i=0;i<tree->current->npoints;++i){
    			  for(j=0,radius=0.0;j<2;++j) radius+=pow(tree->pointlist->current->x[j]-ray[j],2);
    			  if( radius < *rmax**rmax ){
       				  if(markpoints == 1){
       					  tree->pointlist->current->in_image=True;
      					  tree->pointlist->current->image->in_image=True;
      				  }else if(markpoints == -1){
      					  tree->pointlist->current->in_image=False;
     					  tree->pointlist->current->image->in_image=False;
     					  tree->pointlist->current->surface_brightness = tree->pointlist->current->image->surface_brightness = 0.0;
       				  }else if(markpoints == 0){
         				  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
    						  ,tree->pointlist->current->id,tree->pointlist->current->image);

         				  MoveDownList(neighborlist);
         				  PointCopyData(neighborlist->current,tree->pointlist->current);
         				  MoveUpList(neighborlist);
      				  }
    			  }
    			  MoveDownList(tree->pointlist);
    		  }
    	  }else{ // put all of points in box into neighborlist
       		  for(i=0;i<tree->current->npoints;++i){
       			  if(markpoints == 1){
       				  tree->pointlist->current->in_image=True;
       				  tree->pointlist->current->image->in_image=True;
       			  }else if(markpoints == -1){
       				  tree->pointlist->current->in_image=False;
       				  tree->pointlist->current->image->in_image=False;
 					  tree->pointlist->current->surface_brightness = tree->pointlist->current->image->surface_brightness = 0.0;
      			  }else if(markpoints == 0){
       				  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
  							  ,tree->pointlist->current->id,tree->pointlist->current->image);

  					  MoveDownList(neighborlist);
  					  PointCopyData(neighborlist->current,tree->pointlist->current);
  					  MoveUpList(neighborlist);
  				  }

       			  MoveDownList(tree->pointlist);
       		  }
    	  }

    	}else{ // keep going down the tree

    	  //printf("moving to child1 from level %i\n",tree->current->level);
    	  if(tree->current->child1 !=NULL){
    		  moveToChild(tree,1);
    		  _PointsWithin(tree,ray,rmax,neighborlist,markpoints);
    		  //printf("moving up from level %i\n",tree->current->level);
    		  moveUp(tree);

    		  incell2=incell;
    	  }

    	  if(tree->current->child2 !=NULL){
    		  //printf("moving to child2 from level %i\n",tree->current->level);
    		  moveToChild(tree,2);
    		  _PointsWithin(tree,ray,rmax,neighborlist,markpoints);
    		  //printf("moving up from level %i\n",tree->current->level);
    		  moveUp(tree);
    	  }

    	  // if ray found in second child go back to first to search for neighbors
    	  if( (incell2==1) && (incell==0) ){
    		  if(tree->current->child1 !=NULL){
    			  //printf("moving to child1 again from level %i\n",tree->current->level);
    			  moveToChild(tree,1);
    			  _PointsWithin(tree,ray,rmax,neighborlist,markpoints);
    			  //printf("moving up from level %i\n",tree->current->level);
    			  moveUp(tree);
    		  }
    	  }
      }
    }  // not in the box

  }else{    // found cell

	  //printf("finding neighboring boxes at level = %i\n",tree->current->level);

	  pass=cutbox(ray,tree->current->boundery_p1,tree->current->boundery_p2,*rmax);
	  // does radius cut into the box
	  if( pass ){

		  if( (tree->current->child1 == NULL)*(tree->current->child2 == NULL)  ){  /* leaf case */

			  tree->pointlist->current=tree->current->points;
			  for(i=0;i<tree->current->npoints;++i){

				  for(j=0,radius=0.0;j<2;++j) radius+=pow(tree->pointlist->current->x[j]-ray[j],2);
				  if( radius < *rmax**rmax ){
					  if(markpoints==1){
						  tree->pointlist->current->in_image=True;
						  tree->pointlist->current->image->in_image=True;
					  }else if(markpoints==-1){
						  tree->pointlist->current->in_image=False;
						  tree->pointlist->current->image->in_image=False;
     					  tree->pointlist->current->surface_brightness = tree->pointlist->current->image->surface_brightness = 0.0;
					  }else if(markpoints==0){
						  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
								  ,tree->pointlist->current->id,tree->pointlist->current->image);

           				  MoveDownList(neighborlist);
           				  PointCopyData(neighborlist->current,tree->pointlist->current);
           				  MoveUpList(neighborlist);
      				  }

				  }
				  MoveDownList(tree->pointlist);
			  }
		  }else if(pass==1){ // whole box is inside radius
			  tree->pointlist->current=tree->current->points;
			  for(i=0;i<tree->current->npoints;++i){
  				  if(markpoints==1){
   					  tree->pointlist->current->in_image=True;
  					  tree->pointlist->current->image->in_image=True;
  				  }else if(markpoints==-1){
  					  tree->pointlist->current->in_image=False;
 					  tree->pointlist->current->image->in_image=False;
 					  tree->pointlist->current->surface_brightness = tree->pointlist->current->image->surface_brightness = 0.0;
   				  }else if(markpoints==0){
   					  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
   							  ,tree->pointlist->current->id,tree->pointlist->current->image);

   					  MoveDownList(neighborlist);
   					  PointCopyData(neighborlist->current,tree->pointlist->current);
   					  MoveUpList(neighborlist);
   				  }

				  MoveDownList(tree->pointlist);
			  }
		  }else{
			  //printf("moving to child1 from level %i\n",tree->current->level);
			  if(tree->current->child1 !=NULL){
				  moveToChild(tree,1);
				  _PointsWithin(tree,ray,rmax,neighborlist,markpoints);
				  //printf("moving up from level %i\n",tree->current->level);
				  moveUp(tree);
			  }

			  if(tree->current->child2 !=NULL){
				  //printf("moving to child2 from level %i\n",tree->current->level);
				  moveToChild(tree,2);
				  _PointsWithin(tree,ray,rmax,neighborlist,markpoints);
				  //printf("moving up from level %i\n",tree->current->level);
				  moveUp(tree);
			  }
		  }

	  }
  }

	  //  printf("end of _PointsWithin incell=%i level=%i p1= %e %e %e\n",incell,tree->current->level
	//	,tree->current->boundery_p1[0],tree->current->boundery_p1[1],tree->current->boundery_p1[2]);/**/
  return;
}

void NeighborsOfNeighbors(ListHndl neighbors,ListHndl wholelist){
	/* finds all the neighbors of neighbors of the
	 * point neighbors->current that are part in the list wholelist
	 * the neighbors are added to neighbors and removed from wholelist
	 */
	short check=0;
	Point *point;

	if(neighbors->Npoints < 1 || wholelist->Npoints < 1) return;

	do{
		MoveToTopList(wholelist);
		do{
			if(check==1){
				MoveUpList(wholelist);
				check=0;
			}
			if(AreBoxNeighbors(neighbors->current,wholelist->current) ){
				if(AtTopList(wholelist)) check=1;
				point=TakeOutCurrent(wholelist);
				InsertPointAfterCurrent(neighbors,point);
			}
		}while(MoveDownList(wholelist));
	}while(MoveDownList(neighbors) && wholelist->Npoints > 0);

	return ;
}

void FriendsOfFriends(TreeHndl tree,double *start_point,float linkinglength,ListHndl neighborlist
		      ,Point *filter,unsigned long Nfilter,unsigned long *filter_place){
  unsigned long Nlocal_filter,i;
  Point *placemark,*local_filter;
  short malloced=0;

/*   printf("entering FOF\n"); */

  if(filter == NULL){ printf("FriendsOfFriends cannot handle no unfiltered points\n"); return;}

  EmptyList(neighborlist);

/*   if(filter != NULL){   */
/*     printf("filter first id %i\n   x start = %e %e\n",filter[*filter_place].id,start_point[0],start_point[1]); */
/*     MoveToTopList(tree->pointlist); */
/*     for(i=0;i<tree->pointlist->Npoints;++i){ */
/*       if(filter[*filter_place].id == tree->pointlist->current->id) printf("x in tree %e %e\n" */
/* 					 ,tree->pointlist->current->x[0],tree->pointlist->current->x[1]); */
/*       MoveDownList(tree->pointlist); */
/*     } */
/*   } */

  incell=1;
  realray[0]=start_point[0];
  realray[1]=start_point[1];

  moveTop(tree);
  _PointsWithin2(tree,start_point,&linkinglength,neighborlist,filter,Nfilter,filter_place,1);

  if(filter == NULL &&  neighborlist->Npoints > 0){
    Nlocal_filter=neighborlist->Npoints;
    //local_filter=(Point *)malloc(Nlocal_filter*sizeof(Point));
	local_filter=NewPointArray(Nlocal_filter,False);

    MoveToTopList(neighborlist);
    for(i=0;i<Nlocal_filter;++i){
      PointCopyData(&local_filter[i],neighborlist->current);
      MoveDownList(neighborlist);
    }
    malloced=1;
  }

  MoveToTopList(neighborlist);
  if( neighborlist->Npoints > 0 ){
    for(;;){

      if(filter==NULL){
	/* find points that are not in filter */
	realray[0]=neighborlist->current->x[0];
	realray[1]=neighborlist->current->x[1];

	_PointsWithin2(tree,neighborlist->current->x,&linkinglength,neighborlist
		      ,local_filter,Nlocal_filter,filter_place,0);

	/* add found points to filter */
	if(Nlocal_filter < neighborlist->Npoints){
	  local_filter=(Point *) realloc(local_filter,neighborlist->Npoints*sizeof(Point));
	  local_filter->head = neighborlist->Npoints;
	  //local_filter = AddPointToArray(local_filter,neighborlist->Npoints,unsigned long Nold);
	  placemark=neighborlist->current;
	  MoveDownList(neighborlist);
	  for(i=Nlocal_filter;i<=neighborlist->Npoints;++i){
	    PointCopyData(&local_filter[i],neighborlist->current);
	    MoveDownList(neighborlist);
	  }
	  neighborlist->current=placemark;
	  Nlocal_filter=neighborlist->Npoints;
	}
      }else{
	realray[0]=neighborlist->current->x[0];
	realray[1]=neighborlist->current->x[1];

	_PointsWithin2(tree,neighborlist->current->x,&linkinglength,neighborlist,filter,Nfilter,filter_place,1);
      }

      if(neighborlist->current == neighborlist->bottom) break;
      else MoveDownList(neighborlist);
    }
  }else{
	ERROR_MESSAGE();
    printf("ERROR: no neighbors in Friends of Friends\n");
    exit(0);
  }

  if(malloced) free(local_filter);

/*   printf("returning from FOF\n"); */
  return;
}

void _PointsWithin2(TreeHndl tree,double *ray,float *rmax,ListHndl neighborlist
		   ,Point *filter,unsigned long Nfilter,unsigned long *filter_place,short compliment){

/*   id numbers in filter[*filterplace ... Nfilter-1] */
/*   reorders filter[] and moves filterplace as points are found */
/*   if filter=NULL no filter is used */
/*   compliment = 0 search only points not in the filter */
/*              = 1 search only points in filter after filter_place */

  int i,j,incell2=1;
  double radius;
  unsigned long k;

  //if(filter != NULL) printf("filter_place=%i Nneighborlist=%i\n",*filter_place,neighborlist->Npoints);

  if((filter != NULL) && (*filter_place >= Nfilter) ) return; // finished filter

  //printf("**************************************\nlevel %i\n",tree->current->level);
  //   printf("   %i incell=%i\n",tree->current->points->id,incell);

  if(incell){  // not found cell yet

    if( inbox(ray,tree->current->boundery_p1,tree->current->boundery_p2) ){

      // found the box small enough
      if( (tree->current->child1 == NULL)*(tree->current->child2 == NULL)){  // leaf case
    	  incell=0;
    	  //printf("found box with %i points\n",tree->current->npoints);

    	  // this sets ray back to real value once closest leaf bax is found
    	  if( (ray[0]!=realray[0])*(ray[1]!=realray[1]) ){ printf("ray != realray _PointsWithin\n"); exit(0);}

    	  ray[0]=realray[0];
    	  ray[1]=realray[1];

    	  // calculate the distance to all the points in cell

    	  tree->pointlist->current=tree->current->points;
    	  for(i=0;i<tree->current->npoints;++i){

    		  for(j=0,radius=0.0;j<2;++j) radius+=pow(tree->pointlist->current->x[j]-ray[j],2);
    		  if( radius < *rmax**rmax ){
    			  if(filter == NULL){

    				  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
    						  ,tree->pointlist->current->id,tree->pointlist->current->image);

    				  MoveDownList(neighborlist);
    				  PointCopyData(neighborlist->current,tree->pointlist->current);
    				  MoveUpList(neighborlist);

    			  }else if(compliment){

    				  for(k=*filter_place;k<Nfilter;++k){
    					  if(tree->pointlist->current->id == filter[k].id){
    						  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
    								  ,tree->pointlist->current->id,tree->pointlist->current->image);

    						  MoveDownList(neighborlist);
    						  PointCopyData(neighborlist->current,tree->pointlist->current);
    						  MoveUpList(neighborlist);

    						  // move point in filter to filter_place and increment filter_place
    						  PointCopyData(&tmp,&filter[k]);
    						  for(;k>*filter_place;k--) PointCopyData(&filter[k],&filter[k-1]);
    						  PointCopyData(&filter[*filter_place],&tmp);
    						  ++*filter_place;
    						  break;
    					  }
    				  }
    			  }else{ // do compliment of filter

    				  for(k=0;k<Nfilter;++k) if(tree->pointlist->current->id == filter[k].id) break;
    				  if(k==Nfilter){
    					  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
    							  ,tree->pointlist->current->id,tree->pointlist->current->image);

    					  MoveDownList(neighborlist);
    					  PointCopyData(neighborlist->current,tree->pointlist->current);
    					  MoveUpList(neighborlist);
    				  }
    			  }

    		  }
    		  MoveDownList(tree->pointlist);
    	  }
      }else{ // keep going down the tree

    	  //printf("moving to child1 from level %i\n",tree->current->level);
    	  if(tree->current->child1 !=NULL){
    		  moveToChild(tree,1);
    		  _PointsWithin2(tree,ray,rmax,neighborlist,filter,Nfilter,filter_place,1);
    		  //printf("moving up from level %i\n",tree->current->level);
    		  moveUp(tree);

    		  incell2=incell;
    	  }

    	  if(tree->current->child2 !=NULL){
    		  //printf("moving to child2 from level %i\n",tree->current->level);
    		  moveToChild(tree,2);
    		  _PointsWithin2(tree,ray,rmax,neighborlist,filter,Nfilter,filter_place,1);
    		  //printf("moving up from level %i\n",tree->current->level);
    		  moveUp(tree);
    	  }

    	  // if ray found in second child go back to first to search for neighbors
    	  if( (incell2==1) && (incell==0) ){
    		  if(tree->current->child1 !=NULL){
    			  //printf("moving to child1 again from level %i\n",tree->current->level);
    			  moveToChild(tree,1);
    			  _PointsWithin2(tree,ray,rmax,neighborlist,filter,Nfilter,filter_place,1);
    			  //printf("moving up from level %i\n",tree->current->level);
    			  moveUp(tree);
    		  }
    	  }
      }
    }  // not in the box
    //printf("not in box \n");}
  }else{    // found cell

	  //printf("finding neighboring boxes at level = %i\n",tree->current->level);

	  // does radius cut into the box
	  if( cutbox(ray,tree->current->boundery_p1,tree->current->boundery_p2,*rmax) ){

		  if( (tree->current->child1 == NULL)*(tree->current->child2 == NULL)){  /* leaf case */

			  tree->pointlist->current=tree->current->points;
			  for(i=0;i<tree->current->npoints;++i){

				  for(j=0,radius=0.0;j<2;++j) radius+=pow(tree->pointlist->current->x[j]-ray[j],2);
				  if( radius < *rmax**rmax ){
					  if(filter == NULL){

						  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
								  ,tree->pointlist->current->id,tree->pointlist->current->image);

						  MoveDownList(neighborlist);
						  PointCopyData(neighborlist->current,tree->pointlist->current);
						  MoveUpList(neighborlist);

					  }else if(compliment){

						  for(k=*filter_place;k<Nfilter;++k){
							  if(tree->pointlist->current->id == filter[k].id){
								  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
										  ,tree->pointlist->current->id,tree->pointlist->current->image);

								  MoveDownList(neighborlist);
								  PointCopyData(neighborlist->current,tree->pointlist->current);
								  MoveUpList(neighborlist);

								  PointCopyData(&tmp,&filter[k]);
								  for(;k>*filter_place;k--) PointCopyData(&filter[k],&filter[k-1]);
								  PointCopyData(&filter[*filter_place],&tmp);
								  ++*filter_place;
								  break;
							  }
						  }

					  }else{ // do compliment of filter

						  for(k=0;k<Nfilter;++k) if(tree->pointlist->current->id == filter[k].id) break;
						  if(k==Nfilter){
							  InsertAfterCurrent(neighborlist,tree->pointlist->current->x
									  ,tree->pointlist->current->id,tree->pointlist->current->image);

							  MoveDownList(neighborlist);
							  PointCopyData(neighborlist->current,tree->pointlist->current);
							  MoveUpList(neighborlist);

						  }
					  }
				  }
				  MoveDownList(tree->pointlist);
			  }

		  }else{
			  //printf("moving to child1 from level %i\n",tree->current->level);
			  if(tree->current->child1 !=NULL){
				  moveToChild(tree,1);
				  _PointsWithin2(tree,ray,rmax,neighborlist,filter,Nfilter,filter_place,1);
				  //printf("moving up from level %i\n",tree->current->level);
				  moveUp(tree);
			  }

			  if(tree->current->child2 !=NULL){
				  //printf("moving to child2 from level %i\n",tree->current->level);
				  moveToChild(tree,2);
				  _PointsWithin2(tree,ray,rmax,neighborlist,filter,Nfilter,filter_place,1);
				  //printf("moving up from level %i\n",tree->current->level);
				  moveUp(tree);
			  }
		  }

	  }//else{printf("box too distant at level %i\n",tree->current->level);}
  }

	  //  printf("end of _PointsWithin incell=%i level=%i p1= %e %e %e\n",incell,tree->current->level
	//	,tree->current->boundery_p1[0],tree->current->boundery_p1[1],tree->current->boundery_p1[2]);/**/
  return;
}

void TransformPoints(ListHndl listout,ListHndl listin){
  /* tranforms between image and source planes */
  /* the image pointers in listin must be set properly */
  /* the image pointers of listout will be set to listin members - it would be better to fix this */

  unsigned long i;
  Point *placemarker;

  placemarker=listin->current;
  EmptyList(listout);
  MoveToTopList(listin);
  for(i=0;i<listin->Npoints;++i){
    InsertAfterCurrent(listout,listin->current->image->x,listin->current->image->id,listin->current);
    MoveDownList(listin);
  }
  listin->current=placemarker;
}

Boolean ArePointsUniqueList(ListHndl list){
	long i,j;
	Point *point,*init;

	init=list->current;
	MoveToTopList(list);
	for(i=0;i<list->Npoints-1;++i){
		point=list->current->next;
		for(j=i+1;j<list->Npoints;++j){
			if(list->current->id==point->id){
				list->current=init;
				return False;
			}
			point=point->next;
		}
		MoveDownList(list);
	}

	list->current=init;
	return True;
}

Boolean IntersectionList(ListHndl list1,ListHndl list2,ListHndl intersection){
	long i,j;
	Point *init1,*init2;

	EmptyList(intersection);
	init1=list1->current;
	init2=list2->current;

	MoveToTopList(list1);
	for(i=0;i<list1->Npoints;++i){
		MoveToTopList(list2);
		for(j=0;j<list2->Npoints;++j){
			if(list1->current->id==list2->current->id){
				InsertAfterCurrent(intersection,list1->current->x,list1->current->id
						,list1->current->image);
				MoveDownList(intersection);
				PointCopyData(intersection->current,list1->current);
			}
			MoveDownList(list2);
		}
		MoveDownList(list1);
	}

	list1->current=init1;
	list2->current=init2;

	if(intersection->Npoints < 1) return False;
	return True;
}

void UnionList(ListHndl list1,ListHndl list2){
	// on exit list1 will be the union of the two lists
	// and list 2 will still be the subset
	list1->bottom->next=list2->top;
	list2->top->prev=list1->bottom;
	list1->bottom=list2->bottom;
	list1->Npoints=list1->Npoints + list2->Npoints;

	return ;
}

void ClearAllMarks(TreeHndl tree){
	unsigned long i;

	MoveToTopList(tree->pointlist);
	for(i=0;i<tree->pointlist->Npoints;++i){
		tree->pointlist->current->in_image=False;
		tree->pointlist->current->image->in_image=False;
		MoveDownList(tree->pointlist);
	}
}

void FindAllBoxNeighbors(TreeHndl tree,Point *point,ListHndl neighbors){
	// finds all the leaves that are neighboring branch
	// points outside of grid have no box neighbors
	void _FindAllBoxNeighbors(TreeHndl tree,Branch *leaf,ListHndl neighbors);
	static int count=0;

	++count;
	EmptyList(neighbors);

	// point is outside of initial region
	if(!inbox(point->x,tree->top->boundery_p1,tree->top->boundery_p2)) return;

	tree->current=point->leaf;

	// find smallest box that surrounds box and its neighbors
	moveUp(tree);
	while( (tree->current->boundery_p1[0]==point->leaf->boundery_p1[0] && point->leaf->boundery_p1[0] != tree->top->boundery_p1[0] )
			|| (tree->current->boundery_p1[1]==point->leaf->boundery_p1[1] && point->leaf->boundery_p1[1] != tree->top->boundery_p1[1])
			|| (tree->current->boundery_p2[0]==point->leaf->boundery_p2[0] && point->leaf->boundery_p2[0] != tree->top->boundery_p2[0])
			|| (tree->current->boundery_p2[1]==point->leaf->boundery_p2[1] && point->leaf->boundery_p2[1] != tree->top->boundery_p2[1]) ){
		moveUp(tree);
	}

	_FindAllBoxNeighbors(tree,point->leaf,neighbors);

	return;
}

// this should be made into a loop instead of a recursion
void _FindAllBoxNeighbors(TreeHndl tree,Branch *leaf,ListHndl neighbors){

	if(  leaf->boundery_p1[0]<=tree->current->boundery_p2[0]
	  && leaf->boundery_p2[0]>=tree->current->boundery_p1[0]
	  && leaf->boundery_p1[1]<=tree->current->boundery_p2[1]
	  && leaf->boundery_p2[1]>=tree->current->boundery_p1[1]){

		if(tree->current->npoints == Nbucket){
			if(tree->current->number != leaf->number){

				//InsertAfterCurrentKist(neighbors,tree->current->points);

				InsertAfterCurrent(neighbors,tree->current->points->x
					,tree->current->points->id
					,tree->current->points->image);
				MoveDownList(neighbors);
				PointCopyData(neighbors->current,tree->current->points);

			}
			return;
		}

		if(tree->current->child1 !=NULL){
			moveToChild(tree,1);
			_FindAllBoxNeighbors(tree,leaf,neighbors);
			moveUp(tree);
		}

		if(tree->current->child2 !=NULL){
			moveToChild(tree,2);
			_FindAllBoxNeighbors(tree,leaf,neighbors);
			moveUp(tree);
		}
	}

	return;
}

void PrintImages(ImageInfo *images,long Nimages){
	long i,j;

	printf("%li",Nimages);
	for(i=0;i<Nimages;++i){
		printf("%li\n",images[i].Npoints);
		for(j=0;j<images[i].Npoints;++j)
			printf("%e %e  %e\n",images[i].points[j].x[0],images[i].points[j].x[1]
			                         ,images[i].points[j].gridsize);
	}
}

void PrintImageInfo(ImageInfo *image){

	printf(" PrintImageInfo\n");
	printf("  Npoints = %li  area = %e +/- %e\n",image->Npoints,image->area,image->area_error);
	printf("  gridrange = %e %e %e\n",image->gridrange[0],image->gridrange[1],image->gridrange[2]);
	printf("  borders inner N = %li  outer N = %li\n",image->innerborder->Nunits,image->outerborder->Nunits);
}
