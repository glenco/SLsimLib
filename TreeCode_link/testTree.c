#include <math.h>
#include <stdio.h>
#include <time.h>
#include "../../Library/Recipes/nrutil.h"
#include "../../Library/Recipes/nrutil.c"
/*#include "../Include/nr.h"*/
#include "../../Library/Recipes/ran2.c"
#include "../../Library/RecipesD/locateD.c"
#include "../../Library/cosmo.h"

#include "Tree.h"
#include "Tree.c"
#include "TreeDriver.c"
#include "image_finder.c"
#include "rayshooter_stars.c"

int main(int arg,char **argv){
 
  TreeHndl tree;
  ListHndl points,neighborlist;
  int i;
  unsigned long Npoints,Nfilter,filter_place,*filter;
  double length,*ray;
  long seed=0;
  time_t to,t1;
  clock_t tcpu;
  Point *i_points,*s_points;

  int Nsph;
  Point **neighbors;
  float tmp;
  FILE *treefile;

  Npoints=600;

  /*********************************/

  if(atoi(argv[1])){
  /*** make fake data ***/

    i_points=NewPointArray(Npoints,True);
    s_points=NewPointArray(Npoints,True);

    /* make some random points */
    length=1.0;
    for(i=0;i<Npoints;++i){
      i_points[i].id=i;

      i_points[i].x[0]=length*(ran2(&seed)-0.5);
      i_points[i].x[1]=length*(ran2(&seed)-0.5);

      /* link images and source points */
      i_points[i].image=&s_points[i];
      s_points[i].image=&i_points[i];

      /*xp[i][0]=length*(ran2(&seed)-0.5);
	xp[i][1]=length*(ran2(&seed)-0.5);*/

      /*printf("   %e  %e\n",i_points[i].x[0],i_points[i].x[1]);/**/
    }

    printf("\n***********************************************\n");
    printf("********* building tree ***********************\n\n");

    to=clock();
    tree=BuildTree(i_points,Npoints);
    t1=clock();
    printf("%f sec to build tree of %i points\n"
	   ,(float)(clock()-to)/CLOCKS_PER_SEC,Npoints);

  }else{

    if(arg == 3){
      printf("\n***********************************************\n");
      printf("********* reading tree *************************\n\n");

      tree=readTree(argv[2]);

      PrintList(tree->pointlist);

      seed=18298;
      length=1.0;
/*       for(i=0;i<Npoints;++i){ */
/* 	xp[i][0]=length*(ran2(&seed)-0.5); */
/* 	xp[i][1]=length*(ran2(&seed)-0.5); */
/* 	printf("   %e  %e\n",xp[i][0],xp[i][1]);/\**\/ */
/*       } */
    }else{
      printf("need input file name\n");
      exit(0);
    }

  }


  printf("\n***********************************************\n");
  printf("********* insert points *******************\n\n");

  /* make some more points to add */

  i_points=AddPointToArray(i_points,Npoints+8,Npoints);
  s_points=AddPointToArray(s_points,Npoints+8,Npoints);

  for(i=Npoints;i<(Npoints+8);++i){
    i_points[i].id=i;

    i_points[i].x[0]=length*(ran2(&seed)-0.5);
    i_points[i].x[1]=length*(ran2(&seed)-0.5);

      /* link images and source points */
    i_points[i].image=&s_points[i];
    s_points[i].image=&i_points[i];
   
    AddPointsToTree(tree,&i_points[i]);/**/
  }
  PrintList(tree->pointlist);/**/

  neighborlist=NewList();
  printf("\n***********************************************\n");
  printf("********* finding neighbors *******************\n\n");

  printf("neighbors to %f %f\n",i_points[10].x[0],i_points[10].x[1]);
  NearestNeighbor(tree,i_points[10].x,7,neighborlist);
  PrintList(neighborlist);/**/

/*   for(i=0;i<7;++i) printf("%i    %f %f  %f\n" */
/* 			  ,neighbors[i]->id,neighbors[i]->x[0] */
/* 			  ,neighbors[i]->x[1] */
/* 			  ,sqrt( pow(neighbors[i]->x[0]-i_points[10].x[0],2) */
/* 				 + pow(neighbors[i]->x[1]-i_points[10].x[1],2) ) ); */

  printf("\n******************************************************\n");
  printf("********* find points within radius of point ***********\n\n");

  PointsWithin(tree,i_points[10].x,0.2,neighborlist);
  PrintList(neighborlist);/**/

  printf("\n******************************************************\n");
  printf("********* find freinds-of-freinds ***********\n\n");

  /* make filter */

  filter=(unsigned long *) malloc(neighborlist->Npoints*sizeof(unsigned long));
  MoveToTopList(neighborlist);
  for(i=0;i<neighborlist->Npoints;++i){
    filter[i]=neighborlist->current->id;
    MoveDownList(neighborlist);
  }
  Nfilter=neighborlist->Npoints;
  filter_place=0;

  EmptyList(neighborlist);
  FriendsOfFriends(tree,i_points[10].x,0.05,neighborlist,filter,Nfilter,&filter_place);

  PrintList(neighborlist);/**/

  exit(0);

  if(arg == 3 ){
    printf("\n***********************************************\n");
    printf("********* saving tree *************************\n\n");

    saveTree(tree,argv[2]);

    printf("tree saved to %s\n",argv[2]);
  }

  return 0;
}
