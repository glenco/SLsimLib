#include <math.h>
#include <stdio.h>
#include <time.h>
#include "../../Library/Recipes/nrutil.h"
#include "../../Library/Recipes/nrutil.c"
/*#include "../Include/nr.h"*/
#include "../../Library/Recipes/ran2.c"
#include "../../Library/RecipesD/locateD.c"

#include "../../Library/cosmo.h"

#include "TreeNB.h"
#include "TreeNB.c"
#include "TreeNBDriver.c"
#include "TreeNBForce.c"
#include "minimization.c"

int main(int arg,char **argv){
 
  TreeNBHndl tree;
  int i;
  unsigned long Nparticles,*particles;
  double **xp;
  double length,*ray;
  long seed=0;
  time_t to,t1;
  clock_t tcpu;

  int Nsph;
  unsigned long *neighbors;
  float *rsph,tmp;
  double **coord,theta,alpha[2],kappa,gamma[2];
  FILE *treefile;

  Nparticles=1000000;

  /*printf("%i %i %i %i %i %i\n",sizeof(short),sizeof(int),sizeof(long),sizeof(unsigned long),sizeof(float),sizeof(double));*/

  /*for(i=0;i<=18;++i) printf("%i %i \n",i,i % 3);*/


  /*** make fake data ***/
  xp=PosTypeMatrix(0,Nparticles-1,0,2);

  length=1.0;
  for(i=0;i<Nparticles;++i){
    xp[i][0]=length*(ran2(&seed)-0.5);
    xp[i][1]=length*(ran2(&seed)-0.5);
    xp[i][2]=length*(ran2(&seed)-0.5);

    /*printf("   %f  %f  %f\n",xp[i][0],xp[i][1],xp[i][2]);*/
  }

  /*********************************/

  particles=(unsigned long *)malloc(Nparticles*sizeof(unsigned long));

  if(atoi(argv[1])){
    printf("\n***********************************************\n");
    printf("********* building tree ***********************\n\n");

    to=clock();
    tree=BuildTreeNB(xp,Nparticles,particles);
    t1=clock();
    printf("%f sec to build tree of %i particles\n",(float)(clock()-to)/CLOCKS_PER_SEC,Nparticles);


    printf("\n***********************************************\n");
    printf("********* calculate SPH smoothing *************\n\n");

    Nsph=64;

    to=clock();
    rsph=FindRSPH(tree,xp,Nsph);
    time(&t1);
    printf("%f sec to find %i particle SPH smoothing for %i particles\n",(float)(clock()-to)/CLOCKS_PER_SEC,Nsph,Nparticles);
    /*printTreeNB(tree,xp);*/


    printf("\n***********************************************\n");
    printf("********* saving tree *************************\n\n");

    saveTreeNB(tree,particles,rsph,argv[2]);

    printf("tree saved to %s\n",argv[2]);
  }else{

    printf("\n***********************************************\n");
    printf("********* reading tree *************************\n\n");

    rsph=(float *)malloc(Nparticles*sizeof(float));
    tree=readTreeNB(particles,rsph,Nparticles,argv[2]);
  }

  /* coordinates defining the line of sight projection */
  coord=dmatrix(0,1,0,2);
  coord[0][0]=1; coord[0][1]=0; coord[0][2]=0;
  coord[1][0]=0; coord[1][1]=1; coord[1][2]=0;
  theta=0.1;

  printf("\n***********************************************\n");
  printf("********* calculate lensing quantities ************\n\n");


  ray=(double *)malloc(3*sizeof(double));

  time(&to);
  for(i=0;i<=5;++i){
    ray[0]=length*(ran2(&seed)-0.5);
    ray[1]=length*(ran2(&seed)-0.5);
    ray[2]=length*(ran2(&seed)-0.5);

    tcpu=clock();
    TreeNBForce2D(tree,xp,rsph,ray,coord,theta,alpha,&kappa,gamma,False);
    /*printf("%f sec per ray\n",(float)(clock()-tcpu)/CLOCKS_PER_SEC);*/

    printf("alpha= %e %e kappa=%e gammma= %e %e\n",alpha[0],alpha[1],kappa,gamma[0],gamma[1]);

  }
  time(&t1);
  printf("%f sec of wall time find force at 100 rays \n",difftime(t1,to),Nparticles);

  printf("alpha= %e %e kappa=%e gammma= %e %e\n",alpha[0],alpha[1],kappa,gamma[0],gamma[1]);
  printf("at ray = %e %e %e particle smoothing\n",ray[0],ray[1],ray[2]);


/*   printf("\n***********************************************\n"); */
/*   printf("********* find nearest neighbors **************\n\n"); */


/*   Nsph=5; */
/*   ray[0]=0; ray[1]=0.25; ray[2]=0; */
/*   neighbors=NearestNeighborNB(tree,xp,ray,Nsph,&tmp); */
/*   printf("closest %i particles to ray = %e %e %e\n",Nsph,ray[0],ray[1],ray[2]); */
/*   for(i=0;i<Nsph;++i) printf("   %i   %e %e %e\n",neighbors[i] */
/* 				   ,xp[neighbors[i]][0],xp[neighbors[i]][1],xp[neighbors[i]][2]); */


/*   printf("Out of NearestNeighborNB\n"); */
  /*  for(i=0;i<Nsph;++i) printf("%e %e %e\n",part[neighbors[i]][0],part[neighbors[i]][1],part[neighbors[i]][2]);*/

  return 0;
}
