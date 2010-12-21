/*#include "/Applications/itt/idl/external/include/idl_export.h"*/
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "../../Library/Recipes/nrutil.h"
#include "../../Library/Recipes/nrutil.c"
#include "../../Library/RecipesD/locateD.c"
#include "../../Library/Recipes/ran2.c"

#include "../../Library/cosmo.h"

#include "../TreeCode_link/utilities.c"
#include "../TreeCode_link/Tree.h"
#include "../TreeCode_link/double_sort.c"
#include "../TreeCode_link/Tree.c"
#include "../TreeCode_link/TreeDriver.c"

#include "rayshooter.c"


int main(){
  double *ray,*alpha,*gamma,*kappa,*invmag;
  long Nrays,seed=-1;
  char paramfile[30];
  time_t to,t1,t2,t3;
  int i;

  Nrays=1;
  ray=dvector(0,2*Nrays-1);
  ray[0]=ray[1]=1.0e-3;

  alpha=dvector(0,2*Nrays-1);
  gamma=dvector(0,2*Nrays-1);
  kappa=dvector(0,Nrays-1);
  invmag=dvector(0,Nrays-1);

  sprintf(paramfile,"paramfileSIE1M");
  printf("hi\n");

  for(i=0;i<25;++i){
    if(i==1) to=clock();
    ray[0]=1.0e-2*ran2D(&seed);
    ray[1]=1.0e-2*ran2D(&seed);

    t1=clock();
    printf("before shoot\n");
    rayshooter(ray,Nrays,alpha,gamma,kappa,invmag,paramfile);
    printf("     %f sec\n",(float)(clock()-t1)/CLOCKS_PER_SEC);
  }
  printf("total = %f sec\n",(float)(clock()-to)/CLOCKS_PER_SEC/24);
}
