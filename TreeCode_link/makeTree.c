#include <math.h>
#include <stdio.h>
#include <time.h>
#include "../../Library/Recipes/nr.h"
#include "../../Library/Recipes/nrutil.h"
#include "../../Library/Recipes/nrutil.c"
#include "../../Library/RecipesD/nrD.h"
#include "../../Library/RecipesD/locateD.c"
#include "../../Library/RecipesD/odeintD.c"
#include "../../Library/RecipesD/bsstepD.c"
#include "../../Library/RecipesD/mmidD.c"
#include "../../Library/RecipesD/pzextrD.c"
#include "../../Library/RecipesD/polintD.c"
#include "../../Library/cosmo.h"
#include "../../Library/cosmo.c"

#include "Tree.h"
#include "Tree.c"
#include "TreeNBDriver.c"
#include "TreeNBForce.c"
#include "readfiles.c"

struct cosmology cosmo;
SimLens lens;
int kmax,kount; 
double *xp,**yp,dxsav;

int main(int arg,char **argv){
 
  int i;
  time_t to,t1;
  clock_t tcpu;
  float tmp;

  /* read parameter file and simulation */

  printf("READING PARAMETERS FROM %s\n\n",argv[1]);
  readparams(argv[1],&lens,&cosmo);
  printf("\nreading particle positions from %s\n\n",lens.simfilename);
  readpositions(&lens);
  printf("number of particles=%i\n",lens.Nparticles);

  /*for(i=0;i<Nparticles;++i) printf("%e %e %e\n",xp[i][0],xp[i][1],xp[i][2]);*/

  /*********************************/

  printf("\n***********************************************\n");
  printf("********* building tree ***********************\n\n");

  to=clock();
  lens.tree=BuildTree(lens.xp,lens.Nparticles,lens.particles);
  t1=clock();
  printf("%f sec to build tree of %i particles\n",(float)(clock()-to)/CLOCKS_PER_SEC,lens.Nparticles);


  printf("\n***********************************************\n");
  printf("********* calculate SPH smoothing *************\n\n");

  to=clock();
  lens.rsph=FindRSPH(lens.tree,lens.xp,lens.Nsph);
  time(&t1);
  printf("%f sec to find %i particle SPH smoothing for %i particles\n",(float)(clock()-to)/CLOCKS_PER_SEC,lens.Nsph,lens.Nparticles);


  printf("\n***********************************************\n");
  printf("********* saving tree *************************\n\n");

  saveTree(lens.tree,lens.particles,lens.rsph,lens.treefilenames);

  printf("tree saved to %s\n\n",lens.treefilenames);

  return 0;
}

