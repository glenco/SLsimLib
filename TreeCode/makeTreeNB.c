#include <math.h>
#include <stdio.h>
#include <time.h>
#include "../../Library/Recipes/nr.h"
#include "../../Library/Recipes/nrutil.h"
#include "../../Library/Recipes/nrutil.c"
#include "../../Library/Recipes/ran2.c"

#include "../../Library/RecipesD/nrD.h"
#include "../../Library/RecipesD/locateD.c"
#include "../../Library/RecipesD/odeintD.c"
#include "../../Library/RecipesD/bsstepD.c"
#include "../../Library/RecipesD/mmidD.c"
#include "../../Library/RecipesD/pzextrD.c"
#include "../../Library/RecipesD/polintD.c"
#include "../../Library/RecipesD/ran2D.c"
#include "../../Library/RecipesD/powellD.c"
#include "../../Library/RecipesD/dfridrD.c"
#include "../../Library/RecipesD/integratorD.c"
#include "../../Library/RecipesD/rootfinder2D.c"

#include "../../Library/Recipes/gasdev.c"
#include "../../Library/Recipes/poidev.c"
#include "../../Library/Recipes/gammln.c"

#include "../../Library/cosmo.h"
#include "../../Library/cosmo.c"

#include "../../Library/powerCDM.c"
#include "TreeNB.h"

struct cosmology cosmo;
SimLens *lens;
int kmax,kount; 
double *xp,**yp,dxsav;

char *paramfile,*outputfile;

int main(int arg,char **argv){
  int dimension=3;
  time_t to,t1;

  /* read parameter file and simulation */
  printf("READING PARAMETERS FROM %s\n\n",argv[1]);
  lens=readparams(argv[1],&cosmo);
  if(arg>2) dimension=atof(argv[2]);
  printf("\nreading particle positions from %s\n\n",lens->simfilename);
  readpositions(lens);
  printf("number of particles=%li\n",lens->Nparticles);
  printf("critical surface density = %e M/Mpc^2\n",lens->Sigma_crit);

  /*for(i=0;i<Nparticles;++i) printf("%e %e %e\n",xp[i][0],xp[i][1],xp[i][2]);*/

  /*********************************/

  printf("\n***********************************************\n");
  printf("********* building tree ***********************\n\n");

  to=clock();
  lens->tree=BuildTreeNB(lens->xp,lens->rsph,lens->Nparticles,lens->particles,dimension);

  t1=clock();
  printf("%f sec to build tree of %i particles\n     dimension: %li",
		  (float)(clock()-to)/CLOCKS_PER_SEC,lens->Nparticles,lens->tree->Ndimensions);


  printf("\n***********************************************\n");
  printf("********* calculate SPH smoothing *************\n\n");

  to=clock();
  lens->rsph=FindRSPH(lens->tree,lens->xp,lens->Nsph);
  time(&t1);
  printf("%f sec to find %i particle SPH smoothing for %i particles\n",(float)(clock()-to)/CLOCKS_PER_SEC,lens->Nsph,lens->Nparticles);


  printf("\n***********************************************\n");
  printf("********* saving tree *************************\n\n");

  //saveTreeNB(lens->tree,lens->particles,lens->rsph,lens->treefilenames);
  saveSPHsmoothing(lens->tree,lens->particles,lens->rsph,lens->treefilenames);

  printf("tree saved to %s\n\n",lens->treefilenames);
  
  return 0;
}

