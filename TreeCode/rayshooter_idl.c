#include <math.h>
#include <stdio.h>
#include <time.h>
#include "../../Library/Recipes/nrutil.h"
//#include "../../Library/Recipes/nrutil.c"
#include <nr.h>
#include <nrD.h>
#include "../../Library/RecipesD/locateD.c"
#include "../../Library/RecipesD/odeintD.c"
#include "../../Library/RecipesD/bsstepD.c"
#include "../../Library/RecipesD/mmidD.c"
#include "../../Library/RecipesD/pzextrD.c"
#include "../../Library/RecipesD/polintD.c"
#include "../../Library/cosmo.h"
#include "../../Library/cosmo.c"

#include "TreeNB.h"
//#include "TreeNB.c"
//#include "TreeNBDriver.c"
//#include "TreeNBForce.c"
//#include "readfiles.c"

/*#include "minimization.c"*/
//#include "idl_export.h"


struct cosmology cosmo;
int kmax,kount; 
double *xp,**yp,dxsav;

void rayshooter(double *ray,IDL_LONG *Nrays,double *alpha,double *gamma
		,double *kappa,double *invmag,IDL_STRING *paramfile){
 
  int i;
  static SimLens lens;
  static int counter=0;
  static double convert_factor=0;

  for(i=0;i<*Nrays;++i){
    printf("%i ray = %e %e \n",i,ray[2*i],ray[2*i+1]);
    printf("ray 1= %e %e\n",(&ray[2*i])[0],(&ray[2*i])[1]);
  }

  ++counter;

  /*printf("%s  %s\n",teststring,paramfile[0].s);*/

  /* read in simulation data on first use */
  if(counter == 1){

    if(paramfile[0].slen == 0){ ERROR_MESSAGE(); print("ERROR: in rayshooter, no paramfile\n"); exit(0);}
    /* read parameter file and simulation */

    printf("READING PARAMETERS FROM %s\n\n",paramfile[0].s);
    readparams(paramfile[0].s,&lens,&cosmo);

    /*printf("READING PARAMETERS FROM %s\n\n",teststring);
      readparams(teststring,&lens,&cosmo);*/

    printf("\nreading particle positions from %s\n\n",lens.simfilename);
    readpositions(&lens);
    printf("*********    loaded    *********************\n\n",lens.treefilenames);
    printf("number of particles=%i\n",lens.Nparticles);

    /*for(i=0;i<Nparticles;++i) printf("%e %e %e\n",xp[i][0],xp[i][1],xp[i][2]);*/

    /*********************************/

    printf("\n***********************************************\n");
    printf("********* tree from %s *********************\n\n",lens.treefilenames);

    lens.tree=readTreeNB(lens.particles,lens.rsph,lens.Nparticles,lens.treefilenames);
    printf("*********    loaded    *********************\n\n",lens.treefilenames);

    convert_factor=lens.mass_units/lens.Sigma_crit;
  }

  /*** loop through rays ***/
  for(i=0;i<*Nrays;++i){

    /*printf("ray = %i     %e  %e\n",i,ray[2*i],ray[2*i+1]);*/

    /* [x,y] = x + x_size*y */
    TreeNBForce2D(lens.tree,lens.xp,lens.rsph,&ray[2*i],lens.coord,lens.theta,&alpha[2*i],&kappa[i],&gamma[2*i]);

    kappa[i]*=convert_factor;
    gamma[2*i]*=convert_factor;
    gamma[2*i+1]*=convert_factor;
    alpha[2*i]*=convert_factor;
    alpha[2*i+1]*=convert_factor;

    invmag[i]=pow(1-kappa[i],2) - pow(gamma[2*i],2) - pow(gamma[2*i+1],2);

    /*printf("   alpha = %e %e\n   gamma=%e %e\n   kappa = %e\n   invmag = %e\n"
	   ,alpha[2*i],alpha[2*i+1],gamma[2*i],gamma[2*i+1],kappa[i],invmag[i]);/**/
  }


}
