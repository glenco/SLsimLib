/*
 * internal_rayshooter_analytic.c
 *
 *  Created on: Dec 8, 2009
 *      Author: R.B. Metcalf
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <assert.h>
#include "../../Library/Recipes/nr.h"
#include "../../Library/Recipes/nrutil.h"
#include "../../Library/RecipesD/nrD.h"

#include "../TreeCode_link/Tree.h"
#include "../../Library/cosmo.h"
#include "../AnalyticNSIE/analytic_lens.h"
#include "../TreeCode/TreeNB.h"

extern char *paramfile,*outputfile;
struct cosmology cosmo;
extern SimLens *lens;

/*
 * this version is for particles
 */
void rayshooterInternal(unsigned long Npoints,Point *i_points,TreeHndl i_tree,Boolean kappa_off){
  static char *paramfile=NULL;

  if(lens==NULL || !lens->set){
	  ERROR_MESSAGE();
	  printf("ERROR: rayshooterInternal  lens not set!\n");
	  exit(0);
  }

  if(paramfile==NULL){
	  paramfile=(char *)malloc(60*sizeof(char));
	  //sprintf(paramfile,"Example/paramfileParticles");
	  sprintf(paramfile,"Example/paramfileStefan");
  }

  rayshooterNB(Npoints,i_points,i_tree,paramfile,kappa_off);
}


// this is used for making non-circular sources
void in_source(double *y_source,ListHndl sourcelist){
  return;
}

