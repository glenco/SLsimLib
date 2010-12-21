/*
 * internal_rayshooter.c
 *
 *  Created on: Nov 2, 2009
 *      Author: R.B. Metcalf
 */
#include <stdio.h>
#include "../TreeCode_link/Tree.h"
#include "../TreeCode/TreeNB.h"


void rayshooterInternal(unsigned long Npoints,Point *i_points,TreeHndl i_tree,Boolean kappa_off){
//void rayshooterInternal(double *x,double *alpha,double *gamma,double *kappa,double *invmag){
  char paramfile[50];
  long i;

  sprintf(paramfile,"ParameterFiles/param_stars");
  rayshooter(Npoints,i_points,i_tree,paramfile);

 // rayshooter(x,1,alpha,gamma,kappa,invmag,paramfile);
}

void in_source(double *y_source,ListHndl sourcelist){
  return;
}
