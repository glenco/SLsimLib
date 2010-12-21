/*
 * internal_rayshooter.c
 *
 *  Created on: Sep 13, 2009
 *      Author: R.B. Metcalf
 */

#include <stdio.h>
#include "../TreeCode_link/Tree.h"
#include "../TreeCode/TreeNB.h"

extern char *paramfile;

void rayshooterInternal(unsigned long Npoints,Point *i_points,TreeHndl i_tree,Boolean kappa_off){
  /* i_points need to be already linked to s_points */
//  static char paramfile[50];
  static short i=0;

  if(i==0){
	  if(paramfile==NULL){
		  paramfile=(char *)malloc(60*sizeof(char));
		  //   	  sprintf(paramfile,"ParameterFiles/paramfileSIE1M_vol_512");
		  //	  sprintf(paramfile,"ParameterFiles/paramfileSIE1M_2");
		  //	  sprintf(paramfile,"ParameterFiles/paramfileSIE1M_vol");
		  sprintf(paramfile,"ParameterFiles/paramfileSIS1M_256");
		  //	  sprintf(paramfile,"ParameterFiles/paramfileSIE1M_vol_2D");
	  }
	  ++i;
  }
  /*rayshooter(x,1,alpha,gamma,kappa,invmag,paramfile);*/
  rayshooter(Npoints,i_points,i_tree,paramfile);
}

void in_source(double *y_source,ListHndl sourcelist){
  return;
}

