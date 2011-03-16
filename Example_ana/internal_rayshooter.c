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
#include <nr.h>
#include "../../Library/Recipes/nrutil.h"
#include <nrD.h>

#include "../TreeCode_link/Tree.h"
#include "../../Library/cosmo.h"
#include "../AnalyticNSIE/analytic_lens.h"
#include "../TreeCode/TreeNB.h"

extern char *paramfile,*outputfile;
struct cosmology cosmo;
extern AnaLens *lens;

/*
 * this function is for calculating the deflection using an analytic SIE model
 */

void rayshooterInternal(unsigned long Npoints,Point *i_points,TreeHndl i_tree,Boolean kappa_off){
  /* i_points need to be already linked to s_points */
  static double zs_old=-1;
  long i;
  static short init=1;
  double x_rescale[2];

  if(lens==NULL || !lens->set){
	  ERROR_MESSAGE();
	  printf("ERROR: rayshooterInternal  lens not set!\n");
	  exit(0);
  }

/*  if(init){

	  // read in lens parameters
	  if(paramfile==NULL){
		  paramfile=(char *)malloc(60*sizeof(char));
		  sprintf(paramfile,"Example/paramfileAna");
	  }

	  lens=(AnaLens *)malloc(sizeof(AnaLens));
	  readparams_ana(paramfile,&cosmo,lens);

	  //PrintAnaLens(lens,False,False);
	  outputfile=lens->outputfile;



	  init=0;
  }*/

  if( lens->zsource != zs_old){
	  lens->ro=4*pi*pow(lens->sigma/2.99792e5,2)*angDist(0,lens->zlens)
    				*angDist(lens->zlens,lens->zsource)/angDist(0,lens->zsource)/(1+lens->zlens);
	  zs_old=lens->zsource;
  }

//#pragma omp parallel for private(x_rescale)
  for(i=0;i<Npoints;++i){

	if(isnan(i_points[i].x[0]*i_points[i].x[1])){
		printf("x nan in internal_rayshooter\n    i=%li x= %e %e\n",
				i,i_points[i].x[0],i_points[i].x[1]);

	}

	// if source is in front of lens no lensing
    if(lens->zsource <= lens->zlens){
      i_points[i].image->x[0]=i_points[i].x[0];
      i_points[i].image->x[1]=i_points[i].x[1];
      i_points[i].kappa=0.0;
      i_points[i].gamma[0]=0.0; i_points[i].gamma[1]=0.0;
      i_points[i].invmag=1.0;

    }else{

    	//if(lens->ro < 1.0e-10 || isnan(lens->ro))

      x_rescale[0]=i_points[i].x[0]/lens->ro;
      x_rescale[1]=i_points[i].x[1]/lens->ro;

      ////////////////////////////////////////////////////////////

      if(lens->ro > 0){
    	  //PrintAnaLens(lens,False,False);
    	  alphaNSIE(i_points[i].image->x,x_rescale,lens->axis_ratio,lens->core/lens->ro,lens->theta);
    	  if(!kappa_off){
    		  gammaNSIE(i_points[i].gamma,x_rescale,lens->axis_ratio,lens->core/lens->ro,lens->theta);
    		  i_points[i].kappa=kappaNSIE(x_rescale,lens->axis_ratio,lens->core/lens->ro,lens->theta);
    	  }else{
    		  i_points[i].kappa=0;
    		  i_points[i].gamma[0]=i_points[i].gamma[1]=0.0;
    	  }

          i_points[i].image->x[0] *= lens->ro;
          i_points[i].image->x[1] *= lens->ro;

      }else{
    	  i_points[i].kappa=0;
    	  i_points[i].gamma[0]=i_points[i].gamma[1]=0;
      }

      //rayshooter(Npoints,i_points,i_tree,paramfile);

      i_points[i].image->x[0]=i_points[i].x[0] - i_points[i].image->x[0];
      i_points[i].image->x[1]=i_points[i].x[1] - i_points[i].image->x[1];

      i_points[i].invmag=(1-i_points[i].kappa)*(1-i_points[i].kappa)
			- i_points[i].gamma[0]*i_points[i].gamma[0] - i_points[i].gamma[1]*i_points[i].gamma[1];
    }

    i_points[i].image->invmag=i_points[i].invmag;
    i_points[i].image->kappa=i_points[i].kappa;
    i_points[i].image->gamma[0]=i_points[i].gamma[0];
    i_points[i].image->gamma[1]=i_points[i].gamma[1];
  }

  return ;
}

// this is used for making non-circular sources
void in_source(double *y_source,ListHndl sourcelist){
  return;
}

