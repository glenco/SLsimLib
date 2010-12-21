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
#include "analytic_lens.h"

extern char *paramfile,*outputfile;
struct cosmology cosmo;
extern AnaLens *lens;

void rayshooterInternal(unsigned long Npoints,Point *i_points,TreeHndl i_tree,Boolean kappa_off){
  /* i_points need to be already linked to s_points */
  double x_rescale[2],alpha[2],gamma[2],tmp;
  static double zs_old=-1;
  long i,j;
  static short init=1;
  //static AnaLens *lens;
  static long seed=-17283;

  if(lens==NULL || !lens->set){
	  printf("ERROR: rayshooterInternal  lens not set!\n");
	  exit(0);
  }

  /*
  if(init){

	  // read in lens parameters
	  if(paramfile==NULL){
		  paramfile=(char *)malloc(60*sizeof(char));
		  sprintf(paramfile,"AnalyticNSIE/paramfile");
	  }

	  lens=(AnaLens *)malloc(sizeof(AnaLens));
	  readparams_ana(paramfile,&cosmo,lens);
	  // construct random host for initialization of all lens variables
	  RandomizeHost(lens,1.0,&seed,True);
	  RandomizeSubstructure(lens,2,&seed);
	  //printf("Fraction within 2 Re: %e\n",FractionWithinRe(lens,2.0));

	  //PrintAnaLens(lens,False,False);
	  outputfile=lens->outputfile;
	  init=0;
  }
*/

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

      /*
       * For putting in a specific substructure
       *
      if(lens->NSubstruct==1){

    	  // position 1
    	  lens->xSubstruct[0][0]=1.5e-3;
    	  lens->xSubstruct[0][1]=-5.5e-3;

    	  // position 2
    	  lens->xSubstruct[0][0]=-3.0e-3;
    	  lens->xSubstruct[0][1]=-4.2e-3;

    	  // position 3
    	  lens->xSubstruct[0][0]=-2.5e-3*1.1;
    	  lens->xSubstruct[0][1]=4.0e-3*1.1;

    	  // position 4
    	  lens->xSubstruct[0][0]=0;
    	  lens->xSubstruct[0][1]=-6.0e-3;

    	  // position 5
    	  lens->xSubstruct[0][0]=-2.0e-3;
    	  lens->xSubstruct[0][1]=-5.25e-3;
      }
       */
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


          if(lens->Nmodes>0){
        	  i_points[i].kappa += lens_expand(lens->beta_perturb,lens->modes,lens->Nmodes,i_points[i].x,alpha,gamma);
        	  /*tmp=lens_expand(lens->beta_perturb,lens->modes,lens->Nmodes,i_points[i].x,alpha,gamma);

        	  if(alpha[0]/i_points[i].image->x[0]  < 0 ||
        	     alpha[1]/i_points[i].image->x[1]  < 0 ||
        	     gamma[0]/i_points[i].gamma[0]  < 0 ||
        		 gamma[1]/i_points[i].gamma[1]  < 0 ){
        	    printf("%e %e  %e  %e %e\n",i_points[i].image->x[0],i_points[i].image->x[1]
        	              ,i_points[i].kappa,i_points[i].gamma[0],i_points[i].gamma[1]);
        	 	  printf("kappa =%e alpha = %e %e gamma =%e %e\n",tmp,alpha[0],alpha[1]
        	        	                                     ,gamma[0],gamma[1]);
        	 	  //exit(0);
        	  }
	*/
        	  i_points[i].image->x[0] += alpha[0];
           	  i_points[i].image->x[1] += alpha[1];
           	  if(!kappa_off){
            	  i_points[i].gamma[0] += gamma[0];
            	  i_points[i].gamma[1] += gamma[1];
           	  } else i_points[i].kappa = 0;
         }

          // add substructure

    	  //PrintAnaLens(lens,False,False);
          for(j=0;j<lens->NSubstruct;++j){
//        	  printf("%e %e %e  %e %e\n",lens->RcutSubstruct[j],lens->massSubstruct[j]
//        	             ,lens->betaSubstruct,lens->xSubstruct[j][0],lens->xSubstruct[j][1]);

        	  alphaPowLaw(alpha,i_points[i].x,lens->RcutSubstruct[j],lens->massSubstruct[j]
        	             ,lens->betaSubstruct,lens->xSubstruct[j],lens->Sigma_crit);
        	  i_points[i].image->x[0] += alpha[0];
        	  i_points[i].image->x[1] += alpha[1];

           	  if(!kappa_off){
           		  i_points[i].kappa+=kappaPowLaw(i_points[i].x,lens->RcutSubstruct[j],lens->massSubstruct[j]
        	          ,lens->betaSubstruct,lens->xSubstruct[j],lens->Sigma_crit);
           		  gammaPowLaw(gamma,i_points[i].x,lens->RcutSubstruct[j],lens->massSubstruct[j]
        	              ,lens->betaSubstruct,lens->xSubstruct[j],lens->Sigma_crit);
           		  i_points[i].gamma[0] += gamma[0];
           		  i_points[i].gamma[1] += gamma[1];
           	  }
          }
      }else{
    	  i_points[i].kappa=0;
    	  i_points[i].gamma[0]=i_points[i].gamma[1]=0;
      }

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


void in_source(double *y_source,ListHndl sourcelist){
  return;
}

