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
#include <nrutil.h>
#include <nrD.h>
#include <cosmo.h>

#include <Tree.h>
#include <point.h>
#include <analytic_lens.h>
#include <nsie.h>
#include <TreeNB.h>
#include <source_model.h>

//extern char *paramfile,*outputfile;
extern COSMOLOGY cosmo;
extern AnaLens *lens;

//const float Concentration=0.0776;  // ratio between truncation radius and scale length

void rayshooterInternal(unsigned long Npoints,Point *i_points,Boolean kappa_off){
  /* i_points need to be already linked to s_points */
  double x_rescale[2],alpha[2],gamma[2],tmp,dt=0;
  static double zs_old=-1,convert_factor=0;
  long i,j;
  //TreeNBHndl star_tree;

  if(lens==NULL || !lens->set){
	  ERROR_MESSAGE();
 	  printf("ERROR: rayshooterInternal  lens not set!\n");
 	  exit(0);
   }


  if( lens->zsource != zs_old){
	  lens->host_ro = 4*pi*pow(lens->host_sigma/2.99792e5,2)*angDist(0,lens->zlens,&cosmo)
    				*angDist(lens->zlens,lens->zsource,&cosmo)
    				/angDist(0,lens->zsource,&cosmo)/(1+lens->zlens);
	  zs_old=lens->zsource;
  }

  convert_factor = lens->star_massscale/lens->Sigma_crit;
  //star_tree = lens->star_tree;
  //#pragma omp parallel for private(x_rescale,star_tree,alpha,gamma,tmp)
  for(i=0;i<Npoints;++i){

	i_points[i].dt = 0;

    if(isnan(i_points[i].x[0]*i_points[i].x[1])){
		ERROR_MESSAGE();
		printf("x nan in internal_rayshooter\n    i=%li x= %e %e\n",
				i,i_points[i].x[0],i_points[i].x[1]);

	}
    if(lens->zsource <= lens->zlens){
      i_points[i].image->x[0]=i_points[i].x[0];
      i_points[i].image->x[1]=i_points[i].x[1];
      i_points[i].kappa=0.0;
      i_points[i].gamma[0]=0.0; i_points[i].gamma[1]=0.0;
      i_points[i].invmag=1.0;
      i_points[i].dt = 0;

    }else{

        // host lens
        if(lens->host_ro > 0){

        	x_rescale[0]=i_points[i].x[0]/lens->host_ro;
        	x_rescale[1]=i_points[i].x[1]/lens->host_ro;

        	alphaNSIE(i_points[i].image->x,x_rescale,lens->host_axis_ratio
    			  ,lens->host_core/lens->host_ro,lens->host_pos_angle);

        	if(!kappa_off){
        		gammaNSIE(i_points[i].gamma,x_rescale,lens->host_axis_ratio
    				  ,lens->host_core/lens->host_ro,lens->host_pos_angle);
        		i_points[i].kappa=kappaNSIE(x_rescale,lens->host_axis_ratio
    				  ,lens->host_core/lens->host_ro,lens->host_pos_angle);
        		i_points[i].dt = phiNSIE(x_rescale,lens->host_axis_ratio
      				  ,lens->host_core/lens->host_ro,lens->host_pos_angle);
        	}else{
        		i_points[i].kappa=0;
        		i_points[i].gamma[0]=i_points[i].gamma[1]=0.0;
        		i_points[i].dt = 0.0;
        	}

        	i_points[i].image->x[0] *= lens->host_ro;
        	i_points[i].image->x[1] *= lens->host_ro;

        }else{
        	i_points[i].image->x[0] = 0.0;
        	i_points[i].image->x[1] = 0.0;
        	i_points[i].kappa=0;
        	i_points[i].gamma[0]=i_points[i].gamma[1]=0;
       		i_points[i].dt = 0.0;
        }


          // perturbations of host lens
        if(lens->perturb_Nmodes > 0){
        	i_points[i].kappa += lens_expand(lens->perturb_beta,lens->perturb_modes
    			  ,lens->perturb_Nmodes,i_points[i].x,alpha,gamma,&dt);

        	i_points[i].image->x[0] += alpha[0];
        	i_points[i].image->x[1] += alpha[1];

        	if(!kappa_off){
        		i_points[i].gamma[0] += gamma[0];
        		i_points[i].gamma[1] += gamma[1];
        		i_points[i].dt += dt;
        	} else i_points[i].kappa = 0;
        }

          // add substructure
         if(lens->substruct_implanted){

  /*      	  if(lens->sub_tree){
        		  // do substructures with tree code
        		  TreeNBForce2D(lens->sub_tree,i_points[i].x,alpha,&tmp,gamma,True);

        		  //printf("alpha . r = %e substructure\n"
        		  //		  , sqrt(alpha[0]*alpha[0] + alpha[1]*alpha[1])/lens->Sigma_crit  );
        		  i_points[i].image->x[0] += alpha[0]/lens->Sigma_crit;
        		  i_points[i].image->x[1] += alpha[1]/lens->Sigma_crit;

        		  if(!kappa_off){
        			  i_points[i].kappa += tmp/lens->Sigma_crit;
        			  i_points[i].gamma[0] += gamma[0]/lens->Sigma_crit;
        			  i_points[i].gamma[1] += gamma[1]/lens->Sigma_crit;

        			  **** dt
        		  }
        	  }else{
*/
        	 for(j=0;j<lens->sub_N;++j){

        		 lens->sub_alpha_func(alpha,i_points[i].x,lens->sub_Rcut[j]
        			       ,lens->sub_mass[j],lens->sub_beta
        			       ,lens->sub_x[j],lens->Sigma_crit);

        		 assert(fabs(alpha[0]) > 0 && fabs(alpha[1]) >0 );

        		 //*** alphaNFW does not have the same input format

        		 //alphaPowLaw(alpha,i_points[i].x,lens->RcutSubstruct[j],lens->massSubstruct[j]
        		 //   ,Concentration*lens->RcutSubstruct[j],lens->xSubstruct[j],lens->Sigma_crit);
        		 //alphaNFW(alpha,i_points[i].x,lens->RcutSubstruct[j],lens->massSubstruct[j]
        		 //   ,Concentration*lens->RcutSubstruct[j],lens->xSubstruct[j],lens->Sigma_crit);
        		 i_points[i].image->x[0] += alpha[0];
        		 i_points[i].image->x[1] += alpha[1];

        		 if(!kappa_off){
        			 i_points[i].kappa += lens->sub_kappa_func(i_points[i].x,lens->sub_Rcut[j]
           		            ,lens->sub_mass[j],lens->sub_beta,lens->sub_x[j]
           		            ,lens->Sigma_crit);

        				  //i_points[i].kappa+=kappaNFW(i_points[i].x,lens->RcutSubstruct[j]
           		          //  ,lens->massSubstruct[j],Concentration*lens->RcutSubstruct[j],lens->xSubstruct[j],lens->Sigma_crit);

        			 lens->sub_gamma_func(gamma,i_points[i].x,lens->sub_Rcut[j],lens->sub_mass[j]
        			                      ,lens->sub_beta,lens->sub_x[j],lens->Sigma_crit);
        				  //gammaNFW(gamma,i_points[i].x,lens->RcutSubstruct[j],lens->massSubstruct[j]
        				  //        ,Concentration*lens->RcutSubstruct[j],lens->xSubstruct[j],lens->Sigma_crit);
        			 i_points[i].gamma[0] += gamma[0];
        			 i_points[i].gamma[1] += gamma[1];

        			 i_points[i].dt += lens->sub_phi_func(i_points[i].x,lens->sub_Rcut[j],lens->sub_mass[j]
        			                         ,lens->sub_beta,lens->sub_x[j],lens->Sigma_crit);
        			 //printf(" dt = %e\n",lens->sub_phi_func(i_points[i].x,lens->sub_Rcut[j],lens->sub_mass[j]
        			 //                        ,lens->sub_beta,lens->sub_x[j],lens->Sigma_crit));
        			 //if(isinf(i_points[i].dt)) exit(0);
        		 }
        	 }
          }

          // add stars for microlensing
          if(lens->stars_N > 0 && lens->stars_implanted){

        	  substract_stars_disks(lens,i_points[i].x,i_points[i].image->x,
        			  &(i_points[i].kappa),i_points[i].gamma);

        	  // do stars with tree code
        	  TreeNBForce2D(lens->star_tree,i_points[i].x,alpha,&tmp,gamma,True);

        	  i_points[i].image->x[0] += convert_factor*alpha[0];
        	  i_points[i].image->x[1] += convert_factor*alpha[1];
        	  //printf("alpha = %e %e\n",alpha[0]*convert_factor,alpha[1]*convert_factor);

        	  if(!kappa_off){
        		  i_points[i].kappa += convert_factor*tmp;
        		  i_points[i].gamma[0] += convert_factor*gamma[0];
        		  i_points[i].gamma[1] += convert_factor*gamma[1];
        	  }

          }

    	  if(!kappa_off){
    		  i_points[i].dt = 0.5*(i_points[i].image->x[0]*i_points[i].image->x[0]
                                    + i_points[i].image->x[1]*i_points[i].image->x[1])
                                 - i_points[i].dt;
    		  i_points[i].dt *= lens->to;
    	  }

          i_points[i].image->x[0] = i_points[i].x[0] - i_points[i].image->x[0];
          i_points[i].image->x[1] = i_points[i].x[1] - i_points[i].image->x[1];

          i_points[i].invmag = (1-i_points[i].kappa)*(1-i_points[i].kappa)
			- i_points[i].gamma[0]*i_points[i].gamma[0] - i_points[i].gamma[1]*i_points[i].gamma[1];
    }

    i_points[i].image->invmag=i_points[i].invmag;
    i_points[i].image->kappa=i_points[i].kappa;
    i_points[i].image->gamma[0]=i_points[i].gamma[0];
    i_points[i].image->gamma[1]=i_points[i].gamma[1];
    i_points[i].image->dt = i_points[i].dt;

/*    // register surface brightness at point
    x_rescale[0]=i_points[i].image->x[0] - lens->source_x[0];
    x_rescale[1]=i_points[i].image->x[1] - lens->source_x[1];

    j = (lens->source_r*lens->source_r > (x_rescale[0]*x_rescale[0]+x_rescale[1]*x_rescale[1]) );
    i_points[i].surface_brightness = j;
    if(j && lens->source_sb_type != Uniform) i_points[i].surface_brightness = lens->source_sb_func(x_rescale);
    */
   }


  return ;
}

double uniform_SB(double *y){
	return (double)( (y[0]*y[0] + y[1]*y[1]) < lens->source_r*lens->source_r );
}

double gaussian_SB(double *y){
	return exp( -(y[0]*y[0] + y[1]*y[1])/lens->source_gauss_r2 );
}

double BLR_SB(double *y){
	return blr_surface_brightness_disk(y,lens,&cosmo);
	return blr_surface_brightness_spherical_random_motions(sqrt(y[0]*y[0] + y[1]*y[1]),lens,&cosmo);
	return blr_surface_brightness_spherical_circular_motions(sqrt(y[0]*y[0] + y[1]*y[1]),lens,&cosmo);
}

void in_source(double *y_source,ListHndl sourcelist){
  return;
}
