/*
 * internal_rayshooter_nfw.c
 *
 *  Created on: Dec 8, 2009
 *      Author: R.B. Metcalf
 */

#include <slsimlib.h>

extern CosmoHndl cosmo;
extern AnaLens *lens;

struct temp_data
{
  double alpha[2], gamma[2];
}
  *temp;

/** \ingroup DeflectionL2
 *
 * \brief Routine for calculating the deflection and other lensing quantities for
 * a analytic one plane lens (AnaLens).
 *
 * Can be switched with rayshooterNB() to change
 * to particle lens model.  This transition needs to be made more automatic and
 * fail safe.
 */
  void AnaLens::rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off){
    /* i_points need to be already linked to s_points */
    double x_rescale[2], tmp, dt = 0;
    static double zs_old=-1,convert_factor=0;
    double convert_fac=0;
    long i,j;

    temp = (struct temp_data *) malloc(Npoints * sizeof(struct temp_data));

    if(lens == NULL || !set)
      {
        ERROR_MESSAGE();
        std::printf("ERROR: rayshooterInternal  lens not set!\n");
        exit(0);
      }

    if(zsource != zs_old)
      {
        host_ro = 4*pi*pow(host_sigma/2.99792e5,2) * cosmo->angDist(0,zlens)
  	*cosmo->angDist(zlens,zsource)
  	/cosmo->angDist(0,zsource)/(1+zlens);

        zs_old=zsource;
      }

    convert_factor = star_massscale / Sigma_crit;

  #pragma omp parallel for private(x_rescale, dt, i)
    for(i = 0; i < Npoints; i++)
      {
        i_points[i].dt = 0;

        if(isnan(i_points[i].x[0]*i_points[i].x[1]))
  	{
  	  ERROR_MESSAGE();
  	  std::printf("x nan in internal_rayshooter\n    i=%li x= %e %e\n",
  		 i,i_points[i].x[0],i_points[i].x[1]);
  	}

        if(zsource <= zlens)
  	{
  	  i_points[i].image->x[0]=i_points[i].x[0];
  	  i_points[i].image->x[1]=i_points[i].x[1];
  	  i_points[i].kappa=0.0;
  	  i_points[i].gamma[0]=0.0; i_points[i].gamma[1]=0.0;
  	  i_points[i].invmag=1.0;
  	  i_points[i].dt = 0;
  	}
        else
  	{
  	  // host lens
  	  if(host_ro > 0)
  	    {
  	      x_rescale[0] = i_points[i].x[0] / host_ro;
  	      x_rescale[1] = i_points[i].x[1] / host_ro;

  	      alphaNSIE(i_points[i].image->x, x_rescale, host_axis_ratio,
  			host_core / host_ro, host_pos_angle);

  	      if(!kappa_off)
  		{
  		  gammaNSIE(i_points[i].gamma,x_rescale,host_axis_ratio
  			    ,host_core/host_ro,host_pos_angle);
  		  i_points[i].kappa=kappaNSIE(x_rescale,host_axis_ratio
  					      ,host_core/host_ro,host_pos_angle);
  		  i_points[i].dt = phiNSIE(x_rescale,host_axis_ratio
  					   ,host_core/host_ro,host_pos_angle);
  		}
  	      else
  		{
  		  i_points[i].kappa=0;
  		  i_points[i].gamma[0]=i_points[i].gamma[1]=0.0;
  		  i_points[i].dt = 0.0;
  		}

  	      i_points[i].image->x[0] *= host_ro;
  	      i_points[i].image->x[1] *= host_ro;

  	    }
  	  else
  	    {
  	      i_points[i].image->x[0] = 0.0;
  	      i_points[i].image->x[1] = 0.0;
  	      i_points[i].kappa=0;
  	      i_points[i].gamma[0]=i_points[i].gamma[1]=0;
  	      i_points[i].dt = 0.0;
  	    }


            // perturbations of host lens
  	  if(perturb_Nmodes > 0)
  	    {
  	      i_points[i].kappa += lens_expand(perturb_beta,perturb_modes
  					       ,perturb_Nmodes,i_points[i].x,temp[i].alpha,temp[i].gamma,&dt);

  	      i_points[i].image->x[0] += temp[i].alpha[0];
  	      i_points[i].image->x[1] += temp[i].alpha[1];

  	      if(!kappa_off)
  		{
  		  i_points[i].gamma[0] += temp[i].gamma[0];
  		  i_points[i].gamma[1] += temp[i].gamma[1];
  		  i_points[i].dt += dt;
  		}
  	      else
  		i_points[i].kappa = 0;
  	    }

            // add substructure
           if(substruct_implanted)
  	   {
  	     for(j=0;j<sub_N;++j)
  	       {
  		 sub_alpha_func(temp[i].alpha,i_points[i].x,sub_Rcut[j]
  				      ,sub_mass[j],sub_beta
  				      ,sub_x[j],Sigma_crit);

  		 //assert(fabs(alpha[0]) > 0 && fabs(alpha[1]) >0 );

  		 //*** alphaNFW does not have the same input format

  		 i_points[i].image->x[0] += temp[i].alpha[0];
  		 i_points[i].image->x[1] += temp[i].alpha[1];

  		 if(!kappa_off)
  		   {
  		     i_points[i].kappa += sub_kappa_func(i_points[i].x,sub_Rcut[j]
  							       ,sub_mass[j],sub_beta,sub_x[j]
  							       ,Sigma_crit);


  		     sub_gamma_func(temp[i].gamma,i_points[i].x,sub_Rcut[j],sub_mass[j]
  					  ,sub_beta,sub_x[j],Sigma_crit);

  		     i_points[i].gamma[0] += temp[i].gamma[0];
  		     i_points[i].gamma[1] += temp[i].gamma[1];

  		     i_points[i].dt += sub_phi_func(i_points[i].x,sub_Rcut[j],sub_mass[j]
  							  ,sub_beta,sub_x[j],Sigma_crit);

  		   }
  	       }
  	   }



      	  if(!kappa_off)
  	    {
  	      i_points[i].dt = 0.5*(i_points[i].image->x[0]*i_points[i].image->x[0]
                                      + i_points[i].image->x[1]*i_points[i].image->x[1])
  		- i_points[i].dt;
  	      i_points[i].dt *= to;
  	    }

            i_points[i].image->x[0] = i_points[i].x[0] - i_points[i].image->x[0];
            i_points[i].image->x[1] = i_points[i].x[1] - i_points[i].image->x[1];

            i_points[i].invmag = (1-i_points[i].kappa)*(1-i_points[i].kappa)
  	    - i_points[i].gamma[0]*i_points[i].gamma[0] - i_points[i].gamma[1]*i_points[i].gamma[1];
  	}
      }

    for(i = 0; i< Npoints; i++)
      {
  	 // add stars for microlensing
        if(stars_N > 0 && stars_implanted)
  	{
  	  substract_stars_disks(lens,i_points[i].x,i_points[i].image->x,
  				&(i_points[i].kappa),i_points[i].gamma);

  	  // do stars with tree code
  	  star_tree->force2D(i_points[i].x,temp[i].alpha,&tmp,temp[i].gamma,true);

  	  i_points[i].image->x[0] += convert_factor*temp[i].alpha[0];
  	  i_points[i].image->x[1] += convert_factor*temp[i].alpha[1];

  	  if(!kappa_off)
  	    {
  		 i_points[i].kappa += convert_factor*tmp;
  		 i_points[i].gamma[0] += convert_factor*temp[i].gamma[0];
  		 i_points[i].gamma[1] += convert_factor*temp[i].gamma[1];
  	    }

  	}

      }

  #pragma omp parallel for private(i)
      for(i = 0; i < Npoints; i++)
        {
  	i_points[i].image->invmag=i_points[i].invmag;
  	i_points[i].image->kappa=i_points[i].kappa;
  	i_points[i].image->gamma[0]=i_points[i].gamma[0];
  	i_points[i].image->gamma[1]=i_points[i].gamma[1];
  	i_points[i].image->dt = i_points[i].dt;
        }

      free(temp);

    return ;
  }

double uniform_SB(double *y){
	return (double)( (y[0]*y[0] + y[1]*y[1]) < lens->source_r*lens->source_r );
}

double gaussian_SB(double *y){
	return exp( -(y[0]*y[0] + y[1]*y[1])/lens->source_gauss_r2 );
}

// surface brightness for models of the Broad Line Region
double BLR_Disk_SB(double *y){
	return blr_surface_brightness_disk(y,lens,cosmo);
}

double BLR_Sph1_SB(double *y){
	return blr_surface_brightness_spherical_circular_motions(sqrt(y[0]*y[0] + y[1]*y[1]),lens,cosmo);
}
double BLR_Sph2_SB(double *y){
	return blr_surface_brightness_spherical_random_motions(sqrt(y[0]*y[0] + y[1]*y[1]),lens,cosmo);
}

void in_source(double *y_source,ListHndl sourcelist){
  return;
}
