/*
 * internal_rayshooter_nfw.c
 *
 *  Created on: Dec 8, 2009
 *      Author: R.B. Metcalf
 */

#include "slsimlib.h"

/** \ingroup DeflectionL2
 *
 * \brief Routine for calculating the deflection and other lensing quantities for
 * a analytic one plane lens (AnaLens).
 *
 * Can be switched with rayshooterNB() to change
 * to particle lens model.  This transition needs to be made more automatic and
 * fail safe.
 */

#ifndef N_THREADS
#define N_THREADS 1
#endif

void *compute_rays_parallel_nfw(void *_p);

/**
 * A data structure for temporarily distribute the data amongst threads.
 */
struct TmpParams{
	Point *i_points;
	bool kappa_off;
	int tid;
	int start;
	int size;
	BaseAnaLens *lens;
};

void BaseAnaLens::rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off){

    if(this == NULL || !set){
    	ERROR_MESSAGE();
    	std::printf("ERROR: rayshooterInternal  lens not set!\n");
    	exit(0);
    }

    int nthreads, rc;

    nthreads = N_THREADS;

    int chunk_size;
    do{
      chunk_size = (int)Npoints/nthreads;
      if(chunk_size == 0) nthreads /= 2;
    }while(chunk_size == 0);

    pthread_t threads[nthreads];
    TmpParams *thread_params = new TmpParams[nthreads];

    for(int i=0; i<nthreads;i++){
      thread_params[i].i_points = i_points;
      thread_params[i].kappa_off = kappa_off;
      thread_params[i].size = chunk_size;
      if(i == nthreads-1)
        thread_params[i].size = (int)Npoints - (nthreads-1)*chunk_size;
      thread_params[i].start = i*chunk_size;
      thread_params[i].tid = i;
      thread_params[i].lens = this;
      rc = pthread_create(&threads[i], NULL, compute_rays_parallel_nfw, (void*) &thread_params[i]);
      assert(rc==0);
    }

    for(int i = 0; i < nthreads; i++){
      rc = pthread_join(threads[i], NULL);
      assert(rc==0);
    }

    delete[] thread_params;
}

void *compute_rays_parallel_nfw(void *_p){
	TmpParams *p = (TmpParams *) _p;
	bool kappa_off = p->kappa_off;
	BaseAnaLens *lens = p->lens;
	int tid        = p->tid;
	int chunk_size = p->size;
	int start      = p->start;
	int end        = start + chunk_size;

	/* i_points need to be already linked to s_points */

	double convert_factor = lens->star_massscale / lens->Sigma_crit;
	long i;

	double x_rescale[2];
    long j;
    float dt,kappa;
	double alpha[2];
	float gamma[3];

    for(i = start; i < end; i++){
    	alpha[0]=alpha[1]=0.0;
    	gamma[0]=gamma[1]=gamma[2]=0.0;
    	dt = kappa = 0.0;

    	p->i_points[i].dt = 0.0;
    	p->i_points[i].gamma[2] = 0.0;

    	// host lens
    	if(lens->host_ro > 0){
    		x_rescale[0] = p->i_points[i].x[0] / lens->host_ro;
    		x_rescale[1] = p->i_points[i].x[1] / lens->host_ro;

    		alphaNSIE(p->i_points[i].image->x, x_rescale, lens->host_axis_ratio,
    				lens->host_core/ lens->host_ro, lens->host_pos_angle);

    		if(!kappa_off){
    			gammaNSIE(p->i_points[i].gamma,x_rescale,lens->host_axis_ratio
    					,lens->host_core/lens->host_ro,lens->host_pos_angle);
    			p->i_points[i].kappa=kappaNSIE(x_rescale,lens->host_axis_ratio
    					,lens->host_core/lens->host_ro,lens->host_pos_angle);
    			p->i_points[i].dt = phiNSIE(x_rescale,lens->host_axis_ratio
    					,lens->host_core/lens->host_ro,lens->host_pos_angle);
    		}
    		else{
    			p->i_points[i].kappa=0;
    			p->i_points[i].gamma[0]=p->i_points[i].gamma[1]=0.0;
    			p->i_points[i].dt = 0.0;
    		}

    		p->i_points[i].image->x[0] *= lens->host_ro;
    		p->i_points[i].image->x[1] *= lens->host_ro;

    	}
    	else{
    		p->i_points[i].image->x[0] = 0.0;
    		p->i_points[i].image->x[1] = 0.0;
    		p->i_points[i].kappa=0;
    		p->i_points[i].gamma[0]=p->i_points[i].gamma[1]=0;
    		p->i_points[i].dt = 0.0;
    	}

    	// perturbations of host lens
    	if(lens->perturb_Nmodes > 0){
    		p->i_points[i].kappa += lens_expand(lens->perturb_beta,lens->perturb_modes
    				,lens->perturb_Nmodes,p->i_points[i].x,alpha,gamma,&dt);

    		p->i_points[i].image->x[0] += alpha[0];
    		p->i_points[i].image->x[1] += alpha[1];

    		if(!kappa_off){
    			p->i_points[i].gamma[0] += gamma[0];
    			p->i_points[i].gamma[1] += gamma[1];
    			p->i_points[i].dt += dt;
    		}
    		else
    			p->i_points[i].kappa = 0;
    	} // end of perturb modes

    	alpha[0]=alpha[1]=0.0;
    	gamma[0]=gamma[1]=gamma[2]=0.0;

    	// add substructure
    	if(lens->AreSubStructImaplated()){
    		for(j=0;j<lens->sub_N;++j){
    			lens->sub_alpha_func(alpha,p->i_points[i].x,lens->sub_Rcut[j],lens->sub_mass[j],lens->sub_beta,lens->sub_x[j],lens->Sigma_crit);

    			p->i_points[i].image->x[0] += alpha[0];
    			p->i_points[i].image->x[1] += alpha[1];

    			if(!kappa_off){
    				p->i_points[i].kappa += lens->sub_kappa_func(p->i_points[i].x,lens->sub_Rcut[j],lens->sub_mass[j],lens->sub_beta,lens->sub_x[j],lens->Sigma_crit);
    				lens->sub_gamma_func(gamma,p->i_points[i].x,lens->sub_Rcut[j],lens->sub_mass[j],lens->sub_beta,lens->sub_x[j],lens->Sigma_crit);
    				p->i_points[i].gamma[0] += gamma[0];
    				p->i_points[i].gamma[1] += gamma[1];
    				p->i_points[i].dt += lens->sub_phi_func(p->i_points[i].x,lens->sub_Rcut[j],lens->sub_mass[j],lens->sub_beta,lens->sub_x[j],lens->Sigma_crit);

    			}
    		}
    	} // end of substructure

    	if(!kappa_off){
    		p->i_points[i].dt = 0.5*(p->i_points[i].image->x[0]*p->i_points[i].image->x[0]
    		                     + p->i_points[i].image->x[1]*p->i_points[i].image->x[1])
      		                     - p->i_points[i].dt;
    		p->i_points[i].dt *= lens->to;
    	}

    	p->i_points[i].image->x[0] = p->i_points[i].x[0] - p->i_points[i].image->x[0];
    	p->i_points[i].image->x[1] = p->i_points[i].x[1] - p->i_points[i].image->x[1];

    	alpha[0]=alpha[1]=0.0;
    	gamma[0]=gamma[1]=gamma[2]=0.0;


    	if(lens->stars_N > 0 && lens->AreStarsImaplated()){
    		// add stars for microlensing
    		lens->substract_stars_disks(p->i_points[i].x,p->i_points[i].image->x,
    				&(p->i_points[i].kappa),p->i_points[i].gamma);

    		// do stars with tree code
    		lens->star_tree->force2D_recur(p->i_points[i].x,alpha,&kappa,gamma,true);

    		p->i_points[i].image->x[0] -= convert_factor*alpha[0];
    		p->i_points[i].image->x[1] -= convert_factor*alpha[1];

    		if(!kappa_off){
    			p->i_points[i].kappa += convert_factor*kappa;
    			p->i_points[i].gamma[0] += convert_factor*gamma[0];
    			p->i_points[i].gamma[1] += convert_factor*gamma[1];
    		}

    	} // end of stars


		// final operations to get the inverse magnification
		p->i_points[i].invmag = (1-p->i_points[i].kappa)*(1-p->i_points[i].kappa)
				- p->i_points[i].gamma[0]*p->i_points[i].gamma[0] - p->i_points[i].gamma[1]*p->i_points[i].gamma[1];

    	p->i_points[i].image->invmag=p->i_points[i].invmag;
    	p->i_points[i].image->kappa=p->i_points[i].kappa;
    	p->i_points[i].image->gamma[0]=p->i_points[i].gamma[0];
    	p->i_points[i].image->gamma[1]=p->i_points[i].gamma[1];
    	p->i_points[i].image->dt = p->i_points[i].dt;


    }

    return 0;
}


/** \ingroup DeflectionL2
   *
   * \brief Routine for calculating the deflection and other lensing quantities for
   * a analytic one plane lens (AnaLens), for just one ray!!
   *
*/
void BaseAnaLens::rayshooterInternal(double *ray, double *alpha, float *gamma, float *kappa, bool kappa_off){
     double x_rescale[2];
     long j;
     double alpha_tmp[2];
     float gamma_tmp[3], dt = 0,tmp = 0;

     gamma_tmp[0] = gamma_tmp[1] = gamma_tmp[2] = 0.0;
     alpha_tmp[0] = alpha_tmp[1] = 0.0;

     alpha[0] = alpha[1] = 0.0;
     gamma[0] = gamma[1] = gamma[2] = 0.0;
     *kappa = 0.0;

     double convert_factor = star_massscale/Sigma_crit;

     if(host_ro > 0){
    	 x_rescale[0] = ray[0]/host_ro;
    	 x_rescale[1] = ray[1]/host_ro;

    	 alphaNSIE(alpha,x_rescale,host_axis_ratio,host_core/host_ro,host_pos_angle);

    	 if(!kappa_off){
    		 gammaNSIE(gamma,x_rescale,host_axis_ratio,host_core/host_ro,host_pos_angle);
    		 *kappa=kappaNSIE(x_rescale,host_axis_ratio,host_core/host_ro,host_pos_angle);
    	 }

    	 alpha[0] *= host_ro;
    	 alpha[1] *= host_ro;
     }

  // perturbations of host lens
     if(perturb_Nmodes > 0){
    	 *kappa += lens_expand(perturb_beta,perturb_modes
    			 ,perturb_Nmodes,ray,alpha_tmp,gamma_tmp,&dt);

    	 alpha[0] += alpha_tmp[0];
    	 alpha[1] += alpha_tmp[1];

   	      if(!kappa_off){
   	    	  gamma[0] += gamma_tmp[0];
   	    	  gamma[1] += gamma_tmp[1];
   	      }

   	     gamma_tmp[0] = gamma_tmp[1] = gamma_tmp[2] = 0.0;
   	     alpha_tmp[0] = alpha_tmp[1] = 0.0;
     }

     // add substructure
     if(substruct_implanted){
    	 for(j=0;j<sub_N;++j){
    		 sub_alpha_func(alpha_tmp,ray,sub_Rcut[j],sub_mass[j],sub_beta,sub_x[j],Sigma_crit);

    		 alpha[0] += alpha_tmp[0];
    		 alpha[1] += alpha_tmp[1];

    		 if(!kappa_off){
    			 *kappa += sub_kappa_func(ray,sub_Rcut[j],sub_mass[j],sub_beta,sub_x[j],Sigma_crit);
    			 sub_gamma_func(gamma_tmp,ray,sub_Rcut[j],sub_mass[j],sub_beta,sub_x[j],Sigma_crit);

    			 gamma[0] += gamma_tmp[0];
    			 gamma[1] += gamma_tmp[1];
    		 }
    	 }

         gamma_tmp[0] = gamma_tmp[1] = gamma_tmp[2] = 0.0;
         alpha_tmp[0] = alpha_tmp[1] = 0.0;
     }

     // add stars for microlensing
     if(stars_N > 0 && stars_implanted){

    	 substract_stars_disks(ray,alpha,kappa,gamma);

    	 // do stars with tree code
    	 star_tree->force2D_recur(ray,alpha_tmp,&tmp,gamma_tmp,true);

    	 alpha[0] -= convert_factor*alpha_tmp[0];
    	 alpha[1] -= convert_factor*alpha_tmp[1];

    	 if(!kappa_off){
    		 *kappa += convert_factor*tmp;
    		 gamma[0] += convert_factor*gamma_tmp[0];
    		 gamma[1] += convert_factor*gamma_tmp[1];
    	 }
     }

     // final operations on results
     convert_factor = 4*pi*Grav*Sigma_crit;

     // convert from physical distance on the lens plane to an angle
	 alpha[0] *= convert_factor;
	 alpha[1] *= convert_factor;

	 // in the multi-plane formalism G^i=partial deflection_angle^i / partial x^i
	 // therefore the quantities need to be in units (1/physical_distance)
	 // --> convert from unitless quantity to (1/physical_distance)
	 *kappa *= convert_factor;
	 gamma[0] *= convert_factor;
	 gamma[1] *= convert_factor;
	 gamma[2] *= convert_factor;

     return ;
}
