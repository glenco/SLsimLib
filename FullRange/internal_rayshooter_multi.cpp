/*
 * internal_rayshooter_multi.cpp
 *
 *  Created on: Jan 19, 2012
 *      Author: mpetkova, bmetcalf
 */

#include "slsimlib.h"

/** \ingroup DeflectionL2
 *
 * \brief This is the function that does the deflection calculation with multiple lens planes.
 *
 * The convergence, shear and rotation will be calculated if kappa_off == false
 *
 * Needs to be calculated before calling:
 *
 *
 * Dl[j = 0...Nplanes-1] - The angular size distance between the observer and the jth plane not counting the observer plane.
 * 	                    Dl[0] is the first plane with mass on it and Dl[Nplane-1] is the distance to the source plane.
 *
 * dDl[j = 0...Nplanes-1] - The angular size distance between the (j-1)th and jth planes counting the observer plane as j = -1.
 *                      dDl[0] = Dl[0], dDl[Nplane-1] is between the last plane with mass and the source plane.
 *
 * charge = 4*G*mass_scale/c^2 in units of Mpc
 *
 * i_points[].x[] is in angular units.
 *
 * Warning: Is not valid for a non-flat universe.
 */

#ifndef N_THREADS
#define N_THREADS 1
#endif

void *compute_rays_parallel(void *_p);

/**
 * \brief A data structure for temporarily distribute the data amongst threads.
 */
struct TmpParams{
  Point *i_points;
  bool kappa_off;
  int tid;
  int start;
  int size;
  int NPlanes;
  bool flag_switch_deflection_off;
  double charge;
  LensPlane** lensing_planes;
  double* plane_redshifts;
  double* Dl;
  double* dDl;
};

void Lens::rayshooterInternal(
		unsigned long Npoints   /// number of points to be shot
		,Point *i_points        /// point on the image plane
		,bool kappa_off         /// turns calculation of convergence and shear off to save time.
		){
    
    kappa_off =false;
    
  int NLastPlane;
  double tmpDs,tmpdDs,tmpZs;

  if(Npoints == 0) return;

// If a lower redshift source is being used
  if(toggle_source_plane){
    NLastPlane = index_of_new_sourceplane;

    assert(NLastPlane <= lensing_planes.size());
    tmpDs = Dl[index_of_new_sourceplane];
    tmpdDs = dDl[index_of_new_sourceplane];
    tmpZs = plane_redshifts[index_of_new_sourceplane];

    Dl[index_of_new_sourceplane] = Ds_implant;
    dDl[index_of_new_sourceplane] = dDs_implant;
    plane_redshifts[index_of_new_sourceplane] = zs_implant;
  }else{
    NLastPlane = lensing_planes.size();
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
    thread_params[i].flag_switch_deflection_off = flag_switch_deflection_off;
    thread_params[i].charge = charge;
    thread_params[i].lensing_planes = &lensing_planes[0];
    thread_params[i].plane_redshifts = &plane_redshifts[0];
    thread_params[i].Dl = &Dl[0];
    thread_params[i].dDl = &dDl[0];
    thread_params[i].NPlanes = NLastPlane;
    rc = pthread_create(&threads[i], NULL, compute_rays_parallel, (void*) &thread_params[i]);
    assert(rc==0);
  }
    
  for(int i = 0; i < nthreads; i++){
    rc = pthread_join(threads[i], NULL);
    assert(rc==0);
  }

  delete[] thread_params;

  if(toggle_source_plane){
    // The initial values for the plane are reset here
    Dl[index_of_new_sourceplane] = tmpDs;
    dDl[index_of_new_sourceplane] = tmpdDs;
    plane_redshifts[index_of_new_sourceplane] = tmpZs;
  }

}

void *compute_rays_parallel(void *_p)
{
  TmpParams *p = (TmpParams *) _p;
  bool kappa_off = p->kappa_off;
  int chunk_size = p->size;
  int start      = p->start;
  int end        = start + chunk_size;
    
  int i, j;
  
  double xx[2],fac;
  double aa,bb,cc;
  double alpha[2];
  KappaType kappa,gamma[3];
  double xminus[2],xplus[2];
  double kappa_minus,gamma_minus[3],kappa_plus,gamma_plus[3];
  
  for(i = start; i< end; i++){
    
	  if(p->i_points[i].in_image == MAYBE)
			 continue;

    // find position on first lens plane in comoving units
    p->i_points[i].image->x[0] = p->i_points[i].x[0]*p->Dl[0];
    p->i_points[i].image->x[1] = p->i_points[i].x[1]*p->Dl[0];
    
    xminus[0] = 0;
    xminus[1] = 0;
    
    // Set magnification matrix on first plane.  Also the default if kappa_off == false
    kappa_minus = 0;
    gamma_minus[0] = 0;
    gamma_minus[1] = 0;
    gamma_minus[2] = 0;
    
    p->i_points[i].kappa = 1;  // This is actually 1-kappa until after the loop through the planes.
    p->i_points[i].gamma[0] = 0;
    p->i_points[i].gamma[1] = 0;
    p->i_points[i].gamma[2] = 0;
    
    for(j = 0; j < p->NPlanes ; ++j){  // each iteration leaves i_point[i].image on plane (j+1)
      
      // convert to physical coordinates on the plane j
      
      xx[0] = p->i_points[i].image->x[0]/(1+p->plane_redshifts[j]);
      xx[1] = p->i_points[i].image->x[1]/(1+p->plane_redshifts[j]);
      
      assert(xx[0] == xx[0] && xx[1] == xx[1]);

      p->lensing_planes[j]->force(alpha,&kappa,gamma,xx,kappa_off);

      cc = p->charge*p->dDl[j+1];

      assert(alpha[0] == alpha[0] && alpha[1] == alpha[1]);

      if(!kappa_off){
    	  fac = 1/(1+p->plane_redshifts[j]);
    	  /* multiply by fac to obtain 1/comoving_distance/physical_distance
    	   * such that a multiplication with the charge (in units of physical distance)
    	   * will result in a 1/comoving_distance quantity */
    	  kappa*=fac;
    	  gamma[0]*=fac;
    	  gamma[1]*=fac;
    	  gamma[2]*=fac;
	
    	  assert(gamma[0] == gamma[0] && gamma[1] == gamma[1]);
    	  assert(kappa == kappa);
      }
      
      if(p->flag_switch_deflection_off)
    	  alpha[0] = alpha[1] = 0.0;
            
      aa = (p->dDl[j+1]+p->dDl[j])/p->dDl[j];
      bb = p->dDl[j+1]/p->dDl[j];
      
      xplus[0] = aa*p->i_points[i].image->x[0] - bb*xminus[0] - cc*alpha[0];
      xplus[1] = aa*p->i_points[i].image->x[1] - bb*xminus[1] - cc*alpha[1];
      
      xminus[0] = p->i_points[i].image->x[0];
      xminus[1] = p->i_points[i].image->x[1];
      
      p->i_points[i].image->x[0] = xplus[0];
      p->i_points[i].image->x[1] = xplus[1];
      
      if(!kappa_off){
	
    	  aa = (p->dDl[j+1]+p->dDl[j])*p->Dl[j]/p->dDl[j]/p->Dl[j+1];
	
    	  if(j>0){
    		  bb = p->dDl[j+1]*p->Dl[j-1]/p->dDl[j]/p->Dl[j+1];
    	  }
    	  else
    		  bb = 0;
	

    	  cc = p->charge*p->dDl[j+1]*p->Dl[j]/p->Dl[j+1];
	
    	  kappa_plus = aa*p->i_points[i].kappa - bb*kappa_minus
    			  - cc*(kappa*p->i_points[i].kappa + gamma[0]*p->i_points[i].gamma[0] + gamma[1]*p->i_points[i].gamma[1]);
	
    	  gamma_plus[0] = aa*p->i_points[i].gamma[0] - bb*gamma_minus[0]
    	          - cc*(gamma[0]*p->i_points[i].kappa + kappa*p->i_points[i].gamma[0] - gamma[1]*p->i_points[i].gamma[2]);
	
    	  gamma_plus[1] = aa*p->i_points[i].gamma[1] - bb*gamma_minus[1]
    	          - cc*(gamma[1]*p->i_points[i].kappa + kappa*p->i_points[i].gamma[1] + gamma[0]*p->i_points[i].gamma[2]);
	
    	  gamma_plus[2] = aa*p->i_points[i].gamma[2] - bb*gamma_minus[2]
    	          - cc*(kappa*p->i_points[i].gamma[2] - gamma[1]*p->i_points[i].gamma[0] + gamma[0]*p->i_points[i].gamma[1]);
	
    	  kappa_minus = p->i_points[i].kappa;
    	  gamma_minus[0] = p->i_points[i].gamma[0];
    	  gamma_minus[1] = p->i_points[i].gamma[1];
    	  gamma_minus[2] = p->i_points[i].gamma[2];
	
    	  assert(kappa_plus==kappa_plus && gamma_minus[0]==gamma_minus[0] && gamma_minus[1]==gamma_minus[1] && gamma_minus[2]==gamma_minus[2]);
	
    	  p->i_points[i].kappa = kappa_plus;
    	  p->i_points[i].gamma[0] = gamma_plus[0];
    	  p->i_points[i].gamma[1] = gamma_plus[1];
    	  p->i_points[i].gamma[2] = gamma_plus[2];
	
      }
    }
    
    // Convert units back to angles.
    p->i_points[i].image->x[0] /= p->Dl[p->NPlanes];
    p->i_points[i].image->x[1] /= p->Dl[p->NPlanes];
    
    p->i_points[i].kappa = 1 - p->i_points[i].kappa;

    if(!kappa_off) p->i_points[i].invmag = (1-p->i_points[i].kappa)*(1-p->i_points[i].kappa)
		     - p->i_points[i].gamma[0]*p->i_points[i].gamma[0]
		     - p->i_points[i].gamma[1]*p->i_points[i].gamma[1]
		     + p->i_points[i].gamma[2]*p->i_points[i].gamma[2];
    else p->i_points[i].invmag = 0.0;

    p->i_points[i].image->invmag=p->i_points[i].invmag;
    p->i_points[i].image->kappa=p->i_points[i].kappa;
    p->i_points[i].image->gamma[0]=p->i_points[i].gamma[0];
    p->i_points[i].image->gamma[1]=p->i_points[i].gamma[1];
    p->i_points[i].image->gamma[2]=p->i_points[i].gamma[2];

    if(p->i_points[i].image->x[0] != p->i_points[i].image->x[0] ||
       p->i_points[i].image->x[1] != p->i_points[i].image->x[1] ||
       p->i_points[i].invmag != p->i_points[i].invmag){
      ERROR_MESSAGE();
      std::cout << p->i_points[i].image->x[0] << "  " << p->i_points[i].image->x[1] << "  " << p->i_points[i].invmag << std::endl;
      std::cout << p->i_points[i].gamma[0] << "  " << p->i_points[i].gamma[1] << "  " << p->i_points[i].gamma[2] << "  " <<
    		  p->i_points[i].kappa << "  "  << kappa_off << std::endl;
      //	assert(0);
      exit(1);
    }
    
  }
  
  return 0;
  
}

/**
 * Calculates the source plane position x_source, gamma, and kappa for a single ray (in radians)
 *
 */
void Lens::rayshooterInternal(
		double *ray 	  /// position on the image plane in radians
		,double *x_source /// position on the source plane in radians
		,KappaType *gamma
		,KappaType *kappa
		,bool kappa_off		/// if true -- no kappa and gamma values are calculated
		){

    kappa_off = false;
    
  std::size_t j;
  std::size_t NPlanes = lensing_planes.size();

  double xx[2],fac;
  double aa,bb,cc;
  double alpha_tmp[2];
  KappaType kappa_tmp,gamma_tmp[3];
  double xminus[2],xplus[2];
  double kappa_tmp_minus,gamma_tmp_minus[3],kappa_tmp_plus,gamma_tmp_plus[3];

  assert(NPlanes > 0);
  // If a lower redshift source is being used

    // find position on first lens plane in comoving units
    x_source[0] = ray[0]*Dl[0];
    x_source[1] = ray[1]*Dl[0];

    xminus[0] = 0;
    xminus[1] = 0;

    // Set magnification matrix on first plane.  Also the default if kappa_off == false
    kappa_tmp_minus = 0;
    gamma_tmp_minus[0] = 0;
    gamma_tmp_minus[1] = 0;
    gamma_tmp_minus[2] = 0;

    *kappa = 1;  // This is actually 1-kappa_tmp until after the loop through the planes.
    gamma[0] = 0;
    gamma[1] = 0;
    gamma[2] = 0;

    for(j = 0; j < NPlanes; ++j){  // each iteration leaves i_point[i].image on plane (j+1)

      // convert to physical coordinates on the plane j

      xx[0] = x_source[0]/(1+plane_redshifts[j]);
      xx[1] = x_source[1]/(1+plane_redshifts[j]);

      assert(xx[0] == xx[0] && xx[1] == xx[1]);

      lensing_planes[j]->force(alpha_tmp,&kappa_tmp,gamma_tmp,xx,kappa_off);

      assert(alpha_tmp[0] == alpha_tmp[0] && alpha_tmp[1] == alpha_tmp[1]);

      cc = charge*dDl[j+1];

      if(!kappa_off){
    	  fac = 1/(1+plane_redshifts[j]);
    	  /* multiply by fac to obtain 1/comoving_distance/physical_distance
    	   * such that a multiplication with the charge (in units of physical distance)
    	   * will result in a 1/comoving_distance quantity */
    	  kappa_tmp*=fac;
    	  gamma_tmp[0]*=fac;
    	  gamma_tmp[1]*=fac;
    	  gamma_tmp[2]*=fac;

    	  assert(gamma_tmp[0] == gamma_tmp[0] && gamma_tmp[1] == gamma_tmp[1]);
    	  assert(kappa_tmp == kappa_tmp);
      }

      if(flag_switch_deflection_off)
    	  alpha_tmp[0] = alpha_tmp[1] = 0.0;

      aa = (dDl[j+1]+dDl[j])/dDl[j];
      bb = dDl[j+1]/dDl[j];

      xplus[0] = aa*x_source[0] - bb*xminus[0] - cc*alpha_tmp[0];
      xplus[1] = aa*x_source[1] - bb*xminus[1] - cc*alpha_tmp[1];

      xminus[0] = x_source[0];
      xminus[1] = x_source[1];

      x_source[0] = xplus[0];
      x_source[1] = xplus[1];

      if(!kappa_off){

    	  aa = (dDl[j+1]+dDl[j])*Dl[j]/dDl[j]/Dl[j+1];

    	  if(j>0){
    		  bb = dDl[j+1]*Dl[j-1]/dDl[j]/Dl[j+1];
    	  }
    	  else
    		  bb = 0;

    	  cc = charge*dDl[j+1]*Dl[j]/Dl[j+1];

    	  kappa_tmp_plus = aa*(*kappa) - bb*kappa_tmp_minus
    			  - cc*(kappa_tmp*(*kappa) + gamma_tmp[0]*gamma[0] + gamma_tmp[1]*gamma[1]);

    	  gamma_tmp_plus[0] = aa*gamma[0] - bb*gamma_tmp_minus[0]
    	          - cc*(gamma_tmp[0]*(*kappa) + kappa_tmp*gamma[0] - gamma_tmp[1]*gamma[2]);

    	  gamma_tmp_plus[1] = aa*gamma[1] - bb*gamma_tmp_minus[1]
    	          - cc*(gamma_tmp[1]*(*kappa) + kappa_tmp*gamma[1] + gamma_tmp[0]*gamma[2]);

    	  gamma_tmp_plus[2] = aa*gamma[2] - bb*gamma_tmp_minus[2]
    	          - cc*(kappa_tmp*gamma[2] - gamma_tmp[1]*gamma[0] + gamma_tmp[0]*gamma[1]);

    	  kappa_tmp_minus = *kappa;
    	  gamma_tmp_minus[0] = gamma[0];
    	  gamma_tmp_minus[1] = gamma[1];
    	  gamma_tmp_minus[2] = gamma[2];

    	  assert(kappa_tmp_plus==kappa_tmp_plus && gamma_tmp_minus[0]==gamma_tmp_minus[0] &&gamma_tmp_minus[1]==gamma_tmp_minus[1] && gamma_tmp_minus[2]==gamma_tmp_minus[2]);

    	  *kappa = kappa_tmp_plus;
    	  gamma[0] = gamma_tmp_plus[0];
    	  gamma[1] = gamma_tmp_plus[1];
    	  gamma[2] = gamma_tmp_plus[2];

      }
    }

    // Convert units back to angles.
    x_source[0] /= Dl[NPlanes];
    x_source[1] /= Dl[NPlanes];

    *kappa = 1 - *kappa;

    if(x_source[0] != x_source[0] ||
       x_source[1] != x_source[1]){
      ERROR_MESSAGE();
      std::cout << x_source[0] << "  " << x_source[1] << "  " << std::endl;
      std::cout << gamma[0] << "  " << gamma[1] << "  " << gamma[2] << "  " <<
    		  *kappa << "  "  << kappa_off << std::endl;
      //	assert(0);
      exit(1);
    }


}
