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
 * charge = 4*pi*G*mass_scale/c^2 in units of Mpc // Fabien : we should remove this ?
 * charge = 4*pi*G in units of physical Mpc
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
  int tid;
  int start;
  int size;
  int NPlanes;
  bool flag_switch_deflection_off;
  bool flag_switch_lensing_off;
  PosType charge;
  LensPlane** lensing_planes;
  PosType* plane_redshifts;
  PosType* Dl;
  PosType* dDl;
};

/** \brief This function calculates the deflection, shear, convergence, rotation
 and time-delay of rays in parallel.
*/
void Lens::rayshooterInternal(
		unsigned long Npoints   /// number of points to be shot
		,Point *i_points        /// point on the image plane
		){
    
  // To force the computation of convergence, shear... -----
  // -------------------------------------------------------
            
  int NLastPlane;
  PosType tmpDs,tmpdDs,tmpZs;

  // If there are no points to shoot, then we quit.
  if(Npoints == 0) return;

  // If a lower redshift source (compared to the farthest lens plane) is being used
  if(toggle_source_plane)
  {
    NLastPlane = index_of_new_sourceplane;

    assert(NLastPlane <= lensing_planes.size());
    tmpDs = Dl[index_of_new_sourceplane];
    tmpdDs = dDl[index_of_new_sourceplane];
    tmpZs = plane_redshifts[index_of_new_sourceplane];

    Dl[index_of_new_sourceplane] = Ds_implant;
    dDl[index_of_new_sourceplane] = dDs_implant;
    plane_redshifts[index_of_new_sourceplane] = zs_implant;
  }
  else{ NLastPlane = lensing_planes.size(); }


  // For refining the grid and shoot new rays.
  int nthreads, rc;
  nthreads = N_THREADS;

  int chunk_size;
  do{
    chunk_size = (int)Npoints/nthreads;
    if(chunk_size == 0) nthreads /= 2;
  }while(chunk_size == 0);

  pthread_t threads[nthreads];
  TmpParams *thread_params = new TmpParams[nthreads];
            
  // This is for multi-threading :
  for(int i=0; i<nthreads;i++)
  {
    thread_params[i].i_points = i_points;
    thread_params[i].size = chunk_size;
    if(i == nthreads-1)
      thread_params[i].size = (int)Npoints - (nthreads-1)*chunk_size;
    thread_params[i].start = i*chunk_size;
    thread_params[i].tid = i;
    thread_params[i].flag_switch_deflection_off = flag_switch_deflection_off;
    thread_params[i].flag_switch_lensing_off = flag_switch_lensing_off;
    thread_params[i].charge = charge;
    thread_params[i].lensing_planes = &lensing_planes[0];
    thread_params[i].plane_redshifts = &plane_redshifts[0];
    thread_params[i].Dl = &Dl[0];
    thread_params[i].dDl = &dDl[0];
    thread_params[i].NPlanes = NLastPlane;
    rc = pthread_create(&threads[i], NULL, compute_rays_parallel, (void*) &thread_params[i]);
    assert(rc==0);
  }
    
  for(int i = 0; i < nthreads; i++)
  {
    rc = pthread_join(threads[i], NULL);
    assert(rc==0);
  }

  delete[] thread_params;

  if(toggle_source_plane)
  {
    // The initial values for the plane are reset here
    Dl[index_of_new_sourceplane] = tmpDs;
    dDl[index_of_new_sourceplane] = tmpdDs;
    plane_redshifts[index_of_new_sourceplane] = tmpZs;
  }

}



void *compute_rays_parallel(void *_p)
{
  TmpParams *p = (TmpParams *) _p;
  int chunk_size = p->size;
  int start      = p->start;
  int end        = start + chunk_size;
    
  int i, j;
  
  PosType xx[2],fac;
  PosType aa,bb,cc;
  PosType alpha[2];
    
  KappaType kappa,gamma[3];
    
  PosType xminus[2],xplus[2];
  PosType kappa_minus,gamma_minus[3],kappa_plus,gamma_plus[3];
    
  KappaType phi;
    
// Main loop : loop over the points of the image
for(i = start; i < end; i++)
  {
    
    // In case e.g. a temporary point is outside of the grid.
    if(p->i_points[i].in_image == MAYBE) continue;
      
    // find position on first lens plane in comoving units
    p->i_points[i].image->x[0] = p->i_points[i].x[0] * p->Dl[0];
    p->i_points[i].image->x[1] = p->i_points[i].x[1] * p->Dl[0];
      
    xminus[0] = 0;
    xminus[1] = 0;
    
    // Set magnification matrix on first plane.
    kappa_minus = 0;
    gamma_minus[0] = 0;
    gamma_minus[1] = 0;
    gamma_minus[2] = 0;
      
    // Setting phi on the first plane.
    phi = 0.0;
      
    // Default values :
    p->i_points[i].kappa = 1;  // This is actually 1-kappa until after the loop through the planes.
    p->i_points[i].gamma[0] = 0;
    p->i_points[i].gamma[1] = 0;
    p->i_points[i].gamma[2] = 0;
    p->i_points[i].dt = 0;
    
    // In case we don't want to compute the values :
    if(p->flag_switch_lensing_off)
    {
      p->i_points[i].image->x[0] /= p->Dl[0];
      p->i_points[i].image->x[1] /= p->Dl[0];
      p->i_points[i].kappa = p->i_points[i].image->kappa = 0.0;
      p->i_points[i].invmag = 1.0;
      p->i_points[i].dt = 0.0;
      
      continue;
    }
      
      
    /* ************************************************************************************
    We compute deflection, shear, convergence, rotation
    and time-delay of rays in parallel.
    ************************************************************************************ */
      
    // Time delay at first plane : position on the observer plane is (0,0) => no need to take difference of positions.
    p->i_points[i].dt = 0.5*( p->i_points[i].image->x[0]*p->i_points[i].image->x[0] + p->i_points[i].image->x[1]*p->i_points[i].image->x[1] )/ p->dDl[0] ;

      
    //  std::cout << "x1 = " << p->i_points[i].image->x[0] << "  ;  x2 = " << p->i_points[i].image->x[1] << std::endl ;
    //  std::cout << "d_A(lense) = " << p->Dl[0] << std::endl ;
      
    // Begining of the loop through the planes :
    // Each iteration leaves i_point[i].image on plane (j+1)
    for(j = 0; j < p->NPlanes ; ++j)
    {

      // convert to physical coordinates on the plane j
      xx[0] = p->i_points[i].image->x[0]/(1+p->plane_redshifts[j]);
      xx[1] = p->i_points[i].image->x[1]/(1+p->plane_redshifts[j]);
      
      assert(xx[0] == xx[0] && xx[1] == xx[1]);

      // std::cout << "xx[0] = " << xx[0] << "  ;  xx[1] = " << xx[1] << std::endl ;
        
      p->lensing_planes[j]->force(alpha,&kappa,gamma,&phi,xx); // Computed in physical coordinates.

        assert(alpha[0] == alpha[0] && alpha[1] == alpha[1]);
        assert(gamma[0] == gamma[0] && gamma[1] == gamma[1]);
        assert(kappa == kappa);
        assert(!isinf(kappa));

        fac = 1/(1+p->plane_redshifts[j]);
    	  /* multiply by fac to obtain 1/comoving_distance/physical_distance
    	   * such that a multiplication with the charge (in units of physical distance)
    	   * will result in a 1/comoving_distance quantity */ // 1 / comoving_distance squared ?
    	  kappa *= fac;
    	  gamma[0] *= fac;
    	  gamma[1] *= fac;
    	  gamma[2] *= fac;
          // dt *= fac ;
	
        assert(gamma[0] == gamma[0] && gamma[1] == gamma[1]);
        assert(kappa == kappa);
        assert(phi == phi);
        
      if(p->flag_switch_deflection_off){ alpha[0] = alpha[1] = 0.0; }
      
        
      // This computes \vec{x}^{j+1} in terms of \vec{x}^{j} , \vec{x}^{j-1} and \vec{\alpha}^{j}
      // according to Eq. (19) of paper GLAMER II -----------------------------------------------

      aa = (p->dDl[j+1] + p->dDl[j])/p->dDl[j];
      bb = p->dDl[j+1]/p->dDl[j];
      cc = p->charge * p->dDl[j+1];
        
      xplus[0] = aa*p->i_points[i].image->x[0] - bb*xminus[0] - cc*alpha[0];
      xplus[1] = aa*p->i_points[i].image->x[1] - bb*xminus[1] - cc*alpha[1];
      
      xminus[0] = p->i_points[i].image->x[0];
      xminus[1] = p->i_points[i].image->x[1];
      

      // Change in the value of the position.
      p->i_points[i].image->x[0] = xplus[0];
      p->i_points[i].image->x[1] = xplus[1];

        
      // ----------------------------------------------------------------------------------------

        
        // This computes (\kappa^{j+1}, \gamma_1^{j+1}, \gamma_2^{j+1}, \gamma_3^{j+1})
        // in terms of the j-plane quantities and according to Eq. (22) of GLAMER II.
        
            // Here the coefficients aa, bb and cc are used for a completely different calculation,
            // they are not the same as they were defined above. ----------------------------------
            aa = (p->dDl[j+1] + p->dDl[j]) * p->Dl[j] / p->dDl[j] / p->Dl[j+1];
            if(j>0)
            {
                bb = p->dDl[j+1] * p->Dl[j-1] / p->dDl[j] / p->Dl[j+1];
            }
            else bb = 0;
            cc = p->charge * p->dDl[j+1] * p->Dl[j] / p->Dl[j+1];
            // ------------------------------------------------------------------------------------
            
            
        // Computation of the "plus quantities", i.e. the  next plane quantities --------------------
        kappa_plus = aa*p->i_points[i].kappa - bb*kappa_minus
    			  - cc*(kappa*p->i_points[i].kappa + gamma[0]*p->i_points[i].gamma[0] + gamma[1]*p->i_points[i].gamma[1]);
            
        gamma_plus[0] = aa*p->i_points[i].gamma[0] - bb*gamma_minus[0]
    	          - cc*(gamma[0]*p->i_points[i].kappa + kappa*p->i_points[i].gamma[0] - gamma[1]*p->i_points[i].gamma[2]);
	
        gamma_plus[1] = aa*p->i_points[i].gamma[1] - bb*gamma_minus[1]
    	          - cc*(gamma[1]*p->i_points[i].kappa + kappa*p->i_points[i].gamma[1] + gamma[0]*p->i_points[i].gamma[2]);
	
        gamma_plus[2] = aa*p->i_points[i].gamma[2] - bb*gamma_minus[2]
    	          - cc*(kappa*p->i_points[i].gamma[2] - gamma[1]*p->i_points[i].gamma[0] + gamma[0]*p->i_points[i].gamma[1]);
        // ------------------------------------------------------------------------------------------
            
            
        // Assigning them to the "minus quantities" for next plane occurence of the loop ------------
        kappa_minus = p->i_points[i].kappa;
        gamma_minus[0] = p->i_points[i].gamma[0];
        gamma_minus[1] = p->i_points[i].gamma[1];
        gamma_minus[2] = p->i_points[i].gamma[2];
        // ------------------------------------------------------------------------------------------
            
        assert(kappa_plus==kappa_plus && gamma_minus[0]==gamma_minus[0] && gamma_minus[1]==gamma_minus[1] && gamma_minus[2]==gamma_minus[2]);

            
        // Updating the point quantities ----------------
        p->i_points[i].kappa = kappa_plus;
        p->i_points[i].gamma[0] = gamma_plus[0];
        p->i_points[i].gamma[1] = gamma_plus[1];
        p->i_points[i].gamma[2] = gamma_plus[2];
        // ----------------------------------------------

            
        
        // Geometric time delay with added potential
            p->i_points[i].dt += 0.5*( (xplus[0] - xminus[0])*(xplus[0] - xminus[0]) + (xplus[1] - xminus[1])*(xplus[1] - xminus[1]) )/p->dDl[j+1] - (1 + p->plane_redshifts[j]) * phi * p->charge ;
      
        // Check that the 1+z factor must indeed be there (because the x positions have been rescaled, so it may be different compared to the draft).

    } // End of the loop going through the planes

      
      
    // Subtracting off a term that makes the unperturbed ray to have zero time delay
    p->i_points[i].dt -= 0.5*( p->i_points[i].image->x[0]*p->i_points[i].image->x[0] + p->i_points[i].image->x[1]*p->i_points[i].image->x[1] ) / p->Dl[p->NPlanes];

      
    // Convert units back to angles.
    // Fabien : be careful ! These angles are not the same as those computed after the comment 'find position on first lens plane in comoving units' above, namely the angles we start with in this function. Values are close but still different. The change occurs after the comment 'Change in the value of the position.' above and by the fact that below we divide by Dl[p->NPlanes] and not Dl[0].
    p->i_points[i].image->x[0] /= p->Dl[p->NPlanes];
    p->i_points[i].image->x[1] /= p->Dl[p->NPlanes];
      
      
    // We go from kappa denoting 1-kappa to kappa denoting kappa
    p->i_points[i].kappa = 1 - p->i_points[i].kappa;

      
    // Computation of the inverse magnitude --------------------------------------------------------
    p->i_points[i].invmag = (1-p->i_points[i].kappa)*(1-p->i_points[i].kappa)
                                            - p->i_points[i].gamma[0]*p->i_points[i].gamma[0]
                                            - p->i_points[i].gamma[1]*p->i_points[i].gamma[1]
                                            + p->i_points[i].gamma[2]*p->i_points[i].gamma[2];
    // ---------------------------------------------------------------------------------------------
      
    
    // Putting the final values of the quantities in the real image point -----
    p->i_points[i].image->invmag = p->i_points[i].invmag;
    p->i_points[i].image->kappa = p->i_points[i].kappa;
    p->i_points[i].image->gamma[0] = p->i_points[i].gamma[0];
    p->i_points[i].image->gamma[1] = p->i_points[i].gamma[1];
    p->i_points[i].image->gamma[2] = p->i_points[i].gamma[2];
    p->i_points[i].image->dt = p->i_points[i].dt;
    // ------------------------------------------------------------------------
      
      

      
/*
// TODO: check
    if(p->i_points[i].image->x[0] != p->i_points[i].image->x[0] ||
       p->i_points[i].image->x[1] != p->i_points[i].image->x[1] ||
       p->i_points[i].invmag != p->i_points[i].invmag)
    {
      ERROR_MESSAGE();
      std::cout << p->i_points[i].image->x[0] << "  " << p->i_points[i].image->x[1] << "  " << p->i_points[i].invmag << std::endl;
      std::cout << p->i_points[i].gamma[0] << "  " << p->i_points[i].gamma[1] << "  " << p->i_points[i].gamma[2] << "  " <<
    		  p->i_points[i].kappa << std::endl;
      //	assert(0);
      exit(1);
    }
*/
      
} // End of the main loop.
  
  return 0;
  
}
