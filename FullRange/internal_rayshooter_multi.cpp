/*
 * internal_rayshooter_multi.cpp
 *
 *  Created on: Jan 19, 2012
 *      Author: mpetkova, bmetcalf
 */

#include "slsimlib.h"
#include <powell.h>

/** 
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
 * charge = 4*pi*G*mass_scale/c^2 in units of Mpc // we should remove this ?
 * charge = 4*pi*G in units of physical Mpc
 *
 * i_points[].x[] is in angular units.
 *
 * Warning: Is not valid for a non-flat universe.
 */

#ifndef N_THREADS
#define N_THREADS 1
#endif


//void *compute_rays_parallel(void *_p);


/** \brief This function calculates the deflection, shear, convergence, rotation
 and time-delay of rays in parallel.
 */
void Lens::rayshooter(RAY &ray){
  
  LinkedPoint point;
  point.x[0] = ray.x[0];
  point.x[1] = ray.x[1];
  
  rayshooterInternal(1,&point);
  ray = point;
};

void Lens::rayshooterInternal(
                                unsigned long Npoints   /// number of points to be shot
                              , Point *i_points         /// point on the image plane
                              , bool verbose         /// verbose option
){
  
  // To force the computation of convergence, shear... -----
  // -------------------------------------------------------
  
  
  assert(plane_redshifts.size() > 0);
  if(plane_redshifts.size() == 1){  // case of no lens plane
    
    for(int ii = 0; ii < Npoints; ++ii){
    
      i_points[ii].image->x[0] = i_points[ii].x[0];
      i_points[ii].image->x[1] = i_points[ii].x[1];
      i_points[ii].dt = 0;
      i_points[ii].image->dt = i_points[ii].dt;
      
      i_points[ii].A.setToI();
      i_points[ii].image->A.setToI();
    }
    
    return;
  }
  
  
  //int NLastPlane;
  //PosType tmpDs,tmpdDs,tmpdTs,tmpZs;
  
  // If there are no points to shoot, then we quit.
  if(Npoints == 0) return;
  
  double source_z;
  
  // If a lower redshift source (compared to the farthest lens plane) is being used
  if(toggle_source_plane)
  {
    source_z = zs_implant;
  }
  else{ source_z = plane_redshifts.back(); }
  
  
  // For refining the grid and shoot new rays.
  int nthreads, rc;
  nthreads = Utilities::GetNThreads();
  
  int chunk_size;
  do{
    chunk_size = (int)Npoints/nthreads;
    if(chunk_size == 0) nthreads /= 2;
  }while(chunk_size == 0);
  
  std::thread thr[nthreads];
  
  // This is for multi-threading :
  for(int i=0; i<nthreads;i++)
  {
    
    int size = chunk_size;
    if(i == nthreads-1)
      size = (int)Npoints - (nthreads-1)*chunk_size;
    int start = i*chunk_size;
    
    thr[i] = std::thread(&Lens::compute_rays_parallel,this,start,size,i_points,&source_z,false,false);
  }
 
  for(int i = 0; i < nthreads; i++) thr[i].join();
}


// NEW VERSION OF RAYSHOOTER USING ONLY PLANE i TO COMPUTE THE RAY POSITIONS AND LENSING QUANTITIES ON PLANE i+1.

void Lens::compute_rays_parallel(int start
                                 ,int chunk_size
                                 ,Point *i_points
                                 ,double *source_zs
                                 ,bool multiZs
                                 ,bool verbose)
{
  //TmpParams *p = (TmpParams *) _p;
  //int chunk_size = p->size;
  //int start      = p->start;
  int end        = start + chunk_size;
  
  int i, j;
  
  PosType xx[2];
  PosType aa,bb;
  PosType alpha[2];
  
  KappaType kappa,gamma[3];
  KappaType phi;
  
  Matrix2x2<PosType> G;

  PosType SumPrevAlphas[2];
  Matrix2x2<PosType> SumPrevAG;
  
  PosType *theta;
  
  long jmax = lensing_planes.size();
  double Dls_Ds; // this is the ratio between of the distance between the last lens plane and the source to the distance to the source
  double D_Ds; // this is the ratio between of the distance to the last lens plane and the source to the distance to the source

  if(!multiZs){
    if(source_zs[0] == plane_redshifts.back() ){
      Dls_Ds = dDl.back() / Dl.back();
      D_Ds = Dl[Dl.size() - 2] / Dl.back();
    }else{
      PosType Dls,Ds;
      FindSourcePlane(source_zs[0],jmax,Dls,Ds);
      Dls_Ds = Dls / Ds;
      if(jmax > 0) D_Ds = Dl[jmax-1] / Ds;
    }
  }
  
  // Main loop : loop over the points of the image
  for(i = start; i < end; i++)
  {
    // In case e.g. a temporary point is outside of the grid.
    if(i_points[i].in_image == MAYBE) continue;
    
    theta = i_points[i].image->x;
    theta[0] = i_points[i].x[0];
    theta[1] = i_points[i].x[1];

    // Initializing SumPrevAlphas :
    SumPrevAlphas[0] = theta[0];
    SumPrevAlphas[1] = theta[1];

    // Initializing SumPrevAG :
    SumPrevAG.setToI();
    
    // Setting phi on the first plane.
    phi = 0.0;
    
    // Default values :
    i_points[i].A.setToI();
    i_points[i].dt = 0;
    
    // In case we don't want to compute the values :
    if(flag_switch_lensing_off)
    {
      i_points[i].image->A.setToI();
      i_points[i].dt = 0.0;
      
      continue;
    }
    
    // Time delay at first plane : position on the observer plane is (0,0) => no need to take difference of positions.
    i_points[i].dt = 0;
    
    //0.5*( p->i_points[i].image->x[0]*p->i_points[i].image->x[0] + p->i_points[i].image->x[1]*p->i_points[i].image->x[1] )/ p->dDl[0] ;
    
    if(multiZs){
      PosType Dls,Ds;
      FindSourcePlane(source_zs[i],jmax,Dls,Ds);
      Dls_Ds = Dls / Ds;
      if(jmax > 0) D_Ds = Dl[jmax-1]/Ds;
    }
    
    // Begining of the loop through the planes :
    // Each iteration leaves i_point[i].image on plane (j+1)

    for(j = 0; j < jmax ; ++j)
      {
      
      double Dphysical = Dl[j]/(1 + plane_redshifts[j]);
      // convert to physical coordinates on the plane j, just for force calculation
      xx[0] = theta[0] *  Dphysical;
      xx[1] = theta[1] *  Dphysical;
      // PhysMpc = ComMpc / (1+z)
      
      assert(xx[0] == xx[0] && xx[1] == xx[1]);
      
      ////////////////////////////////////////////////////////
      
      lensing_planes[j]->force(alpha,&kappa,gamma,&phi,xx);
      // Computed in physical coordinates, xx is in PhysMpc.
      
      ////////////////////////////////////////////////////////
      
      assert(alpha[0] == alpha[0] && alpha[1] == alpha[1]);
      assert(gamma[0] == gamma[0] && gamma[1] == gamma[1]);
      assert(kappa == kappa);
      if(std::isinf(kappa)) { std::cout << "xx = " << xx[0] << " " << xx[1] << std::endl ;}
      assert(!std::isinf(kappa));
      
      G[0] = kappa + gamma[0];    G[1] = gamma[1];
      G[2] = gamma[1]; G[3] = kappa - gamma[0];
  
      /* multiply by fac to obtain 1/comoving_distance/physical_distance
       * such that a multiplication with the charge (in units of physical distance)
       * will result in a 1/comoving_distance quantity */
      
      G *= charge * Dl[j] / (1 + plane_redshifts[j]);
        
      assert(gamma[0] == gamma[0] && gamma[1] == gamma[1]);
      assert(kappa == kappa);
      assert(phi == phi);
      
      if(flag_switch_deflection_off){ alpha[0] = alpha[1] = 0.0; }
      
      // This computes \vec{x}^{j+1} in terms of \vec{x}^{j}
      // according to the corrected Eq. (18) of paper GLAMER II ---------------------------------
      
      // Adding the j-plane alpha contribution to the sum \Sum_{k=1}^{j} \vec{alpha_j} :
      SumPrevAlphas[0] -= charge * alpha[0] ;
      SumPrevAlphas[1] -= charge * alpha[1] ;
      
      if(j < jmax-1 ){
        aa = dDl[j+1] / Dl[j+1];
        bb = Dl[j] / Dl[j+1];
      }else{
        aa = Dls_Ds;
        bb = D_Ds;
      }
      
      theta[0] = bb * theta[0] + aa * SumPrevAlphas[0];
      theta[1] = bb * theta[1] + aa * SumPrevAlphas[1];
      
      // ----------------------------------------------------------------------------------------
            
      // Sum_{k=1}^{j} Dl[k] A^k.G^k
      SumPrevAG -= (G * (i_points[i].A)) ;
      
      // Computation of the "plus quantities", i.e. the  next plane quantities :
      i_points[i].A = i_points[i].A * bb + SumPrevAG * aa;
      
      // ----------------------------------------------
      
      // Geometric time delay with added potential
      //p->i_points[i].dt += 0.5*( (xplus[0] - xminus[0])*(xplus[0] - xminus[0]) + (xplus[1] - xminus[1])*(xplus[1] - xminus[1]) ) * p->dTl[j+1] /p->dDl[j+1] /p->dDl[j+1] - phi * p->charge ; /// in Mpc  ???
      
      // Check that the 1+z factor must indeed be there (because the x positions have been rescaled, so it may be different compared to the draft).
      // Remark : Here the true lensing potential is not "phi" but "phi * p->charge = phi * 4 pi G".
      
      
    } // End of the loop going through the planes
    
    
    // Subtracting off a term that makes the unperturbed ray to have zero time delay
    //p->i_points[i].dt -= 0.5*( p->i_points[i].image->x[0]*p->i_points[i].image->x[0] + p->i_points[i].image->x[1]*p->i_points[i].image->x[1] ) / p->Dl[NLastPlane];
    
    // Conversion of dt from Mpc (physical Mpc) to Years -----------------
    i_points[i].dt *= MpcToSeconds * SecondToYears ;
    
    // ---------------------------------------------------------------------------------------------
    
    
    // Putting the final values of the quantities in the image point -----
    i_points[i].image->A = i_points[i].A;
    i_points[i].image->dt = i_points[i].dt;
    // ------------------------------------------------------------------------
    
    
    // TEST : showing final quantities
    // ------------------------------=
    if(verbose) std::cout << "RSI final : X X | " << Dl[jmax] << " | " << i_points[i].kappa() << " " << i_points[i].gamma1() << " " << i_points[i].gamma2() << " " << i_points[i].gamma3() << " " << i_points[i].invmag() << " | " << i_points[i].dt << std::endl ;
    
  } // End of the main loop.
  
  return 0;
}




// FORMER VERSION OF RAYSHOOTER USING PLANES i AND i-1 TO COMPUTE THE RAY POSITIONS AND LENSING QUANTITIES ON PLANE i+1.
/*
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
  KappaType phi;
  
  PosType xminus[2],xplus[2];
  PosType kappa_minus,gamma_minus[3],kappa_plus,gamma_plus[3];

  
  // Main loop : loop over the points of the image
  for(i = start; i < end; i++)
  {
    
    // In case e.g. a temporary point is outside of the grid.
    if(p->i_points[i].in_image == MAYBE) continue;
    
    // find position on first lens plane in comoving units
    p->i_points[i].image->x[0] = p->i_points[i].x[0] * p->Dl[0]; // x^1 = \theta * D_1
    p->i_points[i].image->x[1] = p->i_points[i].x[1] * p->Dl[0];
    
    xminus[0] = 0; // x^0 = 0
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
    
    
    // Time delay at first plane : position on the observer plane is (0,0) => no need to take difference of positions.
    p->i_points[i].dt = 0.5*( p->i_points[i].image->x[0]*p->i_points[i].image->x[0] + p->i_points[i].image->x[1]*p->i_points[i].image->x[1] )/ p->dDl[0] ;


    // TEST : showing initial quantities
    // ---------------------------------
    std::cout << "RSI initial : X X | X | " << p->i_points[i].image->x[0] << " " << p->i_points[i].image->x[1] << " | " << p->i_points[i].kappa << " " << p->i_points[i].gamma[0] << " " << p->i_points[i].gamma[1] << " " << p->i_points[i].gamma[2] << " X | " << p->i_points[i].dt << std::endl ;
    
    // Begining of the loop through the planes :
    // Each iteration leaves i_point[i].image on plane (j+1)
    for(j = 0; j < p->NPlanes ; ++j)
    {

      // For Test :
      // std::cout << "> p->plane_redshifts = " << p->plane_redshifts[j] << std::endl ;
      
      // convert to physical coordinates on the plane j, just for force calculation
      xx[0] = p->i_points[i].image->x[0]/(1+p->plane_redshifts[j]);
      xx[1] = p->i_points[i].image->x[1]/(1+p->plane_redshifts[j]);
      // PhysMpc = ComMpc / (1+z)
      
      assert(xx[0] == xx[0] && xx[1] == xx[1]);
      
      // std::cout << "p->i_points[i].image->x[0] = " << p->i_points[i].image->x[0] << " , p->i_points[i].image->x[1] = " << p->i_points[i].image->x[1] << std::endl ;
      // std::cout << "dDl[" << j << "] = " << p->dDl[j] << " , dDl[" << j+1 << "] = " << p->dDl[j+1] << std::endl;
      

      ////////////////////////////////////////////////////////
      
      p->lensing_planes[j]->force(alpha,&kappa,gamma,&phi,xx);
      // Computed in physical coordinates, xx is in PhysMpc.
      
      ////////////////////////////////////////////////////////
      
      
      assert(alpha[0] == alpha[0] && alpha[1] == alpha[1]);
      assert(gamma[0] == gamma[0] && gamma[1] == gamma[1]);
      assert(kappa == kappa);
      if(std::isinf(kappa)) { std::cout << "xx = " << xx[0] << " " << xx[1] << std::endl ;}
      assert(!std::isinf(kappa));
      
      fac = 1/(1+p->plane_redshifts[j]);
      // multiply by fac to obtain 1/comoving_distance/physical_distance
      // such that a multiplication with the charge (in units of physical distance)
      // will result in a 1/comoving_distance quantity
      // 1 / comoving_distance squared ?
      kappa *= fac;
      gamma[0] *= fac;
      gamma[1] *= fac;
      gamma[2] *= fac;
      
      assert(gamma[0] == gamma[0] && gamma[1] == gamma[1]);
      assert(kappa == kappa);
      assert(phi == phi);
      
      if(p->flag_switch_deflection_off){ alpha[0] = alpha[1] = 0.0; }
      
      
      // This computes \vec{x}^{j+1} in terms of \vec{x}^{j} , \vec{x}^{j-1} and \vec{\alpha}^{j}
      // according to Eq. (19) of paper GLAMER II -----------------------------------------------
      
      aa = (p->dDl[j+1] + p->dDl[j])/p->dDl[j]; // (Dls*(1+zs) + Dl*(1+zl)) / Dl*(1+zl) = 1 + (Dls/Dl)*((1+zs)/(1+zl)) ;
      bb = p->dDl[j+1]/p->dDl[j];               // (Dls/Dl)*((1+zs)/(1+zl)) ;
      cc = p->charge * p->dDl[j+1];             // 4 pi G * Dls * (1+zs) ;
           
      assert(aa == aa);
      assert(bb == bb);
      assert(cc == cc);
      // std::cout << "RayshooterInternal : aa = " << aa << " , bb = " << bb << " , cc = " << cc << std::endl ;
      
      xplus[0] = aa*p->i_points[i].image->x[0] - bb*xminus[0] - cc*alpha[0];
      xplus[1] = aa*p->i_points[i].image->x[1] - bb*xminus[1] - cc*alpha[1];
      // x is in ComMpc, cc*alpha in [(PhysMpc/mass)*ComMpc]*(mass/PhysMpc) = ComMpc.
        
      // For the test with a source at the position (0,0) centered on the lens :
      // std::cout << "alpha we have at the end of force_halo : " << alpha[0] << " " << alpha[1] << std::endl ;
      // std::cout << "alpha we should have at the end of force_halo : " << aa*p->i_points[i].image->x[0] / cc << " " << aa*p->i_points[i].image->x[1] / cc << std::endl << std::endl ;

      // x^j becomes x^{j-1}
      xminus[0] = p->i_points[i].image->x[0];
      xminus[1] = p->i_points[i].image->x[1];
      
      // x^{j+1} becomes x^j
      p->i_points[i].image->x[0] = xplus[0];
      p->i_points[i].image->x[1] = xplus[1];
      
      // std::cout << "In rayshooter : !!! " << p->i_points[i].image->x[0] / p->Dl[p->NPlanes] << " " << p->i_points[i].image->x[1] / p->Dl[p->NPlanes] << " !!!" << std::endl ; // This step is done after !
      // ----------------------------------------------------------------------------------------
      
      
      // This computes (\kappa^{j+1}, \gamma_1^{j+1}, \gamma_2^{j+1}, \gamma_3^{j+1})
      // in terms of the j-plane quantities and according to Eq. (22) of GLAMER II ----------------
      
      // Here the coefficients aa, bb and cc are used for a completely different calculation,
      // they are not the same as they were defined above.
      
      aa = (p->dDl[j+1] + p->dDl[j]) * p->Dl[j] / p->dDl[j] / p->Dl[j+1];
      if(j>0)
      {
        bb = p->dDl[j+1] * p->Dl[j-1] / p->dDl[j] / p->Dl[j+1];
      }
      else bb = 0;
      cc = p->charge * p->dDl[j+1] * p->Dl[j] / p->Dl[j+1];
      

      // Computation of the "plus quantities", i.e. the  next plane quantities :
      kappa_plus = aa*p->i_points[i].kappa - bb*kappa_minus
      - cc*(kappa*p->i_points[i].kappa + gamma[0]*p->i_points[i].gamma[0] + gamma[1]*p->i_points[i].gamma[1]);
      
      gamma_plus[0] = aa*p->i_points[i].gamma[0] - bb*gamma_minus[0]
      - cc*(gamma[0]*p->i_points[i].kappa + kappa*p->i_points[i].gamma[0] - gamma[1]*p->i_points[i].gamma[2]);
      
      gamma_plus[1] = aa*p->i_points[i].gamma[1] - bb*gamma_minus[1]
      - cc*(gamma[1]*p->i_points[i].kappa + kappa*p->i_points[i].gamma[1] + gamma[0]*p->i_points[i].gamma[2]);
      
      gamma_plus[2] = aa*p->i_points[i].gamma[2] - bb*gamma_minus[2]
      - cc*(kappa*p->i_points[i].gamma[2] - gamma[1]*p->i_points[i].gamma[0] + gamma[0]*p->i_points[i].gamma[1]);
      
      // ------------------------------------------------------------------------------------------

      // std::cout << std::endl ;
      // std::cout << "RayshooterInternal : kappa_minus = " << kappa_minus << " , kappa = " << kappa << " , gamma[0] = " << gamma[0] << " , gamma[1] = " << gamma[1] << std::endl ;
      // std::cout << "RayshooterInternal : aa = " << aa << " , bb = " << bb << " , cc = " << cc << std::endl ;
      // std::cout << "RayshooterInternal : p->points[i] : kappa = " << p->i_points[i].kappa << " , gamma[0] = " << p->i_points[i].gamma[0] << " , gamma[1] = " << p->i_points[i].gamma[1] << std::endl ;

      // std::cout << "RayshooterInternal : plane " << j << " , z = " << p->plane_redshifts[j] << " , 1-kappa_plus = " << 1-kappa_plus << std::endl ;
        
      // Assigning them to the "minus quantities" for next plane occurence of the loop ------------
      kappa_minus = p->i_points[i].kappa;
      gamma_minus[0] = p->i_points[i].gamma[0];
      gamma_minus[1] = p->i_points[i].gamma[1];
      gamma_minus[2] = p->i_points[i].gamma[2];
      // ------------------------------------------------------------------------------------------
      
      
      if(!(kappa_plus==kappa_plus && gamma_minus[0]==gamma_minus[0] && gamma_minus[1]==gamma_minus[1] && gamma_minus[2]==gamma_minus[2])){
        std::cout << "p->dDl[j-1]" << "\t" << "p->dDl[j]" << "\t" << "p->dDl[j+1]" << std::endl;
        std::cout << p->dDl[j-1] << "\t" << p->dDl[j] << "\t" << p->dDl[j+1] << std::endl;
        std::cout << "aa" << "\t" << "bb" << "\t" << "cc" << std::endl;
        std::cout << aa << "\t" << bb << "\t" << cc << std::endl;
        
        
        std::cout << "kappa" << "\t" << "gamma[0]" << "\t" << "gamma[1]" << "\t" << "gamma[2]" << std::endl ;
        std::cout << kappa << "\t" << gamma[0] << "\t" << gamma[1] << "\t" << gamma[2] << std::endl ;
        
        
        std::cout << "p->i_points[i].kappa" << "\t" << "p->i_points[i].gamma[0]" << "\t" << "p->i_points[i].gamma[1]" << "\t" << "p->i_points[i].gamma[2]" << std::endl ;
        std::cout << p->i_points[i].kappa << "\t" << p->i_points[i].gamma[0] << "\t" << p->i_points[i].gamma[1] << "\t" << p->i_points[i].gamma[2] << std::endl ;
        
        std::cout << "x" << "\t" << "y" << std::endl;
        std::cout << p->i_points[i].image->x[0] << "\t" << p->i_points[i].image->x[1] << std::endl;
        std::cout << "kappa_minus" << "\t" << "gamma_minus[0]" << "\t" << "gamma_minus[1]" << "\t" << "gamma_minus[2]" << "\t" << std::endl ;
        std::cout << kappa_minus << "\t" << gamma_minus[0] << "\t" << gamma_minus[1] << "\t" << gamma_minus[2] << "\t" << std::endl ;
        std::cout << "kappa_plus" << "\t" << "gamma_plus[0]" << "\t" << "gamma_plus[1]" << "\t" << "gamma_plus[2]" << "\t" << std::endl ;
        std::cout << kappa_plus << "\t" << gamma_plus[0] << "\t" << gamma_plus[1] << "\t" << gamma_plus[2] << "\t" << std::endl ;
      }
      
      assert(kappa_plus==kappa_plus && gamma_minus[0]==gamma_minus[0] && gamma_minus[1]==gamma_minus[1] && gamma_minus[2]==gamma_minus[2]);
      
      
      // Updating the point quantities ----------------
      p->i_points[i].kappa = kappa_plus;
      p->i_points[i].gamma[0] = gamma_plus[0];
      p->i_points[i].gamma[1] = gamma_plus[1];
      p->i_points[i].gamma[2] = gamma_plus[2];
      // ----------------------------------------------
      
      
      // Geometric time delay with added potential
      p->i_points[i].dt += 0.5*( (xplus[0] - xminus[0])*(xplus[0] - xminus[0]) + (xplus[1] - xminus[1])*(xplus[1] - xminus[1]) )/p->dDl[j+1] - (1 + p->plane_redshifts[j]) * phi * p->charge ; /// in Mpc

      
      // Check that the 1+z factor must indeed be there (because the x positions have been rescaled, so it may be different compared to the draft).
      // Remark : Here the true lensing potential is not "phi" but "phi * p->charge = phi * 4 pi G".
      
      
      // TEST : showing plus quantities
      // ------------------------------
      std::cout << "RSI plane " << j << " : " << p->i_points[i].image->x[0] << " " << p->i_points[i].image->x[1] << " | " << p->Dl[j] << " | " << p->i_points[i].kappa << " " << p->i_points[i].gamma[0] << " " << p->i_points[i].gamma[1] << " " << p->i_points[i].gamma[2] << " X | " << p->i_points[i].dt << std::endl ;
      
      
    } // End of the loop going through the planes
    
    
    // Subtracting off a term that makes the unperturbed ray to have zero time delay
    p->i_points[i].dt -= 0.5*( p->i_points[i].image->x[0]*p->i_points[i].image->x[0] + p->i_points[i].image->x[1]*p->i_points[i].image->x[1] ) / p->Dl[p->NPlanes];
    
    
    // Convert units back to angles.
    p->i_points[i].image->x[0] /= p->Dl[p->NPlanes];
    p->i_points[i].image->x[1] /= p->Dl[p->NPlanes];
    
    // std::cout << "   Ds = " << p->Dl[p->NPlanes] << std::endl;
    // std::cout << "   dx = " << p->i_points[i].image->x[0]-p->i_points[i].x[0] << "  " << p->i_points[i].image->x[1]-p->i_points[i].x[1] << std::endl;
    
    // We go from kappa denoting 1-kappa to kappa denoting kappa
    p->i_points[i].kappa = 1 - p->i_points[i].kappa;
    
    // std::cout << "RayshooterInternal : kappa final = " << p->i_points[i].kappa << std::endl ;
    // std::cout << std::endl ;
      
    // Computation of the inverse magnitude --------------------------------------------------------
    p->i_points[i].invmag = (1-p->i_points[i].kappa)*(1-p->i_points[i].kappa)
    - p->i_points[i].gamma[0]*p->i_points[i].gamma[0]
    - p->i_points[i].gamma[1]*p->i_points[i].gamma[1]
    + p->i_points[i].gamma[2]*p->i_points[i].gamma[2];
    // ---------------------------------------------------------------------------------------------
    
    
    // Conversion of dt from Mpc (physical Mpc) to Years (also possible into days) -----------------
    
    p->i_points[i].dt *= MpcToSeconds * SecondToYears ;
    
    // ---------------------------------------------------------------------------------------------
    
    
    
    // Putting the final values of the quantities in the real image point -----
    p->i_points[i].image->invmag = p->i_points[i].invmag;
    p->i_points[i].image->kappa = p->i_points[i].kappa;
    p->i_points[i].image->gamma[0] = p->i_points[i].gamma[0];
    p->i_points[i].image->gamma[1] = p->i_points[i].gamma[1];
    p->i_points[i].image->gamma[2] = p->i_points[i].gamma[2];
    p->i_points[i].image->dt = p->i_points[i].dt;
    // ------------------------------------------------------------------------
    
    
    // TEST : showing final quantities
    // ------------------------------=
    std::cout << "RSI final : X X | " << p->Dl[p->NPlanes] << " | " << p->i_points[i].kappa << " " << p->i_points[i].gamma[0] << " " << p->i_points[i].gamma[1] << " " << p->i_points[i].gamma[2] << " " << p->i_points[i].invmag << " | " << p->i_points[i].dt << std::endl ;
    
 
     // TODO: check
     // if(p->i_points[i].image->x[0] != p->i_points[i].image->x[0] ||
     // p->i_points[i].image->x[1] != p->i_points[i].image->x[1] ||
     // p->i_points[i].invmag != p->i_points[i].invmag)
     // {
     // ERROR_MESSAGE();
     // std::cout << p->i_points[i].image->x[0] << "  " << p->i_points[i].image->x[1] << "  " << p->i_points[i].invmag << std::endl;
     // std::cout << p->i_points[i].gamma[0] << "  " << p->i_points[i].gamma[1] << "  " << p->i_points[i].gamma[2] << "  " <<
     // p->i_points[i].kappa << std::endl;
     //	assert(0);
     // exit(1);
     // }
    
  } // End of the main loop.
  
  return 0;
}
*/





/** \brief Collects information about the halos and kappa contributions along the light path
 
 Information on the nearest halos is collected where nearest is defined by rmax and mode.
 When mode == 2 the unlensed angular coordinates are used to evaluate proximity not the lensed ones.
 
 */
void Lens::info_rayshooter(
                           Point *i_point     /// point to be shot, must have image point linked
                           ,std::vector<Point_2d> & ang_positions  /// angular positions on each plane
                           ,std::vector<KappaType> & kappa_on_planes          /// convergence on each plane
                           ,std::vector<std::vector<LensHalo*>> & halo_neighbors  /// neighboring halos within rmax of ray on each plane
                           ,LensHalo **halo_max
                           ,KappaType &kappa_max
                           ,KappaType gamma_max[]
                           ,PosType rmax  /// distance from ray on each plane, units depend on mode parameter
                           ,short mode  /// 0:physical distance (Mpc), 1: comoving distance (Mpc), 2: angular distance (rad)
                           ,bool verbose
)
{
  
  // !!! would like to find the maximum contributing halo so that its contribution can be subtracted from the total
  
  std::cout << " Warning : Lens::info_rayshooter() has not been fully tested. " << std::endl;
  
  double source_z;
  
  // If a lower redshift source (compared to the farthest lens plane) is being used
  if(toggle_source_plane)
  {
    source_z = zs_implant;
  }
  else{ source_z = plane_redshifts.back(); }

  long jmax = lensing_planes.size();
  double Dls_Ds; // this is the ratio between of the distance between the last lens plane and the source to the distance to the source
  double D_Ds; // this is the ratio between of the distance to the last lens plane and the source to the distance to the source

  if(source_z == plane_redshifts.back() ){
    Dls_Ds = dDl.back() / Dl.back();
    D_Ds = Dl[Dl.size() - 2] / Dl.back();
  }else{
    PosType Dls,Ds;
    FindSourcePlane(source_z,jmax,Dls,Ds);
    Dls_Ds = Dls / Ds;
    if(jmax > 0) D_Ds = Dl[jmax-1] / Ds;
  }

  
  // If there are no points to shoot, then we quit.
  
  ang_positions.resize(jmax);
  halo_neighbors.resize(jmax);
  kappa_on_planes.resize(jmax);
  
  int j;
  
  PosType tmp_r,x_tmp[2];
  KappaType phi_tmp,kappa_tmp,gamma_tmp[3];
  
  PosType xx[2];
  PosType aa,bb;
  PosType alpha[2],alpha_tmp[2];
  
  KappaType kappa,gamma[3];
  KappaType phi;
  
  Matrix2x2<PosType> G;

  PosType SumPrevAlphas[2];
  Matrix2x2<PosType> SumPrevAG;
  
  PosType *theta;
 
  kappa_max = -1.0;
  
  // find angular position on first lens plane
  i_point->image->x[0] = i_point->x[0] ;
  i_point->image->x[1] = i_point->x[1] ;
  
  theta = i_point->image->x;
  theta[0] = i_point->x[0];
  theta[1] = i_point->x[1];

  // Initializing SumPrevAlphas :
  SumPrevAlphas[0] = theta[0];
  SumPrevAlphas[1] = theta[1];

  // Initializing SumPrevAG :
  SumPrevAG.setToI();
  
  // Setting phi on the first plane.
  phi = 0.0;
  
  // Default values :
  i_point->A.setToI();
  i_point->dt = 0;
  
  // Begining of the loop through the planes :
  // Each iteration leaves i_point[i].image on plane (j+1)
  
  for(j = 0; j < jmax ; ++j)
  {
    
    // convert to physical coordinates on the plane j
    double Dphysical = Dl[j]/(1 + plane_redshifts[j]);
    // convert to physical coordinates on the plane j, just for force calculation
    xx[0] = theta[0] *  Dphysical;
    xx[1] = theta[1] *  Dphysical;
    
    lensing_planes[j]->force(alpha,&kappa_on_planes[j],gamma,&phi,xx); // Computed in physical coordinates.
    if(flag_switch_deflection_off){ alpha[0] = alpha[1] = 0.0; }
    
    G[0] = kappa_on_planes[j] + gamma[0];    G[1] = gamma[1];
    G[2] = gamma[1]; G[3] = kappa_on_planes[j] - gamma[0];
    
    G *= charge * Dl[j] / (1 + plane_redshifts[j]);
 
    PosType SigmaCrit = cosmo.SigmaCrit(plane_redshifts[j]
                                        ,source_z);
    
    //kappa_on_planes[j] *= 1/(1+plane_redshifts[j]);
    kappa_on_planes[j] *= charge / SigmaCrit;
    
    tmp_r = rmax;
    if(mode == 1) tmp_r /= (1+plane_redshifts[j]);
    if(mode == 2) tmp_r *= Dl[j];
    
    lensing_planes[j]->getNeighborHalos(xx,tmp_r,halo_neighbors[j]);
    
    ang_positions[j][0] = i_point->image->x[0];
    ang_positions[j][1] = i_point->image->x[1];
 
   
    // Find the halo with the largest kappa
    for(int ii=0;ii<halo_neighbors[j].size();++ii){
      alpha_tmp[0] = alpha_tmp[1] = 0.0;
      kappa_tmp = 0.0;
      gamma_tmp[0] = gamma_tmp[1] = gamma_tmp[2] = 0.0;
      phi_tmp = 0.0;
      
      // Getting the halo position (in physical Mpc) :
      halo_neighbors[j][ii]->getX(x_tmp);
      
      // Taking the shift into account :
      x_tmp[0] = xx[0] - x_tmp[0];
      x_tmp[1] = xx[1] - x_tmp[1];
      
      halo_neighbors[j][ii]->force_halo(alpha_tmp,&kappa_tmp,gamma_tmp,&phi_tmp,x_tmp,false);
      kappa_tmp /= SigmaCrit;
      if(kappa_tmp > kappa_max){
        kappa_max = kappa_tmp;
        *halo_max = halo_neighbors[j][ii];
        gamma_max[0] = gamma_tmp[0]/SigmaCrit;
        gamma_max[1] = gamma_tmp[1]/SigmaCrit;
      }
    }
  
    SumPrevAlphas[0] -= charge * alpha[0] ;
    SumPrevAlphas[1] -= charge * alpha[1] ;
    
    if(j < jmax-1 ){
      aa = dDl[j+1] / Dl[j+1];
      bb = Dl[j] / Dl[j+1];
    }else{
      aa = Dls_Ds;
      bb = D_Ds;
    }
    
    theta[0] = bb * theta[0] + aa * SumPrevAlphas[0];
    theta[1] = bb * theta[1] + aa * SumPrevAlphas[1];

    // Sum_{k=1}^{j} Dl[k] A^k.G^k
    SumPrevAG -= (G * (i_point->A)) ;
    
    // Computation of the "plus quantities", i.e. the  next plane quantities :
    i_point->A = i_point->A * bb + SumPrevAG * aa;
    
  } // End of the loop going through the planes
  
  return;
  
}

/**  \brief Find an image position for a source position.
 
 This routine finds the image position by minimizing the seporation on the source plane with Powell's method of minimization.  This will not find all images.  For that you must use another routine.  In the weak lensing regiam this should be sufficient.
 */
/*
RAY Lens::find_image(
          Point_2d y_source     /// input position of source (radians)
          ,Point_2d &x_image    /// initialized with guess for image postion (radians)
          ,PosType z_source     /// redshift of source
          ,PosType ytol2        /// target tolerance in source position squared
          ,PosType &dy2        /// final value of Delta y ^2
          ,int sign             /// sign of magnification
){
  
  PosType tmp_zs = getSourceZ();
  if(tmp_zs != zsource) ResetSourcePlane(z_source,false);
  
  LinkedPoint p;
  if(sign==0){
    // get sign of magnification at inital point
    p[0]=x_image[0];
    p[1]=x_image[1];
    rayshooterInternal(1,&p);
    sign = sgn(p.invmag());
  }
  
  MINyFunction minfunc(*this,y_source,sign);
  int iter;
  
  powell_tp(x_image.x,2,ytol2,&iter,&dy2,minfunc);
  
  p[0]=x_image[0];
  p[1]=x_image[1];
  rayshooterInternal(1,&p);
  
  ResetSourcePlane(tmp_zs,false);
  
  return RAY(p);
}
*/


RAY Lens::find_image(
          Point &p       /// p[] is
          ,PosType ytol2        /// target tolerance in source position squared
          ,PosType &dy2        /// final value of Delta y ^2
          ,bool use_image_guess
){

  if(p.image == nullptr){
    std::cerr << " point in Lens::find_image() must have an attached image point" << std::endl;
  }

  int MaxSteps = 10;
  
  Point_2d yo = *(p.image),dx;

  std::cout << " yo = " << yo << std::endl;
  if(!use_image_guess) p = yo;   // in this case the input image position is not used
  
  rayshooterInternal(1,&p);
  
  double invmago = p.invmag();
  
  Point_2d dy = *(p.image) - yo,dy_tmp;

  if(!use_image_guess){
    // in this case the input image position is not used
    
    p.x[0] = p.x[0] + dy[0];
    p.x[1] = p.x[1] + dy[1];
  
    rayshooterInternal(1,&p);

    dy = *(p.image) - yo;
  }
  std::cout << " dy = " << dy / arcsecTOradians << std::endl;
  
  dy2 = dy.length_sqr();

  int steps = 0;
  double f = 1 , dy2_tmp;
  while(dy2 > ytol2 && steps < MaxSteps){
    
    dx = p.A.inverse() * dy * f;

    p.x[0] = p.x[0] - dx[0];
    p.x[1] = p.x[1] - dx[1];

    rayshooterInternal(1,&p);
    
    if(use_image_guess && p.invmag()*invmago < 0){ // prevents image from changin pairity
      f /= 2;
    }else{

      dy_tmp = *(p.image) - yo;
      std::cout << " dy = " << dy / arcsecTOradians << std::endl;
      
      dy2_tmp = dy.length_sqr();
      
      if(dy2_tmp > dy2){
        f /= 2;
      }else{
        f = 1;
        dy2= dy2_tmp;
        dy = dy_tmp;
      }
      
    }
  }
  
  return RAY(p);
}

RAY Lens::find_image(
          RAY &ray       /// p[] is
          ,PosType ytol2        /// target tolerance in source position squared
          ,PosType &dy2        /// final value of Delta y ^2
          ,bool use_image_guess
                     ){
  LinkedPoint pp;
  pp[0] = ray.x[0];
  pp[1] = ray.x[1];
  pp.image->x[0] = ray.y[0];
  pp.image->x[1] = ray.y[1];
 
  return find_image(pp,ytol2,dy2,use_image_guess);
}
