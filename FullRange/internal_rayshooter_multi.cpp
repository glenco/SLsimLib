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

/** \brief This function calculates the deflection, shear, convergence, rotation
 and time-delay of rays in parallel.  The source redshift must be set for the ray.
 */
void Lens::rayshooter(RAY &ray){
  
  LinkedPoint point;
  point.x[0] = ray.x[0];
  point.x[1] = ray.x[1];
  
  rayshooterInternal(1,&point);
  ray = point;

};

void Lens::rayshooterInternal(
                                unsigned long Npoints
                              , Point *i_points
                              , bool verbose
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
    
    thr[i] = std::thread(&Lens::compute_points_parallel<Point>,this,start,size,i_points,&source_z,false,false);
  }
 
  for(int i = 0; i < nthreads; i++) thr[i].join();
}

void Lens::rayshooterInternal(
                                unsigned long Npoints
                              , LinkedPoint *i_points
                              , std::vector<double> &source_zs
                              , bool verbose
){
  
  // To force the computation of convergence, shear... -----
  // -------------------------------------------------------
  
  if(source_zs.size() != Npoints){
    std::cerr << " Lens::rayshooterInternal() - number of source redshifts must be the same as the number of rays." << std::endl;
    throw std::invalid_argument("missing redshifts");
  }
  
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
  
  // If there are no points to shoot, then we quit.
  if(Npoints == 0) return;
  
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
    
    thr[i] = std::thread(&Lens::compute_points_parallel<LinkedPoint>,this,start,size,i_points
                         ,source_zs.data(),true,false);
  }
 
  for(int i = 0; i < nthreads; i++) thr[i].join();
}

void Lens::rayshooterInternal(
                                unsigned long Npoints
                              , RAY *rays
){
  
  // To force the computation of convergence, shear... -----
  // -------------------------------------------------------
  
  assert(plane_redshifts.size() > 0);
  if(plane_redshifts.size() == 1){  // case of no lens plane
    
    for(int ii = 0; ii < Npoints; ++ii){
    
      rays[ii].y[0] = rays[ii].x[0];
      rays[ii].y[1] = rays[ii].x[1];
      rays[ii].dt = 0;
      rays[ii].A.setToI();
    }
    
    return;
  }
  
  // If there are no points to shoot, then we quit.
  if(Npoints == 0) return;
  
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
    
    thr[i] = std::thread(&Lens::compute_rays_parallel,this,start,size,rays
                         ,false);
  }
 
  for(int i = 0; i < nthreads; i++) thr[i].join();
}

/** \brief Collects information about the halos and kappa contributions along the light path
 
 Information on the nearest halos is collected where nearest is defined by rmax and mode.
 When mode == 2 the unlensed angular coordinates are used to evaluate proximity not the lensed ones.
 
 */
void Lens::info_rayshooter(
                           RAY &ray     /// ray, ray.x needs to be set, all other quantities will be calculated
                           ,std::vector<Point_2d> &ang_positions  /// angular positions on each plane
                           ,std::vector<KappaType> &kappa_on_planes          /// convergence on each plane
                           ,std::vector<std::vector<LensHalo*>> & halo_neighbors  /// neighboring halos within rmax of ray on each plane
                           ,LensHalo &halo_max                /// halo with the larges kappa
                           ,KappaType &kappa_max              /// the kappa from that halo
                           ,KappaType gamma_max[]             /// shear from that halo
                           ,PosType rmax
                           ,int tag
                           ,short mode  /// 0:physical distance (Mpc), 1: comoving distance (Mpc), 2: angular distance (rad)
                           ,bool verbose
)
{
 
  // !!! would like to find the maximum contributing halo so that its contribution can be subtracted from the total
  
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
 
  kappa_max = -100.0;
  gamma_max[0] = 0;
  gamma_max[1] = 0;
  
  // find angular position on first lens plane
  ray.y = ray.x;
  //i_point.image->x[0] = i_point.x[0] ;
  //i_point.image->x[1] = i_point.x[1] ;
  
  theta = ray.y.x;
  theta[0] = ray.x[0];
  theta[1] = ray.x[1];

  // Initializing SumPrevAlphas :
  SumPrevAlphas[0] = theta[0];
  SumPrevAlphas[1] = theta[1];

  // Initializing SumPrevAG :
  SumPrevAG.setToI();
  
  // Setting phi on the first plane.
  phi = 0.0;
  
  // Default values :
  ray.A.setToI();
  ray.dt = 0;
  
 
  // Begining of the loop through the planes :
  // Each iteration leaves i_point[i].image on plane (j+1)
  
  long imax_halo=-1,jmax_halo=-1;
  
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
    
    ang_positions[j][0] = ray.y[0];
    ang_positions[j][1] = ray.y[1];
 
   
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
      
      if(tag==0 || halo_neighbors[j][ii]->tag == tag ){
        if(kappa_tmp > kappa_max){
          kappa_max = kappa_tmp;
          imax_halo=ii;
          jmax_halo=j;
          gamma_max[0] = gamma_tmp[0]/SigmaCrit;
          gamma_max[1] = gamma_tmp[1]/SigmaCrit;
        }
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
    SumPrevAG -= (G * (ray.A)) ;
    
    // Computation of the "plus quantities", i.e. the  next plane quantities :
    ray.A = ray.A * bb + SumPrevAG * aa;
    
  } // End of the loop going through the planes
  
  if(imax_halo !=-1 && jmax_halo != -1)
    halo_max = *(halo_neighbors[jmax_halo][imax_halo]);

  return;
}


void Lens::compute_rays_parallel(int start
                                 ,int chunk_size
                                 ,RAY *rays
                                 ,bool verbose)
{
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
  {
    PosType Dls,Ds;
    FindSourcePlane(rays[0].z,jmax,Dls,Ds);
    Dls_Ds = Dls / Ds;
    if(jmax > 0) D_Ds = Dl[jmax-1] / Ds;
  }
  
  // Main loop : loop over the points of the image
  for(i = start; i < end; i++)
  {
    //theta = i_points[i].image->x;
    theta = rays[i].ptr_y();
    
    theta[0] = rays[i].x[0];
    theta[1] = rays[i].x[1];

    // Initializing SumPrevAlphas :
    SumPrevAlphas[0] = theta[0];
    SumPrevAlphas[1] = theta[1];

    // Initializing SumPrevAG :
    SumPrevAG.setToI();
    
    // Setting phi on the first plane.
    phi = 0.0;
    
    // Default values :
    rays[i].A.setToI();
    rays[i].dt = 0;
    
    // Time delay at first plane : position on the observer plane is (0,0) => no need to take difference of positions.
    rays[i].dt = 0;
    
    //0.5*( p->i_points[i].image->x[0]*p->i_points[i].image->x[0] + p->i_points[i].image->x[1]*p->i_points[i].image->x[1] )/ p->dDl[0] ;
    
    {
      PosType Dls,Ds;
      FindSourcePlane(rays[i].z,jmax,Dls,Ds);
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
      
        if(!flag_switch_deflection_off){
          theta[0] = bb * theta[0] + aa * SumPrevAlphas[0];
          theta[1] = bb * theta[1] + aa * SumPrevAlphas[1];
        }
   
      // ----------------------------------------------------------------------------------------
            
      // Sum_{k=1}^{j} Dl[k] A^k.G^k
      SumPrevAG -= (G * (rays[i].A)) ;
      
      // Computation of the "plus quantities", i.e. the  next plane quantities :
      rays[i].A = rays[i].A * bb + SumPrevAG * aa;
      
      // ----------------------------------------------
      
      // Geometric time delay with added potential
      //p->i_points[i].dt += 0.5*( (xplus[0] - xminus[0])*(xplus[0] - xminus[0]) + (xplus[1] - xminus[1])*(xplus[1] - xminus[1]) ) * p->dTl[j+1] /p->dDl[j+1] /p->dDl[j+1] - phi * p->charge ; /// in Mpc  ???
      
      // Check that the 1+z factor must indeed be there (because the x positions have been rescaled, so it may be different compared to the draft).
      // Remark : Here the true lensing potential is not "phi" but "phi * p->charge = phi * 4 pi G".
      
      
    } // End of the loop going through the planes
    
    if(flag_switch_deflection_off){
      rays[i].A = Matrix2x2<KappaType>::I() - SumPrevAG;
    }
    
    // Subtracting off a term that makes the unperturbed ray to have zero time delay
    //p->i_points[i].dt -= 0.5*( p->i_points[i].image->x[0]*p->i_points[i].image->x[0] + p->i_points[i].image->x[1]*p->i_points[i].image->x[1] ) / p->Dl[NLastPlane];
    
    // Conversion of dt from Mpc (physical Mpc) to Years -----------------
    rays[i].dt *= MpcToSeconds * SecondToYears ;
        
    
    // TEST : showing final quantities
    // ------------------------------=
    if(verbose)  std::cout << "RSI final : X X | " << i << "  " << rays[i].z << " | " << rays[i].kappa() << " " << rays[i].gamma1() << " " << rays[i].gamma2() << " " << rays[i].gamma3() << " " << rays[i].invmag() << " | " << rays[i].dt << std::endl ;
    
  } // End of the main loop.

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
          Point &p              /// p[] is
          ,double zs           ///
          ,PosType ytol2        /// target tolerance in source position squared
          ,PosType &dy2         /// final value of Delta y ^2
          ,bool use_image_guess
){

  if(p.image == nullptr){
    std::cerr << " point in Lens::find_image() must have an attached image point" << std::endl;
  }

  int MaxSteps = 100;
  
  Point_2d yo = *(p.image),dx;

  //std::cout << " yo = " << yo / arcsecTOradians << std::endl;
  if(!use_image_guess) p = yo;   // in this case the input image position is not used
  
  //rayshooterInternal(1,&p);
  compute_points_parallel<Point>(0,1,&p,&zs,true);
  
  double invmago = p.invmag();
  
  Point_2d dy = *(p.image) - yo,dy_tmp;
  if(dy.length_sqr() == 0) return RAY(p,zs);

  if(!use_image_guess){
    // in this case the input image position is not used
    
    p.x[0] = p.x[0] + dy[0];
    p.x[1] = p.x[1] + dy[1];
  
    //rayshooterInternal(1,&p);
    compute_points_parallel<Point>(0,1,&p,&zs,true);

    dy = *(p.image) - yo;
  }
  
  dy2 = dy.length_sqr();
 
  //std::cout << " dy = " << dy / arcsecTOradians << " |dy2| = " << sqrt(dy2) << std::endl;
 
  LinkedPoint pt;
  
  int steps = 0;
  double f = 1 , dy2_tmp,ddy2=ytol2;
  while(dy2 > ytol2 && steps < MaxSteps){
    
    ++steps;
    dx = p.A.inverse() * dy * f;
    while(dx.length_sqr() > dy2*1.0e2 ){ // don't take too big a step at once
      dx /= 2;
      f /= 2;
    }

    pt.x[0] = p.x[0] - dx[0];
    pt.x[1] = p.x[1] - dx[1];

    //rayshooterInternal(1,&pt);
    compute_points_parallel<Point>(0,1,&pt,&zs,true);
    
    if(use_image_guess && pt.invmag()*invmago < 0){ // prevents image from changin pairity
      f /= 2;
    }else{

      dy_tmp = *(pt.image) - yo;
      dy2_tmp = dy.length_sqr();
            
      if(dy2_tmp > dy2){
        f /= 2;
        if(f < 0.01) break;
      }else{
        f = 1;
        dy2= dy2_tmp;
        dy = dy_tmp;
        
        ddy2 = (*(pt.image) -  *(p.image)).length_sqr();
        p = pt;
        //if(ddy2 < ytol2*1.0e-9) break;
      }
      
      //std::cout << " dy = " << dy / arcsecTOradians << " |dy2| = " << sqrt(dy2) << std::endl;
    }
  }
  
  return RAY(p,zs);
}

RAY Lens::find_image(
          RAY &ray             /// p[] is
          ,PosType ytol2       /// target tolerance in source position squared
          ,PosType &dy2        /// final value of Delta y ^2
          ,bool use_image_guess
                     ){
  LinkedPoint pp;
  pp[0] = ray.x[0];
  pp[1] = ray.x[1];
  pp.image->x[0] = ray.y[0];
  pp.image->x[1] = ray.y[1];
 
  return find_image(pp,ray.z,ytol2,dy2,use_image_guess);
}

RAY Lens::find_image(
          Point &p             /// p[] is
          ,double zs           /// source redshift
          ,PosType ytol2        /// target tolerance in source position squared
          ,PosType &dy2        /// final value of Delta y ^2
          ,std::vector<Point_2d> &boundary   /// image will be limited to within this boundary
){

  if(p.image == nullptr){
    std::cerr << " point in Lens::find_image() must have an attached image point" << std::endl;
  }

  int MaxSteps = 100;
  
  Point_2d yo = *(p.image),dx;

  if(! Utilities::inhull(p.x , boundary) ){
    p.x[0] = boundary[0][0];
    p.x[1] = boundary[0][1];
  }

  //rayshooterInternal(1,&p);
  compute_points_parallel<Point>(0,1,&p,&zs,true);
  
  double invmago = p.invmag();
  
  Point_2d dy = *(p.image) - yo,dy_tmp;

  dy2 = dy.length_sqr();
 
  //std::cout << " dy = " << dy / arcsecTOradians << " |dy2| = " << sqrt(dy2) << std::endl;
 
  LinkedPoint pt;
  
  int steps = 0;
  double f = 1 , dy2_tmp,ddy2=ytol2;
  while(dy2 > ytol2 && f > 1.0e-4 && steps < MaxSteps){
    
    ++steps;
    dx = p.A.inverse() * dy * f;
    while(dx.length_sqr() > dy2*1.0e2 ){ // don't take too big a step at once
      dx /= 2;
      f /= 2;
    }

    pt.x[0] = p.x[0] - dx[0];
    pt.x[1] = p.x[1] - dx[1];

    //rayshooterInternal(1,&pt);
    compute_points_parallel<Point>(0,1,&pt,&zs,true);
  
    if(pt.invmag()*invmago < 0){ // prevents image from changing pairity
      f /= 2;
    }else if( !Utilities::inhull(pt.x , boundary) ){
      f /= 2;
    }else{

      dy_tmp = *(pt.image) - yo;
      dy2_tmp = dy_tmp.length_sqr();
            
      if(dy2_tmp > dy2){
        f /= 2;
        //if(f < 0.01) break;
      }else{
        f = 1;
        dy2= dy2_tmp;
        dy = dy_tmp;
        
        ddy2 = (*(pt.image) -  *(p.image)).length_sqr();
        p = pt;
        //if(ddy2 < ytol2*1.0e-9) break;
      }
      
      //std::cout << " dy = " << dy / arcsecTOradians << " |dy2| = " << sqrt(dy2) << std::endl;
    }
  }
  
  return RAY(p,zs);
}

RAY Lens::find_image(
          RAY &ray       /// p[] is
          ,PosType ytol2        /// target tolerance in source position squared
          ,PosType &dy2        /// final value of Delta y ^2
          ,std::vector<Point_2d> &boundary   /// image will be limited to within this boundary
                     ){
  LinkedPoint pp;
  pp[0] = ray.x[0];
  pp[1] = ray.x[1];
  pp.image->x[0] = ray.y[0];
  pp.image->x[1] = ray.y[1];
 
  return find_image(pp,ray.z,ytol2,dy2,boundary);
}
