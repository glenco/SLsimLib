#include "slsimlib.h"
#include <nrutil.h>
#include <iomanip>      // std::setprecision

using namespace std;

static double betaT,*modT,**xobT,**dx_subT,sigGT,*modTT,*modoT,**vT,x_centerT[2],**xgT,**dx_subTt;
static int NmodT,NsourcesT,NimagesT,*pairingT,degenT,Nmin;
static double oldsm;//,tang[2],length,yot[2],radsourceT;

/** 
 *
 *  \brief Wrapper that allows simple lens to be found with a single
 * lens with a single source and translates result into data structures used in the other code.
 *
 * The lens is centered on [0,0] source position in lens is updated along with all the modes.
 *
 */

void LensHaloFit::FindLensSimple(
                                 int Nimages               /// Number of images to be fit
                                 ,Point *image_positions   /// Array of points with point[i].x set to the image positions
                                 ,double *y                /// output source position
                                 ,double **dx_sub          /// dx_sub[Nimages][2] pre-calculated deflections caused by substructures or external masses at each image
){
  
  ImageInfo* imageinfo = new ImageInfo[Nimages];
  
  for(int i=0;i<Nimages;++i){
    imageinfo[i].centroid[0] = image_positions[i].x[0];
    imageinfo[i].centroid[1] = image_positions[i].x[1];
  }
  
  FindLensSimple(imageinfo,Nimages,y,dx_sub);
  
  std::cout << "perturbation modes (in LensHaloFit::FindLensSimple)" << std::endl;
  for(int i=0;i<perturb_Nmodes;++i) std::cout << perturb_modes[i] << " " ;
  std::cout << std::endl;
  
  delete[] imageinfo;
}



/** 
 *
 *  \brief Same as FindLensSimple but with some tests in it.
 *
 *
 */

bool LensHaloFit::SafeFindLensSimple(
                                 int Nimages                      /// Number of images to be fit
                                 ,Point *image_positions          /// Array of points with point[i].x set to the image positions
                                 ,double *y                       /// output source position
                                 ,double **dx_sub                 /// dx_sub[Nimages][2] pre-calculated deflections caused by substructures or external masses at each image
                                 ,int SafetyNum                   /// integer of the number of time you want to check the modes
                                 ,PosType PixelSizeRad            /// Pixel size in radians used for the maps
                                 ,std::vector<std::vector<PosType>> & PrecisionBackTracedPos  /// Table that will receive the back-traced images uncertainty (in units of PixelSizeRad).
                                 ,std::vector<std::vector<PosType>> & alphaTab                /// Table that will receive the deviation angle in radians
                                 ,bool verbose                    /// verbose mode switch
){
  // Array to store the modes over the different calls of FindLensSimple :
  PosType ModesMin [perturb_Nmodes] ;         // Minimal value for the modes (among all calls)
  PosType ModesMax [perturb_Nmodes] ;         // Maximal value for the modes (among all calls)
  PosType ModesAve [perturb_Nmodes] ;         // Average value for the modes (among all calls)
  const double ToleranceModes = 0.01 ;        // Tolerance on the ratio (Max-Min)/Average for the modes
  // way 1 :
  // const double ToleranceSourcePos = 0.1 ;  // Tolerance on the ratio (y-x)/alpha on the source position reconstruction
  // way 2 :
  const double ToleranceSourcePos = 1. ;      // Tolerance on the difference y-x+alpha , with respect to the pixel size, on the source position reconstruction
  bool ReturnCode = true ;                    // Boolean that returns whether or not the Fit was a success.
  
  
  // Doing the proper initialisation of these quantities :
  
  for(int k=0;k<perturb_Nmodes;++k)
  {
    ModesMin[k] = 1.e200 ;    // So the modes have to be less than that !
    ModesMax[k] = -1.e200 ;   // So the modes have to be more than that !
    ModesAve[k] = 0. ;
  }
  
  
  // We compute the modes SafetyNum times and retain the min, the max, and the sum :
  // ------------------------------------------------------------------------------=
  
  // We call FindLensSimple many times :
  for(int i=0;i<SafetyNum;i++)
  {
    // Calling FindLensSimple (the one that really computes the modes) :
    ////////////////////////////////////////////////
    FindLensSimple(Nimages,image_positions,y,dx_sub);
    ////////////////////////////////////////////////

    // Test that no 'nan' occurs :
    for(int k=0;k<perturb_Nmodes;++k)
    {
      assert(perturb_modes[k] == perturb_modes[k]) ;
    }
    
    // Filling the check tables :
    for(int k=0;k<perturb_Nmodes;++k)
    {
      if(perturb_modes[k]<ModesMin[k]) ModesMin[k] = perturb_modes[k];
      if(perturb_modes[k]>ModesMax[k]) ModesMax[k] = perturb_modes[k];
      ModesAve[k] += perturb_modes[k]; // Summing the modes
    }
  }
  // Dividing by the number of values to get the average :
  for(int k=0;k<perturb_Nmodes;++k) ModesAve[k] /= perturb_Nmodes ;
  
  // We check that the modes are not > 1e50 (which is clearly wrong) :
  PosType UpperBoundForModes = 1.e50 ;
  for(int k=0;k<perturb_Nmodes;++k)
  {
    if(abs(ModesAve[k])>UpperBoundForModes)
    {
      std::cout << "Mode k = " << k << " : Av = " << abs(ModesAve[k]) << " > " << UpperBoundForModes << std::endl ;
      ReturnCode = false ;
    }
  }
  if(ReturnCode == false)
  {
    ERROR_MESSAGE();
    std::cout << "Error of unstability in SafeFindLensSimple !" << std::endl ;
    // exit(0);
  }
  
  // We compare the min and max values with the average :
  if(ReturnCode == true)
  {
    for(int k=0;k<perturb_Nmodes;++k)
    {
      if(abs(ModesMax[k]-ModesMin[k]) > abs(ToleranceModes*ModesAve[k]))
      {
        std::cout << "Mode k = " << k << " : Max-Min = " << abs(ModesMax[k]-ModesMin[k]) << " < " << "Tol * Av = " << abs(ToleranceModes*ModesAve[k]) << std::endl ;
        ReturnCode = false ;
      }
    }
    if(ReturnCode == false)
    {
      ERROR_MESSAGE();
      std::cout << "Error of unstability in SafeFindLensSimple !" << std::endl ;
      // exit(0);
    }
  }

  // If the modes are crazy then we return the false value here :
  if(ReturnCode == false)
  {
    std::cout << "Not doing the check of the back-traced images because of the unstability in SafeFindLensSimple !" << std::endl ;
    return ReturnCode ;
  }
  
  // Printing values otherwise :
  if(verbose)
  {
    std::cout << "Values of the modes :" << std::endl ;
    std::cout << "Ave : " ;
    for(int k=0;k<perturb_Nmodes;++k) std::cout << std::setprecision(7) << perturb_modes[k] << " " ;
    std::cout << std::endl << "Min : " ;
    for(int k=0;k<perturb_Nmodes;++k) std::cout << std::setprecision(7) << ModesMin[k] << " " ;
    std::cout << std::endl << "Max : " ;
    for(int k=0;k<perturb_Nmodes;++k) std::cout << std::setprecision(7) << ModesMax[k] << " " ;
    std::cout << std::endl ;
    
    std::cout << std::endl << "Estimation made with " << SafetyNum << " calls of FindLensSimple, with a tolerance of " << ToleranceModes*100. << " % on the modes." << std::endl;
    
    // Caption :
    // q[1] : ellipticity of nearest elliptical lens
    // q[2] : position angle of nearest elliptical lens
    // if sigG > 0 :
    // q[3,4] : center of lens, copied to x_center
    
    cout << "Ellipticity : " << qpriv[1] << " , position angle : " << qpriv[2] << std::endl ;
    if (sigGT > 0) cout << "Center of the lens : " << qpriv[3] << " , " << qpriv[4] << std::endl ;
  }
  

  // Testing that the image positions are consistent when traced back to the source plane :
  // ------------------------------------------------------------------------------------==
  
  std::cout << std::endl << "Evaluated position of the source from each image :" << std::endl ;
  
  // alpha that we are going to use after to get the source position form the back-traced image position :
  PosType * alphaTMP = new PosType [2];
  alphaTMP[0] = 0. ; alphaTMP[1] = 0. ;
  // Quantities not used for the moment :
  double betaTMP = perturb_beta ;
  KappaType * gammaTMP = new KappaType [3];
  KappaType * phiTMP = new KappaType ;
  KappaType kappaTMP;
  gammaTMP[0] = 0. ; gammaTMP[1] = 0. ; gammaTMP[2] = 0. ;
  *phiTMP = 0. ;
  kappaTMP = 0. ;
  
  // Ratios for (y-x)/alpha :
  PosType ratioSourcePos [2] ;
  ratioSourcePos[0] = ratioSourcePos[1] = 0. ;

  // Displaying the back-traced positions :
  for(int i=0;i<Nimages;++i)
  {
    Point pointTMP ;
    pointTMP.x[0] = image_positions[i].x[0] * Dl ;
    pointTMP.x[1] = image_positions[i].x[1] * Dl ;
    
    kappaTMP = lens_expand(betaTMP, perturb_modes, perturb_Nmodes, pointTMP.x, alphaTMP, gammaTMP, phiTMP); // We removed the -1 after perturb_Nmodes AS IT SHOULD NOT BE THERE but it is still used this way in LensHaloBaseNSIE::force_halo !
    
    // Applying the factors of LensHaloBaseNSIE::force_halo :
    alphaTMP[0] *= -1. ;
    alphaTMP[1] *= -1. ;
    
    // Applying the extra factors of RayShooter :
    alphaTMP[0] *= -4*PI*Grav ;
    alphaTMP[1] *= -4*PI*Grav ;
    alphaTMP[0] *= Dls / Ds ;
    alphaTMP[1] *= Dls / Ds ;
    
    // Saving alpha in the argument :
    alphaTab[i][0] = alphaTMP[0];
    alphaTab[i][1] = alphaTMP[1];
    
    // Computing ratios :
    ratioSourcePos[0] = abs(alphaTMP[0]) / abs(y[0] - image_positions[i].x[0]) ;
    ratioSourcePos[1] = abs(alphaTMP[1]) / abs(y[1] - image_positions[i].x[1]) ;
    
    // Displaying the values :
    if(verbose)
    {
      std::cout << "Image : " << image_positions[i].x[0] << " " << image_positions[i].x[1] << std::endl ;
      std::cout << "y - x : " << y[0] - image_positions[i].x[0] << " " << y[1] - image_positions[i].x[1] << std::endl ;
      std::cout << "alpha : " << alphaTMP[0] << " " << alphaTMP[1] << std::endl ;
      std::cout << "ratios : " << ratioSourcePos[0] << " " << ratioSourcePos[1] << std::endl ;
    }
    
    // Saving the relative precision of the back-traced images in the source plane (in units of pixels of the map) :
    PrecisionBackTracedPos[i][0] = abs(y[0] - image_positions[i].x[0] + alphaTMP[0]) / PixelSizeRad ;
    PrecisionBackTracedPos[i][1] = abs(y[1] - image_positions[i].x[1] + alphaTMP[1]) / PixelSizeRad ;
    
    // Deciding if the test is sufficient to keep on with the rest :
    if(verbose)
    {
      std::cout << "(y-x+alpha)/pixelsize : " << PrecisionBackTracedPos[i][0] << " , " << PrecisionBackTracedPos[i][1] << " pixels." << std::endl ;
      std::cout << "! y = x - alpha : " << image_positions[i].x[0] - alphaTMP[0] << " , " << image_positions[i].x[1] - alphaTMP[1] << " !" << std::endl << std::endl ;
    }
    for(int k=0;k<2;k++)
    {
      // In the case where we are above the tolerance on the source position & the real source position is not too close to (0,0) :
      // way 1 :
      // if( abs(ratioSourcePos[k]-1) > ToleranceSourcePos && abs(y[k] - image_positions[i].x[k]) > 1.e-15 ) // assuming positions to be in radians.
      // way 2 :
      if( abs(y[k] - image_positions[i].x[k] + alphaTMP[k]) > PixelSizeRad * ToleranceSourcePos && abs(y[k] - image_positions[i].x[k]) > 1.e-15 ) // assuming positions to be in radians.
      {
        ERROR_MESSAGE();
        std::cout << "Error of precision in source-position reconstruction in SafeFindLensSimple !" << std::endl ;
        ReturnCode = false ;
        // exit(0);
      }
      // Else we can continue !
    }
    if(!verbose) std::cout << "Back-traced source position : " << image_positions[i].x[0] - alphaTMP[0] << " " << image_positions[i].x[1] - alphaTMP[1] << std::endl << std::endl ;
    
  }

  // OTHER TESTS ?
  // ------------=
  
  // Otherwise we keep the last computed modes and display them :
  std::cout << std::endl << "Perturbation modes (in LensHaloFit::FindLensSimple) :" << std::endl;
  for(int i=0;i<perturb_Nmodes;++i) std::cout << perturb_modes[i] << " " ;
  std::cout << std::endl;
  
  return ReturnCode ;
}


/** 
 *
 *  \brief Wrapper that allows simple lens to be found with a single
 * lens with a single source and translates result into data structures used in the other code.
 *
 * The lens is centered on [0,0] source position in lens is updated along with all the modes.
 *
 */

void LensHaloFit::FindLensSimple(
                                 ImageInfo *imageinfo    /// Positions of images relative to center of lens.  Only imageinfo[].centoid[] is used. Centroids must be in radians.
                                 ,int Nimages             /// input number of images
                                 ,double *y               /// output source position
                                 ,double **dx_sub         /// dx_sub[Nimages][2] pre-calculated deflections caused by substructures or external masses at each image (in radians)
){
  
  assert(Nimages < 200);
  
  if(perturb_Nmodes <= 0 || Nimages <= 0){
    ERROR_MESSAGE();
    std::printf("must set perturb_Nmodes lens->perturb_Nmodes = %i Nimages = %i \n"
                ,perturb_Nmodes,Nimages);
    exit(1);
  }
  
  int i;
  
  //PrintAnaLens(lens,false,false);
  //for(i=0;i<Nimages;++i) std::printf("  x = %e %e\n",image_positions[i].x[0],image_positions[i].x[1]);
  
  if(Nimages == 1){
    for(i=1;i<perturb_Nmodes;++i) perturb_modes[i] = 0.0;
    y[0] = imageinfo[0].centroid[0]; // (radians)
    y[1] = imageinfo[0].centroid[1]; // (radians)
    
    return ;
  }
  
  int pairing[Nimages],Nsources = 1;
  double **xob,**xg,q[6],*mods;
  double re2 = 0,x_center[2],scale;
  //for(int i=0;i<6;i++) q[i] = 0.;
  //  !!! need to set initial guess to something
  
  xob = dmatrix(0,Nimages-1,0,1); // For the rescaled positions of the images
  xg = dmatrix(0,1,0,1);
  mods = dvector(0,perturb_Nmodes + 2*Nsources + 1 );
  
  xg[0][0] = xg[0][1] = 0.0;
  x_center[0] = x_center[1] = 0.0;
  
  // calculate scale to re-normalize (in radians). Otherwise the linear algebra routines will fail.
  for(i=0,scale=0;i<Nimages;++i) scale = DMAX(scale,sqrt( pow(imageinfo[0].centroid[0] - imageinfo[i].centroid[0],2) + pow(imageinfo[0].centroid[1] - imageinfo[i].centroid[1],2) ) );
  
  for(i=0;i<Nimages;++i){
    pairing[i] = 1;
    xob[i][0] = imageinfo[i].centroid[0]/scale; // xob is rescaled here (i.e. values = xob / xobmax ~ 1)
    xob[i][1] = imageinfo[i].centroid[1]/scale; // i.e. xob is dimensionless
    
    dx_sub[i][0] /= scale; // same with dx_sub now in units of scale
    dx_sub[i][1] /= scale; // in units of scale
    //std::printf("xob = %e %e  dx = %e %e\n",xob[i][0],xob[i][1],dx_sub[i][0],dx_sub[i][1]);
  }
  
  x_center[0] /= scale; // x_center now in units of scale
  x_center[1] /= scale;
  
  //ERROR_MESSAGE();
  ElliptisizeLens(Nimages,Nsources,1,pairing,xob,x_center,xg,0,perturb_beta,perturb_Nmodes
                  ,mods,dx_sub,&re2,q); // The -1 in after perturb_Nmodes WAS MAKING FINDLENSSIMPLE UNSTABLE !
  
  // At this point we should have :
  // mod[0] = 0
  // mod[1] and mod[2] : in units of scale (i.e. dimensionless)
  // mod[3] and further : in units of scale^(beta-1) = 1 (i.e. no units) for beta = 1
    
  // Assigning the modes :
  for(int j=1;j<perturb_Nmodes;++j) perturb_modes[j] = mods[j];
  cout << "# " ;
  for(int j=1;j<perturb_Nmodes;++j) cout << perturb_modes[j] << " " ;
  cout << endl ;
  
  perturb_modes[0] = 0.0;
  perturb_modes[1] *= -1;  // checked
  perturb_modes[2] *= -1;  // checked
  
  // source position :
  y[0] = mods[perturb_Nmodes+1]*scale;
  y[1] = mods[perturb_Nmodes+2]*scale;
  // y[0] and y[1] are now in radians.
    
  // std::cout << "i = " << i << std::endl ;
  std::cout << "scale = " << scale << std::endl;
  std::cout << "source : y[0] = " << y[0] << " , y[1] = " << y[1] << std::endl;
  
  // For convenience :
  //PosType zl = LensHalo::getZlens() ;
  //PosType zs = zsource_reference ;
    
  // Converting source position to physical angle :
  // y[0] *= Dl * (1+zl) / (Ds * (1+zs)) ;
  // y[1] *= Dl * (1+zl) / (Ds * (1+zs)) ;
  // y[0] and y[1] are still in radians.
  
  // dx_sub and x_center back in radians :
  for(i=0;i<Nimages;++i)
  {
    dx_sub[i][0] *= scale;
    dx_sub[i][1] *= scale;
  }
  x_center[0] *= scale;
  x_center[1] *= scale;
  
  //Einstein_ro = 0.0; // the monople is now included in the modes
  //sigma = 0.0;
    
    
  // ----- Applying factors to the modes ------------------------------------------------------
    
  // Multiplying the first 3 modes by scale :
  for(i=3;i<perturb_Nmodes;i++) perturb_modes[i] *= scale ; // Important step !
  // mod[0,1,2] are already in radians (see in ElliptisizeLens).
  // mod[3,4,5,...] are now in radians.
    
  for(i=3;i<perturb_Nmodes;i++) perturb_modes[i] *= Dl ;
  // mod[0,1,2] are in radians.
  // mod[3,4,5,...] are now in PhysMpc.

  // Check that Dl*(1+zl) + Dls*(1+zs) = Ds*(1+zs), i.e. that D*(1+z) are here the comoving distances :
  //assert(Dl*(1+zl) + Dls*(1+zs) - Ds*(1+zs) == 0.);
  // std::cout << "Dl (1+zl) + Dls (1+zs) = " << Dl*(1+zl) + Dls*(1+zs) << " , Ds (1+zs) = " << Ds*(1+zs) << std::endl ;
  
  for(i=0;i<perturb_Nmodes;i++) perturb_modes[i] /= (4*PI*Grav * Dls * Dl / Ds) ;
  // mod[0,1,2] are now in radians / ((PhysMpc / mass) * PhysMpc) = mass / PhysMpc^2.
  // mod[3,4,5,...] are now in PhysMpc / ((PhysMpc / mass) * PhysMpc) = mass / PhysMpc for beta = 1.
    
  // ------------------------------------------------------------------------------------------
    
    
  // ------------------------------------------------
  // So when the modes go out they are in the units :
  // mod[0,1,2] in mass / PhysMpc^2.
  // mod[3,4,5,...] in mass / PhysMpc.
  // ------------------------------------------------
    
    
  // Print modes :
  // for(i=0;i<perturb_Nmodes;i++) std::cout << perturb_modes[i] << " " ;
  // std::cout << std::endl ;
    
  free_dmatrix(xob,0,Nimages-1,0,1);
  free_dmatrix(xg,0,1,0,1);
  free_dvector(mods,0,perturb_Nmodes + 2*Nsources + 1);
  
  // storing q :
  for(int i = 0 ; i < 7 ; i++) qpriv[i] = q[i] ;
  
  return ;
}

// OPERATIONS ON THE MODES BEFORE I DISCOVER THE PROBLEM :
// for(i=3;i<perturb_Nmodes;i++) perturb_modes[i] *= scale ;
// for(i=3;i<perturb_Nmodes;i++) perturb_modes[i] *= Dl ;
// for(i=0;i<perturb_Nmodes;i++) perturb_modes[i] /= (4*PI*Grav * Dls * (1+zs)) ;





// test solutions (removed from the code above) :
/*
 {
 
 PosType alpha[2];
 KappaType gamma[2],phi;
 
 for(int i=0 ; i<4 ; i++)
 {
 xob[i][0] *= scale ;
 xob[i][1] *= scale ;
 }
 
 // std::cout << "/// In HaloFit ///" << std::endl ;
 // std::cout << "test FindLensSimple solution" << std::endl;
 // std::cout << "perturb_beta = " << perturb_beta << " , perturb_Nmodes = " << perturb_Nmodes << std::endl ;
 // std::cout << "Modes in FindLensSimple before calling lens_expand :" << std::endl ;
 // for(int i=0 ; i<perturb_Nmodes ; i++) std::cout << perturb_modes[i] << " " ;
 // std::cout << std::endl ;
 
 
 // Here mod[0,1,2] are dimensionless and mod[4,5,...] are in radians so that alpha is in radians !
 for(int i=0;i<Nimages;++i){
 lens_expand(perturb_beta,perturb_modes,perturb_Nmodes-1,xob[i],alpha,gamma,&phi); // The -1 is needed here too !
 
 
 std::cout << "xob : @@@ " << xob[i][0] << " " << xob[i][1] << " @@@" << std::endl ;
 std::cout << "alpha in FindLensSimple : ??? " << alpha[0] << "  " << alpha[1] << " ???" << std::endl;
 // std::cout << "source in FindLensSimple : !!! " << xob[i][0] - alpha[0] << "  " << xob[i][1] - alpha[1] << " !!!" << std::endl;
 
 }
 }
 */




/** L2
 * \brief Find most elliptical lens
 *
 ******
 ***  Output:
 ***
 *   - mod[]     - [1...Nmod+2*Nsources]  modes of final model
 *   - re2       - Einstein radius of additional lens
 *      - q[1]     - ellipticity of nearest elliptical lens
 *      - q[2]     - position angle of nearest elliptical lens
 *   - if sigG == 0
 *      - q[3]     - re2 of additional lens if Nlenses = 2, copied to re2
 *   - if sigG > 0
 *      - q[3,4]   - center of lens, copied to x_center
 *      - q[5]     - re2 of additional lens if Nlenses = 2, copied to re2
 *      - q[6]  ???
 *
 *************************************
 **************************************/
double LensHaloFit::ElliptisizeLens(
                                    int Nimages   /// number of images
                                    ,int Nsources /// number of sources
                                    ,int Nlenses  /// number of lens centers
                                    ,int *pairing  /// [0...Nimages-1] which source each image belongs to, indexed 1 through Nsources
                                    ,double **xob  /// [0...Nimages-1][0...1]  observed image positions (in units of scale when used in LensHaloFit::FindLensSimple)
                                    ,double *x_center  /// x_center[][] - ??? expected center of lens (in units of scale when used in LensHaloFit::FindLensSimple)
                                    ,double **xg  /// ??? [0...Nlenses-1][0...1] centers of additional lenses
                                    ,double sigG  /// amount by which the center of the lens is allowed to vary
                                    ,double beta  /// - slope of density profile  kappa propto r^-beta
                                    ,int Nmod    /// number of modes used in lens fitting, this must be greater then ???
                                    ,double *mod  /// Axial modes of lens model (to be computed here)
                                    ,double **dx_sub  /// [0,Nimages-1][0,1] deflection offset for each image caused by any masses not included in the host model, ex substructure (in units of scale when used in LensHaloFit::FindLensSimple)
                                    ,double *re2     /// Einstein radius of additional lens (to be computed here)
                                    ,double *q       /// output, see comment (to be computed here)
){
  
  int iter,i;
  double **xi,sm;
  static int count=0;
  double s;
  
  assert(Nimages > 1);
  assert(Nsources > 0);
  assert(Nlenses > 0);
  for(i=0;i<Nimages;++i) assert(pairing[i] > 0 && pairing[i] <= Nimages);
  
  //std::printf("x_center = %e %e\n",x_center[0],x_center[1]);
  //std::printf("sigG = %e\n",sigG);
  
  ++count;
  
  // allocate memory
  pairingT=ivector(0,Nimages-1);
  xobT=dmatrix(0,Nimages-1,0,1);
  dx_subT=dmatrix(0,Nimages-1,0,1);
  modT=dvector(1,Nmod+2*Nsources+1);
  modTT=dvector(1,Nmod+2*Nsources+1);
  modoT=dvector(1,Nmod+2*Nsources+1);
  vT=dmatrix(1,Nmod+2*Nsources+1,1,Nmod+2*Nsources+1);
  xgT=dmatrix(0,Nlenses-1,0,1);
  if(Nlenses>1) dx_subTt=dmatrix(0,Nimages-1,0,1);
  
  // copy input parameters into global variables
  NimagesT=Nimages;
  NsourcesT=Nsources;
  betaT=beta;
  NmodT=Nmod;
  sigGT=sigG;
  
  // clean modes to zero (in case it was not done before calling the function)
  for(i=1;i<=Nmod+2*Nsources+1;++i) mod[i]=0.0;
  
  // copying the input objects
  for(i=0;i<Nimages;++i){
    xobT[i][0] = xob[i][0];
    xobT[i][1] = xob[i][1];
    pairingT[i] = pairing[i];
    dx_subT[i][0] = dx_sub[i][0];
    dx_subT[i][1] = dx_sub[i][1];
    if(Nlenses>1 && count>1 ){
      s=atan2(xg[1][1]-xob[i][1],xg[1][0]-xob[i][0]);
      dx_subT[i][0] = dx_sub[i][0] + *re2*cos(s);
      dx_subT[i][1] = dx_sub[i][1] + *re2*sin(s);
    }
  }
  for(i=0;i<Nlenses;++i){
    xgT[i][0] = xg[i][0];
    xgT[i][1] = xg[i][1];
  }
  
  // find lens with degeneracy information
  find_lens(NimagesT,NsourcesT,pairingT,xobT,x_center,betaT,NmodT,&degenT,modT,vT,dx_subT);
  
  // Looking inside find_lens, we have :
  // mod[0] : not touched, it is zero
  // mod[1] and mod[2] : in units of scale when used in LensHaloFit::FindLensSimple
  // mod[3] and further : in units of scale^(beta-1) = 1 (i.e. no units) when used in LensHaloFit::FindLensSimple and for beta = 1
    
  //std::printf("found model\n");
  for(i=1;i<=Nmod+2*Nsources+1;++i) mod[i] = modT[i];
  
  /*
   * find the most elliptical model amongst the degenerate models that
   * fit the image positions
   */
  //return sm;
  //if(count == 1){
		q[1] = 0.9;
		q[2] = 0.0; /*find_axis(modT,NmodT);*/
		//}
  
  x_centerT[0] = x_center[0];
  x_centerT[1] = x_center[1];
  
  xi=dmatrix(1,5,1,5);
  
  xi[1][1]=0.01;  xi[1][2]=0.00;    xi[1][3]=0.0;       xi[1][4]=0.0;        xi[1][5]=0.0;
  xi[2][1]=0.00;  xi[2][2]=1.0e-3;  xi[2][3]=0.0;       xi[2][4]=0.0;        xi[2][5]=0.0;
  xi[3][1]=0.00;  xi[3][2]=0.0;     xi[3][3]=1.0e-3;    xi[3][4]=0.0;        xi[3][5]=0.0;
  xi[4][1]=0.00;  xi[4][2]=0.0;     xi[4][3]=0.0;       xi[4][4]=1.0e-3;     xi[4][5]=0.0;
  xi[5][1]=0.00;  xi[5][2]=0.0;     xi[5][3]=0.0;       xi[5][4]=0.0;        xi[5][5]=1.0e-3;
  
  oldsm=-1;
  if(sigG == 0.0){   /// keep center of lens fixed
    if(Nlenses==1 || count>1){
      Nmin=2;
      powellD(q,xi,Nmin,1.0e-18,&iter,&sm,minEllip); // This should not affect the units of the modes.
    }else{
      Nmin=3;
      if(count == 1) q[3]=1.0e-5;
      powellD(q,xi,Nmin,1.0e-18,&iter,&sm,minEllip);
      *re2=q[3];
    }
  }else{
    q[3] = x_center[0];
    q[4] = x_center[1];
    if(Nlenses==1){
      Nmin=4;
      powellD(q,xi,Nmin,1.0e-18,&iter,&sm,minEllip);
      *re2=0.0;
    }else{
      Nmin=5;
      if(count == 1) q[5]=0.001;
      powellD(q,xi,Nmin,1.0e-18,&iter,&sm,minEllip);
      *re2=q[5];
    }
    x_center[0] = q[3];
    x_center[1] = q[4];
  }
  for(i=1;i<=Nmod+2*Nsources+1;++i) mod[i]=modTT[i];
  
  //std::printf("iter = %i\n",iter);
  free_ivector(pairingT,0,Nimages-1);
  free_dmatrix(xobT,0,Nimages-1,0,1);
  free_dmatrix(dx_subT,0,Nimages-1,0,1);
  free_dvector(modT,1,Nmod+2*Nsources+1);
  free_dvector(modTT,1,Nmod+2*Nsources+1);
  free_dvector(modoT,1,Nmod+2*Nsources+1);
  free_dmatrix(vT,1,Nmod+2*Nsources+1,1,Nmod+2*Nsources+1);
  free_dmatrix(xgT,0,Nlenses-1,0,1);
  if(Nlenses>1) free_dmatrix(dx_subTt,0,Nimages-1,0,1);
  free_dmatrix(xi,1,5,1,5);
  
  return sm;
}

/// minimized to find model closest to elliptical
double minEllip(double *par){
  double K,E,q,theta,sm,r,s;
  int i;
  static double x_center[2];
  
  q=par[1];   // ellipticity
  theta=par[2];  // position angle
  
  // set modo to elliptical model
  for(i=1;i<=NmodT;++i){
    modoT[i]=0;
  }
  
  // elliptical integrals
  K = rfD(0,1./q/q,1);
  E = K - (1-1./q/q)*rdD(0,1./q/q,1)/3;
  
  // fill in modes with their values for an elliptical lens
  modoT[3]=4*K/PI;
  if(q != 1.0){
    if(NmodT>3) modoT[4] =4*( (1+q*q)*K-2*q*q*E )/(1-q*q)/PI/(1-4);
    if(NmodT>7) modoT[8] =4*( (3*q*q+1)*(q*q+3)*K-8*q*q*(1+q*q)*E )
      /( 3*PI*pow(1-q*q,2) )/(1-16);
    if(NmodT>11) modoT[12] =4*( (1+q*q)*(15+98*q*q+15*q*q*q*q)*K-2*q*q*(23+82*q*q+23*q*q*q*q)*E )
      /( 15*PI*pow(1-q*q,3) )/(1-36);
    if(NmodT>15) modoT[16]=4*( -32*q*q*(1+q*q)*(11+74*q*q+11*q*q*q*q)*E
                              +(105+1436*q*q+3062*q*q*q*q+1436*pow(q,6)+105*pow(q,8))*K )
      /(105*PI*pow(1-q*q,4))/(1-64);
  }
  
  // rotate model
  RotateModel(theta,modoT,NmodT,NsourcesT);   // xobT is used as a dumby variable
  
  x_center[0]=x_centerT[0];
  x_center[1]=x_centerT[1];
  
  if(Nmin == 3){
    for(i=1;i<=NimagesT;++i){  // add in second lens
      s=atan2(xgT[1][1]-xobT[i-1][1],xgT[1][0]-xobT[i-1][0]);
      dx_subTt[i-1][0]=dx_subT[i-1][0] + par[3]*cos(s);
      dx_subTt[i-1][1]=dx_subT[i-1][1] + par[3]*sin(s);
    }
    find_lens(NimagesT,NsourcesT,pairingT,xobT,x_centerT,betaT,NmodT,&degenT,modT,vT,dx_subTt);
  }
  if(Nmin == 4){
    if(par[3] != x_center[0] || par[4] != x_center[1]){
      x_center[0]=par[3];
      x_center[1]=par[4];
      find_lens(NimagesT,NsourcesT,pairingT,xobT,x_center,betaT,NmodT,&degenT,modT,vT,dx_subT);
    }
  }
  if(Nmin == 5){  /** add in second lens **/
    if(par[3] != x_center[0] || par[4] != x_center[1]){
      x_center[0]=par[3];
      x_center[1]=par[4];
    }
    for(i=1;i<=NimagesT;++i){  /** add in second lens **/
      s=atan2(xgT[1][1]-xobT[i-1][1],xgT[1][0]-xobT[i-1][0]);
      dx_subTt[i-1][0]=dx_subT[i-1][0] + par[5]*cos(s);
      dx_subTt[i-1][1]=dx_subT[i-1][1] + par[5]*sin(s);
    }
    find_lens(NimagesT,NsourcesT,pairingT,xobT,x_center,betaT,NmodT,&degenT,modT,vT,dx_subTt);
  }
  
  // find most elliptical model
  sm=regularize(NmodT,3,NmodT,NsourcesT,degenT,modT,vT,modoT);
  
  if( sm < oldsm || oldsm < 0.0){
    oldsm=sm;
    for(i=1;i<=NmodT+2*NsourcesT+1;++i) modTT[i]=modT[i];
  }
  //std::printf("q=%e   theta=%e gamma=%e %e\n",q,theta*180/PI,modT[1],modT[2]);
  //std::printf("%e %e %e %e %e %e\n",modoT[3],modoT[4],modoT[5],modoT[6],modoT[7],modoT[8]);
  //std::printf("%e %e %e %e %e %e\n\n",modT[3],modT[4],modT[5],modT[6],modT[7],modT[8]);
  
  // keep orientation within range
  if(theta > PI/2) sm *= exp( (theta - PI/2)/1.0e-5 );
  if(theta < -PI/2 ) sm *= exp( -(theta - PI/2)/1.0e-5 );
  // keep center of lens within sigG
  if(Nmin == 4 || Nmin == 5){
    r=pow(x_center[0]-x_centerT[0],2) + pow(x_center[1]-x_centerT[1],2);
    if(sqrt(r) > sigGT ) sm *= exp( 10*(r-sigGT*sigGT)/sigGT/sigGT );
  }
  if(fabs(log10(q)) > 2) sm *=exp( pow( (fabs(log10(q))- 2)/0.0001   ,2) );
  
  //sm *= exp( (modT[1]*modT[1] + modT[2]*modT[2])/pow(0.01,2) );
  return sm;
}

/** L2
 *
 * \brief Calculates a lens that fits the image positions
 *
 *  dx_sub[] is a perturbation to the deflection angle calculated elsewhere
 *
 * - [1]=gamma1
 * - [2]=gamma2
 * - [3]=ao
 *******************************************************/
//void AnaNSIELensHalo::find_lens(int Nimages,int Nsources,int *pairing,double **xob,double *x_center,double beta
void find_lens(int Nimages,int Nsources,int *pairing,double **xob,double *x_center,double beta,int Nmodes,int *degen,double *mod,double **v,double **dx_sub){
  double **c,*b,*w,r,theta,wmax,**a,*y,*temp,**x;
  int i,k,j;
  
  if(Nmodes+2*Nsources < 2*Nimages){ std::printf("ERROR: too few parameters in find_lens\n"); exit(0);}
  
  y=dvector(0,1);
  temp=dvector(0,1);
  
  
   x=dmatrix(0,Nimages-1,0,1);
   c=dmatrix(1,2*Nimages+1,1,Nmodes+2*Nsources+1);
   a=dmatrix(1,2*Nimages+1,1,Nmodes+2*Nsources+1);
   
  b=dvector(1,2*Nimages+1);
  w=dvector(1,Nmodes+2*Nsources+1);
  /*
  x= Utilities::PosTypeMatrix(0,Nimages-1,0,1);
  c= Utilities::PosTypeMatrix(1,2*Nimages+1,1,Nmodes+2*Nsources+1);
  a= Utilities::PosTypeMatrix(1,2*Nimages+1,1,Nmodes+2*Nsources+1);
  */
  betaT = beta;
  NmodT = Nmodes;
  
  // recenter on lens center
  for(i=0;i<Nimages;++i){
    x[i][0]=xob[i][0]-x_center[0]; x[i][1]=xob[i][1]-x_center[1];
    //std::printf("%e  %e\n",x[i][0],x[i][1]);
  }
  
  // fill in data matrix
  for(i=1;i<=Nimages;++i){
    r = sqrt( x[i-1][0]*x[i-1][0] + x[i-1][1]*x[i-1][1] );
    theta = atan2(x[i-1][1],x[i-1][0]);
    
    // shear from mod[1] and mod[2]
    c[2*i-1][1]= -x[i-1][0]; c[2*i-1][2]= -x[i-1][1];
    c[2*i][1]  =  x[i-1][1]; c[2*i][2]  = -x[i-1][0];
    
    // monopole from mod[3]
    c[2*i-1][3]= 0.5*beta*pow(r,beta-1)*cos(theta);
    c[2*i][3]  = 0.5*beta*pow(r,beta-1)*sin(theta);
    
    // contributions from higher modes
    for(j=4;j<=Nmodes-1;j+=2){
      k=j/2;
      c[2*i-1][j]  = (beta*cos(theta)*cos(k*theta)+k*sin(theta)*sin(k*theta))*pow(r,beta-1);
      c[2*i-1][j+1]= (beta*cos(theta)*sin(k*theta)-k*sin(theta)*cos(k*theta))*pow(r,beta-1);
      
      c[2*i][j]  = (beta*sin(theta)*cos(k*theta)-k*cos(theta)*sin(k*theta))*pow(r,beta-1);
      c[2*i][j+1]= (beta*sin(theta)*sin(k*theta)+k*cos(theta)*cos(k*theta))*pow(r,beta-1);
    }
    
    // assign images to sources
    for(j=0;j<Nsources;++j){
      if(pairing[i-1]==j+1){
        c[2*i-1][2*j+Nmodes+1] = 1; c[2*i-1][2*j+Nmodes+2] = 0;
        c[2*i][2*j+Nmodes+1]   = 0; c[2*i][2*j+Nmodes+2]   = 1;
      }else{
        c[2*i-1][2*j+Nmodes+1] = 0; c[2*i-1][2*j+Nmodes+2] = 0;
        c[2*i][2*j+Nmodes+1]   = 0; c[2*i][2*j+Nmodes+2]   = 0;
      }
    }
    b[2*i-1] = x[i-1][0] + dx_sub[i-1][0];
    b[2*i]   = x[i-1][1] + dx_sub[i-1][1];
    /*std::printf("dx_sub = %e %e\n",dx_sub[i-1][0],dx_sub[i-1][1]);*/
  }
  
  // preserve image data matrix
  for(i=1;i<=2*Nimages;++i){
    for(j=1;j<=Nmodes+2*Nsources;++j){
      a[i][j] = c[i][j];
    }
  }
  
  //**** fit largest modes to first four images  ***
  // SVD decomposition
  svdcmpD(c,2*Nimages,Nmodes+2*Nsources,w,v);
  // weed out small w's
  wmax=0.0;
  *degen=0;
  for(i=1;i<=Nmodes+2*Nsources;++i) if (w[i] > wmax) wmax=w[i];
  wmax=1.0e-6*wmax;
  for(i=1;i<=Nmodes+2*Nsources;++i) if (w[i] < wmax){ w[i]=0.0; ++*degen;}  // removes degenerate modes
  //ERROR_MESSAGE();
  svbksbD(c,w,v,2*Nimages,Nmodes+2*Nsources,b,mod);
  
  //for(i=1;i<=Nmodes + 2*Nsources;++i) std::printf("mod[%i]=%e\n",i,mod[i]);
  //for(i=1;i<=Nmodes + 2*Nsources;++i) std::printf("w[%i]=%e\n",i,w[i]);
  
  /* find degeneracy vectors
   * and make them the first *degen columns of v
   * mod[i] + a_1 v[i][1] + ... + a_degen v[i][degen] is a solution
   **/
  
  for(i=1,j=0;i<=Nmodes+2*Nsources+1;++i){
    if (w[i] == 0.0){
      ++j;
      for(k=1;k<=Nmodes+2*Nsources+1;++k){
        v[k][j]=v[k][i];
      }
    }
  }
  
  //********************************************************************
  //for(i=1;i<=N+2*Nsources;++i) std::printf("%e\n",modo[i]);
  
  free_dvector(y,0,1);
  free_dvector(temp,0,1);
  
   free_dmatrix(x,0,Nimages-1,0,1);
   free_dmatrix(c,1,2*Nimages+1,1,Nmodes+2*Nsources+1);
   free_dmatrix(a,1,2*Nimages+1,1,Nmodes+2*Nsources+1);
   
  /*Utilities::free_PosTypeMatrix(x,0,Nimages-1,0,1);
  Utilities::free_PosTypeMatrix(c,1,2*Nimages+1,1,Nmodes+2*Nsources+1);
  Utilities::free_PosTypeMatrix(a,1,2*Nimages+1,1,Nmodes+2*Nsources+1);
  */
  free_dvector(b,1,2*Nimages+1);
  free_dvector(w,1,Nmodes+2*Nsources+1);
  
  return ;
}


/**
 *   \brief Sets the perturbation modes in the LensHaloFit
 *
 */
void LensHaloFit::set_perturbmodes(PosType * ListModes, const int Nmodes)
{
  perturb_Nmodes = Nmodes ;
  for(int i=0; i < Nmodes; i++) perturb_modes[i] = ListModes[i];
  return ;
}

/**
 *   \brief Gets the perturbation modes in the LensHaloFit
 *
 */
void LensHaloFit::get_perturbmodes(PosType * ListModes, const int Nmodes)
{
  assert(Nmodes == perturb_Nmodes);
  
  for(int i=0 ; i < perturb_Nmodes ; i++)
  {
    ListModes[i] = perturb_modes[i] ;
  }
  return ;
};


/**
 *   \brief Gets the perturbation modes in the LensHaloFit
 *
 */
void LensHaloFit::get_perturbmodes(std::vector<PosType> & ListModes)
{
  assert(ListModes.size() == perturb_Nmodes);
  
  for(int i=0 ; i < perturb_Nmodes ; i++)
  {
    ListModes[i] = perturb_modes[i] ;
  }
  return ;
};



/** L2
 *
 * \brief  calculate the sources position, surface density and magnification at x
 * given lens model
 *        - mod given
 *        - Re2 - Einstein radius of second lens
 *        - x2 - position of second lens relative to center of lens
 ***************************************************************/

double LensHaloAnaNSIE::deflect_translated(double beta,double *mod,double *x,double *y,double *mag,int Nmodes
                                           ,int Nlenses,double Re2,double *x2){
  KappaType kappa,gamma[2];
  KappaType *phi = new KappaType ;
  
  if(mod[0] != 0.0){ERROR_MESSAGE(); std::printf("background kappa should be zero\n"); exit(0);}
  assert(Nlenses < 3);
  
  // use deflection calculator to reduce code duplication
  kappa = lens_expand(beta,mod,Nmodes,x,y,gamma,phi); // Shouldn't it be Nmodes-1 ?
  
  // translate result to convention used here
  
  /// changed from alpha to y,  also convention on shear is opposite
  y[0] = x[0] + 2*(x[0]*mod[1] + x[1]*mod[2]) - y[0];
  y[1] = x[1] - 2*(x[1]*mod[1] + x[0]*mod[2]) - y[1];
  
  mag[0] = 1 - kappa - gamma[0];
  mag[1] = 1 - kappa + gamma[0];
  mag[3] = -gamma[1];
  
  /*
   theta=atan2(x[1],x[0]);
   r=sqrt(x[0]*x[0] + x[1]*x[1]);
   cosx=x[0]/r;
   sinx=x[1]/r;
   */
  
  /*
   F=0.5*mod[3];
   F1=0;
   F2=0;
   for(i=4;i<=Nmodes-1;i+=2){
   k=i/2;
   F += mod[i]*cos(k*theta)     + mod[i+1]*sin(k*theta);
   F1+=-mod[i]*k*sin(k*theta)   + mod[i+1]*k*cos(k*theta);
   F2+=-mod[i]*k*k*cos(k*theta) - mod[i+1]*k*k*sin(k*theta);
   }
   
   
   y[0] = -pow(r,beta-1)*(beta*cosx*F - sinx*F1);
   y[1] = -pow(r,beta-1)*(beta*sinx*F + cosx*F1);
   
   // add shear
   y[0] += x[0] + x[0]*mod[1] + x[1]*mod[2];
   y[1] += x[1] - x[1]*mod[1] + x[0]*mod[2];
   */
  
  
  /*
   // magnification matrix in polar coordinates
   dxdr= (1+mod[1])*cosx + mod[2]*sinx
   - (beta-1)*pow(r,beta-2)*(beta*cosx*F-sinx*F1);
   dxda=-(1+mod[1])*r*sinx + mod[2]*r*cosx
   + pow(r,beta-1)*(beta*sinx*F + (1-beta)*cosx*F1 + sinx*F2);
   dydr= (1-mod[1])*sinx + mod[2]*cosx
   - (beta-1)*pow(r,beta-2)*(beta*sinx*F+cosx*F1);
   dyda= (1-mod[1])*r*cosx - mod[2]*r*sinx
   + pow(r,beta-1)*(-beta*cosx*F + (1-beta)*sinx*F1 - cosx*F2);
   */
  
  
  /*   if(Nlenses>1){ */
  /*     dxdr= cosx -Re2*( sinx2*sinx2*cosx - sinx2*cosx2*sinx )/r2; */
  /*     dydr= sinx -Re2*( cosx2*cosx2*sinx - sinx2*cosx2*cosx )/r2; */
  /*     dxda= -r*sinx + r*Re2*( sinx2*sinx2*sinx + sinx2*cosx2*cosx )/r2; */
  /*     dyda=  r*cosx - r*Re2*( cosx2*cosx2*cosx + sinx2*cosx2*sinx )/r2; */
  /*   } */
  
  /** actually the inverse magnification **/
  /* *mag=(dxdr*dyda - dxda*dydr)/r; */
  
  // add additional lens
  if(Nlenses > 1){
    double dxdr,dxda,dydr,dyda;
    double theta2,r2,cosx2,sinx2,r,cosx,sinx;
    
    //theta=atan2(x[1],x[0]);
    r=sqrt(x[0]*x[0] + x[1]*x[1]);
    cosx=x[0]/r;
    sinx=x[1]/r;
    
    theta2=atan2(x2[1]-x[1],x2[0]-x[0]);
    r2=sqrt( pow(x[0]-x2[0],2) + pow(x[1]-x2[1],2));
    cosx2=cos(theta2);
    sinx2=sin(theta2);
    
    y[0] += Re2*cos(theta2);
    y[1] += Re2*sin(theta2);
    
    dxdr =  -Re2*( sinx2*sinx2*cosx - sinx2*cosx2*sinx )/r2;
    dydr =  -Re2*( cosx2*cosx2*sinx - sinx2*cosx2*cosx )/r2;
    dxda = r*Re2*( sinx2*sinx2*sinx + sinx2*cosx2*cosx )/r2;
    dyda =-r*Re2*( cosx2*cosx2*cosx + sinx2*cosx2*sinx )/r2;
    
    mag[0] += ( -x[1]*dxda + r*x[0]*dxdr )/r/r;  /** xx **/
    mag[1] += (  x[0]*dyda + r*x[1]*dydr )/r/r;  /** yy **/
    mag[2] += (  x[0]*dxda + r*x[1]*dxdr )/r/r;  /** xy **/
    
    return kappa + 0.5*Re2/r2;
  }
  
  return kappa;
}
/** L2
 * \brief find degenerate model most like modo modulo a normalization
 */
double regularize(int Nmax,int Nmin,int N,int Nsources,int degen
                  ,double *mod,double **v,double *modo){
  double Dsum,sum=0,sumold,aa,*weights;
  

  int i,j;
  
  /*
   for(i=Nmin,sum=0.0;i<=Nmax;i+=1){
   std::printf("%i %e %e\n",i,mod[i],modo[i]);
   }
   */
  
  weights=dvector(1,degen);
  
  /*   aa=mod[3]/modo[3]; */
  /*   for(i=3;i<=Nmax;++i) modo[i]*=aa; */
  
  for(i=Nmin,Dsum=0;i<=Nmax;i+=1){
    if(i<=3) Dsum += pow( mod[i]-modo[i],2);
    else Dsum += pow(1-pow(i/2,2.),2)*pow( mod[i]-modo[i],2);
  }
  sumold=Dsum;
  
  while(Dsum > 1.0e-6*sumold){
    
    if(modo[3] != 0.0){
      /** renormalize model **/
      for(i=3,sum=0.0;i<=Nmax;++i){
        if(i<=3) sum += mod[i]*modo[i];
        else  sum += pow(1-pow(i/2,2.),2)*mod[i]*modo[i];
      }
      aa=sum;
      for(i=3,sum=0.0;i<=Nmax;++i){
        if(i<=3) sum += modo[i]*modo[i];
        else  sum += pow(1-pow(i/2,2.),2.)*modo[i]*modo[i];
      }
      aa/=sum;
      
      if(aa > 0.0) for(i=3;i<=Nmax;++i) modo[i]*=aa;
    }
    /** move in degenerate space to find best model **/
    for(j=1;j<=degen;++j){
      for(i=Nmin,sum=0.0;i<=Nmax;i+=1){
        if(i<=3) sum += ( mod[i] - modo[i]  )*v[i][j];
        else sum += pow(1-pow(i/2,2.),2.)*( mod[i] - modo[i]  )*v[i][j];
      }
      weights[j]=-sum;
      for(i=Nmin,sum=0.0;i<=Nmax;i+=1){
        if(i<=3) sum += v[i][j]*v[i][j];
        else sum += pow(1-pow(i/2,2.),2)*v[i][j]*v[i][j];
      }
      weights[j] /= sum;
      
      for(i=1;i<=N+2*Nsources+1;++i){
        //std::printf("i=%i j=%i w=%e v=%e\n",i,j,weights[j],v[i][j]);
        if(weights[j] == weights[j]) mod[i] += weights[j]*v[i][j];
      }
      if(sum < 0) std::printf("max found\n");
      /*std::printf("weights[%i]=%e\n",j,weights[j]);*/
    }
    
    for(i=Nmin,sum=0.0;i<=Nmax;i+=1){
      if(i<=3) sum += pow(mod[i]-modo[i],2);
      else sum += pow(1-pow(i/2,2.),2.)*pow(mod[i]-modo[i],2);
    }
    Dsum=fabs(sumold-sum);
    /*std::printf("%e\n",sumold-sum);*/
    sumold=sum;
    /*std::printf("Dsum=%e %e\n",Dsum,sumold);*/
  }
  
  /* test result */
  free_dvector(weights,1,degen);
  
  return sum;
}

/*
 double critmin(double r){
 double kappa,mag[3],x[2],y[2],AA[4],angle[2],ang_lens[2];
 
 ang_lens[0]=0;
 ang_lens[1]=0;
 x[0]=r*cos(thetaT)+xt[0];
 x[1]=r*sin(thetaT)+xt[1];
 
 kappa=deflect_translated(1.0,modT,x,y,mag,NmodT,NlensesT,Re2T,x2T);
 if(withsub==1){
 deflection_total(x,angle,1,1,ang_lens,1,AA);
 return (mag[0]+AA[0]-1)*(mag[1]+AA[1]-1)-(mag[2]+AA[2])*(mag[2]+AA[3]);
 }else{ return mag[0]*mag[1]-mag[2]*mag[2];}
 }
 */

/**************************************************
 ************** rotate lens model *****************
 **************************************************/
void RotateModel(double thetaX,double *mod,int N,int Nsources){
  double temp[2];
  int i,k;
  /*  std::printf("rotate: theta = %e\n",thetaX);*/
  
  temp[0]=mod[1];  temp[1]=mod[2];
  mod[1]= temp[0]*cos(2*thetaX)+temp[1]*sin(2*thetaX);
  mod[2]=-temp[0]*sin(2*thetaX)+temp[1]*cos(2*thetaX);
  for(i=4;i<=N-1;i+=2){
    k=i/2;
    temp[0]=mod[i]; temp[1]=mod[i+1];
    mod[i]  = temp[0]*cos(k*thetaX)+temp[1]*sin(k*thetaX);
    mod[i+1]=-temp[0]*sin(k*thetaX)+temp[1]*cos(k*thetaX);
  }
  
  /** source positions **/
  for(i=0;i<Nsources;++i){
    temp[0]=mod[N+1+2*i];  temp[1]=mod[N+2+2*i];
    mod[N+1+2*i]= temp[0]*cos(thetaX)+temp[1]*sin(thetaX);
    mod[N+2+2*i]=-temp[0]*sin(thetaX)+temp[1]*cos(thetaX);
  }
}


double LensHaloAnaNSIE::find_axis(double *mod,int Nmod){
  /********************************************
   ************* find rotation axis ***********
   ********************************************/
  double theta,ans;
  int i,k;
  
  NmodT=Nmod;
  for(i=1;i<=NmodT;++i) modT[i] = mod[i];
  theta = zriddrD(minaxis,0,PI/4.,1.0e-9);
  
  /* calc 2nd deriv */
  for(i=4,ans=0.0;i<=5;i+=2){
    k=i/2;
    /*std::printf("ans=%e   %e %e\n",ans,modT[i],modT[i+1]);*/
    ans += -4*k*k*modT[i]*modT[i+1]*sin(2*k*theta)
    - 2*k*k*(modT[i]*modT[i] - modT[i+1]*modT[i+1])*cos(2*k*theta);
  }
  
  if(ans < 0) theta += PI/4;  /* make answer a minimum */
  return theta;
}

double minaxis(double thetaX){
  double ans;
  int i,k;
  
  for(i=4,ans=0.0;i<=5;i+=2){
    k=i/2;
    /*std::printf("ans=%e   %e %e\n",ans,modT[i],modT[i+1]);*/
    ans += 2*k*modT[i]*modT[i+1]*cos(2*k*thetaX)
    - k*(modT[i]*modT[i] - modT[i+1]*modT[i+1])*sin(2*k*thetaX);
  }
  return ans;
}


int LensHaloAnaNSIE::check_model(int Nimages,int Nsources,int Nlenses,int *pairing,double **xob,double *x_center
                                 ,int Nmod,double *mod,double **xg,double Re2,double **dx_sub,double **Amag,double ytol){
  int i;
  double dy,dymax,tmp[2],kappa,mag[3],x2[2],par;
  
  
  // check parity of images
  for(i=0,par=0;i<Nimages;++i) par+=(Amag[i][0]*Amag[i][1]-Amag[i][2]*Amag[i][2])
    /fabs(Amag[i][0]*Amag[i][1]-Amag[i][2]*Amag[i][2]);
  
  //std::printf("par=%e\n",par);
  if(fabs(par) > 1.0e-6 ) return 1;
  
  if(Nlenses>1){
    x2[0]=xg[1][0]-x_center[0]; x2[1]=xg[1][1]-x_center[1];
  }else{ Re2=0.0;}
  
  // check that there is only one source
  double **y = dmatrix(0,Nimages-1,0,1);
  for(i=0;i<Nimages;++i){
    tmp[0]=xob[i][0]-x_center[0];
    tmp[1]=xob[i][1]-x_center[1];
    kappa=deflect_translated(1.0,mod,tmp,y[i],mag,Nmod,Nlenses,Re2,x2);
    y[i][0]+=dx_sub[i][0];
    y[i][1]+=dx_sub[i][1];
    //std::printf("y= %e %e    %e %e  %i\n",y[i][0],y[i][1],kappa,1.0/(mag[0]*mag[1]-mag[2]*mag[2]),pairing[i-1]);
  }
  
  for(i=1,dymax=0;i<Nimages;++i){
    dy=sqrt( pow(y[0][0]-y[i][0],2)+pow(y[0][1]-y[i][1],2) );
    if(dy > dymax) dymax=dy;
  }
  
  free_dmatrix(y,0,Nimages-1,0,1);
  
  if(dymax > ytol) return 1;
  
  // check number of images
  tmp[0]=xob[1][0]-x_center[0];
  tmp[1]=xob[1][1]-x_center[1];
  //std::printf("number of images: %i\n",find_image_number(y[0],tmp,mod,Nmod,Nlenses,Re2,x2) );
  //std::printf("number of images: %i\n",magandcrit(y[0],tmp,mod,Nmod,Nlenses,Re2,x2) );
  //std::printf("finite mag: %e\n",finiteMag(1.0e-3,xob[0],mod,Nmod,Nlenses,Re2,x2) );
  //std::printf("finite mag: %e\n",finiteMag(1.0e-3,xob[1],mod,Nmod,Nlenses,Re2,x2) );
  //std::printf("finite mag: %e\n",finiteMag(1.0e-3,xob[2],mod,Nmod,Nlenses,Re2,x2) );
  //std::printf("finite mag: %e\n",finiteMag(1.0e-3,xob[3],mod,Nmod,Nlenses,Re2,x2) );
  
  return 0;
}


// Copied from lens_expand.c
PosType LensHaloFit::lens_expand(PosType beta,PosType *mod,int Nmodes,PosType const *x,PosType *alpha,KappaType *gamma,KappaType *phi)
{
  PosType F,F1,F2,theta,r,cosx,sinx,cos2theta,sin2theta,gt,gx;
  int i,k;
  
  if(Nmodes<=0){
    alpha[0]=alpha[1]=0;
    gamma[0]=gamma[1]=0;
    return 0.0;
  }
  
  r=sqrt(x[0]*x[0] + x[1]*x[1]);
  theta=atan2(x[1],x[0]);
  cosx=x[0]/r;
  sinx=x[1]/r;
  
  if(Nmodes > 3) F=0.5*mod[3];
  else F = 0;
  F1=0;
  F2=0;
  for(i=4;i<Nmodes;i+=2){
    k=i/2;
    F += mod[i]*cos(k*theta)     + mod[i+1]*sin(k*theta);
    F1+=-mod[i]*k*sin(k*theta)   + mod[i+1]*k*cos(k*theta);
    F2+=-mod[i]*k*k*cos(k*theta) - mod[i+1]*k*k*sin(k*theta);
  }
  
  alpha[0] = pow(r,beta-1)*(beta*cosx*F - sinx*F1);
  alpha[1] = pow(r,beta-1)*(beta*sinx*F + cosx*F1);
  
  // add shear
  alpha[0] +=  x[0]*mod[1] + x[1]*mod[2];
  alpha[1] += -x[1]*mod[1] + x[0]*mod[2];
  
  // add flat kappa
  alpha[0] += -1.0*x[0]*mod[0];
  alpha[1] += -1.0*x[1]*mod[0];
  
  gt=-0.5*pow(r,beta-2)*(beta*(beta-2)*F-F2);
  gx=pow(r,beta-2)*(beta-1)*F1;
  
  cos2theta=2*cosx*cosx-1;
  sin2theta=2*cosx*sinx;
  
  gamma[0]=-(gt*cos2theta+gx*sin2theta) + mod[1];
  gamma[1]=-gt*sin2theta+gx*cos2theta  + mod[2];
  
  //gamma[0]=-0.5*( x[0]*(r*dxdr - dyda) - x[1]*(r*dydr + dxda) )/r/r;
  //gamma[1]=-(  x[0]*dxda + r*x[1]*dxdr )/r/r;
  
  // potential
  *phi = F*pow(r,beta) + r*r*(mod[0] + mod[1]*cos2theta + mod[2]*sin2theta)/2;
  
  //printf("  lens_expand *phi = %e\n",*phi);
  return 0.5*(beta*beta*F+F2)*pow(r,beta-2) + mod[0];
}


