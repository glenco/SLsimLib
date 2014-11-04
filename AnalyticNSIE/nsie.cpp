/** \file
 * \brief routines for calculating the lensing properties of a non-singular isothermal ellipsoid
*  written by R.B. Metcalf, March 18, 2009
*  based on the analytic solutions of Kormann et al. 1993
*  convention here is gamma_1 = -(Axx-Ayy)/2 and gamma_2= -Axy
*/

#include "base_analens.h"
#include <complex>

/** \ingroup DeflectionL2 \ingroup function
 * \brief Deflection angle for non-singular isothermal ellipsoid in units of Einstein radii
 */
void alphaNSIE(
               PosType *alpha /// Deflection angle in units of the Einstein radius
               ,PosType const *xt   /// position on the image plane in Einstein radius units
               ,PosType f     /// axis ratio of mass
               ,PosType bc    /// core size in same units as alpha
               ,PosType theta /// position angle of ellipsoid
               ){
  
  PosType r=sqrt(xt[0]*xt[0]+xt[1]*xt[1]);
  
  if(f>1.) throw std::runtime_error("f should not be greater than 1 !") ;
                 
  if( r < 1.0e-20 || r > 1.0e20){
	  alpha[0]=0.0;
	  alpha[1]=0.0;
	  return;
  }
  
  // deflection angle alpha has opposite sign with respect to ray position xt
  if( f==1.0 ){
    if(bc == 0.0){
      alpha[0]=-1.0*xt[0]/r;
      alpha[1]=-1.0*xt[1]/r;
    }else{
      alpha[0]=-1.0*(sqrt(r*r+bc*bc) - bc)*xt[0]/r/r;
      alpha[1]=-1.0*(sqrt(r*r+bc*bc) - bc)*xt[1]/r/r;
    }
    return;
  }
  
  PosType x[2],angle[2],fp,b2;
  int s=0;
  
  Utilities::rotation(x,xt,theta);
  
  if(x[0] > 0){ s = -1; x[0] *= -1;}
  else s = 1;
  
  fp=sqrt(1-f*f);
  b2=x[0]*x[0]+f*f*x[1]*x[1];
  
  std::complex<double> c_alpha(0,0),xi(0,0),c_x;
  
  c_x.real(x[0]);
  c_x.imag(x[1]);
  
  std::complex<double> tmp = f*f*c_x*c_x - fp*fp*bc*bc;
  xi = fp/sqrt(tmp);
  
  c_alpha = sqrt(f)*( asinh(xi*sqrt(b2 + bc*bc)) - asinh(xi*bc) )/fp;

  alpha[0] = s*c_alpha.real();
  alpha[1] = -c_alpha.imag();
  
  Utilities::rotation(alpha,alpha,-theta);
  
  //alpha[0] *= -1.0;
  //alpha[1] *= -1.0;
  
  if(alpha[0] != alpha[0] || alpha[1] != alpha[1] ){
	  printf("alpha is %e %e in nsie.c \n fp=%e b2=%e r=%e bc=%e f=%e theta=%e\n x = %e %e xt= %e %e\n"
           ,alpha[0],alpha[1],fp,b2,r,bc,f,theta,x[0],x[1],xt[0],xt[1]);
	  printf("angle=%e %e\n",angle[0],angle[1]);
	  ERROR_MESSAGE();
	  throw std::runtime_error("Invalid input to alphaNSIE");
	  alpha[0]=alpha[1]=0;
  }
  return;
}

/**\ingroup DeflectionL2 \ingroup function
 * \brief Convergence for non-singular isothermal ellipsoid, units \f$ \frac{r_{einstein}}{units(x)} \f$
 * or \f$ \frac{\sigma^2}{\Sigma_{crit}G\, units(xt) } \f$
 */
KappaType kappaNSIE(
		PosType const *xt     /// position on the image plane in Einstein radius units
		,PosType f      /// axis ratio of mass
		,PosType bc     /// core size in units of Einstein radius
		,PosType theta  /// position angle of ellipsoid
		)
{
  PosType x[2],b2;

  Utilities::rotation(x,xt,theta);

  b2=x[0]*x[0]+f*f*x[1]*x[1];
  if( (b2 < 1.0e-20)*(bc < 1.0e-20)){
	  return 1.0e10;
  }
  if(b2>1.0e20 ) return 0.0;
  return 0.5*sqrt(f/(b2+bc*bc));
}

/**\ingroup DeflectionL2 \ingroup function
 * \brief Shear for non-singular isothermal ellipsoid, units \f$ \frac{r_{einstein}}{units(x)} \f$
* or \f$ \frac{\sigma^2}{\Sigma_{crit}G\, units(xt) } \f$
* */
void gammaNSIE(
		KappaType *gam    /// output shear
		,PosType const *xt     /// position on the image plane in Einstein radius units
		,PosType f       /// axis ratio of mass
		,PosType bc      /// core size in units of Einstein radius
		,PosType theta   /// position angle of ellipsoid
		){
  PosType r=sqrt(xt[0]*xt[0]+xt[1]*xt[1]);

  if(r < 1.0e-20 || r > 1.0e20){
	  gam[0]=gam[1]=0.0;
 	  return;
  }

  PosType x[2],fp,P,b2;
  
  Utilities::rotation(x,xt,theta);

  fp=sqrt(1-f*f);

  b2=x[0]*x[0]+f*f*x[1]*x[1];
  //r=sqrt(x[0]*x[0]+x[1]*x[1]);

  P=sqrt(f)*( kappaNSIE(x,f,bc,0.0)*(x[0]*x[0] + f*f*f*f*x[1]*x[1])/sqrt(f)
              - 0.5*(1+f*f)*sqrt(b2+bc*bc)+f*bc ) 
    /( pow(f*r,4) - 2*f*f*fp*fp*bc*bc*(x[0]*x[0]-x[1]*x[1]) + pow(fp*bc,4) );

  gam[0]=(f*f*(x[0]*x[0]-x[1]*x[1])-fp*fp*bc*bc)*P;
  gam[1]=2*f*f*x[0]*x[1]*P;

  Utilities::rotation(gam,gam,-2*theta);


  return;
}
/** \ingroup function
 *  \brief Elliptical radius \f$ R^2 = x^2 + f^2 y^2 \f$ of a NonSingular Isothermal Ellipsoid
 */

// TODO: BEN Check pi factor against notes.
PosType rmaxNSIE(
		PosType sigma    /// velocity dispersion in km/s
		,PosType mass    /// mass in Msun
		,PosType f       /// axis ratio
		,PosType rc      /// core radius Mpc
		){
	return sqrt( pow(mass*Grav*lightspeed*lightspeed*f/pi/sigma/sigma + rc,2) - rc*rc );
}
/** \ingroup function
 *  \brief Elliptical radius \f$ R^2 = x^2 + f^2 y^2 \f$ given f and position angle of model
 */
PosType ellipticRadiusNSIE(
		PosType const *x
		,PosType f
		,PosType pa
		){
	PosType y[2];

	Utilities::rotation(y,x,pa);

	return sqrt(y[0]*y[0] + f*f*y[1]*y[1]);
}
     /* inverse magnification */
KappaType invmagNSIE(
		PosType *x
		,PosType f
		,PosType bc
		,PosType theta
		,KappaType *gam
		,KappaType kap
		){

  gammaNSIE(gam,x,f,bc,theta);
  kap=kappaNSIE(x,f,bc,theta);
  return pow(1-kap,2) - gam[0]*gam[0] - gam[1]*gam[1];
}

namespace Utilities{
	/// Rotates 2 dimensional point without changing input point
	void rotation(float *xout,float *xin,PosType theta){
    
    float tmp = xin[0];  // to make it work if xout == xin
		xout[0]=tmp*cos(theta)-xin[1]*sin(theta);
		xout[1]=xin[1]*cos(theta)+tmp*sin(theta);
	}
	/// Rotates 2 dimensional point without changing input point
	void rotation(PosType *xout,PosType const *xin,PosType theta){
    
    PosType tmp = xin[0];  // to make it work if xout == xin
		xout[0]=tmp*cos(theta)-xin[1]*sin(theta);
		xout[1]=xin[1]*cos(theta)+tmp*sin(theta);
	}
}



/**\ingroup function
 *
 * Structure that does allow the integration of alphaNSIE in phiNSIE.
 *
 */
struct alphaForInt {
  
  alphaForInt(PosType f, PosType bc, PosType theta) : f(f), bc(bc), theta(theta) {}

  KappaType operator()(PosType r) /// PosType r : |position| on the image plane in Einstein radius units
  {
    PosType alpha[2] ;
    PosType rmin = 1.e-7 ;
    PosType x[2] = {rmin,r} ;  // rmin is arbitrary, any other value should pass the assertion as long as we do not get nan.
    
    alphaNSIE(alpha,x,f,bc,theta) ;
    
    if(f==1.)     // Spherical case
    {
      if(r < rmin) { r = rmin ;
        // And jump over assertion !
      }
      else { // std::cout << "alpha[0]/x[0] = " << alpha[0]/x[0] << " , alpha[1]/x[1] = " << alpha[1]/x[1] << std::endl ;
      assert( (alpha[0]/x[0])/(alpha[1]/x[1]) - 1. < 1.e-7 );   // Works only in symmetric case !!!
      }
    }
    else          // Elliptical case
    {
      // std::cout << "We have to write the code for phi in the NSIE elliptical case." << std::endl ;
      // exit(0);
    }
    
    return (alpha[0]/x[0]) * r ;
  }
  
  private :
  PosType f;       /// axis ratio of mass
  PosType bc;      /// core size in units of Einstein radius
  PosType theta;   /// position angle of ellipsoid
};



/**\ingroup function
 *
 * Compute the potential for the NSIE (in physical Mpc) by integration of alphaNSIE.
 *
 */
KappaType LensHaloBaseNSIE::phiNSIE(PosType const *xt    /// position on the image plane in Einstein radius units
                  ,PosType f      /// axis ratio of mass
                  ,PosType bc     /// core size in units of Einstein radius
                  ,PosType theta  /// position angle of ellipsoid
                )
{
  // Here f, bc, and theta are fixed and known.

  PosType r = sqrt(xt[0]*xt[0]+xt[1]*xt[1]) ;

  // std::cout << "xt = (" << xt[0] << " , " << xt[1] << " ) " << std::endl ;
  // std::cout << fratio << "  " << beta << "  " << pa << std::endl ;
  // std::cout << f << "  " << bc << "  " << theta << std::endl ;
  
  
  set_slope(1.); // sets the slope beta. // HAVE TO GENERALISE THAT !
  
  
  if(fratio==1) // Spherical case
  {
  
    struct alphaForInt alphaForIntFunc(f,bc,theta);
    
    // Returning phi from the integration of alphaForIntFunc :
    // =======================================================
    return Utilities::nintegrate<alphaForInt,PosType>(alphaForIntFunc, 1.0e-7, r, 1.0e-7);
    
  }
  else          // Elliptical case
  {
    
    // Calculating the modes necessary to ellipticize phi : OLD METHOD -- DOES NOT WORK !
    // ====================================================

    //for(int i=0;i<perturb_Nmodes;i++) { std::cout << perturb_modes[i] << "  " ; }
    //std::cout << std::endl ;
    
    // PosType tmp_pert_modes[Nmod];
    // calcModes(f, 2-beta, theta, tmp_pert_modes);
    
    // main_perturbation beta  -> perturb_beta
    // main_perturbation kappa -> perturb_rms[0]
    // main_perturbation gamma -> perturb_rms[1]
    // main_perturb_monopole   -> perturb_rms[2]
    // main_perturb_quadrapole -> perturb_rms[3]
    // main_perturb_hexopole   -> perturb_rms[4]
    // main_perturb_octopole   -> perturb_rms[5]
    
    // perturb_modes[0] = 1. ;  // perturbation kappa
    // perturb_modes[1] = 0.03;  // perturbation gamma1
    // perturb_modes[2] = 0.03;  // perturbation gamma2
    // perturb_modes[3] = 0.0;   // Monopole
    // perturb_modes[4] = 0.005; // Quadrupole
    // perturb_modes[5] = 0.005; // Hexapole
    // perturb_modes[6] = 0.01;  // Octupole
    
    // std::cout << "perturbation modes extracted from phiNSIE : " << std::endl ;
    // for(int i=0;i<perturb_Nmodes;i++) { std::cout << perturb_modes[i] << "  " ; }
    // std::cout << std::endl << std::endl ;
    
    

    
    // Computing phi from the method in lens_expand :
    // ==============================================
    
    PosType tmp_theta, cosx, sinx, cos2theta, sin2theta ;
    PosType F, F1, F2;
    int i, k;
    
    tmp_theta = atan2(xt[1],xt[0]);
    cosx = xt[0]/r;
    sinx = xt[1]/r;

    if(perturb_Nmodes > 3) F=0.5*perturb_modes[3];
    else F = 0;
    F1=0;
    F2=0;
    for(i=4; i<perturb_Nmodes; i+=2)
    {
      k=i/2;
      F += perturb_modes[i]*cos(k*theta)     + perturb_modes[i+1]*sin(k*theta);
      F1+=-perturb_modes[i]*k*sin(k*theta)   + perturb_modes[i+1]*k*cos(k*theta);
      F2+=-perturb_modes[i]*k*k*cos(k*theta) - perturb_modes[i+1]*k*k*sin(k*theta);
    }
    cos2theta=2*cosx*cosx-1;
    sin2theta=2*cosx*sinx;
    
    return F*pow(r,beta) + r*r*(perturb_modes[0] + perturb_modes[1]*cos2theta + perturb_modes[2]*sin2theta)/2;
    
    
    // Test :
    // Compare kappa, gamma, alpha as computed with lens_expand with kappaNSIE, gammaNSIE, alphaNSIE
    // in order to validate the elliptisation of phi.
    /*
    bool test = false ;
    
    KappaType kappa_tmp, gamma_tmp, phi_tmp;
    kappa_tmp = gamma_tmp = phi_tmp = 0.;
    PosType * alpha_tmp = new PosType ;
    
    // Computing all quantities with lens_expand :
    kappa_tmp = lens_expand(beta, tmp_pert_modes, perturb_Nmodes, xt, alpha_tmp, &gamma_tmp, &phi_tmp);
    
    
    // Computing alpha with NSIE function :
    PosType * alpha_tmpNew = new PosType ;
    alphaNSIE(alpha_tmpNew, xt, f, beta, theta);
    
    // Computing kappa with NSIE function :
    KappaType kappa_tmpNew = kappaNSIE(xt, f, beta, theta);
    
    // Showing comparison :
    std::cout << "Comparison alpha1 : " << alpha_tmp[0] << "  " << alpha_tmpNew[0] << std::endl;
    std::cout << "Comparison alpha2 : " << alpha_tmp[1] << "  " << alpha_tmpNew[1] << std::endl;
    std::cout << "Comparison kappa : " << kappa_tmp << "  " << kappa_tmpNew << std::endl;
    
    // if(alpha_tmp[0] == alpha_tmp[1] && alpha_tmp[1] == alpha_tmpNew[1]) test = true ;
    
    
    // Returning phi given by lens_expand :
    // ====================================
    
    if (test == true) {
      return F*pow(r,beta) + r*r*(mod[0] + mod[1]*cos2theta + mod[2]*sin2theta)/2;
    }
    else {
      return 0. ;
    }
    */
    
  }

}



/**\ingroup function
 *
 * Quadropole moment of an elliptically truncated NSIE
 * Units are unit[mass]*unit[Rmax]^2
 */
void quadMomNSIE(
		float mass     /// total mass
		,float Rmax    /// elliptical maximum radius
		,float f       /// axis ratio of mass
		,float rc      /// core size in same units as Rmax
		,float theta   /// position angle of ellipsoid
		,PosType *quad   /// output
	){

	PosType m3,b;
	b = rc/Rmax;
	m3 = f*Rmax*Rmax*mass*(1-f*f)*( (1-2*b*b)*sqrt(1+b*b) +2*b*b*b)/(sqrt(1+b*b)-b)/6/f/f;

	quad[0] = m3*cos(-2*theta);
	quad[1] = -m3*cos(-2*theta);
	quad[2] = m3*sin(-2*theta);

	return;
}
