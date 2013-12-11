/** \file
 * \brief routines for calculating the lensing properties of a non-singular isothermal ellipsoid
*  written by R.B. Metcalf, March 18, 2009
*  based on the analytic solutions of Kormann et al. 1993
*  convention here is gamma_1 = -(Axx-Ayy)/2 and gamma_2= -Axy
*/

#include "base_analens.h"

/** \ingroup DeflectionL2 \ingroup function
 * \brief Deflection angle for non-singular isothermal ellipsoid in units of Einstein radii
 */
void alphaNSIE(
		PosType *alpha /// Deflection angle in units of the Einstein radius
		,PosType *xt   /// position on the image plane in Einstein radius units
		,PosType f     /// axis ratio of mass
		,PosType bc    /// core size in same units as alpha
		,PosType theta /// position angle of ellipsoid
		){

  PosType r=sqrt(xt[0]*xt[0]+xt[1]*xt[1]);

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

  PosType x[2],angle[2],fp,b2,RCphase,SCphase;

  Utilities::rotation(x,xt,theta);

  fp=sqrt(1-f*f);
  b2=x[0]*x[0]+f*f*x[1]*x[1];

  PosType Qp=( pow(fp*sqrt(b2+bc*bc)+x[0],2) + f*f*f*f*x[1]*x[1] )
    /( pow(f*r*r+fp*bc*x[0],2)+fp*fp*bc*bc*x[1]*x[1] );

  PosType Qm=( pow(fp*sqrt(b2+bc*bc)-x[0],2) + f*f*f*f*x[1]*x[1] )
    /( pow(f*r*r-fp*bc*x[0],2)+fp*fp*bc*bc*x[1]*x[1] );

  //printf("Q = %e %e\n",Qp,Qm);

  //PosType QpQm = ( pow(f*r*r-fp*bc*x[0],2)+fp*fp*bc*bc*x[1]*x[1] )/( pow(f*r*r+fp*bc*x[0],2)+fp*fp*bc*bc*x[1]*x[1] );

  //printf(" fp=%e f=%e \n",fp,f);
  //printf("Qp=%e Qm=%e\n",Qp,Qm);
  //printf("Qp/Qm=%e log(Qp/Qm)=%e %e\n",Qp/Qm,log((float)(Qp/Qm)),log(6.853450e-02));
  angle[0]=0.25*sqrt(f)*log(Qp/Qm)/fp;

  //printf("angle[0]=%e Qp=%e Qm=%e\n",angle[0],Qp,Qm);
  //exit(0);
  RCphase=atan2(-2*f*f*fp*sqrt(b2+bc*bc)*x[1],x[0]*x[0]+pow(f*f*x[1],2)-fp*fp*(b2+bc*bc));
  SCphase=atan2(-2*f*fp*bc*x[1],f*f*r*r-fp*fp*bc*bc);

  angle[1]= -0.5*sqrt(f)*(RCphase-SCphase)/fp;

  Utilities::rotation(alpha,angle,-theta);

  alpha[0] = -1.0*alpha[0];
  alpha[1] = -1.0*alpha[1];

  if(alpha[0] != alpha[0] || alpha[1] != alpha[1] ){
	  printf("alpha is %e %e in nsie.c \n fp=%e b2=%e r=%e bc=%e f=%e theta=%e\n x = %e %e xt= %e %e\n"
			  ,alpha[0],alpha[1],fp,b2,r,bc,f,theta,x[0],x[1],xt[0],xt[1]);
	  printf("angle=%e %e Qp=%e Qm=%e RCphase=%e SCphase=%e\n",angle[0],angle[1],Qp,Qm,RCphase,SCphase);
	  ERROR_MESSAGE();
	  exit(0);
	  alpha[0]=alpha[1]=0;
  }
  return;
}

/**\ingroup DeflectionL2 \ingroup function
 * \brief Convergence for non-singular isothermal ellipsoid, units \f$ \frac{r_{einstein}}{units(x)} \f$
 * or \f$ \frac{\sigma^2}{\Sigma_{crit}G\, units(xt) } \f$
 */
KappaType kappaNSIE(
		PosType *xt     /// position on the image plane in Einstein radius units
		,PosType f      /// axis ratio of mass
		,PosType bc     /// core size in units of Einstein radius
		,PosType theta  /// position angle of ellipsoid
		){
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
		,PosType *xt     /// position on the image plane in Einstein radius units
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
		PosType *x
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
	void rotation(PosType *xout,PosType *xin,PosType theta){
    
    PosType tmp = xin[0];  // to make it work if xout == xin
		xout[0]=tmp*cos(theta)-xin[1]*sin(theta);
		xout[1]=xin[1]*cos(theta)+tmp*sin(theta);
	}
}
/* potential in Mpc^2 */
KappaType phiNSIE(PosType *xt,PosType f,PosType bc,PosType theta){

	return 0.0;
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
