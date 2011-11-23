/*  routines for calculating the lensing properties of  */
/*   a non-singular isothermal ellipsoid */
/*  written by R.B. Metcalf, March 18, 2009 */
/*  based on the analytic solutions of Kormann et al. 1993
 *  convention here is gamma_1 = -(Axx-Ayy)/2 and gamma_2= -Axy
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "analytic_lens.h"

void alphaNSIE(double *alpha,double *xt,double f,double bc,double theta){
  double x[2],angle[2],fp,b2,r,Qp,Qm,RCphase,SCphase;
  void rotation(double *xout,double *xin,double theta);

  r=sqrt(xt[0]*xt[0]+xt[1]*xt[1]);

  if( r < 1.0e-20 || r > 1.0e20){
	  alpha[0]=0.0;
	  alpha[1]=0.0;
	  return;
  }

  if( f==1.0 ){
    if(bc == 0.0){
      alpha[0]=xt[0]/r;
      alpha[1]=xt[1]/r;
    }else{
      alpha[0]=(sqrt(r*r+bc*bc) - bc)*xt[0]/r/r;
      alpha[1]=(sqrt(r*r+bc*bc) - bc)*xt[1]/r/r;
    }
    return;
  }

  rotation(x,xt,theta);

  fp=sqrt(1-f*f);
  b2=x[0]*x[0]+f*f*x[1]*x[1];
  r=sqrt(x[0]*x[0]+x[1]*x[1]);

  Qp=( pow(fp*sqrt(b2+bc*bc)+x[0],2) + f*f*f*f*x[1]*x[1] ) 
    /( pow(f*r*r+fp*bc*x[0],2)+fp*fp*bc*bc*x[1]*x[1] );
  Qm=( pow(fp*sqrt(b2+bc*bc)-x[0],2) + f*f*f*f*x[1]*x[1] ) 
    /( pow(f*r*r-fp*bc*x[0],2)+fp*fp*bc*bc*x[1]*x[1] );

  /*  printf("Q = %e %e\n",Qp,Qm);*/

  //printf(" fp=%e f=%e \n",fp,f);
  //printf("Qp=%e Qm=%e\n",Qp,Qm);
  //printf("Qp/Qm=%e log(Qp/Qm)=%e %e\n",Qp/Qm,log((float)(Qp/Qm)),log(6.853450e-02));
  angle[0]=0.25*sqrt(f)*log(Qp/Qm)/fp;

  //printf("angle[0]=%e Qp=%e Qm=%e\n",angle[0],Qp,Qm);
//exit(0);
  RCphase=atan2(-2*f*f*fp*sqrt(b2+bc*bc)*x[1],x[0]*x[0]+pow(f*f*x[1],2)-fp*fp*(b2+bc*bc));
  SCphase=atan2(-2*f*fp*bc*x[1],f*f*r*r-fp*fp*bc*bc);

  angle[1]= -0.5*sqrt(f)*(RCphase-SCphase)/fp;

  rotation(alpha,angle,-theta);

  if(isnan(alpha[0]) || isnan(alpha[1]) ){
	  printf("alpha is %e %e in nsie.c \n fp=%e b2=%e r=%e bc=%e f=%e theta=%e\n x = %e %e xt= %e %e\n"
			  ,alpha[0],alpha[1],fp,b2,r,bc,f,theta,x[0],x[1],xt[0],xt[1]);
	  printf("angle=%e %e Qp=%e Qm=%e RCphase=%e SCphase=%e\n",angle[0],angle[1],Qp,Qm
			  ,RCphase,SCphase);
	  exit(0);
	  alpha[0]=alpha[1]=0;
  }
}

/* surface density */
double kappaNSIE(double *xt,double f,double bc,double theta){
  double x[2],b2;
  void rotation(double *xout,double *xin,double theta);

  rotation(x,xt,theta);

  b2=x[0]*x[0]+f*f*x[1]*x[1];
  if( (b2 < 1.0e-20)*(bc < 1.0e-20)){
	  return 1.0e10;
  }
  if(b2>1.0e20 ) return 0.0;
  return 0.5*sqrt(f/(b2+bc*bc));
}
     /* shear */
void gammaNSIE(double gam[2],double *xt,double f,double bc,double theta){
  double x[2],fp,P,b2,r;
  void rotation(double *xout,double *xin,double theta);

  r=sqrt(xt[0]*xt[0]+xt[1]*xt[1]);

  if(r < 1.0e-20 || r > 1.0e20){
	  gam[0]=gam[1]=0.0;
 	  return;
  }

  rotation(x,xt,theta);

  fp=sqrt(1-f*f);

  b2=x[0]*x[0]+f*f*x[1]*x[1];
  //r=sqrt(x[0]*x[0]+x[1]*x[1]);

  P=sqrt(f)*( kappaNSIE(x,f,bc,0.0)*(x[0]*x[0]+f*f*f*f*x[1]*x[1])/sqrt(f) 
              - 0.5*(1+f*f)*sqrt(b2+bc*bc)+f*bc ) 
    /( pow(f*r,4) -2*f*f*fp*fp*bc*bc*(x[0]*x[0]-x[1]*x[1])+pow(fp*bc,4) );

  gam[0]=(f*f*(x[0]*x[0]-x[1]*x[1])-fp*fp*bc*bc)*P;
  gam[1]=2*f*f*x[0]*x[1]*P;

  rotation(gam,gam,-2*theta);
  return;
}

     /* inverse magnification */
double invmagNSIE(double *x,double f,double bc,double theta
		     ,double *gam,double kap){

  gammaNSIE(gam,x,f,bc,theta);
  kap=kappaNSIE(x,f,bc,theta);
  return pow(1-kap,2) - gam[0]*gam[0] - gam[1]*gam[1];
}

void rotation(double *xout,double *xin,double theta){

  xout[0]=xin[0]*cos(theta)-xin[1]*sin(theta);
  xout[1]=xin[1]*cos(theta)+xin[0]*sin(theta);
}

/* potential in Mpc^2 */
double phiNSIE(double *xt,double f,double bc,double theta){

	return 0.0;
}
