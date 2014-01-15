/** \ingroup FitLensL2
 * lens_expand.c
 *
 *  Created on: Feb 22, 2010
 *      Author: R.B. Metcalf
 *

\brief  calculate the sources position, surface density and magnification at x
        given lens model mod
        mod[0] uniform kappa sheet
        mod[1] first component of background shear
        mod[2] second component of background shear
        mod[3] monopole mode  ( = 2*sigma^2/G/Sigma_crit for isothermal )
        mod[4-5] dipolepole n=2
        mod[6-7] quadrapole n=3
        mod[8-9]  octopole  n=4
        mode[10 ... ] higher order modes of frequency mode/2
        beta - radial power-law index, Note kappa propto r^-beta !!!!
        returns kappa


 kappa = 0.5 sum_n (beta^2 - n^2) [ a_n cos(n theta) + b_n sin(n theta) ] r^(beta-2)

 *  convention here is gamma_1 = -(Axx-Ayy)/2 and gamma_2= -Axy

        came from previous program lensexpand1.1.c
***************************************************************/

#include "lens_halos.h"

/// removed declaration of Grav -- caused compiler problems and was not needed

PosType lens_expand(PosType beta,PosType *mod,int Nmodes,PosType *x,PosType *alpha,KappaType *gamma,KappaType *phi){
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

