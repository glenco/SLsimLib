/*
 * expanded_powerlaw.c
 *
 *  Created on: Feb 22, 2010
 *      Author: R.B. Metcalf
 *

  calculate the sources position, surface density and magnification at x
        given lens model mod
        mod[1] first component of shear
        mod[2] second component of shear
        mod[3...Nmodes-1]
        beta - radial power-law index
        returns kappa

        came from previous program lensexpand1.1.c
***************************************************************/
#include "analytic_lens.h"

double deflection_model(double beta,double *mod,int Nmodes,double *x,double *y,double *mag){
  double F,F1,F2,theta,r,cosx,sinx,dxdr,dxda,dydr,dyda;
  int i,k;

  theta=atan2(x[1],x[0]);
  r=sqrt(x[0]*x[0] + x[1]*x[1]);
  cosx=x[0]/r;
  sinx=x[1]/r;

  F=0.5*mod[3];
  F1=0;
  F2=0;
  for(i=4;i<Nmodes;i+=2){
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

  // magnification matrix in polar coordinates
  dxdr= (1+mod[1])*cosx + mod[2]*sinx
    - (beta-1)*pow(r,beta-2)*(beta*cosx*F-sinx*F1);
  dxda=-(1+mod[1])*r*sinx + mod[2]*r*cosx
    + pow(r,beta-1)*(beta*sinx*F + (1-beta)*cosx*F1 + sinx*F2);
  dydr= (1-mod[1])*sinx + mod[2]*cosx
    - (beta-1)*pow(r,beta-2)*(beta*sinx*F+cosx*F1);
  dyda= (1-mod[1])*r*cosx - mod[2]*r*sinx
    + pow(r,beta-1)*(-beta*cosx*F + (1-beta)*sinx*F1 - cosx*F2);

  // actually the inverse magnification
  // *mag=(dxdr*dyda - dxda*dydr)/r;

  // convert magnification matrix to Cartesian coordinates
  mag[0]=( -x[1]*dxda + r*x[0]*dxdr )/r/r;  // xx
  mag[1]=(  x[0]*dyda + r*x[1]*dydr )/r/r;  // yy
  mag[2]=(  x[0]*dxda + r*x[1]*dxdr )/r/r;  // xy

  return 0.5*(beta*beta*F+F2)*pow(r,beta-2);
}

