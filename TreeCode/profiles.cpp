/*
 * profiles.cpp
 *
 *  Created on: Mar 6, 2012
 *      Author: bmetcalf
 */
#include <slsimlib.h>

/// Gaussian particles
double alpha_o(double r2,float sigma){
  if(r2==0) return 0.0;
  if(sigma == 0.0 ) return -1.0/r2/pi;
  if(sigma < 0.0)  return -exp(-r2/sigma/sigma/2)/r2/pi;  // screened potential
  return ( exp(-r2/sigma/sigma/2) - 1.0 )/r2/pi;
}
double kappa_o(double r2,float sigma){
  if(sigma == 0.0) return 0.0;
  if(sigma < 0.0) return -exp(-r2/sigma/sigma/2)/2/pi/sigma/sigma;
  return exp(-r2/sigma/sigma/2)/2/pi/sigma/sigma;
}
double gamma_o(double r2,float sigma){
  if(r2==0) return 0.0;
  if(sigma == 0.0) return -2.0/pi/pow(r2,2);
  if(sigma < 0.0) return -(-2.0 + (2.0+(r2/sigma/sigma))*exp(-r2/sigma/sigma/2) )/pi/pow(r2,2);
  return (-2.0 + (2.0 + r2/sigma/sigma)*exp(-r2/sigma/sigma/2) )/pi/pow(r2,2);
}

/** power-law profiles
 *  All quantities should be divided by Sigma_crit to get the usual
 */

double alphaPowLaw(double r,double *alpha,HaloInternal &par){
	double b=0;

	if(r==0) return 0.0;

	b=par.mass/r/pi;
	if(r<par.Rmax) b *= pow(r/par.Rmax,par.beta+2);

	return b;
}
double kappaPowLaw(double r,HaloInternal &par){

	if(r>par.Rmax) return 0.0;
	if(r < 1.0-20) r=1.0e-20;
	return (par.beta+2)*par.mass*pow(r/par.Rmax,par.beta)/(2*pi*pow(par.Rmax,2));
}
double gammaPowLaw(double r,HaloInternal &par){
	double gt=0;

	if( r<= 0.0) return 0.0;
	gt=par.mass/pi/pow(r,2);
	if(r<par.Rmax) gt *= -par.beta*pow(r/par.Rmax,par.beta+2)/2;

	return gt;
}
double phiPowLaw(double r,HaloInternal &par){
	double b;

	b=par.mass/pi;
	if(r <= par.Rmax) return b*pow(r/par.Rmax,par.beta+2);
	return b*(log(r/par.Rmax) + 1);
}
