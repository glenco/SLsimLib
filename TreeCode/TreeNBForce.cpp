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
