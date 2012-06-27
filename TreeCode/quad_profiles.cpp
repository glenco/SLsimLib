/*
 * profiles.cpp
 *
 *  Created on: Mar 6, 2012
 *      Author: bmetcalf
 */
#include <slsimlib.h>

/// Gaussian particles
double QuadTree::alpha_o(double r2,float sigma){
  if(r2==0) return 0.0;
  if(sigma <= 0.0 ) return -1.0/r2/pi;
  return ( exp(-r2/sigma/sigma/2) - 1.0 )/r2/pi;
}
double QuadTree::kappa_o(double r2,float sigma){
  if(sigma <= 0.0) return 0.0;
  return exp(-r2/sigma/sigma/2)/2/pi/sigma/sigma;
}
double QuadTree::gamma_o(double r2,float sigma){
  if(r2==0) return 0.0;
  if(sigma <= 0.0) return -2.0/pi/pow(r2,2);
  return (-2.0 + (2.0 + r2/sigma/sigma)*exp(-r2/sigma/sigma/2) )/pi/pow(r2,2);
}
double QuadTree::phi_o(double r2,float sigma){
	ERROR_MESSAGE();  // not yet written
	exit(1);
	return 0;
}

/** power-law profiles
 *  All quantities should be divided by Sigma_crit to get the usual
 */

double QuadTreePowerLaw::alpha_h(double r2,HaloStructure &par){

	if(r2==0) return 0.0;

	double b=0,r = sqrt(r2);

	b = -par.mass/r2/pi;
	if(r < par.Rmax) b *= pow(r/par.Rmax,beta+2);

	return b;
}
double QuadTreePowerLaw::kappa_h(double r2,HaloStructure &par){

	if(r2 > par.Rmax*par.Rmax) return 0.0;
	if(r2 < 1.0-20) r2=1.0e-20;
	return (beta+2)*par.mass*pow(r2/par.Rmax/par.Rmax,beta/2)/(pi*pow(par.Rmax,2));
}
double QuadTreePowerLaw::gamma_h(double r2,HaloStructure &par){
	double gt=0;

	if( r2 <= 0.0) return 0.0;
	gt = -2.0*par.mass/pi/pow(r2,2);
	if(r2 < par.Rmax*par.Rmax) gt *= -beta*pow(r2/par.Rmax/par.Rmax,beta/2+1)/2;

	return gt;
}
double QuadTreePowerLaw::phi_h(double r2,HaloStructure &par){
	double b;

	b=par.mass/pi;
	double r = sqrt(r2);
	if(r <= par.Rmax) return b*pow(r/par.Rmax,beta+2);
	return b*(log(r/par.Rmax) + 1);
}

// NFW profile

// Gives the magnitude of the deflection angle modulo Sigma_crit
double QuadTreeNFW::alpha_h(double r2,HaloStructure &par){
	double b=0;

	if(r2 <= 0) return 0.0;

	b = -par.mass/r2/pi;

	double r = sqrt(r2);
	if(r < par.Rmax){
		double y;

		y = par.Rmax/par.rscale;
		b = rhos(y);
		y = r/par.rscale;
		b*= gfunction(y)/par.rscale;
	}

	return b;
}
/// Convergence for an NFW halo
double QuadTreeNFW::kappa_h(double r2,HaloStructure &par){
	double r = sqrt(r2);

	if(r >= par.Rmax) return 0.0;
	if(r < 1.0-20) r=1.0e-20;

	double y,b;

	b = par.mass/2/pi/par.rscale;
	y = par.Rmax/par.rscale;
	b *= rhos(y);
	y = r/par.rscale;
	b*= ffunction(y);

	return b;
}

/// Shear for and NFW halo. this might have a flaw in it
double QuadTreeNFW::gamma_h(double r2,HaloStructure &par){
	double gt=0;

	if(r2 <= 0.0) return 0.0;

	double r = sqrt(r2);

	if(r > par.Rmax)
		gt = -2.0*par.mass/pi/pow(r2,2);
	else{
		double y;
		gt = par.mass/2/pi/par.rscale;
		y = par.Rmax/par.rscale;
		gt *= rhos(y);
		y = r/par.rscale;
		gt *= g2function(y);
	}

	return gt;
}
double QuadTreeNFW::phi_h(double r2,HaloStructure &par){
	ERROR_MESSAGE();
	std::cout << "time delay has not been fixed for NFW profile yet." << std::endl;
	exit(1);
	return 0.0;
}

double QuadTreeNFW::gfunction(double x){
	double ans;

	ans=log(x/2);
	if(x==1.0){ ans += 1.0; return ans;}
	if(x>1.0){  ans +=  2*atan(sqrt((x-1)/(x+1)))/sqrt(x*x-1);; return ans;}
	if(x<1.0){  ans += 2*atanh(sqrt((1-x)/(x+1)))/sqrt(1-x*x);; return ans;}
	return 0.0;
}
double QuadTreeNFW::ffunction(double x){
	double ans;

	if(x==1.0){ return 1.0;}
	if(x>1.0){  ans = (1-2*atan(sqrt((x-1)/(x+1)))/sqrt(x*x-1))/(x*x-1); return ans;}
	if(x<1.0){  ans = (1-2*atanh(sqrt((1-x)/(x+1)))/sqrt(1-x*x))/(x*x-1); return ans;}
	return 0.0;
}

double QuadTreeNFW::g2function(double x){
	double ans;

	if(x==1) return 10/3. + 4*log(0.5);

	ans=4*log(x/2.0)/x/x - 2/(x*x-1);
	if(x<1.0){ ans+= 4*(2/x/x+1/(x*x-1))*atanh(sqrt((1-x)/(x+1)))/sqrt(1-x*x); return ans;}
	if(x>1.0){ ans+= 4*(2/x/x+1/(x*x-1))*atan(sqrt((x-1)/(x+1)))/sqrt(x*x-1); return ans;}

	return 0.0;
}

double QuadTreeNFW::rhos(double x){
	double ans;

	ans = 1.0;
	ans /= log(1+x) - x/(1+x);

	return ans;
}
/* Profiles with a flat interior and a power-law drop after rscale
 *  All quantities should be divided by Sigma_crit to get the usual lens plane units.
 */

double QuadTreePseudoNFW::alpha_h(double r2,HaloStructure &par){

	if(r2<=0.0) return 0.0;

	double b=0,r = sqrt(r2);

	b = -par.mass/r2/pi;
	if(r < par.Rmax) b *= mhat(r/par.rscale)/mhat(par.Rmax/par.rscale);

	return b;
}
double QuadTreePseudoNFW::kappa_h(double r2,HaloStructure &par){

	if(r2 > par.Rmax*par.Rmax) return 0.0;
	if(r2 < 1.0-20) r2=1.0e-20;
	return par.mass/(2*pi*pow(par.rscale,2)*mhat(par.Rmax/par.rscale))/pow(1+sqrt(r2)/par.rscale,beta);
}

double QuadTreePseudoNFW::gamma_h(double r2,HaloStructure &par){
	double gt=0;

	if( r2 <= 0.0) return 0.0;
	gt = -2.0*par.mass/pi/r2/r2;
	if(r2 < par.Rmax*par.Rmax){
		double y = sqrt(r2)/par.rscale;
		//gt *= (1+2*y)*r2/pow(par.rscale,2)/mhat(par.Rmax/par.rscale)/pow(1+y,beta);
		gt *= -(y*y/pow(1+y,beta) - 2*mhat(y))/mhat(par.Rmax/par.rscale)/2;
	}

	return gt;
}

double QuadTreePseudoNFW::phi_h(double r2,HaloStructure &par){

	ERROR_MESSAGE();
	std::cout << "time delay has not been fixed for NFW profile yet." << std::endl;
	exit(1);
	return 0.0;
}

double QuadTreePseudoNFW::mhat(double y){
	switch(beta){
	case 1:
		return y - log(1+y);
		break;
	case 2:
		return log(1+y) - y/(1+y);
		break;
	default:
		return - y/(beta-2)/pow(1+y,beta-1) - 1/(beta-2)/(beta-1)/pow(1+y,beta-1);
		break;
	}
}

/*
 *
double QuadTreePseudoNFW::alpha_h(double r2,HaloStructure &par){

	if(r2<=0.0) return 0.0;

	double b=0,r = sqrt(r2);

	b = -par.mass/r2/pi;
	if(r < par.Rmax) b *= mhat(r/par.rscale)/mhat(par.Rmax/par.rscale);

	return b;
}
double QuadTreePseudoNFW::kappa_h(double r2,HaloStructure &par){

	if(r2 > par.Rmax*par.Rmax) return 0.0;
	if(r2 < 1.0-20) r2=1.0e-20;
	return par.mass/(2*pi*pow(par.rscale,2)*mhat(par.Rmax/par.rscale))/pow(1+sqrt(r2)/par.rscale,beta);
}
double QuadTreePseudoNFW::gamma_h(double r2,HaloStructure &par){
	double gt=0;

	if( r2 <= 0.0) return 0.0;
	gt = -par.mass/pi/pow(r2,2);
	if(r2 < par.Rmax*par.Rmax){
		double y = sqrt(r2)/par.rscale;
		gt *= (mhat(y) - 0.5*y*y/pow(1+y,beta))/mhat(sqrt(r2)/par.Rmax);
	}

	return gt;
}
double QuadTreePseudoNFW::phi_h(double r2,HaloStructure &par){

	ERROR_MESSAGE();
	std::cout << "time delay has not been fixed for NFW profile yet." << std::endl;
	exit(1);
	return 0.0;
}

double QuadTreePseudoNFW::mhat(double y){
	return (1 + ( 1 + (beta-1)*y)/pow(y+1,beta-1)  )/(beta-2)/(beta-1);
}
*/
