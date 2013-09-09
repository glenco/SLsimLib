/** 
 *  Created on: July 30, 2013
 *      Author: D. Leier
 */

#include "lens_halos.h"

double LensHaloJaffe::gfunction(double x){
	double ans;
	ans=pi;
	if(x==0) x=1e-5;
	if(x==1.0) return (ans-2.)*x;
	if(x<1.0){ ans -= 2.*x*acosh(1./x)/sqrt(1.-x*x) ; return ans*x;}
	if(x>1.0){ ans -= 2.*x*acos(1./x)/sqrt(x*x-1.)  ; return ans*x;}
	return 0.0;
}

double LensHaloJaffe::ffunction(double x){
	double ans;
	if(x==0) x=1e-5;
	ans=pi/x;
	if(x==1.0){ return ans-(2.+2./3.);}
	if(x>1.0){  ans += 2./(1.-x*x)*(1.-(2.-x*x)*acos(1./x)/sqrt(x*x-1.)); return ans;}
	if(x<1.0){  ans += 2./(1.-x*x)*(1.-(2.-x*x)*acosh(1./x)/sqrt(1.-x*x)); return ans;}
	return 0.0;
}


double LensHaloJaffe::g2function(double x){
	double ans;
	if(x==0) x=1e-5;
	ans=pi/x;
	if(x==1.0){ ans -=4./3. ; return ans*x/3.;}
	if(x>1.0){  ans += 2./x/x*(-1.0*(2.*x*x*acos(1./x))/sqrt(x*x-1.)+(sqrt(1.-1./x)*sqrt(1.+1./x)*x*(log(x-1.0)+log(1.+x)))/sqrt(x*x-1.)-log(x*x-1.))-(2./(1.-x*x)*(1.-(2.-x*x)*acos(1./x)/sqrt(x*x-1.))); return ans*x/3.;}
	if(x<1.0){  ans += 2./x/x*(-(2.*x*x*acosh(1./x))/sqrt(1.-x*x)+(sqrt(-1.+1./x)*sqrt(1.+1./x)*x*(log(1.-x)+log(1.+x)))/sqrt(1.-x*x)-log(1.-x*x))-(2./(1.-x*x)*(1.-(2.-x*x)*acosh(1./x)/sqrt(1.-x*x))); return ans*x/3.;}
	return 0.0;
}

double LensHaloJaffe::hfunction(double x){
	ERROR_MESSAGE();
	std::cout << "no analytic expression" << std::endl;
	exit(1);
}
