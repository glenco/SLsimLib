/** 
 *  Created on: July 30, 2013
 *      Author: D. Leier
 */

#include "lens_halos.h"

PosType LensHaloJaffe::gfunction(PosType x) const{
	PosType ans;
	ans=pi;
	if(x==0) x=1e-5;
	if(x==1.0) return (ans-2.)*x;
	if(x<1.0){ ans -= 2.*x*acosh(1./x)/sqrt(1.-x*x) ; return ans*x;}
	if(x>1.0){ ans -= 2.*x*acos(1./x)/sqrt(x*x-1.)  ; return ans*x;}
	return 0.0;
}

PosType LensHaloJaffe::ffunction(PosType x) const{
	PosType ans;
	if(x==0) x=1e-5;
	ans=pi/x;
	if(x==1.0){ return ans-(2.+2./3.);}
	if(x>1.0){  ans += 2./(1.-x*x)*(1.-(2.-x*x)*acos(1./x)/sqrt(x*x-1.)); return ans;}
	if(x<1.0){  ans += 2./(1.-x*x)*(1.-(2.-x*x)*acosh(1./x)/sqrt(1.-x*x)); return ans;}
	return 0.0;
}


PosType LensHaloJaffe::g2function(PosType x) const{
	PosType ans;
	if(x==0) x=1e-5;
	ans=pi/x;
	if(x==1.0){ ans -=4./3. ; return ans*x/3.;}
	if(x>1.0){  ans += 2./x/x*(-1.0*(2.*x*x*acos(1./x))/sqrt(x*x-1.)+(sqrt(1.-1./x)*sqrt(1.+1./x)*x*(log(x-1.0)+log(1.+x)))/sqrt(x*x-1.)-log(x*x-1.))-(2./(1.-x*x)*(1.-(2.-x*x)*acos(1./x)/sqrt(x*x-1.))); return ans*x/3.;}
	if(x<1.0){  ans += 2./x/x*(-(2.*x*x*acosh(1./x))/sqrt(1.-x*x)+(sqrt(-1.+1./x)*sqrt(1.+1./x)*x*(log(1.-x)+log(1.+x)))/sqrt(1.-x*x)-log(1.-x*x))-(2./(1.-x*x)*(1.-(2.-x*x)*acosh(1./x)/sqrt(1.-x*x))); return ans*x/3.;}
	return 0.0;
}

PosType LensHaloJaffe::hfunction(PosType x) const{
	ERROR_MESSAGE();
	std::cout << "there is yet no analytic expression" << std::endl;
	exit(1);
}


PosType LensHaloJaffe::bfunction(PosType fx){
    PosType ans;
    PosType fac=-1.0;
    PosType x=fx;
    if(x==0) x=1e-5;
    if(x==1){return -2.1231*fac;} // (x/ffunction(x))*(2.+2./15.)
    PosType aux=sqrt(1.-x*x);
    if(x<1){ans=x*((-(pi/x/x)+(2.*(((2.-x*x))/(sqrt(-1+1./x)*sqrt(1+1./x) *x*x* aux)+(2*x*acosh(1./x))/aux-(x*(2.-x*x)*acosh(1./x))/pow(aux,3)))/(1.-x*x)+(4.*x*(1.-((2.-x*x)*acosh(1./x))/aux))/(1.-x*x)/(1.-x*x))/(pi/x+(2.*(1.-((2.-x*x)*acosh(1./x))/aux))/(1.-x*x))); return fac*ans;}
    if(x>1){ans=(x*(-1.0*(pi/x/x)+(2.*(-(((2.-x*x))/(sqrt(1-1./x/x)*x*x*sqrt(-1.+x*x)))+(x*(2.-x*x)*acos(1./x))/(pow(-1.+x*x,3./2.))+(2.*x*acos(1./x))/sqrt(-1.+x*x)))/(1.-x*x)+(4.*x*(1.-((2.-x*x)*acos(1./x))/sqrt(-1.+x*x)))/(1.-x*x)/(1.-x*x)))/(pi/x+(2.*(1.-((2.-x*x)*acos(1./x))/sqrt(-1.+x*x)))/(1.-x*x)); return fac*ans;}
//if(x>1.0){  ans =(x/ffunction(x))*((4.-6.*x*x+2.*pow(x,4)+2.*aux/x*pow(x,5)*acos(1./x))/(aux*x*(1.-x*x)*pow(aux,3))+(4.*x*(1.+((-2.+x*x)*acos(1/x))/aux))/(1.-x*x)/(1.-x*x)) ; return fac*ans;}
//	if(x<1.0){  ans = (x/ffunction(x))*(((4.-6.*x*x+2.*pow(x,4)-2.*sqrt(-1+1./x)*pow(x,5)*sqrt((1.+x)/x)*acosh(1./x))/(sqrt(-1+1./x)*x*x*sqrt((1.+x)/x)*sqrt(1.-x*x))+4.*x*(1.+((-2.+x*x)*acosh(1./x))/sqrt(1. -x*x)))/(1. -x*x)/(1. -x*x))
//        ; return fac*ans;}
    return 0.0;
}

PosType LensHaloJaffe::dbfunction(PosType x){
    PosType h=1e-5;
    if(x==0) x=1e-5;
    //if(x>0){return ((1./12.)*log(bfunction(x-2*h))-(2./3.)*log(bfunction(x-h))+(2./3.)*log(bfunction(x+h))-(1./12.)*log(bfunction(x+2*h)))/h;};
    if(x>0){return ((1./12.)*(bfunction(x-2*h))-(2./3.)*(bfunction(x-h))+(2./3.)*(bfunction(x+h))-(1./12.)*(bfunction(x+2*h)))/h;};
    return 0.0;
}

PosType LensHaloJaffe::ddbfunction(PosType x){
    PosType h=1e-5;
    if(x==0) x=1e-5;
    if(x>0){return (bfunction(x-h)-2.*bfunction(x)+bfunction(x+h))/(h*h);};
    //if(x>0){return ((-1./12.)*bfunction(x-2*h)+(4./3.)*bfunction(x-h)-(5./2.)*bfunction(x)+(4./3.)*bfunction(x+h)-(1./12.)*bfunction(x+2*h))/(h*h);};
    //if(x>0){return ((1./90.)*bfunction(x-3*h)+(-3./20.)*bfunction(x-2*h)+(3./2.)*bfunction(x-h)-(49./18.)*bfunction(x)+(3./2.)*bfunction(x+h)-(3./20.)*bfunction(x+2*h)+(1./90.)*bfunction(x+3*h))/(h*h);};
    return 0.0;
}




