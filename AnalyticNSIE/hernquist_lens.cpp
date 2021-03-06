/** 
 *  Created on: June 5, 2013
 *      Author: D. Leier
 */

#include "image_processing.h"
#include "lens_halos.h"


PosType LensHaloHernquist::gfunction(PosType x) const{
	PosType ans;
	if(x==0) x=1e-5;
	if(x==1.0) return 2./3. ;//2./3.;
	if(x<1.0){ ans =  2.*x*(1.-(atanh(sqrt(1.-x*x))/sqrt(1.-x*x)))/(x*x-1.) ;; return ans*x;} // ans
    //if(x<1.0){ ans =  2.*(1.+1./(x*x-1.)+sqrt(1.-x*x)*x*x*(0.5*(log((1.+sqrt(1.-x*x))/(1.-sqrt(1.-x*x)))))/((x*x-1.)*(x*x-1.))) ;; return ans/x;}
	if(x>1.0){ ans =  2.*x*(1.-(atan(sqrt(x*x-1.))/sqrt(x*x-1.)))/(x*x-1.) ;; return ans*x;}
	return 0.0;
}

PosType LensHaloHernquist::ffunction(PosType x) const{
	PosType ans;

	if(x==0) x=1e-5;
	if(x==1.0){ return (0.6-(1./3.));}
	if(x>1.0){  ans = (-3.+(2.+x*x)*(atan(sqrt(x*x-1.))/sqrt(x*x-1.)))/(x*x-1.)/(x*x-1.); return ans;}
	if(x<1.0){  ans = (-3+(2+x*x)*(atanh(sqrt(1.-x*x))/(sqrt(1.-x*x))))/(x*x-1.)/(x*x-1.); return ans;}
	return 0.0;
}

PosType LensHaloHernquist::g2function(PosType x) const{
	PosType ans,ax;
	if(x==0) x=1e-5;
	if(x==1) return 2./3./(x*x) - ffunction(x);
	ax=sqrt((1.0-x)*(x+1.0));
	if(x<1.0){ ans= 2*(1+1/(x*x-1)+ax*x*x*(0.5*(log(1+ax)-log(1-ax)))/((x*x-1)*(x*x-1)))/(x*x) ; return ans-ffunction(x);}
	ax=sqrt((x-1.0)*(x+1.0));
	if(x>1.0){ ans= 2*(1+1/(x*x-1)-x*x*atan(ax)/pow(x*x-1,3./2.))/(x*x); return ans-ffunction(x);}
	return 0.0;
}

PosType LensHaloHernquist::hfunction(PosType x) const{
	PosType ans,ax;
	if(x==0) x=1e-5;
	ans=log(0.25*x*x);
	ax=sqrt(x*x-1.);
	if(x>=1.0){ ans+= 2*atan(ax)/ax ; return ans;}
	ax=sqrt(1.-x*x);
	if(x<1.0){ ans+= 2*atanh(ax)/ax; return ans;}
	return 0.0;
}

