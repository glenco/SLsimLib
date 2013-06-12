/** 
 *  Created on: June 5, 2013
 *      Author: D. Leier
 */

#include "lens_halos.h"

/// deflection caused by NFW halo
void alphaHern(double *alpha,double *x,double Rtrunc,double mass,double r_scale
		,double *center,double Sigma_crit){
	double r,b=0;

	r=sqrt(pow(x[0]-center[0],2) + pow(x[1]-center[1],2));
	if(r==0){
		alpha[0]=alpha[1]=0.0;
		return ;
	}
	b=mass/pow(r,2)/pi/Sigma_crit;

	if(r < Rtrunc){
		double y;
		double ghernfunction(double);

		y = Rtrunc/r_scale;
		b/= ghernfunction(y);
		y = r/r_scale;
		b*= ghernfunction(y);
	}

	alpha[0]=b*(x[0]-center[0]);
	alpha[1]=b*(x[1]-center[1]);

	return ;
}
/// Convergence for an NFW halo
KappaType kappaHern(double *x,double Rtrunc,double mass,double r_scale
		,double *center,double Sigma_crit){
	double r;
	double ghernfunction(double),fhernfunction(double);

	r=sqrt(pow(x[0]-center[0],2) + pow(x[1]-center[1],2));
	if(r>=Rtrunc) return 0.0;
	if(r < 1.0e-20) r=1.0e-20;

	double y,b;

	b=1.0;
	y = Rtrunc/r_scale;
	b/= ghernfunction(y);
	y = r/r_scale;
	b*= fhernfunction(y);

	return b*mass/(pi*pow(r_scale,2)*Sigma_crit);
}

/// Shear for and NFW halo. this might have a flaw in it
void gammaHern(KappaType *gamma,double *x,double Rtrunc,double mass,double r_scale
		,double *center,double Sigma_crit){
	double r,gt=0;
	double g2hernfunction(double x);

	r=sqrt(pow(x[0]-center[0],2) + pow(x[1]-center[1],2));
	if(r==0.0){
		gamma[0]=gamma[1]=0.0;
		return ;
	}

	gt=mass/pi/Sigma_crit/pow(r_scale,2);
	if(r<Rtrunc){
		double y;
		double ghernfunction(double);

		y = Rtrunc/r_scale;
		gt /= ghernfunction(y);
		y = r/r_scale;
		gt *= g2hernfunction(y);
		//gt -=kappaHern(&y,Rtrunc,mass,r_scale,center,Sigma_crit);
	}

	gamma[0]=-gt*(pow(x[0]-center[0],2)-pow(x[1]-center[1],2))/r/r;
	gamma[1]=-2*gt*(x[0]-center[0])*(x[1]-center[1])/r/r;

	return ;
}

double ghernfunction(double x){
	double ans;

	if(x==0) x=1e-5;
	if(x<=1.0){ ans =  (x*x)*(((log((1.+sqrt(1.-x*x))/x))/sqrt(1.-x*x))-1.)/(1.-x*x) ;; return ans;}
	if(x>1.0){  ans =  (x*x)*(((acos(1./x))/sqrt(x*x-1.))-1.)/(1.-x*x) ;; return ans;}
	return 0.0;
}
double fhernfunction(double x){
	double ans;

	if(x==0) x=1e-5;
	if(x==1.0){ return (0.6-(1./3.));}
	if(x>1.0){  ans = ((2.+x*x)*(atan(sqrt(x*x-1.))/(sqrt(x*x-1.0)))-3.)/((x*x-1.)*(x*x-1.)); return ans;}
	if(x<1.0){  ans = ((2.+x*x)*(0.5*(log(1+sqrt(1.0-x*x))-log(1-sqrt(1.0-x*x)))/sqrt(1.0-x*x))-3.)/((x*x-1.)*(x*x-1.)); return ans;}
	return 0.0;
}

double g2hernfunction(double x){
	double ans,ax;

	if(x==0) x=1e-5;
	if(x==1) return 2./3./(x*x) - fhernfunction(x);
	ax=sqrt((1.0-x)*(x+1.0));
	if(x<1.0){ ans= 2*(1+1/(x*x-1)+ax*x*x*(0.5*(log(1+ax)-log(1-ax)))/((x*x-1)*(x*x-1)))/(x*x) ; return ans-fhernfunction(x);}
	ax=sqrt((x-1.0)*(x+1.0));
	if(x>1.0){ ans= 2*(1+1/(x*x-1)-x*x*atan(ax)/pow(x*x-1,3./2.))/(x*x); return ans-fhernfunction(x);}

	return 0.0;
}
