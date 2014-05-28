/** 
 *  Created on: Aug 26, 2010
 *      Author: R.B. Metcalf
 */

#include "lens_halos.h"

/// deflection caused by NFW halo
void LensHaloNFW::alphaNFW(PosType *alpha,PosType *x,PosType Rtrunc,PosType mass,PosType r_scale
		,PosType *center,PosType Sigma_crit){
	PosType r,b=0;

	r=sqrt(pow(x[0]-center[0],2) + pow(x[1]-center[1],2));
	if(r==0){
		alpha[0]=alpha[1]=0.0;
		return ;
	}
	b=mass/pow(r,2)/pi/Sigma_crit;

	if(r < Rtrunc){
		PosType y;
		//PosType gfunction(PosType);

		y = Rtrunc/r_scale;
		b/= gfunction(y);
		y = r/r_scale;
		b*= gfunction(y);
	}

	alpha[0]=b*(x[0]-center[0]);
	alpha[1]=b*(x[1]-center[1]);

	return ;
}
/// Convergence for an NFW halo
KappaType LensHaloNFW::kappaNFW(PosType *x,PosType Rtrunc,PosType mass,PosType r_scale
		,PosType *center,PosType Sigma_crit){
	PosType r;
	//PosType gfunction(PosType),ffunction(PosType);

	r=sqrt(pow(x[0]-center[0],2) + pow(x[1]-center[1],2));
	if(r>=Rtrunc) return 0.0;
	if(r < 1.0e-20) r=1.0e-20;

	PosType y,b;

	b=1.0;
	y = Rtrunc/r_scale;
	b/= gfunction(y);
	y = r/r_scale;
	b*= ffunction(y);

	return b*mass/(2*pi*pow(r_scale,2)*Sigma_crit);
}

/// Shear for and NFW halo. this might have a flaw in it
void LensHaloNFW::gammaNFW(KappaType *gamma,PosType *x,PosType Rtrunc,PosType mass,PosType r_scale
		,PosType *center,PosType Sigma_crit){
	PosType r,gt=0;
	//PosType g2function(PosType x);

	r=sqrt(pow(x[0]-center[0],2) + pow(x[1]-center[1],2));
	if(r==0.0){
		gamma[0]=gamma[1]=0.0;
		return ;
	}

	gt=mass/pi/Sigma_crit/pow(r,2);
	if(r<Rtrunc){
		PosType y;
		//PosType gfunction(PosType),ffunction(PosType);

		y = Rtrunc/r_scale;
		gt /= gfunction(y);
		y = r/r_scale;
		gt *= g2function(y)*pow(y,2);
	}

	gamma[0]=-gt*(pow(x[0]-center[0],2)-pow(x[1]-center[1],2))/r/r;
	gamma[1]=-2*gt*(x[0]-center[0])*(x[1]-center[1])/r/r;

	return ;
}

PosType LensHaloNFW::gfunction(PosType x){
	PosType ans;

	if(x==0) x=1e-5;
	ans=log(x/2);
	if(x==1.0){ ans += 1.0; return ans;}
	if(x>1.0){  ans +=  2*atan(sqrt((x-1)/(x+1)))/sqrt(x*x-1);; return ans;}
	if(x<1.0){  ans += 2*atanh(sqrt((1-x)/(x+1)))/sqrt(1-x*x);; return ans;}
	return 0.0;
}

PosType LensHaloNFW::ffunction(PosType x){
	PosType ans;

	if(x==0) x=1e-5;
	if(x==1.0){ return 1.0/3.0;}
	if(x>1.0){  ans = (1-2*atan(sqrt((x-1)/(x+1)))/sqrt(x*x-1))/(x*x-1); return ans;}
	if(x<1.0){  ans = (1-2*atanh(sqrt((1-x)/(x+1)))/sqrt(1-x*x))/(x*x-1); return ans;}
	return 0.0;
}

PosType LensHaloNFW::g2function(PosType x){
	PosType ans,y;

	if(x==0) x=1e-5;
	if(x==1) return 10/3. + 4*log(0.5);
	y=x*x-1;
	ans=4*log(x/2)/x/x -2/y;
	if(x<1.0){ ans+= 4*(2/x/x+1/y)*atanh(sqrt((1-x)/(x+1)))/sqrt(-y); return ans;}
	if(x>1.0){ ans+= 4*(2/x/x+1/y)*atan(sqrt((x-1)/(x+1)))/sqrt(y); return ans;}

	return 0.0;
}

PosType LensHaloNFW::hfunction(PosType x){
	if(x==0) x=1e-5;
	if(x<=1.0) return 2.*((log(x/2.)*log(x/2.))-atanh(sqrt(1.-x*x))*atanh(sqrt(1.-x*x)));
	if(x>1.0) return 2.* ((log(x/2.)*log(x/2.))+ atan(sqrt(x*x-1.))*atan(sqrt(x*x-1.)));
	return 0.0;
}

PosType LensHaloNFW::bfunction(PosType fx){
    PosType ans;
    PosType fac=1.0;
    PosType x=1.0*fx;
    if(x==0) x=1e-5;
    if(x==1){return fac*1.2;}
    PosType aux=sqrt(x*x-1);
    if(x>1.0){  ans = (x/ffunction(x)*((aux+2.*x*x*aux-6.0*x*x*atan(sqrt((x-1)/(x+1))))/(x*pow(aux,5)))); return fac*ans;}
	if(x<1.0){  ans = (x/ffunction(x))*
        ( sqrt(1./(1.+x))*
         ( (x-1.)*(1.+2.*x*x*sqrt(1.0/(1.0+x))*sqrt(1.+x))+6.0*x*x*sqrt(-1.0+2.0/(1.+x))*atanh( sqrt(-1.0+2.0/(1.0+x))) )
        )/(pow(-1.0+x,3)*x*pow(1.0+x,3./2.))
                       ; return fac*ans;}
    return 0.0;
}


/// dbfunction and ddbfunction are approximation formulae for dbeta/dr and d^2beta/dr^2 whereas dbnum and ddbnum calculate those derivates numerically using the 5-point rule, i.e. 4th order accuracy.

PosType LensHaloNFW::dbnum(PosType x){
    PosType h=1e-5;
    if(x==0) x=1e-5;
    //if(x>0){return ((1./12.)*log(bfunction(x-2*h))-(2./3.)*log(bfunction(x-h))+(2./3.)*log(bfunction(x+h))-(1./12.)*log(bfunction(x+2*h)))/h;};
    if(x>0){return ((1./12.)*(bfunction(x-2*h))-(2./3.)*(bfunction(x-h))+(2./3.)*(bfunction(x+h))-(1./12.)*(bfunction(x+2*h)))/h;};
    return 0.0;
}

PosType LensHaloNFW::ddbnum(PosType x){
    PosType h=1e-5;
    if(x==0) x=1e-5;
    if(x>0){return (bfunction(x-h)-2.*bfunction(x)+bfunction(x+h))/(h*h);};
    //if(x>0){return ((-1./12.)*bfunction(x-2*h)+(4./3.)*bfunction(x-h)-(5./2.)*bfunction(x)+(4./3.)*bfunction(x+h)-(1./12.)*bfunction(x+2*h))/(h*h);};
    //if(x>0){return ((1./90.)*bfunction(x-3*h)+(-3./20.)*bfunction(x-2*h)+(3./2.)*bfunction(x-h)-(49./18.)*bfunction(x)+(3./2.)*bfunction(x+h)-(3./20.)*bfunction(x+2*h)+(1./90.)*bfunction(x+3*h))/(h*h);};
    return 0.0;
}


PosType LensHaloNFW::dbfunction(PosType x){
    PosType p[]={0.72708293,0.22082787,1.41708596,-0.22678577,-0.00238127};
    if(x==0) x=1e-5;
    if(x>0){return p[0]/(pow(x,p[1])*(1+pow(x,p[2]))+p[3])+p[4];};
    return 0.0;
}

PosType LensHaloNFW::ddbfunction(PosType x){
    PosType p[]={0.72708293,0.22082787,1.41708596,-0.22678577,-0.00238127};
    if(x==0) x=1e-5;
    if(x>0){return (-1.0*p[0]*p[1]*pow(x,(p[1]-1.))*(1+(1+p[2]/p[1])*pow(x,p[2])))/pow(p[3]+pow(x,p[1])*(1+pow(x,p[2])),2.0);};
    return 0.0;
}






