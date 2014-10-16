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

	r=sqrt(pow(x[0]-center[0],2) + pow(x[1]-center[1],2));
	if(r==0.0){
		gamma[0]=gamma[1]=0.0;
		return ;
	}

	gt=mass/pi/Sigma_crit/pow(r,2);
	if(r<Rtrunc){
		PosType y;

		y = Rtrunc/r_scale;
		gt /= gfunction(y);
		y = r/r_scale;
		gt *= g2function(y)*pow(y,2);
	}

	gamma[0]=-gt*(pow(x[0]-center[0],2)-pow(x[1]-center[1],2))/r/r;
	gamma[1]=-2*gt*(x[0]-center[0])*(x[1]-center[1])/r/r;

	return ;
}

PosType LensHaloNFW::gfunction(PosType x) const{
	PosType ans;

	if(x==0) x=1e-5;
	ans=log(x/2);
	if(x==1.0){ ans += 1.0; return ans;}
	if(x>1.0){  ans +=  2*atan(sqrt((x-1)/(x+1)))/sqrt(x*x-1);; return ans;}
	if(x<1.0){  ans += 2*atanh(sqrt((1-x)/(x+1)))/sqrt(1-x*x);; return ans;}
	return 0.0;
}

PosType LensHaloNFW::dgfunctiondx(PosType x){
	PosType ans;
	if(x==0) x=1e-5;
	ans=1./x;
	if(x==1.0){ return 1./3.;}
	if(x>1.0){  ans += 1/(x*(x*x-1)) -  2*x*atan(sqrt((x-1)/(x+1)))/pow(x*x-1,1.5); return ans;}
	if(x<1.0){  ans += 1/(x*(x*x-1)) +  2*x*atanh(sqrt((1-x)/(x+1)))/pow(1-x*x,1.5); return ans;}
	return 0.0;
}


PosType LensHaloNFW::ffunction(PosType x) const{
	PosType ans;

	if(x==0) x=1e-5;
	if(x==1.0){ return 1.0/3.0;}
	if(x>1.0){  ans = (1-2*atan(sqrt((x-1)/(x+1)))/sqrt(x*x-1))/(x*x-1); return ans;}
	if(x<1.0){  ans = (1-2*atanh(sqrt((1-x)/(x+1)))/sqrt(1-x*x))/(x*x-1); return ans;}
	return 0.0;
}

PosType LensHaloNFW::g2function(PosType x) const{
	PosType ans,y;

	if(x==0) x=1e-5;
	if(x==1) return 10/3. + 4*log(0.5);
	y=x*x-1;
	ans=4*log(x/2)/x/x -2/y;
	if(x<1.0){ ans+= 4*(2/x/x+1/y)*atanh(sqrt((1-x)/(x+1)))/sqrt(-y); return ans;}
	if(x>1.0){ ans+= 4*(2/x/x+1/y)*atan(sqrt((x-1)/(x+1)))/sqrt(y); return ans;}

	return 0.0;
}

PosType LensHaloNFW::hfunction(PosType x) const{
	if(x==0) x=1e-5;
	if(x<=1.0) return 2.*((log(x/2.)*log(x/2.))-atanh(sqrt(1.-x*x))*atanh(sqrt(1.-x*x)));
	if(x>1.0) return 2.* ((log(x/2.)*log(x/2.))+ atan(sqrt(x*x-1.))*atan(sqrt(x*x-1.)));
	return 0.0;
}

PosType LensHaloNFW::dhfunction(PosType x) const { // d ln phi(x) / d ln(x) = beta
    if(x==0) x=1e-5;
    if(x<1.0)  return x*((((2*x*atanh(sqrt(1.-x*x)))/(sqrt(1.-x*x)*x*x)+(2.*log(x/2))/x))/(-pow(atanh(sqrt(1.-x*x)),2)+log(0.5*x)*log(0.5*x)));
    if(x==1.0) return 2.0*(log(0.5)+1.0)/(log(0.5)*log(0.5));
    if(x>1.0)  return x*((((2*x*atan(sqrt(-1.+x*x)))/(sqrt(-1.+x*x)*(x*x))+(2.*log(x/2))/x))/(atan(sqrt(-1.+x*x))*atan(sqrt(-1.+x*x))+log(x/2)*log(x/2)));
    return 0.0;
}

PosType LensHaloNFW::ddhfunction(PosType x, bool numerical){ // d dhfunction(x) / dx i.e. dbeta(x)/dx
    if(numerical==false){
        if(x==0) x=1e-5;
        PosType dh=dhfunction(x);
        PosType aux=-1.0*dh*dh/x;
        if(x<1.0)  return aux + (2./x-2./(x*(1-x*x))+2*x*atanh(sqrt(1.-x*x))/(pow(1.-x*x,3./2.)))/(-1.0*atanh(sqrt(1.-x*x))*atanh(sqrt(1.-x*x))+log(x/2.)*log(x/2.));
        if(x==1.0) return aux + 1.38758;
            if(x>1.0)  return aux + (2./x+2./(x*(x*x-1.))-2*x*atan(sqrt(x*x-1.))/(pow(x*x-1.,3./2.)))/(atan(sqrt(x*x-1.))*atan(sqrt(x*x-1.))+log(x/2.)*log(x/2.));
    }
    if(numerical==true){
        PosType h=1e-5;
        if(x==0) x=2e-5;
        if(x>0){return ((1./12.)*(dhfunction(x-2*h))-(2./3.)*(dhfunction(x-h))+(2./3.)*(dhfunction(x+h))-(1./12.)*(dhfunction(x+2*h)))/h;};
    }
    return 0.0;
}


PosType LensHaloNFW::dddhfunction(PosType x, bool numerical){ // d^2 dhfunction(x) / dx^2 i.e. d^2 beta(x) / dx^2
    if(numerical==false){
        PosType b=dhfunction(x);
        PosType db=ddhfunction(x,false);
        PosType aux=-b*db/x+(b*b/x/x)*(1-b);
        if(x==0) x=1e-5;
        if(x<1.0){
            PosType aux3=sqrt(1.-x*x);
            PosType aux2=(pow(atanh(aux3),2)-log(x/2)*log(x/2));
            return aux+(-(2/x/x)-6/(1-x*x)/(1-x*x)+2/(x*x*(1-x*x))+(6*x*x*atanh(aux3))/ (pow(1-x*x,5./2.)) +(2*atanh(aux3))/pow(aux3,3))/(-1.0*aux2);
        }
        if(x==1.0) return aux + 1.38758;
        if(x>1.0){
            PosType aux3=sqrt(-1.+x*x);
            PosType aux2=(pow(atan(aux3),2)+log(x/2)*log(x/2));
            return aux+(-(2/x/x)-6/(-1+x*x)/(1-x*x)-2/(x*x*(-1+x*x))+(6*x*x*atan(aux3))/ (pow(-1+x*x,5./2.)) -(2*atan(aux3))/pow(aux3,3))/(aux2);
        }
    }
    if(numerical==true){
        PosType h=1e-5;
        if(x==0) x=2e-5;
        if(x>0){return ((-1./12.)*dhfunction(x-2*h)+(4./3.)*dhfunction(x-h)-(5./2.)*dhfunction(x)+(4./3.)*dhfunction(x+h)-(1./12.)*dhfunction(x+2*h))/(h*h);};
    }
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

PosType LensHaloNFW::dmoddb(int whichmod, PosType q, PosType b){
    PosType h=1e-5;
    if(b==0) b=2e-5;
    if(b>0){return ((1./12.)*(InterpolateModes(whichmod, q, b-2*h))-(2./3.)*(InterpolateModes(whichmod, q, b-h))+(2./3.)*(InterpolateModes(whichmod, q, b+h))-(1./12.)*(InterpolateModes(whichmod, q, b+2*h)))/h;};
    return 0.0;
}

PosType LensHaloNFW::dmoddq(int whichmod, PosType q, PosType b){
    PosType h=1e-5;
    if(q==0) q=2e-5;
    if(q>0){return ((1./12.)*(InterpolateModes(whichmod, q-2*h, b))-(2./3.)*(InterpolateModes(whichmod, q-h, b))+(2./3.)*(InterpolateModes(whichmod, q+h, b))-(1./12.)*(InterpolateModes(whichmod, q+2*h, b)))/h;};
    return 0.0;
}

PosType LensHaloNFW::ddmoddb(int whichmod, PosType q, PosType b){
 PosType h=1e-5;
 if(b==0) b=2e-5;
 if(b>0){return ((-1./12.)*InterpolateModes(whichmod, q, b-2*h)+(4./3.)*InterpolateModes(whichmod, q, b-h)-(5./2.)*InterpolateModes(whichmod, q, b)+(4./3.)*InterpolateModes(whichmod, q, b+h)-(1./12.)*InterpolateModes(whichmod, q, b+2*h))/(h*h);};
 //if(x>0){return (bfunction(x-h)-2.*bfunction(x)+bfunction(x+h))/(h*h);};
 //if(x>0){return ((-1./12.)*bfunction(x-2*h)+(4./3.)*bfunction(x-h)-(5./2.)*bfunction(x)+(4./3.)*bfunction(x+h)-(1./12.)*bfunction(x+2*h))/(h*h);};
 //if(x>0){return ((1./90.)*bfunction(x-3*h)+(-3./20.)*bfunction(x-2*h)+(3./2.)*bfunction(x-h)-(49./18.)*bfunction(x)+(3./2.)*bfunction(x+h)-(3./20.)*bfunction(x+2*h)+(1./90.)*bfunction(x+3*h))/(h*h);};
 return 0.0;
 }

PosType LensHaloNFW::ddmoddq(int whichmod, PosType q, PosType b){
 PosType h=1e-5;
 if(q==0) q=2e-5;
 if(q>0){return ((-1./12.)*InterpolateModes(whichmod, q-2*h, b)+(4./3.)*InterpolateModes(whichmod, q-h, b)-(5./2.)*InterpolateModes(whichmod, q, b)+(4./3.)*InterpolateModes(whichmod, q+h, b)-(1./12.)*InterpolateModes(whichmod, q+2*h, b))/(h*h);};
 //if(x>0){return (bfunction(x-h)-2.*bfunction(x)+bfunction(x+h))/(h*h);};
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




