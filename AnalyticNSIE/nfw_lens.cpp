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
