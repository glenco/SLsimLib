/*
 * powerlaw.c
 *
 *  Created on: Mar 8, 2010
 *      Author: R.B. Metcalf
 *
 *      powerlaw gives the deflection etc. for a truncated power-law
 *      density profile
 *     convention here is gamma_1 = -(Axx-Ayy)/2 and gamma_2= -Axy
 *
 */

#include "lens_halos.h"

///
void alphaPowLaw(PosType *alpha,PosType *x,PosType R,PosType mass,PosType beta,PosType *center,PosType Sigma_crit){
	PosType r,b=0;

	r=sqrt(pow(x[0]-center[0],2) + pow(x[1]-center[1],2));
	if(r==0){
		alpha[0]=alpha[1]=0.0;
		return ;
	}
	b=mass/pow(r,2)/PI/Sigma_crit;
	if(r<R) b *= pow(r/R,beta+2);

	alpha[0]=b*(x[0]-center[0]);
	alpha[1]=b*(x[1]-center[1]);

	return ;
}
///
KappaType kappaPowLaw(PosType *x,PosType R,PosType mass,PosType beta,PosType *center,PosType Sigma_crit){
	PosType r;

	r=sqrt(pow(x[0]-center[0],2) + pow(x[1]-center[1],2));
	if(r>R) return 0.0;
	if(r < 1.0-20) r=1.0e-20;
	return (beta+2)*mass*pow(r/R,beta)/(2*PI*pow(R,2)*Sigma_crit);
}
///
void gammaPowLaw(KappaType *gamma,PosType *x,PosType R,PosType mass,PosType beta
		,PosType *center,PosType Sigma_crit){
	PosType r,gt=0;

	r=sqrt(pow(x[0]-center[0],2) + pow(x[1]-center[1],2));
	if(r==0.0){
		gamma[0]=gamma[1]=0.0;
		return ;
	}
	gt=mass/PI/Sigma_crit/pow(r,2);
	if(r<R) gt *= -beta*pow(r/R,beta+2)/2;

	gamma[0]=-gt*(pow(x[0]-center[0],2)-pow(x[1]-center[1],2))/r/r;
	gamma[1]=-2*gt*(x[0]-center[0])*(x[1]-center[1])/r/r;

	return ;
}
///
KappaType phiPowLaw(PosType *x,PosType R,PosType mass,PosType beta
		,PosType *center,PosType Sigma_crit){
	PosType b,r;

	r=sqrt(pow(x[0]-center[0],2) + pow(x[1]-center[1],2));

	b=mass/PI/Sigma_crit;
	if(r<=R) return b*pow(r/R,beta+2);
	return b*(log(r/R) + 1);
}

