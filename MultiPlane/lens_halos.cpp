/*
 * lens_halos.cpp
 *
 *  Created on: 06.05.2013
 *      Author: mpetkova
 */

#include "slsimlib.h"

LensHalo::LensHalo(){

}

LensHalo::~LensHalo(){

}

NFWLensHalo::NFWLensHalo() : LensHalo(){
}

NFWLensHalo::~NFWLensHalo(){

}

void NFWLensHalo::set_internal(long *seed, float vmax, float r_halfmass){

	NFW_Utility nfw_util;

	// Find the NFW profile with the same mass, Vmax and R_halfmass
	nfw_util.match_nfw(vmax,r_halfmass,mass,&rscale,&Rmax);

	rscale = Rmax/rscale; // Was the concentration
}

PseudoNFWLensHalo::PseudoNFWLensHalo() : LensHalo(){

}

PseudoNFWLensHalo::~PseudoNFWLensHalo(){

}



PowerLawLensHalo::PowerLawLensHalo() : LensHalo(){
	rscale = 1.0;
}

PowerLawLensHalo::~PowerLawLensHalo(){

}


NSIELensHalo::NSIELensHalo() : LensHalo(){
	rscale = 1.0;
	rcore = 0.0;
}


NSIELensHalo::~NSIELensHalo(){

}

void NSIELensHalo::set_internal(long *seed, float vmax, float r_halfmass){
	sigma = 126*pow(mass/1.0e10,0.25); // From Tully-Fisher and Bell & de Jong 2001
	fratio = (ran2(seed)+1)*0.5;  //TODO This is a kluge.
	pa = 2*pi*ran2(seed);  //TODO This is a kluge.
	Rsize = rmaxNSIE(sigma,mass,fratio,rcore);
	Rmax = MAX(1.0,1.0/fratio)*Rsize;  // redefine

	assert(Rmax >= Rsize);
}

void LensHalo::force_halo(
		double *alpha     /// mass/Mpc
		,KappaType *kappa
		,KappaType *gamma
		,double *xcm
		,bool no_kappa
		){

	double rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
	if(rcm2 < 1e-20) rcm2 = 1e-20;

	/// intersecting, subtract the point particle
	if(rcm2 < Rmax*Rmax){
		double prefac = mass/rcm2/pi;
		double arg1 = sqrt(rcm2)/rscale;
		double arg2 = Rmax/rscale;

		double tmp = (alpha_h(arg1,arg2) + 1.0)*prefac;
		alpha[0] += tmp*xcm[0];
		alpha[1] += tmp*xcm[1];

		// can turn off kappa and gamma calculations to save times
		if(!no_kappa){
			*kappa += kappa_h(arg1,arg2)*prefac;

			tmp = (gamma_h(arg1,arg2) + 2.0)*prefac/rcm2;

			gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
			gamma[1] += xcm[0]*xcm[1]*tmp;
		}
	}

	return;
}

void NSIELensHalo::force_halo(
		double *alpha
		,KappaType *kappa
		,KappaType *gamma
		,double *xcm
		,bool no_kappa
){

	double rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
	if(rcm2 < 1e-20) rcm2 = 1e-20;

	double ellipR = ellipticRadiusNSIE(xcm,fratio,pa);
	if(ellipR > Rsize){
		double rout = Rsize*MAX(1.0,1.0/fratio);
		// This is the case when the ray is within the NSIE's circular region of influence but outside its elliptical truncation
		// if the ray misses the halo treat it as a point mass
		double prefac = -1.0*mass/rcm2/pi;

		if(rcm2 > rout*rout){
			alpha[0] += prefac*xcm[0];
			alpha[1] += prefac*xcm[1];
		}else{
			double alpha_out[2],alpha_in[2],rin,x_in[2];
			double prefac = -1.0*mass/rout/pi;
			double r = sqrt(rcm2);

			alpha_out[0] = prefac*xcm[0]/r;
			alpha_out[1] = prefac*xcm[1]/r;

			Utilities::rotation(x_in,xcm,pa);
			rin = r*Rsize
					/sqrt( x_in[0]*x_in[0] + pow(fratio*x_in[1],2) );
			//rin = Rsize;

			x_in[0] = rin*xcm[0]/r;
			x_in[1] = rin*xcm[1]/r;

			alpha_in[0] = alpha_in[1] = 0;
			float units = pow(sigma/lightspeed,2)/Grav/sqrt(fratio); // mass/distance(physical)
			alphaNSIE(alpha_in,x_in,fratio,rcore,pa);
			alpha_in[0] *= -units;
			alpha_in[1] *= -units;

			alpha[0] += (r - rin)*(alpha_out[0] - alpha_in[0])/(rout - rin) + alpha_in[0];
			alpha[1] += (r - rin)*(alpha_out[1] - alpha_in[1])/(rout - rin) + alpha_in[1];
		}

		// can turn off kappa and gamma calculations to save times
		if(!no_kappa){
			prefac *= 2.0/rcm2;

			gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*prefac;
			gamma[1] += xcm[0]*xcm[1]*prefac;
		}

	}else{
		double xt[2]={0,0},tmp[2]={0,0};
		float units = pow(sigma/lightspeed,2)/Grav/sqrt(fratio); // mass/distance(physical)
		xt[0]=xcm[0];
		xt[1]=xcm[1];
		alphaNSIE(tmp,xt,fratio,rcore,pa);
		alpha[0] -= units*tmp[0];
		alpha[1] -= units*tmp[1];
		if(!no_kappa){
			KappaType tmp[2]={0,0};
			*kappa += units*kappaNSIE(xt,fratio,rcore,pa);
			gammaNSIE(tmp,xt,fratio,rcore,pa);
			gamma[0] += units*tmp[0];
			gamma[1] += units*tmp[1];
		}
	}
	return;
}
