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


long NFWLensHalo::ob_count = 0;
double *NFWLensHalo::ftable = NULL,*NFWLensHalo::gtable = NULL,*NFWLensHalo::g2table = NULL,*NFWLensHalo::xtable = NULL;

NFWLensHalo::NFWLensHalo() : LensHalo(){
	// make halo profile lookup tables if this is the first instance of a NFWLensHalo
	if(ob_count == 0) make_tables();
	++ob_count;
}

NFWLensHalo::~NFWLensHalo(){
	--ob_count;
	if(ob_count == 0){
		delete[] xtable;
		delete[] gtable;
		delete[] ftable;
		delete[] g2table;
	}
}

void NFWLensHalo::set_internal(long *seed, float vmax, float r_halfmass){

	NFW_Utility nfw_util;

	// Find the NFW profile with the same mass, Vmax and R_halfmass
	nfw_util.match_nfw(vmax,r_halfmass,mass,&rscale,&Rmax);

	rscale = Rmax/rscale; // Was the concentration
}

long PseudoNFWLensHalo::ob_count = 0;
double * PseudoNFWLensHalo::mhattable = NULL,*PseudoNFWLensHalo::xtable = NULL;

PseudoNFWLensHalo::PseudoNFWLensHalo() : NFWLensHalo(){
	// make halo profile lookup tables if this is the first instance of a QuadTreePseudoNFW
	if(ob_count == 0) make_tables();
	++ob_count;
}

PseudoNFWLensHalo::~PseudoNFWLensHalo(){
	--ob_count;
	if(ob_count == 0){
		delete[] xtable;
		delete[] mhattable;
	}
}



PowerLawLensHalo::PowerLawLensHalo() : LensHalo(){
	rscale = 1.0;
}

PowerLawLensHalo::~PowerLawLensHalo(){

}


NSIELensHalo::NSIELensHalo() : LensHalo(){

}


NSIELensHalo::~NSIELensHalo(){

}

void NSIELensHalo::set_internal(long *seed, float vmax, float r_halfmass){
	rscale = 1.0;
	rcore = 0.0;

	sigma = 126*pow(mass/1.0e10,0.25); // From Tully-Fisher and Bell & de Jong 2001
	fratio = (ran2(seed)+1)*0.5;  //TODO This is a kluge.
	pa = 2*pi*ran2(seed);  //TODO This is a kluge.
	Rsize = rmaxNSIE(sigma,mass,fratio,rcore);
	Rmax = MAX(1.0,1.0/fratio)*Rsize;  // redefine

	assert(Rmax >= Rsize);
}


NFW_NSIELensHalo::NFW_NSIELensHalo() : NFWLensHalo(){
	nsie_halo = new NSIELensHalo;
}


NFW_NSIELensHalo::~NFW_NSIELensHalo(){
	delete nsie_halo;
}

void NFW_NSIELensHalo::set_internal(long *seed, float vmax, float r_halfmass){
	double mo=7.3113e10,M1=2.8575e10,gam1=7.17,gam2=0.201,beta=0.557;

	NFWLensHalo::set_internal(seed,vmax,r_halfmass);

	double gmf = 2*mo*pow(mass/M1,gam1)
	  /pow(1+pow(mass/M1,beta),(gam1-gam2)/beta)/mass;
	if(gmf > 1.0) gmf = 1;

	nsie_halo->set_mass(mass*gmf);   //TODO This is a kluge. A mass dependent ratio would be better

	mass *= (1-gmf);

	nsie_halo->set_internal(seed,vmax,r_halfmass);
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

void NFW_NSIELensHalo::force_halo(
		double *alpha
		,KappaType *kappa
		,KappaType *gamma
		,double *xcm
		,bool no_kappa
){
	double alpha_tmp[2];
	KappaType kappa_tmp;
	KappaType gamma_tmp[3];

	LensHalo::force_halo(alpha,kappa,gamma,xcm,no_kappa);
	nsie_halo->force_halo(&alpha_tmp[0],&kappa_tmp,&gamma_tmp[0],xcm,no_kappa);

	alpha[0]+=alpha_tmp[0];
	alpha[1]+=alpha_tmp[1];
	if(!no_kappa){
		*kappa+=kappa_tmp;
		gamma[0]+=gamma_tmp[0];
		gamma[1]+=gamma_tmp[1];
		gamma[2]+=gamma_tmp[2];
	}
}
