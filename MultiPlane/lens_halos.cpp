/*
 * lens_halos.cpp
 *
 *  Created on: 06.05.2013
 *      Author: mpetkova
 */

#include "slsimlib.h"

LensHalo::LensHalo(){
	rscale = 1.0;
	mass = Rmax = 0.0;
}

LensHalo::LensHalo(InputParams& params){
	assignParams(params);
}

void LensHalo::initFromMassFunc(float my_mass, double mass_scale, float my_Rmax, float my_rscale, double my_slope, long *seed){
	mass = my_mass/mass_scale;
	Rmax = my_Rmax;
	rscale = my_rscale;
}

void LensHalo::error_message1(std::string parameter,std::string file){
	ERROR_MESSAGE();
	std::cout << "Parameter " << parameter << " is needed to construct a LensHalo.  It needs to be set in parameter file " << file << "!" << std::endl;
	exit(0);
}

void LensHalo::assignParams(InputParams& params){
	if(!params.get("mass",mass)) error_message1("mass",params.filename());
	if(!params.get("Rmax",Rmax)) error_message1("Rmax",params.filename());
	if(!params.get("z_lens",zlens)) error_message1("z_lens",params.filename());
}

LensHalo::~LensHalo(){

}

int NFWLensHalo::count = 0;
double *NFWLensHalo::xtable = NULL,*NFWLensHalo::ftable = NULL,*NFWLensHalo::gtable = NULL,*NFWLensHalo::g2table = NULL;
NFWLensHalo::NFWLensHalo() : LensHalo(){
	make_tables();
}

NFWLensHalo::NFWLensHalo(InputParams& params){
	assignParams(params);
	make_tables();
}

void NFWLensHalo::make_tables(){
	if(count == 0){
		int i;
		double x, dx = maxrm/(double)NTABLE;

		xtable = new double[NTABLE];
		ftable = new double[NTABLE];
		gtable = new double[NTABLE];
		g2table = new double[NTABLE];

		for(i = 0 ; i< NTABLE; i++){
			x = i*dx;
			xtable[i] = x;
			ftable[i] = ffunction(x);
			gtable[i] = gfunction(x);
			g2table[i] = g2function(x);
		}

		count++;
	}
}

double NFWLensHalo::InterpolateFromTable(double *table, double y){
	int j;
	j=(int)(y/maxrm*NTABLE);

	assert(y>=xtable[j] && y<=xtable[j+1]);
	return (table[j+1]-table[j])/(xtable[j+1]-xtable[j])*(y-xtable[j]) + table[j];
}

void NFWLensHalo::assignParams(InputParams& params){
	if(!params.get("mass",mass)) error_message1("mass",params.filename());
	if(!params.get("Rmax",Rmax)) error_message1("Rmax",params.filename());
	if(!params.get("z_lens",zlens)) error_message1("z_lens",params.filename());
	if(!params.get("concentration",rscale)) error_message1("concentration",params.filename());
	rscale = rscale*Rmax; // was the concentration
}

NFWLensHalo::~NFWLensHalo(){
	--count;
	if(count == 0){
		delete[] xtable;
		delete[] gtable;
		delete[] ftable;
		delete[] g2table;
	}
}

void NFWLensHalo::initFromFile(float my_mass, double mass_scale, long *seed, float vmax, float r_halfmass){

	mass = my_mass;

	NFW_Utility nfw_util;

	// Find the NFW profile with the same mass, Vmax and R_halfmass
	nfw_util.match_nfw(vmax,r_halfmass,mass,&rscale,&Rmax);

	rscale = Rmax/rscale; // Was the concentration
	mass /= mass_scale;
}

int PseudoNFWLensHalo::count = 0;
double *PseudoNFWLensHalo::xtable = NULL,*PseudoNFWLensHalo::mhattable = NULL;
PseudoNFWLensHalo::PseudoNFWLensHalo() : LensHalo(){
}

PseudoNFWLensHalo::PseudoNFWLensHalo(InputParams& params){
	assignParams(params);
	make_tables();
}

void PseudoNFWLensHalo::make_tables(){
	if(count == 0){
		int i;
		double x, dx = maxrm/(double)NTABLE;

		xtable = new double[NTABLE];
		mhattable = new double[NTABLE];

		for(i = 0 ; i< NTABLE; i++){
			x = i*dx;
			xtable[i] = x;
			mhattable[i] = mhat(x,beta);
		}

		count++;
	}
}

double PseudoNFWLensHalo::InterpolateFromTable(double y){
	int j;
	j=(int)(y/maxrm*NTABLE);

	assert(y>=xtable[j] && y<=xtable[j+1]);
	return (mhattable[j+1]-mhattable[j])/(xtable[j+1]-xtable[j])*(y-xtable[j]) + mhattable[j];
}

void PseudoNFWLensHalo::initFromMassFunc(float my_mass, double mass_scale, float my_Rmax, float my_rscale, double my_slope, long *seed){
	LensHalo::initFromMassFunc(my_mass,mass_scale,my_Rmax,my_rscale,my_slope,seed);
	beta = my_slope;
	make_tables();
}

void PseudoNFWLensHalo::assignParams(InputParams& params){
	if(!params.get("mass",mass)) error_message1("mass",params.filename());
	if(!params.get("Rmax",Rmax)) error_message1("Rmax",params.filename());
	if(!params.get("z_lens",zlens)) error_message1("z_lens",params.filename());
	if(!params.get("concentration",rscale)) error_message1("concentration",params.filename());
	if(!params.get("slope_pnfw",beta)) error_message1("slope_pnfw",params.filename());
	rscale = rscale*Rmax; // was the concentration
}

PseudoNFWLensHalo::~PseudoNFWLensHalo(){
    --count;
    if(count == 0){
    	delete[] xtable;
    	delete[] mhattable;
    }
}



PowerLawLensHalo::PowerLawLensHalo() : LensHalo(){
	rscale = 1.0;
}

PowerLawLensHalo::PowerLawLensHalo(InputParams& params){
	assignParams(params);
}

void PowerLawLensHalo::initFromMassFunc(float my_mass, double mass_scale, float my_Rmax, float my_rscale, double my_slope, long *seed){
	LensHalo::initFromMassFunc(my_mass,mass_scale,my_Rmax,my_rscale,my_slope,seed);
	beta = my_slope;
}

void PowerLawLensHalo::assignParams(InputParams& params){
	if(!params.get("mass",mass)) error_message1("mass",params.filename());
	if(!params.get("Rmax",Rmax)) error_message1("Rmax",params.filename());
	if(!params.get("z_lens",zlens)) error_message1("z_lens",params.filename());
	if(!params.get("slope_pl",beta)) error_message1("slope_pl",params.filename());
}

PowerLawLensHalo::~PowerLawLensHalo(){

}


SimpleNSIELensHalo::SimpleNSIELensHalo() : LensHalo(){
	rscale = 1.0;
}

SimpleNSIELensHalo::SimpleNSIELensHalo(InputParams& params){
	assignParams(params);
}

void SimpleNSIELensHalo::assignParams(InputParams& params){
	if(!params.get("mass",mass)) error_message1("mass",params.filename());
	if(!params.get("Rmax",Rmax)) error_message1("Rmax",params.filename());
	if(!params.get("z_lens",zlens)) error_message1("z_lens",params.filename());

	if(!params.get("sigma",sigma)) error_message1("sigma",params.filename());
	if(!params.get("core",rcore)) error_message1("core",params.filename());
	if(!params.get("axis_ratio",fratio)) error_message1("axis_ratio",params.filename());
	if(!params.get("pos_angle",pa)) error_message1("pos_angle",params.filename());

	Rsize = rmaxNSIE(sigma,mass,fratio,rcore);
	Rmax = MAX(1.0,1.0/fratio)*Rsize;  // redefine

	assert(Rmax >= Rsize);
}

SimpleNSIELensHalo::~SimpleNSIELensHalo(){

}
void SimpleNSIELensHalo::initFromMass(float my_mass, double mass_scale, long *seed){
	mass = my_mass;
	rcore = 0.0;
	sigma = 126*pow(mass/1.0e10,0.25); // From Tully-Fisher and Bell & de Jong 2001
	fratio = (ran2(seed)+1)*0.5;  //TODO This is a kluge.
	pa = 2*pi*ran2(seed);  //TODO This is a kluge.
	Rsize = rmaxNSIE(sigma,mass,fratio,rcore);
	Rmax = MAX(1.0,1.0/fratio)*Rsize;  // redefine

	mass /= mass_scale;
	assert(Rmax >= Rsize);
}

void SimpleNSIELensHalo::initFromFile(float my_mass, double mass_scale, long *seed, float vmax, float r_halfmass){
	initFromMass(my_mass,mass_scale,seed);
}

void SimpleNSIELensHalo::initFromMassFunc(float my_mass, double mass_scale, float my_Rmax, float my_rscale, double my_slope, long *seed){
	initFromMass(my_mass,mass_scale,seed);
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
		double x = sqrt(rcm2)/rscale;
		double xmax = Rmax/rscale;

		double tmp = (alpha_h(x,xmax) + 1.0)*prefac;
		alpha[0] += tmp*xcm[0];
		alpha[1] += tmp*xcm[1];

		// can turn off kappa and gamma calculations to save times
		if(!no_kappa){
			*kappa += kappa_h(x,xmax)*prefac;

			tmp = (gamma_h(x,xmax) + 2.0)*prefac/rcm2;

			gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
			gamma[1] += xcm[0]*xcm[1]*tmp;
		}
	}

	return;
}

void SimpleNSIELensHalo::force_halo(
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
