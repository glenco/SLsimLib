/*
 * lens_halos.cpp
 *
 *  Created on: 06.05.2013
 *      Author: mpetkova
 */

#include "slsimlib.h"

LensHalo::LensHalo(){
	rscale = 1.0;
	mass = Rmax = xmax = 0.0;
}

LensHalo::LensHalo(InputParams& params){
	assignParams(params);
}

void LensHalo::initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, double my_slope, long *seed){
	mass = my_mass;
	Rmax = my_Rmax;
	rscale = my_rscale;
  xmax = Rmax/rscale;
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

int LensHaloNFW::count = 0;
double *LensHaloNFW::xtable = NULL,*LensHaloNFW::ftable = NULL,*LensHaloNFW::gtable = NULL,*LensHaloNFW::g2table = NULL;

LensHaloNFW::LensHaloNFW() : LensHalo(){
  gmax=0;
	make_tables();
  gmax = InterpolateFromTable(gtable,xmax);
}

LensHaloNFW::LensHaloNFW(InputParams& params){
	assignParams(params);
	make_tables();
  gmax = InterpolateFromTable(gtable,xmax);
}

void LensHaloNFW::make_tables(){
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
  }
  count++;
}

double LensHaloNFW::InterpolateFromTable(double *table, double y){
	int j;
	j=(int)(y/maxrm*NTABLE);

	assert(y>=xtable[j] && y<=xtable[j+1]);
	return (table[j+1]-table[j])/(xtable[j+1]-xtable[j])*(y-xtable[j]) + table[j];
}

void LensHaloNFW::assignParams(InputParams& params){
	if(!params.get("mass_nfw",mass)) error_message1("mass_nfw",params.filename());
	if(!params.get("Rmax_nfw",Rmax)) error_message1("Rmax_nfw",params.filename());
	if(!params.get("zlens_nfw",zlens)) error_message1("lens_nfw",params.filename());
	if(!params.get("concentration_nfw",rscale)) error_message1("concentration_nfw",params.filename());
	rscale = Rmax/rscale; // was the concentration
  xmax = Rmax/rscale;

}

LensHaloNFW::~LensHaloNFW(){
	--count;
	if(count == 0){
		delete[] xtable;
		delete[] gtable;
		delete[] ftable;
		delete[] g2table;
	}
}

/// Sets the profile to match the mass, Vmax and R_halfmass
void LensHaloNFW::initFromFile(float my_mass, long *seed, float vmax, float r_halfmass){

	mass = my_mass;

	NFW_Utility nfw_util;

	// Find the NFW profile with the same mass, Vmax and R_halfmass
	nfw_util.match_nfw(vmax,r_halfmass,mass,&rscale,&Rmax);
	rscale = Rmax/rscale; // Was the concentration
  xmax = Rmax/rscale;
  gmax = InterpolateFromTable(gtable,xmax);
}

int LensHaloPseudoNFW::count = 0;
double *LensHaloPseudoNFW::xtable = NULL,*LensHaloPseudoNFW::mhattable = NULL;
LensHaloPseudoNFW::LensHaloPseudoNFW() : LensHalo(){
}

LensHaloPseudoNFW::LensHaloPseudoNFW(InputParams& params){
	assignParams(params);
	make_tables();
}

/// Auxiliary function for PseudoNFW profile
// previously defined in tables.cpp
double LensHaloPseudoNFW::mhat(double y, double beta){
  if(y==0) y=1e-5;
	if(beta == 1.0) return y - log(1+y);
	if(beta == 2.0) return log(1+y) - y/(1+y);
	if(beta>=3.0) return ( (1 - beta)*y + pow(1+y,beta-1) - 1)/(beta-2)/(beta-1)/pow(1+y,beta-1);

	ERROR_MESSAGE();
	std::cout << "Only beta ==1, ==2 and >=3 are valid" << std::endl;
	exit(1);
	return 0.0;
}


void LensHaloPseudoNFW::make_tables(){
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

double LensHaloPseudoNFW::InterpolateFromTable(double y){
	int j;
	j=(int)(y/maxrm*NTABLE);

	assert(y>=xtable[j] && y<=xtable[j+1]);
	return (mhattable[j+1]-mhattable[j])/(xtable[j+1]-xtable[j])*(y-xtable[j]) + mhattable[j];
}

void LensHaloPseudoNFW::initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, double my_slope, long *seed){
	LensHalo::initFromMassFunc(my_mass,my_Rmax,my_rscale,my_slope,seed);
	beta = my_slope;
  xmax = my_Rmax/my_rscale;
	make_tables();
}

void LensHaloPseudoNFW::assignParams(InputParams& params){
	if(!params.get("mass_pnfw",mass)) error_message1("mass_pnfw",params.filename());
	if(!params.get("Rmax_pnfw",Rmax)) error_message1("Rmax_pnfw",params.filename());
	if(!params.get("zlens_pnfw",zlens)) error_message1("zlens_pnfw",params.filename());
	if(!params.get("concentration_pnfw",rscale)) error_message1("concentration_pnfw",params.filename());
	if(!params.get("slope_pnfw",beta)) error_message1("slope_pnfw",params.filename());
	rscale = Rmax/rscale; // was the concentration
  xmax = Rmax/rscale;
}

LensHaloPseudoNFW::~LensHaloPseudoNFW(){
    --count;
    if(count == 0){
    	delete[] xtable;
    	delete[] mhattable;
    }
}



LensHaloPowerLaw::LensHaloPowerLaw() : LensHalo(){
	rscale = 1.0;
}

LensHaloPowerLaw::LensHaloPowerLaw(InputParams& params){
	assignParams(params);
}

void LensHaloPowerLaw::initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, double my_slope, long *seed){
	LensHalo::initFromMassFunc(my_mass,my_Rmax,my_rscale,my_slope,seed);
	beta = my_slope;
  xmax = my_Rmax/my_rscale;
}

void LensHaloPowerLaw::assignParams(InputParams& params){
	if(!params.get("mass_pl",mass)) error_message1("mass_pl",params.filename());
	if(!params.get("Rmax_pl",Rmax)) error_message1("Rmax_pl",params.filename());
	if(!params.get("zlens_pl",zlens)) error_message1("zlens_pl",params.filename());
	if(!params.get("slope_pl",beta)) error_message1("slope_pl",params.filename());
}

LensHaloPowerLaw::~LensHaloPowerLaw(){

}


LensHaloSimpleNSIE::LensHaloSimpleNSIE() : LensHalo(){
	rscale = 1.0;
}

LensHaloSimpleNSIE::LensHaloSimpleNSIE(InputParams& params){
	assignParams(params);
}

void LensHaloSimpleNSIE::assignParams(InputParams& params){
	if(!params.get("mass_nsie",mass)) error_message1("mass_nsie",params.filename());
	if(!params.get("zlens_nsie",zlens)) error_message1("zlens_nsie",params.filename());

	if(!params.get("sigma",sigma)) error_message1("sigma",params.filename());
	if(!params.get("core",rcore)) error_message1("core",params.filename());
	if(!params.get("axis_ratio",fratio)) error_message1("axis_ratio",params.filename());
	if(!params.get("pos_angle",pa)) error_message1("pos_angle",params.filename());

	Rsize = rmaxNSIE(sigma,mass,fratio,rcore);
	Rmax = MAX(1.0,1.0/fratio)*Rsize;  // redefine

	assert(Rmax >= Rsize);
}

LensHaloSimpleNSIE::~LensHaloSimpleNSIE(){

}
void LensHaloSimpleNSIE::initFromMass(float my_mass, long *seed){
	mass = my_mass;
	rcore = 0.0;
	sigma = 126*pow(mass/1.0e10,0.25); // From Tully-Fisher and Bell & de Jong 2001
	fratio = (ran2(seed)+1)*0.5;  //TODO This is a kluge.
	pa = 2*pi*ran2(seed);  //TODO This is a kluge.
	Rsize = rmaxNSIE(sigma,mass,fratio,rcore);
	Rmax = MAX(1.0,1.0/fratio)*Rsize;  // redefine

	assert(Rmax >= Rsize);
}

void LensHaloSimpleNSIE::initFromFile(float my_mass, long *seed, float vmax, float r_halfmass){
	initFromMass(my_mass,seed);
}

void LensHaloSimpleNSIE::initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, double my_slope, long *seed){
	initFromMass(my_mass,seed);
}

void LensHalo::force_halo(
		double *alpha     /// mass/Mpc
		,KappaType *kappa
		,KappaType *gamma
		,double *xcm
		,bool no_kappa
		,bool subtract_point /// if true contribution from a point mass is subtracted
		){

	double rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
	if(rcm2 < 1e-20) rcm2 = 1e-20;

	/// intersecting, subtract the point particle
	if(rcm2 < Rmax*Rmax){
		double prefac = mass/rcm2/pi;
		double x = sqrt(rcm2)/rscale;
		//double xmax = Rmax/rscale;

    //double tmp = (alpha_h(x,xmax) + 1.0*subtract_point)*prefac;
		double tmp = (alpha_h(x) + 1.0*subtract_point)*prefac;
		alpha[0] += tmp*xcm[0];
		alpha[1] += tmp*xcm[1];

		// can turn off kappa and gamma calculations to save times
		if(!no_kappa){
			*kappa += kappa_h(x)*prefac;

			tmp = (gamma_h(x) + 2.0*subtract_point)*prefac/rcm2;

			gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
			gamma[1] += xcm[0]*xcm[1]*tmp;
		}

	}
	else
	{
		if (subtract_point == false)
		{
			double prefac = mass/rcm2/pi;
			alpha[0] += -1.0*prefac*xcm[0];
			alpha[1] += -1.0*prefac*xcm[1];

			// can turn off kappa and gamma calculations to save times
			if(!no_kappa){
				double tmp = -2.0*prefac/rcm2;

				gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
				gamma[1] += xcm[0]*xcm[1]*tmp;
			}
		}
	}

	return;
}

void LensHaloSimpleNSIE::force_halo(
		double *alpha
		,KappaType *kappa
		,KappaType *gamma
		,double *xcm
		,bool no_kappa
		,bool subtract_point /// if true contribution from a point mass is subtracted
){

	double rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
	if(rcm2 < 1e-20) rcm2 = 1e-20;

  //**** test line

	if(rcm2 < Rmax*Rmax){
		double ellipR = ellipticRadiusNSIE(xcm,fratio,pa);
		if(ellipR > Rsize){
			// This is the case when the ray is within the NSIE's circular region of influence but outside its elliptical truncation

			double alpha_out[2],alpha_in[2],rin,x_in[2];
			double prefac = -1.0*mass/Rmax/pi;
			double r = sqrt(rcm2);

			alpha_out[0] = prefac*xcm[0]/r;
			alpha_out[1] = prefac*xcm[1]/r;

			rin = r*Rsize/ellipR;

			x_in[0] = rin*xcm[0]/r;
			x_in[1] = rin*xcm[1]/r;

			alpha_in[0] = alpha_in[1] = 0;
			float units = pow(sigma/lightspeed,2)/Grav/sqrt(fratio); // mass/distance(physical)
			alphaNSIE(alpha_in,x_in,fratio,rcore,pa);
			alpha_in[0] *= units;  // minus sign removed because already included in alphaNSIE
			alpha_in[1] *= units;
      
			alpha[0] += (r - rin)*(alpha_out[0] - alpha_in[0])/(Rmax - rin) + alpha_in[0];
			alpha[1] += (r - rin)*(alpha_out[1] - alpha_in[1])/(Rmax - rin) + alpha_in[1];
			//alpha[0] -= (r - rin)*(alpha_out[0] - alpha_in[0])/(Rmax - rin) + alpha_in[0];
			//alpha[1] -= (r - rin)*(alpha_out[1] - alpha_in[1])/(Rmax - rin) + alpha_in[1];

		}else{
			double xt[2]={0,0},tmp[2]={0,0};
			float units = pow(sigma/lightspeed,2)/Grav/sqrt(fratio); // mass/distance(physical)
			xt[0]=xcm[0];
			xt[1]=xcm[1];
			alphaNSIE(tmp,xt,fratio,rcore,pa);
           
			alpha[0] = units*tmp[0];  // minus sign removed because already included in alphaNSIE
			alpha[1] = units*tmp[1];
			//alpha[0] += units*tmp[0];
			//alpha[1] += units*tmp[1];
			if(!no_kappa){
				KappaType tmp[2]={0,0};
				*kappa += units*kappaNSIE(xt,fratio,rcore,pa);
				gammaNSIE(tmp,xt,fratio,rcore,pa);
				gamma[0] += units*tmp[0];
				gamma[1] += units*tmp[1];
			}
		}
    
		if(subtract_point){
			double fac = mass/rcm2/pi;
			alpha[0] += fac*xcm[0];
			alpha[1] += fac*xcm[1];

			// can turn off kappa and gamma calculations to save times
			if(!no_kappa){
				fac = 2.0*fac/rcm2;

				gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*fac;
				gamma[1] += xcm[0]*xcm[1]*fac;
			}

		}

	}
	return;
}


int LensHaloHernquist::count = 0;
double *LensHaloHernquist::xtable = NULL,*LensHaloHernquist::ftable = NULL,*LensHaloHernquist::gtable = NULL,*LensHaloHernquist::g2table = NULL;

LensHaloHernquist::LensHaloHernquist() : LensHalo(){
  gmax=0;
	make_tables();
  gmax = InterpolateFromTable(gtable,xmax);
}

LensHaloHernquist::LensHaloHernquist(InputParams& params){
	assignParams(params);
	make_tables();
  gmax = InterpolateFromTable(gtable,xmax);
}

void LensHaloHernquist::make_tables(){
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
  }
  count++;
}

double LensHaloHernquist::InterpolateFromTable(double *table, double y){
	int j;
	j=(int)(y/maxrm*NTABLE);

	assert(y>=xtable[j] && y<=xtable[j+1]);
	return (table[j+1]-table[j])/(xtable[j+1]-xtable[j])*(y-xtable[j]) + table[j];
}

void LensHaloHernquist::assignParams(InputParams& params){
	if(!params.get("mass_hernquist",mass)) error_message1("mass_hernquist",params.filename());
	if(!params.get("Rmax_hernquist",Rmax)) error_message1("Rmax_hernquist",params.filename());
	if(!params.get("zlens_hernquist",zlens)) error_message1("zlens_hernquist",params.filename());
	if(!params.get("rscale_hernquist",rscale)) error_message1("rscale_hernquist",params.filename());
  xmax = Rmax/rscale;
}

LensHaloHernquist::~LensHaloHernquist(){
	--count;
	if(count == 0){
		delete[] xtable;
		delete[] gtable;
		delete[] ftable;
		delete[] g2table;
	}
}


LensHaloDummy::LensHaloDummy()
: LensHalo()
{
//	mass = 0.;
}

LensHaloDummy::LensHaloDummy(InputParams& params)
: LensHalo()
{
	assignParams(params);
//	mass = 0.;
}

void LensHaloDummy::initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, double my_slope, long *seed){
	mass = 1.e-10;
	Rmax = my_Rmax;
	rscale = my_rscale;
  xmax = Rmax/rscale;
}


void LensHaloDummy::force_halo(double *alpha,KappaType *kappa,KappaType *gamma,double *xcm,bool no_kappa,bool subtract_point)
{
	double rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
	double prefac = mass/rcm2/pi;
	double tmp = subtract_point*prefac;
	alpha[0] += tmp*xcm[0];
	alpha[1] += tmp*xcm[1];
	double x = sqrt(rcm2)/rscale;

	// can turn off kappa and gamma calculations to save times
	if(!no_kappa){
		*kappa += kappa_h(x)*prefac;

		tmp = (gamma_h(x) + 2.0*subtract_point)*prefac/rcm2;

		gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
		gamma[1] += xcm[0]*xcm[1]*tmp;
	}

}

void LensHaloDummy::assignParams(InputParams& params)
{
	if(!params.get("z_lens",zlens)) error_message1("z_lens",params.filename());
}
