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
	stars_implanted = false;
}

LensHalo::LensHalo(InputParams& params){
	assignParams(params);
	stars_implanted = false;
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

void LensHalo::PrintStars(bool show_stars) const
{
std::cout << std::endl << "Nstars "<<stars_N << std::endl << std::endl;
if(stars_N>0){
	if(star_Nregions > 0)
		std::cout << "stars_Nregions "<<star_Nregions << std::endl;
	std::cout << "stars_massscale "<<star_massscale << std::endl;
	std::cout << "stars_fstars "<<star_fstars << std::endl;
	std::cout << "stars_theta_force "<<star_theta_force << std::endl;
	if(show_stars){
		if(stars_implanted){
		  for(int i=0 ; i < stars_N ; ++i) std::cout << "    x["<<i<<"]="
						    << stars_xp[i][0] << " " << stars_xp[i][1] << std::endl;
		}else std::cout << "stars are not implanted yet" << std::endl;
	}
}
}

void LensHalo::force_stars(
		double *alpha     /// mass/Mpc
		,KappaType *kappa
		,KappaType *gamma
		,double *xcm
		,bool no_kappa
		)
{
    double alpha_tmp[2];
    KappaType gamma_tmp[3], tmp = 0;

    gamma_tmp[0] = gamma_tmp[1] = gamma_tmp[2] = 0.0;
    alpha_tmp[0] = alpha_tmp[1] = 0.0;

 	 substract_stars_disks(xcm,alpha,kappa,gamma);

	 // do stars with tree code
    star_tree->force2D_recur(xcm,alpha_tmp,&tmp,gamma_tmp,no_kappa);

    alpha[0] -= star_massscale*alpha_tmp[0];
    alpha[1] -= star_massscale*alpha_tmp[1];

    if(!no_kappa){
      *kappa += star_massscale*tmp;
      gamma[0] -= star_massscale*gamma_tmp[0];
      gamma[1] -= star_massscale*gamma_tmp[1];
    }
  
}

LensHalo::~LensHalo()
{
}

const long LensHaloNFW::NTABLE = 10000;
const double LensHaloNFW::maxrm = 100.0;
int LensHaloNFW::count = 0;

double* LensHaloNFW::xtable = NULL;
double* LensHaloNFW::ftable = NULL;
double* LensHaloNFW::gtable = NULL;
double* LensHaloNFW::g2table = NULL;
double* LensHaloNFW::htable = NULL;

LensHaloNFW::LensHaloNFW()
: LensHalo(), gmax(0)
{
	make_tables();
	gmax = InterpolateFromTable(gtable, xmax);
}

LensHaloNFW::LensHaloNFW(InputParams& params)
{
	assignParams(params);
	make_tables();
	gmax = InterpolateFromTable(gtable, xmax);
}

void LensHaloNFW::make_tables(){
	if(count == 0){
		int i;
		double x, dx = maxrm/(double)NTABLE;

		xtable = new double[NTABLE];
		ftable = new double[NTABLE];
		gtable = new double[NTABLE];
		g2table = new double[NTABLE];
		htable = new double[NTABLE];

		for(i = 0 ; i< NTABLE; i++){
			x = i*dx;
			xtable[i] = x;
			ftable[i] = ffunction(x);
			gtable[i] = gfunction(x);
			g2table[i] = g2function(x);
			htable[i] = hfunction(x);
		}
  }
  count++;
}

double LensHaloNFW::InterpolateFromTable(double *table, double y){
	int j;
	j=(int)(y/maxrm*NTABLE);

	assert(y>=xtable[j] && y<=xtable[j+1]);
	if (j==0)
		{
		if (table==ftable) return ffunction(y);
		if (table==gtable) return gfunction(y);
		if (table==g2table) return g2function(y);
		if (table==htable) return hfunction(y);
		}
	return (table[j+1]-table[j])/(xtable[j+1]-xtable[j])*(y-xtable[j]) + table[j];
}

void LensHaloNFW::assignParams(InputParams& params){
	if(!params.get("main_mass",mass)) error_message1("main_mass",params.filename());
	if(!params.get("main_Rmax",Rmax)) error_message1("main_Rmax",params.filename());
	if(!params.get("main_zlens",zlens)) error_message1("main_zlens",params.filename());
	if(!params.get("main_concentration",rscale)) error_message1("main_concentration",params.filename());
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
		delete[] htable;
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

void LensHaloNFW::initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, double my_slope, long* seed)
{
	LensHalo::initFromMassFunc(my_mass, my_Rmax, my_rscale, my_slope, seed);
    gmax = InterpolateFromTable(gtable,xmax);
}

const long LensHaloPseudoNFW::NTABLE = 10000;
const double LensHaloPseudoNFW::maxrm = 100.0;
int LensHaloPseudoNFW::count = 0;

double* LensHaloPseudoNFW::xtable = NULL;
double* LensHaloPseudoNFW::mhattable = NULL;

LensHaloPseudoNFW::LensHaloPseudoNFW()
: LensHalo()
{
}

LensHaloPseudoNFW::LensHaloPseudoNFW(InputParams& params)
{
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
	if (j==0) return mhat(y,beta);
	return (mhattable[j+1]-mhattable[j])/(xtable[j+1]-xtable[j])*(y-xtable[j]) + mhattable[j];
}

void LensHaloPseudoNFW::initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, double my_slope, long *seed){
	LensHalo::initFromMassFunc(my_mass,my_Rmax,my_rscale,my_slope,seed);
	beta = my_slope;
  xmax = my_Rmax/my_rscale;
	make_tables();
}

void LensHaloPseudoNFW::assignParams(InputParams& params){
	if(!params.get("main_mass",mass)) error_message1("main_mass",params.filename());
	if(!params.get("main_Rmax",Rmax)) error_message1("main_Rmax",params.filename());
	if(!params.get("main_zlens",zlens)) error_message1("main_zlens",params.filename());
	if(!params.get("main_concentration",rscale)) error_message1("main_concentration",params.filename());
	if(!params.get("main_slope",beta)) error_message1("main_slope",params.filename());
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
  beta = -2;
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
	if(!params.get("main_mass",mass)) error_message1("main_mass",params.filename());
	if(!params.get("main_Rmax",Rmax)) error_message1("main_Rmax",params.filename());
	if(!params.get("main_zlens",zlens)) error_message1("main_zlens",params.filename());
	if(!params.get("main_slope",beta)) error_message1("main_slope",params.filename());

	if(!params.get("main_stars_N",stars_N)) error_message1("main_stars_N",params.filename());
    else if(stars_N){
    	assignParams_stars(params);
    }
}

LensHaloPowerLaw::~LensHaloPowerLaw(){

}


LensHaloSimpleNSIE::LensHaloSimpleNSIE() : LensHalo(){
	sigma = 0.;
	zlens = 0.;
	fratio = 0.;
	pa = 0.;
	rcore = 0.;

}

LensHaloSimpleNSIE::LensHaloSimpleNSIE(InputParams& params){
	sigma = 0.;
	zlens = 0.;
	fratio = 0.;
	pa = 0.;
	rcore = 0.;

	assignParams(params);
}

void LensHaloSimpleNSIE::assignParams(InputParams& params){
	if(!params.get("main_mass",mass)) error_message1("main_mass",params.filename());
	if(!params.get("main_zlens",zlens)) error_message1("main_zlens",params.filename());

	if(!params.get("main_sigma",sigma)) error_message1("main_sigma",params.filename());
	if(!params.get("main_core",rcore)) error_message1("main_core",params.filename());
	if(!params.get("main_axis_ratio",fratio)) error_message1("main_axis_ratio",params.filename());
	if(!params.get("main_pos_angle",pa)) error_message1("main_pos_angle",params.filename());

	Rsize = rmaxNSIE(sigma,mass,fratio,rcore);
	Rmax = MAX(1.0,1.0/fratio)*Rsize;  // redefine

	assert(Rmax >= Rsize);

	if(!params.get("main_stars_N",stars_N)) error_message1("main_stars_N",params.filename());
    else if(stars_N){
    	assignParams_stars(params);
    }

}

LensHaloSimpleNSIE::~LensHaloSimpleNSIE(){

}
void LensHaloSimpleNSIE::initFromMass(float my_mass, long *seed){
	mass = my_mass;
	rcore = 0.0;
	sigma = 126*pow(mass/1.0e10,0.25); // From Tully-Fisher and Bell & de Jong 2001
  //std::cout << "Warning: All galaxies are spherical" << std::endl;
	fratio = (ran2(seed)+1)*0.5;  //TODO: Ben change this!  This is a kluge.
	pa = 2*pi*ran2(seed);  //TODO: This is a kluge.
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
		,bool kappa_off
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
		if(!kappa_off){
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
			if(!kappa_off){
				double tmp = -2.0*prefac/rcm2;

				gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
				gamma[1] += xcm[0]*xcm[1]*tmp;
			}
		}
	}

    // add stars for microlensing
    if(stars_N > 0 && stars_implanted){
   	 force_stars(alpha,kappa,gamma,xcm,kappa_off);
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
			if(!no_kappa){
        // TODO: this makes the kappa and gamma disagree with the alpha as calculated above
				KappaType tmp[2]={0,0};
        double xt[2]={0,0};
        float units = pow(sigma/lightspeed,2)/Grav/sqrt(fratio); // mass/distance(physical)
        xt[0]=xcm[0];
        xt[1]=xcm[1];
        
				*kappa += units*kappaNSIE(xt,fratio,rcore,pa);
				gammaNSIE(tmp,xt,fratio,rcore,pa);
				gamma[0] += units*tmp[0];
				gamma[1] += units*tmp[1];
			}

		}else{
			double xt[2]={0,0},tmp[2]={0,0};
			float units = pow(sigma/lightspeed,2)/Grav/sqrt(fratio); // mass/distance(physical)
			xt[0]=xcm[0];
			xt[1]=xcm[1];
			alphaNSIE(tmp,xt,fratio,rcore,pa);
      
			//alpha[0] = units*tmp[0];  // minus sign removed because already included in alphaNSIE
			//alpha[1] = units*tmp[1];  // Why was the "+=" removed?
			alpha[0] += units*tmp[0];
			alpha[1] += units*tmp[1];

			if(!no_kappa){
				KappaType tmp[2]={0,0};
				*kappa += units*kappaNSIE(xt,fratio,rcore,pa);
				gammaNSIE(tmp,xt,fratio,rcore,pa);
				gamma[0] += units*tmp[0];
				gamma[1] += units*tmp[1];
			}
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

    // add stars for microlensing
    if(stars_N > 0 && stars_implanted){
   	 force_stars(alpha,kappa,gamma,xcm,no_kappa);
    }

       return;
}

const long LensHaloHernquist::NTABLE = 10000;
const double LensHaloHernquist::maxrm = 100.0;
int LensHaloHernquist::count = 0;

double* LensHaloHernquist::xtable = NULL;
double* LensHaloHernquist::ftable = NULL;
double* LensHaloHernquist::gtable = NULL;
double* LensHaloHernquist::g2table = NULL;
double* LensHaloHernquist::htable = NULL;

LensHaloHernquist::LensHaloHernquist()
: LensHalo(), gmax(0)
{
	make_tables();
	gmax = InterpolateFromTable(gtable,xmax);
}

LensHaloHernquist::LensHaloHernquist(InputParams& params)
{
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
		htable = new double[NTABLE];
		g2table = new double[NTABLE];

		for(i = 0 ; i< NTABLE; i++){
			x = i*dx;
			xtable[i] = x;
			ftable[i] = ffunction(x);
			gtable[i] = gfunction(x);
			htable[i] = hfunction(x);
			g2table[i] = g2function(x);

		}
  }
  count++;
}

double LensHaloHernquist::InterpolateFromTable(double *table, double y){
	int j;
	j=(int)(y/maxrm*NTABLE);

	assert(y>=xtable[j] && y<=xtable[j+1]);
	if (j==0)
		{
		if (table==ftable) return ffunction(y);
		if (table==gtable) return gfunction(y);
		if (table==g2table) return g2function(y);
		if (table==htable) return hfunction(y);
		}
	return (table[j+1]-table[j])/(xtable[j+1]-xtable[j])*(y-xtable[j]) + table[j];
}

void LensHaloHernquist::assignParams(InputParams& params){
	if(!params.get("main_mass",mass)) error_message1("main_mass",params.filename());
	if(!params.get("main_Rmax",Rmax)) error_message1("main_Rmax",params.filename());
	if(!params.get("main_zlens",zlens)) error_message1("main_zlens",params.filename());
	if(!params.get("main_rscale",rscale)) error_message1("main_rscale",params.filename());
  xmax = Rmax/rscale;

	if(!params.get("main_stars_N",stars_N)) error_message1("main_stars_N",params.filename());
  else if(stars_N){
  	assignParams_stars(params);

  	std::cout << "Rmax " << Rmax <<std::endl;
  }

}

LensHaloHernquist::~LensHaloHernquist(){
	--count;
	if(count == 0){
		delete[] xtable;
		delete[] gtable;
		delete[] ftable;
		delete[] htable;
		delete[] g2table;
	}
}

const long LensHaloJaffe::NTABLE = 10000;
const double LensHaloJaffe::maxrm = 100.0;
int LensHaloJaffe::count = 0;

double* LensHaloJaffe::xtable = NULL;
double* LensHaloJaffe::ftable = NULL;
double* LensHaloJaffe::gtable = NULL;
double* LensHaloJaffe::g2table = NULL;
//double* LensHaloJaffe::htable = NULL;

LensHaloJaffe::LensHaloJaffe()
: LensHalo(), gmax(0)
{
	make_tables();
	gmax = InterpolateFromTable(gtable,xmax);
}

LensHaloJaffe::LensHaloJaffe(InputParams& params)
{
	assignParams(params);
	make_tables();
	gmax = InterpolateFromTable(gtable,xmax);
}

void LensHaloJaffe::make_tables(){
	if(count == 0){
		int i;
		double x, dx = maxrm/(double)NTABLE;

		xtable = new double[NTABLE];
		ftable = new double[NTABLE];
		gtable = new double[NTABLE];
		g2table = new double[NTABLE];
		//htable = new double[NTABLE];

		for(i = 0 ; i< NTABLE; i++){
			x = i*dx;
			xtable[i] = x;
			ftable[i] = ffunction(x);
			gtable[i] = gfunction(x);
			g2table[i] = g2function(x);
			//htable[i] = hfunction(x);
		}
  }
  count++;
}

double LensHaloJaffe::InterpolateFromTable(double *table, double y){
	int j;
	j=(int)(y/maxrm*NTABLE);

	assert(y>=xtable[j] && y<=xtable[j+1]);
	if (j==0)
		{
		if (table==ftable) return ffunction(y);
		if (table==gtable) return gfunction(y);
		if (table==g2table) return g2function(y);
//		if (table==htable) return hfunction(y);
		}
	return (table[j+1]-table[j])/(xtable[j+1]-xtable[j])*(y-xtable[j]) + table[j];
}

void LensHaloJaffe::assignParams(InputParams& params){
	if(!params.get("main_mass",mass)) error_message1("main_mass",params.filename());
	if(!params.get("main_Rmax",Rmax)) error_message1("main_Rmax",params.filename());
	if(!params.get("main_zlens",zlens)) error_message1("main_zlens",params.filename());
	if(!params.get("main_rscale",rscale)) error_message1("main_rscale",params.filename());
  xmax = Rmax/rscale;

	if(!params.get("main_stars_N",stars_N)) error_message1("main_stars_N",params.filename());
  else if(stars_N){
  	assignParams_stars(params);
  }

}

LensHaloJaffe::~LensHaloJaffe(){
	--count;
	if(count == 0){
		delete[] xtable;
		delete[] gtable;
		delete[] ftable;
		delete[] g2table;
		//delete[] htable;
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
  if(subtract_point){
    double x = sqrt(rcm2)/rscale;

    // can turn off kappa and gamma calculations to save times
    if(!no_kappa){
      *kappa += kappa_h(x)*prefac;

      tmp = (gamma_h(x) + 2.0*subtract_point)*prefac/rcm2;

      gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
      gamma[1] += xcm[0]*xcm[1]*tmp;
    }
  }
  // add stars for microlensing
  if(stars_N > 0 && stars_implanted){
 	 force_stars(alpha,kappa,gamma,xcm,no_kappa);
  }

}

void LensHaloDummy::assignParams(InputParams& params)
{
	if(!params.get("main_zlens",zlens)) error_message1("main_zlens",params.filename());
}

void LensHalo::serialize(RawData& d) const
{
	d << mass << Rmax << rscale << zlens << xmax;
	
	for(std::size_t i = 0; i < Nmod; ++i)
		d << mod[i];
	
	d << r_eps;
}

void LensHalo::unserialize(RawData& d)
{
	d >> mass >> Rmax >> rscale >> zlens >> xmax;
	
	for(std::size_t i = 0; i < Nmod; ++i)
		d >> mod[i];
	
	d >> r_eps;
}

void LensHalo::randomize(double step, long* seed)
{
	const std::type_info& type = typeid(*this);
	std::cerr << "Error: " << type.name() << "::randomize() not implemented!" << std::endl;
	exit(1);
}
