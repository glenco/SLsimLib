/*
 * lens_halos.cpp
 *
 *  Created on: 06.05.2013
 *      Author: mpetkova
 */

#include "slsimlib.h"

using namespace std;

LensHalo::LensHalo(){
	rscale = 1.0;
	mass = Rmax = xmax = posHalo[0] = posHalo[1] = 0.0;
	stars_implanted = false;
}

LensHalo::LensHalo(InputParams& params){
	assignParams(params);
	stars_implanted = false;
    posHalo[0] = posHalo[1] = 0.0;
}


void LensHalo::initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, PosType my_slope, long *seed){
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
		PosType *alpha     /// mass/Mpc
		,KappaType *kappa
		,KappaType *gamma
		,PosType *xcm     /// physical position on lens plane
		,bool no_kappa
		)
{
    PosType alpha_tmp[2];
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
const PosType LensHaloNFW::maxrm = 100.0;
int LensHaloNFW::count = 0;

PosType* LensHaloNFW::xtable = NULL;
PosType* LensHaloNFW::ftable = NULL;
PosType* LensHaloNFW::gtable = NULL;
PosType* LensHaloNFW::g2table = NULL;
PosType* LensHaloNFW::htable = NULL;

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
		PosType x, dx = maxrm/(PosType)NTABLE;

		xtable = new PosType[NTABLE];
		ftable = new PosType[NTABLE];
		gtable = new PosType[NTABLE];
		g2table = new PosType[NTABLE];
		htable = new PosType[NTABLE];

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

PosType LensHaloNFW::InterpolateFromTable(PosType *table, PosType y){
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

void LensHaloNFW::initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, PosType my_slope, long* seed)
{
	LensHalo::initFromMassFunc(my_mass, my_Rmax, my_rscale, my_slope, seed);
    gmax = InterpolateFromTable(gtable,xmax);
}

const long LensHaloPseudoNFW::NTABLE = 10000;
const PosType LensHaloPseudoNFW::maxrm = 100.0;
int LensHaloPseudoNFW::count = 0;

PosType* LensHaloPseudoNFW::xtable = NULL;
PosType* LensHaloPseudoNFW::mhattable = NULL;

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
PosType LensHaloPseudoNFW::mhat(PosType y, PosType beta){
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
		PosType x, dx = maxrm/(PosType)NTABLE;

		xtable = new PosType[NTABLE];
		mhattable = new PosType[NTABLE];

		for(i = 0 ; i< NTABLE; i++){
			x = i*dx;
			xtable[i] = x;
			mhattable[i] = mhat(x,beta);
		}

		count++;
	}
}

PosType LensHaloPseudoNFW::InterpolateFromTable(PosType y){
	int j;
	j=(int)(y/maxrm*NTABLE);

	assert(y>=xtable[j] && y<=xtable[j+1]);
	if (j==0) return mhat(y,beta);
	return (mhattable[j+1]-mhattable[j])/(xtable[j+1]-xtable[j])*(y-xtable[j]) + mhattable[j];
}

void LensHaloPseudoNFW::initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, PosType my_slope, long *seed){
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
  beta = -1;
  fratio = 1;
  rscale = xmax = 1.0;
}

LensHaloPowerLaw::LensHaloPowerLaw(InputParams& params){
	assignParams(params);
    calcModes(fratio, beta, pa, mod);
    for(int i=1;i<Nmod;i++){
        //std::cout << mod[i] << std::endl;
        if(mod[i]!=0){set_flag_elliptical(true);};
    }
    if(get_flag_elliptical()==false){beta=-beta;}; /// TODO the beta used for calculating kappa,gamma,alpha in the symmetric cases is a factor of -1 different from the asymmetric case
    std::cout << "if mods!=0 this must be 1: " << get_flag_elliptical() << std::endl;
    rscale = xmax = 1.0;
}

void LensHaloPowerLaw::initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, PosType my_slope, long *seed){
	LensHalo::initFromMassFunc(my_mass,my_Rmax,my_rscale,my_slope,seed);
	beta = my_slope;
    xmax = my_Rmax/my_rscale;
}

void LensHaloPowerLaw::assignParams(InputParams& params){
	if(!params.get("main_mass",mass)) error_message1("main_mass",params.filename());
	if(!params.get("main_Rmax",Rmax)) error_message1("main_Rmax",params.filename());
	if(!params.get("main_zlens",zlens)) error_message1("main_zlens",params.filename());
	if(!params.get("main_slope",beta)) error_message1("main_slope, example -1",params.filename());
    if(beta>=2.0) error_message1("main_slope < 2",params.filename());
    if(!params.get("main_axis_ratio",fratio)) error_message1("main_axis_ratio, example 1",params.filename());
    if(!params.get("main_pos_angle",pa)) error_message1("main_pos_angle",params.filename());


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
  else if(fratio > 1){
    ERROR_MESSAGE();
    std::cout << "parameter main_axis_ratio must be < 1 in file " << params.filename() << ". Use main_pos_angle to rotate the halo." << std::endl;
    exit(1);
  }

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

void LensHaloSimpleNSIE::initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, PosType my_slope, long *seed){
	initFromMass(my_mass,seed);
}

void LensHalo::force_halo(
	PosType *alpha     /// mass/Mpc
    ,KappaType *kappa
    ,KappaType *gamma
    ,double *xcm
    ,bool kappa_off
    ,bool subtract_point /// if true contribution from a point mass is subtracted
    ){
    //std::cout << "In lens_halo.cpp force_halo " << elliptical << std::endl;
    bool IsElliptical=get_flag_elliptical();
	if (IsElliptical==true){
        force_halo_asym(alpha,kappa,gamma,xcm,kappa_off,subtract_point);
    }else{
        force_halo_sym(alpha,kappa,gamma,xcm,kappa_off,subtract_point);
	}
}


void LensHalo::force_halo_sym(
		PosType *alpha     /// mass/Mpc
		,KappaType *kappa
		,KappaType *gamma
		,PosType *xcm
		,bool kappa_off
		,bool subtract_point /// if true contribution from a point mass is subtracted
		){

	PosType rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
	if(rcm2 < 1e-20) rcm2 = 1e-20;

	/// intersecting, subtract the point particle
	if(rcm2 < Rmax*Rmax){
		PosType prefac = mass/rcm2/pi;
		PosType x = sqrt(rcm2)/rscale;
		//PosType xmax = Rmax/rscale;

    //PosType tmp = (alpha_h(x,xmax) + 1.0*subtract_point)*prefac;
		PosType tmp = (alpha_h(x) + 1.0*subtract_point)*prefac;
		alpha[0] += tmp*xcm[0];
		alpha[1] += tmp*xcm[1];

		// can turn off kappa and gamma calculations to save times
		if(!kappa_off){
			*kappa += kappa_h(x)*prefac;
			//cout << x << endl;

			tmp = (gamma_h(x) + 2.0*subtract_point)*prefac/rcm2;

			gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
			gamma[1] += xcm[0]*xcm[1]*tmp;
		}

	}
	else
	{
		if (subtract_point == false)
		{
			PosType prefac = mass/rcm2/pi;
			alpha[0] += -1.0*prefac*xcm[0];
			alpha[1] += -1.0*prefac*xcm[1];

			// can turn off kappa and gamma calculations to save times
			if(!kappa_off){
				PosType tmp = -2.0*prefac/rcm2;

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

void LensHalo::force_halo_asym(
		PosType *alpha     /// mass/Mpc
		,KappaType *kappa
		,KappaType *gamma
		,double *xcm
		,bool kappa_off
		,bool subtract_point /// if true contribution from a point mass is subtracted
		){
	std::ofstream dfunc;
	dfunc.open( "dfunc.dat", ios::out | ios::app );
            
	double rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
    

	if(rcm2 < 1e-20) rcm2 = 1e-20;

	/// intersecting, subtract the point particle
	if(rcm2 < Rmax*Rmax){
		double x = sqrt(rcm2)/rscale;
		double theta;
        
        if (xcm[0]==0){
			xcm[0]=7.0710678e-21;
		}
        
		if (xcm[1]==0){
			xcm[1]=7.0710678e-21;
		}
        
        theta=atan2(xcm[1],xcm[0]);
        
		/*if (xcm[0]>0){
			theta=atan(xcm[1]/xcm[0]);
		}
		if (xcm[0]<0 && xcm[1]>=0){
			theta=atan(xcm[1]/xcm[0])+pi;
		}
		if (xcm[0]<0 && xcm[1]<0){
		   theta=atan(xcm[1]/xcm[0])-pi;
		}
		if (xcm[0]==0 && xcm[1]>0){
			theta=pi/2;
		}
		if (xcm[0]==0 && xcm[1]>0){
			theta=-1*pi/2;
		}*/
		//double prefac = mass/rcm2/pi;

		//double ellr2=rcm2*0.5*(cos(theta)*cos(theta)+(1./0.5/0.5)*sin(theta)*sin(theta))*pi;
		//double prefac = mass/ellr2;
		double prefac = mass/rscale/rscale/pi;

		//double xmax = Rmax/rscale;
		//double tmp = (alpha_h(x,xmax) + 1.0*subtract_point)*prefac;
        PosType alpha_tmp[2];
        alpha_asym(x,theta, alpha_tmp);
		double tmp =  1.0*subtract_point*prefac;
		alpha[0] +=  alpha_tmp[0]*prefac + tmp*xcm[0];
		alpha[1] +=  alpha_tmp[1]*prefac + tmp*xcm[1];

		// can turn off kappa and gamma calculations to save times
		if(!kappa_off){
			*kappa += kappa_asym(x,theta)*prefac;
			//if(kappa_asym(x,theta) < 0.01){
			//	dfunc << x << " " << theta << " " << kappa_asym(x,theta) << " " << xcm[0] << " " << xcm[1] << std::endl;
			//}
            //std::cout << x << std::endl;
            PosType gamma_tmp[2];
            gamma_asym(x,theta,gamma_tmp);
			tmp = (2.0*subtract_point)*prefac/rcm2;
			//gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*(tmp+gamma_tmp[0]*prefac/rcm2);
            gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp+gamma_tmp[0];
            gamma[1] += xcm[0]*xcm[1]*tmp+gamma_tmp[1];
			//gamma[1] += xcm[0]*xcm[1]*(tmp+gamma_tmp[0]*prefac/rcm2);
            

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



/*

void LensHaloSimpleNSIE::force_halo(
		PosType *alpha
		,KappaType *kappa
		,KappaType *gamma
		,PosType *xcm
		,bool no_kappa
		,bool subtract_point /// if true contribution from a point mass is subtracted
){

	PosType rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
	if(rcm2 < 1e-20) rcm2 = 1e-20;

  //**** test line

	if(rcm2 < Rmax*Rmax){
		PosType ellipR = ellipticRadiusNSIE(xcm,fratio,pa);
		if(ellipR > Rsize){
			// This is the case when the ray is within the NSIE's circular region of influence but outside its elliptical truncation

			PosType alpha_out[2],alpha_in[2],rin,x_in[2];
			PosType prefac = -1.0*mass/Rmax/pi;
			PosType r = sqrt(rcm2);
			float units = pow(sigma/lightspeed,2)/Grav/sqrt(fratio); // mass/distance(physical)

 			alpha_in[0] = alpha_in[1] = 0;
     
      rin = r*Rsize/ellipR;
     
			alpha_out[0] = prefac*xcm[0]/r;
			alpha_out[1] = prefac*xcm[1]/r;


			x_in[0] = rin*xcm[0]/r;
			x_in[1] = rin*xcm[1]/r;

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
        PosType xt[2]={0,0};
        float units = pow(sigma/lightspeed,2)/Grav/sqrt(fratio); // mass/distance(physical)
        xt[0]=xcm[0];
        xt[1]=xcm[1];
        
				*kappa += units*kappaNSIE(xt,fratio,rcore,pa);
				gammaNSIE(tmp,xt,fratio,rcore,pa);
				gamma[0] += units*tmp[0];
				gamma[1] += units*tmp[1];
			}

		}else{
			PosType xt[2]={0,0},tmp[2]={0,0};
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
			PosType prefac = mass/rcm2/pi;
			alpha[0] += -1.0*prefac*xcm[0];
			alpha[1] += -1.0*prefac*xcm[1];

			// can turn off kappa and gamma calculations to save times
			if(!no_kappa){
				PosType tmp = -2.0*prefac/rcm2;

				gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
				gamma[1] += xcm[0]*xcm[1]*tmp;
			}
		}
	}


	if(subtract_point){
		PosType fac = mass/rcm2/pi;
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
*/
void LensHaloSimpleNSIE::force_halo(
                                    PosType *alpha
                                    ,KappaType *kappa
                                    ,KappaType *gamma
                                    ,PosType *xcm
                                    ,bool no_kappa
                                    ,bool subtract_point /// if true contribution from a point mass is subtracted
                                    ){
  
	PosType rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
	if(rcm2 < 1e-20) rcm2 = 1e-20;
  
  //**** test line
  
	if(rcm2 < Rmax*Rmax){
		PosType ellipR = ellipticRadiusNSIE(xcm,fratio,pa);
		if(rcm2 > Rsize*Rsize){
			// This is the case when the ray is within the NSIE's circular region of influence but outside its elliptical truncation
      
			PosType alpha_iso[2],alpha_ellip[2];
			PosType prefac = -1.0*mass/Rmax/pi;
			PosType r = sqrt(rcm2);
			float units = pow(sigma/lightspeed,2)/Grav/sqrt(fratio); // mass/distance(physical)
      
      double f1 = (Rmax - r)/(Rmax - Rsize),f2 = (r - Rsize)/(Rmax - Rsize);
    
      // SIE solution
 			alpha_ellip[0] = alpha_ellip[1] = 0;
      alphaNSIE(alpha_ellip,xcm,fratio,rcore,pa);
      alpha_ellip[0] *= units;
      alpha_ellip[1] *= units;
      
      /*/ SIS solution
      alpha_iso[0] = alpha_iso[1] = 0;
      alphaNSIE(alpha_iso,xcm,1,rcore,pa);
      alpha_iso[0] *= units;
      alpha_iso[1] *= units;
      //*/
      
      // point mass solution
      //PosType tmp = mass/rcm2/pi;
      PosType tmp = mass/Rmax/pi/r;
			alpha_iso[0] = -1.0*tmp*xcm[0];
			alpha_iso[1] = -1.0*tmp*xcm[1];

			alpha[0] += alpha_iso[0]*f2 + alpha_ellip[0]*f1;
			alpha[1] += alpha_iso[1]*f2 + alpha_ellip[1]*f1;
      
			// can turn off kappa and gamma calculations to save times
			if(!no_kappa){
				PosType tmp = -2.0*prefac/rcm2;
        
				gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp*f2;
				gamma[1] += xcm[0]*xcm[1]*tmp*f2;

        KappaType tmp_k[2]={0,0};
				*kappa += units*kappaNSIE(xcm,fratio,rcore,pa)*f1;
				gammaNSIE(tmp_k,xcm,fratio,rcore,pa);
				gamma[0] += units*tmp_k[0]*f1;
				gamma[1] += units*tmp_k[1]*f1;
			}
      
		}else{
			PosType xt[2]={0,0},tmp[2]={0,0};
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

    if(subtract_point){
      PosType fac = mass/rcm2/pi;
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
	else
	{  // outside of the halo
		if (subtract_point == false)
		{
			PosType prefac = mass/rcm2/pi;
			alpha[0] += -1.0*prefac*xcm[0];
			alpha[1] += -1.0*prefac*xcm[1];
      
			// can turn off kappa and gamma calculations to save times
			if(!no_kappa){
				PosType tmp = -2.0*prefac/rcm2;
        
				gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
				gamma[1] += xcm[0]*xcm[1]*tmp;
			}
		}
	}
  
  
  // add stars for microlensing
  if(stars_N > 0 && stars_implanted){
    force_stars(alpha,kappa,gamma,xcm,no_kappa);
  }
  
  return;
}

const long LensHaloHernquist::NTABLE = 10000;
const PosType LensHaloHernquist::maxrm = 100.0;
int LensHaloHernquist::count = 0;

PosType* LensHaloHernquist::xtable = NULL;
PosType* LensHaloHernquist::ftable = NULL;
PosType* LensHaloHernquist::gtable = NULL;
PosType* LensHaloHernquist::g2table = NULL;
PosType* LensHaloHernquist::htable = NULL;

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
		PosType x, dx = maxrm/(PosType)NTABLE;

		xtable = new PosType[NTABLE];
		ftable = new PosType[NTABLE];
		gtable = new PosType[NTABLE];
		htable = new PosType[NTABLE];
		g2table = new PosType[NTABLE];

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

PosType LensHaloHernquist::InterpolateFromTable(PosType *table, PosType y){
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
const PosType LensHaloJaffe::maxrm = 100.0;
int LensHaloJaffe::count = 0;

PosType* LensHaloJaffe::xtable = NULL;
PosType* LensHaloJaffe::ftable = NULL;
PosType* LensHaloJaffe::gtable = NULL;
PosType* LensHaloJaffe::g2table = NULL;
//PosType* LensHaloJaffe::htable = NULL;

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
		PosType x, dx = maxrm/(PosType)NTABLE;

		xtable = new PosType[NTABLE];
		ftable = new PosType[NTABLE];
		gtable = new PosType[NTABLE];
		g2table = new PosType[NTABLE];
		//htable = new PosType[NTABLE];

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

PosType LensHaloJaffe::InterpolateFromTable(PosType *table, PosType y){
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

void LensHaloDummy::initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, PosType my_slope, long *seed){
	mass = 1.e-10;
	Rmax = my_Rmax;
	rscale = my_rscale;
  xmax = Rmax/rscale;
}


void LensHaloDummy::force_halo(PosType *alpha,KappaType *kappa,KappaType *gamma,PosType *xcm,bool no_kappa,bool subtract_point)
{
	PosType rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
	PosType prefac = mass/rcm2/pi;
	PosType tmp = subtract_point*prefac;
	alpha[0] += tmp*xcm[0];
	alpha[1] += tmp*xcm[1];
  if(subtract_point){
    PosType x = sqrt(rcm2)/rscale;

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

std::size_t LensHalo::Nparams() const
{
	return 0;
}

PosType LensHalo::getParam(std::size_t p) const
{
	switch(p)
	{
		default:
			throw std::invalid_argument("bad parameter index for getParam()");
	}
}

PosType LensHalo::setParam(std::size_t p, PosType val)
{
	switch(p)
	{
		default:
			throw std::invalid_argument("bad parameter index for setParam()");
	}
}

void LensHalo::printCSV(std::ostream&, bool header) const
{
	const std::type_info& type = typeid(*this);
	std::cerr << "LensHalo subclass " << type.name() << " does not implement printCSV()" << std::endl;
	std::exit(1);
}
