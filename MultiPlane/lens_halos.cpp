/*
 * lens_halos.cpp
 *
 *  Created on: 06.05.2013
 *      Author: mpetkova
 */

#include "slsimlib.h"

using namespace std;

/// Shell constructor
LensHalo::LensHalo(){
	rscale = 1.0;
	mass = Rmax = xmax = posHalo[0] = posHalo[1] = 0.0;
	stars_implanted = false;
  elliptical_flag = false;
}

LensHalo::LensHalo(InputParams& params){
	assignParams(params);
	stars_implanted = false;
  posHalo[0] = posHalo[1] = 0.0;
  elliptical_flag = false;
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

/// calculates the deflection etc. caused by stars alone
void LensHalo::force_stars(
		PosType *alpha     /// mass/Mpc
		,KappaType *kappa
		,KappaType *gamma
		,PosType const *xcm     /// physical position on lens plane
		)
{
    PosType alpha_tmp[2];
    KappaType gamma_tmp[3], tmp = 0;
    KappaType phi;

    gamma_tmp[0] = gamma_tmp[1] = gamma_tmp[2] = 0.0;
    alpha_tmp[0] = alpha_tmp[1] = 0.0;

 	 substract_stars_disks(xcm,alpha,kappa,gamma);

	 // do stars with tree code
    star_tree->force2D_recur(xcm,alpha_tmp,&tmp,gamma_tmp,&phi);

    alpha[0] -= star_massscale*alpha_tmp[0];
    alpha[1] -= star_massscale*alpha_tmp[1];

    {
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
PosType* LensHaloNFW::xgtable = NULL;
PosType*** LensHaloNFW::modtable= NULL; // was used for Ansatz IV



LensHaloNFW::LensHaloNFW()
: LensHalo(), gmax(0)
{
	make_tables();
	gmax = InterpolateFromTable(gtable, xmax);
}

LensHaloNFW::LensHaloNFW(float my_mass,float my_Rmax,PosType my_zlens,float my_concentration,float my_fratio,float my_pa,int my_stars_N){
    mass=my_mass, Rmax=my_Rmax, zlens=my_zlens, rscale=my_concentration;
    fratio=my_fratio, pa=my_pa, stars_N=my_stars_N;
    stars_implanted = false;

	rscale = Rmax/rscale; // TODO make use of rscale/concentration in NFW clearer
  xmax = Rmax/rscale;
    
  make_tables();
	gmax = InterpolateFromTable(gtable, xmax);
    
  set_slope(1);
  /// If the axis ratio given in the parameter file is set to 1 all ellipticizing routines are skipped.
  if(fratio!=1){
    std::cout << "NFW constructor: slope set to " << get_slope() << std::endl;
    calcModes(fratio, get_slope(), pa, mod); // to ellipticize potential instead of  kappa take calcModes(fratio, 2-get_slope(), pa, mod);
    for(int i=1;i<Nmod;i++){
      if(mod[i]!=0){set_flag_elliptical(true);};
    }

  }else set_flag_elliptical(false);
    
    // std::cout << mass << " " << rscale << std::endl;
    
 }

/* LensHalo::LensHalo(mass,Rmax,zlens, // base
 rscale,fratio,pa,stars_N, //NFW,Hernquist, Jaffe
 rscale,fratio,pa,beta // Pseudo NFW
 rscale,fratio,pa,sigma,rcore // NSIE
 zlens,stars_N // dummy
 ){
 
 stars_implanted = false;
 posHalo[0] = posHalo[1] = 0.0;
 }*/

LensHaloNFW::LensHaloNFW(InputParams& params)
{
	assignParams(params);
	make_tables();
	gmax = InterpolateFromTable(gtable, xmax);

    mnorm = renormalization(get_Rmax());
    std::cout << "mass normalization: " << mnorm << std::endl;
    
    // If the 2nd argument in calcModes(fratio, slope, pa, mod), the slope, is set to 1 it yields an elliptical kappa contour of given axis ratio (fratio) at the radius where the slope of the 3D density profile is -2, which is defined as the scale radius for the NFW profile. To ellipticize the potential instead of the convergence use calcModes(fratio, 2-get_slope(), pa, mod), this produces also an ellipse in the convergence map, but at the radius where the slope is 2-get_slope().
    set_slope(1);
    // If the axis ratio given in the parameter file is set to 1 all ellipticizing routines are skipped.
    if(fratio!=1){
        std::cout << "NFW constructor: slope set to " << get_slope() << std::endl;
        //for(int i=1;i<20;i++){
         //  calcModes(fratio, 0.1*i, pa, mod);
        //}
        //calcModes(fratio, get_slope()-0.5, pa, mod); // to ellipticize potential instead of  kappa take calcModes(fratio, 2-get_slope(), pa, mod);
        calcModes(fratio, get_slope(), pa, mod); // to ellipticize potential instead of  kappa take calcModes(fratio, 2-get_slope(), pa, mod);
        //calcModes(fratio, get_slope()+0.5, pa, mod2); // to ellipticize potential instead of  kappa take calcModes(fratio, 2-get_slope(), pa, mod);
        
        for(int i=1;i<Nmod;i++){
            if(mod[i]!=0){set_flag_elliptical(true);};
        }
        
        //std::cout << "if mods!=0 this must be 1: " << get_flag_elliptical() << std::endl;
    }
    std::cout << get_flag_elliptical() << std::endl;
}

void LensHaloNFW::make_tables(){
	if(count == 0){
		int i;
        //struct Ig_func g(*this);
		PosType x, dx = maxrm/(PosType)NTABLE;

		xtable = new PosType[NTABLE];
		ftable = new PosType[NTABLE];
		gtable = new PosType[NTABLE];
		g2table = new PosType[NTABLE];
		htable = new PosType[NTABLE];
        xgtable = new PosType[NTABLE];
        
        
		for(i = 0 ; i< NTABLE; i++){
			x = i*dx;
			xtable[i] = x;
			ftable[i] = ffunction(x);
			gtable[i] = gfunction(x);
			g2table[i] = g2function(x);
			htable[i] = hfunction(x);
            if(i==0){xgtable[i]=0;}
            if(i!=0){
                xgtable[i] = alpha_int(x);
                //Utilities::nintegrate<Ig_func>(g,1E-4,x,dx/10.);
            }
        }
        
        // modtable[axis ratio 100][potential slope beta 1000][Nmods 32] for Ansatz IV
        
        int j;
        modtable = new PosType**[100];
 
         
        for(i = 0; i < 100; i++){
            modtable[i] = new PosType*[200];
            for(j = 0; j< 200; j++){
                modtable[i][j] = new PosType[Nmod];
            }
        }
        /*
         
        for(i = 0; i<99; i++){
            std::cout<< i << std::endl;
            PosType iq=0.01*(i+1);
            for(j = 0; j< 200; j++){
                PosType beta_r=0.01*(j+1);
                calcModesC(beta_r, iq, pa, mod);
                for(k=0;k<Nmod;k++){
                    modtable[i][j][k]=mod[k];
                }
            }
        }
        */
        
  }
  count++;
}

// InterpolateModes was used for Ansatz IV and is an efficient way to calculate the Fourier modes used for elliptisizing the isotropic profiles before the program starts

PosType LensHaloNFW::InterpolateModes(int whichmod, PosType q, PosType b){
    PosType x1,x2,y1,y2,f11,f12,f21,f22;
    int i,j,k;
    int NTABLEB=200;
    int NTABLEQ=99;
    PosType const maxb=2.0;
    PosType const maxq=0.99;
    k=whichmod;
    j=(int)(b/maxb*NTABLEB);
    i=(int)(q/maxq*NTABLEQ);
    f11=modtable[i][j][k];
    f12=modtable[i][j+1][k];
    f21=modtable[i+1][j][k];
    f22=modtable[i+1][j+1][k];
    x1=i*maxq/NTABLEQ;
    x2=(i+1)*maxq/NTABLEQ;
    y1=j*maxb/NTABLEB;
    y2=(j+1)*maxb/NTABLEB;
    //std::cout << "x12y12: " << q << " " << x1 << " " << x2 << " "<< b << " " << y1 << " " << y2 << " " << std::endl;
    //std::cout << "IM: " << f11 << " " << f12 << " " << f21 << " " << f22 << " " << res << std::endl;
    return 1.0/(x2-x1)/(y2-y1)*(f11*(x2-q)*(y2-b)+f21*(q-x1)*(y2-b)+f12*(x2-q)*(b-y1)+f22*(q-x1)*(b-y1));
}


PosType LensHaloNFW::InterpolateFromTable(PosType *table, PosType y) const{
	int j;
	j=(int)(y/maxrm*NTABLE);

	assert(y>=xtable[j] && y<=xtable[j+1]);
	if (j==0)
		{
		if (table==ftable) return ffunction(y);
		if (table==gtable) return gfunction(y);
		if (table==g2table) return g2function(y);
		if (table==htable) return hfunction(y);
    if (table==xgtable) return alpha_int(y);
		}
	return (table[j+1]-table[j])/(xtable[j+1]-xtable[j])*(y-xtable[j]) + table[j];
}


void LensHaloNFW::assignParams(InputParams& params){
	if(!params.get("main_mass",mass)) error_message1("main_mass",params.filename());
	if(!params.get("main_Rmax",Rmax)) error_message1("main_Rmax",params.filename());
	if(!params.get("main_zlens",zlens)) error_message1("main_zlens",params.filename());
	if(!params.get("main_concentration",rscale)) error_message1("main_concentration",params.filename());
    if(!params.get("main_axis_ratio",fratio)){fratio=1; std::cout << "main_axis_ratio not defined in file " << params.filename() << ", hence set to 1." << std::endl;};
    if(!params.get("main_pos_angle",pa)){pa=0; std::cout << "main_pos_angle not defined in file " << params.filename() << ", hence set to 0." << std::endl;};
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
        delete[] xgtable;
        
        // was used for Ansatz IV
        
        for(int i=0; i<99; i++){
             for(int j=0; j<200; j++){
                 delete[] modtable[i][j];
             }
            delete[] modtable[i];
        }
        delete[] modtable;
        
        
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
/// constructor
LensHaloPseudoNFW::LensHaloPseudoNFW(
                                     float my_mass             /// mass in solar masses
                                     ,float my_Rmax            /// maximum radius in Mpc
                                     ,PosType my_zlens         /// redshift
                                     ,float my_concentration   /// Rmax/rscale
                                     ,PosType my_beta          /// large r slope, see class description
                                     ,float my_fratio          /// axis ratio
                                     ,float my_pa              /// position angle
                                     ,int my_stars_N           /// number of stars, not yet implanted
                                     )
  {
    mass = my_mass;
    Rmax = my_Rmax;
    zlens = my_zlens;
    beta = my_beta;
    fratio = my_fratio;
    pa = my_pa;
    stars_N = my_stars_N;
  
    stars_implanted = false;
	  rscale = Rmax/my_concentration;
    xmax = Rmax/rscale;
    
    make_tables();
    if(fratio!=1){
        std::cout << "Note: Fourier modes set to ellipticize kappa at slope main_slope+0.5, i.e. "<< get_slope()+0.5 << std::endl;
        calcModes(fratio, get_slope()+0.5, pa, mod);
        for(int i=1;i<Nmod;i++){
            if(mod[i]!=0){set_flag_elliptical(true);};
        }
    }

}

/// The Fourier modes set to ellipticize kappa at slope main_slope+0.5, i.e. e.g. 1.5 for main_slope = 1. Note that set_slope is overridden for PseudoNFW to recalculate tables for different beta. But only fixed values of beta, i.e. 1,2 and >=3 are allowed!
LensHaloPseudoNFW::LensHaloPseudoNFW(InputParams& params)
{
	assignParams(params);
	make_tables();
    if(fratio!=1){
        std::cout << "Note: Fourier modes set to ellipticize kappa at slope main_slope+0.5, i.e. "<< get_slope()+0.5 << std::endl;
        calcModes(fratio, get_slope()+0.5, pa, mod);
        for(int i=1;i<Nmod;i++){
            //std::cout << mod[i] << std::endl;
            if(mod[i]!=0){set_flag_elliptical(true);};
        }
    }

}

/// Auxiliary function for PseudoNFW profile
// previously defined in tables.cpp
PosType LensHaloPseudoNFW::mhat(PosType y, PosType beta) const{
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

PosType LensHaloPseudoNFW::gfunction(PosType y) const{
    int j;
    j=(int)(y/maxrm*NTABLE);
	assert(y>=xtable[j] && y<=xtable[j+1]);
	if (j==0) return mhat(y,beta);
    return ((mhattable[j+1]-mhattable[j])/(xtable[j+1]-xtable[j])*(y-xtable[j]) + mhattable[j]);
}

PosType LensHaloPseudoNFW::InterpolateFromTable(PosType y) const{
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
    if(!params.get("main_axis_ratio",fratio)){fratio=1; std::cout << "main_axis_ratio not defined in file " << params.filename() << ", hence set to 1." << std::endl;};
    if(!params.get("main_pos_angle",pa)){pa=0; std::cout << "main_pos_angle not defined in file " << params.filename() << ", hence set to 0." << std::endl;};
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
  beta = 1;
  fratio = 1;
  rscale = xmax = 1.0;
}

/// constructor
LensHaloPowerLaw::LensHaloPowerLaw(
                                   float my_mass       /// mass of halo in solar masses
                                   ,float my_Rmax      /// maximum radius of halo in Mpc
                                   ,PosType my_zlens   /// redshift of halo
                                   ,float my_rscale    /// dummy scale,  not acutally used
                                   ,PosType my_beta    /// logarithmic slop of surface density, kappa \propto r^{-beta}
                                   ,float my_fratio    /// axis ratio in asymetric case
                                   ,float my_pa        /// position angle
                                   ,int my_stars_N     /// number of stars, not yet implanted
                                   ){
    mass=my_mass, Rmax=my_Rmax, zlens=my_zlens, rscale=my_rscale;
    beta=my_beta;
    fratio=my_fratio, pa=my_pa, stars_N=my_stars_N;
    stars_implanted = false;
    if(fratio!=1){
        calcModes(fratio, beta, pa, mod);
        for(int i=1;i<Nmod;i++){
            //std::cout << i << " " << mod[i] << std::endl;
            if(mod[i]!=0){set_flag_elliptical(true);};
        }
    }
    // rscale = xmax = 1.0;
    // Fabien : replacing it by :
    rscale = 1;
    xmax = Rmax/rscale ;
}

LensHaloPowerLaw::LensHaloPowerLaw(InputParams& params){
	assignParams(params);
    /// If the 2nd argument in calcModes(fratio, slope, pa, mod), the slope, is set to 1 it yields an elliptical kappa contour of given axis ratio (fratio) at the radius where the slope of the 3D density profile is -2, which is defined as the scale radius for the NFW profile. To ellipticize the potential instead of the convergence use calcModes(fratio, 2-get_slope(), pa, mod), this produces also an ellipse in the convergence map, but at the radius where the slope is 2-get_slope().
    /// If the axis ratio given in the parameter file is set to 1 all ellipticizing routines are skipped.
    if(fratio!=1){
        //for(int islope=1;islope<20;islope++){
        //beta=islope*0.1;
        /*
            for(int islope=1;islope<20;islope++){
                for(int ifratio=1;ifratio<10;ifratio++){
                    calcModes(ifratio*0.1, islope*0.1, pa, mod);
                    for(int i=4;i<Nmod;i=i+4){
                        std::cout << i << " " << islope*0.1 << " " << ifratio*0.1 << " "<< mod[i] << " " << modfunc(i, islope*0.1, ifratio*0.1)<< std::endl;
                    }
                }
            }
    
        calcModes(fratio, beta, pa, mod);
        
        for(int i=0;i<Nmod;i++){
            std::cout<< "before: " << mod[i] << std::endl;
            mod[i]=modfunc(i, beta, fratio);
            std::cout<< "after: " << mod[i] << std::endl;
        }

        */
        calcModes(fratio, beta, pa, mod);
        //    std::cout << mod[4] << " " << modfunc(4, 1, 0.5) << std::endl;
            for(int i=1;i<Nmod;i++){
                //std::cout << i << " " << mod[i] << std::endl;
                if(mod[i]!=0){set_flag_elliptical(true);};
            }
        //}

    }else elliptical_flag = false;

  // rscale = xmax = 1.0;
  // mnorm = renormalization(get_Rmax());
  mnorm = 1. ;
  std::cout << "mass normalization: " << mnorm << std::endl;

    // rscale = xmax = 1.0; // Commented by Fabien in order to have a correct computation of the potential term in the time delay.
    // Fabien : replacing it by :
    rscale = 1;
    xmax = Rmax/rscale;
    // Be careful : the other constructors have not been changed !
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
    //if(beta>=2.0) error_message1("main_slope < 2",params.filename());
    if(!params.get("main_axis_ratio",fratio)){fratio=1; std::cout << "main_axis_ratio not defined in file " << params.filename() << ", hence set to 1." << std::endl;};
    if(!params.get("main_pos_angle",pa)){pa=0.0; std::cout << "main_pos_angle not defined in file " << params.filename() << ", hence set to 0." << std::endl;};


	if(!params.get("main_stars_N",stars_N)) error_message1("main_stars_N",params.filename());
    else if(stars_N){
    	assignParams_stars(params);
    }
}

LensHaloPowerLaw::~LensHaloPowerLaw(){

}


/*
LensHaloSimpleNSIE::LensHaloSimpleNSIE(float my_mass,float my_Rmax,PosType my_zlens,float my_rscale,float my_sigma
                                      , float my_rcore,float my_fratio,float my_pa,int my_stars_N){
    mass=my_mass, Rmax=my_Rmax, zlens=my_zlens, rscale=my_rscale;
    sigma=my_sigma, rcore=my_rcore;
    fratio=my_fratio, pa=my_pa, stars_N=my_stars_N;
    stars_implanted = false;
	Rsize = rmaxNSIE(sigma,mass,fratio,rcore);
	Rmax = MAX(1.0,1.0/fratio)*Rsize;  // redefine
	assert(Rmax >= Rsize);
}
*/
LensHaloSimpleNSIE::LensHaloSimpleNSIE(float my_mass,PosType my_zlens,float my_sigma
                                       , float my_rcore,float my_fratio,float my_pa,int my_stars_N):LensHalo(){
  mass=my_mass, zlens=my_zlens, rscale=1.0;
  
  sigma=my_sigma, rcore=my_rcore;
  fratio=my_fratio, pa=my_pa, stars_N=my_stars_N;
  stars_implanted = false;
	Rsize = rmaxNSIE(sigma,mass,fratio,rcore);
	Rmax = MAX(1.0,1.0/fratio)*Rsize;  // redefine
	assert(Rmax >= Rsize);
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
/*
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
*/

/** \brief returns the lensing quantaties of a ray in center of mass coordinates.
 *
 *  Warning: This adds to input value of alpha, kappa, gamma, and phi.  They need 
 *  to be veroed out if the contribution of just this halo is wanted.
 */
void LensHalo::force_halo(
    PosType *alpha          /// deflection solar mass/Mpc
    ,KappaType *kappa     /// surface density in Msun/Mpc^2 (?)
    ,KappaType *gamma     /// three components of shear
    ,KappaType *phi       /// potential in solar masses
    ,PosType const *xcm   /// position relative to center (Mpc?)
    ,bool subtract_point /// if true contribution from a point mass is subtracted
    )
{
	if (elliptical_flag){
        force_halo_asym(alpha,kappa,gamma,xcm,subtract_point);
        //assert(!isinf(*kappa) );
    }else{
        force_halo_sym(alpha,kappa,gamma,phi,xcm,subtract_point);
        //assert(!isinf(*kappa) );
	}
}

void LensHalo::force_halo_sym(
		PosType *alpha     /// solar mass/Mpc
		,KappaType *kappa  /// convergence
		,KappaType *gamma  /// three components of shear
    ,KappaType *phi      /// potential solar masses
		,PosType const *xcm     /// position relative to center (Mpc)
		,bool subtract_point /// if true contribution from a point mass is subtracted
		)
{
            
	PosType rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
            
	if(rcm2 < 1e-20) rcm2 = 1e-20;
            
	/// intersecting, subtract the point particle
	if(rcm2 < Rmax*Rmax)
  {
    PosType prefac = mass/rcm2/pi;
		PosType x = sqrt(rcm2)/rscale;
		// PosType xmax = Rmax/rscale;
    
    PosType tmp = (alpha_h(x) + 1.0*subtract_point)*prefac;
		alpha[0] += tmp*xcm[0];
		alpha[1] += tmp*xcm[1];

    *kappa += kappa_h(x)*prefac;
    
    tmp = (gamma_h(x) + 2.0*subtract_point) * prefac / rcm2;
		gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
		gamma[1] += xcm[0]*xcm[1]*tmp;
    
    *phi += phi_h(x) * mass / pi ;
    // *phi += phi_h(x,Rmax) * mass / pi ;
	}
	else // the point particle is not subtracted
	{
		if (subtract_point == false)
		{
			PosType prefac = mass/rcm2/pi;
			alpha[0] += -1.0 * prefac * xcm[0];
			alpha[1] += -1.0 * prefac * xcm[1];
      
      //std::cout << "rcm2  = " << rcm2 << std::endl;
      //std::cout << "prefac  = " << prefac << std::endl;
      //std::cout << "xcm  = " << xcm[0] << " " << xcm[1] << std::endl;

			PosType tmp = -2.0*prefac/rcm2;
                
      // kappa is equal to 0 in the point mass case.

			gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
      gamma[1] += xcm[0]*xcm[1]*tmp;
      
      *phi += 0.5 * log(rcm2) * mass / pi ;
		}
	}

    /// add stars for microlensing
    if(stars_N > 0 && stars_implanted)
    {
   	 force_stars(alpha,kappa,gamma,xcm);
    }

	return;
}
// TODO: put in some comments about the units used
void LensHalo::force_halo_asym(
		PosType *alpha     /// mass/Mpc
		,KappaType *kappa
		,KappaType *gamma
		,PosType const *xcm
		,bool subtract_point /// if true contribution from a point mass is subtracted
		){
	//std::ofstream dfunc;
	//dfunc.open( "dfunc.dat", ios::out | ios::app );
            
	double rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
    

	if(rcm2 < 1e-20) rcm2 = 1e-20;

	/// intersecting, subtract the point particle
	if(rcm2 < Rmax*Rmax)
  {
		double r = sqrt(rcm2);///rscale;
        //std::cout << sqrt(rcm2) << " " << rscale << " " << Rmax << std::endl;
		double theta;
        
    if(xcm[0] == 0.0 && xcm[1] == 0.0) theta = 0.0;
    else theta=atan2(xcm[1],xcm[0]);

		// double xmax = Rmax/rscale;
    PosType alpha_tmp[2],kappa_tmp,gamma_tmp[2];
        
    alphakappagamma_asym(r,theta, alpha_tmp,&kappa_tmp,gamma_tmp);
    
		//alpha[0] +=  alpha_tmp[0]*prefac*xcm[0] + tmp*xcm[0];
    //alpha[1] +=  alpha_tmp[1]*prefac*xcm[1] + tmp*xcm[1];
    
		alpha[0] +=  alpha_tmp[0];
    alpha[1] +=  alpha_tmp[1];

    *kappa += kappa_tmp;
    gamma[0] += gamma_tmp[0];
    gamma[1] += gamma_tmp[1];

    if(subtract_point){
      double tmp =  mass/pi/rcm2;
      alpha[0] +=  tmp*xcm[0];
      alpha[1] +=  tmp*xcm[1];

      tmp = 2.0*tmp/rcm2;
      gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
      gamma[1] += xcm[0]*xcm[1]*tmp;
    }

    /*
		*kappa += kappa_asym(r,theta);

    //dfunc << x << " " << theta << " " << kappa_asym(x,theta) << " " << bfunction(x) << " " << xcm[0] << " " << xcm[1] << std::endl;

    //std::cout << x << " " << rscale << " " << xmax << std::endl;
    //PosType gamma_tmp[2];
    //gamma_asym(r,theta,gamma_tmp);
		double tmp = (2.0*subtract_point)*prefac/rcm2;
		//gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*(tmp+gamma_tmp[0]*prefac/rcm2);
    //std::cout << prefac << std::endl;
    gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp + gamma_tmp[0]*prefac/rcm2;
    gamma[1] += xcm[0]*xcm[1]*tmp + gamma_tmp[1]*prefac/rcm2;
		//gamma[1] += xcm[0]*xcm[1]*(tmp+gamma_tmp[0]*prefac/rcm2);
     */
  }
	else
	{
		if (subtract_point == false)
		{
			double prefac = mass/rcm2/pi;
			alpha[0] += -1.0*prefac*xcm[0];
			alpha[1] += -1.0*prefac*xcm[1];

      double tmp = -2.0*prefac/rcm2;
      gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
      gamma[1] += xcm[0]*xcm[1]*tmp;

		}
	}

    // add stars for microlensing
    if(stars_N > 0 && stars_implanted){
   	 force_stars(alpha,kappa,gamma,xcm);
    }


	return;
}



/*

void LensHaloSimpleNSIE::force_halo(
		PosType *alpha
		,KappaType *kappa
		,KappaType *gamma
		,PosType *xcm
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
      
			{
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

			{
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
			{
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
		{
			fac = 2.0*fac/rcm2;

			gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*fac;
			gamma[1] += xcm[0]*xcm[1]*fac;
		}

	}

    // add stars for microlensing
    if(stars_N > 0 && stars_implanted){
   	 force_stars(alpha,kappa,gamma,xcm);
    }

       return;
}
*/
void LensHaloSimpleNSIE::force_halo(
                                    PosType *alpha
                                    ,KappaType *kappa
                                    ,KappaType *gamma
                                    ,KappaType *phi
                                    ,PosType const *xcm
                                    ,bool subtract_point /// if true contribution from a point mass is subtracted
                                    )
{
  
	PosType rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
	if(rcm2 < 1e-20) rcm2 = 1e-20;
  
	if(rcm2 < Rmax*Rmax)
    {
		PosType ellipR = ellipticRadiusNSIE(xcm,fratio,pa);
        
		if(rcm2 > Rsize*Rsize)
        {
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
        // PosType tmp = mass/rcm2/pi;
        PosType tmp = mass/Rmax/pi/r;
        alpha_iso[0] = -1.0*tmp*xcm[0];
        alpha_iso[1] = -1.0*tmp*xcm[1];

        alpha[0] += alpha_iso[0]*f2 + alpha_ellip[0]*f1;
        alpha[1] += alpha_iso[1]*f2 + alpha_ellip[1]*f1;
      
            // can turn off kappa and gamma calculations to save times
            {
            PosType tmp = -2.0*prefac/rcm2;
        
            gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp*f2;
            gamma[1] += xcm[0]*xcm[1]*tmp*f2;

            KappaType tmp_k[2]={0,0};
            *kappa += units*kappaNSIE(xcm,fratio,rcore,pa)*f1;
            gammaNSIE(tmp_k,xcm,fratio,rcore,pa);
            gamma[0] += units*tmp_k[0]*f1;
            gamma[1] += units*tmp_k[1]*f1;
            }
      
        }
        else
        {
        PosType xt[2]={0,0},tmp[2]={0,0};
        float units = pow(sigma/lightspeed,2)/Grav/sqrt(fratio); // mass/distance(physical)
        xt[0]=xcm[0];
        xt[1]=xcm[1];
        alphaNSIE(tmp,xt,fratio,rcore,pa);
      
        //alpha[0] = units*tmp[0];  // minus sign removed because already included in alphaNSIE
        //alpha[1] = units*tmp[1];  // Why was the "+=" removed?
        alpha[0] += units*tmp[0];
        alpha[1] += units*tmp[1];
      
            {
            KappaType tmp[2]={0,0};
            *kappa += units*kappaNSIE(xt,fratio,rcore,pa);
            gammaNSIE(tmp,xt,fratio,rcore,pa);
            gamma[0] += units*tmp[0];
            gamma[1] += units*tmp[1];
			}
		}

        
        if(subtract_point)
        {
            PosType fac = mass/rcm2/pi;
            alpha[0] += fac*xcm[0];
            alpha[1] += fac*xcm[1];
      
            // can turn off kappa and gamma calculations to save times
            {
                fac = 2.0*fac/rcm2;
                
                gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*fac;
                gamma[1] += xcm[0]*xcm[1]*fac;
            }
        }
    
	}
	else
	{
        // outside of the halo
		if (subtract_point == false)
		{
			PosType prefac = mass/rcm2/pi;
			alpha[0] += -1.0*prefac*xcm[0];
			alpha[1] += -1.0*prefac*xcm[1];
      
			// can turn off kappa and gamma calculations to save times
            {
				PosType tmp = -2.0*prefac/rcm2;
        
				gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
				gamma[1] += xcm[0]*xcm[1]*tmp;
			}
		}
	}
  
  
  // add stars for microlensing
  if(stars_N > 0 && stars_implanted)
  {
    force_stars(alpha,kappa,gamma,xcm);
  }
  
  return;
}




const long LensHaloHernquist::NTABLE = 100000;
const PosType LensHaloHernquist::maxrm = 100.0;
int LensHaloHernquist::count = 0;

PosType* LensHaloHernquist::xtable = NULL;
PosType* LensHaloHernquist::ftable = NULL;
PosType* LensHaloHernquist::gtable = NULL;
PosType* LensHaloHernquist::g2table = NULL;
PosType* LensHaloHernquist::htable = NULL;
PosType* LensHaloHernquist::xgtable = NULL;

/*
LensHaloHernquist::LensHaloHernquist()
: LensHalo(), gmax(0)
{
	make_tables();
	gmax = InterpolateFromTable(gtable,xmax);
}
*/
LensHaloHernquist::LensHaloHernquist(float my_mass,float my_Rmax,PosType my_zlens,float my_rscale,float my_fratio,float my_pa,int my_stars_N){
    
    mass=my_mass, Rmax=my_Rmax, zlens=my_zlens, rscale=my_rscale;
    fratio=my_fratio, pa=my_pa, stars_N=my_stars_N;
    stars_implanted = false;

    xmax = Rmax/rscale;
    make_tables();
	gmax = InterpolateFromTable(gtable,xmax);
    
    set_slope(1);
    /// If the axis ratio given in the parameter file is set to 1 all ellipticizing routines are skipped.
    if(fratio!=1){
        std::cout << "Hernquist constructor: slope set to " << get_slope() << std::endl;
        calcModes(fratio, get_slope(), pa, mod); // to ellipticize potential instead of kappa use (fratio, get_slope()-2, pa, mod)
        for(int i=1;i<Nmod;i++){
            if(mod[i]!=0){set_flag_elliptical(true);};
        }
    }
    
}

LensHaloHernquist::LensHaloHernquist(InputParams& params)
{
	assignParams(params);
	make_tables();
	gmax = InterpolateFromTable(gtable,xmax);
    
    set_slope(1);
    /// If the axis ratio given in the parameter file is set to 1 all ellipticizing routines are skipped.
    if(fratio!=1){
        std::cout << "Hernquist constructor: slope set to " << get_slope() << std::endl;
        calcModes(fratio, get_slope(), pa, mod); // to ellipticize potential instead of kappa use (fratio, get_slope()-2, pa, mod)
        for(int i=1;i<Nmod;i++){
            if(mod[i]!=0){set_flag_elliptical(true);};
        }
    }
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
        xgtable = new PosType[NTABLE];
        
		for(i = 0 ; i< NTABLE; i++){
			x = i*dx;
			xtable[i] = x;
			ftable[i] = ffunction(x);
			gtable[i] = gfunction(x);
			htable[i] = hfunction(x);
			g2table[i] = g2function(x);
            if(i==0){xgtable[i]=0;}
            if(i!=0){
                xgtable[i] = alpha_int(x);
            }
		}
  }
  count++;
}



PosType LensHaloHernquist::InterpolateFromTable(PosType *table, PosType y) const{
	int j;
	j=(int)(y/maxrm*NTABLE);

	assert(y>=xtable[j] && y<=xtable[j+1]);
	if (j==0)
		{
		if (table==ftable) return ffunction(y);
		if (table==gtable) return gfunction(y);
		if (table==g2table) return g2function(y);
		if (table==htable) return hfunction(y);
        if (table==xgtable) return alpha_int(y);
		}
	return (table[j+1]-table[j])/(xtable[j+1]-xtable[j])*(y-xtable[j]) + table[j];
}

void LensHaloHernquist::assignParams(InputParams& params){
	if(!params.get("main_mass",mass)) error_message1("main_mass",params.filename());
	if(!params.get("main_Rmax",Rmax)) error_message1("main_Rmax",params.filename());
	if(!params.get("main_zlens",zlens)) error_message1("main_zlens",params.filename());
	if(!params.get("main_rscale",rscale)) error_message1("main_rscale",params.filename());
  xmax = Rmax/rscale;
    if(!params.get("main_axis_ratio",fratio)){fratio=1; std::cout << "main_axis_ratio not defined in file " << params.filename() << ", hence set to 1." << std::endl;};
    if(!params.get("main_pos_angle",pa)){pa=0; std::cout << "main_pos_angle not defined in file " << params.filename() << ", hence set to 0." << std::endl;};
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
        delete[] xgtable;
	}
}

const long LensHaloJaffe::NTABLE = 100000;
const PosType LensHaloJaffe::maxrm = 100.0;
int LensHaloJaffe::count = 0;

PosType* LensHaloJaffe::xtable = NULL;
PosType* LensHaloJaffe::ftable = NULL;
PosType* LensHaloJaffe::gtable = NULL;
PosType* LensHaloJaffe::g2table = NULL;
PosType* LensHaloJaffe::xgtable = NULL;


//PosType* LensHaloJaffe::htable = NULL;
/*
LensHaloJaffe::LensHaloJaffe()
: LensHalo(), gmax(0)
{
	make_tables();
	gmax = InterpolateFromTable(gtable,xmax);
}
*/
LensHaloJaffe::LensHaloJaffe(float my_mass,float my_Rmax,PosType my_zlens,float my_rscale,float my_fratio,float my_pa,int my_stars_N){
    
    mass=my_mass, Rmax=my_Rmax, zlens=my_zlens, rscale=my_rscale;
    fratio=my_fratio, pa=my_pa, stars_N=my_stars_N;
    stars_implanted = false;
    xmax = Rmax/rscale;
    make_tables();
	gmax = InterpolateFromTable(gtable,xmax);
    
    set_slope(1);
    if(fratio!=1){
        std::cout << "Jaffe constructor: slope set to " << get_slope() << std::endl;
        calcModes(fratio, get_slope(), pa, mod);
        for(int i=1;i<Nmod;i++){
            if(mod[i]!=0){set_flag_elliptical(true);};
        }
    }

    
}

LensHaloJaffe::LensHaloJaffe(InputParams& params)
{
	assignParams(params);
	make_tables();
	gmax = InterpolateFromTable(gtable,xmax);
    
    set_slope(1);
    
    /// If the axis ratio given in the parameter file is set to 1 all ellipticizing routines are skipped.
    if(fratio!=1){
        std::cout << "Jaffe constructor: slope set to " << get_slope() << std::endl;
        calcModes(fratio, get_slope(), pa, mod);
        for(int i=1;i<Nmod;i++){
            if(mod[i]!=0){set_flag_elliptical(true);};
        }
        //std::cout << "if mods!=0 this must be 1: " << get_flag_elliptical() << std::endl;
    }
}

void LensHaloJaffe::make_tables(){
	if(count == 0){
		int i;
		PosType x, dx = maxrm/(PosType)NTABLE;

		xtable = new PosType[NTABLE];
		ftable = new PosType[NTABLE];
		gtable = new PosType[NTABLE];
		g2table = new PosType[NTABLE];
		xgtable = new PosType[NTABLE];

		for(i = 0 ; i< NTABLE; i++){
			x = i*dx;
			xtable[i] = x;
			ftable[i] = ffunction(x);
			gtable[i] = gfunction(x);
			g2table[i] = g2function(x);
            if(i==0){xgtable[i]=0;}
            if(i!=0){
                xgtable[i] = alpha_int(x);
            }
		}
  }
  count++;
}

PosType LensHaloJaffe::InterpolateFromTable(PosType *table, PosType y) const{
	int j;
	j=(int)(y/maxrm*NTABLE);

	assert(y>=xtable[j] && y<=xtable[j+1]);
	if (j==0)
		{
		if (table==ftable) return ffunction(y);
		if (table==gtable) return gfunction(y);
		if (table==g2table) return g2function(y);
 		if (table==xgtable) return alpha_int(y);
		}
	return (table[j+1]-table[j])/(xtable[j+1]-xtable[j])*(y-xtable[j]) + table[j];
}

void LensHaloJaffe::assignParams(InputParams& params){
	if(!params.get("main_mass",mass)) error_message1("main_mass",params.filename());
	if(!params.get("main_Rmax",Rmax)) error_message1("main_Rmax",params.filename());
	if(!params.get("main_zlens",zlens)) error_message1("main_zlens",params.filename());
	if(!params.get("main_rscale",rscale)) error_message1("main_rscale",params.filename());
  xmax = Rmax/rscale;
    if(!params.get("main_axis_ratio",fratio)){fratio=1; std::cout << "main_axis_ratio not defined in file " << params.filename() << ", hence set to 1." << std::endl;};
    if(!params.get("main_pos_angle",pa)){pa=0; std::cout << "main_pos_angle not defined in file " << params.filename() << ", hence set to 0." << std::endl;};
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
        delete[] xgtable;
	}
}





LensHaloDummy::LensHaloDummy()
: LensHalo()
{
//	mass = 0.;
}

LensHaloDummy::LensHaloDummy(float my_mass,float my_Rmax,PosType my_zlens,float my_rscale, int my_stars_N){
    mass=my_mass, Rmax=my_Rmax, zlens=my_zlens, rscale=my_rscale;
    stars_N=my_stars_N;
    stars_implanted = false;
    posHalo[0] = posHalo[1] = 0.0;
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


void LensHaloDummy::force_halo(PosType *alpha
                               ,KappaType *kappa
                               ,KappaType *gamma
                               ,KappaType *phi
                               ,PosType const *xcm
                               ,bool subtract_point)
{
	PosType rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
	PosType prefac = mass/rcm2/pi;
	PosType tmp = subtract_point*prefac;
	alpha[0] += tmp*xcm[0];
	alpha[1] += tmp*xcm[1];

    // intersecting, subtract the point particle
    if(subtract_point)
    {
        PosType x = sqrt(rcm2)/rscale;
       
        // can turn off kappa and gamma calculations to save times
        {
            *kappa += kappa_h(x)*prefac;

            tmp = (gamma_h(x) + 2.0*subtract_point)*prefac/rcm2;

            gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
            gamma[1] += xcm[0]*xcm[1]*tmp;
            
            *phi += phi_h(x);
        }
    }
    
    // add stars for microlensing
    if(stars_N > 0 && stars_implanted)
    {
        force_stars(alpha,kappa,gamma,xcm);
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

/// calculates the mass within radius R by integating kappa in theta and R, used only for testing
PosType LensHalo::MassBy2DIntegation(PosType R){
  LensHalo::DMDR dmdr(this);
  
  return Utilities::nintegrate<LensHalo::DMDR,PosType>(dmdr,-12,log(R),1.0e-5);
}
/// calculates the mass within radius R by integating alpha on a ring and using Gauss' law, used only for testing
PosType LensHalo::MassBy1DIntegation(PosType R){
  LensHalo::DMDTHETA dmdtheta(R,this);
  
  return R*Utilities::nintegrate<LensHalo::DMDTHETA,PosType>(dmdtheta, 0, 2*pi, 1.0e-6)/2;
}

