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
  mass = Rsize = Rmax = xmax = posHalo[0] = posHalo[1] = 0.0;
  stars_implanted = false;
  elliptical_flag = false;
  Dist = -1;
  stars_N =0.0;
  zlens = 0.0;
}

LensHalo::LensHalo(PosType z,const COSMOLOGY &cosmo){
  rscale = 1.0;
  mass = Rsize = Rmax = xmax = posHalo[0] = posHalo[1] = 0.0;
  stars_implanted = false;
  elliptical_flag = false;
  Dist = cosmo.angDist(z);
  stars_N =0.0;
  zlens = z;
}

LensHalo::LensHalo(InputParams& params,COSMOLOGY &cosmo,bool needRsize){
  Dist = -1;
  assignParams(params,needRsize);
  stars_implanted = false;
  posHalo[0] = posHalo[1] = 0.0;
  elliptical_flag = false;
  setDist(cosmo);
}
LensHalo::LensHalo(InputParams& params,bool needRsize){
  Dist = -1;
  assignParams(params,needRsize);
  stars_implanted = false;
  posHalo[0] = posHalo[1] = 0.0;
  elliptical_flag = false;
}

void LensHalo::initFromMassFunc(float my_mass, float my_Rsize, float my_rscale
                                , PosType my_slope, long *seed){
  mass = my_mass;
  Rsize = my_Rsize;
  rscale = my_rscale;
  xmax = Rsize/rscale;
}

void LensHalo::error_message1(std::string parameter,std::string file){
  ERROR_MESSAGE();
  std::cout << "Parameter " << parameter << " is needed to construct a LensHalo.  It needs to be set in parameter file " << file << "!" << std::endl;
  throw std::runtime_error(parameter);
}

void LensHalo::assignParams(InputParams& params,bool needRsize){
  if(!params.get("main_mass",mass)) error_message1("main_mass",params.filename());
  if(needRsize){  // this might not be required for lenses like NSIE
    if(!params.get("main_Rsize",Rsize)) error_message1("main_Rsize",params.filename());
  }
  if(!params.get("main_zlens",zlens)) error_message1("main_zlens",params.filename());
}

void LensHalo::PrintStars(bool show_stars)
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
  
  LensHalo::setRsize(1.0);
  Rmax = LensHalo::getRsize();

  LensHalo::setMass(0.0);
  LensHalo::setZlens(0);
  
  fratio=1;
  pa = stars_N = 0;
  stars_implanted = false;
  rscale = LensHalo::getRsize()/5;
  xmax = LensHalo::getRsize()/rscale;

  make_tables();
  gmax = InterpolateFromTable(gtable, xmax);
  set_flag_elliptical(false);
}

LensHaloNFW::LensHaloNFW(float my_mass,float my_Rsize,PosType my_zlens,float my_concentration,float my_fratio,float my_pa,int my_stars_N, EllipMethod my_ellip_method){

  LensHalo::setRsize(my_Rsize);
  LensHalo::setMass(my_mass);
  LensHalo::setZlens(my_zlens);


  fratio=my_fratio, pa=my_pa, stars_N=my_stars_N, main_ellip_method=my_ellip_method;
  stars_implanted = false;
  rscale = LensHalo::getRsize()/my_concentration;
  xmax = LensHalo::getRsize()/rscale;
  make_tables();
  gmax = InterpolateFromTable(gtable, xmax);
  
  
  set_slope(1);
  /// If the axis ratio given in the parameter file is set to 1 all ellipticizing routines are skipped.
  if(fratio!=1){
    Rmax = Rmax_to_Rsize_ratio*LensHalo::getRsize();
    //std::cout << getEllipMethod() << " method to ellipticise" << std::endl;
    if(getEllipMethod()==Fourier){
      std::cout << "NFW constructor: slope set to " << get_slope() << std::endl;
      calcModes(fratio, get_slope(), pa, mod); // to ellipticize potential instead of  kappa take calcModes(fratio, 2-get_slope(), pa, mod);
      std::cout << "NFW constructor: Fourier modes computed" << std::endl;
      for(int i=1;i<Nmod;i++){
        if(mod[i]!=0){set_flag_elliptical(true);};
      }
    }else set_flag_elliptical(true);
    if (getEllipMethod()==Pseudo or getEllipMethod()==Fourier){
      set_norm_factor();
      set_switch_flag(true);
    }
    
  }else{
    set_flag_elliptical(false);
    Rmax = LensHalo::getRsize();
  }
  
}

/* LensHalo::LensHalo(mass,Rsize,zlens, // base
 rscale,fratio,pa,stars_N, //NFW,Hernquist, Jaffe
 rscale,fratio,pa,beta // Pseudo NFW
 rscale,fratio,pa,sigma,rcore // NSIE
 zlens,stars_N // dummy
 ){
 
 stars_implanted = false;
 posHalo[0] = posHalo[1] = 0.0;
 }*/

LensHaloNFW::LensHaloNFW(InputParams& params):LensHalo(params)
{
  assignParams(params);
  make_tables();
  gmax = InterpolateFromTable(gtable, xmax);
  
  mnorm = renormalization(LensHalo::getRsize());
  
  std::cout << "mass normalization: " << mnorm << std::endl;
  
  // If the 2nd argument in calcModes(fratio, slope, pa, mod), the slope, is set to 1 it yields an elliptical kappa contour of given axis ratio (fratio) at the radius where the slope of the 3D density profile is -2, which is defined as the scale radius for the NFW profile. To ellipticize the potential instead of the convergence use calcModes(fratio, 2-get_slope(), pa, mod), this produces also an ellipse in the convergence map, but at the radius where the slope is 2-get_slope().
  set_slope(1);
  // If the axis ratio given in the parameter file is set to 1 all ellipticizing routines are skipped.
  if(fratio!=1){
    Rmax = Rmax_to_Rsize_ratio*LensHalo::getRsize();
    //std::cout << getEllipMethod() << " method to ellipticise" << std::endl;
    if(getEllipMethod()==Fourier){
      std::cout << "NFW constructor: slope set to " << get_slope() << std::endl;
      //for(int i=1;i<20;i++){
      //  calcModes(fratio, 0.1*i, pa, mod);
      //}
      //calcModes(fratio, get_slope()-0.5, pa, mod); // to ellipticize potential instead of  kappa take calcModes(fratio, 2-get_slope(), pa, mod);
      calcModes(fratio, get_slope(), pa, mod); // to ellipticize potential instead of  kappa take calcModes(fratio, 2-get_slope(), pa, mod);
      //calcModes(fratio, get_slope()+0.5, pa, mod); // to ellipticize potential instead of  kappa take calcModes(fratio, 2-get_slope(), pa, mod);
      
      for(int i=1;i<Nmod;i++){
        if(mod[i]!=0){set_flag_elliptical(true);};
      }
    }else set_flag_elliptical(true);
    if (getEllipMethod()==Pseudo){
      set_norm_factor();
    }
  }else{
    set_flag_elliptical(false);
    Rmax = LensHalo::getRsize();
  }
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
  //std::cout << "Interp: " << std::setprecision(7) << y-0.95 << " " << std::setprecision(7) << xtable[j]-0.95 << " " << xtable[j+1] <<std::endl;
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
  PosType tmp;
  if(!params.get("main_zlens",tmp)) error_message1("main_zlens",params.filename());
  if(!params.get("main_concentration",rscale)) error_message1("main_concentration",params.filename());
  if(!params.get("main_axis_ratio",fratio)){fratio=1; std::cout << "main_axis_ratio not defined in file " << params.filename() << ", hence set to 1." << std::endl;};
  if(!params.get("main_pos_angle",pa)){pa=0; std::cout << "main_pos_angle not defined in file " << params.filename() << ", hence set to 0." << std::endl;};
  if(!params.get("main_ellip_method",main_ellip_method)){if(fratio!=1){main_ellip_method=Pseudo;std::cout << "main_ellip_method is not defined in file " << params.filename() << ", hence set to Pseudo." << endl;};};
  
  rscale = LensHalo::getRsize()/rscale; // was the concentration
  xmax = LensHalo::getRsize()/rscale;
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
  
  LensHalo::setMass(my_mass);
  
  NFW_Utility nfw_util;
  // Find the NFW profile with the same mass, Vmax and R_halfmass
  nfw_util.match_nfw(vmax,r_halfmass,LensHalo::get_mass(),&rscale,&Rmax);
  LensHalo::setRsize(Rmax);
  rscale = LensHalo::getRsize()/rscale; // Was the concentration
  xmax = LensHalo::getRsize()/rscale;
  gmax = InterpolateFromTable(gtable,xmax);
 // std::cout << Rmax_halo << " " << LensHalo::getRsize() << std::endl;
  
}

void LensHaloNFW::initFromMassFunc(float my_mass, float my_Rsize, float my_rscale, PosType my_slope, long* seed)
{
  LensHalo::initFromMassFunc(my_mass, my_Rsize, my_rscale, my_slope, seed);
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
                                     ,float my_Rsize            /// maximum radius in Mpc
                                     ,PosType my_zlens         /// redshift
                                     ,float my_concentration   /// Rsize/rscale
                                     ,PosType my_beta          /// large r slope, see class description
                                     ,float my_fratio          /// axis ratio
                                     ,float my_pa              /// position angle
                                     ,int my_stars_N           /// number of stars, not yet implanted
                                     ,EllipMethod my_ellip_method /// ellipticizing method
)
{
  LensHalo::setMass(my_mass);
  LensHalo::setZlens(my_zlens);
  LensHalo::setRsize(my_Rsize);

  beta = my_beta;
  fratio = my_fratio;
  pa = my_pa;
  stars_N = my_stars_N;
  
  stars_implanted = false;
  rscale = LensHalo::getRsize()/my_concentration;
  xmax = LensHalo::getRsize()/rscale;
  
  make_tables();
  
  if(fratio!=1){
    Rmax = Rmax_to_Rsize_ratio*LensHalo::getRsize();
    //std::cout << getEllipMethod() << " method to ellipticise" << std::endl;
    if(getEllipMethod()==Fourier){
      std::cout << "Note: Fourier modes set to ellipticize kappa at slope main_slope+0.5, i.e. "<< get_slope()+0.5 << std::endl;
      calcModes(fratio, get_slope()+0.5, pa, mod);
      for(int i=1;i<Nmod;i++){
        if(mod[i]!=0){set_flag_elliptical(true);};
      }
    }else set_flag_elliptical(true);
    if (getEllipMethod()==Pseudo){
      set_norm_factor();
    }
  }else{
    set_flag_elliptical(false);
    Rmax = LensHalo::getRsize();
  }
}

/// The Fourier modes set to ellipticize kappa at slope main_slope+0.5, i.e. e.g. 1.5 for main_slope = 1. Note that set_slope is overridden for PseudoNFW to recalculate tables for different beta. But only fixed values of beta, i.e. 1,2 and >=3 are allowed!
LensHaloPseudoNFW::LensHaloPseudoNFW(InputParams& params):LensHalo(params)
{
  assignParams(params);
  make_tables();
  if(fratio!=1){
    Rmax = Rmax_to_Rsize_ratio*LensHalo::getRsize();
    //std::cout << getEllipMethod() << " method to ellipticise" << std::endl;
    if(getEllipMethod()==Fourier){
      std::cout << "Note: Fourier modes set to ellipticize kappa at slope main_slope+0.5, i.e. "<< get_slope()+0.5 << std::endl;
      calcModes(fratio, get_slope()+0.5, pa, mod);
      for(int i=1;i<Nmod;i++){
        //std::cout << mod[i] << std::endl;
        if(mod[i]!=0){set_flag_elliptical(true);};
      }
    }else set_flag_elliptical(true);
    if (getEllipMethod()==Pseudo){
      set_norm_factor();
    }
  }else{
    set_flag_elliptical(false);
    Rmax = LensHalo::getRsize();
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

void LensHaloPseudoNFW::initFromMassFunc(float my_mass, float my_Rsize, float my_rscale, PosType my_slope, long *seed){
  LensHalo::initFromMassFunc(my_mass,my_Rsize,my_rscale,my_slope,seed);
  beta = my_slope;
  xmax = my_Rsize/my_rscale;
  make_tables();
}

void LensHaloPseudoNFW::assignParams(InputParams& params){
  if(!params.get("main_concentration",rscale)) error_message1("main_concentration",params.filename());
  if(!params.get("main_slope",beta)) error_message1("main_slope",params.filename());
  if(!params.get("main_axis_ratio",fratio)){fratio=1; std::cout << "main_axis_ratio not defined in file " << params.filename() << ", hence set to 1." << std::endl;};
  if(!params.get("main_pos_angle",pa)){pa=0; std::cout << "main_pos_angle not defined in file " << params.filename() << ", hence set to 0." << std::endl;};
  if(!params.get("main_ellip_method",main_ellip_method)){if(fratio!=1){main_ellip_method=Pseudo;std::cout << "main_ellip_method is not defined in file " << params.filename() << ", hence set to Pseudo. CAUTION: Ellipticizing methods have not been tested with PseudoNFW yet!" << endl;};};
  
  rscale = LensHalo::getRsize()/rscale; // was the concentration
  xmax = LensHalo::getRsize()/rscale;
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
                                   ,float my_Rsize      /// maximum radius of halo in Mpc
                                   ,PosType my_zlens   /// redshift of halo
                                   ,PosType my_beta    /// logarithmic slop of surface density, kappa \propto r^{-beta}
                                   ,float my_fratio    /// axis ratio in asymetric case
                                   ,float my_pa        /// position angle
                                   ,int my_stars_N     /// number of stars, not yet implanted
                                   ,EllipMethod my_ellip_method /// ellipticizing method
){

  LensHalo::setMass(my_mass);
  LensHalo::setZlens(my_zlens);
  LensHalo::setRsize(my_Rsize);
  
  beta=my_beta;
  fratio=my_fratio, pa=my_pa, main_ellip_method=my_ellip_method, stars_N=my_stars_N;
  stars_implanted = false;
  rscale = 1;
  xmax = LensHalo::getRsize()/rscale ; /// xmax needs to be in initialized before the mass_norm_factor for Pseudo ellip method is calculated via  set_norm_factor()
  //mnorm = renormalization(get_Rmax());
  //std::cout << "PA in PowerLawConstructor: " << pa << std::endl;
  
  
  if(fratio!=1){
    Rmax = Rmax_to_Rsize_ratio*LensHalo::getRsize();
    //std::cout << getEllipMethod() << " method to ellipticise" << std::endl;
    if(getEllipMethod()==Fourier){
      calcModes(fratio, beta, pa, mod);
      for(int i=1;i<Nmod;i++){
        //std::cout << i << " " << mod[i] << std::endl;
        if(mod[i]!=0){set_flag_elliptical(true);};
      }
    }else set_flag_elliptical(true);
    if (getEllipMethod()==Pseudo){
      fratio=0.00890632+0.99209115*pow(fratio,0.33697702); /// translate fratio's needed for kappa into fratio used for potential
      //std::cout << "Pseudo-elliptical method requires fratio transformation, new fratio=" << fratio <<   std::endl;
      set_norm_factor();
    }
  }else{
    set_flag_elliptical(false);
    Rmax = LensHalo::getRsize();
  }
}

LensHaloPowerLaw::LensHaloPowerLaw(InputParams& params):LensHalo(params)
{
  assignParams(params);
  /// If the 2nd argument in calcModes(fratio, slope, pa, mod), the slope, is set to 1 it yields an elliptical kappa contour of given axis ratio (fratio) at the radius where the slope of the 3D density profile is -2, which is defined as the scale radius for the NFW profile. To ellipticize the potential instead of the convergence use calcModes(fratio, 2-get_slope(), pa, mod), this produces also an ellipse in the convergence map, but at the radius where the slope is 2-get_slope().
  /// If the axis ratio given in the parameter file is set to 1 all ellipticizing routines are skipped.
  
  // rscale = xmax = 1.0; // Commented in order to have a correct computation of the potential term in the time delay.
  // Replacing it by :
  rscale = 1;
  xmax = LensHalo::getRsize()/rscale;
  
  if(fratio!=1){
     Rmax = Rmax_to_Rsize_ratio*LensHalo::getRsize();

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
    //std::cout << getEllipMethod() << " method to ellipticise" << std::endl;
    if(getEllipMethod()==Fourier){
      calcModes(fratio, beta, pa, mod);
      for(int i=1;i<Nmod;i++){
        //std::cout << i << " " << mod[i] << " " << fratio << " " << beta <<  " " << pa << " " <<  std::endl;
        if(mod[i]!=0){set_flag_elliptical(true);};
      }
    }else set_flag_elliptical(true);
    if (getEllipMethod()==Pseudo){
      set_norm_factor();
    }
  }else{
    set_flag_elliptical(false);
    Rmax = LensHalo::getRsize();
  }
  // rscale = xmax = 1.0;
  // mnorm = renormalization(get_Rmax());
  mnorm = 1.;
  
  
}


void LensHaloPowerLaw::initFromMassFunc(float my_mass, float my_Rsize, float my_rscale, PosType my_slope, long *seed){
  LensHalo::initFromMassFunc(my_mass,my_Rsize,my_rscale,my_slope,seed);
  beta = my_slope;
  xmax = my_Rsize/my_rscale;
}

void LensHaloPowerLaw::assignParams(InputParams& params){
  if(!params.get("main_slope",beta)) error_message1("main_slope, example 1",params.filename());
  //if(beta>=2.0) error_message1("main_slope < 2",params.filename());
  if(!params.get("main_axis_ratio",fratio)){fratio=1; std::cout << "main_axis_ratio not defined in file " << params.filename() << ", hence set to 1." << std::endl;};
  if(!params.get("main_pos_angle",pa)){pa=0.0; std::cout << "main_pos_angle not defined in file " << params.filename() << ", hence set to 0." << std::endl;};
  if(!params.get("main_ellip_method",main_ellip_method)){if(fratio!=1){main_ellip_method=Pseudo;std::cout << "main_ellip_method is not defined in file " << params.filename() << ", hence set to Pseudo." << endl;};};
  
  if(!params.get("main_stars_N",stars_N)) error_message1("main_stars_N",params.filename());
  else if(stars_N){
    assignParams_stars(params);
  }
}

LensHaloPowerLaw::~LensHaloPowerLaw(){
  
}


/*
 LensHaloRealNSIE::LensHaloRealNSIE(float my_mass,float my_Rsize,PosType my_zlens,float my_rscale,float my_sigma
 , float my_rcore,float my_fratio,float my_pa,int my_stars_N){
 mass=my_mass, Rsize=my_Rsize, zlens=my_zlens, rscale=my_rscale;
 sigma=my_sigma, rcore=my_rcore;
 fratio=my_fratio, pa=my_pa, stars_N=my_stars_N;
 stars_implanted = false;
	Rsize = rmaxNSIE(sigma,mass,fratio,rcore);
	Rsize = MAX(1.0,1.0/fratio)*Rsize;  // redefine
	assert(Rmax >= Rsize);
 }
 */

size_t LensHaloRealNSIE::objectCount = 0;
std::vector<double> LensHaloRealNSIE::q_table;
std::vector<double> LensHaloRealNSIE::Fofq_table;

LensHaloRealNSIE::LensHaloRealNSIE(
                                   float my_mass
                                   ,PosType my_zlens
                                   ,float my_sigma   /// in km/s
                                   ,float my_rcore   /// in units of R_einstein
                                   ,float my_fratio
                                   ,float my_pa
                                   ,int my_stars_N)
:LensHalo(){
  rscale=1.0;
  LensHalo::setMass(my_mass);
  LensHalo::setZlens(my_zlens);


  sigma=my_sigma, rcore=my_rcore;
  fratio=my_fratio, pa=my_pa, stars_N=my_stars_N;
  stars_implanted = false;
  
  if(fratio  != 1.0) elliptical_flag = true;
  else elliptical_flag = false;
  ++objectCount;
  if(objectCount == 1){   // make table for calculating elliptical integrale
    construct_ellip_tables();
  }
  
  LensHalo::setRsize(rmax_calc());
  //std::cout << "NSIE " << Rsize << std::endl;
  Rmax = Rmax_to_Rsize_ratio*LensHalo::getRsize();
  if(fratio > 1.0 || fratio < 0.01) throw std::invalid_argument("invalid fratio");
  
  if(rcore > 0.0){
    LensHalo::setMass(MassBy1DIntegation(LensHalo::getRsize()));
  }
  
}

LensHaloRealNSIE::LensHaloRealNSIE(InputParams& params):LensHalo(params,false){
  sigma = 0.;
  fratio = 0.;
  pa = 0.;
  rcore = 0.;
  
  assignParams(params);

  if(fratio  != 1.0) elliptical_flag = true;
  else elliptical_flag = false;
  ++objectCount;
  if(objectCount == 1){   // make table for calculating elliptical integrale
    construct_ellip_tables();
  }
  //Rsize = rmaxNSIE(sigma,mass,fratio,rcore);
  //Rmax = MAX(1.0,1.0/fratio)*Rsize;  // redefine
  LensHalo::setRsize(rmax_calc());
  Rmax = Rmax_to_Rsize_ratio*LensHalo::getRsize();

  if(fratio > 1.0 || fratio < 0.01) throw std::invalid_argument("invalid fratio");
  
  if(rcore > 0.0){
    LensHalo::setMass(MassBy1DIntegation(LensHalo::getRsize()) );
  }

}

void LensHaloRealNSIE::assignParams(InputParams& params){
  if(!params.get("main_sigma",sigma)) error_message1("main_sigma",params.filename());
  if(!params.get("main_core",rcore)) error_message1("main_core",params.filename());
  if(!params.get("main_axis_ratio",fratio)) error_message1("main_axis_ratio",params.filename());
  else if(fratio > 1){
    ERROR_MESSAGE();
    std::cout << "parameter main_axis_ratio must be < 1 in file " << params.filename() << ". Use main_pos_angle to rotate the halo." << std::endl;
    exit(1);
  }
  
  if(!params.get("main_pos_angle",pa)) error_message1("main_pos_angle",params.filename());
  
  if(params.get("main_ellip_method",main_ellip_method)){std::cout << "main_ellip_method is NOT needed in file " << params.filename() << ". RealNSIE produces parametric ellipses!" << endl;};
  
  
  if(!params.get("main_stars_N",stars_N)) error_message1("main_stars_N",params.filename());
  else if(stars_N){
    assignParams_stars(params);
  }
  
}

LensHaloRealNSIE::~LensHaloRealNSIE(){
  --objectCount;
  if(objectCount == 0){
    q_table.resize(0);
    Fofq_table.resize(0);
  }
}

void LensHaloRealNSIE::construct_ellip_tables(){
  int N = 200;
  q_table.resize(N);
  Fofq_table.resize(N);
  q_table[0] = 0.01;
  NormFuncer funcer(q_table[0]);
  Fofq_table[0] = Utilities::nintegrate<NormFuncer,double>(funcer,0.0,PI/2,1.0e-6)*2/PI;
  for(int i=1 ; i < N-1 ; ++i){
    q_table[i] = i*1.0/(N-1) + 0.;
    NormFuncer funcer(q_table[i]);
    Fofq_table[i] = Utilities::nintegrate<NormFuncer,double>(funcer,0.0,PI/2,1.0e-6)*2/PI;
  }
  q_table.back() = 1.0;
  Fofq_table.back() = 1.0;
}

PosType LensHaloRealNSIE::rmax_calc(){
 
  if(fratio == 1.0 || rcore > 0.0)
    return sqrt( pow( LensHalo::get_mass()*Grav*lightspeed*lightspeed
                     *sqrt(fratio)/PI/sigma/sigma + rcore,2)
                                - rcore*rcore );
  
  // This is because there is no easy way of finding the circular Rmax for a fixed mass when rcore != 0
  //if(rcore > 0.0) throw std::runtime_error("rcore must be zero for this constructor");
  
  // asymmetric case
  //NormFuncer funcer(fratio);
  //double ellipticint = Utilities::nintegrate<NormFuncer,double>(funcer,0.0,PI/2,1.0e-6)*2/PI;
  
  double ellipticint = Utilities::InterpolateYvec(q_table,Fofq_table,fratio)*sqrt(fratio);
  
  return LensHalo::get_mass()*Grav*lightspeed*lightspeed/PI/sigma/sigma/ellipticint;
}

/*
 void LensHaloRealNSIE::initFromMass(float my_mass, long *seed){
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
 
 void LensHaloRealNSIE::initFromFile(float my_mass, long *seed, float vmax, float r_halfmass){
	initFromMass(my_mass,seed);
 }
 
 void LensHaloRealNSIE::initFromMassFunc(float my_mass, float my_Rmax, float my_rscale, PosType my_slope, long *seed){
	initFromMass(my_mass,seed);
 }
 */

/** \brief returns the lensing quantities of a ray in center of mass coordinates.
 *
 *  Warning: This adds to input value of alpha, kappa, gamma, and phi.  They need
 *  to be zeroed out if the contribution of just this halo is wanted.
 */
void LensHalo::force_halo(
                          PosType *alpha          /// deflection solar mass/Mpc
                          ,KappaType *kappa     /// surface density in Msun/Mpc^2 (?)
                          ,KappaType *gamma     /// three components of shear
                          ,KappaType *phi       /// potential in solar masses
                          ,PosType const *xcm   /// position relative to center (Mpc?)
                          ,bool subtract_point /// if true contribution from a point mass is subtracted
                          ,PosType screening   /// the factor by which to scale the mass for screening of the point mass subtraction
)
{
  if (elliptical_flag){
    force_halo_asym(alpha,kappa,gamma,phi,xcm,subtract_point,screening);
    //assert(!isinf(*kappa) );
  }else{
    force_halo_sym(alpha,kappa,gamma,phi,xcm,subtract_point,screening);
    //assert(!isinf(*kappa) );
  }
}

/** \brief returns the lensing quantities of a ray in center of mass coordinates for a symmetric halo
 *
 *  phi is defined here in such a way that it differs from alpha by a sign (as we should have alpha = \nabla_x phi)
 *  but alpha agrees with the rest of the lensing quantities (kappa and gammas).
 *  Warning : Be careful, the sign of alpha is changed in LensPlaneSingular::force !
 *
 */
void LensHalo::force_halo_sym(
                              PosType *alpha     /// solar mass/Mpc
                              ,KappaType *kappa  /// convergence
                              ,KappaType *gamma  /// three components of shear
                              ,KappaType *phi      /// potential solar masses
                              ,PosType const *xcm     /// position relative to center (Mpc)
                              ,bool subtract_point /// if true contribution from a point mass is subtracted
                              ,PosType screening   /// the factor by which to scale the mass for screening of the point mass subtraction
)
{
  
  PosType rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
  if(rcm2 < 1e-20) rcm2 = 1e-20;
  
  /// intersecting, subtract the point particle
  if(rcm2 < Rmax*Rmax)
  {
    PosType prefac = mass/rcm2/PI;
    PosType x = sqrt(rcm2)/rscale;
    // PosType xmax = Rmax/rscale;
    PosType tmp = (alpha_h(x) + 1.0*subtract_point)*prefac;
    alpha[0] += tmp*xcm[0];
    alpha[1] += tmp*xcm[1];
    
    *kappa += kappa_h(x)*prefac;
    
    tmp = (gamma_h(x) + 2.0*subtract_point) * prefac / rcm2; // ;
    gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
    gamma[1] += xcm[0]*xcm[1]*tmp;
 
    *phi += phi_h(x) * mass / PI ;
    if(subtract_point) *phi -= 0.5 * log(rcm2) * mass / PI;
  }
  else // the point particle is not subtracted
  {
    if (subtract_point == false)
    {
      PosType prefac = screening*mass/rcm2/PI;
      alpha[0] += -1.0 * prefac * xcm[0];
      alpha[1] += -1.0 * prefac * xcm[1];
      
      //std::cout << "rcm2  = " << rcm2 << std::endl;
      //std::cout << "prefac  = " << prefac << std::endl;
      //std::cout << "xcm  = " << xcm[0] << " " << xcm[1] << std::endl;
      
      PosType tmp = -2.0*prefac/rcm2;
      
      // kappa is equal to 0 in the point mass case.
      
      gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
      gamma[1] += xcm[0]*xcm[1]*tmp;
      *phi += 0.5 * log(rcm2) * mass / PI ;
    }
  }
  
  /// add stars for microlensing
  if(stars_N > 0 && stars_implanted)
  {
    force_stars(alpha,kappa,gamma,xcm);
  }
  
  //(alpha[0] == alpha[0] && alpha[1] == alpha[1]);

  return;
}
// TODO: put in some comments about the units used
void LensHalo::force_halo_asym(
                               PosType *alpha     /// mass/Mpc
                               ,KappaType *kappa
                               ,KappaType *gamma
                               ,KappaType *phi      /// potential solar masses
                               ,PosType const *xcm
                               ,bool subtract_point /// if true contribution from a point mass is subtracted
                               ,PosType screening   /// the factor by which to scale the mass for screening of the point mass subtraction
){
  
  //float r_size=get_rsize()*Rmax;
  //Rmax=r_size*1.2;

  PosType rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
  PosType alpha_tmp[2],kappa_tmp,gamma_tmp[2],phi_tmp;
  
  if(rcm2 < 1e-20) rcm2 = 1e-20;
  //std::cout << "rsize , rmax,  mass_norm =" << Rsize << " , " << Rmax << " , " << mass_norm_factor << std::endl;
  //std::cout << subtract_point << std::endl;
  
  /// intersecting, subtract the point particle
  if(rcm2 < Rmax*Rmax){
    
    double r = sqrt(rcm2); // devision by rscale for isotropic halos (see above) here taken out because not used for any halo. it should be taken out for the isotropic case too, if others not affected / make use of it ;
    double theta;
    if(xcm[0] == 0.0 && xcm[1] == 0.0) theta = 0.0;
    else theta=atan2(xcm[1],xcm[0]);
    if(rcm2 > LensHalo::getRsize()*LensHalo::getRsize()){
      PosType alpha_iso[2],alpha_ellip[2];
      alpha_ellip[0] = alpha_ellip[1] = 0;
      if(main_ellip_method==Pseudo){alphakappagamma_asym(LensHalo::getRsize(),theta, alpha_tmp,&kappa_tmp,gamma_tmp,&phi_tmp);}
      if(main_ellip_method==Fourier){alphakappagamma1asym(LensHalo::getRsize(),theta, alpha_tmp,&kappa_tmp,gamma_tmp,&phi_tmp);}
      if(main_ellip_method==Schramm){alphakappagamma2asym(LensHalo::getRsize(),theta, alpha_tmp,&kappa_tmp,gamma_tmp,&phi_tmp);}
      if(main_ellip_method==Keeton){alphakappagamma3asym(LensHalo::getRsize(),theta, alpha_tmp,&kappa_tmp,gamma_tmp,&phi_tmp);}
      alpha_ellip[0]=alpha_tmp[0]*mass_norm_factor;
      alpha_ellip[1]=alpha_tmp[1]*mass_norm_factor;
      double f1 = (Rmax - r)/(Rmax - LensHalo::getRsize()),f2 = (r - LensHalo::getRsize())/(Rmax - LensHalo::getRsize());

     // PosType tmp = mass/Rmax/PI/r;
      PosType tmp = mass/rcm2/PI;
      alpha_iso[0] = -1.0*tmp*xcm[0];
      alpha_iso[1] = -1.0*tmp*xcm[1];
      alpha[0] += alpha_iso[0]*f2 + alpha_ellip[0]*f1;
      alpha[1] += alpha_iso[1]*f2 + alpha_ellip[1]*f1;
      
      {
        PosType tmp = -2.0*mass/rcm2/PI/rcm2;
        gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
        gamma[1] += xcm[0]*xcm[1]*tmp;
        //gamma[0] += 0.5*gamma_tmp[0]*mass_norm_factor;
        //gamma[1] += 0.5*gamma_tmp[1]*mass_norm_factor;
        
        *phi += phi_tmp;
      }
      
    }else{
      if(main_ellip_method==Pseudo){alphakappagamma_asym(r,theta, alpha_tmp,&kappa_tmp,gamma_tmp,&phi_tmp);}
      if(main_ellip_method==Fourier){alphakappagamma1asym(r,theta, alpha_tmp,&kappa_tmp,gamma_tmp,&phi_tmp);}
      if(main_ellip_method==Schramm){alphakappagamma2asym(r,theta, alpha_tmp,&kappa_tmp,gamma_tmp,&phi_tmp);}
      if(main_ellip_method==Keeton){alphakappagamma3asym(r,theta, alpha_tmp,&kappa_tmp,gamma_tmp,&phi_tmp);}
      
      alpha[0] +=  alpha_tmp[0]*mass_norm_factor;//-1.0*subtract_point*mass/rcm2/PI*xcm[0];
      alpha[1] +=  alpha_tmp[1]*mass_norm_factor;//-1.0*subtract_point*mass/rcm2/PI*xcm[1];

      if(get_switch_flag()==true){  /// case distinction used for elliptical NFWs (true) only (get_switch_flag==true)
        *kappa += kappa_tmp*mass_norm_factor*mass_norm_factor;
        gamma[0] += 0.5*gamma_tmp[0]*mass_norm_factor*mass_norm_factor;
        gamma[1] += 0.5*gamma_tmp[1]*mass_norm_factor*mass_norm_factor;
      }else{
        *kappa += kappa_tmp*mass_norm_factor;

        gamma[0] += 0.5*gamma_tmp[0]*mass_norm_factor;//+1.0*subtract_point*mass/rcm2/PI/rcm2*0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1]);
        gamma[1] += 0.5*gamma_tmp[1]*mass_norm_factor;//-1.0*subtract_point*mass/rcm2/PI/rcm2*(xcm[0]*xcm[1]);
        
        //if(theta < 0.660 && theta> 0.659){
        //assert(mass_norm_factor==1);
        //std::cout << theta << " , " <<  0.5*gamma_tmp[0]*mass_norm_factor << " , " << 0.5*gamma_tmp[1]*mass_norm_factor << " , " << kappa_tmp*mass_norm_factor << " , " << alpha_tmp[0]*mass_norm_factor << " , " << alpha_tmp[1]*mass_norm_factor  <<  std::endl;
        //}
      }
      
      /*if (rcm2 < 1E-6){
        std::cout << kappa_tmp*mass_norm_factor << " " << 0.5*gamma_tmp[0]*mass_norm_factor<< " "  << 0.5*gamma_tmp[1]*mass_norm_factor << " " <<rcm2 << " " << alpha_tmp[0]*mass_norm_factor << " " << alpha_tmp[1]*mass_norm_factor << std::endl;
      }
      */
      *phi += phi_tmp;

    }
    
    
    if(subtract_point){
      //std::cout << "DO WE EVEN GET HERE??" << std::endl;
      PosType tmp =  screening*mass/PI/rcm2; // *mass_norm_factor
      alpha[0] +=  tmp*xcm[0];
      alpha[1] +=  tmp*xcm[1];
      
      tmp = 2.0*tmp/rcm2;
      gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
      gamma[1] += xcm[0]*xcm[1]*tmp;
      
      *phi -= 0.5 * log(rcm2) * mass / PI ; // *mass_norm_factor
    }
    
  }
  else // the point particle is not subtracted
  {
    if (subtract_point == false)
    {
      PosType prefac = mass/rcm2/PI;
      alpha[0] += -1.0 * prefac * xcm[0];
      alpha[1] += -1.0 * prefac * xcm[1];
      
      //if(rcm2==1.125){
      //std::cout << "rcm2  = " << rcm2 << std::endl;
      //std::cout << "prefac  = " << prefac << std::endl;
      //std::cout << "xcm  = " << xcm[0] << " " << xcm[1] << std::endl;
      //}
      
      PosType tmp = -2.0*prefac/rcm2;
      
      // kappa is equal to 0 in the point mass case.
      
      gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
      gamma[1] += xcm[0]*xcm[1]*tmp;
      
      *phi += 0.5 * log(rcm2) * mass / PI ;
    }
  }
  
  /// add stars for microlensing
  if(stars_N > 0 && stars_implanted)
  {
    force_stars(alpha,kappa,gamma,xcm);
  }
  
  //assert(alpha[0] == alpha[0] && alpha[1] == alpha[1]);

  return;
}



/* *
 void LensHaloRealNSIE::force_halo(
 PosType *alpha
 ,KappaType *kappa
 ,KappaType *gamma
 ,KappaType *phi
 ,PosType const *xcm
 ,bool subtract_point /// if true contribution from a point mass is subtracted
 ,PosType screening   /// the factor by which to scale the mass for screening of the point mass subtraction
 ){
 
	PosType rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
	if(rcm2 < 1e-20) rcm2 = 1e-20;
 
 // **** test line
 
	if(rcm2 < Rmax*Rmax){
 PosType ellipR = ellipticRadiusNSIE(xcm,fratio,pa);
 if(ellipR > LensHalo::getRsize()){
 // This is the case when the ray is within the NSIE's circular region of influence but outside its elliptical truncation
 
 PosType alpha_out[2],alpha_in[2],rin,x_in[2];
 PosType prefac = -1.0*mass/Rmax/PI;
 PosType r = sqrt(rcm2);
 float units = pow(sigma/lightspeed,2)/Grav/sqrt(fratio); // mass/distance(physical)
 
 alpha_in[0] = alpha_in[1] = 0;
 
 rin = r*LensHalo::getRsize()/ellipR;
 
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
 PosType prefac = mass/rcm2/PI;
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
 PosType fac = mass/rcm2/PI;
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
void LensHaloRealNSIE::force_halo(
                                  PosType *alpha
                                  ,KappaType *kappa
                                  ,KappaType *gamma
                                  ,KappaType *phi
                                  ,PosType const *xcm
                                  ,bool subtract_point /// if true contribution from a point mass is subtracted
                                  ,PosType screening   /// the factor by which to scale the mass for screening of the point mass subtraction
)
{
  PosType rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
 
  
  if(rcm2 < 1e-20) rcm2 = 1e-20;
  if(rcm2 < Rmax*Rmax){
    //PosType ellipR = ellipticRadiusNSIE(xcm,fratio,pa);
    float units = pow(sigma/lightspeed,2)/Grav;///sqrt(fratio); // mass/distance(physical)
   // std::cout << "rsize , rmax,  mass_norm =" << LensHalo::getRsize() << " , " << Rmax << " , " << mass_norm_factor << std::endl;
    if(rcm2 > LensHalo::getRsize()*LensHalo::getRsize())
      //if(ellipR > LensHalo::getRsize()*LensHalo::getRsize())
    {
      // This is the case when the ray is within the NSIE's circular region of influence but outside its elliptical truncation
      
      PosType alpha_iso[2],alpha_ellip[2];
      //PosType prefac = -1.0*mass/Rmax/PI;
      PosType r = sqrt(rcm2);
      
      //double Rin = sqrt(rcm2)*LensHalo::getRsize()/ellipR;
      //std::cout << Rmax << " " << LensHalo::getRsize() << " " << Rmax/LensHalo::getRsize() << std::endl;
      double f1 = (Rmax - r)/(Rmax - LensHalo::getRsize()),f2 = (r - LensHalo::getRsize())/(Rmax - LensHalo::getRsize());
      //double f1 = (Rmax - r)/(Rmax - Rin),f2 = (r - Rin)/(Rmax - Rin);
      
      // SIE solution
      alpha_ellip[0] = alpha_ellip[1] = 0;
      alphaNSIE(alpha_ellip,xcm,fratio,rcore,pa);
      alpha_ellip[0] *= units;
      alpha_ellip[1] *= units;
      
      // SIS solution
      //alpha_iso[0] = alpha_iso[1] = 0;
      //alphaNSIE(alpha_iso,xcm,1,rcore,pa);
      //alpha_iso[0] *= units;
      //alpha_iso[1] *= units;
      //
      
      // point mass solution
      // PosType tmp = mass/rcm2/PI;
      PosType tmp = LensHalo::get_mass()/rcm2/PI;
      alpha_iso[0] = -1.0*tmp*xcm[0];
      alpha_iso[1] = -1.0*tmp*xcm[1];
      
      alpha[0] += alpha_iso[0]*f2 + alpha_ellip[0]*f1;
      alpha[1] += alpha_iso[1]*f2 + alpha_ellip[1]*f1;
      {
        PosType tmp = -2.0*LensHalo::get_mass()/rcm2/PI/rcm2;
        
        gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
        gamma[1] += xcm[0]*xcm[1]*tmp;
        
      }
      
    }else{
      PosType xt[2]={0,0},tmp[2]={0,0};
      xt[0]=xcm[0];
      xt[1]=xcm[1];
      alphaNSIE(tmp,xt,fratio,rcore,pa);
      //alpha[0] = units*tmp[0];  // minus sign removed because already included in alphaNSIE
      //alpha[1] = units*tmp[1];  // Why was the "+=" removed?
      alpha[0] += units*tmp[0];//*sqrt(fratio);
      alpha[1] += units*tmp[1];//*sqrt(fratio);
      {
        KappaType tmp[2]={0,0};
        *kappa += units*kappaNSIE(xt,fratio,rcore,pa);///sqrt(fratio);
        gammaNSIE(tmp,xt,fratio,rcore,pa);
        gamma[0] += units*tmp[0];
        gamma[1] += units*tmp[1];
      }
    }
    
    
    if(subtract_point)
    {
      PosType fac = screening*LensHalo::get_mass()/rcm2/PI;
      alpha[0] += fac*xcm[0];
      alpha[1] += fac*xcm[1];
      
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
      PosType prefac = LensHalo::get_mass()/rcm2/PI;
      alpha[0] += -1.0*prefac*xcm[0];
      alpha[1] += -1.0*prefac*xcm[1];
      
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
  
  //assert(alpha[0] == alpha[0] && alpha[1] == alpha[1]);

  return;
}
/**/



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
LensHaloHernquist::LensHaloHernquist(float my_mass,float my_Rsize,PosType my_zlens,float my_rscale,float my_fratio,float my_pa,int my_stars_N, EllipMethod my_ellip_method){
  
  LensHalo::setMass(my_mass);
  LensHalo::setZlens(my_zlens);
  LensHalo::setRsize(my_Rsize);

  rscale=my_rscale;
  fratio=my_fratio, pa=my_pa, stars_N=my_stars_N;
  stars_implanted = false;
  
  xmax = LensHalo::getRsize()/rscale;
  make_tables();
  gmax = InterpolateFromTable(gtable,xmax);
  
  set_slope(1);
  /// If the axis ratio given in the parameter file is set to 1 all ellipticizing routines are skipped.
  if(fratio!=1){
    Rmax = Rmax_to_Rsize_ratio*LensHalo::getRsize();
    //std::cout << getEllipMethod() << " method to ellipticise" << std::endl;
    if(getEllipMethod()==Fourier){
      //std::cout << "Hernquist constructor: slope set to " << get_slope() << std::endl;
      calcModes(fratio, get_slope(), pa, mod); // to ellipticize potential instead of kappa use (fratio, get_slope()-2, pa, mod)
      for(int i=1;i<Nmod;i++){
        if(mod[i]!=0){set_flag_elliptical(true);};
      }
    }else set_flag_elliptical(true);
    if (getEllipMethod()==Pseudo){
      fratio=0.00890632+0.99209115*pow(fratio,0.33697702);
      set_norm_factor();
    }
  }else{
    set_flag_elliptical(false);
    Rmax = LensHalo::getRsize();
  }
}

LensHaloHernquist::LensHaloHernquist(InputParams& params): LensHalo(params)
{
  assignParams(params);
  make_tables();
  gmax = InterpolateFromTable(gtable,xmax);
  
  set_slope(1);
  /// If the axis ratio given in the parameter file is set to 1 all ellipticizing routines are skipped.
  if(fratio!=1){
    //std::cout << getEllipMethod() << " method to ellipticise" << std::endl;
    if(getEllipMethod()==Fourier){
      std::cout << "Hernquist constructor: slope set to " << get_slope() << std::endl;
      calcModes(fratio, get_slope(), pa, mod); // to ellipticize potential instead of kappa use (fratio, get_slope()-2, pa, mod)
      for(int i=1;i<Nmod;i++){
        if(mod[i]!=0){set_flag_elliptical(true);};
      }
    }else set_flag_elliptical(true);
    if (getEllipMethod()==Pseudo){
      set_norm_factor();
    }
  }else set_flag_elliptical(false);
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
  if(!params.get("main_rscale",rscale)) error_message1("main_rscale",params.filename());
  xmax = LensHalo::getRsize()/rscale;
  if(!params.get("main_axis_ratio",fratio)){fratio=1; std::cout << "main_axis_ratio not defined in file " << params.filename() << ", hence set to 1." << std::endl;};
  if(!params.get("main_pos_angle",pa)){pa=0; std::cout << "main_pos_angle not defined in file " << params.filename() << ", hence set to 0." << std::endl;};
  if(!params.get("main_ellip_method",main_ellip_method)){if(fratio!=1){main_ellip_method=Pseudo;std::cout << "main_ellip_method is not defined in file " << params.filename() << ", hence set to Pseudo." << endl;};};
  if(!params.get("main_stars_N",stars_N)) error_message1("main_stars_N",params.filename());
  else if(stars_N){
    assignParams_stars(params);
    
    std::cout << "LensHalo::getRsize() " << LensHalo::getRsize() <<std::endl;
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
LensHaloJaffe::LensHaloJaffe(float my_mass,float my_Rsize,PosType my_zlens,float my_rscale,float my_fratio,float my_pa,int my_stars_N, EllipMethod my_ellip_method){
  
  LensHalo::setMass(my_mass);
  LensHalo::setZlens(my_zlens);
  LensHalo::setRsize(my_Rsize);

  rscale=my_rscale;
  fratio=my_fratio, pa=my_pa, stars_N=my_stars_N;
  stars_implanted = false;
  xmax = LensHalo::getRsize()/rscale;
  make_tables();
  gmax = InterpolateFromTable(gtable,xmax);
  
  set_slope(1);
  if(fratio!=1){
    Rmax = Rmax_to_Rsize_ratio*LensHalo::getRsize();
    //std::cout << getEllipMethod() << " method to ellipticise" << std::endl;
    if(getEllipMethod()==Fourier){
      //std::cout << "Jaffe constructor: slope set to " << get_slope() << std::endl;
      calcModes(fratio, get_slope(), pa, mod);
      for(int i=1;i<Nmod;i++){
        if(mod[i]!=0){set_flag_elliptical(true);};
      }
    }else set_flag_elliptical(true);
    if (getEllipMethod()==Pseudo){
      fratio=0.00890632+0.99209115*pow(fratio,0.33697702);
      set_norm_factor();
    }
  }else{
    set_flag_elliptical(false);
    Rmax = LensHalo::getRsize();
  }
}

LensHaloJaffe::LensHaloJaffe(InputParams& params): LensHalo(params)
{
  assignParams(params);
  make_tables();
  gmax = InterpolateFromTable(gtable,xmax);
  
  set_slope(1);
  
  /// If the axis ratio given in the parameter file is set to 1 all ellipticizing routines are skipped.
  if(fratio!=1){
    //std::cout << getEllipMethod() << " method to ellipticise" << std::endl;
    if(getEllipMethod()==Fourier){
      std::cout << "Jaffe constructor: slope set to " << get_slope() << std::endl;
      calcModes(fratio, get_slope(), pa, mod);
      for(int i=1;i<Nmod;i++){
        if(mod[i]!=0){set_flag_elliptical(true);};
      }
    }else set_flag_elliptical(true);
    if (getEllipMethod()==Pseudo){
      set_norm_factor();
    }
  }else set_flag_elliptical(false);
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
  if(!params.get("main_rscale",rscale)) error_message1("main_rscale",params.filename());
  xmax = LensHalo::getRsize()/rscale;
  if(!params.get("main_axis_ratio",fratio)){fratio=1; std::cout << "main_axis_ratio not defined in file " << params.filename() << ", hence set to 1." << std::endl;};
  if(!params.get("main_pos_angle",pa)){pa=0; std::cout << "main_pos_angle not defined in file " << params.filename() << ", hence set to 0." << std::endl;};
  
  if(!params.get("main_ellip_method",main_ellip_method)){if(fratio!=1){main_ellip_method=Pseudo;std::cout << "main_ellip_method is not defined in file " << params.filename() << ", hence set to Pseudo." << endl;};};
  
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

LensHaloDummy::LensHaloDummy(float my_mass,float my_Rsize,PosType my_zlens,float my_rscale, int my_stars_N){
  LensHalo::setMass(my_mass);
  LensHalo::setZlens(my_zlens);
  LensHalo::setRsize(my_Rsize);

  rscale=my_rscale;
  stars_N=my_stars_N;
  stars_implanted = false;
  setTheta(0.0,0.0);
}

LensHaloDummy::LensHaloDummy(InputParams& params): LensHalo(params)
{
  assignParams(params);
  //	mass = 0.;
}

void LensHaloDummy::initFromMassFunc(float my_mass, float my_Rsize, float my_rscale, PosType my_slope, long *seed){
  LensHalo::setMass(1.e-10);
  LensHalo::setRsize(my_Rsize);
  
  rscale = my_rscale;
  xmax = LensHalo::getRsize()/rscale;
  Rmax = LensHalo::getRsize();
}


void LensHaloDummy::force_halo(PosType *alpha
                               ,KappaType *kappa
                               ,KappaType *gamma
                               ,KappaType *phi
                               ,PosType const *xcm
                               ,bool subtract_point
                               ,PosType screening   /// the factor by which to scale the mass for screening of the point mass subtraction
)
{
  PosType rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
  PosType prefac = LensHalo::get_mass()/rcm2/PI;
  PosType tmp = subtract_point*prefac;
  alpha[0] += tmp*xcm[0];
  alpha[1] += tmp*xcm[1];
  
  // intersecting, subtract the point particle
  if(subtract_point)
  {
    PosType x = screening*sqrt(rcm2)/rscale;
    
    *kappa += kappa_h(x)*prefac;
      
    tmp = (gamma_h(x) + 2.0*subtract_point)*prefac/rcm2;
      
    gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
    gamma[1] += xcm[0]*xcm[1]*tmp;
      
    *phi += phi_h(x);
    
  }
  
  // add stars for microlensing
  if(stars_N > 0 && stars_implanted)
  {
    force_stars(alpha,kappa,gamma,xcm);
  }
  
  //assert(alpha[0] == alpha[0] && alpha[1] == alpha[1]);

}

void LensHaloDummy::assignParams(InputParams& params)
{  
  if(params.get("main_ellip_method",main_ellip_method)){std::cout << "main_ellip_method is NOT needed in file " << params.filename() << ". LensHaloDummy does not require ellipticity!" << endl;};
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
  
  return R*Utilities::nintegrate<LensHalo::DMDTHETA,PosType>(dmdtheta, 0, 2*PI, 1.0e-6)/2;
}

/// calculates the average gamma_t for LensHalo::test()
double LensHalo::test_average_gt(PosType R){
  struct test_gt_func f(R,this);
  return Utilities::nintegrate<test_gt_func>(f,0.0,2.*PI,1.0e-3);
}

// returns <kappa> x 2PI on a ring at R
double LensHalo::test_average_kappa(PosType R){
  struct test_kappa_func f(R,this);
  return Utilities::nintegrate<test_kappa_func>(f,0.0,2.*PI,1.0e-3);
}

/// Three tests: 1st - Mass via 1D integration vs mass via 2D integration. 2nd: gamma_t=alpha/r - kappa(R) which can be used for spherical distributions. Deviations are expected for axis ratios <1. For the latter case we use the next test. 3rd: The average along a circular aperture of gamma_t should be equal to <kappa(<R)> minus the average along a circular aperture over kappa. Note that also  alpha/r - kappa is checked for consistency with kappa(<R)-<kappa(R)>. For axis ratios < 1 the factor between the two is expected to be of order O(10%).
bool LensHalo::test(){
  std::cout << "test alpha's consistance with kappa by comparing mass interior to a radius by 1D integration and Gauss' law and by 2D integration" << std::endl << "  The total internal mass is " << mass << std::endl;
  
  std::cout << "R/Rmax     R/Rsize     Mass 1 D (from alpha)     Mass 2 D         (m1 - m2)/m1       m2/m1" << std::endl;
  
  int N=25;
  PosType m1,m2;
  for(int i=1;i<N;++i){
    m1 = MassBy1DIntegation(LensHalo::getRsize()*i/(N-4));
    m2 = MassBy2DIntegation(LensHalo::getRsize()*i/(N-4));
    std::cout <<  i*1./(N-4) << "      " << i/(N-4) << "      " << m1 << "       "
    << m2 << "        "<< (m1-m2)/m1 << "         " << m2/m1  << std::endl;
    
  }
  
  
  PosType r;
  
  std::cout << "test gamma_t's consistance with kappa and alpha by comparing gamma_t to alpha/r - kappa along the x-axis" << std::endl
  << "Not expected to be equal for asymmetric cases."<< std::endl;
  std::cout << std::endl <<"R/Rmax         R/Rsize         gamma_t       alpha/r - kappa          alpha/r           kappa        delta/gt "  << std::endl;
  for(int i=1;i<N;++i){
    r = Rsize*i/(N-2);
    
    PosType alpha[2] = {0,0},x[2] = {0,0};
    KappaType kappa = 0,gamma[3] = {0,0,0} ,phi=0;
    
    x[0] = r;
    x[1] = 0;
    
    force_halo(alpha,&kappa,gamma,&phi,x);
    
    std::cout << r/Rmax << "      " << r/Rsize << "       " <<  -gamma[0]  << "         " << -alpha[0]/r - kappa << "         " << -alpha[0]/r << "         " << kappa << "      " <<  (alpha[0]/r + kappa)/gamma[0]   <<std::endl;
  }
  
  
  std::cout << "test average tangential shear's, gamma_t's, consistance with the average convergence at a radius, kappa(r) and average kappa within a radius calculated using alpha and Gauss' law.  gamma_t should be equal to <kappa>_R - kappa(R)" << std::endl;
  std::cout << std::endl <<"R/Rmax         R/Rsize        gamma_t       <kappa>_R-kappa(R)       <kappa>_R            kappa(R)     [<kappa>_R-kappa(R)]/gt  " << std::endl;
  
  
  for(int i=1;i<N;++i){
    r = Rsize*i/(N-2);
    
    //integrate over t
    PosType average_gt, average_kappa;
    average_gt=test_average_gt(r)/2/PI;
    average_kappa=test_average_kappa(r)/2/PI;
    m1 = MassBy1DIntegation(r)/PI/r/r;
    
    
    std::cout << r/Rmax << "       "  << r/Rsize << "       " << -1.0*average_gt << "         " << m1-average_kappa << "         " <<  m1 << "          " <<  average_kappa << "         "  << -1.0*(m1-average_kappa)/average_gt << std::endl;
    
    
    PosType alpha[2] = {0,0},x[2] = {0,0};
    KappaType kappa = 0,gamma[3] = {0,0,0} ,phi=0;
    x[0] = r;
    x[1] = 0;
    force_halo(alpha,&kappa,gamma,&phi,x);
    /*
     assert( abs(-1.0*(m1-average_kappa)/average_gt-1.) < 1e-2 ); // <g_t> = <k(<R)>-<kappa(R)> test
     if(!elliptical_flag){
     assert( abs(abs(-alpha[0]/r)/m1-1.) < 1e-1 ); // alpha/r ~ <kappa(R)>
     assert( abs(abs(alpha[0]/r + kappa)/gamma[0])-1.0 < 1e-2); // g_t = alpha/r - kappa test
     }
     */
    
    // This is a list of possible assertion test that can be made
    //if(!elliptical_flag){
    //std::cout << abs(abs(-alpha[0]/r)/m1-1.)  << std::endl;
    //assert( abs(abs(-alpha[0]/r)/m1-1.) < 1e-1 ); // alpha/r ~ <kappa(R)>
    //std::cout << abs(abs(alpha[0]/r + kappa)/gamma[0])-1.0 << std::endl;
    //assert( abs(abs(alpha[0]/r + kappa)/gamma[0])-1.0 < 1e-2); // g_t = alpha/r - kappa test
    //}
    //std::cout << abs(-1.0*(m1-average_kappa)/average_gt-1.) << std::endl;
    //std::cout << abs( -alpha[0]/r - kappa ) / abs(m1-average_kappa  ) -1.  << std::endl;
    //assert( abs( -alpha[0]/r - kappa ) / abs(m1-average_kappa  ) -1.  < 1 ); // alpha/r ~ <kappa(R)>
    
    
  }
  
  return true;
};

/// The following functions calculate the integrands of the Schramm 1990 method to obtain elliptical halos
PosType LensHalo::DALPHAXDM::operator()(PosType m){
  
  double ap = m*m*a2 + lambda,bp = m*m*b2 + lambda;
  double p2 = x[0]*x[0]/ap/ap/ap/ap + x[1]*x[1]/bp/bp/bp/bp;  // actually the inverse of equation (5) in Schramm 1990
  PosType tmp = m*(isohalo->getRsize());
  KappaType kappa=0;
  
  double xiso=tmp/isohalo->rscale;
  //PosType alpha[2]={0,0},tm[2] = {m*(isohalo->getRsize()),0};
  //KappaType kappam=0,gamma[2]={0,0},phi;
  
  kappa=isohalo->kappa_h(xiso)/PI/xiso/xiso*isohalo->mass;
  
  //isohalo->force_halo_sym(alpha,&kappam,gamma,&phi,tm);
  //std::cout << "kappa: " << kappa << " " << kappam << " " << kappa/kappam << std::endl;
  assert(kappa >= 0.0);
  std::cout << "output x: " << m << " " << m*kappa/(ap*ap*ap*bp*p2) << std::endl;
  
  return m*kappa/(ap*ap*ap*bp*p2); // integrand of equation (28) in Schramm 1990
}

PosType LensHalo::DALPHAYDM::operator()(PosType m){
  
  double ap = m*m*a2 + lambda,bp = m*m*b2 + lambda;
  double p2 = x[0]*x[0]/ap/ap/ap/ap + x[1]*x[1]/bp/bp/bp/bp;  // actually the inverse of equation (5) in Schramm 1990
  PosType tmp = m*(isohalo->getRsize());
  KappaType kappa=0;
  double xiso=tmp/isohalo->rscale;
  //PosType alpha[2]={0,0},tm[2] = {m*(isohalo->getRsize()),0};
  //KappaType kappam=0,gamma[2]={0,0},phi;
  
  kappa=isohalo->kappa_h(xiso)/PI/xiso/xiso*isohalo->mass;
  //isohalo->force_halo_sym(alpha,&kappam,gamma,&phi,tm);
  //std::cout << "kappa: " << kappa << " " << kappam << " " << kappa/kappam <<   std::endl;
  assert(kappa >= 0.0);
  return m*kappa/(ap*bp*bp*bp*p2); // integrand of equation (29) in Schramm 1990
}

//bool LensHaloZcompare(LensHalo *lh1,LensHalo *lh2){return (lh1->getZlens() < lh2->getZlens());}
//bool compare(LensHalo *lh1,LensHalo *lh2){return (lh1->getZlens() < lh2->getZlens());}
