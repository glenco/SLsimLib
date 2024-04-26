/*
 * lens_halos.cpp
 *
 *  Created on: 06.05.2013
 *      Author: mpetkova
 */

#include "slsimlib.h"
#include "lens_halos.h"

using namespace std;

/// Shell constructor
LensHalo::LensHalo(){
  rscale = 1.0;
  mass = Rsize = Rmax = xmax = posHalo[0] = posHalo[1] = 0.0;
  //stars_implanted = false;
  elliptical_flag = false;
  Dist = -1;
  //stars_N =0.0;
  zlens = 0.0;
  tag=0;
}

LensHalo::LensHalo(PosType z,const COSMOLOGY &cosmo){
  rscale = 1.0;
  mass = Rsize = Rmax = xmax = posHalo[0] = posHalo[1] = 0.0;
  //stars_implanted = false;
  elliptical_flag = false;
  Dist = cosmo.angDist(z);
  //stars_N =0.0;
  zlens = z;
  tag=0;
}

//LensHalo::LensHalo(InputParams& params,COSMOLOGY &cosmo,bool needRsize){
//  Dist = -1;
//  assignParams(params,needRsize);
//  stars_implanted = false;
//  posHalo[0] = posHalo[1] = 0.0;
//  elliptical_flag = false;
//  setDist(cosmo);
//}
//LensHalo::LensHalo(InputParams& params,bool needRsize){
//  Dist = -1;
//  assignParams(params,needRsize);
//  stars_implanted = false;
//  posHalo[0] = posHalo[1] = 0.0;
//  elliptical_flag = false;
//}

LensHalo::LensHalo(const LensHalo &h){
  
  idnumber = h.idnumber; /// Identification number of halo.  It is not always used.
  Dist = h.Dist;
  posHalo[0] = h.posHalo[0]; posHalo[1] = h.posHalo[1];
  zlens = h. zlens;
  mass = h.mass;
  Rsize = h.Rsize;
  mnorm = h.mnorm;
  Rmax = h.Rmax;
  
//  stars_index = h.stars_index;
//  stars_xp = h.stars_xp;
//
//  stars_N = h.stars_N;
//  star_theta_force = h.star_theta_force;
//
//  if(stars_N > 1){
//    star_tree = new TreeQuadParticles<StarType>(stars_xp.data(),stars_N
//                                              ,false,false,0,4,star_theta_force);
//  }else{
//    star_tree = nullptr;
//  }
//  star_massscale = h.star_massscale;
//  star_fstars = h.star_fstars;
//
//  star_Nregions = h.star_Nregions;
//  star_region = h.star_region;
  
  beta = h.beta;
  
  Rmax_to_Rsize_ratio = h.Rmax_to_Rsize_ratio;
  rscale = h.rscale;
  
//  stars_implanted = h.stars_implanted;
//  main_stars_imf_type = h.main_stars_imf_type;
//  main_stars_min_mass = h. main_stars_min_mass;
//  main_stars_max_mass = h.main_stars_max_mass;
  main_ellip_method = h.main_ellip_method;
//  bend_mstar =h.bend_mstar;
//  lo_mass_slope = h.lo_mass_slope;
//  hi_mass_slope = h.hi_mass_slope;
//
//  star_Sigma = h.star_Sigma;
//  star_xdisk = h.star_xdisk;
  
  xmax = h.xmax;
  mass_norm_factor = h.mass_norm_factor;
  pa = h.pa;
  fratio = h.fratio;
  elliptical_flag = h.elliptical_flag;
  switch_flag = h.switch_flag;
  //Nmod = h.Nmod;
  
  tag = h.tag;
};

LensHalo & LensHalo::operator=(LensHalo &&h){
  
  if (this != &h){
  
    idnumber = h.idnumber; /// Identification number of halo.  It is not always used.
    Dist = h.Dist;
    posHalo[0] = h.posHalo[0]; posHalo[1] = h.posHalo[1];
    zlens = h. zlens;
    mass = h.mass;
    Rsize = h.Rsize;
    mnorm = h.mnorm;
    Rmax = h.Rmax;
  
//    stars_index = h.stars_index;
//    stars_xp = h.stars_xp;
//
//    stars_N = h.stars_N;
//    star_theta_force = h.star_theta_force;
//
//    delete star_tree;
//    star_tree = h.star_tree;
//    h.star_tree = nullptr;
//
//    star_massscale = h.star_massscale;
//    star_fstars = h.star_fstars;
//
//    star_Nregions = h.star_Nregions;
//    star_region = h.star_region;
  
    beta = h.beta;
  
    Rmax_to_Rsize_ratio = h.Rmax_to_Rsize_ratio;
    rscale = h.rscale;
  
//    stars_implanted = h.stars_implanted;
//    main_stars_imf_type = h.main_stars_imf_type;
//    main_stars_min_mass = h. main_stars_min_mass;
//    main_stars_max_mass = h.main_stars_max_mass;
    main_ellip_method = h.main_ellip_method;
//    bend_mstar =h.bend_mstar;
//    lo_mass_slope = h.lo_mass_slope;
//    hi_mass_slope = h.hi_mass_slope;
//
//    star_Sigma = h.star_Sigma;
//    star_xdisk = h.star_xdisk;
  
    xmax = h.xmax;
    mass_norm_factor = h.mass_norm_factor;
    pa = h.pa;
    fratio = h.fratio;
    elliptical_flag = h.elliptical_flag;
    switch_flag = h.switch_flag;
    //Nmod = h.Nmod;
    
    tag = h.tag;
  }
  return *this;
};


LensHalo & LensHalo::operator=(const LensHalo &h){
  
  if(this == &h) return *this;
  
  idnumber = h.idnumber; /// Identification number of halo.  It is not always used.
  Dist = h.Dist;
  posHalo[0] = h.posHalo[0]; posHalo[1] = h.posHalo[1];
  zlens = h. zlens;
  mass = h.mass;
  Rsize = h.Rsize;
  mnorm = h.mnorm;
  Rmax = h.Rmax;
  
//  stars_index = h.stars_index;
//  stars_xp = h.stars_xp;
//
//  stars_N = h.stars_N;
//  star_theta_force = h.star_theta_force;
//
//  if(stars_N > 1){
//    star_tree = new TreeQuadParticles<StarType>(stars_xp.data(),stars_N
//                                              ,false,false,0,4,star_theta_force);
//  }else{
//    star_tree = nullptr;
//  }
//  star_massscale = h.star_massscale;
//  star_fstars = h.star_fstars;
//
//  star_Nregions = h.star_Nregions;
//  star_region = h.star_region;
  
  beta = h.beta;
  
  Rmax_to_Rsize_ratio = h.Rmax_to_Rsize_ratio;
  rscale = h.rscale;
  
//  stars_implanted = h.stars_implanted;
//  main_stars_imf_type = h.main_stars_imf_type;
//  main_stars_min_mass = h. main_stars_min_mass;
//  main_stars_max_mass = h.main_stars_max_mass;
  main_ellip_method = h.main_ellip_method;
//  bend_mstar =h.bend_mstar;
//  lo_mass_slope = h.lo_mass_slope;
//  hi_mass_slope = h.hi_mass_slope;
//
//  star_Sigma = h.star_Sigma;
//  star_xdisk = h.star_xdisk;
  
  xmax = h.xmax;
  mass_norm_factor = h.mass_norm_factor;
  pa = h.pa;
  fratio = h.fratio;
  elliptical_flag = h.elliptical_flag;
  switch_flag = h.switch_flag;
  //Nmod = h.Nmod;
  tag = h.tag;
  
  return *this;
};

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

//void LensHalo::PrintStars(bool show_stars)
//{
//  std::cout << std::endl << "Nstars "<<stars_N << std::endl << std::endl;
//  if(stars_N>0){
//    if(star_Nregions > 0)
//      std::cout << "stars_Nregions "<<star_Nregions << std::endl;
//    std::cout << "stars_massscale "<<star_massscale << std::endl;
//    std::cout << "stars_fstars "<<star_fstars << std::endl;
//    std::cout << "stars_theta_force "<<star_theta_force << std::endl;
//    if(show_stars){
//      if(stars_implanted){
//        for(int i=0 ; i < stars_N ; ++i) std::cout << "    x["<<i<<"]="
//						    << stars_xp[i][0] << " " << stars_xp[i][1] << std::endl;
//      }else std::cout << "stars are not implanted yet" << std::endl;
//    }
//  }
//}

PixelMap<double> LensHalo::map_variables(
                       LensingVariable lensvar /// lensing variable - KAPPA, ALPHA1, ALPHA2, GAMMA1, GAMMA2 or PHI  in units of sigma crit
                       ,size_t Nx
                       ,size_t Ny
                       ,double res
){
  
  Point_2d center;
  PixelMap<double> map(center.data(),Nx,Ny,res);
  Point_2d x,alpha;
  size_t N = Nx*Ny;
  KappaType kappa,phi,gamma[3];
  for(size_t i = 0 ; i < N ; ++i){
    map.find_position(x.data(),i);
    force_halo(alpha.data(),&kappa,gamma,&phi,x.data());
    
    switch (lensvar) {
      case LensingVariable::ALPHA:
        map[i] = sqrt(alpha[0]*alpha[0] + alpha[1]*alpha[1]);
        break;
      case LensingVariable::ALPHA1:
        map[i] = alpha[0];
        break;
      case LensingVariable::ALPHA2:
        map[i] = alpha[1];
        break;
      case LensingVariable::GAMMA1:
        map[i] = gamma[0];
        break;
      case LensingVariable::GAMMA2:
        map[i] = gamma[1];
        break;
      case LensingVariable::KAPPA:
        map[i] = kappa;
        break;
      case LensingVariable::PHI:
        map[i] = phi;
        break;
        
      default:
        std::cerr << "Error :  LensHalo::map_variables - lensing variable not acceptable " << std::endl;
        throw std::invalid_argument("bad lensing variable");
    }
  }
  return map;
}


// calculates the deflection etc. caused by stars alone
//void LensHalo::force_stars(
//                           PosType *alpha     /// mass/Mpc
//                           ,KappaType *kappa
//                           ,KappaType *gamma
//                           ,PosType const *xcm     /// physical position on lens plane
//)
//{
//  PosType alpha_tmp[2];
//  KappaType gamma_tmp[3], tmp = 0;
//  KappaType phi;
//
//  gamma_tmp[0] = gamma_tmp[1] = gamma_tmp[2] = 0.0;
//  alpha_tmp[0] = alpha_tmp[1] = 0.0;
//
//  substract_stars_disks(xcm,alpha,kappa,gamma);
//
//	 // do stars with tree code
//  star_tree->force2D_recur(xcm,alpha_tmp,&tmp,gamma_tmp,&phi);
//
//  alpha[0] -= star_massscale*alpha_tmp[0];
//  alpha[1] -= star_massscale*alpha_tmp[1];
//
//  {
//    *kappa += star_massscale*tmp;
//    gamma[0] -= star_massscale*gamma_tmp[0];
//    gamma[1] -= star_massscale*gamma_tmp[1];
//  }
//
//}

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
  LensHalo::setZlens(0,COSMOLOGY(CosmoParamSet::Planck18));
  LensHalo::Dist = -1;  // to be set later
  
  fratio=1;
  pa = 0;
  //stars_N = 0;
  //stars_implanted = false;
  rscale = LensHalo::getRsize()/5;
  xmax = LensHalo::getRsize()/rscale;

  make_tables();
  gmax = InterpolateFromTable(gtable, xmax);
  set_flag_elliptical(false);
  ++count;
}

LensHaloNFW::LensHaloNFW(float my_mass,float my_Rsize,PosType my_zlens,float my_concentration
                         ,float my_fratio,float my_pa,const COSMOLOGY &cosmo
                         ,EllipMethod my_ellip_method
                         ):LensHalo(my_zlens,cosmo)
{

  LensHalo::setRsize(my_Rsize);
  LensHalo::setMass(my_mass);

  fratio=my_fratio;
  pa= PI/2 - my_pa;
  main_ellip_method=my_ellip_method;
  rscale = LensHalo::getRsize()/my_concentration;
  xmax = LensHalo::getRsize()/rscale;
  make_tables();
  gmax = InterpolateFromTable(gtable, xmax);
  
  
  set_slope(1);
  /// If the axis ratio given in the parameter file is set to 1 all ellipticizing routines are skipped.
  if(fratio!=1){
    Rmax = Rmax_to_Rsize_ratio*LensHalo::getRsize();
    //std::cout << getEllipMethod() << " method to ellipticise" << std::endl;
    if(getEllipMethod()==EllipMethod::Fourier){
      std::cout << "NFW constructor: slope set to " << get_slope() << std::endl;
      calcModes(fratio, get_slope(), pa, mod); // to ellipticize potential instead of  kappa take calcModes(fratio, 2-get_slope(), pa, mod);
      std::cout << "NFW constructor: Fourier modes computed" << std::endl;
      for(int i=1;i<Nmod;i++){
        if(mod[i]!=0){set_flag_elliptical(true);};
      }
    }else set_flag_elliptical(true);
    if (getEllipMethod()==EllipMethod::Pseudo or getEllipMethod()==EllipMethod::Fourier){
      set_norm_factor();
      set_switch_flag(true);
    }
    
  }else{
    set_flag_elliptical(false);
    Rmax = LensHalo::getRsize();
  }
  
  ++count;
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

//LensHaloNFW::LensHaloNFW(InputParams& params):LensHalo(params)
//{
//  assignParams(params);
//  make_tables();
//  gmax = InterpolateFromTable(gtable, xmax);
//  
//  mnorm = renormalization(LensHalo::getRsize());
//  
//  std::cout << "mass normalization: " << mnorm << std::endl;
//  
//  // If the 2nd argument in calcModes(fratio, slope, pa, mod), the slope, is set to 1 it yields an elliptical kappa contour of given axis ratio (fratio) at the radius where the slope of the 3D density profile is -2, which is defined as the scale radius for the NFW profile. To ellipticize the potential instead of the convergence use calcModes(fratio, 2-get_slope(), pa, mod), this produces also an ellipse in the convergence map, but at the radius where the slope is 2-get_slope().
//  set_slope(1);
//  // If the axis ratio given in the parameter file is set to 1 all ellipticizing routines are skipped.
//  if(fratio!=1){
//    Rmax = Rmax_to_Rsize_ratio*LensHalo::getRsize();
//    //std::cout << getEllipMethod() << " method to ellipticise" << std::endl;
//    if(getEllipMethod()==Fourier){
//      std::cout << "NFW constructor: slope set to " << get_slope() << std::endl;
//      //for(int i=1;i<20;i++){
//      //  calcModes(fratio, 0.1*i, pa, mod);
//      //}
//      //calcModes(fratio, get_slope()-0.5, pa, mod); // to ellipticize potential instead of  kappa take calcModes(fratio, 2-get_slope(), pa, mod);
//      calcModes(fratio, get_slope(), pa, mod); // to ellipticize potential instead of  kappa take calcModes(fratio, 2-get_slope(), pa, mod);
//      //calcModes(fratio, get_slope()+0.5, pa, mod); // to ellipticize potential instead of  kappa take calcModes(fratio, 2-get_slope(), pa, mod);
//      
//      for(int i=1;i<Nmod;i++){
//        if(mod[i]!=0){set_flag_elliptical(true);};
//      }
//    }else set_flag_elliptical(true);
//    if (getEllipMethod()==Pseudo){
//      set_norm_factor();
//    }
//  }else{
//    set_flag_elliptical(false);
//    Rmax = LensHalo::getRsize();
//  }
//}

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
  if(!params.get("main_ellip_method",main_ellip_method)){if(fratio!=1){main_ellip_method=EllipMethod::Pseudo;std::cout << "main_ellip_method is not defined in file " << params.filename() << ", hence set to Pseudo." << endl;};};
  
  rscale = LensHalo::getRsize()/rscale; // was the concentration
  xmax = LensHalo::getRsize()/rscale;
}

LensHaloNFW::~LensHaloNFW(){
  --LensHaloNFW::count;
  if(LensHaloNFW::count == 0){
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
  make_tables();
  ++count;
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
                                     ,const COSMOLOGY &cosmo
                                     ,EllipMethod my_ellip_method /// ellipticizing method
)
{
  LensHalo::setMass(my_mass);
  LensHalo::setZlens(my_zlens,cosmo);
  LensHalo::setRsize(my_Rsize);

  beta = my_beta;
  fratio = my_fratio;
  pa = my_pa;
  rscale = LensHalo::getRsize()/my_concentration;
  xmax = LensHalo::getRsize()/rscale;
  
  make_tables();
  
  if(fratio!=1){
    Rmax = Rmax_to_Rsize_ratio*LensHalo::getRsize();
    //std::cout << getEllipMethod() << " method to ellipticise" << std::endl;
    if(getEllipMethod()==EllipMethod::Fourier){
      std::cout << "Note: Fourier modes set to ellipticize kappa at slope main_slope+0.5, i.e. "<< get_slope()+0.5 << std::endl;
      calcModes(fratio, get_slope()+0.5, pa, mod);
      for(int i=1;i<Nmod;i++){
        if(mod[i]!=0){set_flag_elliptical(true);};
      }
    }else set_flag_elliptical(true);
    if (getEllipMethod()==EllipMethod::Pseudo){
      set_norm_factor();
    }
  }else{
    set_flag_elliptical(false);
    Rmax = LensHalo::getRsize();
  }
  ++count;
}

// The Fourier modes set to ellipticize kappa at slope main_slope+0.5, i.e. e.g. 1.5 for main_slope = 1. Note that set_slope is overridden for PseudoNFW to recalculate tables for different beta. But only fixed values of beta, i.e. 1,2 and >=3 are allowed!
//LensHaloPseudoNFW::LensHaloPseudoNFW(InputParams& params):LensHalo(params)
//{
//  assignParams(params);
//  make_tables();
//  if(fratio!=1){
//    Rmax = Rmax_to_Rsize_ratio*LensHalo::getRsize();
//    //std::cout << getEllipMethod() << " method to ellipticise" << std::endl;
//    if(getEllipMethod()==Fourier){
//      std::cout << "Note: Fourier modes set to ellipticize kappa at slope main_slope+0.5, i.e. "<< get_slope()+0.5 << std::endl;
//      calcModes(fratio, get_slope()+0.5, pa, mod);
//      for(int i=1;i<Nmod;i++){
//        //std::cout << mod[i] << std::endl;
//        if(mod[i]!=0){set_flag_elliptical(true);};
//      }
//    }else set_flag_elliptical(true);
//    if (getEllipMethod()==Pseudo){
//      set_norm_factor();
//    }
//  }else{
//    set_flag_elliptical(false);
//    Rmax = LensHalo::getRsize();
//  }
//}

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
  ++count;
}

void LensHaloPseudoNFW::assignParams(InputParams& params){
  if(!params.get("main_concentration",rscale)) error_message1("main_concentration",params.filename());
  if(!params.get("main_slope",beta)) error_message1("main_slope",params.filename());
  if(!params.get("main_axis_ratio",fratio)){fratio=1; std::cout << "main_axis_ratio not defined in file " << params.filename() << ", hence set to 1." << std::endl;};
  if(!params.get("main_pos_angle",pa)){pa=0; std::cout << "main_pos_angle not defined in file " << params.filename() << ", hence set to 0." << std::endl;};
  if(!params.get("main_ellip_method",main_ellip_method)){if(fratio!=1){main_ellip_method=EllipMethod::Pseudo;std::cout << "main_ellip_method is not defined in file " << params.filename() << ", hence set to Pseudo. CAUTION: Ellipticizing methods have not been tested with PseudoNFW yet!" << endl;};};
  
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
                                   ,const COSMOLOGY &cosmo
                                   ,EllipMethod my_ellip_method /// ellipticizing method
){

  if(my_beta >= 2.0){
    std::cerr << "ERROR: Power-law index in LensHaloPowerLaw cannot be >= 2" << std::endl;
    throw std::invalid_argument("bad gamma");
  }
  if(my_Rsize <= 0.0 || std::isnan(my_Rsize)){
    std::cerr << "ERROR : my_Rsize = " << my_Rsize << std::endl;
    throw std::runtime_error("Bad halo");
  }
  
  LensHalo::setMass(my_mass);
  LensHalo::setZlens(my_zlens,cosmo);
  LensHalo::setRsize(my_Rsize);
  
  beta=my_beta;
  fratio=my_fratio, pa=my_pa, main_ellip_method=my_ellip_method;
  rscale = 1;
  xmax = LensHalo::getRsize()/rscale ; /// xmax needs to be in initialized before the mass_norm_factor for Pseudo ellip method is calculated via  set_norm_factor()
  //mnorm = renormalization(get_Rmax());
  //std::cout << "PA in PowerLawConstructor: " << pa << std::endl;
  
  
  
  
  if(fratio!=1){
    Rmax = Rmax_to_Rsize_ratio*LensHalo::getRsize();
    //std::cout << getEllipMethod() << " method to ellipticise" << std::endl;
    if(getEllipMethod()==EllipMethod::Fourier){
      calcModes(fratio, beta, pa, mod);
      for(int i=1;i<Nmod;i++){
        //std::cout << i << " " << mod[i] << std::endl;
        if(mod[i]!=0){set_flag_elliptical(true);};
      }
    }else set_flag_elliptical(true);
    if (getEllipMethod()==EllipMethod::Pseudo){
      fratio=0.00890632+0.99209115*pow(fratio,0.33697702); /// translate fratio's needed for kappa into fratio used for potential
      //std::cout << "Pseudo-elliptical method requires fratio transformation, new fratio=" << fratio <<   std::endl;
      set_norm_factor();
    }
  }else{
    set_flag_elliptical(false);
    Rmax = LensHalo::getRsize();
  }
  
  if(xmax <= 0.0 || std::isnan(xmax)){
    std::cerr << "ERROR : xmax = " << xmax << std::endl;
    std::cerr << "        rscale = " << rscale << std::endl;
    std::cerr << "        Rmax = " << Rmax << std::endl;

    throw std::runtime_error("Bad halo");
  }
}

//LensHaloPowerLaw::LensHaloPowerLaw(InputParams& params):LensHalo(params)
//{
//  assignParams(params);
//  /// If the 2nd argument in calcModes(fratio, slope, pa, mod), the slope, is set to 1 it yields an elliptical kappa contour of given axis ratio (fratio) at the radius where the slope of the 3D density profile is -2, which is defined as the scale radius for the NFW profile. To ellipticize the potential instead of the convergence use calcModes(fratio, 2-get_slope(), pa, mod), this produces also an ellipse in the convergence map, but at the radius where the slope is 2-get_slope().
//  /// If the axis ratio given in the parameter file is set to 1 all ellipticizing routines are skipped.
//  
//  // rscale = xmax = 1.0; // Commented in order to have a correct computation of the potential term in the time delay.
//  // Replacing it by :
//  rscale = 1;
//  xmax = LensHalo::getRsize()/rscale;
//  
//  if(fratio!=1){
//     Rmax = Rmax_to_Rsize_ratio*LensHalo::getRsize();
//
//     //std::cout << getEllipMethod() << " method to ellipticise" << std::endl;
//    if(getEllipMethod()==Fourier){
//      calcModes(fratio, beta, pa, mod);
//      for(int i=1;i<Nmod;i++){
//        //std::cout << i << " " << mod[i] << " " << fratio << " " << beta <<  " " << pa << " " <<  std::endl;
//        if(mod[i]!=0){set_flag_elliptical(true);};
//      }
//    }else set_flag_elliptical(true);
//    if (getEllipMethod()==Pseudo){
//      set_norm_factor();
//    }
//  }else{
//    set_flag_elliptical(false);
//    Rmax = LensHalo::getRsize();
//  }
//  // rscale = xmax = 1.0;
//  // mnorm = renormalization(get_Rmax());
//  mnorm = 1.;
//}


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
  if(!params.get("main_ellip_method",main_ellip_method)){if(fratio!=1){main_ellip_method=EllipMethod::Pseudo;std::cout << "main_ellip_method is not defined in file " << params.filename() << ", hence set to Pseudo." << endl;};};
  
//  if(!params.get("main_stars_N",stars_N)) error_message1("main_stars_N",params.filename());
//  else if(stars_N){
//    assignParams_stars(params);
//  }
}

LensHaloPowerLaw::~LensHaloPowerLaw(){
  
}


LensHaloTNSIE::LensHaloTNSIE(
                                   float my_mass
                                   ,PosType my_zlens
                                   ,float my_sigma
                                   ,float my_rcore
                                   ,float my_fratio
                                   ,float my_pa
                                   ,const COSMOLOGY &cosmo
                                   ,float f)
:LensHalo(),sigma(my_sigma),fratio(my_fratio)
,pa(PI/2 - my_pa),rcore(my_rcore)
{
  rscale=1.0;
  LensHalo::setMass(my_mass);
  LensHalo::setZlens(my_zlens,cosmo);

  if(fratio > 1.0 || fratio < 0.01) throw std::invalid_argument("invalid fratio");
  
  units = sigma*sigma/lightspeed/lightspeed/Grav;///sqrt(fratio); // mass/distance(physical);
  
  rtrunc = my_mass*sqrt(my_fratio)/units/PI + rcore;
  Rmax = f * rtrunc;
  LensHalo::setRsize(Rmax);
}

void LensHaloTNSIE::force_halo(
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
 
  if(force_point(alpha,kappa,gamma,phi,xcm,rcm2
                 ,subtract_point,screening)) return;

  if(rcm2 < 1.0e-5) rcm2 = 1e-5;
  if(rcm2 < Rmax*Rmax){
 
    PosType tmp[2]={0,0};
    alphaNSIE(tmp,xcm,fratio,rcore,pa);
    alpha[0] += units*tmp[0];
    alpha[1] += units*tmp[1];

    alphaNSIE(tmp,xcm,fratio,rtrunc,pa);
    alpha[0] -= units*tmp[0];
    alpha[1] -= units*tmp[1];
    {
      KappaType tmp[2]={0,0};
      *kappa += units*(kappaNSIE(xcm,fratio,rcore,pa) - kappaNSIE(xcm,fratio,rtrunc,pa) );
                       
      gammaNSIE(tmp,xcm,fratio,rcore,pa);
      gamma[0] += units*tmp[0];
      gamma[1] += units*tmp[1];

      gammaNSIE(tmp,xcm,fratio,rtrunc,pa);
      gamma[0] -= units*tmp[0];
      gamma[1] -= units*tmp[1];
    }
    
//    if(subtract_point)
//    {
//      PosType fac = screening*LensHalo::get_mass()/rcm2/PI;
//      alpha[0] += fac*xcm[0];
//      alpha[1] += fac*xcm[1];
//
//      {
//        fac = 2.0*fac/rcm2;
//
//        gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*fac;
//        gamma[1] += xcm[0]*xcm[1]*fac;
//      }
//    }
  }
//  else
//  {
//    // outside of the halo
//    if (subtract_point == false)
//    {
//      PosType prefac = LensHalo::get_mass()/rcm2/PI;
//      alpha[0] += -1.0*prefac*xcm[0];
//      alpha[1] += -1.0*prefac*xcm[1];
//
//      {
//        PosType tmp = -2.0*prefac/rcm2;
//
//        gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
//        gamma[1] += xcm[0]*xcm[1]*tmp;
//      }
//    }
//  }
//
  return;
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
                                   float my_mass     /// mass, sets truncation radius
                                   ,PosType my_zlens /// redshift
                                   ,float my_sigma   /// in km/s
                                   ,float my_rcore   /// core radius
                                   ,float my_fratio  /// axis ratio
                                   ,float my_pa      /// postion angle
                                   ,const COSMOLOGY &cosmo)
:LensHalo(){
  rscale=1.0;
  LensHalo::setMass(my_mass);
  LensHalo::setZlens(my_zlens,cosmo);


  sigma=my_sigma, rcore=my_rcore;
  fratio=my_fratio, pa = PI/2 - my_pa;//, stars_N=my_stars_N;
  //stars_implanted = false;
  
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
  
  units = pow(sigma/lightspeed,2)/Grav;///sqrt(fratio); // mass/distance(physical);
}

//LensHaloRealNSIE::LensHaloRealNSIE(InputParams& params):LensHalo(params,false){
//  sigma = 0.;
//  fratio = 0.;
//  pa = 0.;
//  rcore = 0.;
//
//  assignParams(params);
//
//  if(fratio  != 1.0) elliptical_flag = true;
//  else elliptical_flag = false;
//  ++objectCount;
//  if(objectCount == 1){   // make table for calculating elliptical integrale
//    construct_ellip_tables();
//  }
//  //Rsize = rmaxNSIE(sigma,mass,fratio,rcore);
//  //Rmax = MAX(1.0,1.0/fratio)*Rsize;  // redefine
//  LensHalo::setRsize(rmax_calc());
//  Rmax = Rmax_to_Rsize_ratio*LensHalo::getRsize();
//
//  if(fratio > 1.0 || fratio < 0.01) throw std::invalid_argument("invalid fratio");
//
//  if(rcore > 0.0){
//    LensHalo::setMass(MassBy1DIntegation(LensHalo::getRsize()) );
//  }
//
//   units = pow(sigma/lightspeed,2)/Grav;///sqrt(fratio); // mass/distance(physical)
//}

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
  
  
//  if(!params.get("main_stars_N",stars_N)) error_message1("main_stars_N",params.filename());
//  else if(stars_N){
//    assignParams_stars(params);
//  }
  
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

LensHaloTEPL::LensHaloTEPL(
              float my_mass
              ,PosType my_zlens
              ,PosType r_trunc
              ,PosType gamma
              ,float my_fratio
              ,float my_pa
              ,const COSMOLOGY &cosmo
              ,float f
):LensHalo(),tt(gamma),x_T(r_trunc),q(my_fratio),pa(my_pa)
{
  if(tt >= 2){
    std::cerr << "LensHaloTEPL : power-law index cannot be >= 2 because the mass will be unbounded." << std::endl;
    throw std::invalid_argument("bad gamma");
  }
 
  LensHalo::setMass(my_mass);
  LensHalo::setZlens(my_zlens,cosmo);
  if(q > 1) q = 1/q;
  
  if( q > 0.99 || q < 0.3){
    std::cerr << "LensHaloTEPL : axis ratio cannot be too large or 1 ." << std::endl;
    throw std::invalid_argument("bad gamma");
  }
  
  q_prime = (1-q*q)/q/q;
  
  R = std::complex<double>(cos(pa),sin(pa));
  
  SigmaT = my_mass * q * (2-tt) / (2 * PI * x_T * x_T);
  mass_pi = mass / PI;
  
  LensHalo::setRsize(x_T*f);
  Rmax = Rsize;
}

void LensHaloTEPL::force_halo(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType const *xcm,bool subtract_point,PosType screening){

  std::complex<double> z(xcm[0],xcm[1]);
  z = z * std::conj(R);

  double r2=std::norm(z);
  if(force_point(alpha,kappa,gamma,phi,xcm,r2
                 ,subtract_point,screening)) return;

  if(r2 < 1.0e-12){
    *kappa = SigmaT * pow(x_T/1.0e-6,tt);
    alpha[0] = alpha[1] = gamma[0] = gamma[1] = 0;
    return;
  }
  
  std::complex<double> a=0;
  std::complex<double> g=0;
  double k=0;
  
  deflection(z,a,g,k);

  a = - std::conj(a) * R;
  g = std::conj(g) * R * R;
  
  *kappa += k;
  alpha[0] += a.real();
  alpha[1] += a.imag();
  gamma[0] += g.real();
  gamma[1] += g.imag();
  
  assert(alpha[0] == alpha[0] && alpha[1] == alpha[1]);
  assert(gamma[0] == gamma[0] && gamma[1] == gamma[1]);
}

void LensHaloTEPL::deflection(std::complex<double> &z
                ,std::complex<double> &a
                ,std::complex<double> &g
                ,KappaType &sigma) const{

  double x_e = sqrt(q*q*z.real()*z.real() + z.imag()*z.imag());
  
  if(x_e <= 1.0e-5){
    sigma = SigmaT * pow(x_T/1.0e-5,tt);
    a = 0;
    g = 0;
    return;
  }
  if(x_e <= x_T){
    sigma = SigmaT * pow(x_T/x_e,tt);
    a = mass_pi * pow(x_T/x_e,tt-2) * F(x_e,tt,z) / z;
  }else{
    sigma = 0;
    a = mass_pi * F(x_T,tt,z) / z;
   }

  g =  (1 - tt) * a / z ;
  
  if(x_e <= x_T){
    g -= sigma * std::conj(z) / z;
  }else{
    g -= (2-tt) * mass_pi / sqrt(1. - q_prime * x_T*x_T / z / z ) / z / z;
  }
}

// for calculating exterior solution
std::complex<double> LensHaloTEPL::F(double r_e,double t,std::complex<double> z) const{

  std::complex<double> u = (1. - sqrt(1. - q_prime * r_e * r_e / z / z ) )/2.;
  //assert(std::norm(u) < 1);
  if(std::norm(u) > 1 ) u = u/std::norm(u);
  double a = 1;
  std::complex<double> sum(1,0);
  for(int n=1 ; n < 20 ; ++n){
    a *= 2*(n+1-t)/(2*n+2-t);
    sum += a * u;
    u *= u;
  }
  
  return sum;
}

LensHaloTEBPL::LensHaloTEBPL(
              float my_mass
              ,PosType my_zlens
              ,PosType r_break
              ,PosType r_trunc
              ,PosType t1
              ,PosType t2
              ,float my_fratio
              ,float my_pa
              ,const COSMOLOGY &cosmo
              ,float f
):LensHalo(),
q(my_fratio),rb(r_break),rt(r_trunc)
,m2(my_mass/(1 + (t1-t2)/(2-t1)*pow(rb/rt,2-t2)))
,m3(m2*pow(rb/rt,2-t2))
,m1(m2*(2-t2)/(2-t1)*pow(rb/rt,2-t2))
,h1(m1,my_zlens,rb,t1,q,0,cosmo)
,h2(m2,my_zlens,rt,t2,q,0,cosmo)
,h3(m3,my_zlens,rb,t2,q,0,cosmo)
{
  
  LensHalo::setMass(my_mass);
  LensHalo::setZlens(my_zlens,cosmo);
     
  R = std::complex<double>(cos(my_pa),sin(my_pa));
  
  Rmax = Rsize = f*rt;
}

void LensHaloTEBPL::force_halo(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType const *xcm,bool subtract_point,PosType screening){
  
  std::complex<PosType> z(xcm[0],xcm[1]);
   
  if(force_point(alpha,kappa,gamma,phi,xcm,std::norm(z)
                 ,subtract_point,screening)) return;
  
  z = z * std::conj(R);

  std::complex<PosType> a=0,g=0;
  std::complex<PosType> at,gt;
  double kappat;

  double xe2 = q*q*z.real()*z.real() + z.imag()*z.imag();
  if(xe2 > rb*rb){
    h3.deflection(z,at,gt,kappat);
    a -= at;
    g -= gt;
    h2.deflection(z,at,gt,kappat);
    a += at;
    g += gt;
    *kappa += kappat;
  }
  h1.deflection(z,at,gt,kappat);
  a += at;
  g += gt;
  *kappa += kappat;

  a = - std::conj(a) * R;
  g = std::conj(g) * R * R;
  
  alpha[0] += a.real();
  alpha[1] += a.imag();
  gamma[0] += g.real();
  gamma[1] += g.imag();
  
  assert(alpha[0] == alpha[0] && alpha[1] == alpha[1]);
  assert(gamma[0] == gamma[0] && gamma[1] == gamma[1]);

}

#ifdef ENABLE_CERF

LensHaloGaussian::LensHaloGaussian(
              float my_mass
              ,PosType my_zlens
              ,PosType r_scale
              ,float my_fratio
              ,float my_pa
              ,const COSMOLOGY &cosmo
              ,float f
):LensHalo(),q(my_fratio),pa(my_pa),I(0,1)
{

  I_sqpi = I / sqrt(PI);
  one_sqpi = 1.0 / sqrt(PI);
  
  LensHalo::setMass(my_mass);
  LensHalo::setZlens(my_zlens,cosmo);
  q = abs(q);
  if(q > 1){
    q = 1/q;
  }
  Rhight = r_scale*sqrt(q);
  q_prime = (1-q*q)/q/q;
  
  R = std::complex<double>(cos(pa),sin(pa));
  
  SigmaO = mass * q / (2* PI * Rhight * Rhight);
  
  ss = sqrt(2)*Rhight;
  if(q != 1.0){
    norm = SigmaO * sqrt(2*PI / q_prime) * Rhight / q;
    norm_g = norm /ss /sqrt(q_prime);  // norm x dzz/dz
  }else{
    // normalization
    norm = -mass / PI;
    norm_g = 0;
  }
  
  LensHalo::setRsize(Rhight*f);
  Rmax = Rsize;
}

void LensHaloGaussian::force_halo(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType const *xcm,bool subtract_point,PosType screening){
  
  std::complex<double> z(xcm[0],xcm[1]);
  double r2=std::norm(z);
  if(force_point(alpha,kappa,gamma,phi,xcm,r2
                 ,subtract_point,screening)) return;

  z = z * std::conj(R);
  
  std::complex<double> a=0;
  std::complex<double> g=0;
  double k=0;
  
  deflection(z,a,g,k);

  a = std::conj(a) * R;
  g = std::conj(g) * R * R;
  
  *kappa += k;
  alpha[0] += a.real();
  alpha[1] += a.imag();
  gamma[0] += g.real();
  gamma[1] += g.imag();
  
  assert(alpha[0] == alpha[0] && alpha[1] == alpha[1]);
  assert(gamma[0] == gamma[0] && gamma[1] == gamma[1]);
}

void LensHaloGaussian::deflection(std::complex<double> &z
                ,std::complex<double> &a
                ,std::complex<double> &g
                ,KappaType &sigma) const {

  std::complex<double> zz = z / ss / sqrt(q_prime);

  zz.real(abs(zz.real()));
  zz.imag(abs(zz.imag()));

  double x_e2 = (q*q*z.real()*z.real() + z.imag()*z.imag()) / ss / ss;
  double exp_x_e2 = exp(-x_e2);

  sigma = SigmaO * exp(- x_e2 );
  
  if(x_e2 < 1.0e-8){
    a = 0.0;
    g = 0.0;
    return;
  }
   
  if(q==1.0){
    a = norm / z * ( 1 - exp_x_e2 );
    g = a / z  - norm / z * exp_x_e2 * (z.real() - I*z.imag())/ss/ss;
  }else{
    
    std::complex<double> b = sqrt(x_e2 - zz*zz);
    
    std::complex<double> wfzz = wF(zz);
    std::complex<double> wfib = wF(I*b);
 

    // this is the correct result in all quadrants
    //a =  norm * sqrt(z*z)/ z
    //* ( sqrt(-zz2)/sqrt(zz2)*wF(I*sqrt(-zz2)) - b/(sqrt(zz2 - x_e2))*wF(I*b) * exp(- x_e2 ));
    a = norm * I * ( wfzz - wfib * exp_x_e2);
  
    std::complex<double> dx_e2dzz = (q*q*abs(z.real()) - I*abs(z.imag()))*sqrt(q_prime)/ss;
  
    g = -norm_g * I * ( 2.0*(I_sqpi - zz*wfzz)
                    + ( (1/sqrt(PI) - b*wfib) * (dx_e2dzz - 2.0*zz)/b
                    + wfib*dx_e2dzz ) * exp_x_e2
                    );

    if(z.real() < 0 && z.imag() < 0){
      a *= -1.0;
    }
    if(z.real() > 0 && z.imag() < 0){
      a = std::conj(a);
      g = std::conj(g);
    }
    if(z.real() < 0 && z.imag() > 0){
      a = -1.0*std::conj(a);
      g = std::conj(g);
    }
  }
  
}

std::complex<double> LensHaloGaussian::wF(std::complex<double> z) const{
  
  double _Complex w = w_of_z(reinterpret_cast<double _Complex(&)>(z));
  assert( w == w);
  return reinterpret_cast<std::complex<double>(&)>(w);
}
std::complex<double> LensHaloGaussian::my_erfc(std::complex<double> z) const{
  double _Complex w = cerfc(reinterpret_cast<double _Complex(&)>(z));
  assert( w == w);
  return reinterpret_cast<std::complex<double>(&)>(w);
}

#endif

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

void LensHalo::force_halo(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType const *xcm,bool subtract_point,PosType screening)
{
  PosType rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
 
  if(force_point(alpha,kappa,gamma,phi,xcm,rcm2
                 ,subtract_point,screening)) return;

  
  if (elliptical_flag){
    force_halo_asym(alpha,kappa,gamma,phi,xcm,subtract_point,screening);
    //assert(!isinf(*kappa) );
  }else{
    force_halo_sym(alpha,kappa,gamma,phi,xcm,subtract_point,screening);
    //assert(!isinf(*kappa) );
  }
}

/*
 Used in derived classes to subtract the point mass if necescary and take care of quantities beyond Rmax
 
 Returns true when no further calculation of the quantities should be done.
 
 This should be put in the beginning of all of the force_halo() functions.
 
 rcm2 needs to be calculated first
 
 */
bool LensHalo::force_point(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi
    ,PosType const *xcm,PosType rcm2
    ,bool subtract_point,PosType screening)
{
  
  if(rcm2 < Rmax*Rmax){
    if(rcm2==0) return false;
    if(subtract_point){
    
      PosType fac = screening*mass/rcm2/PI;
    
      alpha[0] += fac*xcm[0];
      alpha[1] += fac*xcm[1];
  
      fac = 2.0*fac/rcm2;
          
      gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*fac;
      gamma[1] += xcm[0]*xcm[1]*fac;
      
      *phi += -0.5*mass*log(rcm2) / PI;
    }
    return false;
  }
  if(!subtract_point){
    
    PosType prefac = mass/rcm2/PI;
    alpha[0] += -1.0*prefac*xcm[0];
    alpha[1] += -1.0*prefac*xcm[1];
      
    PosType tmp = -2.0*prefac/rcm2;
    gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
    gamma[1] += xcm[0]*xcm[1]*tmp;
    
    *phi += 0.5*mass*log(rcm2) / PI;
  }
//  else{
//    alpha[0] = alpha[1] = 0.0;
//    gamma[0] = gamma[1] = gamma[2] = 0.0;
//    *kappa = 0.0;
//    *phi = 0.0 ;
//  }
  
  return true;
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
//  if(stars_N > 0 && stars_implanted)
//  {
//    force_stars(alpha,kappa,gamma,xcm);
//  }
  
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
      if(main_ellip_method==EllipMethod::Pseudo){alphakappagamma_asym(LensHalo::getRsize(),theta, alpha_tmp,&kappa_tmp,gamma_tmp,&phi_tmp);}
      if(main_ellip_method==EllipMethod::Fourier){alphakappagamma1asym(LensHalo::getRsize(),theta, alpha_tmp,&kappa_tmp,gamma_tmp,&phi_tmp);}
      if(main_ellip_method==EllipMethod::Schramm){alphakappagamma2asym(LensHalo::getRsize(),theta, alpha_tmp,&kappa_tmp,gamma_tmp,&phi_tmp);}
      if(main_ellip_method==EllipMethod::Keeton){alphakappagamma3asym(LensHalo::getRsize(),theta, alpha_tmp,&kappa_tmp,gamma_tmp,&phi_tmp);}
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
      if(main_ellip_method==EllipMethod::Pseudo){alphakappagamma_asym(r,theta, alpha_tmp,&kappa_tmp,gamma_tmp,&phi_tmp);}
      if(main_ellip_method==EllipMethod::Fourier){alphakappagamma1asym(r,theta, alpha_tmp,&kappa_tmp,gamma_tmp,&phi_tmp);}
      if(main_ellip_method==EllipMethod::Schramm){alphakappagamma2asym(r,theta, alpha_tmp,&kappa_tmp,gamma_tmp,&phi_tmp);}
      if(main_ellip_method==EllipMethod::Keeton){alphakappagamma3asym(r,theta, alpha_tmp,&kappa_tmp,gamma_tmp,&phi_tmp);}
      
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
//  if(stars_N > 0 && stars_implanted)
//  {
//    force_stars(alpha,kappa,gamma,xcm);
//  }
  
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
 
  if(force_point(alpha,kappa,gamma,phi,xcm,rcm2
                 ,subtract_point,screening)) return;

  
  if(rcm2 < 1e-20) rcm2 = 1e-20;
  if(rcm2 < Rmax*Rmax){
    //PosType ellipR = ellipticRadiusNSIE(xcm,fratio,pa);
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
    
    
//    if(subtract_point)
//    {
//      PosType fac = screening*LensHalo::get_mass()/rcm2/PI;
//      alpha[0] += fac*xcm[0];
//      alpha[1] += fac*xcm[1];
//
//      {
//        fac = 2.0*fac/rcm2;
//
//        gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*fac;
//        gamma[1] += xcm[0]*xcm[1]*fac;
//      }
//    }
    
  }
//  else
//  {
//    // outside of the halo
//    if (subtract_point == false)
//    {
//      PosType prefac = LensHalo::get_mass()/rcm2/PI;
//      alpha[0] += -1.0*prefac*xcm[0];
//      alpha[1] += -1.0*prefac*xcm[1];
//
//      {
//        PosType tmp = -2.0*prefac/rcm2;
//
//        gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
//        gamma[1] += xcm[0]*xcm[1]*tmp;
//      }
//    }
//  }
  
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
LensHaloHernquist::LensHaloHernquist(float my_mass,float my_Rsize,PosType my_zlens,float my_rscale,float my_fratio,float my_pa,const COSMOLOGY &cosmo, EllipMethod my_ellip_method){
  
  LensHalo::setMass(my_mass);
  LensHalo::setZlens(my_zlens,cosmo);
  LensHalo::setRsize(my_Rsize);

  rscale=my_rscale;
  fratio=my_fratio, pa=my_pa;
  
  xmax = LensHalo::getRsize()/rscale;
  make_tables();
  gmax = InterpolateFromTable(gtable,xmax);
  
  set_slope(1);
  /// If the axis ratio given in the parameter file is set to 1 all ellipticizing routines are skipped.
  if(fratio!=1){
    Rmax = Rmax_to_Rsize_ratio*LensHalo::getRsize();
    //std::cout << getEllipMethod() << " method to ellipticise" << std::endl;
    if(getEllipMethod()==EllipMethod::Fourier){
      //std::cout << "Hernquist constructor: slope set to " << get_slope() << std::endl;
      calcModes(fratio, get_slope(), pa, mod); // to ellipticize potential instead of kappa use (fratio, get_slope()-2, pa, mod)
      for(int i=1;i<Nmod;i++){
        if(mod[i]!=0){set_flag_elliptical(true);};
      }
    }else set_flag_elliptical(true);
    if (getEllipMethod()==EllipMethod::Pseudo){
      fratio=0.00890632+0.99209115*pow(fratio,0.33697702);
      set_norm_factor();
    }
  }else{
    set_flag_elliptical(false);
    Rmax = LensHalo::getRsize();
  }
}

//LensHaloHernquist::LensHaloHernquist(InputParams& params): LensHalo(params)
//{
//  assignParams(params);
//  make_tables();
//  gmax = InterpolateFromTable(gtable,xmax);
//
//  set_slope(1);
//  /// If the axis ratio given in the parameter file is set to 1 all ellipticizing routines are skipped.
//  if(fratio!=1){
//    //std::cout << getEllipMethod() << " method to ellipticise" << std::endl;
//    if(getEllipMethod()==Fourier){
//      std::cout << "Hernquist constructor: slope set to " << get_slope() << std::endl;
//      calcModes(fratio, get_slope(), pa, mod); // to ellipticize potential instead of kappa use (fratio, get_slope()-2, pa, mod)
//      for(int i=1;i<Nmod;i++){
//        if(mod[i]!=0){set_flag_elliptical(true);};
//      }
//    }else set_flag_elliptical(true);
//    if (getEllipMethod()==Pseudo){
//      set_norm_factor();
//    }
//  }else set_flag_elliptical(false);
//}

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
  if(!params.get("main_ellip_method",main_ellip_method)){if(fratio!=1){main_ellip_method=EllipMethod::Pseudo;std::cout << "main_ellip_method is not defined in file " << params.filename() << ", hence set to Pseudo." << endl;};};
//  if(!params.get("main_stars_N",stars_N)) error_message1("main_stars_N",params.filename());
//  else if(stars_N){
//    assignParams_stars(params);
//
//    std::cout << "LensHalo::getRsize() " << LensHalo::getRsize() <<std::endl;
//  }
  
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
LensHaloJaffe::LensHaloJaffe(float my_mass,float my_Rsize,PosType my_zlens,float my_rscale,float my_fratio,float my_pa,const COSMOLOGY &cosmo, EllipMethod my_ellip_method){
  
  LensHalo::setMass(my_mass);
  LensHalo::setZlens(my_zlens,cosmo);
  LensHalo::setRsize(my_Rsize);

  rscale=my_rscale;
  fratio=my_fratio, pa=my_pa;
  xmax = LensHalo::getRsize()/rscale;
  make_tables();
  gmax = InterpolateFromTable(gtable,xmax);
  
  set_slope(1);
  if(fratio!=1){
    Rmax = Rmax_to_Rsize_ratio*LensHalo::getRsize();
    //std::cout << getEllipMethod() << " method to ellipticise" << std::endl;
    if(getEllipMethod()==EllipMethod::Fourier){
      //std::cout << "Jaffe constructor: slope set to " << get_slope() << std::endl;
      calcModes(fratio, get_slope(), pa, mod);
      for(int i=1;i<Nmod;i++){
        if(mod[i]!=0){set_flag_elliptical(true);};
      }
    }else set_flag_elliptical(true);
    if (getEllipMethod()==EllipMethod::Pseudo){
      fratio=0.00890632+0.99209115*pow(fratio,0.33697702);
      set_norm_factor();
    }
  }else{
    set_flag_elliptical(false);
    Rmax = LensHalo::getRsize();
  }
}

//LensHaloJaffe::LensHaloJaffe(InputParams& params): LensHalo(params)
//{
//  assignParams(params);
//  make_tables();
//  gmax = InterpolateFromTable(gtable,xmax);
//
//  set_slope(1);
//
//  /// If the axis ratio given in the parameter file is set to 1 all ellipticizing routines are skipped.
//  if(fratio!=1){
//    //std::cout << getEllipMethod() << " method to ellipticise" << std::endl;
//    if(getEllipMethod()==Fourier){
//      std::cout << "Jaffe constructor: slope set to " << get_slope() << std::endl;
//      calcModes(fratio, get_slope(), pa, mod);
//      for(int i=1;i<Nmod;i++){
//        if(mod[i]!=0){set_flag_elliptical(true);};
//      }
//    }else set_flag_elliptical(true);
//    if (getEllipMethod()==Pseudo){
//      set_norm_factor();
//    }
//  }else set_flag_elliptical(false);
//}

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
  
  if(!params.get("main_ellip_method",main_ellip_method)){if(fratio!=1){main_ellip_method=EllipMethod::Pseudo;std::cout << "main_ellip_method is not defined in file " << params.filename() << ", hence set to Pseudo." << endl;};};
  
//  if(!params.get("main_stars_N",stars_N)) error_message1("main_stars_N",params.filename());
//  else if(stars_N){
//    assignParams_stars(params);
//  }
  
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

LensHaloDummy::LensHaloDummy(float my_mass,float my_Rsize,PosType my_zlens,float my_rscale,const COSMOLOGY &cosmo){
  LensHalo::setMass(my_mass);
  LensHalo::setZlens(my_zlens,cosmo);
  LensHalo::setRsize(my_Rsize);

  rscale=my_rscale;
  setTheta(0.0,0.0);
}

//LensHaloDummy::LensHaloDummy(InputParams& params): LensHalo(params)
//{
//  assignParams(params);
//  //	mass = 0.;
//}

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
 
  if(force_point(alpha,kappa,gamma,phi,xcm,rcm2
                 ,subtract_point,screening)) return;

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
  std::cout << " LensHalo test :" << std::endl;
  std::cout << "test alpha's consistance with kappa by comparing mass interior to a radius by 1D integration and Gauss' law and by 2D integration" << std::endl << "  The total internal mass is " << mass << std::endl;
  
  std::cout << "R            R/Rmax     R/Rsize     Mass 1 D (from alpha)     Mass 2 D         (m1 - m2)/m1       m2/m1" << std::endl;
  
  int N=100;
  PosType m1,m2;
  for(int i=1;i<N;++i){
    m1 = MassBy1DIntegation(LensHalo::getRsize()*i/(N-4));
    m2 = MassBy2DIntegation(LensHalo::getRsize()*i/(N-4));
    std::cout << LensHalo::getRsize()*i*1./(N-4) << "      " <<  i*1./(N-4) << "      " << i*1./(N-4) << "      " << m1 << "       "
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

#ifdef ENABLE_EIGEN
#ifdef ENABLE_CERF

//int LensHaloMultiGauss::count = 0;

LensHaloMultiGauss::LensHaloMultiGauss(
                   double mass_norm
                   ,double Rnorm
                   ,MultiGauss::PROFILE &profile
                   ,int Ngaussians
                   ,int Nradii
                   ,PosType r_min
                   ,PosType r_max
                   ,PosType my_zlens
                   ,float my_fratio
                   ,float my_pa
                   ,const COSMOLOGY &cosmo
                   ,float f
                   ,bool verbose
):LensHalo(my_zlens,cosmo),nn(Ngaussians),mm(Nradii),q(my_fratio),pa(my_pa)
,mass_norm(mass_norm),r_norm(Rnorm)
{
 //++LensHaloMultiGauss::count;
 if(q <= 0){throw std::runtime_error("bad axis ratio"); }
 if(nn > mm){
    std::cerr << "LensHaloMultiGauss : nn must be less than mm." << std::endl;
    throw std::runtime_error("");
  }
  if( mass_norm <= 0){
     std::cerr << "LensHaloMultiGauss : mass_norm <= 0." << std::endl;
     throw std::runtime_error("");
   }

  if(q > 1){
    q = 1/q;
  }
  
  double totalmass=0;
  calc_masses_scales(profile,nn,mm,r_min,r_max,mass_norm,r_norm
                     ,totalmass,sigmas,A,rms_error,verbose);

  // construct Gaussian components
  for(int i=0 ; i<nn ; ++i){
    gaussians.emplace_back(A[i],my_zlens,sigmas[i],q,0,cosmo,f);
  }
  Rotation = std::complex<double>(cos(my_pa),sin(my_pa));
  Rmax = Rsize = f*r_max;
  
  totalmass = 0;
  for(int n=0 ; n<nn ; ++n){
    totalmass += A[n] * ( 1 - exp(-Rmax*Rmax / 2 / (sigmas[n]*sigmas[n]) ) );
  }
  LensHalo::setMass(totalmass);
}

LensHaloMultiGauss::LensHaloMultiGauss(
                   double my_mass_norm
                   ,double Rnorm
                   ,double my_scale   // radial scale in units of the scale that was used to produce relative_scales
                   ,const std::vector<double> &relative_scales
                   ,const std::vector<double> &relative_masses
                   ,PosType my_zlens /// redshift
                   ,float my_fratio /// axis ratio
                   ,float my_pa     /// position angle, 0 has long axis along the veritical axis and goes clockwise
                   ,const COSMOLOGY &cosmo  /// cosmology
                   ,float f        /// cuttoff radius in units of the larges scale
                   ,bool verbose
):LensHalo(my_zlens,cosmo),q(my_fratio),pa(my_pa),mass_norm(my_mass_norm),r_norm(Rnorm)

{
  //++LensHaloMultiGauss::count;
  
  if(q <= 0){throw std::runtime_error("bad axis ratio"); }
  if(relative_scales.size() != relative_masses.size()){
    throw std::runtime_error("arrays are wrong size.");
  }
  if( mass_norm <= 0){
     std::cerr << "LensHaloMultiGauss : mass_norm <= 0." << std::endl;
     throw std::runtime_error("");
   }

  nn = relative_scales.size();
  q = abs(q);
  if(q > 1){
    q = 1/q;
  }
  mm=0;
  
  sigmas = relative_scales;
  for(double &s : sigmas) s *= my_scale;

  
  A = relative_masses;
  double tmp_mass = 0;
  
  for(int n=0 ; n<nn ; ++n){
    tmp_mass += relative_masses[n] * ( 1 - exp(-r_norm*r_norm / 2 / (sigmas[n]*sigmas[n]) ) );
  }
  for(auto &a : A) a *= my_mass_norm/tmp_mass;

  Rotation = std::complex<double>(cos(my_pa),sin(my_pa));
  Rmax = Rsize = f*sigmas.back();

  double totalmass = 0;
  for(int n=0 ; n<nn ; ++n){
    totalmass += A[n] * ( 1 - exp(-Rmax*Rmax / 2 / (sigmas[n]*sigmas[n]) ) );
  }
  LensHalo::setMass(totalmass);
  
  // construct Gaussian components
  for(int i=0 ; i<nn ; ++i){
    gaussians.emplace_back(A[i],my_zlens,sigmas[i],q,0,cosmo
                           ,Rmax/sigmas[i]);
  }

  if(verbose){
    double mass=0;
    for(int n=0 ; n<nn ; ++n){
      mass += A[n] * ( 1 - exp(-r_norm*r_norm / 2 / (sigmas[n]*sigmas[n]) ) );
    }
    std::cout << " total mass at r_norm : " << mass << " mass / mass_in = " << mass/my_mass_norm << std::endl;
    
    mass=0;
    for(int n=0 ; n<nn ; ++n){
      mass += A[n] * ( 1 - exp(-Rmax*Rmax / 2 / (sigmas[n]*sigmas[n]) ) );
    }
    std::cout << " total mass at Rmax : " << mass << " mass / mass_in = " << mass/my_mass_norm << std::endl;

  }

}

LensHaloMultiGauss::LensHaloMultiGauss(LensHaloMultiGauss &&halo):
LensHalo(std::move(halo))
{
  //++LensHaloMultiGauss::count;
  nn=halo.nn; // number of gaussians
  mm=halo.mm; // number of fit radii
  q=halo.q; // axis ratio
  pa=halo.pa; // position angle
  mass_norm=halo.mass_norm;
  r_norm=halo.r_norm;
  sigmas=std::move(halo.sigmas);
  radii=std::move(halo.radii);
  Rotation=std::move(halo.Rotation);
  gaussians=std::move(halo.gaussians);
  A=std::move(halo.A);
  rms_error=halo.rms_error;
}
LensHaloMultiGauss::LensHaloMultiGauss(const LensHaloMultiGauss &halo):
LensHalo(halo)
{
  //++LensHaloMultiGauss::count;
  nn=halo.nn; // number of gaussians
  mm=halo.mm; // number of fit radii
  q=halo.q; // axis ratio
  pa=halo.pa; // position angle
  mass_norm=halo.mass_norm;
  r_norm=halo.r_norm;
  sigmas=halo.sigmas;
  radii=halo.radii;
  Rotation=halo.Rotation;
  gaussians=halo.gaussians;
  A=halo.A;
  rms_error=halo.rms_error;
}

LensHaloMultiGauss& LensHaloMultiGauss::operator=(const LensHaloMultiGauss &&halo){
 
  LensHalo::operator= (std::move(halo));
  nn=halo.nn; // number of gaussians
  mm=halo.mm; // number of fit radii
  q=halo.q; // axis ratio
  pa=halo.pa; // position angle
  mass_norm=halo.mass_norm;
  r_norm=halo.r_norm;
  sigmas=std::move(halo.sigmas);
  radii=std::move(halo.radii);
  Rotation=std::move(halo.Rotation);
  gaussians=std::move(halo.gaussians);
  A=std::move(halo.A);
  rms_error=halo.rms_error;
  
  return *this;
}
LensHaloMultiGauss& LensHaloMultiGauss::operator=(const LensHaloMultiGauss &halo){
  if(this==&halo) return *this;
  
  LensHalo::operator= (halo);
  nn=halo.nn; // number of gaussians
  mm=halo.mm; // number of fit radii
  q=halo.q; // axis ratio
  pa=halo.pa; // position angle
  mass_norm=halo.mass_norm;
  r_norm=halo.r_norm;
  sigmas=halo.sigmas;
  radii=halo.radii;
  Rotation=halo.Rotation;
  gaussians=halo.gaussians;
  A=halo.A;
  rms_error=halo.rms_error;
  
  return *this;
}

void LensHaloMultiGauss::set_pa(double my_pa){
  pa = my_pa;
  Rotation = std::complex<double>(cos(my_pa),sin(my_pa));
}

double LensHaloMultiGauss::profile(double r){
  double sigma=0;
  for(int n=0 ; n < nn ; ++n){
    sigma += A[n] * exp(-r*r / 2 / sigmas[n] /sigmas[n]) /2/PI/sigmas[n] /sigmas[n];
  }
  
  return sigma;
}

/// mass within elliptical radius
double LensHaloMultiGauss::mass_cum(double r){
  double mass_tmp=0;
  for(int n=0 ; n<nn ; ++n) mass_tmp += A[n]*( 1 - exp(-r*r / 2 / sigmas[n] /sigmas[n]) );
  
  return mass_tmp;
}

/// returns the RMS relative error for the fit points in the profile
float LensHaloMultiGauss::error(){
  return rms_error;
}

void LensHaloMultiGauss::force_halo(PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,PosType const *xcm,bool subtract_point,PosType screening){
  
  std::complex<PosType> z(xcm[0],xcm[1]);
  
  if(force_point(alpha,kappa,gamma,phi,xcm,std::norm(z)
                 ,subtract_point,screening)) return;
  
  z = z * std::conj(Rotation);
  
  std::complex<PosType> a=0,g=0;
  std::complex<PosType> at,gt;
  double kappat;
  
  for(auto &h : gaussians){
    h.deflection(z,at,gt,kappat);
    a += at;
    g += gt;
    *kappa += kappat;
  }
  
  a = std::conj(a) * Rotation;
  g = std::conj(g) * Rotation * Rotation;
  
  alpha[0] += a.real();
  alpha[1] += a.imag();
  gamma[0] += g.real();
  gamma[1] += g.imag();
  
  assert(alpha[0] == alpha[0] && alpha[1] == alpha[1]);
  assert(gamma[0] == gamma[0] && gamma[1] == gamma[1]);
}

#endif // eigen
#endif // libcerf

//bool LensHaloZcompare(LensHalo *lh1,LensHalo *lh2){return (lh1->getZlens() < lh2->getZlens());}
//bool compare(LensHalo *lh1,LensHalo *lh2){return (lh1->getZlens() < lh2->getZlens());}
