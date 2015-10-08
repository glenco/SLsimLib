/*
 * Source/overzier.cpp
 *
 *  Created on: Apr 12, 2012
 *      Author: bmetcalf
 *
 *      This surface brightness model comes from Roderik Overzier in private communication although it
 *      is probably written into a paper somewhere.
 */
#include "slsimlib.h"

/// sets everything to zero
SourceOverzier::SourceOverzier()
: haloID(0), Reff(0), Rh(0),  PA(0), inclination(0),
  cxx(0), cyy(0), cxy(0), sbDo(0), sbSo(0), mag(0), mag_bulge(0)
{
}

SourceOverzier::SourceOverzier(
		PosType my_mag              /// Total magnitude
		,double my_mag_bulge            /// Bulge to total ratio
		,double my_Reff         /// Bulge half light radius (arcs)
		,double my_Rh           /// disk scale hight (arcs)
		,double my_PA           /// Position angle (radians)
		,double my_inclination  /// inclination of disk (radians)
		,unsigned long my_id          ///          id number
		,double my_z            /// optional redshift
		,const double *my_theta          /// optional angular position on the sky
		){
	setInternals(my_mag,my_mag_bulge,my_Reff,my_Rh,my_PA,my_inclination,my_id,my_z,my_theta);
}

SourceOverzier::~SourceOverzier()
{
}

/// Sets internal variables.  If default constructor is used this must be called before the surface brightness function.
void SourceOverzier::setInternals(double my_mag,double my_mag_bulge,double my_Reff,double my_Rh,double my_PA,double incl,unsigned long my_id,double my_z,const double *my_theta){

	haloID = my_id;

	Reff = my_Reff*arcsecTOradians;
	Rh = my_Rh*arcsecTOradians;
	mag = my_mag;
  mag_bulge = my_mag_bulge;
	PA = my_PA;
	inclination = incl;

	if(Rh > 0.0){
		cxx = ( pow(cos(PA),2) + pow(sin(PA)/cos(incl),2) )/Rh/Rh;
		cyy = ( pow(sin(PA),2) + pow(cos(PA)/cos(incl),2) )/Rh/Rh;
		cxy = ( 2*cos(PA)*sin(PA)*(1-pow(1/cos(incl),2)) )/Rh/Rh;
	}else{
		cxx = cyy = cxy = 0.0;
	}
	
  renormalize();
  
	// redshift
	setZ(my_z);
	
	// radius
	// weighted mean between the radii that enclose 99% of the flux
	// in the pure De Vacouleur/exponential disk case
	// 6.670 = 3.975*Re = 3.975*1.678*Rh
  float BtoT = getBtoT();
	setRadius(6.670*Rh*(1-BtoT)+18.936*Reff*BtoT);
	
	// position
	if(my_theta != NULL)
		setX(my_theta[0], my_theta[1]);
	else
		setX(0, 0);
}

/// Surface brightness in erg/cm^2/sec/rad^2/Hz
PosType SourceOverzier::SurfaceBrightness(
		PosType *y  /// position in radians
		){
	// position relative to center
	PosType x[2];
	x[0] = y[0]-getX()[0];
	x[1] = y[1]-getX()[1];
	
	PosType R = cxx*x[0]*x[0] + cyy*x[1]*x[1] + cxy*x[0]*x[1],sb;
	R = sqrt(R);

	//sb = sbDo*exp(-(R)) + sbSo*exp(-7.6693*pow(R/Reff,0.25));
	sb = sbDo*exp(-R);
      
	if(Reff > 0.0) sb += sbSo*exp(-7.6693*pow((x[0]*x[0] + x[1]*x[1])/Reff/Reff,0.125));
  //	if(sb < 1.0e-4*(sbDo + sbSo) ) return 0.0;
	sb *= pow(10,-0.4*48.6)*inv_hplanck;
	
	if(sb< sb_limit)
		return 0.;
	
	return sb;
}

PosType SourceOverzier::getTotalFlux() const{
	return pow(10,-(48.6+mag)/2.5);
}

void SourceOverzier::printSource(){
	std::cout << "bulge half light radius: " << Reff << " arcs   disk scale hight: " << Rh << " arcs" << std::endl;
}

void SourceOverzier::assignParams(InputParams& /* params */)
{
}

PosType SourceOverzier::getMag(Band band) const {
  
  switch(band){
    case SDSS_U:
      return mag_u;
    case SDSS_G:
      return mag_g;
    case SDSS_R:
      return mag_r;
    case SDSS_I:
      return mag_i;
    case SDSS_Z:
      return mag_z;
    case J:
      return mag_J;
    case Ks:
      return mag_Ks;
      
    default:
      throw std::invalid_argument("band not supported");
      return 0.0;
      break;
  }
}
PosType SourceOverzier::getMagBulge(Band band) const {
  
  switch(band){
    case SDSS_U:
      return mag_u_bulge;
    case SDSS_G:
      return mag_g_bulge;
    case SDSS_R:
      return mag_r_bulge;
    case SDSS_I:
      return mag_i_bulge;
    case SDSS_Z:
      return mag_z_bulge;
    case J:
      return mag_J_bulge;
    case Ks:
      return mag_Ks_bulge;
      
    default:
      throw std::invalid_argument("band not supported");
      return 0.0;
      break;
  }
}

void SourceOverzier::setMag(Band band,PosType my_mag){
  switch(band){
    case SDSS_U:
      mag_u = my_mag;
      break;
    case SDSS_G:
      mag_g = my_mag;
      break;
    case SDSS_R:
      mag_r = my_mag;
      break;
    case SDSS_I:
      mag_i = my_mag;
      break;
    case SDSS_Z:
      mag_z = my_mag;
      break;
    case J:
      mag_J = my_mag;
      break;
    case Ks:
      mag_Ks = my_mag;
      break;
    case H:
      mag_H = my_mag;
      break;
    default:
      throw std::invalid_argument("band not supported");
      break;
  }
}

void SourceOverzier::setMagBulge(Band band,PosType my_mag){
  switch(band){
    case SDSS_U:
      mag_u_bulge = my_mag;
      break;
    case SDSS_G:
      mag_g_bulge = my_mag;
      break;
    case SDSS_R:
      mag_r_bulge = my_mag;
      break;
    case SDSS_I:
      mag_i_bulge = my_mag;
      break;
    case SDSS_Z:
      mag_z_bulge = my_mag;
      break;
    case J:
      mag_J_bulge = my_mag;
      break;
    case Ks:
      mag_Ks_bulge = my_mag;
      break;
    case H:
      mag_H_bulge = my_mag;
      break;
    default:
      throw std::invalid_argument("band not supported");
      break;
  }
}

 void SourceOverzier::changeBand(Band band){
  
  switch(band){
    case SDSS_U:
      mag = mag_u;
      mag_bulge = mag_u_bulge;
      break;
    case SDSS_G:
      mag = mag_g;
      mag_bulge = mag_g_bulge;
      break;
    case SDSS_R:
      mag = mag_r;
      mag_bulge = mag_r_bulge;
      break;
    case SDSS_I:
      mag = mag_i;
      mag_bulge = mag_i_bulge;
      break;
    case SDSS_Z:
      mag = mag_z;
      mag_bulge = mag_z_bulge;
      break;
    case J:
      mag = mag_J;
      mag_bulge = mag_J_bulge;
      break;
    case Ks:
      mag = mag_Ks;
      mag_bulge = mag_Ks_bulge;
      break;
    case H:
      mag = mag_H;
      mag_bulge = mag_H_bulge;
      break;
    default:
      throw std::invalid_argument("band not supported");
      break;
  }
  
  renormalize();
}

void SourceOverzier::renormalize(){
  float BtoT = getBtoT();
  if(Rh > 0.0) sbDo = pow(10,-mag/2.5)*0.159148*(1-BtoT)/pow(Rh,2);
  else sbDo = 0.0;
  if(Reff > 0.0) sbSo = pow(10,-mag/2.5)*94.484376*BtoT/pow(Reff,2);
  else sbSo = 0.0;
}


SourceOverzierPlus::SourceOverzierPlus(PosType my_mag,PosType my_mag_bulge,PosType my_Reff,PosType my_Rh,PosType my_PA,PosType inclination,unsigned long my_id,PosType my_z,const PosType *theta,Utilities::RandomNumbers_NR &ran):
SourceOverzier(my_mag,my_mag_bulge,my_Reff,my_Rh,my_PA,inclination,my_id,my_z,theta)
{
  
  float minA = 0.01,maxA = 0.2; // minimum and maximum amplitude of arms
  Narms = (ran() > 0.2) ? 2 : 4;  // number of arms
  arm_alpha = (21 + 10*(ran()-0.5)*2)*degreesTOradians; // arm pitch angle
  mctalpha = Narms/tan(arm_alpha);
  disk_phase = pi*ran(); // add phase of arms
  Ad = minA + (maxA-minA)*ran();
  
  // extra cersic component
  //double index = 4 + 3*(ran()-0.5)*2;
  
  double index = 4*pow(MAX(getBtoT(),0.03),0.4)*pow(10,0.2*(ran()-0.5));
  
  double q = 1 + (0.5-1)*ran();
  
  spheroid = new SourceSersic(my_mag_bulge,my_Reff,-my_PA + 10*(ran() - 0.5)*pi/180,index,q,my_z,theta);
  
  cospa = cos(PA);
  sinpa = sin(PA);
  cosi  = cos(inclination);
  
  modes.resize(6);
  for(PosType &mod : modes){
    mod = 2.0e-2*ran();
  }

}

SourceOverzierPlus::~SourceOverzierPlus(){
  delete spheroid;
}

SourceOverzierPlus::SourceOverzierPlus(const SourceOverzierPlus &p):
SourceOverzier(p),
Narms(p.Narms),Ad(p.Ad),mctalpha(p.mctalpha),arm_alpha(p.arm_alpha)
,disk_phase(p.disk_phase),cospa(p.cospa),sinpa(p.sinpa),cosi(p.cosi)
{
  spheroid = new SourceSersic(p.spheroid->getMag(),p.getReff(),p.spheroid->getPA()
                              ,p.spheroid->getSersicIndex(),p.spheroid->getAxesRatio()
                              ,p.spheroid->getZ(),p.spheroid->getX());
  modes = p.modes;
}

SourceOverzierPlus & SourceOverzierPlus::operator=(const SourceOverzierPlus &p){
  if(this == &p) return *this;
  SourceOverzier::operator=(p);
  
  Narms=p.Narms;
  Ad=p.Ad;
  mctalpha=p.mctalpha;
  arm_alpha=p.arm_alpha;
  disk_phase=p.disk_phase;
  cospa=p.cospa;
  sinpa=p.sinpa;
  cosi=p.cosi;
  
  spheroid = new SourceSersic(p.spheroid->getMag(),p.getReff(),p.spheroid->getPA()
                              ,p.spheroid->getSersicIndex(),p.spheroid->getAxesRatio()
                              ,p.spheroid->getZ(),p.spheroid->getX());
  modes = p.modes;
  
  return *this;
}


PosType SourceOverzierPlus::SurfaceBrightness(PosType *y){
  // position relative to center
  Point_2d x;
  x[0] = y[0]-getX()[0];
  x[1] = y[1]-getX()[1];
  
  PosType xlength = x.length();
  
  Point_2d z;
  PosType sb;

  if(xlength > 0){
    z[0] = cospa*x[0] - sinpa*x[1];
    z[1] = ( sinpa*x[0] + cospa*x[1] )/cosi;
    
    //PosType R = sqrt( cxx*x[0]*x[0] + cyy*x[1]*x[1] + cxy*x[0]*x[1] );
    PosType R = z.length()/Rh;
    PosType theta = atan2(z[1],z[0]);
    
    //disk
    sb = sbDo*exp(-R);
    //spiral arms
    PosType phir = mctalpha*log(R/0.5) + disk_phase;
    //PosType phiro = 0.4*mctalpha*log(R/0.5) + disk_phase;
    //if(R > 0.5 )
    sb *= 1 + Ad*cos(Narms*theta + phir);
    //else disk_sb *= 1 + Ad*cos(Narms*theta + phiro);
  }else{
    sb = sbDo;
  }
  
  sb *= pow(10,-0.4*48.6)*inv_hplanck;
  
  //std::cout << "disk sb " << sb << std::endl;
  
  // bulge
  if(Reff > 0.0){
    PosType perturb = 1.0,tmp;
    if(xlength > 0){
      // bulge perturbations
      int N = modes.size()/2;
      PosType c1 = x[0]/xlength;
      PosType s1 = x[1]/xlength;
      PosType cn=1,sn=0;
      for(int n=0;n<N;++n){
        tmp = cn;
        cn = (cn*c1 - sn*s1);
        sn = (sn*c1 + tmp*s1);
        perturb += modes[2*n]*cn + modes[2*n+1]*sn;
      }
    }
    // spheroid contribution
    sb += spheroid->SurfaceBrightness(y)*perturb;
    
    
    // !!!
    //sb = sbSo*exp(-7.6693*pow((x[0]*x[0] + x[1]*x[1])/Reff/Reff,0.125))
    //*pow(10,-0.4*48.6)*inv_hplanck;

   }
  
  //std::cout << "total sb " << sb << std::endl;

  if(sb< sb_limit)  return 0.;
  
  assert(sb >= 0.0);
  return sb;
}

void SourceOverzierPlus::setBand(Band band){

  switch(band){
    case SDSS_U:
      mag = mag_u;
      mag_bulge = mag_u_bulge;
      break;
    case SDSS_G:
      mag = mag_g;
      mag_bulge = mag_g_bulge;
      break;
    case SDSS_R:
      mag = mag_r;
      mag_bulge = mag_r_bulge;
      break;
    case SDSS_I:
      mag = mag_i;
      mag_bulge = mag_i_bulge;
      break;
    case SDSS_Z:
      mag = mag_z;
      mag_bulge = mag_z_bulge;
      break;
    case J:
      mag = mag_J;
      mag_bulge = mag_J_bulge;
      break;
    case Ks:
      mag = mag_Ks;
      mag_bulge = mag_Ks_bulge;
      break;
    default:
      throw std::invalid_argument("band not supported");
      break;
  }
  
  SourceOverzier::renormalize();
  if(Reff > 0.0){
    spheroid->setMag(mag_bulge);
  }
}

void SourceOverzierPlus::randomize(Utilities::RandomNumbers_NR &ran){
  
  float BtoT;
  { // SourceOverzier variables
    
    Reff *= (1 + 0.2*(2*ran()-1.));
    Rh *= (1 + 0.2*(2*ran()-1.));
    
    PosType tmp = 0.1*(2*ran()-1.);
    
    setUMag(getMag(SDSS_U) + tmp);
     setGMag(getMag(SDSS_G) + tmp);
     setRMag(getMag(SDSS_R) + tmp);
     setIMag(getMag(SDSS_I) + tmp);
     setZMag(getMag(SDSS_Z) + tmp);
     setJMag(getMag(J) + tmp);
     setKMag(getMag(Ks) + tmp);
     mag += tmp;
     
    
    PA = pi*ran();
    inclination = 0.9*pi/2*ran();
    
    cospa = cos(PA);
    sinpa = sin(PA);
    cosi  = cos(inclination);

    if(Rh > 0.0){
      cxx = ( cospa*cospa + pow(sinpa/cosi,2) )/Rh/Rh;
      cyy = ( sinpa*sinpa + pow(cospa/cosi,2) )/Rh/Rh;
      cxy = ( 2*cospa*sinpa*(1-1./cosi/cosi) )/Rh/Rh;
    }else{
      cxx = cyy = cxy = 0.0;
    }
    
    SourceOverzier::renormalize();
    
    // radius
    // weighted mean between the radii that enclose 99% of the flux
    // in the pure De Vacouleur/exponential disk case
    // 6.670 = 3.975*Re = 3.975*1.678*Rh
    BtoT = getBtoT();
    setRadius(6.670*Rh*(1-BtoT)+18.936*Reff*BtoT);
  }
  
  float minA = 0.01,maxA = 0.2; // minimum and maximum amplitude of arms
  Narms = (ran() > 0.2) ? 2 : 4;  // number of arms
  arm_alpha = (21 + 10*(ran()-0.5)*2)*degreesTOradians; // arm pitch angle
  mctalpha = Narms/tan(arm_alpha);
  disk_phase = pi*ran(); // add phase of arms
  Ad = minA + (maxA-minA)*ran();
 
  // spheroid
  
  // extra cersic component
  //double index = 4 + 3*(ran()-0.5)*2;
  
  double index = 4*pow(MAX(BtoT,0.03),0.4)*pow(10,0.2*(ran()-0.5));
  
  double q = 1 + (0.5-1)*ran();
  
  delete spheroid;
  spheroid = new SourceSersic(mag_bulge,Reff/arcsecTOradians,-PA + 10*(ran() - 0.5)*pi/180,index,q,zsource,getX());
  
  
  for(PosType &mod : modes){
    mod = 5.0e-2*ran();
  }
  
}


