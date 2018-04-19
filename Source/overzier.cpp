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
		PosType my_mag          /// Total magnitude
		,double my_mag_bulge    /// Bulge to total ratio
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

SourceOverzier::SourceOverzier(const SourceOverzier &s)
:Source(s){
/// bulge half light radius
  Reff = s.Rh;
  Rh = s.Rh;
  PA = s.PA;
  inclination = s.inclination;

  cxx = s.cxx;
  cyy = s.cyy;
  cxy = s.cxy;
  sbDo = s.sbDo;
  sbSo = s.sbSo;
  mag = s.mag;
  mag_bulge = s.mag_bulge;
  sedtype = s.sedtype;

  mag_map = s.mag_map;
  bulge_mag_map = s.bulge_mag_map;
}
SourceOverzier& SourceOverzier::operator=(const SourceOverzier &s){
  if(this == &s) return *this;
  
  Source::operator=(s);
  /// bulge half light radius
  Reff = s.Reff;
  Rh = s.Rh;
  PA = s.PA;
  inclination = s.inclination;
  
  cxx = s.cxx;
  cyy = s.cyy;
  cxy = s.cxy;
  sbDo = s.sbDo;
  sbSo = s.sbSo;
  mag = s.mag;
  mag_bulge = s.mag_bulge;
  sedtype = s.sedtype;

  mag_map = s.mag_map;
  bulge_mag_map = s.bulge_mag_map;

  return *this;
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
		setTheta(my_theta[0], my_theta[1]);
	else
		setTheta(0, 0);
}

/// Surface brightness in erg/cm^2/sec/rad^2/Hz
PosType SourceOverzier::SurfaceBrightness(
		PosType *y  /// position in radians
		){
	// position relative to center
	PosType x[2];
	x[0] = y[0]-getTheta()[0];
	x[1] = y[1]-getTheta()[1];
	
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
	std::cout << "bulge half light radius: " << Reff*180*60*60/pi << " arcs   disk scale hight: " << Rh*180*60*60/pi << " arcs" << std::endl;
}

void SourceOverzier::assignParams(InputParams& /* params */)
{
}

PosType SourceOverzier::getMag(Band band) const {
  return mag_map.at(band);
}

PosType SourceOverzier::getMagBulge(Band band) const {
  return bulge_mag_map.at(band);
}

void SourceOverzier::setMag(Band band,PosType my_mag){
  mag_map[band] = my_mag;
 }

void SourceOverzier::setMagBulge(Band band,PosType my_mag){
  bulge_mag_map[band] = my_mag;
 }

 void SourceOverzier::changeBand(Band band){
   
   mag = mag_map[band];
   mag_bulge = bulge_mag_map[band];
  
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
  
  //double index = 4*pow(MAX(getBtoT(),0.03),0.4)*pow(10,0.2*(ran()-0.5));
  double index = ran() + 3.5 ;

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
                              ,p.spheroid->getZ(),p.spheroid->getTheta().x);
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
                              ,p.spheroid->getZ(),p.spheroid->getTheta().x);
  modes = p.modes;
  
  return *this;
}


PosType SourceOverzierPlus::SurfaceBrightness(PosType *y){
  // position relative to center
  Point_2d x;
  x[0] = y[0]-getTheta()[0];
  x[1] = y[1]-getTheta()[1];
  
  PosType xlength = x.length();
  
  Point_2d z;
  PosType sb;

  if(xlength > 0 && sbDo > 0){
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

void SourceOverzierPlus::changeBand(Band band){

  SourceOverzier::changeBand(band);
  
  //mag = mag_map.at(band);
  //mag_bulge = bulge_mag_map.at(band);
  //SourceOverzier::renormalize();
  if(Reff > 0.0){
    spheroid->setMag(mag_bulge);
  }
}

void SourceOverzierPlus::randomize(Utilities::RandomNumbers_NR &ran){
  
  float BtoT;
  { // SourceOverzier variables
    
    Reff *= (1 + 0.2*(2*ran()-1.));
    Rh *= (1 + 0.2*(2*ran()-1.));
    
    PosType tmp = 0.01*(2*ran()-1.);
    
    for(auto mag = mag_map.begin() ; mag != mag_map.end() ; ++mag){
      mag->second = mag->second + tmp;
    }
    
    for(auto mag = bulge_mag_map.begin() ; mag != bulge_mag_map.end() ; ++mag){
      mag->second = mag->second + tmp;
    }

    /*
     setUMag(getMag(SDSS_U) + tmp);
     setGMag(getMag(SDSS_G) + tmp);
     setRMag(getMag(SDSS_R) + tmp);
     setIMag(getMag(SDSS_I) + tmp);
     setZMag(getMag(SDSS_Z) + tmp);
     setJMag(getMag(J) + tmp);
     setKMag(getMag(Ks) + tmp);
     */
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
  spheroid = new SourceSersic(mag_bulge,Reff/arcsecTOradians,-PA + 10*(ran() - 0.5)*pi/180,index,q,zsource,getTheta().x);
  
  
  for(PosType &mod : modes){
    mod = 5.0e-2*ran();
  }
  
}


