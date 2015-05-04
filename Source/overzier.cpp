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
: haloID(0), Reff(0), Rh(0), BtoT(0), PA(0), inclination(0),
  cxx(0), cyy(0), cxy(0), sbDo(0), sbSo(0), mag(0)
{
}

SourceOverzier::SourceOverzier(
		PosType my_mag              /// Total magnitude
		,double my_BtoT            /// Bulge to total ratio
		,double my_Reff         /// Bulge half light radius (arcs)
		,double my_Rh           /// disk scale hight (arcs)
		,double my_PA           /// Position angle (radians)
		,double my_inclination  /// inclination of disk (radians)
		,unsigned long my_id          ///          id number
		,double my_z            /// optional redshift
		,const double *my_theta          /// optional angular position on the sky
		){
	setInternals(my_mag,my_BtoT,my_Reff,my_Rh,my_PA,my_inclination,my_id,my_z,my_theta);
}

SourceOverzier::~SourceOverzier()
{
}

/// Sets internal variables.  If default constructor is used this must be called before the surface brightness function.
void SourceOverzier::setInternals(double my_mag,double my_BtoT,double my_Reff,double my_Rh,double my_PA,double incl,unsigned long my_id,double my_z,const double *my_theta){

	haloID = my_id;

	Reff = my_Reff*pi/180/60/60;
	Rh = my_Rh*pi/180/60/60;
	mag = my_mag;
	BtoT = my_BtoT;
	PA = my_PA;
	inclination = incl;

	if(Rh > 0.0){
		cxx = ( pow(cos(PA),2) + pow(sin(PA)/cos(incl),2) )/Rh/Rh;
		cyy = ( pow(sin(PA),2) + pow(cos(PA)/cos(incl),2) )/Rh/Rh;
		cxy = ( 2*cos(PA)*sin(PA)*(1-pow(1/cos(incl),2)) )/Rh/Rh;
	}else{
		cxx = cyy = cxy = 0.0;
	}

	//muDo = my_mag-2.5*log10(1-BtoT)+5*log10(Rh)+1.9955;
	//muSo = my_mag-2.5*log10(BtoT)+5*log10(Reff)-4.9384;

	if(Rh > 0.0) sbDo = pow(10,-my_mag/2.5)*0.159148*(1-BtoT)/pow(Rh,2);
	else sbDo = 0.0;
	if(Reff > 0.0) sbSo = pow(10,-my_mag/2.5)*94.484376*BtoT/pow(Reff,2);
	else sbSo = 0.0;
	
	// redshift
	setZ(my_z);
	
	// radius
	// weighted mean between the radii that enclose 99% of the flux
	// in the pure De Vacouleur/exponential disk case
	// 6.670 = 3.975*Re = 3.975*1.678*Rh
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

PosType SourceOverzier::getTotalFlux(){
	return pow(10,-(48.6+mag)/2.5);
}

void SourceOverzier::printSource(){
	std::cout << "bulge half light radius: " << Reff << " arcs   disk scale hight: " << Rh << " arcs" << std::endl;
}

void SourceOverzier::assignParams(InputParams& /* params */)
{
}

SourceOverzierPlus::SourceOverzierPlus(PosType mag,PosType BtoT,PosType Reff,PosType Rh,PosType PA,PosType inclination,unsigned long my_id,PosType my_z,const PosType *theta,Utilities::RandomNumbers_NR &ran):
SourceOverzier(mag,BtoT,Reff,Rh,PA,inclination,my_id,my_z,theta)
{
  
  float minA = 0.01,maxA = 0.5; // minimum and maximum amplitude of arms
  Narms = 4;  // number of arms
  arm_alpha = (21 + 10*(ran()-0.5)*2)*degreesTOradians; // arm pitch angle
  mctalpha = Narms/tan(arm_alpha);
  disk_phase = pi*ran(); // add phase of arms
  Ad = minA + (maxA-minA)*ran();
  
  // spheroid
  
  // extra cersic component
  //double index = 4 + 3*(ran()-0.5)*2;

  
  double index = 4*pow(MAX(BtoT,0.03),0.4)*pow(10,0.2*(ran()-0.5));
  double q = 1 + (0.7-1)*ran();
  
  //double q = 0.5;
  double muSo = mag-2.5*log10(BtoT);  // this is not correct!
  
  //spheroid.reset(new SourceSersic(muSo,Reff,-PA + 10*(ran() - 0.5)*pi/180,index,q,my_z,theta));
  spheroid = new SourceSersic(muSo,Reff,-PA + 10*(ran() - 0.5)*pi/180,index,q,my_z,theta);
  
  cospa = cos(PA);
  sinpa = sin(PA);
  cosi  = cos(inclination);
  
  modes.resize(6);
  for(PosType mod : modes){
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
  
  Point_2d z;
  
  z[0] = cospa*x[0] - sinpa*x[1];
  z[1] = ( sinpa*x[0] + cospa*x[1] )/cosi;
  
  //PosType R = sqrt( cxx*x[0]*x[0] + cyy*x[1]*x[1] + cxy*x[0]*x[1] );
  PosType R = z.length()/Rh;
  PosType theta = atan2(z[1],z[0]);
  
  PosType disk_sb,bulge_sb;
  
  //disk
  disk_sb = sbDo*exp(-R);
  //spiral arms
  PosType phir = mctalpha*log(R) + disk_phase;
  disk_sb *= 1 + Ad*cos(Narms*theta + phir);
  
  
  PosType sb = disk_sb;
  sb *= pow(10,-0.4*48.6)*inv_hplanck;
  
  // bulge
  if(Reff > 0.0){
    // bulge perturbations
    PosType perturb = 1.0,tmp;
    int N = modes.size()/2;
    PosType c1 = x[0]/x.length();
    PosType s1 = x[1]/x.length();
    PosType cn=1,sn=0;
    for(int n=0;n<N;++n){
      tmp = cn;
      cn = (cn*c1 - sn*s1);
      sn = (sn*c1 + tmp*s1);
      perturb += modes[2*n]*cn + modes[2*n+1]*sn;
    }
    // spheroid contribution
    sb += spheroid->SurfaceBrightness(x.x)*perturb;
   }else{
    bulge_sb = 0.0;
  }

  if(sb< sb_limit)
    return 0.;
  
  return sb;
}

void SourceOverzierPlus::setBand(Band band){
  switch(band){
    case SDSS_U:
      mag = mag_u;
      break;
    case SDSS_G:
      mag = mag_g;
      break;
    case SDSS_R:
      mag = mag_r;
      break;
    case SDSS_I:
      mag = mag_i;
      break;
    case SDSS_Z:
      mag = mag_z;
      break;
    case J:
      mag = mag_J;
      break;
    case Ks:
      mag = mag_Ks;
      break;
    default:
      throw std::invalid_argument("band not supported");
      break;
  }
  if(Rh > 0.0) sbDo = pow(10,(-mag+oldmag)/2.5);
  if(Reff > 0.0){
    spheroid->setMag(mag-2.5*log10(BtoT));
  }
}



