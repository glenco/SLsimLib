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


SourceOverzier::SourceOverzier(
		PosType my_mag          /// Total magnitude
		,double my_mag_bulge    /// magnitude of bulge
		,double my_Reff         /// Bulge half light radius (arcs)
		,double my_Rdisk        /// disk scale hight (arcs)
		,double my_PA           /// Position angle (radians)
		,double my_inclination  /// inclination of disk (radians)
		,unsigned long my_id          ///          id number
		,double my_z            /// optional redshift
		,const double *my_theta          /// optional angular position on the sky
    ,double zeropoint       /// magnitude zero point

		):Source(0,Point_2d(0,0),my_z,-1,zeropoint),
 spheroid(my_mag_bulge,my_Reff,0,4,(1-0.5*sin(my_inclination)),my_z,zeropoint)
{

      //std::cout << "SourceOverzier constructor" << std::endl;
  setInternals(my_mag,my_mag_bulge,my_Reff
               ,my_Rdisk,my_PA,my_inclination
               ,my_id,my_z,my_theta);

  assert(current.Reff != 0 || current.Rdisk !=0 );
  spheroid.setTheta(0,0);
  
  // position
  if(my_theta != NULL)
    setTheta(my_theta[0], my_theta[1]);
  else
    setTheta(0, 0);
}

SourceOverzier::~SourceOverzier()
{
}

SourceOverzier::SourceOverzier(const SourceOverzier &s)
:Source(s),spheroid(s.spheroid){
  //spheroid = s.spheroid;
  current = s.current;
  sedtype = s.sedtype;
}
SourceOverzier& SourceOverzier::operator=(const SourceOverzier &s){
  if(this == &s) return *this;
  
  Source::operator=(s);
  
  spheroid = s.spheroid;
  current = s.current;
  sedtype = s.sedtype;
  return *this;
}

/// Sets internal variables.  If default constructor is used this must be called before the surface brightness function.
void SourceOverzier::setInternals(double my_mag,double my_mag_bulge,double my_Reff,double my_Rdisk,double my_PA,double incl,unsigned long my_id,double my_z,const double *my_theta){

	haloID = my_id;

  current.cosPA = cos(my_PA);
  current.sinPA = sin(my_PA);
  
  if(my_Reff < 0.05) my_Reff = 0.05; // ????
	current.Reff = my_Reff*arcsecTOradians;
	current.Rdisk = my_Rdisk*arcsecTOradians;
	current.mag = my_mag;
  current.mag_bulge = my_mag_bulge;
	current.PA = my_PA;
	current.inclination = incl;

	if(current.Rdisk > 0.0){
		//current.cxx = ( pow(cos(current.PA),2) + pow(sin(current.PA)/cos(incl),2) )/current.Rdisk/current.Rdisk;
		//current.cyy = ( pow(sin(current.PA),2) + pow(cos(current.PA)/cos(incl),2) )/current.Rdisk/current.Rdisk;
		//current.cxy = ( 2*cos(current.PA)*sin(current.PA)*(1-pow(1/cos(incl),2)) )/current.Rdisk/current.Rdisk;
    current.cxx = 1.0/current.Rdisk/current.Rdisk;
    current.cyy = 1.0/cos(incl)/cos(incl)/current.Rdisk/current.Rdisk;

	}
	
  renormalize_current();
  
	// redshift
	setZ(my_z);
	
	// radius
	// weighted mean between the radii that enclose 99% of the flux
	// in the pure De Vacouleur/exponential disk case
	// 6.670 = 3.975*Re = 3.975*1.678*Rdisk
  float BtoT = getBtoT();
	setRadius(6.670*current.Rdisk*(1 - BtoT) + 18.936 * current.Reff * BtoT);
	
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
	x[0] = (y[0]-source_x[0])*current.cosPA +
      (y[1]-source_x[1])*current.sinPA;
  x[1] = -(y[0]-source_x[0])*current.sinPA +
      (y[1]-source_x[1])*current.cosPA;
      
  double sb = spheroid.SurfaceBrightness(x);
      
  PosType R = current.cxx*x[0]*x[0] + current.cyy*x[1]*x[1];
	R = sqrt(R);

	//sb = sbDo*exp(-(R)) + sbSo*exp(-7.6693*pow(R/Reff,0.25));
  sb += current.sbDo*exp(-R);
      
	//if(current.Reff > 0.0) sb += current.sbSo*exp(-7.6693*pow((x[0]*x[0] + x[1]*x[1])
  //                                                          /current.Reff/current.Reff,0.125));

	if(sb< sb_limit)
		return 0.;
	
	return sb;
}

PosType SourceOverzier::getTotalFlux() const{
	return mag_to_flux(current.mag);
}

void SourceOverzier::printSource(){
  current.print();
}

void SourceOverzier::assignParams(InputParams& /* params */)
{
}

PosType SourceOverzier::getMag(Band band) const {
  return current.mag_map.at(band);
}

PosType SourceOverzier::getMagBulge(Band band) const {
  return current.bulge_mag_map.at(band);
}

 void SourceOverzier::changeBand(Band band){
   
   current.mag = current.mag_map[band];
   current.mag_bulge = current.bulge_mag_map[band];
  
  renormalize_current();
}

void SourceOverzier::renormalize_current(){
  float BtoT = getBtoT();
  if(current.Rdisk > 0.0){
    double det = 2*PI/sqrt(current.cxx*current.cyy);
    current.sbDo = mag_to_flux(current.mag)*(1-BtoT)/det;
  }
  else current.sbDo = 0.0;
}


SourceOverzierPlus::SourceOverzierPlus(PosType my_mag
                                       ,PosType my_mag_bulge
                                       ,PosType my_Reff
                                       ,PosType my_Rdisk
                                       ,PosType my_PA
                                       ,PosType inclination
                                       ,unsigned long my_id
                                       ,PosType my_z
                                       ,const PosType *theta
                                       ,PosType zeropoint
                                       ,Utilities::RandomNumbers_NR &ran
                                       ):
SourceOverzier(my_mag,my_mag_bulge,my_Reff,my_Rdisk,0,inclination,my_id,my_z,0,zeropoint)
,spheroid(my_mag_bulge,my_Reff, my_PA,4,1,my_z,zeropoint)
{
  assert(my_mag_bulge >= my_mag);
  //std::cout << "SourceOverzierPlus constructor" << std::endl;
  original = current;

  float minA = 0.01,maxA = 0.4; // minimum and maximum amplitude of arms
  Narms = (ran() > 0.2) ? 2 : 4;  // number of arms
  arm_alpha = (21 + 10*(ran()-0.5)*2)*degreesTOradians; // arm pitch angle
  mctalpha = Narms/tan(arm_alpha);
  disk_phase = PI*ran(); // add phase of arms
  Ad = minA + (maxA-minA)*ran();
  
  // extra cersic component
  //double index = 4 + 3*(ran()-0.5)*2;
  
  //double index = 4*pow(MAX(getBtoT(),0.03),0.4)*pow(10,0.2*(ran()-0.5));
  //double index = ran() + 3.5 ;
  double index = 3.5 + 0.5*ran();  // ????

  double q = 1 - 0.5*ran();
  spheroid.setSersicIndex(index);
  
  spheroid.ReSet(my_mag_bulge,my_Reff,0,index,q,my_z,0);
  spheroid.setTheta(0, 0);
  
  PA = my_PA;
  cosPA = cos(my_PA );
  sinPA = sin(-my_PA );
  cosi  = cos(current.inclination);
  
  modes.resize(6);
  for(PosType &mod : modes){
    mod = 2.0e-2*ran();
  }

  setTheta(theta);
  
  assert(original.Rdisk != 0 || original.Reff != 0);
}

SourceOverzierPlus::~SourceOverzierPlus(){
  //std::cout << "SourceOverzierPlus destructor" << std::endl;
  //spheroid;
}

SourceOverzierPlus::SourceOverzierPlus(const SourceOverzierPlus &p):
SourceOverzier(p),
Narms(p.Narms),Ad(p.Ad),mctalpha(p.mctalpha),arm_alpha(p.arm_alpha)
,disk_phase(p.disk_phase),PA(p.PA),cosPA(p.cosPA),sinPA(p.sinPA),cosi(p.cosi),spheroid(p.spheroid)
{
  //delete spheroid;
  /*spheroid = new SourceSersic(p.spheroid->getMag(),p.getReff()
                              ,p.spheroid->getPA()
                              ,p.spheroid->getSersicIndex()
                              ,p.spheroid->getAxesRatio()
                              ,p.spheroid->getZ()
                              ,p.spheroid->getTheta().x);
   */
  
  Narms=p.Narms;
  Ad=p.Ad;
  mctalpha=p.mctalpha;
  arm_alpha=p.arm_alpha;
  disk_phase=p.disk_phase;
  cosi=p.cosi;
  
  spheroid = p.spheroid;
  modes = p.modes;
  original = p.original;
}

SourceOverzierPlus & SourceOverzierPlus::operator=(const SourceOverzierPlus &p){
  if(this == &p) return *this;
  SourceOverzier::operator=(p);
  
  Narms=p.Narms;
  Ad=p.Ad;
  mctalpha=p.mctalpha;
  arm_alpha=p.arm_alpha;
  disk_phase=p.disk_phase;
  PA = p.PA;
  cosPA=p.cosPA;
  sinPA=p.sinPA;
  cosi=p.cosi;
  
  spheroid = p.spheroid;
  modes = p.modes;
  original = p.original;
  
  return *this;
}


PosType SourceOverzierPlus::SurfaceBrightness(PosType *y){
  // position relative to center
  Point_2d x;
  x[0] = (y[0]-source_x[0]) * cosPA -
         (y[1]-source_x[1]) * sinPA;
  x[1] = (y[0]-source_x[0]) * sinPA +
         (y[1]-source_x[1]) * cosPA;
  
  PosType xlength = x.length();

  PosType sb = 0;

  if(xlength > 0 && current.sbDo > 0){
    Point_2d z=x;
    z[1] /= cosi; // tranformed into face-on
    
    //PosType R = sqrt( cxx*x[0]*x[0] + cyy*x[1]*x[1] + cxy*x[0]*x[1] );
    PosType R = z.length()/current.Rdisk;
    PosType theta = atan2(z[1],z[0]);
    
    //disk
    sb = current.sbDo*exp(-R);
    //spiral arms
    PosType phir = mctalpha*log(R/0.5) + disk_phase;
    //PosType phiro = 0.4*mctalpha*log(R/0.5) + disk_phase;
    //if(R > 0.5 )
    sb *= 1 + Ad*cos(Narms*theta + phir);
    //else disk_sb *= 1 + Ad*cos(Narms*theta + phiro);
  }else{
    sb = current.sbDo;
  }
  
  assert(sb >= 0.0);
  
  //std::cout << "disk sb " << sb << std::endl;
  
  // bulge
  if(current.Reff > 0.0){
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
    sb += spheroid.SurfaceBrightness(x.x)*perturb;
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
  if(current.Reff > 0.0){
    spheroid.setMag(current.mag_bulge);
  }
}

void SourceOverzierPlus::randomize(Utilities::RandomNumbers_NR &ran){
  
  // reset to original parameters to prevent random walk
  current = original;
  
  float BtoT;
  { // SourceOverzier variables
    
    float minsize = 0.01*arcsecTOradians;
    current.Reff *= MAX( (1 + 0.2*(2*ran()-1.))
                        , minsize );
    current.Rdisk *= MAX( (1 + 0.2*(2*ran()-1.))
                      , minsize);
    
    PosType tmp = 0.01*(2*ran()-1.);
    
    for(auto mag = current.mag_map.begin() ; mag != current.mag_map.end() ; ++mag){
      mag->second = mag->second + tmp;
    }
    
    for(auto mag = current.bulge_mag_map.begin() ; mag != current.bulge_mag_map.end() ; ++mag){
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
    current.mag += tmp;
    
    current.inclination = ran()*0.90*PI/2;
    
    cosi  = cos(current.inclination);

    if(current.Rdisk > 0.0){
      current.cxx = 1.0/current.Rdisk/current.Rdisk;
      current.cyy = 1.0/cosi/cosi/current.Rdisk/current.Rdisk;
    }else{
      current.cxx = current.cyy = 0.0;
    }
    
    SourceOverzier::renormalize_current();
    
    // radius
    // weighted mean between the radii that enclose 99% of the flux
    // in the pure De Vacouleur/exponential disk case
    // 6.670 = 3.975*Re = 3.975*1.678*Rdisk
    BtoT = getBtoT();
    setRadius(6.670*current.Rdisk*(1-BtoT) + 18.936*current.Reff*BtoT);
  }
  
  float minA = 0.01,maxA = 0.2; // minimum and maximum amplitude of arms
  Narms = (ran() > 0.2) ? 2 : 4;  // number of arms
  arm_alpha = (21 + 10*(ran()-0.5)*2)*degreesTOradians; // arm pitch angle
  mctalpha = Narms/tan(arm_alpha);
  disk_phase = PI*ran(); // add phase of arms
  Ad = minA + (maxA-minA)*ran();
 
  // spheroid
  
  // extra cersic component
  //double index = 4 + 3*(ran()-0.5)*2;
  
  
  //double index = 4*pow(MAX(BtoT,0.03),0.4)*pow(10,0.2*(ran()-0.5));
  //double index = 4*pow(10,0.2*(ran()-0.5));
  //double index = ran() + 3.5;  //
  double index = 3.5 + 0.5*ran();  // ???
  assert(index > 0.5 && index < 9);
  double q = 1 + (0.5-1)*ran();
  
  /*delete spheroid;
  spheroid = new SourceSersic(mag_bulge,Reff/arcsecTOradians,PA + 10*(ran() - 0.5)*PI/180,index,q,zsource,getTheta().x);
  */
  
  spheroid.ReSet(current.mag_bulge,current.Reff/arcsecTOradians,0,index,q,zsource,0);
  spheroid.setTheta(0, 0);
  
  PA = PI*ran();
  cosPA = cos(PA);
  sinPA = sin(-PA);

  for(PosType &mod : modes){
    mod = 5.0e-2*ran();
  }
  
  //SourceOverzier::printSource();
}


