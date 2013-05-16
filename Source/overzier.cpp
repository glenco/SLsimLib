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

OverzierSource::OverzierSource()
: haloID(0), Reff(0), Rh(0), BtoT(0), PA(0), inclination(0),
  cxx(0), cyy(0), cxy(0), sbDo(0), sbSo(0), mag(0)
{
}

OverzierSource::OverzierSource(
		double my_mag              /// Total magnitude
		,double my_BtoT            /// Bulge to total ratio
		,double my_Reff         /// Bulge half light radius (arcs)
		,double my_Rh           /// disk scale hight (arcs)
		,double my_PA           /// Position angle (radians)
		,double my_inclination  /// inclination of disk (radians)
		,unsigned long my_id          /// optional angular position on the sky
		,double my_z            /// optional redshift
		,const double *my_theta          /// optional angular position on the sky
		){
	setInternals(my_mag,my_BtoT,my_Reff,my_Rh,my_PA,my_inclination,my_id,my_z,my_theta);
}

OverzierSource::~OverzierSource()
{
}

/// Sets internal variables.  If default constructor is used this must be called before the surface brightness function.
void OverzierSource::setInternals(double my_mag,double my_BtoT,double my_Reff,double my_Rh,double my_PA,double incl,unsigned long my_id,double my_z,const double *my_theta){

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
double OverzierSource::SurfaceBrightness(
		double *y  /// position in radians
		){
	// position relative to center
	double x[2];
	x[0] = y[0]-getX()[0];
	x[1] = y[1]-getX()[1];
	
	double R = cxx*x[0]*x[0] + cyy*x[1]*x[1] + cxy*x[0]*x[1],sb;
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

double OverzierSource::getTotalFlux(){
	return pow(10,-(48.6+mag)/2.5);
}

void OverzierSource::printSource(){
	std::cout << "bulge half light radius: " << Reff << " arcs   disk scale hight: " << Rh << " arcs" << std::endl;
}

void OverzierSource::assignParams(InputParams& /* params */)
{
}