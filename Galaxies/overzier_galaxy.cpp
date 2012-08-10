/*
 * overier_galaxy.cpp
 *
 *  Created on: Apr 12, 2012
 *      Author: bmetcalf
 *
 *      This surface brightness model comes from Roderik Overzier in private communication although it
 *      is probably written into a paper somewhere.
 */
#include <slsimlib.h>

OverGalaxy::OverGalaxy(
		double mag              /// Total magnitude
		,double my_BtoT            /// Bulge to total ratio
		,double my_Reff         /// Bulge half light radius (arcs)
		,double my_Rh           /// disk scale hight (arcs)
		,double my_PA           /// Position angle (radians)
		,double my_inclination  /// inclination of disk (radians)
		){
	setInternals(mag,my_BtoT,my_Reff,my_Rh,my_PA,my_inclination);
}
/// Sets internal variables.  If default constructor is used this must be called before the surface brightness function.
void OverGalaxy::setInternals(double mag,double BtoT,double my_Reff,double my_Rh,double PA,double incl){

	Reff = my_Reff*pi/180/60/60;
	Rh = my_Rh*pi/180/60/60;

	cxx = ( pow(cos(PA),2) + pow(sin(PA)/cos(incl),2) )/Rh/Rh;
	cyy = ( pow(sin(PA),2) + pow(cos(PA)/cos(incl),2) )/Rh/Rh;
	cxy = ( 2*cos(PA)*sin(PA)*(1-pow(1/cos(incl),2)) )/Rh/Rh;

	//muDo = mag-2.5*log10(1-BtoT)+5*log10(Rh)+1.9955;
	//muSo = mag-2.5*log10(BtoT)+5*log10(Reff)-4.9384;

	sbDo = pow(10,-mag/2.5)*0.159148*(1-BtoT)/pow(Rh,2);
	sbSo = pow(10,-mag/2.5)*94.484376*BtoT/pow(Reff,2);
}
/// Surface brightness normalized so that the total luminosity is 10^(-mag/2.5). x should be in units of radians
double OverGalaxy::SurfaceBrightness(
		double *x  /// position in radians relative to center of source
		){

	double R = cxx*x[0]*x[0] + cyy*x[1]*x[1] + cxy*x[0]*x[1],sb;
	R = sqrt(R);

	//sb = sbDo*exp(-(R)) + sbSo*exp(-7.6693*pow(R/Reff,0.25));
	sb = sbDo*exp(-R) + sbSo*exp(-7.6693*pow((x[0]*x[0] + x[1]*x[1])/Reff/Reff,0.125));
	if(sb < 1.0e-3*(sbDo + sbSo) ) return 0.0;
	return sb;
}

void OverGalaxy::print(){
	std::cout << "bulge half light radius: " << Reff << " arcs   disk scale hight: " << Rh << " arcs" << std::endl;
}
