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

	Reff = my_Reff;
	Rh = my_Rh;

	cxx = pow(cos(PA),2) + pow(cos(incl)*sin(PA),2);
	cyy = pow(sin(PA),2) + pow(cos(incl)*cos(PA),2);
	cxy = 2*cos(PA)*sin(PA)*(1-pow(cos(incl),2));

	//muDo = mag-2.5*log10(1-BtoT)+5*log10(Rh)+1.9955;
	//muSo = mag-2.5*log10(BtoT)+5*log10(Reff)-4.9384;

	sbDo = pow(10,-mag/2.5)*0.159148*(1-BtoT)/pow(Rh,2);
	sbSo = pow(10,-mag/2.5)*94.484376*BtoT/pow(Reff,2);
}
/// Surface brightness normalized so that the total luminosity is 10^(-mag/2.5)
double OverGalaxy::SurfaceBrightness(double *x){

	double R = cxx*x[0]*x[0] + cyy*x[1]*x[1] + cxy*x[0]*x[1];
	R = sqrt(R);

	return sbDo*exp(-(R/Rh)) + sbSo*exp(-7.6693*pow(R/Reff,0.25));
}

void OverGalaxy::print(){
	cout << "bulge half light radius: " << Reff << " arcs   disk scale hight: " << Rh << " arcs" << endl;
}
