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

OverGalaxy::OverGalaxy(double mag,double BtoT,double my_Reff,double my_Rh,double my_PA,double my_inclination){
	setInternals(mag,BtoT,my_Reff,my_Rh,my_PA,my_inclination);
}
/// Sets internal variables.  If default constructor is used this must be called before the surface brightness function.
void OverGalaxy::setInternals(double mag,double BtoT,double my_Reff,double my_Rh,double my_PA,double my_inclination){

	Reff = my_Reff;
	Rh = my_Rh;
	PA = my_PA;
	incl = my_inclination;

	muDo = mag-2.5*log10(1-BtoT)+5*log10(Rh)+1.9955;
	muSo = mag-2.5*log10(BtoT)+5*log10(Reff)-4.93884;
}
/// Surface brightness in magnitudes per arcsec^2
double OverGalaxy::SurfaceBrightness(double *x){

	double R = pow(cos(PA)*x[0] + sin(PA)*x[1],2)
			+ pow((-sin(PA)*x[0] + cos(PA)*x[1])/incl,2);
	R = sqrt(R);

	return muDo + muSo + 8.3268*pow(R/Reff,0.25) + 1.0857*(R/Rh);
}

