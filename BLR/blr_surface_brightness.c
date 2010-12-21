/*
 * blr_surface_brightness.c
 *
 *  Created on: Sep 20, 2010
 *      Author: bmetcalf
 */
#include <math.h>

// computes the surface brightness associated with our simplest model
// for the BLR. x is the projected "distance from center" coordinate
// and tau the delay time (i.e. specifies an instant as the observer
// sees it. nu is the photon frequency we are interested in and nu_0
// defines the line center (rest) frequency of the line the
// normalization is arbitrary for now.

// inputs are in Mpc,days,Hz,Hz
double blr_surface_brightness(double x,double tau,double nu
		,double nu_o){

	// some constants which need to be set in the real version but will be
	// plucked from the air here

	static const float clight = 8.39428142e-10; // Mpc / day
	static const float r_in = 2.3884456e-9; // Mpc
	static const float r_out = 4.2992021e-7; // Mpc
	static const float gam = -1.0;
	float r,y;

	r = (x*x + (clight*clight*tau*tau))/(2*clight*tau);

	if ( (r < r_in) || (r > r_out) ) return 0.0;

	y = r - clight*tau;

	//if( y < 0.0 ) return 0.0;  // optically thick accretion disk

	return  x/pow(x*x + (clight*clight*tau*tau),1.5) * pow(r/r_in,gam-0.5);
}
/* unit filled version
double blr_surface_brightness(double x,double tau,double nu
		,double nu_o){

	// some constants which need to be set in the real version but will be
	// plucked from the air here

	static const float rho_o = 1.0;
	static const float r_in = 7.3699794e15; // cm
	static const float r_out = 1.3265963e18 ; //cm
	static const float eta_0 = 1.0;
	static const float E_0 = 1.0;
	static const float gam = -1.0;
	static const float sigma_0 = 6.0e8 ; // cm/s
	static const float clight = 3.0e10 ; // cm/s
	static const float pi = 3.141593 ;

	// convert from Mpc to cm and days to seconds
	x*=3.08568e24;
	tau*= 8.64e4;

	float r,y;

	r = (x*x + (clight*clight*tau*tau))/(2*clight*tau);
	y = r - clight*tau;

	if ((r < r_in) || (r > r_out) || (y < 0) ) return 0.0;

	float rho,sig_nu;

	rho = rho_o * pow(r/r_in,gam);
	sig_nu = sigma_0 * pow(r/r_in,-0.5) * nu_o/clight;

	return  1.0e48*x/pow(x*x + (clight*clight*tau*tau),1.5) * rho/sig_nu;

	return eta_0 * E_0 * 2 * clight * x
			/pow(x*x + (clight*clight*tau*tau),1.5)
			/sqrt(2*pi)*exp( -0.5* pow((nu - nu_o)/sig_nu,2)) * rho/sig_nu;
}
*/
