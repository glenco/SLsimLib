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
double blr_surface_brightness_spherical(double x,double tau,double nu
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

double blr_surface_brightness_disk(double x,double phi,double tau,double nu
		,double nu_o){

	/* computes the surface brightness associated with our simplest model
	// for the BLR. x is the projected "distance from center" coordinate
	// and tau the delay time (i.e. specifies an instant as the observer
	// sees it. nu is the photon frequency we are interested in and nu_0
	// defines the line center (rest) frequency of the line the
	// normalization is arbitrary for now.

	// some constants which need to be set in the real version but will be
	// plucked from the air here
	 *
	 *   Translated from IDL routine written by S. Sim
	*/
//***
	static const float rho_0 = 1.0;
	static const float r_in = 8.8439753e+15; // cm
	static const float r_out = 4.0e18 ; // cm
	static const float eta_0 = 1.0;
	static const float E_0 = 1.0;
	static const float gam = 0.0;
	static const float sigma_0 = 6.0e7; // cm/s
	static const float inc=30. * 3.14159/180 ; // in radians
	static const float opening_angle = 5. * 3.14159/180 ; //in radians
	static float MBH = 1.0e9 ; // black hole mass, in Msun

	static const float CLIGHT = 3.0e10 ; // speed of light in cm/s
	static const float GNEWTON = 6.67259e-8 ; // grav. constant
	static const float MSUN = 1.989e33 ; // Solar mass

	double r,zz,theta,xx,yy,try,rho;
	double rr_prime,zz_prime,theta_prime,xx_prime,yy_prime,phi_prime;
	double sig_nu;

	MBH = MBH *MSUN;


	// first get the "r" and "z" coordinates

	r = (x*x + (CLIGHT*CLIGHT*tau*tau))/(2*CLIGHT*tau);
	zz = r - CLIGHT*tau;

	if ((r < r_in) || (r > r_out)) return 0.0;

	theta = acos(zz/r);

	xx = r * sin(theta) * cos(phi);
	yy = r * sin(theta) * sin(phi);

	//now rotate to coordinates aligned with the disk - this is a rotation
	//                                                  by angle inc about
	//                                                  the x-direction

	yy_prime = yy * cos(inc) + zz * sin(inc);
	zz_prime = zz * cos(inc) - yy * sin(inc);
	xx_prime = xx;

	//now put back to polar coordinates in the prime frame

	rr_prime = r; //radius is preserved
	theta_prime = acos(zz_prime/rr_prime);
	if (zz_prime > 0.999999*rr_prime) theta_prime = 0.0;
	if (zz_prime < (-0.9999999*rr_prime)) theta_prime = 3.14159;

	//print, theta_prime
	try = xx_prime/rr_prime/sin(theta_prime);
	if (try < -0.9999999999999){
		phi_prime = 3.14159;
	} else if (try > 0.999999999999){
		phi_prime = 0.0;
	} else {
		phi_prime = acos(try);
	}


	//phi_prime = acos(xx_prime/rr_prime/sin(theta_prime))
	if (yy_prime < 0.0) phi_prime = (2*3.14159) - phi_prime;

	rho = 0.0;
	if ((((3.14159/2.) - theta_prime) > 0.0) && (((3.14159/2.) - theta_prime) < opening_angle))
		rho = rho_0 * pow(r/r_in,gam);

	sig_nu = sigma_0 * nu_o/CLIGHT;

	double v_shift_xprime,v_shift_yprime,v_shift_zprime;

	//now want to compute the Keplerian velocity vector in the primed frame
	v_shift_xprime = -1. * sqrt(GNEWTON * MBH / r) * sin(phi_prime) * sin(theta_prime);
	v_shift_yprime = sqrt(GNEWTON * MBH / r) * cos(phi_prime) * sin(theta_prime);
	v_shift_zprime = 0.0;
	//and derotate it back to the normal frame

	double v_shift_x,v_shift_y,v_shift_z,nu_shift;

	v_shift_x = v_shift_xprime;
	v_shift_y = v_shift_yprime * cos(-1.*inc) + v_shift_zprime * sin(-1.*inc);
	v_shift_z = v_shift_zprime * cos(-1.*inc) - v_shift_yprime * sin(-1.*inc);

	//the projected shift is only the z-component (i.e. the Doppler motion
	//towards the observer
	nu_shift = v_shift_z * nu_o / CLIGHT;

	return eta_0 * E_0 * 2 * CLIGHT / (x*x + (CLIGHT*CLIGHT*tau*tau)) /sqrt(2.*3.14159)*exp( -1.* pow(nu - nu_o - nu_shift,2)
			/ 2 / sig_nu/sig_nu) * rho/sig_nu;
}
