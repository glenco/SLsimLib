/*
 * blr_surface_brightness.c
 *
 *  Created on: Sep 20, 2010
 *      Author: bmetcalf
 */
#include "slsimlib.h"

// computes the surface brightness associated with our simplest model
// for the BLR. x is the projected "distance from center" coordinate
// and tau the delay time (i.e. specifies an instant as the observer
// sees it. nu is the photon frequency we are interested in and nu_0
// defines the line center (rest) frequency of the line the
// normalization is arbitrary for now.

// inputs are in Mpc,days,Hz,Hz
double blr_surface_brightness_spherical(double x,SourceBLR *source){

	// some constants which need to be set in the real version but will be
	// plucked from the air here

	//static const float clight = 8.39428142e-10; // Mpc / day
	//static const float r_in = 2.3884456e-9; // Mpc
	//static const float r_out = 4.2992021e-7; // Mpc
	//static const float gam = -1.0;
	float r,y,tau;

	tau = source->source_tau*8.39428142e-10;  // convert days to Mpc

	r = (x*x + tau*tau)/(2*tau);

	if ( (r < source->source_r_in ) || (r > source->source_r_out) ) return 0.0;

	y = r - tau;

	//if( y < 0.0 ) return 0.0;  // optically thick accretion disk

	return  x/pow(x*x + tau*tau,1.5) * pow(r/source->source_r_in,source->source_gamma-0.5);
}

double blr_surface_brightness_disk_old(double x[],SourceBLR *source){

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

	//static const float r_in = 8.8439753e+15; // cm
	//static const float r_out = 4.0e18 ; // cm
	//static const float gam = 0.0;
	static float sigma_0 = 2.0e-3; // in units of c
	//static const float inc=30. * 3.14159/180 ; // in radians
	//static const float opening_angle = 5. * 3.14159/180 ; //in radians
	//static float MBH = 1.0e9 ; // black hole mass, in Msun
	//static const float CLIGHT = 3.0e10 ; // speed of light in cm/s
	//static const float GNEWTON = 6.67259e-8 ; // grav. constant
	//static const float MSUN = 1.989e33 ; // Solar mass

	static float DlDs;  //
	static double oldzlens=0,oldzsource=0;
	double R,r,zz,rho,tau;
	double zz_prime,xx_prime,yy_prime;
	double sig_nu;

	DlDs = source->getDlDs();

	//printf("hi from BLR disk\n");

	if(source->source_nuo <= 0.0){ERROR_MESSAGE(); exit(1); }

	tau = source->source_tau*8.39428142e-10/(1+source->getZ());  // convert days to Mpc and account for time dilation

	// first get the "r" and "z" coordinates

	r = sqrt( x[0]*x[0] + x[1]*x[1] )/DlDs;                    // 2d radius
	R = (r*r + tau*tau)/(2*tau);   // 3d radius

	if ((R < source->source_r_in ) || (R > source->source_r_out)) return 0.0;

	zz = R - tau;

	//zz = x[1]/DlDs *tan(source->source_inclination);  // allows whole thing to be mapped

	//now rotate to coordinates aligned with the disk - this is a rotation
	//                                                  by angle inc about
	//                                                  the x-direction

	yy_prime = x[1]/DlDs * cos(source->source_inclination) + zz * sin(source->source_inclination);
	zz_prime =  zz * cos(source->source_inclination) - x[1]/DlDs * sin(source->source_inclination);
	xx_prime = x[0]/DlDs;

	if(fabs(zz_prime) > R*sin(source->source_opening_angle) ) return 0.0;  // outside of disk

	rho = pow(R/source->source_r_in ,source->source_gamma);

	if(source->source_monocrome) return pow(source->source_r_in,2) * rho / (r*r + tau*tau);

	float vr = sqrt(4.7788e-20*source->source_BHmass/R);
	float v_shift_yprime;         // in units of the speed of light

	if(source->source_opening_angle < 0.9999*pi/2){

		//double v_shift_xprime,v_shift_yprime,v_shift_zprime;
		//v_shift_xprime = - vr * yy_prime/sqrt(xx_prime*xx_prime + yy_prime*yy_prime);
		v_shift_yprime = vr * xx_prime/sqrt(xx_prime*xx_prime + yy_prime*yy_prime);
		//v_shift_zprime = 0.0;

		// the projected shift is only the z-component (i.e. the Doppler motion
		// towards the observer and de-rotate it back to the normal frame
		//double v_shift_x,v_shift_y,v_shift_z,nu_shift;
		//v_shift_x = v_shift_xprime;
		//v_shift_y = v_shift_yprime * cos(source->source_inclination) - v_shift_zprime * sin(source->source_inclination);
		//v_shift_z = v_shift_zprime * cos(source->source_inclination) + v_shift_yprime * sin(source->source_inclination);

		double v_shift_z,nu_shift;

		v_shift_z = v_shift_yprime * sin(source->source_inclination);

		//sig_nu = sigma_0 * source->source_nuo;
		sig_nu = MAX(sigma_0, 0.01 * vr) * source->source_nuo;
		//sig_nu = 0.05 * vr * source->source_nuo;

		nu_shift = v_shift_z * source->source_nuo;

		return pow(source->source_r_in,2)*exp( -0.5*pow( (source->source_nu*(1+source->getZ()) - source->source_nuo - nu_shift)/sig_nu ,2) )
					* rho / (r*r + tau*tau)/sig_nu;
	}

	// isotropic circular Keplerian orbits

	v_shift_yprime = ( source->source_nuo/source->source_nu/(1+source->getZ()) - 1 ) * R/r;  // convert to tangent of sphere

	if( fabs(v_shift_yprime) > vr ) return 0.0;

	return  pow(source->source_r_in,2) * rho/(r*r + tau*tau) /sqrt(vr*vr - v_shift_yprime*v_shift_yprime);
}
