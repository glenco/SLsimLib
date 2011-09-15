/*
 * blr_surface_brightness2.0.c
 *
 *  Created on: Sep 20, 2010
 *      Author: bmetcalf
 *      revised by S. Sim Sep 1, 2011
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <analytic_lens.h>
#include <source_model.h>
#include <source_model.h>
#include <Tree.h>

// computes the surface brightness associated with our simplest model
// for the BLR. x is the projected "distance from center" coordinate
// and tau the delay time (i.e. specifies an instant as the observer
// sees it. lens->source_nu is the photon frequency we are interested in and lens->source_nuo
// defines the line center (rest) frequency of the line the
// normalization is arbitrary for now.

// inputs are in Mpc,days,Hz,Hz

double blr_surface_brightness_spherical_random_motions(double x,AnaLens *lens,COSMOLOGY *cosmo){

	// some constants which need to be set in the real version but will be
	// plucked from the air here

	float r,tau, sigma2, eta;
	
	static float DlDs;  //
	static double oldzlens=0,oldzsource=0;
	if(lens->zlens != oldzlens || lens->zsource != oldzsource){
		DlDs = angDist(0,lens->zlens,cosmo)
      		/angDist(0,lens->zsource,cosmo);
		oldzlens = lens->zlens;
		oldzsource = lens->zsource;
	}

	x /= DlDs;


	tau = lens->source_tau*8.39428142e-10/(1+lens->zsource);  // convert days to Mpc

	r = (x*x + tau*tau)/(2*tau);

	if ( (r < lens->source_r_in ) || (r > lens->source_r_out) ) return 0.0;

	if(lens->source_monocrome) return pow(r/lens->source_r_in,lens->source_gamma)*r/tau;

	//calculated the eta function
	//sigma2 = (lens->source_nuo*lens->source_nuo*1.1126501e-21) * ((8.25111e+07 * temp/mass) + (vturb*vturb) + (f_K*f_K*42.9497*lens->source_BHmass/r));
	  //                                         1/c^2              kB/Mh

	sigma2 = pow(lens->source_nuo*lens->source_sigma,2)*1.1126501e-11;
	//                                                  1/c^2 in km/s

	eta = lens->source_nuo *(1+lens->zsource)* exp(-0.5 * pow( lens->source_nu*(1+lens->zsource) - lens->source_nuo ,2) / sigma2) / sqrt(sigma2);
	  //  1/sqrt(2pi)

	return  pow(r/lens->source_r_in,lens->source_gamma)*eta*r/tau;
}

double blr_surface_brightness_spherical_circular_motions(double x,AnaLens *lens,COSMOLOGY *cosmo){

	float r,tau, eta, sin_theta, nu_m;

	static float DlDs;  //
	static double oldzlens=0,oldzsource=0;
	if(lens->zlens != oldzlens || lens->zsource != oldzsource){
		DlDs = angDist(0,lens->zlens,cosmo)
        		/angDist(0,lens->zsource,cosmo);
		oldzlens = lens->zlens;
		oldzsource = lens->zsource;
	}

	x /= DlDs;

	tau = lens->source_tau*8.39428142e-10/(1+lens->zsource);  // convert days to Mpc

	r = (x*x + tau*tau)/(2*tau);

	if ( (r < lens->source_r_in ) || (r > lens->source_r_out) ) return 0.0;

	if(lens->source_monocrome) return pow(r/lens->source_r_in,lens->source_gamma)*r/tau;

	sin_theta = x/r;

	//calculated the eta function

	// maximum frequency at theta
	nu_m = lens->source_nuo * sqrt( 4.7788e-20 * lens->source_BHmass /r ) * sin_theta;

	if ( fabs(lens->source_nu*(1+lens->zsource) - lens->source_nuo) < nu_m ) return 0.0;

	eta = lens->source_nuo *(1+lens->zsource)/sqrt( 1. - pow( (lens->source_nu*(1+lens->zsource) - lens->source_nuo)/nu_m ,2) )/nu_m/pi;

	return  pow(r/lens->source_r_in,lens->source_gamma)*eta*r/tau;
}

double blr_surface_brightness_disk(double x[],AnaLens *lens,COSMOLOGY *cosmo){

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
        //static float sigma_0 = 2.0e-3; // in units of c

  //BEN - I'VE REMOVED THE FIXED FLOOR SIGMA_0 SINCE THE THERMAL/TURBULENT CONTRIBUTION IS MORE ELEGANT AND SOLVES THIS PROBLEM

	//static const float inc=30. * 3.14159/180 ; // in radians
	//static const float opening_angle = 5. * 3.14159/180 ; //in radians
	//static float MBH = 1.0e9 ; // black hole mass, in Msun
	//static const float CLIGHT = 3.0e10 ; // speed of light in cm/s
	//static const float GNEWTON = 6.67259e-8 ; // grav. constant
	//static const float MSUN = 1.989e33 ; // Solar mass

	//EXTRA THINGS THIS ROUTINE NEEDS TO KNOW: 
	// mass = MASS OF EMITTING ION, IN UNITS OF HYDROGEN MASS
        // temp = TEMPERATURE FOR THERMAL BROADENING IN KELVIN
        // turb = TURBULENT BROADENING INSIEDE CLOUD IN CM/S
	//I AM ASSUMING THAT source_nuo is in Hz

	double R,r,Z,tau;
	double zz_prime,xx_prime;//,yy_prime;
	double sigma2,eta;

	static float DlDs;  //
	static double oldzlens=0,oldzsource=0;
	if(lens->zlens != oldzlens || lens->zsource != oldzsource){
		DlDs = angDist(0,lens->zlens,cosmo)
        		/angDist(0,lens->zsource,cosmo);
		oldzlens = lens->zlens;
		oldzsource = lens->zsource;
	}

	//printf("hi from BLR disk\n");

	if(lens->source_nuo <= 0.0){ERROR_MESSAGE(); exit(1); }

	tau = lens->source_tau*8.39428142e-10/(1+lens->zsource);  // convert days to Mpc and account for time dilation

	// first get the "r", "R" and "Z" coordinates

	R = sqrt( x[0]*x[0] + x[1]*x[1] )/DlDs;                    // 2d radius
	r = (R*R + tau*tau)/(2*tau);   // 3d radius

	if ((r < lens->source_r_in ) || (r > lens->source_r_out)) return 0.0;

	Z = r - tau;

	//zz = x[1]/DlDs *tan(lens->source_inclination);  // allows whole thing to be mapped

	//now rotate to coordinates aligned with the disk - this is a rotation
	//                                                  by angle inc about
	//                                                  the x-direction

	//yy_prime = x[1]/DlDs * cos(lens->source_inclination) + Z * sin(lens->source_inclination);
	zz_prime = Z * cos(lens->source_inclination) - x[1]/DlDs * sin(lens->source_inclination);
	xx_prime = x[0]/DlDs;

	if(fabs(zz_prime) > r*sin(lens->source_opening_angle) ) return 0.0;  // outside of disk

	// this is an option to do a monochromatic version
	if(lens->source_monocrome) return pow(lens->source_r_in,lens->source_gamma) * r / tau;

  //now to compute eta

	float v_Kep = sqrt(4.7788e-20*lens->source_BHmass/r);
	float v_shift_yprime;         // in units of the speed of light

	//	if(lens->source_opening_angle < 0.9999*pi/2){
	// BEN -- I'VE REMOVED THIS CONDITION SINCE I THINK THE LAST FEW LINES (THE ALTERNATIVE) MIGHT BE REDUNDANT NOW

	v_shift_yprime = v_Kep * xx_prime/r;
	// xx_prime/r comes about from xx_prime/sqrt(xx_prime^2+yy_prime^2) * sqrt(xx_prime^2+yy_prime^2)/r
	//                              this is cos (phi_prime)                this is sin (theta_prime)
	
	// the projected shift is only the z-component (i.e. the Doppler motion
	// towards the observer and de-rotate it back to the normal frame
	
	double v_shift_z,nu_shift;
	
	v_shift_z = v_shift_yprime * sin(lens->source_inclination);
	
	nu_shift = v_shift_z * lens->source_nuo;
	//BEN - CONFIRM UNITS ARE CORRECT HERE?
	//sigma2 = (lens->source_nuo*lens->source_nuo*1.1126501e-21) * ((8.25111e+07 * temp/mass) + (vturb*vturb));
	//                                         1/c^2              kB/Mh
	
	sigma2 = pow(lens->source_nuo*lens->source_sigma,2)*1.1126501e-11;
	//                                                  1/c^2 in km/s

	eta = lens->source_nuo *(1+lens->zsource) * exp(-0.5 * pow( lens->source_nu*(1+lens->zsource) - lens->source_nuo - nu_shift ,2) / sigma2) / sqrt(sigma2);
	 //  1/sqrt(2pi)

	return pow(r/lens->source_r_in,lens->source_gamma)*eta*r/tau;
	//}

	// isotropic circular Keplerian orbits
	// BEN - AM I CORRECT THAT THIS IS NOW REDUNDANT?
	
		/*
	v_shift_yprime = ( lens->source_nuo/lens->source_nu/(1+lens->zsource) - 1 ) * R/r;  // convert to tangent of sphere

	if( fabs(v_shift_yprime) > vr ) return 0.0;

	return  pow(lens->source_r_in,2) * rho/(r*r + tau*tau) /sqrt(vr*vr - v_shift_yprime*v_shift_yprime);
		*/
}
