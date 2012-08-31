/*
 * blr_surface_brightness2.0.c
 *
 *  Created on: Sep 20, 2010
 *      Author: bmetcalf
 *      revised by S. Sim Sep 1, 2011
 */
/*#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <analytic_lens.h>
#include <Tree.h> */

#include <slsimlib.h>

// computes the surface brightness associated with our simplest model
// for the BLR. x is the projected "distance from center" coordinate
// and tau the delay time (i.e. specifies an instant as the observer
// sees it. source->source_nu is the photon frequency we are interested in and source->source_nuo
// defines the line center (rest) frequency of the line the
// normalization is arbitrary for now.

// inputs are in Mpc,days,Hz,Hz

double blr_surface_brightness_spherical_random_motions(double x,SourceBLR *source){

	// some constants which need to be set in the real version but will be
	// plucked from the air here

	float r,tau, sigma2, eta;
	
	static float DlDs;  //
	static double oldzlens=0,oldzsource=0;
	DlDs = source->getDlDs();

	x /= DlDs;


	tau = source->source_tau*8.39428142e-10/(1+source->getZ());  // convert days to Mpc

	r = (x*x + tau*tau)/(2*tau);

	if ( (r < source->source_r_in ) || (r > source->source_r_out) ) return 0.0;

	if(source->source_monocrome) return pow(r/source->source_r_in,source->source_gamma)*r/tau;

	//calculated the eta function
	//sigma2 = (source->source_nuo*source->source_nuo*1.1126501e-21) * ((8.25111e+07 * temp/mass) + (vturb*vturb) + (f_K*f_K*42.9497*source->source_BHmass/r));
	  //                                         1/c^2              kB/Mh

	//sigma2 = pow(source->source_nuo*source->source_sigma,2)*1.1126501e-11;
	//                                                  1/c^2 in km/s
	float v_Kep = sqrt(4.7788e-20*source->source_BHmass/r);

	sigma2 = pow(source->source_nuo*v_Kep*source->source_fK,2);

	eta = source->source_nuo *(1+source->getZ())* exp(-0.5 * pow( source->source_nu*(1+source->getZ()) - source->source_nuo ,2) / sigma2) / sqrt(sigma2);
	  //  1/sqrt(2pi)

	return  pow(r/source->source_r_in,source->source_gamma)*eta*r/tau;
}

double blr_surface_brightness_spherical_circular_motions(double x,SourceBLR *source){

	float r,tau, eta, sin_theta, nu_m;

	static float DlDs;  //
	static double oldzlens=0,oldzsource=0;
	DlDs=source->getDlDs();

	x /= DlDs;

	tau = source->source_tau*8.39428142e-10/(1+source->getZ());  // convert days to Mpc

	r = (x*x + tau*tau)/(2*tau);

	if ( (r < source->source_r_in ) || (r > source->source_r_out) ) return 0.0;

	if(source->source_monocrome) return pow(r/source->source_r_in,source->source_gamma)*r/tau;

	sin_theta = x/r;

	//calculated the eta function

	// maximum frequency at theta
	nu_m = source->source_nuo * sqrt( 4.7788e-20 * source->source_BHmass /r ) * sin_theta;

	if ( fabs(source->source_nu*(1+source->getZ()) - source->source_nuo) > nu_m ) return 0.0;

	eta = source->source_nuo *(1+source->getZ())/sqrt( 1. - pow( (source->source_nu*(1+source->getZ()) - source->source_nuo)/nu_m ,2) )/nu_m/pi;

	//printf("tau = %e eta =%e r=%e xi=%e\n",tau,eta,r,pow(r/source->source_r_in,source->source_gamma));
	return  pow(r/source->source_r_in,source->source_gamma)*eta*r/tau;
}

double blr_surface_brightness_disk(double x[],SourceBLR *source){

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

	double R,r,Z,tau;
	double zz_prime,xx_prime;//,yy_prime;
	double sigma2,eta;

	static float DlDs;  //
	static double oldzlens=0,oldzsource=0;
	DlDs = source->getDlDs();

	//printf("hi from BLR disk\n");

	if(source->source_nuo <= 0.0){ERROR_MESSAGE(); exit(1); }

	tau = source->source_tau*8.39428142e-10/(1+source->getZ());  // convert days to Mpc and account for time dilation

	// first get the "r", "R" and "Z" coordinates

	R = sqrt( x[0]*x[0] + x[1]*x[1] )/DlDs;                    // 2d radius
	r = (R*R + tau*tau)/(2*tau);   // 3d radius

	if ((r < source->source_r_in ) || (r > source->source_r_out)) return 0.0;

	Z = r - tau;

	//zz = x[1]/DlDs *tan(source->source_inclination);  // allows whole thing to be mapped

	//now rotate to coordinates aligned with the disk - this is a rotation
	//                                                  by angle inc about
	//                                                  the x-direction

	//yy_prime = x[1]/DlDs * cos(source->source_inclination) + Z * sin(source->source_inclination);
	zz_prime = Z * cos(source->source_inclination) - x[1]/DlDs * sin(source->source_inclination);
	xx_prime = x[0]/DlDs;

	if(fabs(zz_prime) > r*sin(source->source_opening_angle) ) return 0.0;  // outside of disk

	assert(pow(r/source->source_r_in,source->source_gamma) * r / tau >= 0.0);
	//printf("hello 2 %e \n",pow(r/source->source_r_in,source->source_gamma) * r / tau);
	// this is an option to do a monochromatic version
	if(source->source_monocrome) return pow(r/source->source_r_in,source->source_gamma) * r / tau;

  //now to compute eta

	float v_Kep = sqrt(4.7788e-20*source->source_BHmass/r);
	float v_shift_yprime;         // in units of the speed of light

	//	if(source->source_opening_angle < 0.9999*pi/2){
	// BEN -- I'VE REMOVED THIS CONDITION SINCE I THINK THE LAST FEW LINES (THE ALTERNATIVE) MIGHT BE REDUNDANT NOW

	v_shift_yprime = v_Kep * xx_prime/r;
	// xx_prime/r comes about from xx_prime/sqrt(xx_prime^2+yy_prime^2) * sqrt(xx_prime^2+yy_prime^2)/r
	//                              this is cos (phi_prime)                this is sin (theta_prime)
	
	// the projected shift is only the z-component (i.e. the Doppler motion
	// towards the observer and de-rotate it back to the normal frame
	
	double v_shift_z,nu_shift;
	
	v_shift_z = v_shift_yprime * sin(source->source_inclination);
	
	nu_shift = v_shift_z * source->source_nuo;
	//BEN - CONFIRM UNITS ARE CORRECT HERE?
	//sigma2 = (source->source_nuo*source->source_nuo*1.1126501e-21) * ((8.25111e+07 * temp/mass) + (vturb*vturb));
	//                                         1/c^2              kB/Mh
	
	//sigma2 = pow(source->source_nuo*source->source_sigma,2)*1.1126501e-11;
	sigma2 = pow(source->source_nuo*v_Kep*source->source_fK,2);
	//                                                  1/c^2 in km/s
	eta = pow( source->source_nu*(1+source->getZ()) - source->source_nuo - nu_shift ,2) / sigma2;
	if(eta > 9) return 0.0;  // prevents frequencies very far off the line from causing unnecessary refinements
	eta = pow(source->source_nuo*(1+source->getZ()),2) * exp(-0.5 * eta) / sqrt(sigma2);


	return pow(r/source->source_r_in,source->source_gamma)*eta*r/tau;
	//}

	// isotropic circular Keplerian orbits
	// BEN - AM I CORRECT THAT THIS IS NOW REDUNDANT?
	
		/*
	v_shift_yprime = ( source->source_nuo/source->source_nu/(1+source->getZ()) - 1 ) * R/r;  // convert to tangent of sphere

	if( fabs(v_shift_yprime) > vr ) return 0.0;

	return  pow(source->source_r_in,2) * rho/(r*r + tau*tau) /sqrt(vr*vr - v_shift_yprime*v_shift_yprime);
		*/
}
