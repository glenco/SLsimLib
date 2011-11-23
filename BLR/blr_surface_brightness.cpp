/*
 * blr_surface_brightness.c
 *
 *  Created on: Sep 20, 2010
 *      Author: bmetcalf
 */
#include <slsimlib.h>

// computes the surface brightness associated with our simplest model
// for the BLR. x is the projected "distance from center" coordinate
// and tau the delay time (i.e. specifies an instant as the observer
// sees it. nu is the photon frequency we are interested in and nu_0
// defines the line center (rest) frequency of the line the
// normalization is arbitrary for now.

// inputs are in Mpc,days,Hz,Hz
double blr_surface_brightness_spherical(double x,AnaLens *lens){

	// some constants which need to be set in the real version but will be
	// plucked from the air here

	//static const float clight = 8.39428142e-10; // Mpc / day
	//static const float r_in = 2.3884456e-9; // Mpc
	//static const float r_out = 4.2992021e-7; // Mpc
	//static const float gam = -1.0;
	float r,y,tau;

	tau = lens->source_tau*8.39428142e-10;  // convert days to Mpc

	r = (x*x + tau*tau)/(2*tau);

	if ( (r < lens->source_r_in ) || (r > lens->source_r_out) ) return 0.0;

	y = r - tau;

	//if( y < 0.0 ) return 0.0;  // optically thick accretion disk

	return  x/pow(x*x + tau*tau,1.5) * pow(r/lens->source_r_in,lens->source_gamma-0.5);
}

double blr_surface_brightness_disk_old(double x[],AnaLens *lens,COSMOLOGY *cosmo){

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

	if(lens->zlens != oldzlens || lens->zsource != oldzsource){
		DlDs = cosmo->angDist(0,lens->zlens)
        		/cosmo->angDist(0,lens->zsource);
		oldzlens = lens->zlens;
		oldzsource = lens->zsource;
	}

	//printf("hi from BLR disk\n");

	if(lens->source_nuo <= 0.0){ERROR_MESSAGE(); exit(1); }

	tau = lens->source_tau*8.39428142e-10/(1+lens->zsource);  // convert days to Mpc and account for time dilation

	// first get the "r" and "z" coordinates

	r = sqrt( x[0]*x[0] + x[1]*x[1] )/DlDs;                    // 2d radius
	R = (r*r + tau*tau)/(2*tau);   // 3d radius

	if ((R < lens->source_r_in ) || (R > lens->source_r_out)) return 0.0;

	zz = R - tau;

	//zz = x[1]/DlDs *tan(lens->source_inclination);  // allows whole thing to be mapped

	//now rotate to coordinates aligned with the disk - this is a rotation
	//                                                  by angle inc about
	//                                                  the x-direction

	yy_prime = x[1]/DlDs * cos(lens->source_inclination) + zz * sin(lens->source_inclination);
	zz_prime =  zz * cos(lens->source_inclination) - x[1]/DlDs * sin(lens->source_inclination);
	xx_prime = x[0]/DlDs;

	if(fabs(zz_prime) > R*sin(lens->source_opening_angle) ) return 0.0;  // outside of disk

	rho = pow(R/lens->source_r_in ,lens->source_gamma);

	if(lens->source_monocrome) return pow(lens->source_r_in,2) * rho / (r*r + tau*tau);

	float vr = sqrt(4.7788e-20*lens->source_BHmass/R);
	float v_shift_yprime;         // in units of the speed of light

	if(lens->source_opening_angle < 0.9999*pi/2){

		//double v_shift_xprime,v_shift_yprime,v_shift_zprime;
		//v_shift_xprime = - vr * yy_prime/sqrt(xx_prime*xx_prime + yy_prime*yy_prime);
		v_shift_yprime = vr * xx_prime/sqrt(xx_prime*xx_prime + yy_prime*yy_prime);
		//v_shift_zprime = 0.0;

		// the projected shift is only the z-component (i.e. the Doppler motion
		// towards the observer and de-rotate it back to the normal frame
		//double v_shift_x,v_shift_y,v_shift_z,nu_shift;
		//v_shift_x = v_shift_xprime;
		//v_shift_y = v_shift_yprime * cos(lens->source_inclination) - v_shift_zprime * sin(lens->source_inclination);
		//v_shift_z = v_shift_zprime * cos(lens->source_inclination) + v_shift_yprime * sin(lens->source_inclination);

		double v_shift_z,nu_shift;

		v_shift_z = v_shift_yprime * sin(lens->source_inclination);

		//sig_nu = sigma_0 * lens->source_nuo;
		sig_nu = MAX(sigma_0, 0.01 * vr) * lens->source_nuo;
		//sig_nu = 0.05 * vr * lens->source_nuo;

		nu_shift = v_shift_z * lens->source_nuo;

		return pow(lens->source_r_in,2)*exp( -0.5*pow( (lens->source_nu*(1+lens->zsource) - lens->source_nuo - nu_shift)/sig_nu ,2) )
					* rho / (r*r + tau*tau)/sig_nu;
	}

	// isotropic circular Keplerian orbits

	v_shift_yprime = ( lens->source_nuo/lens->source_nu/(1+lens->zsource) - 1 ) * R/r;  // convert to tangent of sphere

	if( fabs(v_shift_yprime) > vr ) return 0.0;

	return  pow(lens->source_r_in,2) * rho/(r*r + tau*tau) /sqrt(vr*vr - v_shift_yprime*v_shift_yprime);
}
