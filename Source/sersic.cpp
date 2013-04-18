/*
 * Source/sersic.cpp
 *
 *  Created on: Feb 28, 2010
 *      Author: R.B. Metcalf
 */
#include "slsimlib.h"

SersicSource::SersicSource(
		double my_mag              /// Total magnitude
		,double my_Reff         /// Bulge half light radius (arcs)
		,double my_PA           /// Position angle (radians)
		,double my_index  /// Sersic index
		,double my_q     // axes ratio
		,double my_z            /// optional redshift
		,const double *my_theta          /// optional angular position on the sky
		)
:Source()
{
	setInternals(my_mag,my_Reff,my_PA,my_index,my_q,my_z,my_theta);
	sb_limit = 30.;
}

SersicSource::~SersicSource()
{
}

void SersicSource::setInternals(double my_mag,double my_Reff,double my_PA,double my_index,double my_q,double my_z,const double *my_theta){

	Reff = my_Reff*pi/180/60/60;
	mag = my_mag;
	PA = my_PA;
	index = my_index;
	q = my_q;

	zsource = my_z;
	if(my_theta != NULL){
		source_x[0] = my_theta[0];
		source_x[1] = my_theta[1];
	}else{
		source_x[0] = 0;
		source_x[1] = 0;
	}

	// approximation valid for 0.5 < n < 8
	bn = 1.9992*index - 0.3271;
	flux = pow(10,-0.4*(mag+48.6));
	Ieff = flux/2./pi/Reff/Reff/exp(bn)/index*pow(bn,2*index)/tgamma(2*index)/q;
}

double SersicSource::SurfaceBrightness(
		double *x  /// position in radians relative to center of source
		){

	double x_new[2];
	x_new[0] = (x[0]-source_x[0])*cos(PA)-(x[1]-source_x[1])*sin(PA);
	x_new[1] = (x[0]-source_x[0])*sin(PA)+(x[1]-source_x[1])*cos(PA);

	double r = sqrt(x_new[0]*x_new[0]+x_new[1]*x_new[1]/q/q);

	double sb = Ieff * exp(-bn*(pow(r/Reff,1/index)-1))/hplanck;
	if (sb*hplanck< pow(10,-0.4*(48.6+sb_limit))*pow(180*60*60/pi,2)) return 0.;
	return sb;
}

void SersicSource::printSource(){}
void SersicSource::assignParams(InputParams& /* params */){}
