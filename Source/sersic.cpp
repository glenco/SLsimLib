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
: Source(),
  mag(my_mag), Reff(my_Reff*pi/180/60/60), PA(my_PA), index(my_index), q(my_q)
{
	setZ(my_z);
	
	if(my_theta)
		setX(my_theta[0], my_theta[1]);
	
	setInternals();
}

SersicSource::~SersicSource()
{
}

void SersicSource::getParameters(Parameters& p) const
{
	// base class serialization
	Source::getParameters(p);
	
	p << Reff << mag << PA << index << bn << q << Ieff << flux;
}

void SersicSource::setParameters(Parameters& p)
{
	// base class deserialization
	Source::setParameters(p);
	
	p >> Reff >> mag >> PA >> index >> bn >> q >> Ieff >> flux;
}

void SersicSource::randomize(double step, long* seed)
{
	// half light radius
	Reff += step*pi/180/60/60*gasdev(seed);
	
	// magnitude
	mag += step*gasdev(seed);
	
	// position angle
	//PA += step*pi*gasdev(seed);
	
	// Sersic index
	index += step*gasdev(seed);
	
	// axes ratio
	//q += step*gasdev(seed);
	
	// redshift?
	
	// position
	source_x[0] += step*pi/180/60/60*gasdev(seed);
	source_x[1] += step*pi/180/60/60*gasdev(seed);
	
	// update
	setInternals();
}

void SersicSource::setInternals()
{
	// approximation valid for 0.5 < n < 8
	bn = 1.9992*index - 0.3271;
	flux = pow(10,-0.4*(mag+48.6));
	Ieff = flux/2./pi/Reff/Reff/exp(bn)/index*pow(bn,2*index)/tgamma(2*index)/q;
	
	// radius in Source
	setRadius((3.73 - 0.926*index + 1.164*index*index)*Reff);
}

double SersicSource::SurfaceBrightness(
		double *x  /// position in radians relative to center of source
		){

	double x_new[2];
	x_new[0] = (x[0]-source_x[0])*cos(PA)-(x[1]-source_x[1])*sin(PA);
	x_new[1] = (x[0]-source_x[0])*sin(PA)+(x[1]-source_x[1])*cos(PA);

	double r = sqrt(x_new[0]*x_new[0]+x_new[1]*x_new[1]/q/q);

	double sb = Ieff * exp(-bn*(pow(r/Reff,1./index)-1.))*inv_hplanck;
	if (sb < sb_limit) return 0.;
	return sb;
}

void SersicSource::printSource(){}
void SersicSource::assignParams(InputParams& /* params */){}
