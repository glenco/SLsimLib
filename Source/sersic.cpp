/*
 * Source/sersic.cpp
 *
 *  Created on: Feb 28, 2010
 *      Author: R.B. Metcalf
 */
#include "slsimlib.h"

SourceSersic::SourceSersic(
		double my_mag              /// Total magnitude
		,double my_Reff         /// Bulge half light radius (arcs)
		,double my_PA           /// Position angle (radians)
		,double my_index  /// Sersic index
		,double my_q     // axes ratio
		,double my_z            /// optional redshift
		,const double *my_theta          /// optional angular position on the sky
		)
: Source(),
  Reff(my_Reff*pi/180/60/60), mag(my_mag), PA(my_PA), index(my_index), q(my_q)
{
	setZ(my_z);
	
	if(my_theta)
		setX(my_theta[0], my_theta[1]);
	
	if ( q > 1)
	{
		std::cerr << "Error: q must be < 1!" << std::endl;
		exit(1);
	}
	setInternals();
}

SourceSersic::SourceSersic(InputParams& params)
: Source()
{
}

SourceSersic::SourceSersic(): Source()
{
}

SourceSersic::~SourceSersic()
{
}

std::size_t SourceSersic::Nrandomize() const
{
	return Source::Nrandomize() + 5;
}

Utilities::Any SourceSersic::randomize(std::size_t i, double step, long* seed)
{
	Utilities::Any old;
	
	if(i < Source::Nrandomize())
	{
		old = Source::randomize(i, step, seed);
	}
	else
	{
		switch(i - Source::Nrandomize())
		{
			case 0:
				// half light radius
				old = Reff;
				Reff = std::max(1.e-10, Reff + step*pi/180/60/60*gasdev(seed));
				break;
			case 1:
				// magnitude
				old = mag;
				mag += step*gasdev(seed);
				break;
			case 2:
				// position angle
				old = PA;
				PA = std::fmod(PA + step*pi*gasdev(seed), pi) - pi/2;
				break;
			case 3:
				// Sersic index
				old = index;
				index += step*gasdev(seed);
				break;
			case 4:
				// axes ratio
				old = q;
				q = std::min(1., std::max(1.e-05, q + step*gasdev(seed)));
				break;
			default:
				throw std::invalid_argument("bad parameter index for randomize()");
		}
	}
	
	// update
	setInternals();
	
	return old;
}

void SourceSersic::unrandomize(std::size_t i, const Utilities::Any& old)
{
	if(i < Source::Nrandomize())
	{
		Source::unrandomize(i, old);
	}
	else
	{
		switch(i - Source::Nrandomize())
		{
			case 0:
				// half light radius
				Reff = Utilities::AnyCast<double>(old);
				break;
			case 1:
				// magnitude
				mag = Utilities::AnyCast<double>(old);
				break;
			case 2:
				// position angle
				PA = Utilities::AnyCast<double>(old);
				break;
			case 3:
				// Sersic index
				index = Utilities::AnyCast<double>(old);
				break;
			case 4:
				// axes ratio
				q = Utilities::AnyCast<double>(old);
				break;
			default:
				throw std::invalid_argument("bad parameter index for randomize()");
		}
	}
	
	// update
	setInternals();
}

void SourceSersic::setInternals()
{
	// approximation valid for 0.5 < n < 8
	bn = 1.9992*index - 0.3271;
	flux = pow(10,-0.4*(mag+48.6));
	Ieff = flux/2./pi/Reff/Reff/exp(bn)/index*pow(bn,2*index)/tgamma(2*index)/q;
	
	// radius in Source
	// approximation of the radius that encloses 99% of the flux
	setRadius((3.73 - 0.926*index + 1.164*index*index)*Reff);
}

double SourceSersic::SurfaceBrightness(
		double *x  /// position in radians relative to center of source
		){

	double x_new[2];
	x_new[0] = (x[0]-source_x[0])*cos(PA)+(x[1]-source_x[1])*sin(PA);
	x_new[1] = (x[0]-source_x[0])*sin(PA)-(x[1]-source_x[1])*cos(PA);

	double r = sqrt(x_new[0]*x_new[0]+x_new[1]*x_new[1]/q/q);

	double sb = Ieff * exp(-bn*(pow(r/Reff,1./index)-1.))*inv_hplanck;
	if (sb < sb_limit) return 0.;
	return sb;
}

void SourceSersic::printSource(){}
void SourceSersic::assignParams(InputParams& /* params */){}
