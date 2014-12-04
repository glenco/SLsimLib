/*
 * Source/sersic.cpp
 *
 *  Created on: Feb 28, 2010
 *      Author: R.B. Metcalf
 */
#include "slsimlib.h"

SourceSersic::SourceSersic(
	double my_mag            /// Total magnitude
	,double my_Reff          /// Bulge half light radius (arcs)
	,double my_PA            /// Position angle (radians)
	,double my_index         /// Sersic index
	,double my_q             /// axes ratio
	,double my_z             /// redshift
	,const double *my_theta  /// optional angular position on the sky
)
: Source()
{
	setReff(my_Reff);
	setMag(my_mag);
	setPA(my_PA);
	setSersicIndex(my_index);
	setAxesRatio(my_q);
	
	setZ(my_z);
	
	if(my_theta)
		setX(my_theta[0], my_theta[1]);

	if(q > 1)
		throw std::invalid_argument("Error: q must be < 1!");
  

}

SourceSersic::~SourceSersic()
{
}

PosType SourceSersic::SurfaceBrightness(
	PosType *x  /// position in radians relative to center of source
)
{

	PosType x_new[2];
	x_new[0] = (x[0]-source_x[0])*cos(PA)+(x[1]-source_x[1])*sin(PA);
	x_new[1] = (x[0]-source_x[0])*sin(PA)-(x[1]-source_x[1])*cos(PA);

	PosType r = sqrt(x_new[0]*x_new[0]+x_new[1]*x_new[1]/q/q);

	PosType sb = flux * I_n * I_q * I_r * exp(-bn*pow(r/Reff,1./index))*inv_hplanck;
	if (sb < sb_limit) return 0.;
	return sb;
}

void SourceSersic::printSource(){}
void SourceSersic::assignParams(InputParams& /* params */){}
