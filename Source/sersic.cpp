/*
 * Source/sersic.cpp
 *
 *  Created on: Feb 28, 2010
 *      Author: R.B. Metcalf
 */
#include "slsimlib.h"

SourceSersic::SourceSersic()
: Source()
{
  /// set to values that hopefully will cause an error if it is used
  ReSet(1.0e10,1,0,0,-1,-1,0);
}

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
  assert(my_Reff > 0);
  ReSet(my_mag,my_Reff,my_PA,my_index,my_q,my_z,my_theta);
}
SourceSersic::SourceSersic(const SourceSersic &p){
  Reff = p.Reff;
  mag = p.mag;
  PA = p.PA;
  index = p.index;
  bn = p.bn;
  q = p.q;
  flux = p.flux;
  I_r = p.I_r;
  I_n = p.I_n;
  I_q = p.I_q;
  cosPA = p.cosPA;
  sinPA = p.sinPA;
}
SourceSersic& SourceSersic::operator=(const SourceSersic &p){
  if(this == &p) return *this;
  
  Reff = p.Reff;
  mag = p.mag;
  PA = p.PA;
  index = p.index;
  bn = p.bn;
  q = p.q;
  flux = p.flux;
  I_r = p.I_r;
  I_n = p.I_n;
  I_q = p.I_q;
  cosPA = p.cosPA;
  sinPA = p.sinPA;
  
  return *this;
}

/// Reset all the parameters
void SourceSersic::ReSet(
                               double my_mag            /// Total magnitude
                               ,double my_Reff          /// Bulge half light radius (arcs)
                               ,double my_PA            /// Position angle (radians)
                               ,double my_index         /// Sersic index
                               ,double my_q             /// axes ratio
                               ,double my_z             /// redshift
                               ,const double *my_theta  /// optional angular position on the sky
){
//std::cout << "SourceSersic constructor" << std::endl;
  //assert(my_Reff > 0);
  setReff(my_Reff);
  setMag(my_mag);
  setPA(my_PA);
  setSersicIndex(my_index);
  setAxesRatio(my_q);

  setZ(my_z);

  if(my_theta)
    setTheta(my_theta[0], my_theta[1]);

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

  if(Reff <= 0.0) return 0.0;
	PosType x_new[2];
	x_new[0] = (x[0]-source_x[0])*cosPA+(x[1]-source_x[1])*sinPA;
	x_new[1] = (x[0]-source_x[0])*sinPA-(x[1]-source_x[1])*cosPA;

	PosType r = sqrt(x_new[0]*x_new[0]+x_new[1]*x_new[1]/q/q);

	PosType sb = flux * I_n * I_q * I_r *inv_hplanck * exp(-bn*pow(r/Reff,1./index));
  
  if (sb < sb_limit) return 0.;
  assert(sb >= 0.0);
	return sb;
}

void SourceSersic::printSource(){}
void SourceSersic::assignParams(InputParams& /* params */){}
