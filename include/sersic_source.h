/*
 * sersic_source.h
 *
 *  Created on: Mar 6, 2010
 *      Author: R.B. Metcalf
 */
#ifndef SERSIC_SOURCE_H_
#define SERSIC_SOURCE_H_

#include "source.h"

/** \brief Class for sources described by a Sersic profile
 *
 * The sources are constructed from magnitude, half light radius Reff, Sersic index, axis ratio q and main axis orientation PA.
 * For q = 1, the source is circular and the surface brightness profile is I(r) = Ieff * exp(-bn*((r/Reff)^(1/index)-1))
 * The elliptical model is defined for axis ratio q < 1.
 *
 */
class SourceSersic : public Source
{
public:
  /// sets values to invalid values
  SourceSersic();
  SourceSersic(PosType mag,PosType Reff,PosType PA,PosType my_index,PosType my_q,PosType my_z,const PosType *theta=0);
	~SourceSersic();
	
  void ReSet(PosType mag,PosType Reff,PosType PA,PosType my_index,PosType my_q,PosType my_z,const PosType *theta=0);

	/// calculates radius where the surface brightness drops by a factor f with respect to the central peak
	inline PosType FractionRadius (PosType f) {return Reff*pow(-log (f)/bn,index);}
	
	inline PosType getSersicIndex() const { return index; }
	inline PosType getAxesRatio() const { return q; }
	inline PosType getReff() const { return Reff*180*60*60/pi; }
	inline PosType getMag() const { return mag; }
	inline PosType getPA() const { return PA; }
	
	inline void setSersicIndex(PosType x)
	{
		index = x;
		bn = 1.9992*index - 0.3271; // approximation valid for 0.5 < n < 8
		I_n = 1./index*pow(bn,2*index)/tgamma(2*index);
		updateRadius();
	}
	
	inline void setAxesRatio(PosType x)
	{
		q = x;
		I_q = 1./q;
	}
	
	inline void setReff(PosType x)
	{
		Reff = x*pi/180/60/60;
		I_r = 1./2./pi/Reff/Reff;
		updateRadius();
	}
	
	inline void setMag(PosType x)
	{
		mag = x;
		flux = pow(10, -0.4*(mag+48.6));
	}
	
	inline void setPA(PosType x)
	{
		PA = x;
    cosPA = cos(x);
    sinPA = sin(x);
	}
	
	inline PosType getTotalFlux() const { return flux; }
	
	PosType SurfaceBrightness(PosType *x);
	void printSource();
	
private:
	inline void updateRadius()
	{
		// radius in Source
		// approximation of the radius that encloses 99% of the flux
		setRadius((3.73 - 0.926*index + 1.164*index*index)*Reff);
	}
	
	void assignParams(InputParams& params);
	PosType Reff;
	PosType mag;
	PosType PA;
	PosType index;
	PosType bn;
	PosType q;
	PosType flux;
	PosType I_r, I_n, I_q;
  PosType cosPA,sinPA;
};

#endif /* GALAXIES_SERSIC_H_ */
