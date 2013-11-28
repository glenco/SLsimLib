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
 * The sources are constructed from magnitude, effective radius Reff, Sersic index, axis ratio q and main axis orientation PA.
 * For q = 1, the source is circular and the surface brightness profile is I(r) = Ieff * exp(-bn*((r/Reff)^(1/index)-1))
 * The elliptical model is defined for axis ratio q < 1.
 *
 */
class SourceSersic : public Source
{
public:
	SourceSersic(double mag,double Reff,double PA,double my_index,double my_q,double my_z=0,const double *theta=0);
	~SourceSersic();
	
	/// calculates radius where the surface brightness drops by a factor f with respect to the central peak
	inline double FractionRadius (double f) {return Reff*pow(-log (f)/bn,index);}
	
	inline double getSersicIndex() const { return index; }
	inline double getAxesRatio() const { return q; }
	inline double getReff() const { return Reff*180*60*60/pi; }
	inline double getMag() { return mag; }
	inline double getPA() { return PA; }
	
	inline void setSersicIndex(double x)
	{
		index = x;
		bn = 1.9992*index - 0.3271; // approximation valid for 0.5 < n < 8
		I_n = 1./index*pow(bn,2*index)/tgamma(2*index);
		updateRadius();
	}
	
	inline void setAxesRatio(double x)
	{
		q = x;
		I_q = 1./q;
	}
	
	inline void setReff(double x)
	{
		Reff = x*pi/180/60/60;
		I_r = 1./2./pi/Reff/Reff;
		updateRadius();
	}
	
	inline void setMag(double x)
	{
		mag = x;
		flux = pow(10, -0.4*(mag+48.6));
	}
	
	inline void setPA(double x)
	{
		PA = x;
	}
	
	inline double getTotalFlux() { return flux; }
	
	double SurfaceBrightness(double *x);
	void printSource();
	
private:
	inline void updateRadius()
	{
		// radius in Source
		// approximation of the radius that encloses 99% of the flux
		setRadius((3.73 - 0.926*index + 1.164*index*index)*Reff);
	}
	
	void assignParams(InputParams& params);
	double Reff;
	double mag;
	double PA;
	double index;
	double bn;
	double q;
	double flux;
	double I_r, I_n, I_q;
};

#endif /* GALAXIES_SERSIC_H_ */
