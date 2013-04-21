/*
 * sersic_source.h
 *
 *  Created on: Mar 6, 2010
 *      Author: R.B. Metcalf
 */
#ifndef SERSIC_SOURCE_H_
#define SERSIC_SOURCE_H_

#include "source.h"

class SersicSource : public Source
{
public:
	SOURCE_TYPE(SersicSource)
	
	SersicSource();
	SersicSource(double mag,double Reff,double PA,double my_index,double my_q,double my_z=0,const double *theta=0);
	~SersicSource();
	
	void setParameters(Parameters& p);
	void getParameters(Parameters& p);
	
	void randomize(double step, long* seed);
	
	void setInternals();
	
	inline double getMag() { return mag; }
	// approximation of the radius that encloses 99% of the flux
	inline double getPA() { return PA; }
	
	inline double getReff() const { return Reff*180*60*60/pi; }
	inline double getAxesRatio() const { return q; }
	inline double getSersicIndex() const { return index; }
	
	double SurfaceBrightness(double *x);
	inline double getTotalFlux() {return flux;}
	void printSource();
	
private:
	void assignParams(InputParams& params);
	double Reff;
	double mag;
	double PA;
	double index;
	double bn;
	double q;
	double Ieff;
	double flux;
};

#endif /* GALAXIES_SERSIC_H_ */
