/*
 * sersic_source.h
 *
 *  Created on: Mar 6, 2010
 *      Author: R.B. Metcalf
 */
#ifndef SERSIC_SOURCE_H_
#define SERSIC_SOURCE_H_

#include "source.h"

class SourceSersic : public Source
{
public:
	SOURCE_TYPE(SourceSersic)
	
	SersicSource();
	SersicSource(double mag,double Reff,double PA,double my_index,double my_q,double my_z=0,const double *theta=0);
	SersicSource(InputParams&){};
	~SersicSource();
	
	void getParameters(Parameters& p) const;
	void setParameters(Parameters& p);
	
	void randomize(double step, long* seed);
	
	void setInternals();
	
	inline double getMag() { return mag; }

	// calculates radius where the surface brightness drops by a factor f with respect to the central peak
	inline double FractionRadius (double f) {return Reff*pow(-log (f)/bn,index);}

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
