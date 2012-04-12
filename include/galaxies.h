/*
 * galaxies.h
 *
 *  Created on: Mar 6, 2010
 *      Author: R.B. Metcalf
 */
#ifndef GALAXIES_H_
#define GALAXIES_H_

/**
 *
 * Structure for holding parameters that define the image of galaxy according to
 * the Overzier model.
 */
struct OverGalaxy{

	OverGalaxy(){};
	OverGalaxy(double mag,double BtoT,double Reff,double Rh,double PA,double inclination);
	~OverGalaxy();

	void setInternals(double mag,double BtoT,double Reff,double Rh,double PA,double inclination);
	double SurfaceBrightness(double *x);

private:
	/// bulge half light radius
	double Reff;
	/// disk scale hight
	double Rh;
	/// position angle
	double PA;
	/// inclination
	double incl;
	/// internal valuable mag-2.5*log10(1-BtoT)+5*log10(Rh)+1.9955
	double muDo;
	/// internal valuable mag-2.5*log10(BtoT)+5*log10(Reff)-4.93884
	double muSo;
};

void create_sersic(int n,double Ro,double f,double *center,double theta,double **x,long Nsources);

#endif /* GALAXIES_H_ */
