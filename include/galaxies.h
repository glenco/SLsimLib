/*
 * galaxies.h
 *
 *  Created on: Mar 6, 2010
 *      Author: R.B. Metcalf
 */
#ifndef GALAXIES_H_
#define GALAXIES_H_

/**
 *\brief Structure for holding parameters for one or more galaxy images according to
 * the Overzier model.
 */
struct OverGalaxy{

	OverGalaxy(){};
	OverGalaxy(double mag,double BtoT,double Reff,double Rh,double PA,double inclination,double my_z=0,double *theta=NULL);
	~OverGalaxy(){};

	void setInternals(double mag,double BtoT,double Reff,double Rh,double PA,double inclination,double my_z=0,double *my_theta=NULL);
	double SurfaceBrightness(double *x);
	void print();

	/// redshift
	double z;
	/// position on the sky
	double theta[2];
private:
	/// bulge half light radius
	double Reff;
	/// disk scale hight
	double Rh;

	double cxx,cyy,cxy;
	/// internal valuable mag-2.5*log10(1-BtoT)+5*log10(Rh)+1.9955
	double sbDo;
	/// internal valuable mag-2.5*log10(BtoT)+5*log10(Reff)-4.93884
	double sbSo;

	// optional position variables
};

void create_sersic(int n,double Ro,double f,double *center,double theta,double **x,long Nsources);

#endif /* GALAXIES_H_ */
