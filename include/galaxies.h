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
	OverGalaxy(double mag,double BtoT,double Reff,double Rh,double PA,double inclination,unsigned long my_id,double my_z=0,double *theta=NULL);
	~OverGalaxy(){};

	void setInternals(double mag,double BtoT,double Reff,double Rh,double PA,double inclination,unsigned long my_id,double my_z=0,double *my_theta=NULL);
	double SurfaceBrightness(double *x);
	void print();
	double getMag(){return mag;}
	/// bulge half light radius
	double getReff(){return Reff;}
	/// disk scale height
	double getRh(){return Rh;}
	double getBtoT(){return  Reff*Reff*sbSo*pow(10,mag/2.5)/94.484376;}


	/// haloID
	unsigned long haloID;
	/// redshift
	double z;
	/// position on the sky
	double theta[2];
	/// returns the maximum radius of the source galaxy TODO This needs to be done better.
	double getRadius(){return 6*(Reff > Rh ? Reff : Rh);}
private:
	/// bulge half light radius
	double Reff;
	/// disk scale height
	double Rh;

	double cxx,cyy,cxy;
	double sbDo;
	double sbSo;
	double mag;

	// optional position variables
};

void create_sersic(int n,double Ro,double f,double *center,double theta,double **x,long Nsources);

#endif /* GALAXIES_H_ */
