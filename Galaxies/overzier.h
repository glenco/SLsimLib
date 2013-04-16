/*
 * Galaxies/overzier.h
 *
 *  Created on: Mar 6, 2010
 *      Author: R.B. Metcalf
 */
#ifndef GALAXIES_OVERZIER_H_
#define GALAXIES_OVERZIER_H_

/**
 *\brief Structure for holding parameters for one or more galaxy images according to
 * the Overzier model.
 */
struct OverGalaxy{
	OverGalaxy();
	OverGalaxy(double mag,double BtoT,double Reff,double Rh,double PA,double inclination,unsigned long my_id,double my_z=0,double *theta=0);
	~OverGalaxy();

	void setInternals(double mag,double BtoT,double Reff,double Rh,double PA,double inclination,unsigned long my_id,double my_z=0,double *my_theta=0);
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


	// colors
	double mag_u;
	double mag_g;
	double mag_r;
	double mag_i;
	double mag_z;
	double mag_J;
	double mag_H;
	double mag_Ks;
	double mag_i1;
	double mag_i2;

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

#endif /* GALAXIES_OVERZIER_H_ */
