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
	OverGalaxy(double mag,double BtoT,double Reff,double Rh,double PA,double inclination,unsigned long my_id,double my_z=0,const double *theta=0);
	~OverGalaxy();

	void setInternals(double mag,double BtoT,double Reff,double Rh,double PA,double inclination,unsigned long my_id,double my_z=0,const double *my_theta=0);
	double SurfaceBrightness(double *x);
	void print();
	double getMag() const { return mag; }
	/// bulge half light radius
	double getReff() const { return Reff/(pi/180/60/60); }
	/// disk scale height
	double getRh() const { return Rh/(pi/180/60/60); }
	double getBtoT() const { return BtoT; }
	double getPA() const { return PA; }
	double getInclination() const { return inclination; }

	/// haloID
	unsigned long haloID;
	/// redshift
	double z;
	/// position on the sky
	double theta[2];
	/// returns the maximum radius of the source galaxy TODO This needs to be done better.
	double getRadius(){return /*6*(Reff > Rh ? Reff : Rh)*/
		// weighted mean between the radii that enclose 99% of the flux
		// in the pure De Vacouleur/exponential disk case
		// 6.670 = 3.975*Re = 3.975*1.678*Rh
		6.670*Rh*(1-BtoT)+18.936*Reff*BtoT;}
private:
	/// bulge half light radius
	double Reff;
	/// disk scale height
	double Rh;

	double BtoT;
	double PA;
	double inclination;

	double cxx,cyy,cxy;
	double sbDo;
	double sbSo;
	double mag;

	// optional position variables
};

#endif /* GALAXIES_OVERZIER_H_ */
