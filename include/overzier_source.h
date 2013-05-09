/*
 * overzier_source.h
 *
 *  Created on: Mar 6, 2010
 *      Author: R.B. Metcalf
 */
#ifndef OVERZIER_SOURCE_H_
#define OVERZIER_SOURCE_H_

#include "source.h"

// define pi here if not done via include
#ifndef pi
#define pi 3.141592653589793238462643383279502884
#endif

/**
 *\brief Structure for holding parameters for one or more galaxy images according to
 * the Overzier model.
 */
class OverzierSource : public Source
{
public:
	SOURCE_TYPE(OverzierSource)
	
	OverzierSource();
	OverzierSource(double mag,double BtoT,double Reff,double Rh,double PA,double inclination,unsigned long my_id,double my_z=0,const double *theta=0);
	~OverzierSource();
	
	void setInternals(double mag,double BtoT,double Reff,double Rh,double PA,double inclination,unsigned long my_id,double my_z=0,const double *my_theta=0);
	double SurfaceBrightness(double *x);
	double getTotalFlux();
	void printSource();
	
	/// Halo ID.
	unsigned long getID() { return haloID; }
	
	/// get magnitude of whole galaxy.  Which band this is in depends on which was passed in the constructor
	double getMag() const { return mag; }
	double getUMag() const { return mag_u; }
	double getGMag() const { return mag_g; }
	double getRMag() const { return mag_r; }
	double getIMag() const { return mag_i; }
	double getZMag() const { return mag_z; }
	double getJMag() const { return mag_J; }
	double getHMag() const { return mag_H; }
	double getKMag() const { return mag_Ks; }
	
	/// set u band magnitude
	void setUMag(double m) { mag_u = m; }
	/// set g band magnitude
	void setGMag(double m) { mag_g = m; }
	/// set r band magnitude
	void setRMag(double m) { mag_r = m; }
	/// set i band magnitude
	void setIMag(double m) { mag_i = m; }
	/// set z band magnitude
	void setZMag(double m) { mag_z = m; }
	/// set j band magnitude
	void setJMag(double m) { mag_J = m; }
	/// set h band magnitude
	void setHMag(double m) { mag_H = m; }
	/// set k band magnitude
	void setKMag(double m) { mag_Ks = m; }
	
	/// bulge half light radius
	double getReff() const { return Reff/(pi/180/60/60); }
	/// disk scale height
	double getRh() const { return Rh/(pi/180/60/60); }
	
	double getBtoT() const { return BtoT; }
	double getPA() const { return PA; }
	double getInclination() const { return inclination; }
	
	/** Returns minimum of the radii at which disk and bulge have a surf. brightness equal to a fraction f of the central one
	* TODO Fabio: Needs to be tested and improved (Bulge is so steep in the center that output values are very small)
  */
	inline double getMinSize(double f) {return std::min(1.678*Reff*fabs(cos(inclination))*pow(-log (f)/7.67,4),Rh*(-log (f)/1.67));}


private:
	void assignParams(InputParams& params);
	
	/// haloID
	unsigned long haloID;
	
	/// weighted mean between the radii that enclose 99% of the flux
	/// in the pure De Vacouleur/exponential disk case
	/// 6.670 = 3.975*Re = 3.975*1.678*Rh
	
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
	
	// optional position variables
};

#endif /* GALAXIES_OVERZIER_H_ */
