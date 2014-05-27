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
class SourceOverzier : public Source
{
public:
	SourceOverzier();
	SourceOverzier(PosType mag,PosType BtoT,PosType Reff,PosType Rh,PosType PA,PosType inclination,unsigned long my_id,PosType my_z=0,const PosType *theta=0);
	~SourceOverzier();
	
	void setInternals(PosType mag,PosType BtoT,PosType Reff,PosType Rh,PosType PA,PosType inclination,unsigned long my_id,PosType my_z=0,const PosType *my_theta=0);
	PosType SurfaceBrightness(PosType *x);
	PosType getTotalFlux();
	void printSource();
	
	/// Halo ID.
	unsigned long getID() { return haloID; }
	
	/// get magnitude of whole galaxy.  Which band this is in depends on which was passed in the constructor
	PosType getMag() const { return mag; }
	PosType getUMag() const { return mag_u; }
	PosType getGMag() const { return mag_g; }
	PosType getRMag() const { return mag_r; }
	PosType getIMag() const { return mag_i; }
	PosType getZMag() const { return mag_z; }
	PosType getJMag() const { return mag_J; }
	PosType getHMag() const { return mag_H; }
	PosType getKMag() const { return mag_Ks; }
	
	/// set u band magnitude
	void setUMag(PosType m) { mag_u = m; }
	/// set g band magnitude
	void setGMag(PosType m) { mag_g = m; }
	/// set r band magnitude
	void setRMag(PosType m) { mag_r = m; }
	/// set i band magnitude
	void setIMag(PosType m) { mag_i = m; }
	/// set z band magnitude
	void setZMag(PosType m) { mag_z = m; }
	/// set j band magnitude
	void setJMag(PosType m) { mag_J = m; }
	/// set h band magnitude
	void setHMag(PosType m) { mag_H = m; }
	/// set k band magnitude
	void setKMag(PosType m) { mag_Ks = m; }
	
	/// bulge half light radius in radians
	PosType getReff() const { return Reff/(pi/180/60/60); }
	/// disk scale height in radians
	PosType getRh() const { return Rh/(pi/180/60/60); }
	
	PosType getBtoT() const { return BtoT; }
	PosType getPA() const { return PA; }
	PosType getInclination() const { return inclination; }
	
	/** Returns minimum of the radii at which disk and bulge have a surf. brightness equal to a fraction f of the central one
	* TODO: Fabio: Needs to be tested and improved (Bulge is so steep in the center that output values are very small)
  */
	inline PosType getMinSize(PosType f) {return std::min(1.678*Reff*fabs(cos(inclination))*pow(-log (f)/7.67,4),Rh*(-log (f)/1.67));}

  static PosType *getx(SourceOverzier &sourceo){return sourceo.getX();}

private:
	void assignParams(InputParams& params);
	
	/// haloID
	unsigned long haloID;

	/// bulge half light radius
	PosType Reff;
	/// disk scale height
	PosType Rh;
	
	PosType BtoT;
	PosType PA;
	PosType inclination;
	
	PosType cxx,cyy,cxy;
	PosType sbDo;
	PosType sbSo;
	PosType mag;
	
	// colors
	PosType mag_u;
	PosType mag_g;
	PosType mag_r;
	PosType mag_i;
	PosType mag_z;
	PosType mag_J;
	PosType mag_H;
	PosType mag_Ks;
	PosType mag_i1;
	PosType mag_i2;
	
	// optional position variables
};

#endif /* GALAXIES_OVERZIER_H_ */
