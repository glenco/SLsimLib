/*
 * overzier_source.h
 *
 *  Created on: Mar 6, 2010
 *      Author: R.B. Metcalf
 */
#ifndef OVERZIER_SOURCE_H_
#define OVERZIER_SOURCE_H_

#include "source.h"
#include "sersic_source.h"

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
	SourceOverzier(PosType mag,PosType mag_bulge,PosType Reff,PosType Rh,PosType PA,PosType inclination,unsigned long my_id,PosType my_z=0,const PosType *theta=0);
  
  SourceOverzier(const SourceOverzier &s);
  SourceOverzier& operator=(const SourceOverzier &s);
	virtual ~SourceOverzier();
	
	void setInternals(PosType mag,PosType BtoT,PosType Reff,PosType Rh,PosType PA,PosType inclination,unsigned long my_id,PosType my_z=0,const PosType *my_theta=0);
  virtual PosType SurfaceBrightness(PosType *x);
	PosType getTotalFlux() const;
	void printSource();
	
	/// Halo ID.
	unsigned long getID() { return haloID; }
	
	/// get magnitude of whole galaxy.  Which band this is in depends on which was passed in the constructor
  PosType getMag() const { return mag; }
  PosType getMag(Band band) const ;
  PosType getMagBulge() const { return mag_bulge; }
  PosType getMagBulge(Band band) const;
  
  /*
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
  */
  
  /// magnitude in specific band
  void setMag(Band band,PosType my_mag);
  
  /// magnitude in specific band
  void setMagBulge(Band band,PosType my_mag);
  
	/// bulge half light radius in radians
	PosType getReff() const { return Reff/(pi/180/60/60); }
	/// disk scale height in radians
	PosType getRh() const { return Rh/(pi/180/60/60); }
	
  /// the bulge to total flux ratio
	PosType getBtoT() const { return pow(10,(-mag_bulge + mag)/2.5); }
	PosType getPA() const { return PA; }
	PosType getInclination() const { return inclination;}
  float getSEDtype() const {return sedtype;}
  void setSEDtype(float s){ sedtype = s;}

  /// change the working band
  virtual void changeBand(Band band);
	
	/** Returns minimum of the radii at which disk and bulge have a surf. brightness equal to a fraction f of the central one
	* TODO: Fabio: Needs to be tested and improved (Bulge is so steep in the center that output values are very small)
  */
	inline PosType getMinSize(PosType f) {return std::min(1.678*Reff*fabs(cos(inclination))*pow(-log (f)/7.67,4),Rh*(-log (f)/1.67));}

  static PosType *getx(SourceOverzier &sourceo){return sourceo.source_x.x;}

protected:
  
  float sedtype = -1;
  // renormalize the disk and bulge to agree with current mag and mag_bulge
  void renormalize();
	void assignParams(InputParams& params);
	
	/// haloID
	unsigned long haloID;

	/// bulge half light radius
	PosType Reff;
	/// disk scale height
	PosType Rh;
	
	//PosType BtoT;
	PosType PA;
	PosType inclination;
	
	PosType cxx,cyy,cxy;
	PosType sbDo;
	PosType sbSo;
	PosType mag;
  PosType mag_bulge;
	
  // colors
  std::map<Band,double> mag_map;
  std::map<Band,double> bulge_mag_map;

	// optional position variables
};

class SourceOverzierPlus : public SourceOverzier
{
public:
  //SourceOverzierPlus();
  SourceOverzierPlus(PosType mag,PosType BtoT,PosType Reff,PosType Rh,PosType PA,PosType inclination,unsigned long my_id,PosType my_z,const PosType *theta,Utilities::RandomNumbers_NR &ran);
  ~SourceOverzierPlus();
  
  SourceOverzierPlus(const SourceOverzierPlus &p);
  SourceOverzierPlus & operator=(const SourceOverzierPlus &p);
  
  //*** meed to be able to change band
  //*** put modes and phases into surface brightness
  //*** possible put gaussian texture on disk

  PosType SurfaceBrightness(PosType *y);

  int getNarms() const {return Narms;}
  PosType getArmAmplitude() const {return Ad;}
  PosType getArmAlpha() const {return arm_alpha;}
  PosType getSphIndex() const {return spheroid->getSersicIndex();}
  PosType getSphAxisRatio() const {return spheroid->getAxesRatio();}
  PosType getSphPA() const {return spheroid->getPA();}
  
  void changeBand(Band band);
  static PosType* getx(SourceOverzierPlus &sourceo){return sourceo.source_x.x;}

  /// Reset the position of the source in radians
  virtual inline void setTheta(PosType *xx){
    source_x[0] = xx[0];
    source_x[1] = xx[1];
    spheroid->setTheta(xx);
  }
  virtual void setTheta(PosType my_x,PosType my_y){
    source_x[0] = my_x;
    source_x[1] = my_y;
    spheroid->setTheta(my_x,my_y);
  }
  void setBulgeAxisRatio(PosType q){
    spheroid->setAxesRatio(q);
  }
  /// Randomly change some of the internal paramters and angles of the source
  void randomize(Utilities::RandomNumbers_NR &ran);
private:
  int Narms;
  PosType Ad,mctalpha,arm_alpha;
  SourceSersic *spheroid;
  std::vector<PosType> modes;
  PosType disk_phase;
  PosType cospa,sinpa,cosi;
};
#endif /* GALAXIES_OVERZIER_H_ */
