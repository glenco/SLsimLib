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

// define PI here if not done via include
#ifndef PI
#define PI 3.141592653589793238462643383279502884
#endif

/**
 *\brief Structure for holding parameters for one or more galaxy images according to
 * the Overzier model.
 */
class SourceOverzier : public Source
{
public:
	//SourceOverzier();
	SourceOverzier(PosType mag,PosType mag_bulge,PosType Reff,PosType Rdisk,PosType PA,PosType inclination,unsigned long my_id,PosType my_z=0,const PosType *theta=0);
  
  SourceOverzier(const SourceOverzier &s);
  SourceOverzier& operator=(const SourceOverzier &s);
	virtual ~SourceOverzier();
	
	void setInternals(PosType mag,PosType BtoT,PosType Reff,PosType Rdisk,PosType PA,PosType inclination,unsigned long my_id,PosType my_z=0,const PosType *my_theta=0);
  virtual PosType SurfaceBrightness(PosType *x);
	PosType getTotalFlux() const;
	void printSource();
	
	/// Halo ID.
	unsigned long getID() { return haloID; }
	
	/// get magnitude of whole galaxy.  Which band this is in depends on which was passed in the constructor
  PosType getMag() const { return current.mag; }
  PosType getMag(Band band) const ;
  PosType getMagBulge() const { return current.mag_bulge; }
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
  virtual void setMag(Band band,PosType my_mag){
    current.mag_map[band] = my_mag;
  };
  /// magnitude in specific band
  virtual void setMagBulge(Band band,PosType my_mag){
    current.bulge_mag_map[band] = my_mag;
  }
  
	/// bulge half light radius in arcseconds
	PosType getReff() const { return current.Reff/arcsecTOradians; }
	/// disk scale height in arcseconds
	PosType getRdisk() const { return current.Rdisk/arcsecTOradians; }
	
  /// the bulge to total flux ratio
	PosType getBtoT() const { return pow(10,(-current.mag_bulge + current.mag)/2.5); }
  /// position angle in radians
	PosType getPA() const { return current.PA; }
  /// inclination in radians
	PosType getInclination() const { return current.inclination;}
  float getSEDtype() const {return sedtype;}
  void setSEDtype(float s){ sedtype = s;}

  /// change the working band
  virtual void changeBand(Band band);
	
	/** Returns minimum of the radii at which disk and bulge have a surf. brightness equal to a fraction f of the central one
	* TODO: Fabio: Needs to be tested and improved (Bulge is so steep in the center that output values are very small)
  */
	inline PosType getMinSize(PosType f) {return std::min(1.678*current.Reff*fabs(cos(current.inclination))*pow(-log (f)/7.67,4),current.Rdisk*(-log (f)/1.67));}

  static PosType *getx(SourceOverzier &sourceo){return sourceo.source_x.x;}

protected:
  
  SourceSersic spheroid;
  
  
  float sedtype = -1;
  // renormalize the disk and bulge to agree with current mag and mag_bulge
  void renormalize_current();
	void assignParams(InputParams& params);
	
	/// haloID
	unsigned long haloID;

  struct Params{

    /// bulge half light radius
    PosType Reff=0;
    /// disk scale height
    PosType Rdisk=0;
	
    //PosType BtoT;
    PosType PA=0;
    PosType inclination=0;
    
    double cosPA;
    double sinPA;
	
    PosType cxx=0,cyy=0;
    PosType sbDo=0;
    PosType mag=0;
    PosType mag_bulge=0;
    
    // colors
    std::map<Band,double> mag_map;
    std::map<Band,double> bulge_mag_map;
    
    void print(){
      /// bulge half light radius
      std::cout << "Reff :" << Reff/arcsecTOradians << " arcsec ";
      std::cout << "Rdisk :" << Rdisk/arcsecTOradians << " arcsec ";
      std::cout << "PA :" << PA << " ";
      std::cout << "inclination :" << inclination << " radians";
      std::cout << "sbDo :" << sbDo << " ";
      std::cout << "mag :" << mag << " ";
      std::cout << "mag_bulge :" << mag_bulge << " ";
      
      std::cout << "BtoT :" << pow(10,(-mag_bulge + mag)/2.5) << std::endl;
    }
  };
	
  Params current;
	// optional position variables
};

/** \brief Adds some extra features to the SourceOverzier source like spiral
 * arms, and randomizations.
 *
 */

class SourceOverzierPlus : public SourceOverzier
{
public:
  //SourceOverzierPlus();
  SourceOverzierPlus(
                                         PosType my_mag         /// total magnitude
                                         ,PosType my_mag_bulge  /// magnitude of bulge
                                         ,PosType my_Reff       /// effective radius of bulge
                                         ,PosType my_Rdisk      /// scale hight of disk
                                         ,PosType my_PA         /// position angle
                                         ,PosType inclination   /// inclination in radians
                                         ,unsigned long my_id
                                         ,PosType my_z
                                         ,const PosType *theta
                                         ,Utilities::RandomNumbers_NR &ran
                                         );

  ~SourceOverzierPlus();
  
  SourceOverzierPlus(const SourceOverzierPlus &p);
  SourceOverzierPlus & operator=(const SourceOverzierPlus &p);
  
  //*** meed to be able to change band
  //*** put modes and phases into surface brightness
  //*** possible put gaussian texture on disk

  PosType SurfaceBrightness(PosType *y);

  /// magnitude in specific band
  void setMag(Band band,PosType my_mag){
    current.mag_map[band] = my_mag;
    original.mag_map[band] = my_mag;
  }
  /// magnitude in specific band
  void setMagBulge(Band band,PosType my_mag){
    current.bulge_mag_map[band] = my_mag;
    original.bulge_mag_map[band] = my_mag;
  }
  
  /// position angle in radians
  PosType getPA() const { return PA; }
  
  int getNarms() const {return Narms;}
  PosType getArmAmplitude() const {return Ad;}
  PosType getArmAlpha() const {return arm_alpha;}
  PosType getSphIndex() const {return spheroid.getSersicIndex();}
  PosType getSphAxisRatio() const {return spheroid.getAxesRatio();}
  //PosType getSphPA() const {return spheroid.getPA();}
  
  void changeBand(Band band);
  static PosType* getx(SourceOverzierPlus &sourceo){return sourceo.source_x.x;}

  /// Reset the position of the source in radians
  virtual inline void setTheta(PosType *xx){
    source_x[0] = xx[0];
    source_x[1] = xx[1];
    spheroid.setTheta(xx);
  }
  virtual void setTheta(PosType my_x,PosType my_y){
    source_x[0] = my_x;
    source_x[1] = my_y;
    spheroid.setTheta(my_x,my_y);
  }
  virtual void setTheta(const Point_2d &p){
    source_x = p;
    spheroid.setTheta(p[0],p[1]);
  }
  void setBulgeAxisRatio(PosType q){
    spheroid.setAxesRatio(q);
  }
  /// Randomly change some of the internal paramters and angles of the source
  void randomize(Utilities::RandomNumbers_NR &ran);
private:
  int Narms;
  PosType Ad,mctalpha,arm_alpha;
  SourceSersic spheroid;
  std::vector<PosType> modes;
  PosType disk_phase;
  PosType cosPA,sinPA,cosi,PA;
  
  SourceOverzier::Params original;  // original parameters
};
#endif /* GALAXIES_OVERZIER_H_ */
