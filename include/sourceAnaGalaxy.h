/*
 * souceAnaLens.h
 *
 *  Created on: Aug 13, 2012
 *      Author: bmetcalf
 */

#ifndef SOURCE_ANA_H_
#define SOURCE_ANA_H_

#include "source.h"
#include "overzier_source.h"

/**
 * \brief Source that represents an analytic galaxy surface brightness model.  It encapsulates a
 * OverzierSource which is a model from R.Overzier et al. 2012 with a bulge and a disk.
 *
 *<pre>
 * Input parameters are only needed if the third constructor is used so that an input catalog is read
 * Input Parameters:
 *
 *	source_input_galaxy_file       file with catalog of galaxies
 *	source_band             Band that these sources are to be observed in. Must be one of the following SDSS_U,SDSS_G,SDSS_R,SDSS_I,SDSS_Z,J,H,Ks,i1, or i2
 *	source_band
 *	source_mag_limit        magnitude limit
 *	source_sb_limit         Minimum surface brightness limit for truncating sources.  By default it is 30. mag/sq. arcsec
 *
 *</pre>
 */
class SourceMultiAnaGalaxy: public Source{
public:
	SourceMultiAnaGalaxy(PosType mag, PosType BtoT, PosType Reff, PosType Rh, PosType PA, PosType inclination,PosType my_z,PosType *my_theta);
	SourceMultiAnaGalaxy(SourceOverzier *my_galaxy);
	SourceMultiAnaGalaxy(InputParams& params);
	~SourceMultiAnaGalaxy();
	
	/// Surface brightness of current galaxy.
	PosType SurfaceBrightness(PosType* x) {
		PosType sb = galaxies[index].SurfaceBrightness(x);
		if (sb < sb_limit) return 0.;
		return sb;
  }
	
	/// Total flux coming from the current galaxy in erg/sec/Hz/cm^2
	PosType getTotalFlux(){return pow(10,-(48.6+galaxies[index].getMag())/2.5);}

	void printSource();
	// Add a pre-constructed galaxy to the source collection
	void AddAGalaxy(SourceOverzier *my_galaxy){galaxies.push_back(*my_galaxy);}

	/** Used to change the "current" source that is returned when the surface brightness is subsequently
	 * called.  It also returns a reference to the particular OverzierSource source model.
	 */
	SourceOverzier& setIndex (std::size_t i){
		if(i < galaxies.size())
			index = i;
		return galaxies[index];
	}
	/** The indexing operator can be used to change the "current" source that is returned when the surface brightness is subsequently
	 * called.  It also returns a reference to the particular OverzierSource source model.
	 */
	SourceOverzier& operator[] (std::size_t i){
		if(i < galaxies.size())
			return galaxies[i];
		return galaxies[index];
	}
	
	const SourceOverzier& operator[] (std::size_t i) const {
		if(i < galaxies.size())
			return galaxies[i];
		return galaxies[index];
	}

	SourceOverzier& CurrentGalaxy(){
		return galaxies[index];
	}

	/// Return redshift of current source.
	PosType getZ(){return galaxies[index].getZ();}
	//PosType getRadius(){return max(galaxies[index]->Reff,galaxies[index]->Rh);}
	PosType getRadius(){return galaxies[index].getRadius();}
	/// Set redshift of current source.  Only changes the redshift while leaving position fixed.
	void setZ(PosType my_z){	galaxies[index].setZ(my_z);}

	unsigned long getID(){return galaxies[index].getID();}


	/// Return angular position of current source.
	PosType *getX(){return galaxies[index].getX();}
	/// Set angular position of current source.
	void setX(PosType my_theta[2]){galaxies[index].setX(my_theta);}
	void setX(PosType my_x,PosType my_y){galaxies[index].setX(my_x, my_y);}
	std::size_t getNumberOfGalaxies() const {return galaxies.size();}

	void multiplier(PosType z,PosType mag_cut,int Multiplicity,long *seed);
  void sortInRedshift();
  void sortInMag();
  /// returns field-of-view in deg^2 assuming region is square
  PosType getFOV(){return (rangex[1]-rangex[0])*(rangey[1]-rangey[0])*180*180/pi/pi;}
  
private:
	Band band;
	float mag_limit;
	std::size_t index;

	std::vector<SourceOverzier> galaxies;
	std::string input_gal_file;

	void readDataFile();
	void assignParams(InputParams& params);

  PosType rangex[2],rangey[2];
};

bool redshiftcompare(SourceOverzier s1,SourceOverzier s2);
bool magcompare(SourceOverzier s1,SourceOverzier s2);

/**
 * \brief Class for reading in and handling an array of SourceShapelets, made on the model of SourceMultiAnaGalaxy
 * Galaxies are read from files in a predifined directory and put into a std::vector that allows 
 * sorting in redshift and magnitude. When the index is set, an individual object can
 * be get through operator [].
 *
 *</pre>
 */
class SourceMultiShapelets: public Source{
public:
    SourceMultiShapelets();
    SourceMultiShapelets(InputParams params){};
	~SourceMultiShapelets();
    void sortInRedshift();
    void sortInMag();

    /// Surface brightness of current galaxy.
	PosType SurfaceBrightness(PosType* x) {
		PosType sb = galaxies[index].SurfaceBrightness(x);
		if (sb < sb_limit) return 0.;
		return sb;
    }
    
	void printSource();
	std::size_t getNumberOfGalaxies() const {return galaxies.size();}

    /// Total flux coming from the current galaxy in erg/sec/Hz/cm^2
	PosType getTotalFlux(){return pow(10,-(48.6+galaxies[index].getMag())/2.5);}

    /// Return angular position of current source.
	PosType *getX(){return galaxies[index].getX();}
	/// Set angular position of current source.
	void setX(PosType my_theta[2]){galaxies[index].setX(my_theta);}
	void setX(PosType my_x,PosType my_y){galaxies[index].setX(my_x, my_y);}

    /// Return redshift of current source.
	PosType getZ(){return galaxies[index].getZ();}
	//PosType getRadius(){return max(galaxies[index]->Reff,galaxies[index]->Rh);}
	PosType getRadius(){return galaxies[index].getRadius();}

	/** Used to change the "current" source that is returned when the surface brightness is subsequently
	 * called.  It also returns a reference to the particular OverzierSource source model.
	 */
	SourceShapelets& setIndex (std::size_t i){
		if(i < galaxies.size())
			index = i;
		return galaxies[index];
	}
    /** The indexing operator can be used to change the "current" source that is returned when the surface brightness is subsequently
	 * called.  It also returns a reference to the particular OverzierSource source model.
	 */
	SourceShapelets& operator[] (std::size_t i){
		if(i < galaxies.size())
			return galaxies[i];
		return galaxies[index];
	}
	
	const SourceShapelets& operator[] (std::size_t i) const {
		if(i < galaxies.size())
			return galaxies[i];
		return galaxies[index];
	}
    
	SourceShapelets& CurrentGalaxy(){
		return galaxies[index];
	}
    

private:
	void assignParams(InputParams& params);
 	std::size_t index;
	float mag_limit;
   
    void readCatalog();
 	std::vector<SourceShapelets> galaxies;
    
};

bool redshiftcompare_shap(SourceShapelets s1,SourceShapelets s2);
bool magcompare_shap(SourceShapelets s1,SourceShapelets s2);

#endif /* SOURCE_H_ */
