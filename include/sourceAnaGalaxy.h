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
#include "simpleTreeVec.h"
#include "utilities.h"

/**
 * \brief Source that represents an analytic galaxy surface brightness model.  It encapsulates a
 * OverzierSource which is a model from R.Overzier et al. 2012 with a bulge and a disk.
 *
 *<pre>
 * Input parameters are only needed if the third constructor is used so that an input catalog is read
 * Input Parameters:
 *
 *	source_input_galaxy_file       file with catalog of galaxies
 *	source_band             Band that these sources are to be observed in. Must be one of the following SDSS_U,SDSS_G,SDSS_R,SDSS_I,SDSS_Z,J,H,Ks,IRAC1,IRAC2
 *	source_band
 *	source_mag_limit        magnitude limit
 *	source_sb_limit         Minimum surface brightness limit for truncating sources.  By default it is 30. mag/sq. arcsec
 *
 *</pre>
 */
class SourceMultiAnaGalaxy: public Source{
public:
	SourceMultiAnaGalaxy(PosType mag, PosType mag_bulge, PosType Reff, PosType Rdisk, PosType PA, PosType inclination,PosType my_z,PosType *my_theta,PosType zero_point,Utilities::RandomNumbers_NR &ran);
	SourceMultiAnaGalaxy(SourceOverzierPlus *my_galaxy);
	//SourceMultiAnaGalaxy(InputParams& params,Utilities::RandomNumbers_NR &ran);
	~SourceMultiAnaGalaxy();
	
	/// Surface brightness of current galaxy.
	PosType SurfaceBrightness(PosType* x) {
		PosType sb = galaxies[index].SurfaceBrightness(x);
		if (sb < sb_limit) return 0.;
		return sb;
  }
	
	/// Total flux coming from the current galaxy in erg/sec/Hz/cm^2
	PosType getTotalFlux() const {return mag_to_flux(galaxies[index].getMag());}

	void printSource();
	// Add a pre-constructed galaxy to the source collection
	void AddAGalaxy(SourceOverzierPlus *my_galaxy){galaxies.push_back(*my_galaxy);}

	/** \brief Used to change the "current" source that is returned when the surface brightness is subsequently
	 * called.  It also returns a reference to the particular OverzierSource source model.
	 */
	SourceOverzierPlus& setIndex (std::size_t i){
		if(i < galaxies.size())
			index = i;
		return galaxies[index];
	}
	/** \brief The indexing operator can be used to change the "current" source that is returned when the surface brightness is subsequently
	 * called.  It also returns a reference to the particular OverzierSource source model.
	 */
	SourceOverzierPlus& operator[] (std::size_t i){
		if(i < galaxies.size())
			return galaxies[i];
		return galaxies[index];
	}
	
	const SourceOverzierPlus& operator[] (std::size_t i) const {
		if(i < galaxies.size())
			return galaxies[i];
		return galaxies[index];
	}

	SourceOverzierPlus& CurrentGalaxy(){
		return galaxies[index];
	}

	/// Return redshift of current source.
	//PosType getZ() const {return galaxies[index].getZ();}
  	PosType getZ() const {return galaxies[index].getZ();}
  //PosType getRadius() const {return max(galaxies[index]->Reff,galaxies[index]->Rdisk);}
	PosType getRadius() const {return galaxies[index].getRadius();}
	/// Set redshift of current source.  Only changes the redshift while leaving position fixed.
	void setZ(PosType my_z){	galaxies[index].setZ(my_z);}
  void resetBand(Band my_band){
    for(size_t i=0;i<galaxies.size();++i) galaxies[i].changeBand(my_band);
    band = my_band;
  }

	unsigned long getID(){return galaxies[index].getID();}


	/// Return angular position of current source.
	Point_2d getTheta(){return galaxies[index].getTheta();}
	/// Set angular position of current source.
	void setTheta(PosType my_theta[2]){galaxies[index].setTheta(my_theta);}
  void setTheta(PosType my_x,PosType my_y){galaxies[index].setTheta(my_x, my_y);}
  void setTheta(const Point_2d &p){galaxies[index].setTheta(p);}

	std::size_t getNumberOfGalaxies() const {return galaxies.size();}

	void multiplier(PosType z,PosType mag_cut,int Multiplicity,Utilities::RandomNumbers_NR &ran);
  void sortInRedshift();
  void sortInMag(Band tmp_band);
  void sortInID();
  /// returns field-of-view in deg^2 assuming region is square
  PosType getFOV(){return (rangex[1]-rangex[0])*(rangey[1]-rangey[0])*180*180/PI/PI;}
  
  /** \brief Finds the closest source to the position theta[] on the sky in Cartesian distance.
   *   Returns the index of that source and its distance from theta[].
   */
  std::size_t findclosestonsky(PosType theta[],PosType *radius){
    size_t index;
    searchtree->NearestNeighbor(theta,*radius,index);
    return index;
  }

  /** \brief finds closest sources to theta on umlensed sky.  "radius" should have a size equal ti the number wanted
   */
  void findclosestonsky(PosType theta[],std::vector<PosType> &radius,std::vector<size_t> &indexes){
    indexes.resize(radius.size());
    searchtree->NearestNeighbors(theta,radius.size(),radius,indexes);
    return;
  }
 
  /** \brief finds objects within radios of theta[] on umlensed sky.    */
  void findonsky(PosType theta[],float radius,std::list<size_t> &indexes){
    indexes.clear();
    searchtree->PointsWithinCircle(theta,radius,indexes);
    return;
  }
  
  /** \brief finds objects within radios of theta[] on unlensed sky and within a redshift range.    */
  void findnear(PosType theta[],float radius,std::list<size_t> &indexes,PosType z_range[]){
    indexes.clear();
    searchtree->PointsWithinCircle(theta,radius,indexes);
    
    PosType z;
    // make redshift cut for lens sources
    for(std::list<size_t>::iterator itt = indexes.begin();
        itt != indexes.end(); ++itt){
      z = galaxies[*itt].getZ();
      if(z < z_range[0] || z > z_range[1]){
        indexes.erase(itt);
        if(itt != indexes.begin()) --itt;
      }
    }
    
    return;
  }

private:
  // make it uncopyable
  //SourceMultiAnaGalaxy(SourceMultiAnaGalaxy &s){};
  SourceMultiAnaGalaxy & operator=(SourceMultiAnaGalaxy &s){return s;}
  
	Band band;
	float mag_limit;
	std::size_t index;

	std::vector<SourceOverzierPlus> galaxies;
  TreeSimpleVec<SourceOverzierPlus> *searchtree;
	std::string input_gal_file;

	void readDataFileMillenn(Utilities::RandomNumbers_NR &ran);
	void assignParams(InputParams& params);

  PosType rangex[2],rangey[2];
};

bool redshiftcompare(SourceOverzierPlus s1,SourceOverzierPlus s2);
bool magcompare(SourceOverzierPlus s1,SourceOverzierPlus s2);
bool idcompare(SourceOverzierPlus s1,SourceOverzierPlus s2);

/**
 * \brief Class for reading in and handling an array of SourceShapelets, made on the model of SourceMultiAnaGalaxy
 * Galaxies are read from files in a predifined directory and put into a std::vector that allows 
 * sorting in redshift and magnitude. An individual object can be get via CurrentGalaxy() or the overloaded operator [].
 *
 */
class SourceMultiShapelets: public Source{
public:
  
  //SourceMultiShapelets(InputParams& params);
  
  SourceMultiShapelets(double mag_zero_point):Source(0,Point_2d(0,0),0,-1,mag_zero_point){};
  
  /// Reads in sources from a catalog.
  SourceMultiShapelets(const std::string &my_shapelets_folder  /// directory where shapelet files are located
                       ,Band my_band  /// band that will be used as default
                       ,double my_max_mag_limit  /// magnitude limit in that band
                       ,double my_min_mag_limit  /// magnitude limit in that band
                       ,double my_z_max          /// maximum redshift
                       ,double my_sb_limit       /// surface brightness limit
                       ,double maximum_radius   /// maximum radius (as defined in shapelet expansion) in radians
                       ,double zero_point        /// magnitude zreo point
                       );

  void input(const std::string &my_shapelets_folder  /// directory where shapelet files are located
                       ,Band my_band  /// band that will be used as default
                       ,double my_max_mag_limit  /// magnitude limit in that band
                       ,double my_min_mag_limit  /// magnitude limit in that band
                       ,double my_z_max          /// maximum redshift
                       ,double my_sb_limit     /// surface brightness limit
                       ,double maximum_radius   /// maximum radius (as defined in shapelet expansion) in radians
                       ,double zero_point        /// magnitude zreo point
                     );

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
  /// number of galaxies
  std::size_t size() const {return galaxies.size();}

    /// Total flux coming from the current galaxy in erg/sec/Hz/cm^2
	PosType getTotalFlux() const {return mag_to_flux(galaxies[index].getMag());}

    /// Return angular position of current source.
	Point_2d getTheta(){return galaxies[index].getTheta();}
	/// Set angular position of current source.
	void setTheta(PosType my_theta[2]){galaxies[index].setTheta(my_theta);}
	void setTheta(PosType my_x,PosType my_y){galaxies[index].setTheta(my_x, my_y);}
  void setTheta(const Point_2d &p){galaxies[index].setTheta(p);}

    /// Return redshift of current source.
	PosType getZ() const {return galaxies[index].getZ();}
	//PosType getRadius() const {return max(galaxies[index]->Reff,galaxies[index]->Rdisk);}
	PosType getRadius() const {return galaxies[index].getRadius();}

	/** Used to change the "current" source that is returned when the surface brightness is subsequently
	 * called.
	 */
	SourceShapelets& setIndex (std::size_t i){
		if(i < galaxies.size())
			index = i;
		return galaxies[index];
	}
  size_t getIndex() const {return index;}
  
    /** The indexing operator can be used to change the "current" source that is returned when the surface brightness is subsequently
	 * called.
	 */
  SourceShapelets& operator[] (std::size_t i){
    if(i < galaxies.size())
      return galaxies[i];
    return galaxies[index];
  }

  SourceShapelets& back(){
    return galaxies.back();
  }

	const SourceShapelets& operator[] (std::size_t i) const {
		if(i < galaxies.size())
			return galaxies[i];
		return galaxies[index];
	}
    
	SourceShapelets& CurrentGalaxy(){
		return galaxies[index];
	}
    
    /// Sets the active band for all the objects
	SourceShapelets& setBand (Band b){
        band = b;
        for (int i = 0; i < galaxies.size(); i++)
            galaxies[i].setActiveBand(band);
      return galaxies[index];
  }
  
  /// If the sources are already sorted by redshift this will find the index of the first galaxy with redshift above z
  int locateZ(PosType z) const {
      return Utilities::locate<SourceShapelets,PosType>(galaxies, z, [](PosType z,const SourceShapelets &s){return (z < s.getZ());});
  }
  
  int getCurrentID(){
    return galaxies[index].id;
  }
 
  std::vector<SourceShapelets> galaxies;

private:
	void assignParams(InputParams& params);
 	std::size_t index;
  float max_mag_limit;
  float min_mag_limit;
  float z_max;  // maximum redshift allowed
  Band band;
  double radius_max;
 
  void readCatalog();
	std::string shapelets_folder;
};

bool redshiftcompare_shap(SourceShapelets s1,SourceShapelets s2);
bool magcompare_shap(SourceShapelets s1,SourceShapelets s2);

#endif /* SOURCE_H_ */
