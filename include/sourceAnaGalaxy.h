/*
 * souceAnaLens.h
 *
 *  Created on: Aug 13, 2012
 *      Author: bmetcalf
 */

#ifndef SOURCE_ANA_H_
#define SOURCE_ANA_H_

#include "source.h"

/**
 * \brief Source that represents an analytic galaxy surface brightness model.  It encapsulates a
 * OverGalaxy which is a model from R.Overzier et al. 2012 with a bulge and a disk.
 *
 *<pre>
 * Input parameters are only needed if the third constructor is used so that an input catalog is read
 * Input Parameters:
 *
 *	input_galaxy_file       file with catalog of galaxies
 *	source_band             Band that these sources are to be observed in. Must be one of the following SDSS_U,SDSS_G,SDSS_R,SDSS_I,SDSS_Z,J,H,Ks,i1, or i2
 *	source_band
 *	source_mag_limit        magnitude limit
 *	source_sb_limit         Minimum surface brightness limit for truncating sources.  By default it is 30. // TODO Fabio: What units is this?
 *
 *</pre>
 */
class MultiSourceAnaGalaxy: public Source{
public:
	MultiSourceAnaGalaxy(double mag, double BtoT, double Reff, double Rh, double PA, double inclination,double my_z,double *my_theta);
	MultiSourceAnaGalaxy(OverGalaxy *my_galaxy);
	MultiSourceAnaGalaxy(InputParams& params);
	~MultiSourceAnaGalaxy();

	/// Surface brightness of current galaxy in coordinates not centered on current galaxy.
	double SurfaceBrightness(double *y){
		double x[2] = {y[0]-galaxies[index].theta[0] , y[1]-galaxies[index].theta[1]};
		double s = galaxies[index].SurfaceBrightness(x);
		if (s < pow(10,-0.4*(48.6+sb_limit))*pow(180*60*60/pi,2)) return 0.;
		return s;

	}
	/// Total flux coming from the current galaxy in erg/sec/Hz/cm^2
	double getTotalFlux(){return pow(10,-(48.6+galaxies[index].getMag())/2.5);}

	void printSource();
	// Add a pre-constructed galaxy to the source collection
	void AddAGalaxy(OverGalaxy *my_galaxy){galaxies.push_back(*my_galaxy);}

	/** Used to change the "current" source that is returned when the surface brightness is subsequently
	 * called.  It also returns a reference to the particular OverGalaxy source model.
	 */
	OverGalaxy& setIndex (std::size_t i){
		if(i < galaxies.size())
			index = i;
		return galaxies[index];
	}
	/** The indexing operator can be used to change the "current" source that is returned when the surface brightness is subsequently
	 * called.  It also returns a reference to the particular OverGalaxy source model.
	 */
	OverGalaxy& operator[] (std::size_t i){
		if(i < galaxies.size())
			return galaxies[i];
		return galaxies[index];
	}
	
	const OverGalaxy& operator[] (std::size_t i) const {
		if(i < galaxies.size())
			return galaxies[i];
		return galaxies[index];
	}

	/// Return redshift of current source.
	double getZ(){return galaxies[index].z;}
	//double getRadius(){return max(galaxies[index]->Reff,galaxies[index]->Rh);}
	double getRadius(){return galaxies[index].getRadius();}
	/// Set redshift of current source.  Only changes the redshift while leaving position fixed.
	void setZ(double my_z){	galaxies[index].z = my_z;}

	unsigned long getID(){return galaxies[index].haloID;}


	/// Return angular position of current source.
	double *getX(){return galaxies[index].theta;}
	/// Set angular position of current source.
	void setX(double my_theta[2]){galaxies[index].theta[0] = my_theta[0]; galaxies[index].theta[1] = my_theta[1];}
	void setX(double my_x,double my_y){galaxies[index].theta[0] = my_x; galaxies[index].theta[1] = my_y;}
	std::size_t getNumberOfGalaxies() const {return galaxies.size();}

	void multiplier(double z,double mag_cut,int Multiplicity,long *seed);


private:
	Band band;
	float mag_limit;
	std::size_t index;

	std::vector<OverGalaxy> galaxies;
	std::string input_gal_file;

	void readDataFile(std::string input_gal_file);
	void assignParams(InputParams& params);

};

#endif /* SOURCE_H_ */
