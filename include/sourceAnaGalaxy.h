/*
 * souceAnaLens.h
 *
 *  Created on: Aug 13, 2012
 *      Author: bmetcalf
 */

#ifndef SOURCE_ANA_H_
#define SOURCE_ANA_H_

/**
 * \brief Source that represents an analytic galaxy surface brightness model.  It encapsulates a
 * OverGalaxy which is a model from R.Oversier et al. 2012 with a bulge and a disk.
 */
class MultiSourceAnaGalaxy: public Source{
public:
	MultiSourceAnaGalaxy(double mag, double BtoT, double Reff, double Rh, double PA, double inclination,double my_z,double *my_theta);
	MultiSourceAnaGalaxy(OverGalaxy *my_galaxy);
	MultiSourceAnaGalaxy(std::string input_gal_file,double my_mag_limit = 100);
	~MultiSourceAnaGalaxy();

	/// Surface brightness of current galaxy in coordinates not centered on current galaxy.
	double SurfaceBrightness(double *y){
		double x[2] = {y[0]-galaxies[index]->theta[0] , y[1]-galaxies[index]->theta[1]};
		return galaxies[index]->SurfaceBrightness(x);
	}

	void printSource();
	//TODO Does this need to be public?
	void readParamfile(std::string);
	//TODO Write the definition of this
	void AddAGalaxy(OverGalaxy *my_galaxy);

	/** Used to change the "current" source that is returned when the surface brightness is subsequently
	 * called.  It also returns a reference to the particular OverGalaxy source model.
	 */
	OverGalaxy& setIndex (const unsigned long my_index){
		if(my_index < 0 || my_index > galaxies.size()-1) return *galaxies[0];

		index = my_index;
		return *galaxies[my_index];
	}
	/** The indexing operator can be used to change the "current" source that is returned when the surface brightness is subsequently
	 * called.  It also returns a reference to the particular OverGalaxy source model.
	 */
	OverGalaxy& operator[] (const unsigned long my_index){
		if(my_index < 0 || my_index > galaxies.size()-1) return *galaxies[0];
		index = my_index;

		return *galaxies[my_index];
	}

	/// Return redshift of current source.
	double getZ(){return galaxies[index]->z;}

	/// Set redshift of current source.
	void setZ(double my_z){galaxies[index]->z = my_z;}

	/// Return angular position of current source.
	double *getX(){return galaxies[index]->theta;}
	/// Set angular position of current source.
	void setX(double my_theta[2]){galaxies[index]->theta[0] = my_theta[0]; galaxies[index]->theta[1] = my_theta[1];}
	unsigned long getNumberOfGalaxies(){return galaxies.size();}

private:
	unsigned long index;

	bool mem_allocated;
	std::vector<OverGalaxy*> galaxies;
	std::string input_gal_file;

	void readDataFile(std::string input_gal_file,double my_mag_limit = 100);
};

#endif /* SOURCE_H_ */
