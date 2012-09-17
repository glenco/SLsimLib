/*
 * MOKAlens.h
 *
 *  Created on: Jun 8, 2012
 *      Author: mpetkova
 */


#ifndef MOKALENS_H_
#define MOKALENS_H_

#include <cosmo.h>
#include <utilities.h>
#include <source.h>
#include <grid_maintenance.h>
#include <profile.h>
#include <valarray>
#include <MOKAfits.h>

//TODO Improve this comment with more complete description of what a MOKAmap is used for.
/**
 * \brief the MOKA map structure, containing all quantities that define it
 *
 * Note: To use this class requires setting the ENABLE_FITS compiler flag and linking
 * the cfits library.

#include <MOKAfits.h>

 */
struct MOKAmap{
	/// values for the map
	std::valarray<float> convergence;
	std::valarray<float> alpha1;
	std::valarray<float> alpha2;
	std::valarray<float> gamma1;
	std::valarray<float> gamma2;
	std::valarray<float> gamma3;
	std::valarray<float> Signlambdar;
	std::valarray<float> Signlambdat;
	std:: vector<double> x;	 
	int nx,ny;
	double boxl,boxlMpc,zlens,zsource,omegam,omegal,h,DL;
	double inarcsec;
	double center[2];
};
//TODO Improve this comment.
/**
 *  \brief A class to represents the MOKA lens map
 *
 * Note: To use this class requires setting the ENABLE_FITS compiler flag and linking
 * the cfits library.
 */

//TODO Change to physical length units !!!!
class MOKALens : public Lens{
public:


	MOKALens(std::string);
	MOKALens(std::string paramfile,LensHalo *LH);

	~MOKALens();

	bool set;	/// the name of the MOKA input file
	std::string MOKA_input_file;
	/// if >=1 (true), do analyzis only; if = 0 (false) change units to internal GLAMER units and prepare for ray-shooting
	int flag_MOKA_analyze;

	void readParamfile(std::string);
	void rayshooterInternal(double *ray, double *alpha, float *gamma, float *kappa, bool kappa_off);
	void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off){ERROR_MESSAGE(); exit(1);};
	void setZlens(double zlens);
	double getZlens();
	void setInternalParams(CosmoHndl,SourceHndl);
	void saveImage(GridHndl grid, bool saveprofile=true);
	void saveImage(bool saveprofile=true);
	void saveKappaProfile();
	void saveGammaProfile();
	void saveProfiles(double &RE3);
	void initMap();

	MOKAmap *map;
	LensHalo *LH;

	void estSignLambdas();
	void EinsteinRadii(double &RE1, double &RE2);

};

#endif /* MOKALENS_H_ */

