/*
 * source.h
 *
 *  Created on: Feb 6, 2012
 *      Author: mpetkova
 */

#ifndef SOURCE_H_
#define SOURCE_H_

#include "standard.h"
#include "InputParams.h"
#include "image_processing.h"

/** \brief Base class for all sources.
 *
 */
class Source
{
public:
	Source();
	virtual ~Source();
	
	// in lens.cpp
	/// Surface brightness of source in grid coordinates not source centered coordinates.
	virtual double SurfaceBrightness(double *y) = 0;
	virtual double getTotalFlux() = 0;
	virtual void printSource() = 0;

	/// Gets sb_limit in erg/cm^2/sec/rad^2/Hz
	double getSBlimit(){return sb_limit;}
	/// Gets sb_limit in mag/arcsec^2
	double getSBlimit_magarcsec(){return -2.5*log10(sb_limit*hplanck/pow((180*60*60/pi),2))-48.6;}
	
	// accessor functions that will sometimes be over ridden in class derivatives
	/// Redshift of source
	virtual inline double getZ(){return zsource;}
	virtual void setZ(double my_z){zsource = my_z;}
	/// Radius of source in radians
	virtual inline double getRadius(){return source_r;}
  /// Reset the radius of the source in radians
	virtual void setRadius(double my_radius){source_r = my_radius;}
	/// position of source in radians
	virtual inline double* getX(){return source_x;}
  /// Reset the position of the source in radians
	virtual inline void setX(double *xx){source_x[0] = xx[0]; source_x[1] = xx[1];}
	void setX(double my_x,double my_y){source_x[0] = my_x; source_x[1] = my_y;}
	/// In the case of a single plane lens, the ratio of angular size distances
	virtual inline double getDlDs(){return DlDs;}
	//TODO: BEN I think this need only be in the BLR source models
	virtual void setDlDs(double my_DlDs){DlDs = my_DlDs;}

	/// Sets sb_limit in erg/cm^2/sec/rad^2/Hz
	void setSBlimit(float limit) {sb_limit = limit;}
	/// Sets sb_limit in mag/arcsec^2
	void setSBlimit_magarcsec(float limit) {sb_limit = pow(10,-0.4*(48.6+limit))*pow(180*60*60/pi,2)/hplanck;}

	double changeFilter(std::string filter_in, std::string filter_out, std::string sed);
	double integrateFilter(std::vector<double> wavel_fil, std::vector<double> fil);
	double integrateFilterSed(std::vector<double> wavel_fil, std::vector<double> fil, std::vector<double> wavel_sed, std::vector<double> sed);

protected:
	virtual void assignParams(InputParams& params) = 0;
	
	// source parameters
	/// total source size, ie no flux outside this radius
	double source_r;
	/// center of source
	double source_x[2];
	
	/// redshift of source
	double zsource;
	/// Dl / Ds -- needed for the blr source models
	//TODO: Could this be moved into the BLR classes because they are the only ones that use it.
	double DlDs;
	double sb_limit;
  
};

typedef Source *SourceHndl;

/** \brief Class for sources described by an array of pixels
 *
 *  The sources are created as a square array of pixels of dimension Npixels x Npixels and pixel size equal to resolution.
 *
 */
class SourcePixelled: public Source{
public:
	SourcePixelled(double my_z, double* center, int Npixels, double resolution, double* arr_val);
	SourcePixelled(const PixelMap& gal_map, double z, double factor = 1.);
	SourcePixelled(InputParams& params);
	~SourcePixelled();
	double SurfaceBrightness(double *y);
	void printSource();
	inline double getTotalFlux(){return flux;}
	inline double getRadius(){return source_r;}
	inline double* getEll(){return ell;};
	inline double getQuad(int i, int j){return quad[i][j];};
	inline double getSize(){return size;};
	inline double* getCentroid(){return centroid;};
	inline double getMag(){return -2.5*log10(flux)-48.6;};
private:
	void assignParams(InputParams& params);
	void calcEll();
	void calcSize();
	void calcCentroid();
	void calcTotalFlux();
	double resolution;
	double range;
	long Npixels;
	double flux;
	double quad[2][2];
	double ell[2];
	double size;
	double centroid[2];
	std::valarray<double> values;
};

/** \brief Class for sources described by shapelets.
 *
 *  The sources are created from magnitude, scale radius, and the coefficients of their decomposition into the shapelets basis functions (Refregier et al., 2001).
 *  The coefficients can be read from a fits square array.
 *
 */
class SourceShapelets: public Source{
public:
	SourceShapelets(double my_z, double* my_center, double my_mag, double my_scale, std::valarray<double> my_coeff, double my_ang = 0.);
	SourceShapelets(double my_z, double* my_center, double my_mag, std::string shap_file, double my_ang = 0.);
	SourceShapelets(double* my_center, std::string shap_file, double my_ang = 0.);
	double SurfaceBrightness(double *y);
	void printSource();
	inline double getTotalFlux(){return flux;}
	inline double getRadius(){return source_r*10.;}
	inline double getMag(){return mag;}
	inline void setMag(double my_mag){mag = my_mag; NormalizeFlux();}

private:
	void assignParams(InputParams& params);
	double Hermite(int n, double x);
	void NormalizeFlux();
	std::valarray<double> coeff;
	int n1,n2;
	double flux, mag;
	double ang;
};

/// A uniform surface brightness circular source.
class SourceUniform : public Source{
public:
	SourceUniform(InputParams& params);
	~SourceUniform();

	double SurfaceBrightness(double *y);
	void assignParams(InputParams& params);
	void printSource();
	double getTotalFlux(){return pi*source_r*source_r;}
};

/// A source with a Gaussian surface brightness profile
class SourceGaussian : public Source{
public:
	SourceGaussian(InputParams& params);
	~SourceGaussian();
	
	/// internal scale parameter
	double source_gauss_r2;
	
	double SurfaceBrightness(double *y);
	void assignParams(InputParams& params);
	void printSource();
	double getTotalFlux(){std::cout << "No total flux in SourceGaussian yet" << std::endl; exit(1);}
};

/// Base class for all sources representing the Broad Line Region (BLR) of a AGN/QSO
class SourceBLR : public Source{
public:
	SourceBLR(InputParams& params);
	~SourceBLR();
	
	void printSource();
	double getTotalFlux(){std::cout << "No total flux in SourceBLR yet" << std::endl; exit(1);}
	
	virtual inline double getRadius(){return source_r_out;}
	
	/// lag time
	double source_tau;
	/// frequency
	double source_nu;
	float source_nuo;
	/// inner radius of BLR
	float source_r_in;
	/// outer radius of BLR
	float source_r_out;
	///inclination of BLR in radians, face on is
	float source_inclination;
	float source_opening_angle;
	float source_gamma;
	float source_BHmass;
	/// fraction of Keplerian velocity in random motions
	float source_fK;
	/// set to true to integrate over frequency
	bool source_monocrome;
private:
	void assignParams(InputParams& params);
};

/// A source representing a BLR with a Keplarian disk
class SourceBLRDisk : public SourceBLR{
public:
	double SurfaceBrightness(double *y);
	double getTotalFlux(){std::cout << "No total flux in SourceBLRDisk yet" << std::endl; exit(1);}
	
	SourceBLRDisk(InputParams&);
	~SourceBLRDisk();
};

/// A source representing a BLR with a spherical symmetry and circular orbits
class SourceBLRSph1 : public SourceBLR{
public:
	double SurfaceBrightness(double *y);
	double getTotalFlux(){std::cout << "No total flux in SourceBLRSph1 yet" << std::endl; exit(1);}
	
	SourceBLRSph1(InputParams&);
	~SourceBLRSph1();
};

/// A source representing a BLR with a spherical symmetry and random velocity dispersion
class SourceBLRSph2 : public SourceBLR{
public:
	double SurfaceBrightness(double *y);
	double getTotalFlux(){std::cout << "No total flux in SourceBLRSph2 yet" << std::endl; exit(1);}

	SourceBLRSph2(InputParams&);
	~SourceBLRSph2();
};

/// pointer to surface brightness function
//double (Source::*SurfaceBrightness)(double *y);

#endif /* SOURCE_H_ */
