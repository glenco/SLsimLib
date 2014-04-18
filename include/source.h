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
  /// shell constructor
	Source(){
    source_r = source_x[0] = source_x[1] = zsource = 0;
    setSBlimit_magarcsec(30.);
  };

	virtual ~Source();
	
	// in lens.cpp
	/// Surface brightness of source in grid coordinates not source centered coordinates.
	virtual PosType SurfaceBrightness(PosType *y) = 0;
	virtual PosType getTotalFlux() = 0;
	virtual void printSource() = 0;

	/// Gets sb_limit in erg/cm^2/sec/rad^2/Hz
	PosType getSBlimit(){return sb_limit;}
	/// Gets sb_limit in mag/arcsec^2
	PosType getSBlimit_magarcsec(){return -2.5*log10(sb_limit*hplanck/pow((180*60*60/pi),2))-48.6;}
	
	// accessor functions that will sometimes be over ridden in class derivatives
	/// Redshift of source
	virtual inline PosType getZ(){return zsource;}
	virtual void setZ(PosType my_z){zsource = my_z;}
	/// Radius of source in radians
	virtual inline PosType getRadius(){return source_r;}
  /// Reset the radius of the source in radians
	virtual void setRadius(PosType my_radius){source_r = my_radius;}
	/// position of source in radians
	virtual inline PosType* getX(){return source_x;}
  /// Reset the position of the source in radians
	virtual inline void setX(PosType *xx){source_x[0] = xx[0]; source_x[1] = xx[1];}
	void setX(PosType my_x,PosType my_y){source_x[0] = my_x; source_x[1] = my_y;}
	/// In the case of a single plane lens, the ratio of angular size distances
	virtual inline PosType getDlDs(){return DlDs;}
	//TODO: BEN I think this need only be in the BLR source models
	virtual void setDlDs(PosType my_DlDs){DlDs = my_DlDs;}

	/// Sets sb_limit in erg/cm^2/sec/rad^2/Hz
	void setSBlimit(float limit) {sb_limit = limit;}
	/// Sets sb_limit in mag/arcsec^2
	void setSBlimit_magarcsec(float limit) {sb_limit = pow(10,-0.4*(48.6+limit))*pow(180*60*60/pi,2)/hplanck;}

	PosType changeFilter(std::string filter_in, std::string filter_out, std::string sed);
	PosType integrateFilter(std::vector<PosType> wavel_fil, std::vector<PosType> fil);
	PosType integrateFilterSed(std::vector<PosType> wavel_fil, std::vector<PosType> fil, std::vector<PosType> wavel_sed, std::vector<PosType> sed);

protected:
	virtual void assignParams(InputParams& params) = 0;
	
	// source parameters
	/// total source size, ie no flux outside this radius
	PosType source_r;
	/// center of source
	PosType source_x[2];
	
	/// redshift of source
	PosType zsource;
	/// Dl / Ds -- needed for the blr source models
	//TODO: Could this be moved into the BLR classes because they are the only ones that use it.
	PosType DlDs;
	PosType sb_limit;
  
};

typedef Source *SourceHndl;

/** \brief Class for sources described by an array of pixels
 *
 *  The sources are created as a square array of pixels of dimension Npixels x Npixels and pixel size equal to resolution.
 *
 */
class SourcePixelled: public Source{
public:
	SourcePixelled(PosType my_z, PosType* center, int Npixels, PosType resolution, PosType* arr_val);
	SourcePixelled(const PixelMap& gal_map, PosType z, PosType factor = 1.);
	SourcePixelled(InputParams& params);
	~SourcePixelled();
	PosType SurfaceBrightness(PosType *y);
	void printSource();
	inline PosType getTotalFlux(){return flux;}
	inline PosType getRadius(){return source_r;}
	inline PosType* getEll(){return ell;};
	inline PosType getQuad(int i, int j){return quad[i][j];};
	inline PosType getSize(){return size;};
	inline PosType* getCentroid(){return centroid;};
	inline PosType getMag(){return -2.5*log10(flux)-48.6;};
private:
	void assignParams(InputParams& params);
	void calcEll();
	void calcSize();
	void calcCentroid();
	void calcTotalFlux();
	PosType resolution;
	PosType range;
	long Npixels;
	PosType flux;
	PosType quad[2][2];
	PosType ell[2];
	PosType size;
	PosType centroid[2];
	std::valarray<PosType> values;
};


/** \brief Class for sources described by shapelets.
 *
 *  The sources are created from magnitude, scale radius, and the coefficients of their decomposition into the shapelets basis functions (Refregier et al., 2001).
 *  The coefficients can be read from a fits square array.
 *  In case the magnitude is not given as input, the constructor will read an array of values from the shapelet file header. One can then choose the desired magnitude with setActiveMag.
 *
 */
class SourceShapelets: public Source{
public:
    //SourceShapelets();
	SourceShapelets(PosType my_z, PosType my_mag, PosType my_scale, std::valarray<PosType> my_coeff, PosType* my_center = 0, PosType my_ang = 0.);
	SourceShapelets(PosType my_z, PosType my_mag, std::string shap_file, PosType *my_center = 0, PosType my_ang = 0.);
	SourceShapelets(std::string shap_file, PosType* my_center = 0, PosType my_ang = 0.);
	PosType SurfaceBrightness(PosType *y);
	void printSource();
	inline PosType getTotalFlux(){return flux;}
	inline PosType getRadius(){return source_r*10.;}
	inline PosType getMag(){return mag;}
    inline PosType getID(){return id;}
    inline void setActiveBand(shap_band band){mag = mags[band]; flux = fluxes[band];}

private:
	void assignParams(InputParams& params);
	PosType Hermite(int n, PosType x);
	void NormalizeFlux();
	std::valarray<PosType> coeff;
	int n1,n2;
    int id;
	PosType flux, mag;
	PosType ang;
    PosType mags[10], fluxes[10];
    PosType coeff_flux;
};

/// A uniform surface brightness circular source.
class SourceUniform : public Source{
public:
	SourceUniform(InputParams& params);
	~SourceUniform();

	PosType SurfaceBrightness(PosType *y);
	void assignParams(InputParams& params);
	void printSource();
	PosType getTotalFlux(){return pi*source_r*source_r;}
};

/// A source with a Gaussian surface brightness profile
class SourceGaussian : public Source{
public:
	SourceGaussian(InputParams& params);
	~SourceGaussian();
	
	/// internal scale parameter
	PosType source_gauss_r2;
	
	PosType SurfaceBrightness(PosType *y);
	void assignParams(InputParams& params);
	void printSource();
	PosType getTotalFlux(){std::cout << "No total flux in SourceGaussian yet" << std::endl; exit(1);}
};

/// Base class for all sources representing the Broad Line Region (BLR) of a AGN/QSO
class SourceBLR : public Source{
public:
	SourceBLR(InputParams& params);
	~SourceBLR();
	
	void printSource();
	PosType getTotalFlux(){std::cout << "No total flux in SourceBLR yet" << std::endl; exit(1);}
	
	virtual inline PosType getRadius(){return source_r_out;}
	
	/// lag time
	PosType source_tau;
	/// frequency
	PosType source_nu;
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
	PosType SurfaceBrightness(PosType *y);
	PosType getTotalFlux(){std::cout << "No total flux in SourceBLRDisk yet" << std::endl; exit(1);}
	
	SourceBLRDisk(InputParams&);
	~SourceBLRDisk();
};

/// A source representing a BLR with a spherical symmetry and circular orbits
class SourceBLRSph1 : public SourceBLR{
public:
	PosType SurfaceBrightness(PosType *y);
	PosType getTotalFlux(){std::cout << "No total flux in SourceBLRSph1 yet" << std::endl; exit(1);}
	
	SourceBLRSph1(InputParams&);
	~SourceBLRSph1();
};

/// A source representing a BLR with a spherical symmetry and random velocity dispersion
class SourceBLRSph2 : public SourceBLR{
public:
	PosType SurfaceBrightness(PosType *y);
	PosType getTotalFlux(){std::cout << "No total flux in SourceBLRSph2 yet" << std::endl; exit(1);}

	SourceBLRSph2(InputParams&);
	~SourceBLRSph2();
};

/// pointer to surface brightness function
//PosType (Source::*SurfaceBrightness)(PosType *y);


/**  \brief Class to handle redshift-dependent quasar luminosity functions.
 *
 *   At the moment, only i band is available
 *   QLF from Ross et al. 2013
 *   k-correction from Richards et al. 2006
 */
class QuasarLF{
	public:
		QuasarLF(PosType red, PosType mag_limit, long *seed);
    ~QuasarLF();
		// returns the integral of the luminosity function at redshift red
		PosType getNorm() {return pow(10,log_phi)*norm;}; // in Mpc^(-3)
		PosType getRandomMag();
		PosType getRandomFlux();

	private:
		PosType kcorr;
		PosType red;
		PosType mag_limit;
		PosType mstar;
		PosType log_phi;
		PosType alpha;
		PosType beta;
		long *seed;
		PosType norm;
		PosType mag_max, mag_min;
		int arr_nbin;
		PosType* mag_arr;
		PosType* lf_arr;
		PosType dl;

		typedef PosType (QuasarLF::*pt2MemFunc)(PosType) const;
		PosType nintegrateQLF(pt2MemFunc func, PosType a,PosType b,PosType tols) const;
		PosType trapzQLFlocal(pt2MemFunc func, PosType a, PosType b, int n, PosType *s2) const;
		PosType lf_kernel (PosType mag) const;
};

#endif /* SOURCE_H_ */
