/*
 * image_processing.h
 *
 *  Created on: Feb 28, 2010
 *      Author: R.B. Metcalf
 */

#ifndef IMAGE_PROCESSING_H_
#define IMAGE_PROCESSING_H_

#include "Tree.h"



/** \ingroup Image
 * \brief Takes image structure and pixelizes the flux into regular pixel grid which then
 * can be exported as a fits file, smoothed, etc. like an image.
 *
 */

class PixelMap
{
public:
	PixelMap(const PixelMap& pmap, double res_ratio);
	PixelMap();
	PixelMap(const PixelMap& other);
	PixelMap(const PixelMap& pmap, const double* center, std::size_t Npixels);
	PixelMap(const double* center, std::size_t Npixels, double resolution);
	PixelMap(std::string filename);
	~PixelMap();
	
	PixelMap& operator=(PixelMap other);
	
	inline bool valid() const { return map_size; };
	inline std::size_t size() const { return map_size; };
	
	inline std::size_t getNpixels() const { return Npixels; }
	inline double getRange() const { return range; }
	inline const double* getCenter() const { return center; }
	inline double getResolution() const { return resolution; }
	
	void Clean();

	void AddImages(ImageInfo *imageinfo,int Nimages,float rescale = 1.);
	void AddCurve(ImageInfo *curve,double value);
	void drawline(double x1[],double x2[],double value);

	void Renormalize(double factor);
	void AddValue(std::size_t i, double value);
	void AssignValue(std::size_t i, double value);
	void printASCII() const;
	void printASCIItoFile(std::string filename) const;
	void printFITS(std::string filename, bool verbose = false) const;

	void smooth(double sigma);

	inline double getValue(std::size_t i) const { return map[i]; }
	inline double operator[](std::size_t i) const { return map[i]; };
	
	friend void swap(PixelMap&, PixelMap&);
	friend bool agree(const PixelMap& a, const PixelMap& b);
	
private:
	std::size_t map_size;
	float* map;

	std::size_t Npixels;
	double resolution,range,center[2];
	double map_boundary_p1[2],map_boundary_p2[2];

	double LeafPixelArea(IndexType i,Branch * branch1);
	void PointsWithinLeaf(Branch * branch1, std::list <unsigned long> &neighborlist);
	bool inMapBox(Branch * branch1);
	bool inMapBox(double * branch1);
};


typedef enum {Euclid_VIS,Euclid_Y,Euclid_J,Euclid_H} Telescope;

/** \ingroup Image
 * \brief It creates a realistic image from the output of a ray-tracing simulation.
 *
 * It translates pixel values in observed units (counts/sec), applies PSF and noise.
 * Input must be in photons/(cm^2*Hz).
 */
class Observation
{
public:
	Observation(Telescope tel_name);
	Observation(float diameter, float transmission, float exp_time, int exp_num, float back_mag, float ron, float seeing = 0.);
	Observation(float diameter, float transmission, float exp_time, int exp_num, float back_mag, float ron, std::string psf_file, float oversample);
	float getExpTime(){return exp_time;}
	int getExpNum(){return exp_num;}
	float getBackMag(){return back_mag;}
	float getDiameter(){return diameter;}
	float getTransmission(){return transmission;}
	float getRon(){return ron;}
	float getSeeing(){return seeing;}
	float getZeropoint(){return mag_zeropoint;}
	std::valarray<float> getPSF(){return map_psf;}
	PixelMap Convert (PixelMap &map, bool psf, bool noise);
	PixelMap Convert_back (PixelMap &map);

private:
	float diameter;  // diameter of telescope (in cm)
	float transmission;  // total transmission of the instrument
	float mag_zeropoint;  // magnitude of a source that produces one count/sec in the image
	float exp_time;  // total exposure time (in sec)
	int exp_num;  // number of exposures
	float back_mag;  // sky (or background) magnitude in mag/arcsec^2
	float ron;  // read-out noise in electrons/pixel
	float seeing;  // full-width at half maximum of the gaussian smoothing
	std::valarray<float> map_psf;  // array of the point spread function
	float oversample; // psf oversampling factor
	double pix_size; // pixel size (in rad)
	bool telescope; // was the observation created from a default telescope?

	PixelMap AddNoise(PixelMap &pmap);
	PixelMap PhotonToCounts(PixelMap &pmap);
	PixelMap ApplyPSF(PixelMap &pmap);

};

void pixelize(double *map,long Npixels,double range,double *center
		,ImageInfo *imageinfo,int Nimages,bool constant_sb,bool cleanmap
		,bool write_for_skymaker = false, std::string filename="");
void _SplitFluxIntoPixels(TreeHndl ptree,Branch *leaf,double *leaf_sb);
void smoothmap(double *map_out,double *map_in,long Npixels,double range,double sigma);

#endif
