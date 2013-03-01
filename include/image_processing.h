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
 * \brief Takes images and pixelizes the flux into regular pixel grid.
 *
 * The routine constructs an internal tree structure of the pixel points for fast calculation of location
 * of each point in the images.
 */

class PixelMap
{
public:
	PixelMap(const PixelMap& pmap, double degrading_factor);
	PixelMap();
	PixelMap(const PixelMap& other);
	PixelMap(std::size_t Npixels, double range, const double* center);
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
	void AddImages(ImageInfo *imageinfo,int Nimages,bool constant_sb);
	void AddImages(ImageInfo *imageinfo,int Nimages,double sigma);
	void printASCII();
	void printASCIItoFile(std::string filename);
	void printFITS(std::string filename);
	void smooth(double sigma);

	void ApplyPSF(std::string psf_file, double oversample_n = 1);

	inline double getValue(std::size_t i) const { return map[i]; }
	inline double operator[](std::size_t i) const { return map[i]; };
	
	friend void swap(PixelMap&, PixelMap&);
	
private:
	std::size_t map_size;
	float* map;
	
	std::size_t Npixels;
	double resolution,range,center[2];
	double map_boundary_p1[2],map_boundary_p2[2];
	
	double LeafPixelArea(IndexType i,Branch * branch1);
	void PointsWithinLeaf(Branch * branch1, std::list <unsigned long> &neighborlist);
	bool inMapBox(Branch * branch1);
};

/**
 * \ingroup Image
 * \brief Mask for PixelMap.
 * 
 * This class represents a mask that can be applied to a PixelMap to select
 * only a subset of its pixels.
 */
class PixelMask
{
public:
	/**
	 * \brief Threshold types.
	 * 
	 * These values represent the different types of thresholds that can be
	 * used when constructing a PixelMask from a PixelMap.
	 */
	enum ThresholdType
	{
		Greater,
		GreaterOrEqual,
		Less,
		LessOrEqual
	};
	
	/**
	 * \brief Create an empty PixelMask.
	 * 
	 * This creates an invalid and empty PixelMask
	 */
	PixelMask();
	
	/**
	 * \brief Create an PixelMask for a given size.
	 * 
	 * Create PixelMask for a number of pixels, all unmasked. The created mask
	 * is thus empty and all pixels are visible.
	 */
	PixelMask(std::size_t map_size);
	
	/**
	 * \brief Create a PixelMask from a PixelMap.
	 * 
	 * Create a new PixelMask given the pixel values from a PixelMap. The value
	 * given as threshold determines when a pixel is considered unmaskes (ie.
	 * visible), by applying the method given in threshold_type.
	 * 
	 * By default, all non-zero pixels are considered to be unmasked.
	 * 
	 * \param base The base PixelMap to convert to a mask.
	 * \param threshold Value that determines whether a pixel is unmasked.
	 * \param type The method to compare a pixel and the threshold.
	 */
	PixelMask(const PixelMap& base, double threshold = 0, ThresholdType type = Greater);
	
	/**
	 * \brief Create a PixelMask from FITS file.
	 * 
	 * Create a new PixelMask from a FITS file using an intermediate PixelMap. 
	 * 
	 * \param base The base PixelMap to convert to a mask.
	 * \param threshold Value that determines whether a pixel is unmasked.
	 * \param type The method to compare a pixel and the threshold.
	 */
	PixelMask(std::string file, double threshold = 0, ThresholdType type = Greater);
	
	/**
	 * Assignment operator.
	 */
	PixelMask& operator=(PixelMask other);
	
	/**
	 * Access unmasked pixel indices.
	 */
	std::size_t operator[](std::size_t i) const;
	
	/**
	 * Check if mask is valid.
	 */
	bool valid() const;
	
	/**
	 * Check if mask is empty.
	 * 
	 * An empty mask means that all pixels in a PixelMap are visible.
	 */
	bool empty() const;
	
	/**
	 * Get the size of the mask.
	 * 
	 * \return The number of unmasked pixels.
	 */
	std::size_t size() const;
	
	/**
	 * Get the size of the base PixelMap.
	 * 
	 * \return The total number of pixels.
	 */
	std::size_t base_size() const;
	
	friend void swap(PixelMask&, PixelMask&);
	
private:
	std::size_t map_size, mask_size;
	std::vector<std::size_t> pixels;
};

void pixelize(double *map,long Npixels,double range,double *center
		,ImageInfo *imageinfo,int Nimages,bool constant_sb,bool cleanmap
		,bool write_for_skymaker = false, std::string filename="");
void _SplitFluxIntoPixels(TreeHndl ptree,Branch *leaf,double *leaf_sb);
void smoothmap(double *map_out,double *map_in,long Npixels,double range,double sigma);

#endif
