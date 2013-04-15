#ifndef IMAGE_LIKELIHOOD_H_
#define IMAGE_LIKELIHOOD_H_

#include "image_processing.h"
#include "model.h"

#include <cstddef>

/**
 * \brief Calculate the likelihood for a model based on image data.
 * 
 * This class generates a PixelMap by simulating lensing for a source and lens
 * model, and compares the result to the given PixelData. It uses the data's
 * chi_square function to calculate a likelihood.
 * 
 * To render the image PixelMap from source and lens, an internal grid is
 * created with an initial number of points (default: 64) and a center and
 * range which are taken from the PixelData. The initial number of grid points 
 * can be adjusted with the `gridPoints()` method.
 * 
 * Each source can have a maximum number of lensed images (default: 100).
 * This maximum number can be changed using the `imagesSize()` method.
 */
class ImageLikelihood
{
public:
	/**
	 * Indicate that this is a log-likelihood.
	 */
	static const bool logarithmic = true;
	
	/**
	 * Construct ImageLikelihood for the given data.
	 * 
	 * \param data The PixelData to be compared to simulated images.
	 */
	ImageLikelihood(PixelData data);
	
	/**
	 * \brief Destructor.
	 * 
	 * Destruct the ImageLikelihood. The Model destructor will also destroy
	 * the source and lens.
	 */
	~ImageLikelihood();
	
	/** Get the number of free parameters. */
	unsigned long dimension() const;
	
	/** Set the number of free parameters. */
	void dimension(unsigned long n);
	
	/** Get the data PixelMap. */
	PixelData data() const;
	
	/** Set the data PixelMap. */
	void data(PixelData data);
	
	/** Get the maximum number of lensed images. */
	std::size_t imagesSize() const;
	
	/** Set the maximum number of lensed images. */
	void imagesSize(std::size_t size);
	
	/** Get the number of initial grid points. */
	std::size_t gridPoints() const;
	
	/** Set the number of initial grid points. */
	void gridPoints(std::size_t n);
	
	/**
	 * Calculate the likelihood for the provided model.
	 */
	double operator()(Model& model);
	
private:
	inline void redof();
	
	PixelData dta;
	
	unsigned long dof;
	unsigned long dim;
	
	std::size_t images_size;
	ImageInfo* images;
	
	std::size_t grid_points;
	double grid_center[2];
	double grid_range;
};

#endif
