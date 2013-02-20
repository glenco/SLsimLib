#ifndef _image_chi_square_declare_
#define _image_chi_square_declare_

class PixelMap;

/**
 * \brief Calculate the chi^2 between model PixelMap and data PixelMap.
 * 
 * This function takes a PixelMap for model and data each and computes a
 * chi^2 value for them.
 * 
 * An additional mask can be given to only consider parts of the model and
 * data PixelMap. All pixels which have a zero value in the mask will not be
 * considered.
 * 
 * The PixelMap for data, model and (if supplied) mask must each have the same 
 * number of pixels.
 * 
 * The degree of freedom is the number of (unmasked) pixels. A normalization
 * factor can be given to transform the pixel value from [0, 1] to a count
 * number, which is then used to calculate sigma assuming Poisson errors.
 * 
 * \param data The (observed) data.
 * \param model The (calculated) model to compare to data.
 * \param offset A fixed background value to substract from the data.
 * \param noise A fixed amount of noise to be added to sigma.
 * \param norm The normalization factor to calculate sigma.
 * \param mask A mask to specify the pixels for the calculation.
 * 
 * \return The chi^2 for model and data.
 */
double chi_square(
	const PixelMap& data,
	const PixelMap& model,
	double offset = 0,
	double noise = 0,
	double norm = 1,
	PixelMap mask = PixelMap()
);

#endif
