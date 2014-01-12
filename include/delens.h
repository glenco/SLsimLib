#ifndef SLSIMLIB_DELENS_H
#define SLSIMLIB_DELENS_H

#include "lens.h"
#include "image_processing.h"

/**
 * \brief Use a given lens to de-lens a given PixelMap.
 * 
 * This function shoots rays from the image plane to the source plane and uses
 * the given lensed PixelMap to recreate the PixelMap in the source plane.
 * 
 * The base number of rays per pixel in the image plane is given by the Nrays
 * argument. The actual number of rays shot per pixel is the base number
 * multiplied by the inverse magnification at the center of the pixel. This
 * will lead to an approximately uniform distribution of rays in the source
 * plane.
 * 
 * A higher-resolution version of the delensed image can be produced by setting
 * the res_factor argument to a value less than unity.
 * 
 * The minimum magnification that is considered can be specified. This prevents
 * extremely demagnified areas of the image from degrading the reconstruction.
 */
PixelMap delens(Lens& lens, const PixelMap& lensed, std::size_t Nrays = 10, double res_factor = 1, double mag_min = 1);

#endif
