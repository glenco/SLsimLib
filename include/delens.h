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
 * This is done by distributing a number of rays uniformly over the pixel area
 * and raytracing them from the image plane to the source plane. The j'th ray
 * contributes a fraction of luminosity y_ij to the i'th pixel, where the
 * luminosity x_i in the image plane is
 * 
 * \sum_j |\mu_ij| y_ij = x_i
 * 
 * and \mu_ij is the magnification at the image plane position of the ray.
 * 
 * From this, the luminosity in the source plane can be reconstructed as
 * 
 * y_ij = |\mu_ij|^{-1} x_i/n
 * 
 * where n is the number of rays shot for the pixel.
 * 
 * The base number of rays per pixel in the image plane is given by the Nrays
 * argument. The actual number of rays shot per pixel is the base number
 * multiplied by the inverse magnification at the center of the pixel. This
 * will lead to a more uniform distribution of rays in the source plane.
 * 
 * A higher-resolution version of the delensed image can be produced by setting
 * the res_factor argument to a value less than unity.
 */
PixelMap delens(Lens& lens, const PixelMap& lensed, std::size_t Nrays = 10, double res_factor = 1);

#endif
