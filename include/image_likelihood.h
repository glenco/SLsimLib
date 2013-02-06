#ifndef _image_likelihood_declare_
#define _image_likelihood_declare_

#include <image_processing.h>

double chi_square (PixelMap data, PixelMap model, float background_subtract = 0., float background_noise = 0., float norm = 1., PixelMap* mask = NULL);

#endif
