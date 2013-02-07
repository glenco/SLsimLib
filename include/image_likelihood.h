#ifndef _image_likelihood_declare_
#define _image_likelihood_declare_

class PixelMap;

double image_chi_square(
	const PixelMap& data,
	const PixelMap& model,
	double background_subtract = 0,
	double background_noise = 0,
	double norm = 1
);

double image_chi_square(
	const PixelMap& data,
	const PixelMap& model,
	const PixelMap& mask,
	double background_subtract = 0,
	double background_noise = 0,
	double norm = 1
);

#endif
