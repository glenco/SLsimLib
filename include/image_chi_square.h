#ifndef _image_chi_square_declare_
#define _image_chi_square_declare_

class PixelMap;

double chi_square(
	const PixelMap& data,
	const PixelMap& model,
	double offset = 0,
	double noise = 0,
	double norm = 1
);

double chi_square(
	const PixelMap& data,
	const PixelMap& model,
	const PixelMap& mask,
	double offset = 0,
	double noise = 0,
	double norm = 1
);

#endif
