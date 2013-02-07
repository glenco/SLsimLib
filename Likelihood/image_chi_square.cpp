#include "image_likelihood.h"

#include <cassert>
#include <cmath>

#include "image_processing.h"

#ifndef NDEBUG
#include <iostream>
#endif

inline double pixel_chi_square(
	unsigned long i,
	const PixelMap& data,
	const PixelMap& model,
	double bg,
	double noise,
	double norm
)
{
	double data_pixel = data.getValue(i);
	double model_pixel = model.getValue(i);
	
	double diff = norm * (data_pixel - model_pixel - bg);
	
	double sigma = std::sqrt(norm * (data_pixel + noise));
	
	return diff*diff/sigma/sigma;
}

double image_chi_square(
	const PixelMap& data,
	const PixelMap& model,
	double background_subtract,
	double background_noise,
	double norm
)
{
	const unsigned long data_Npixels = data.getNpixels();
	
	if(data_Npixels != model.getNpixels())
	{
#ifndef NDEBUG
		std::cerr << "Size of data and model do not agree!" << std::endl;
#endif
		return 0;
	}
	
	double chi = 0;
	
	unsigned int n = data_Npixels*data_Npixels;
	
	for(unsigned long i = 0; i < n; ++i)
		chi += pixel_chi_square(i, data, model, background_subtract, background_noise, norm);
	
	return chi/n;
}

double image_chi_square(
	const PixelMap& data,
	const PixelMap& model,
	const PixelMap& mask,
	double background_subtract,
	double background_noise,
	double norm
)
{
	const unsigned long data_Npixels = data.getNpixels();

	if(data_Npixels != model.getNpixels())
	{
#ifndef NDEBUG
		std::cerr << "Size of data and model do not agree!" << std::endl;
#endif
		return 0;
	}
	
	if(data_Npixels != mask.getNpixels())
	{
#ifndef NDEBUG
		std::cerr << "Size of data and mask do not agree!" << std::endl;
#endif
		return 0;
	}
	
	double chi = 0;
	unsigned long dof = 0;
	
	for(unsigned long i = 0, n = data_Npixels*data_Npixels; i < n; ++i)
	{
		if(!mask.getValue(i))
			continue;
		
		chi += pixel_chi_square(i, data, model, background_subtract, background_noise, norm);
		
		++dof;
	}
	
	return chi/dof;
}
