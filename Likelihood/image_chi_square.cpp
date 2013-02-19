#include "slsimlib.h"

#include <cmath>

#ifndef NDEBUG
#include <iostream>
#endif

double chi_square(
	const PixelMap& data,
	const PixelMap& model,
	double offset,
	double noise,
	double norm
)
{
	std::size_t pixels = data.getNpixels();
	
	if(model.getNpixels() != pixels)
	{
#ifndef NDEBUG
		std::cerr << "chi_square: Size of data and model do not agree!" << std::endl;
#endif
		return -1;
	}
	
	double chi = 0;
	std::size_t n = pixels*pixels;
	
	for(std::size_t i = 0; i < n; ++i)
		chi += std::pow(norm * (data[i] - model[i] - offset), 2)/(norm * (data[i] + noise));
	
	return chi/n;
}

double chi_square(
	const PixelMap& data,
	const PixelMap& model,
	const PixelMap& mask,
	double offset,
	double noise,
	double norm
)
{
	if(!mask.valid())
		return chi_square(data, model, offset, noise, norm);
	
	std::size_t pixels = data.getNpixels();
	
	if(pixels != model.getNpixels())
	{
#ifndef NDEBUG
		std::cerr << "chi_square: Size of data and model do not agree!" << std::endl;
#endif
		return -1;
	}
	
	if(pixels != mask.getNpixels())
	{
#ifndef NDEBUG
		std::cerr << "chi_square: Size of data and mask do not agree!" << std::endl;
#endif
		return -1;
	}
	
	double chi = 0;
	std::size_t dof = 0;
	
	for(std::size_t i = 0, n = pixels*pixels; i < n; ++i)
	{
		if(!mask[i])
			continue;
		
		chi += std::pow(norm * (data[i] - model[i] - offset), 2)/(norm * (data[i] + noise));
		
		++dof;
	}
	
	return chi/dof;
}
