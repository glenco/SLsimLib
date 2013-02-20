#include "slsimlib.h"

double chi_square(
	const PixelMap& data,
	const PixelMap& model,
	double offset,
	double noise,
	double norm,
	PixelMap mask
)
{
	using std::pow;
	
	const std::size_t pixels = data.getNpixels();
	const bool masked = mask.valid();
	
	if(pixels != model.getNpixels())
	{
		SLSIMLIB_DEBUG("chi_square: Size of data and model do not agree!");
		return -1;
	}
	
	if(masked && pixels != mask.getNpixels())
	{
		SLSIMLIB_DEBUG("chi_square: Size of data and mask do not agree!");
		return -1;
	}
	
	double chi = 0;
	std::size_t dof = masked ? 0 : pixels*pixels;
	
	if(masked)
	{
		for(std::size_t i = 0, n = pixels*pixels; i < n; ++i)
		{
			if(!mask[i])
				continue;
			
			chi += pow(norm * (data[i] - model[i] - offset), 2)/(norm * (data[i] + noise));
			
			++dof;
		}
	}
	else
	{
		for(std::size_t i = 0, n = pixels*pixels; i < n; ++i)
			chi += pow(norm * (data[i] - model[i] - offset), 2)/(norm * (data[i] + noise));
	}
	
	return chi/dof;
}
