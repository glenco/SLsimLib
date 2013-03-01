#include "slsimlib.h"
#include "debug.h"

double chi_square(
	const PixelMap& data,
	const PixelMap& model,
	double offset,
	double noise,
	double norm,
	PixelMask mask
)
{
	using std::pow;
	
	if(data.getNpixels() != model.getNpixels())
	{
		SLSIMLIB_DEBUG("Size of data and model do not agree!");
		return -1;
	}
	
	if(mask.valid() && data.size() != mask.base_size())
	{
		SLSIMLIB_DEBUG("Size of data and mask do not agree!");
		return -1;
	}
	
	double chi2 = 0;
	
	if(mask.empty())
	{
		for(std::size_t i = 0, n = data.size(); i < n; ++i)
			chi2 += pow(norm * (data[i] - model[i] - offset), 2)/(norm * (data[i] + noise));
	}
	else
	{
		for(std::size_t i = 0, n = mask.size(); i < n; ++i)
			chi2 += pow(norm * (data[mask[i]] - model[mask[i]] - offset), 2)/(norm * (data[mask[i]] + noise));
	}
	
	return chi2;
}
