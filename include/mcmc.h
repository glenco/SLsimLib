#ifndef MCMC_H_
#define MCMC_H_

#include "lens.h"
#include "sky.h"
#include "raw_data.h"
#include "image_processing.h"

#include <vector>
#include <cstddef>

/// Blocked Metropolis MCMC algorithm.
class MCMC
{
public:
	// TODO: do this right
	static const int MAX_N_IMAGES = 100;
	static const int GRID_POINTS = 64;
	
	MCMC(Lens& lens, Sky& sky, const PixelData& data);
	
	std::vector<RawData> run(std::size_t n, double step, long* seed);
	
private:
	Lens& lens;
	Sky& sky;
	PixelData data;
};

#endif
