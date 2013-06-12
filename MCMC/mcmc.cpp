#include "../include/mcmc.h"

#include "../include/source.h"
#include "../include/lens.h"
#include "../include/model.h"
#include "../include/image_processing.h"
#include "../include/grid_maintenance.h"
#include "../include/map_images.h"

#include <cmath>
#include <ctime>

namespace
{
	template<typename It>
	static double likelihood(It first, It last, const PixelData& data)
	{
		const PixelMap& image = data.image();
		const PixelMap& noise = data.noise();
		const PixelMask& mask = data.mask();
		double norm = data.getNormalization();
		
		double chi2 = 0;
		
		if(mask.valid())
		{
			for(std::size_t i = 0, n = mask.size(); i < n; ++i)
			{
				std::size_t k = mask[i];
				double model = 0;
				for(It j = first; j != last; ++j)
					model += (*j)[k];
				
				double diff = norm*(model - image[k]);
				double sigma2 = norm*model + (norm*noise[k])*(norm*noise[k]);
				
				chi2 += (diff*diff)/sigma2 + std::log(sigma2);
			}
		}
		else
		{
			for(std::size_t i = 0, n = image.size(); i < n; ++i)
			{
				double model = 0;
				for(It j = first; j != last; ++j)
					model += (*j)[i];
				
				double diff = norm*(model - image[i]);
				double sigma2 = norm*model + (norm*noise[i])*(norm*noise[i]);
				
				chi2 += (diff*diff)/sigma2 + std::log(sigma2);
			}
		}
		
		return -0.5*chi2;
	}
}

MCMC::MCMC(Lens& l, Sky& s, COSMOLOGY& c, const PixelData& d)
: lens(l), sky(s), cosmo(c), data(d)
{
}

std::vector<Parameters> MCMC::run(std::size_t n, double step, long* seed)
{
	std::size_t num_sources = sky.getNsources();
	
	// the Markov chain
	std::vector<Parameters> chain;
	
	// the number of steps in the chain
	chain.reserve(n+1);
	
	// put initial parameters in chain
	chain.push_back(Parameters());
	//lens.getParameters(chain.back());
	sky.getParameters(chain.back());
	
	// the image properties
	double center[2];
	center[0] = data.getCenter()[0];
	center[1] = data.getCenter()[1];
	double range = data.getRange();
	double resolution = data.getResolution();
	std::size_t npixels = data.getNpixels();
	
	// buffer for images
	ImageInfo images[MAX_N_IMAGES];
	int image_count;
	
	// create pixmaps for the individual sources
	std::vector<PixelMap> pixmaps(num_sources);
	
	// the likelihood for each source
	for(std::size_t s = 0; s < num_sources; ++s)
	{
		// get current source
		Source& source = sky.getSource(s);
		
		// check if source plane has changed
		if(lens.getSourceZ() != source.getZ())
			lens.ResetSourcePlane(&cosmo, source.getZ(), false);
		
		// create a grid
		Grid grid(&lens, GRID_POINTS, center, range);
		
		// render model
		map_images(&lens, &source, &grid, &image_count, images, MAX_N_IMAGES, source.getRadius(), 0.1*source.getRadius(), 0, EachImage, true, false, true);
		
		// create pixmap for images
		pixmaps[s] = PixelMap(center, npixels, resolution);
		pixmaps[s].AddImages(images, image_count);
	}
	
	// calculate current likelihood
	double Lx = likelihood(pixmaps.begin(), pixmaps.end(), data);
	
	// run n complete cycles through all sources
	for(std::size_t i = 0; i < n; ++i)
	{
		// randomize each source in turn
		for(std::size_t s = 0; s < num_sources; ++s)
		{
			// get current source
			Source& source = sky.getSource(s);
			
			// check if source plane has changed
			if(lens.getSourceZ() != source.getZ())
				lens.ResetSourcePlane(&cosmo, source.getZ(), false);
			
			// get source parameters
			Parameters p;
			source.getParameters(p);
			
			// randomize source
			source.randomize(step, seed);
			
			// create a grid
			Grid grid(&lens, GRID_POINTS, center, range);
			
			// render model
			map_images(&lens, &source, &grid, &image_count, images, MAX_N_IMAGES, source.getRadius(), 0.1*source.getRadius(), 0, EachImage, true, false, true);
			
			// create pixmap for images
			PixelMap pixmap(center, npixels, resolution);
			pixmap.AddImages(images, image_count);
			
			// swap pixmap into array to calculate likelihood
			swap(pixmaps[s], pixmap);
			
			// calculate candidate likelihood
			double Ly = likelihood(pixmaps.begin(), pixmaps.end(), data);
			
			// probability of accepting a candidate point
			if(Ly < Lx && std::exp(Ly - Lx) < ran2(seed))
			{
				// not accepted, restore parameters
				source.setParameters(p);
				
				// restore pixmap
				swap(pixmaps[s], pixmap);
			}
			else
			{
				// accepted, candidate likelihood becomes current
				Lx = Ly;
			}
		}
		
		// TODO: randomize lens
		
		// add parameters to chain after complete cycle
		chain.push_back(Parameters());
		//lens.getParameters(chain.back());
		sky.getParameters(chain.back());
	}
	
	// chain is complete
	return chain;
}
