#ifndef MCMC_H_
#define MCMC_H_

#include "source.h"
#include "lens.h"
#include "model.h"
#include "parameters.h"
#include "image_processing.h"

#include <vector>
#include <cstddef>
#include <cmath>
#include <ctime>

/// Blocked Metropolis MCMC algorithm.
class MCMC
{
public:
	// TODO: do this right
	static const int MAX_N_IMAGES = 100;
	static const int GRID_POINTS = 64;
	
	MCMC(double step, long seed = 0)
	: step(step), seed(seed)
	{
		if(!this->seed)
			this->seed = time(0);
	}
	
	template<typename Lens>
	std::vector<Parameters> run(Model<Lens, MultiSource>& model, const PixelData& data, const std::size_t n)
	{
		// number of sources
		std::size_t num_sources = model.source->size();
		
		// the Markov chain
		std::vector<Parameters> chain;
		
		// the number of steps in the chain
		chain.reserve(n+1);
		
		// put initial parameters in chain
		chain.push_back(Parameters());
		model.getParameters(chain.back());
		
		// the image properties
		double center[2];
		center[0] = data.getCenter()[0];
		center[1] = data.getCenter()[1];
		double range = data.getRange();
		double resolution = data.getResolution();
		std::size_t npixels = data.getNpixels();
		
		// create a grid
		Grid grid(model.lens, GRID_POINTS, center, range);
		
		// buffer for images
		ImageInfo images[MAX_N_IMAGES];
		int image_count;
		
		// create pixmaps for the individual sources
		std::vector<PixelMap> pixmaps(num_sources);
		
		// the likelihood for each source
		for(std::size_t s = 0; s < num_sources; ++s)
		{
			// make source current
			model.source->setIndex(s);
			
			// update grid
			grid.RefreshSurfaceBrightnesses(model.source);
			
			// render model
			map_images(model.lens, model.source, &grid, &image_count, images, MAX_N_IMAGES, model.source->getRadius(), 0.1*model.source->getRadius(), 0, EachImage, true, false, true);
			
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
				// set current source
				model.source->setIndex(s);
				
				// get current source
				Source* source = model.source->getCurrent();
				
				// get source parameters
				Parameters p;
				source->getParameters(p);
				
				// randomize source
				source->randomize(step, &seed);
				
				// update grid
				grid.RefreshSurfaceBrightnesses(source);
				
				// render model
				map_images(model.lens, source, &grid, &image_count, images, MAX_N_IMAGES, source->getRadius(), 0.1*source->getRadius(), 0, EachImage, true, false, true);
				
				// create pixmap for images
				PixelMap pixmap(center, npixels, resolution);
				pixmap.AddImages(images, image_count);
				
				// swap pixmap into array to calculate likelihood
				swap(pixmaps[s], pixmap);
				
				// calculate candidate likelihood
				double Ly = likelihood(pixmaps.begin(), pixmaps.end(), data);
				
				// probability of accepting a candidate point
				if(Ly < Lx && std::exp(Ly - Lx) < ran2(&seed))
				{
					// not accepted, restore parameters
					source->setParameters(p);
					
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
			model.getParameters(chain.back());
		}
		
		// chain is complete
		return chain;
	}
	
private:
	template<typename It>
	static double likelihood(It first, It last, const PixelData& data)
	{
		const PixelMap& image = data.image();
		const PixelMap& noise = data.noise();
		double norm = data.getNormalization();
		
		double chi2 = 0;
		
		for(std::size_t i = 0, n = image.size(); i < n; ++i)
		{
			double model = 0;
			for(It j = first; j != last; ++j)
				model += (*j)[i];
			
			double diff = norm*(model - image[i]);
			double sigma2 = norm*model + (norm*noise[i])*(norm*noise[i]);
			
			chi2 += (diff*diff)/sigma2 + std::log(sigma2);
		}
		
		return -0.5*chi2;
	}
	
	double step;
	long seed;
};

#endif
