#ifndef MCMC_H_
#define MCMC_H_

#include "source.h"
#include "lens.h"
#include "model.h"
#include "parameters.h"
#include "image_processing.h"

#include <vector>
#include <stdexcept>
#include <string>
#include <sstream>
#include <cstddef>
#include <cmath>
#include <ctime>

// TODO: put into C++ code
#ifndef MCMC_MAX_ITER
#define MCMC_MAX_ITER 1000
#endif

class IterationError : public std::runtime_error
{
public:
	IterationError(std::size_t max_iter) : val(max_iter), std::runtime_error("Maximum number of " + itos(max_iter) + " iterations exceeded.") {}
	
	std::size_t value() const { return val; }
	
private:
	static std::string itos(std::size_t i) { std::stringstream sstr; sstr << i; return sstr.str(); }
	std::size_t val;
};

/// Blocked Gibbs sampling MCMC algorithm.
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
		
		// create a grid
		double center[2];
		center[0] = data.getCenter()[0];
		center[1] = data.getCenter()[1];
		double range = data.getRange();
		Grid grid(model.lens, GRID_POINTS, center, range);
		
		// buffer for images
		ImageInfo images[MAX_N_IMAGES];
		int image_count;
		
		// create a pixmap for the results
		PixelMap pixmap(data.getCenter(), data.getNpixels(), data.getResolution());
		
		// the likelihood for each source
		std::vector<double> Lx(num_sources);
		for(std::size_t s = 0; s < num_sources; ++s)
		{
			// make source current
			model.source->setIndex(s);
			
			// update grid
			grid.RefreshSurfaceBrightnesses(model.source);
			
			// render model
			map_images(model.lens, model.source, &grid, &image_count, images, MAX_N_IMAGES, model.source->getRadius(), 0.1*model.source->getRadius(), 0, EachImage, true, false, true);
			
			// create pixmap from images
			pixmap.Clean();
			pixmap.AddImages(images, image_count);
			
			// calculate chi^2 from data
			Lx[s] = -0.5*data.chi_square(pixmap);
		}
		
		// keep track of iterations per step
		std::size_t iter;
		
		// run n complete cycles through all sources
		for(std::size_t i = 0; i < n; ++i)
		{
			// reset iterations
			iter = 0;
			
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
				
				// candidate likelihood
				double Ly;
				
				// TODO: treat overlap of sources
				
				// generate points until one is accepted
				while(++iter)
				{
					// maximum number of iterations
					if(iter == MCMC_MAX_ITER)
						throw IterationError(MCMC_MAX_ITER);
					
					// randomize source
					source->randomize(step, &seed);
					
					// update grid
					grid.RefreshSurfaceBrightnesses(source);
					
					// render model
					map_images(model.lens, source, &grid, &image_count, images, MAX_N_IMAGES, source->getRadius(), 0.1*source->getRadius(), 0, EachImage, true, false, true);
					
					// create pixmap from images
					pixmap.Clean();
					pixmap.AddImages(images, image_count);
					
					// calculate candidate chi^2
					Ly = -0.5*data.chi_square(pixmap);
					
					// check if density increased
					if(Ly > Lx[s])
						break;
					
					// step with probability a
					if(std::exp(Ly - Lx[s]) > ran2(&seed))
						break;
					
					// restore parameters
					source->setParameters(p);
					p.reset();
				}
				
				// candidate likelihood becomes current
				Lx[s] = Ly;
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
	double step;
	long seed;
};

#endif
