#ifndef MCMC_H_
#define MCMC_H_

#include "source.h"
#include "lens.h"
#include "model.h"
#include "parameters.h"

#include <vector>
#include <stdexcept>
#include <cstddef>
#include <cmath>
#include <ctime>

// TODO: put into C++ code
#ifndef MCMC_MAX_ITER
#define MCMC_MAX_ITER 1000
#endif

class MCMC
{
public:
	MCMC(double step, long seed = 0)
	: step(step), seed(seed)
	{
		if(!this->seed)
			this->seed = time(0);
	}
	
	template<typename Model, typename Likelihood>
	std::vector<Parameters> run(Model& model, Likelihood& L, const std::size_t n)
	{
		// the Markov chain
		std::vector<Parameters> chain;
		
		// the number of steps in the chain
		chain.reserve(n+1);
		
		// initial parameters
		chain.push_back(Parameters());
		
		// load parameters from model
		model.getParameters(chain.back());
		
		// current likelihood
		double Lx = L(model);
		
		// candidate likelihood
		double Ly;
		
		// run
		for(std::size_t i = 0; i < n; ++i)
		{
			// generate points until one is accepted
			for(std::size_t j = 0; true; ++j)
			{
				// maximum number of iterations
				if(j == MCMC_MAX_ITER)
					throw std::runtime_error("MCMC: maximum iterations reached.");
				
				// load parameters into model
				model.setParameters(chain.back());
				
				// reset consumed parameters
				chain.back().reset();
				
				// randomize model
				model.randomize(step, &seed);
				
				// calculate candidate likelihood
				Ly = L(model);
				
				// acceptance
				double a = Ly - Lx;
				
				// check if density increased
				if(a > 0)
					break;
				
				// step with probability a
				if(std::exp(a) > ran2(&seed))
					break;
			}
			
			// add parameters to chain
			chain.push_back(Parameters());
			
			// get candidate parameters into chain
			model.setParameters(chain.back());
			
			// keep likelihood
			Lx = Ly;
		}
		
		// chain is complete
		return chain;
	}
	
private:
	double step;
	long seed;
};

#endif
