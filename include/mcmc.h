#ifndef MCMC_H_
#define MCMC_H_

#include "parameters.h"
#include "image_likelihood.h"

#include <monaco/markov.hpp>
#include <monaco/metropolis.hpp>

#include <list>
#include <iterator>

template<
	typename Generator,
	typename Likelihood = ImageLikelihood<typename Generator::source_type, typename Generator::lens_type>
>
class MCMC
{
public:
	typedef typename Generator::source_type source_type;
	typedef typename Generator::lens_type lens_type;
	
	typedef SourceLensParameters<source_type, lens_type> parameter_type;
	
	typedef std::list<parameter_type> container_type;
	
	typedef typename container_type::iterator iterator;
	
	MCMC(Generator& generator, Likelihood& likelihood)
	: metrop(generator, likelihood)
	{
	}
	
	void add(const source_type& source, const lens_type& lens)
	{
		SourceLensParameters<source_type, lens_type> p;
		p.source << source;
		p.lens << lens;
		chain.push_back(p);
	}
	
	void add(const Parameters<source_type>& source, const Parameters<lens_type>& lens)
	{
		SourceLensParameters<source_type, lens_type> p;
		p.source = source;
		p.lens = lens;
		chain.push_back(p);
	}
	
	void add(const parameter_type& p)
	{
		chain.push_back(p);
	}
	
	void step()
	{
		chain.push_back(metrop(chain.back()));
	}
	
	template<typename Size>
	void step(Size n)
	{
		monaco::markov_chain(chain, n, metrop);
	}
	
	iterator begin()
	{
		return chain.begin();
	}
	
	iterator end()
	{
		return chain.end();
	}
	
	parameter_type front() const
	{
		return chain.front();
	}
	
	parameter_type back() const
	{
		return chain.back();
	}
	
	std::size_t size() const
	{
		return chain.size();
	}
	
	void clear()
	{
		chain.clear();
	}
	
	void clear(std::ptrdiff_t count)
	{
		iterator end = (count > 0) ? chain.begin() : chain.end();
		std::advance(end, count);
		chain.erase(chain.begin(), end);
	}
	
private:
	container_type chain;
	
	monaco::metropolis<Generator, Likelihood> metrop;
};

#endif
