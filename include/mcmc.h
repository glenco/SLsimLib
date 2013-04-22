#ifndef MCMC_H_
#define MCMC_H_

#include <monaco/likelihood.hpp>
#include <monaco/metropolis.hpp>
#include <cstddef>

// macros to define the spec struct

#define MCMC_SOURCE(s) typedef s source_type
#define MCMC_LENS(l) typedef l lens_type

#define MCMC_SOURCE_PARAMETERS(sp) typedef sp source_parameters
#define MCMC_LENS_PARAMETERS(lp) typedef lp lens_parameters

#define MCMC_SOURCE_GENERATOR(sg) typedef sg source_generator
#define MCMC_LENS_GENERATOR(lg) typedef lg lens_generator

#define MCMC_LIKELIHOOD(l) typedef l likelihood

template<typename Spec>
class MCMC
{
public:
	typedef typename Spec::lens_type lens_type;
	typedef typename Spec::source_type source_type;
	
	typedef typename Spec::lens_parameters lens_parameters;
	typedef typename Spec::source_parameters source_parameters;
	
	typedef typename Spec::lens_generator lens_generator;
	typedef typename Spec::source_generator source_generator;
	
	typedef typename Spec::likelihood likelihood;
	
	class parameters
	{
	public:
		lens_parameters lens;
		source_parameters source;
		
		parameters()
		{
		}
		
		parameters(const lens_type& l, const source_type& s)
		{
			lens.from_lens(l);
			source.from_source(s);
		}
		
		parameters(const lens_parameters& l, const source_parameters& s)
		: lens(l), source(s)
		{
		}
	};
	
public:
	MCMC(Model<lens_type, source_type>& model, lens_generator& lens_gen, source_generator& source_gen, likelihood& lh)
	: model(model), gen(lens_gen, source_gen), lh(model, lh)
	{
	}
	
	template<typename Iterator>
	void run(Iterator outp, const std::size_t n)
	{
		// load current parameters into chain
		*outp = parameters(*model.lens, *model.source);
		
		// run metropolis algorithm
		monaco::metropolis(gen, lh, outp, n);
	}
	
private: /* classes */
	class mcmc_generator
	{
	public:
		mcmc_generator(lens_generator& l, source_generator& s)
		: lgen(l), sgen(s)
		{
		}
		
		inline parameters operator()(parameters p)
		{
			p.lens = lgen(p.lens);
			p.source = sgen(p.source);
			return p;
		}
		
	private:
		lens_generator& lgen;
		source_generator& sgen;
	};
	
	class mcmc_likelihood
	{
	public:
		static const bool logarithmic = monaco::likelihood_traits<likelihood>::logarithmic;
		
		mcmc_likelihood(Model<lens_type, source_type>& model, likelihood& lh)
		: model(model), lh(lh)
		{
		}
		
		inline double operator()(const parameters& p)
		{
			p.lens.to_lens(*model.lens);
			p.source.to_source(*model.source);
			
			return lh(model);
		}
		
	private:
		Model<lens_type, source_type>& model;
		likelihood& lh;
	};
	
private: /* members */
	Model<lens_type, source_type> model;
	mcmc_generator gen;
	mcmc_likelihood lh;
};

#endif
