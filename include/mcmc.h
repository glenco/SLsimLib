#ifndef MCMC_H_
#define MCMC_H_

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
	typedef typename Spec::source_type source_type;
	typedef typename Spec::lens_type lens_type;
	
	typedef typename Spec::source_parameters source_parameters;
	typedef typename Spec::lens_parameters lens_parameters;
	
	typedef typename Spec::source_generator source_generator;
	typedef typename Spec::lens_generator lens_generator;
	
	class generator
	{
	public:
		generator(source_generator& s, lens_generator& l)
		: sgen(s), lgen(l)
		{
		}
		
		template<typename Parameters>
		inline Parameters operator()(Parameters p)
		{
			p.source = sgen(p.source);
			p.lens = lgen(p.lens);
			return p;
		}
		
	private:
		source_generator& sgen;
		lens_generator& lgen;
	};
	
	class parameters
	{
	public:
		source_parameters source;
		lens_parameters lens;
		
		parameters(source_parameters& s, lens_parameters& l)
		: source(s), lens(l)
		{
		}
	};
};

#endif
