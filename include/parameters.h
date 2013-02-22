#ifndef PARAMETERS_H_
#define PARAMETERS_H_

template<typename Source, typename Lens>
struct SourceLensParameters
{
	typedef Source source_type;
	typedef typename Source::parameters source_parameters;
	
	typedef Lens lens_type;
	typedef typename Lens::parameters lens_parameters;
	
	source_parameters source;
	lens_parameters lens;
	
	SourceLensParameters()
	{
	}
	
	SourceLensParameters(const Source& s, const Lens& l)
	{
		get(s, source);
		get(l, lens);
	}
};

#endif
