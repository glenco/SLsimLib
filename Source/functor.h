#pragma once

#include "../include/source.h"
#include <functional>

namespace SLsimLib
{
	/// Functor to turn sources into binary functions
	struct SourceFunctor : public std::binary_function<double, double, double>
	{
		SourceFunctor(Source& source) : source(source) {}
		
		double operator()(double x, double y)
		{
			// TODO: make const double[2] as soon as possible
			double z[2] = {x, y};
			return source.SurfaceBrightness(z);
		}
		
		Source& source;
	};
}
