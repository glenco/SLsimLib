#include "isop_render.h"

#include "../Numerics/isop.h"
#include "../Numerics/quadrature.h"
#include "../Source/functor.h"

namespace SLsimLib
{
	namespace Render
	{
		double isop_render(Source& source, const double nodx[], const double nody[], double a0, double b0, double a1, double b1)
		{
			using Numerics::quadrature;
			using Numerics::quadrature_result;
			using Numerics::Rules::ClenshawCurtis;
			using Numerics::make_isop_map;
			
			// TODO: figure out best integration rule to use
			typedef ClenshawCurtis<5> rule;
			
			// the precision and accuracy goals
			double pg = 10e-2;
			double ag = 10e-2;
			
			// do the quadrature
			quadrature_result result = quadrature<rule>(make_isop_map(SourceFunctor(source), nodx, nody), a0, b0, a1, b1, pg, ag);
			
			// return the result, error should be within the required limits
			return result.value;
		}
	}
}
