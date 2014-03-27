#pragma once

#include "../include/source.h"

namespace SLsimLib
{
	namespace Render
	{
		/**
		 * \brief Integrate source flux using the isoparameterized lens.
		 */
		double isop_render(Source& source, const double nodx[], const double nody[], double a0, double b0, double a1, double b1);
	}
}
