#pragma once

#include "source.h"

namespace SLsimLib
{
	namespace detail
	{
		/// Shape functions of the 8-node quadrilateral element.
		template<int N> double shape_fun(double xi, double eta);
		
		template<> double shape_fun<0>(double xi, double eta) { return -0.25*(1-xi)*(1-eta)*(1+xi+eta); }
		template<> double shape_fun<1>(double xi, double eta) { return 0.5*(1-xi)*(1+xi)*(1-eta); }
		template<> double shape_fun<2>(double xi, double eta) { return -0.25*(1+xi)*(1-eta)*(1-xi+eta); }
		template<> double shape_fun<3>(double xi, double eta) { return 0.5*(1+xi)*(1+eta)*(1-eta); }
		template<> double shape_fun<4>(double xi, double eta) { return -0.25*(1+xi)*(1+eta)*(1-xi-eta); }
		template<> double shape_fun<5>(double xi, double eta) { return 0.5*(1-xi)*(1+xi)*(1+eta); }
		template<> double shape_fun<6>(double xi, double eta) { return -0.25*(1-xi)*(1+eta)*(1+xi-eta); }
		template<> double shape_fun<7>(double xi, double eta) { return 0.5*(1-xi)*(1+eta)*(1-eta); }
		
		/// Static-recurse addition of the shape functions via templates.
		template<int N = 8>
		inline double isop(const double nodes[], double xi, double eta)
		{
			return isop<N-1>(nodes, xi, eta) + nodes[N-1]*shape_fun<N-1>(xi, eta);
		}
		
		template<>
		inline double isop<0>(const double[], double, double)
		{
			return 0;
		}
	}
	
	/**
	 * \brief Do the isoparametric interpolation given a list of 8 nodes.
	 * 
	 * This function takes 8 function values on the border of a [-1,1]x[-1,1]
	 * square region interpolates them over the square.
	 * 
	 * The first node is the lower left corner (-1,-1) and the nodes must
	 * continue in counter-clockwise order.
	 * 
	 * \param nodes The nodes to be interpolated.
	 * \param xi First coordinate (x) in isoparametric space.
	 * \param eta Second coordinate (y) in isoparametric space.
	 * 
	 * \return Interpolation of the nodes at (xi,eta).
	 */
	inline double isop(const double nodes[], double xi, double eta)
	{
		return detail::isop(nodes, xi, eta);
	}
	
	/// Functor to transparently map [-1,1]x[-1,1] to interpolation space.
	template<typename Function>
	struct isop_map
	{
		isop_map(Function f, const double* x, const double* y) : f(f), x(x), y(y) {}
		
		double operator()(double xi, double eta)
		{
			return f(isop(x, xi, eta), isop(y, xi, eta));
		}
		
		Function f;
		const double* x;
		const double* y;
	};
	
	/// Creates an isoparametric mapping functor.
	template<typename Function>
	inline isop_map<Function> make_isop_map(Function f, const double* x, const double* y)
	{
		return isop_map<Function>(f, x, y);
	}
	
	/**
	 * \brief Integrate source flux using the isoparameterized lens.
	 * 
	 * Integrate the surface brightness of a source using the isoparametric
	 * representation of the lens mapping.
	 * 
	 * The isoparametric representation is defined by the 8 nodes on the border
	 * of a [-1,1]x[-1,1] square. The first node (index 0) is in the bottom-left
	 * corner, and then nodes are counted counter-clockwise:
	 * 
	 * \verbatim
	 * 6 -- 5 -- 4
	 * |         |
	 * 7         3
	 * |         |
	 * 0 -- 1 -- 2
	 * \endverbatim
	 * 
	 * The positions of the nodes in the isoparametric space are
	 * 
	 * \verbatim
	 * (-1, 1) -- ( 0, 1) -- ( 1, 1)
	 *    |                     |
	 * (-1, 0)               ( 1, 0)
	 *    |                     |
	 * (-1,-1) -- ( 0,-1) -- ( 1,-1)
	 * \endverbatim
	 * 
	 * The interval to be integrated is given in terms of the isoparametric
	 * coordinates xi and eta.
	 */
	double isop_render(
		Source& source       /// Source to be integrated
		,const double nodx[] /// Source plane x coordinates of 8 nodes.
		,const double nody[] /// Source plane y coordinates of 8 nodes.
		,double a_xi         /// Lower bound of xi integration.
		,double b_xi         /// Upper bound of xi integration.
		,double a_eta        /// Lower bound of eta integration.
		,double b_eta        /// Upper bound of eta integration.
	);
}
