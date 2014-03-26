#pragma once

namespace SLsimLib
{
	namespace Numerics
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
		
		/// Do the isoparametric interpolation given a list of 8 nodes.
		inline double isop(const double nodes[], double xi, double eta)
		{
			return detail::isop(nodes, xi, eta);
		}
		
		/// Functor to transparently map [-1,1]x[-1,1] to interpolation space.
		template<typename Function>
		struct isop_map
		{
			isop_map(const Function& f, const double* x, const double* y) : f(f), x(x), y(y) {}
			
			double operator()(double xi, double eta) const
			{
				return f(isop(x, xi, eta), isop(y, xi, eta));
			}
			
			const Function& f;
			const double* x;
			const double* y;
		};
		
		/// Creates an isoparametric mapping functor.
		template<typename Function>
		isop_map<Function> make_isop_map(const Function& f, const double* x, const double* y)
		{
			return isop_map<Function>(f, x, y);
		}
	}
}
