#pragma once

#include <vector>
#include <algorithm>
#include <cmath>

namespace SLsimLib
{	
	namespace Numerics
	{
		namespace detail
		{
			/// Subdivided area in quadrature with value and error estimate.
			struct quadrature_task
			{
				quadrature_task(double value, double error,
								double a0, double b0, double a1, double b1,
								double x00, double x0c, double x0n,
								double xc0, double xcc, double xcn,
								double xn0, double xnc, double xnn)
				: value(value), error(error),
				  a0(a0), b0(b0), a1(a1), b1(b1),
				  x00(x00), x0c(x0c), x0n(x0n),
				  xc0(xc0), xcc(xcc), xcn(xcn),
				  xn0(xn0), xnc(xnc), xnn(xnn)
				{
				}
				
				bool operator<(const quadrature_task& other) const
				{
					return std::abs(error) < std::abs(other.error);
				}
				
				double value;
				double error;
				
				double a0;
				double b0;
				double a1;
				double b1;
				
				double x00, x0c, x0n;
				double xc0, xcc, xcn;
				double xn0, xnc, xnn;
			};
		}
		
		/// Quadrature rules. Must be of the closed type.
		namespace Rules
		{
			/// Clenshaw-Curtis quadrature data
			template<int N>
			struct ClenshawCurtis
			{
				static const int n = 2*N-1;
				static const double p[n];
				static const double w[n];
				static const double e[n];
			};
			
			/// Lobatto-Kronrod quadrature data
			template<int N>
			struct LobattoKronrod
			{
				static const int n = 2*N-1;
				static const double p[n];
				static const double w[n];
				static const double e[n];
			};
		}
		
		/// Result of numerical integration.
		struct quadrature_result
		{
			double value;
			double error;
		};
		
		template<typename Rule, typename Function>
		inline quadrature_result quadrature(Function f, double a0, double b0, double a1, double b1, double pg, double ag)
		{
			using detail::quadrature_task;
			
			// container type for queue
			typedef std::vector<quadrature_task> queue_container;
			
			// iterator types for queue
			typedef queue_container::const_iterator queue_iterator;
			
			// calculated function values
			double values[Rule::n][Rule::n];
			
			// constants for easier access
			const int n = Rule::n-1;
			const int c = (Rule::n-1)/2;
			
			// do initial quadrature
			for(int i = 0; i < Rule::n; ++i)
				for(int j = 0; j < Rule::n; ++j)
					values[i][j] = f(0.5*((b0-a0)*Rule::p[i]+a0+b0), 0.5*((b1-a1)*Rule::p[j]+a1+b1));
			
			// calculate initial result
			quadrature_result result = {0, 0};
			for(int i = 0; i < Rule::n; ++i)
			{
				for(int j = 0; j < Rule::n; ++j)
				{
					result.value += Rule::w[i]*Rule::w[j]*values[i][j];
					result.error += Rule::e[i]*Rule::e[j]*values[i][j];
				}
			}
			result.value = 0.25*(b0-a0)*(b1-a1)*result.value;
			result.error = 0.25*(b0-a0)*(b1-a1)*std::abs(result.error);
			
			// queue of subdivisions to refine
			queue_container queue;
			std::make_heap(queue.begin(), queue.end());
			
			// push initial task to queue
			queue.push_back(quadrature_task(
				result.value, result.error,
				a0, b0, a1, b1,
				values[0][0], values[0][c], values[0][n],
				values[c][0], values[c][c], values[c][n],
				values[n][0], values[n][c], values[n][n]
			));
			std::push_heap(queue.begin(), queue.end());
			
			// refine until accuracy or precision goal is reached
			while(result.error > std::abs(result.value)*pg && result.error > ag)
			{
				// get top task from queue
				quadrature_task task = queue.front();
				
				// pop task from queue
				std::pop_heap(queue.begin(), queue.end());
				queue.pop_back();
				
				// subdivide area
				for(int k = 0; k < 4; ++k)
				{
					// four quadrants for subdivision
					switch(k)
					{
						case 0:
							a0 = task.a0;
							b0 = 0.5*(task.a0+task.b0);
							a1 = task.a1;
							b1 = 0.5*(task.a1+task.b1);
							values[0][0] = task.x00;
							values[0][n] = task.x0c;
							values[n][0] = task.xc0;
							values[n][n] = task.xcc;
							break;
							
						case 1:
							a0 = 0.5*(task.a0+task.b0);
							b0 = task.b0;
							a1 = task.a1;
							b1 = 0.5*(task.a1+task.b1);
							values[0][0] = task.xc0;
							values[0][n] = task.xcc;
							values[n][0] = task.xn0;
							values[n][n] = task.xnc;
							break;
							
						case 2:
							a0 = task.a0;
							b0 = 0.5*(task.a0+task.b0);
							a1 = 0.5*(task.a1+task.b1);
							b1 = task.b1;
							values[0][0] = task.x0c;
							values[0][n] = task.x0n;
							values[n][0] = task.xcc;
							values[n][n] = task.xcn;
							break;
							
						case 3:
							a0 = 0.5*(task.a0+task.b0);
							b0 = task.b0;
							a1 = 0.5*(task.a1+task.b1);
							b1 = task.b1;
							values[0][0] = task.xcc;
							values[0][n] = task.xcn;
							values[n][0] = task.xnc;
							values[n][n] = task.xnn;
							break;
					}
					
					// do subdivision quadrature
					for(int j = 1; j < n; ++j)
						values[0][j] = f(0.5*((b0-a0)*Rule::p[0]+a0+b0), 0.5*((b1-a1)*Rule::p[j]+a1+b1));
					for(int i = 1; i < n; ++i)
						for(int j = 0; j < Rule::n; ++j)
							values[i][j] = f(0.5*((b0-a0)*Rule::p[i]+a0+b0), 0.5*((b1-a1)*Rule::p[j]+a1+b1));
					for(int j = 1; j < n; ++j)
						values[n][j] = f(0.5*((b0-a0)*Rule::p[n]+a0+b0), 0.5*((b1-a1)*Rule::p[j]+a1+b1));
					
					// calculate subdivision result
					quadrature_result sub = {0, 0};
					for(int i = 0; i < Rule::n; ++i)
					{
						for(int j = 0; j < Rule::n; ++j)
						{
							sub.value += Rule::w[i]*Rule::w[j]*values[i][j];
							sub.error += Rule::e[i]*Rule::e[j]*values[i][j];
						}
					}
					sub.value = 0.25*(b0-a0)*(b1-a1)*sub.value;
					sub.error = 0.25*(b0-a0)*(b1-a1)*sub.error;
					
					// push subdivision task to queue
					queue.push_back(quadrature_task(
						sub.value, sub.error,
						a0, b0, a1, b1,
						values[0][0], values[0][c], values[0][n],
						values[c][0], values[c][c], values[c][n],
						values[n][0], values[n][c], values[n][n]
					));
					std::push_heap(queue.begin(), queue.end());
				}
				
				// update result value and error
				result.value = 0;
				result.error = 0;
				for(queue_iterator it = queue.begin(); it != queue.end(); ++it)
				{
					result.value += it->value;
					result.error += it->error;
				}
				result.error = std::abs(result.error);
			}
			
			// done
			return result;
		}
	}
}
