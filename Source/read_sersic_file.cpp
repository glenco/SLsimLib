#include "../include/sky.h"
#include "../include/sersic_source.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <cmath>

const double PI = 3.141592653589793238462643383279502884;
const double LOG10 = std::log(10.);

namespace
{
	template<typename Number>
	inline std::string to_str(Number n)
	{
		std::stringstream sstr;
		sstr << n;
		return sstr.str();
	}
}

void Sky::readSersicFile(std::string filename, double z)
{
	std::ifstream input(filename.c_str());
	if(!input)
		throw std::runtime_error("SExtractor catalog: could not read " + filename + "");
	
	std::size_t line = 1;
	
	while(!input.eof())
	{
		if(input.peek() == std::istream::traits_type::eof())
			break;
		
		if(input.peek() == '#')
		{
			input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			++line;
			continue;
		}
		
		if(input.peek() == '\n')
		{
			input.get();
			++line;
			continue;
		}
		
		double mag;
		double Reff;
		double theta;
		double n;
		double q;
		double x[2];
		
		input >> mag >> Reff >> theta >> n >> q >> x[0] >> x[1];
		
		if(input.fail() || input.get() != '\n')
			throw std::runtime_error("Sersic catalog: line " + to_str(line) + ": syntax error");
		
		++line;
		
		Reff *= 3600;
		
		theta *= -PI/180.;
		
		n *= 0.90; // TODO: don't cheat
		
		x[0] *= -pi/180.;
		x[1] *= pi/180.;
		
		addSource(SourceSersic(mag, Reff, theta, n, q, z, x));
	}
}
