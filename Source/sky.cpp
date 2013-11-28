#include "../include/sky.h"

#include <iterator>
#include <ctime>

Sky::Sky()
{
}

Sky::Sky(InputParams& params)
{
	std::string filename;
	
	// check if there is a Millenium data file
	if(params.get("input_galaxy_file", filename))
	{
		Band band;
		double mag_limit;
		
		if(!params.get("source_band", band))
		{
			std::cerr << "ERROR: Must assign source_band in parameter file " << params.filename() << std::endl;
			std::cerr << "Could be that specified band is not available " << std::endl;
			exit(1);
		}
		
		if(!params.get("source_mag_limit", mag_limit))
		{
			std::cerr << "ERROR: Must assign source_mag_limit in parameter file " << params.filename() << std::endl;
			exit(1);
		}
		
		// load galaxies from input file
		readGalaxyFile(filename, band, mag_limit);
	}
}

Sky::~Sky()
{
}

double Sky::getFov()
{
	if(srcs.empty())
		return 0;
	
	double rangex[2];
	double rangey[2];
	
	const double* x = srcs[0].getX();
	rangex[0] = x[0];
	rangex[1] = x[0];
	rangey[0] = x[1];
	rangey[1] = x[1];
	
	for(std::size_t i = 1, n = srcs.size(); i < n; ++i)
	{
		x = srcs[i].getX();
		
		if(x[0] < rangex[0])
			rangex[0] = x[0];
		else if(x[0] > rangex[1])
			rangex[1] = x[0];
		
		if(x[1] < rangey[0])
			rangey[0] = x[1];
		else if(x[1] > rangey[1])
			rangey[1] = x[1];
	}
	
	return (rangex[1]-rangex[0])*(rangey[1]-rangey[0])*180*180/pi/pi;
}
