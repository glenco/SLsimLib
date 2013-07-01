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
	
	// check if there is a Sersic catalog file
	if(params.get("input_sersic_file", filename))
	{
		// source redshift
		double z_source;
		if(!params.get("z_source", z_source))
		{
			std::cerr << "ERROR: Must assign z_source in parameter file " << params.filename() << std::endl;
			exit(1);
		}
		
		// load Sersic sources from file
		readSersicFile(filename, z_source);
	}
}

Sky::~Sky()
{
}

void Sky::serialize(RawData& d) const
{
	for(std::size_t i = 0, n = srcs.size(); i < n; ++i)
		srcs[i].serialize(d);
}

void Sky::unserialize(RawData& d)
{
	for(std::size_t i = 0, n = srcs.size(); i < n; ++i)
		srcs[i].unserialize(d);
}

void Sky::randomize(double step, long* seed)
{
	//lns->randomize(step, seed);
	for(std::size_t i = 0, n = srcs.size(); i < n; ++i)
		srcs[i].randomize(step, seed);
}
