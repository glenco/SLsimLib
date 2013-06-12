#include "../include/multi_source.h"
#include "../include/InputParams.h"

#include <algorithm>
#include <utility>
#include <iterator>

void swap(SourceMulti& a, SourceMulti& b)
{
	using std::swap;
	
	swap(a.sources, b.sources);
	swap(a.type_map, b.type_map);
	swap(a.index, b.index);
	
	// Source stuff
	swap(a.source_r, b.source_r);
	swap(a.source_x[0], b.source_x[0]);
	swap(a.source_x[1], b.source_x[1]);
	swap(a.zsource, b.zsource);
	swap(a.DlDs, b.DlDs);
	swap(a.sb_limit, b.sb_limit);
}

SourceMulti::SourceMulti()
: index(0)
{
	setSBlimit_magarcsec(30.);
}

SourceMulti::SourceMulti(InputParams& params)
: index(0)
{
	// check if there is a sb_limit set
	if(!params.get("source_sb_limit", sb_limit))
		setSBlimit_magarcsec(30.);
	else
		sb_limit = pow(10,-0.4*(48.6+sb_limit))*pow(180*60*60/pi,2)/hplanck;

	
	// check if there is a Millennium data file
	std::string input_galaxy_file;
	if(params.get("input_galaxy_file", input_galaxy_file))
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
		readGalaxyFile(input_galaxy_file, band, mag_limit);
	}
}

SourceMulti::SourceMulti(const SourceMulti& other)
: Source(other)
{
	// copy all contained sources
	for(std::size_t i = 0, n = other.sources.size(); i < n; ++i)
		addInternal(other.sources[i]->clone());
	
	// set index
	index = other.index;
}

SourceMulti::~SourceMulti()
{
	// delete sources
	for(std::size_t i = 0, n = sources.size(); i < n; ++i)
		delete sources[i];
}

SourceMulti& SourceMulti::operator=(SourceMulti rhs)
{
	swap(*this, rhs);
	return *this;
}

void SourceMulti::getParameters(Parameters& p) const
{
	Source::getParameters(p);
	for(std::size_t i = 0, n = sources.size(); i < n; ++i)
		sources[i]->getParameters(p);
}

void SourceMulti::setParameters(Parameters& p)
{
	Source::setParameters(p);
	for(std::size_t i = 0, n = sources.size(); i < n; ++i)
		sources[i]->setParameters(p);
}

void SourceMulti::randomize(double step, long* seed)
{
	for(std::size_t i = 0, n = sources.size(); i < n; ++i)
		sources[i]->randomize(step, seed);
}

bool SourceMulti::setCurrent(Source* source)
{
	// find source in list
	std::vector<Source*>::iterator pos = std::find(sources.begin(), sources.end(), source);
	
	// check if source was found
	if(pos == sources.end())
		return false;
	
	// set index
	index = (std::size_t)std::distance(sources.begin(), pos);
	
	// all is well
	return true;
}

bool SourceMulti::setIndex(std::size_t i)
{
	// make sure index is in range
	if(i >= sources.size())
		return false;
	
	// set index
	index = i;
	
	// all is well
	return true;
}

std::size_t SourceMulti::add(const Source& source)
{
	// get the last index
	std::size_t i = sources.size();
	
	// add source to internal lists
	addInternal(source.clone());
	
	// return the added source
	return i;
}

std::size_t SourceMulti::add(const Source* source)
{
	// get the last index
	std::size_t i = sources.size();
	
	// add source to internal lists
	addInternal(source->clone());
	
	// return the added source
	return i;
}

std::vector<Source*> SourceMulti::getAll() const
{
	return sources;
}

std::vector<Source*> SourceMulti::getAll(SourceType type) const
{
	// try to find SourceType in map
	std::map<SourceType, std::vector<Source*> >::const_iterator pos = type_map.find(type);
	
	// check if found
	if(pos == type_map.end())
		return std::vector<Source*>();
	
	// return sources found
	return pos->second;
}

void SourceMulti::printSource()
{
	if(sources.empty())
		std::cout << "Empty MultiSource" << std::endl;
	else
	{
		std::cout << "MultiSource\n-----------" << std::endl;
		sources[index]->printSource();
	}
}

void SourceMulti::assignParams(InputParams& /* params */)
{
}

void SourceMulti::addInternal(Source* source)
{
	// add source
	sources.push_back(source);
	
	// add to type map
	type_map[source->type()].push_back(source);
}
