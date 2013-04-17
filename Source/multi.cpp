#include "../include/multi_source.h"
#include "../include/InputParams.h"

MultiSource::MultiSource()
: index(0)
{
	sb_limit = 30.;
}

MultiSource::MultiSource(InputParams& params)
: index(0)
{
	// check if there is a sb_limit set
	if(!params.get("source_sb_limit", sb_limit))
		sb_limit = 30.;
}

MultiSource::MultiSource(const MultiSource& other)
: Source(other), // base copy constructor
  index(other.index), sources(other.sources), created()
{
}

MultiSource::~MultiSource()
{
	// delete created sources
	for(std::size_t i = 0, n = created.size(); i < n; ++i)
		delete created[i];
}

bool MultiSource::setCurrent(Source* source)
{
	// find source in list
	std::size_t i = 0, n = sources.size();
	while(i < n && sources[i].source != source)
		++i;
	
	// check if source was found
	if(i == n)
		return false;
	
	// set index
	index = i;
	
	// all is well
	return true;
}

bool MultiSource::setIndex(std::size_t i)
{
	// make sure index is in range
	if(i >= sources.size())
		return false;
	
	// set index
	index = i;
	
	// all is well
	return true;
}

std::vector<Source*> MultiSource::getAll() const
{
	// list of results
	std::vector<Source*> matches;
	
	// go through list
	for(std::size_t i = 0, n = sources.size(); i < n; ++i)
		matches.push_back(sources[i].source);
	
	// return sources
	return matches;
}

void MultiSource::printSource()
{
	if(sources.empty())
		std::cout << "Empty MultiSource" << std::endl;
	else
	{
		std::cout << "MultiSource\n-----------" << std::endl;
		sources[index].source->printSource();
	}
}

void MultiSource::assignParams(InputParams& /* params */)
{
}

void MultiSource::addInternal(Source* source, const std::type_info& type, bool owned)
{
	// add source to list of created sources if owned by this MultiSource
	if(owned)
		created.push_back(source);
	
	// add source and type to list
	sources.push_back(SourceTypePair(source, type));
}
