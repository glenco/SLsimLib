#ifndef PARAMETERS_MULTISOURCEANAGALAXY_H_
#define PARAMETERS_MULTISOURCEANAGALAXY_H_

#include "../../include/parameters.h"

#include "../../include/sourceAnaGalaxy.h"

#include "OverGalaxy.h"

#include <vector>
#include <cstddef>
#include <cassert>

template<>
struct Parameters<MultiSourceAnaGalaxy>
{
	std::vector<Parameters<OverGalaxy> > galaxies;
};

void operator<<(Parameters<MultiSourceAnaGalaxy>& p, const MultiSourceAnaGalaxy& s)
{
	p.galaxies.resize(s.getNumberOfGalaxies());
	
	for(std::size_t i = 0, n = p.galaxies.size(); i < n; ++i)
		p.galaxies[i] << s[i];
}

void operator>>(const Parameters<MultiSourceAnaGalaxy>& p, MultiSourceAnaGalaxy& s)
{
	assert(p.galaxies.size() == s.getNumberOfGalaxies());
	
	for(std::size_t i = 0, n = p.galaxies.size(); i < n; ++i)
		p.galaxies[i] >> s[i];
}

#endif
