#ifndef MULTI_SOURCE_H_
#define MULTI_SOURCE_H_

#include <vector>
#include <map>
#include <cstddef>

#include "source.h"

class SourceMulti : public Source
{
public:
	SOURCE_TYPE(SourceMulti)
	
	/// Construct an empty MultiSource
	SourceMulti();
	
	/// Create a MultiSource from input parameters.
	SourceMulti(InputParams& params);
	
	/// Copy a MultiSource and its contained sources.
	SourceMulti(const SourceMulti& other);
	
	/// Destroy the MultiSource.
	~SourceMulti();
	
	/// Copy contents of another MultiSource into this. Replaces all sources, invalidates all pointers.
	SourceMulti& operator=(SourceMulti rhs);
	
	/// Get parameters for all sources.
	void getParameters(Parameters& p) const;
	
	/// Set parameters into all sources.
	void setParameters(Parameters& p);
	
	/// Surface brightness of current source. The limits of both the MultiSource and the current source apply.
	double SurfaceBrightness(double* y);
	
	/// Total flux coming from the current source.
	double getTotalFlux();
	
	/// Output source information.
	void printSource();
	
	/// Return redshift of current source.
	double getZ();
	
	/// Set redshift of current source. Only changes the redshift while leaving position fixed.
	void setZ(double z);
	
	/// Get radius of current source.
	double getRadius();
	
	/// Set radius of current source.
	void setRadius(double r);
	
	/// Return angular position of current source.
	double* getX();
	
	/// Set angular position of current source.
	void setX(double x[2]);
	
	/// Set angular position of current source.
	void setX(double x1, double x2);
	
	/// Randomize all sources.
	void randomize(double step, long* seed);
	
	/// Get number of sources.
	std::size_t size() const;
	
	/// Add a source. A copy is created and stored in the MultiSource.
	std::size_t add(const Source& source);
	
	/// Add a source. A copy is created and stored in the MultiSource.
	std::size_t add(const Source* source);
	
	/// Get the current source
	Source* getCurrent() const;
	
	/// Set the current source. Needs to be in list of sources. Does a costly find.
	bool setCurrent(Source* source);
	
	/// Get the current index
	std::size_t getIndex() const;
	
	/// Set the current index. Needs to be in range.
	bool setIndex(std::size_t i);
	
	/// Get current source in a given type. Returns 0 if types don't match.
	template<typename SourceT>
	SourceT* getCurrent() const
	{
		// use the SourceType to safely upcast the current Source*
		return source_cast<SourceT*>(sources[index]);
	}
	
	/// Get the type of the current source.
	const SourceType getCurrentType() const
	{
		return sources[index]->type();
	}
	
	/// Get all sources.
	std::vector<Source*> getAll() const;
	
	/// Get all sources of a given type.
	std::vector<Source*> getAll(SourceType type) const;
	
	/// Get all sources of a given type in that type.
	template<typename SourceT>
	std::vector<SourceT*> getAll() const
	{
		// list of resulting sources
		std::vector<SourceT*> matches;
		
		// try to find SourceType in map
		std::map<SourceType, std::vector<Source*> >::const_iterator pos = type_map.find(source_type_of<SourceT>());
		
		// check if found
		if(pos == type_map.end())
			return std::vector<SourceT*>();
		
		// add sources to list of results
		for(std::vector<Source*>::const_iterator it = pos->second.begin(); it != pos->second.end(); ++it)
			matches.push_back((SourceT*)(*it));
		
		// return sources of requested type
		return matches;
	}
	
	/// Swap contents of MultiSource with another
	friend void swap(SourceMulti& a, SourceMulti& b);
	
private:
	void assignParams(InputParams& params);
	
	/// add a source to the internal list
	void addInternal(Source* source);
	
	/// read a Millenium galaxy data file
	void readGalaxyFile(std::string filename, Band band, double mag_limit);
	
	/// the current source index
	std::size_t index;
	
	/// list of sources
	std::vector<Source*> sources;
	
	/// map of sources and types
	std::map<SourceType, std::vector<Source*> > type_map;
};

/**** inline functions ****/

inline double SourceMulti::SurfaceBrightness(double *y)
{
	double sb = sources[index]->SurfaceBrightness(y);
	
	if(sb < sb_limit)
		return 0.;
	
	return sb;
}

inline double SourceMulti::getTotalFlux() { return sources[index]->getTotalFlux(); }

inline double SourceMulti::getRadius() { return sources[index]->getRadius(); }

inline void SourceMulti::setRadius(double r) { sources[index]->setRadius(r); }

inline double SourceMulti::getZ(){ return sources[index]->getZ(); }

inline void SourceMulti::setZ(double z){ sources[index]->setZ(z); }

inline double* SourceMulti::getX(){return sources[index]->getX(); }

inline void SourceMulti::setX(double x[2]){ sources[index]->setX(x); }

inline void SourceMulti::setX(double x1, double x2){ sources[index]->setX(x1, x2); }

inline std::size_t SourceMulti::size() const { return sources.size(); }

inline Source* SourceMulti::getCurrent() const { return sources[index]; }

inline std::size_t SourceMulti::getIndex() const { return index; }

#endif
