#ifndef MULTI_SOURCE_H_
#define MULTI_SOURCE_H_

#include <vector>
#include <map>
#include <cstddef>

#include "source.h"

class MultiSource : public Source
{
public:
	SOURCE_TYPE(MultiSource)
	
	/// Construct an empty MultiSource
	MultiSource();
	
	/// Create a MultiSource from input parameters.
	MultiSource(InputParams& params);
	
	/// Copy a MultiSource and its contained sources.
	MultiSource(const MultiSource& other);
	
	/// Destroy the MultiSource.
	~MultiSource();
	
	/// Copy contents of another MultiSource into this. Replaces all sources, invalidates all pointers.
	MultiSource& operator=(MultiSource rhs);
	
	/// Set parameters into all sources.
	void setParameters(Parameters& p);
	
	/// Get parameters for all sources.
	void getParameters(Parameters& p);
	
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
	friend void swap(MultiSource& a, MultiSource& b);
	
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

inline double MultiSource::SurfaceBrightness(double *y)
{
	double sb = sources[index]->SurfaceBrightness(y);
	
	if(sb*hplanck < std::pow(10., -0.4*(48.6+getSBlimit())) * (180*60*60/pi)*(180*60*60/pi))
		return 0.;
	
	return sb;
}

inline double MultiSource::getTotalFlux() { return sources[index]->getTotalFlux(); }

inline double MultiSource::getRadius() { return sources[index]->getRadius(); }

inline void MultiSource::setRadius(double r) { sources[index]->setRadius(r); }

inline double MultiSource::getZ(){ return sources[index]->getZ(); }

inline void MultiSource::setZ(double z){ sources[index]->setZ(z); }

inline double* MultiSource::getX(){return sources[index]->getX(); }

inline void MultiSource::setX(double x[2]){ sources[index]->setX(x); }

inline void MultiSource::setX(double x1, double x2){ sources[index]->setX(x1, x2); }

inline std::size_t MultiSource::size() const { return sources.size(); }

inline Source* MultiSource::getCurrent() const { return sources[index]; }

inline std::size_t MultiSource::getIndex() const { return index; }

#endif
