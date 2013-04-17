#ifndef MULTI_SOURCE_H_
#define MULTI_SOURCE_H_

#include <vector>
#include <utility>
#include <typeinfo>
#include <cstddef>

#include "source.h"

class MultiSource : public Source
{
public:
	/// Construct an empty MultiSource
	MultiSource();
	
	/// Create a MultiSource from input parameters.
	MultiSource(InputParams& params);
	
	/// Copy a MultiSource. The new MultiSource does not own any sources.
	MultiSource(const MultiSource& other);
	
	/// Destroy the MultiSource and free created sources.
	~MultiSource();
	
	/// Surface brightness of current source
	double SurfaceBrightness(double* y);
	
	/// Total flux coming from the current source
	double getTotalFlux();
	
	/// Get radius of current source.
	double getRadius();
	
	/// Return redshift of current source.
	double getZ();
	
	/// Set redshift of current source. Only changes the redshift while leaving position fixed.
	void setZ(double z);
	
	/// Return angular position of current source.
	double* getX();
	
	/// Set angular position of current source.
	void setX(double theta[2]);
	
	/// Set angular position of current source.
	void setX(double x1, double x2);
	
	/// Get number of sources.
	std::size_t size() const;
	
	/// Add a source. Ownership is with the MultiSource.
	template<typename SourceT>
	SourceT* add()
	{
		// create the new source of the given type
		SourceT* source = new SourceT();
		
		// add source to internal lists
		addInternal(source, typeid(SourceT), true);
		
		// return the created source
		return source;
	}
	
	/// Add a source using input parameters. Ownership is with the MultiSource.
	template<typename SourceT>
	SourceT* add(InputParams& params)
	{
		// create the new source of the given type
		SourceT* source = new SourceT(params);
		
		// add source to internal lists
		addInternal(source, typeid(SourceT), true);
		
		// return the created source
		return source;
	}
	
	/// Add a source. Ownership is not transferred to the MultiSource.
	template<typename SourceT>
	SourceT* add(SourceT* source)
	{
		// add source to internal lists
		addInternal(source, typeid(SourceT), false);
		
		// return the added source
		return source;
	}
	
	/// Add a source. Ownership is not transferred to the MultiSource.
	template<typename SourceT>
	SourceT* add(SourceT& source)
	{
		// add source to internal lists
		addInternal(&source, typeid(SourceT), false);
		
		// return the added source
		return &source;
	}
	
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
		// use the RTTI to safely upcase Source* current
		return dynamic_cast<SourceT*>(sources[index].source);
	}
	
	/// Get the type of the current source.
	const std::type_info& getCurrentType() const
	{
		return sources[index].type;
	}
	
	/// Get all sources.
	std::vector<Source*> getAll() const;
	
	/// Get all sources of a given type.
	template<typename SourceT>
	std::vector<SourceT*> getAll() const
	{
		// list of results
		std::vector<SourceT*> matches;
		
		// go through list
		for(std::size_t i = 0, n = sources.size(); i < n; ++i)
		{
			// check if type of source matches
			if(sources[i].type == typeid(SourceT))
				matches.push_back(dynamic_cast<SourceT*>(sources[i].source));
		}
		
		// return sources of requested type
		return matches;
	}
	
	void printSource();
	
private:
	void assignParams(InputParams& params);
	
	/// add a source to the internal list
	void addInternal(Source* source, const std::type_info& type, bool owned);
	
	/// read a Millenium galaxy data file
	void readGalaxyFile(std::string filename, Band band, double mag_limit);
	
	/// the current source index
	std::size_t index;
	
	// TODO: make this a map of type_index when we move to C++11
	struct SourceTypePair
	{
		Source* source;
		const std::type_info& type;
		
		SourceTypePair(Source* s, const std::type_info& t) : source(s), type(t) {}
	};
	
	/// list of sources and associated types
	std::vector<SourceTypePair> sources;
	
	// TODO: handle creation better
	std::vector<Source*> created;
};

/**** inline functions ****/

inline double MultiSource::SurfaceBrightness(double *y) { return sources[index].source->SurfaceBrightness(y); }

inline double MultiSource::getTotalFlux() { return sources[index].source->getTotalFlux(); }

inline double MultiSource::getRadius() { return sources[index].source->getRadius(); }

inline double MultiSource::getZ(){ return sources[index].source->getZ(); }

inline void MultiSource::setZ(double z){ sources[index].source->setZ(z); }

inline double* MultiSource::getX(){return sources[index].source->getX(); }

inline void MultiSource::setX(double x[2]){ sources[index].source->setX(x); }

inline void MultiSource::setX(double x1, double x2){ sources[index].source->setX(x1, x2); }

inline std::size_t MultiSource::size() const { return sources.size(); }

inline Source* MultiSource::getCurrent() const { return sources[index].source; }

inline std::size_t MultiSource::getIndex() const { return index; }

#endif
