#ifndef SKY_H_
#define SKY_H_

#include "source.h"
#include "utilities_slsim.h"

#include <vector>
#include <cstddef>

class InputParams;
class Parameters;

/**
 * \brief Container for all sources that make up one sky.
 * 
 * This class contains all the sources belonging to the sky of a given model.
 */
class Sky
{
public: /* iterators */
	template<typename SourceT = Source>
	class iterator : public Utilities::MixedVector<Source>::iterator<SourceT>
	{
	public:
		iterator(const Utilities::MixedVector<Source>::iterator<SourceT>& i) : Utilities::MixedVector<Source>::iterator<SourceT>(i) {}
	};
	
public:
	/// Construct an empty Sky.
	Sky();
	
	/// Construct a new Sky from input parameters.
	Sky(InputParams& params);
	
	/// Delete the Sky and its contents.
	~Sky();
	
	/// Add a source. It is copied into the Sky.
	template<typename SourceT>
	void addSource(const SourceT& source)
	{
		// make sure type is right
		checkType(source);
		
		// add to internal list
		srcs.push_back(source);
	}
	
	/// Get source by index.
	Source& getSource(std::size_t i)
	{
		return srcs.get(i);
	}
	
	/// Get source by index (const).
	const Source& getSource(std::size_t i) const
	{
		return srcs.get(i);
	}
	
	/// Get source of type SourceT by index.
	template<typename SourceT>
	SourceT& getSource(std::size_t i)
	{
		return srcs.get<SourceT>(i);
	}
	
	/// Get source of type SourceT by index (const).
	template<typename SourceT>
	const SourceT& getSource(std::size_t i) const
	{
		return srcs.get<SourceT>(i);
	}
	
	/// Set active source by index.
	void setSource(std::size_t i);
	
	/// Get iterator to first of all sources
	iterator<> begin()
	{
		return srcs.begin();
	}
	
	/// Get iterator to last of all sources
	iterator<> end()
	{
		return srcs.end();
	}
	
	/// Get iterator to first of sources of type SourceT
	template<typename SourceT>
	iterator<SourceT> begin()
	{
		return srcs.begin<SourceT>();
	}
	
	/// Get iterator to last of sources of type SourceT
	template<typename SourceT>
	iterator<SourceT> end()
	{
		return srcs.end<SourceT>();
	}
	
	/// Get number of sources of all types.
	std::size_t getNsources() const
	{
		return srcs.size();
	}
	
	/// Get number of sources of type SourceT.
	template<typename SourceT>
	std::size_t getNsources() const
	{
		return srcs.size<SourceT>();
	}
	
	/// Set parameters for all sources.
	void setParameters(Parameters& p);
	
	/// Get parameters from all sources.
	void getParameters(Parameters& p);
	
	/// Randomize all sources and lenses.
	void randomize(double step, long* seed);
	
private:
	Utilities::MixedVector<Source> srcs;
	
	void checkType(const Source&) {}
	
private: /* input methods */
	void readGalaxyFile(std::string filename, Band band, double mag_limit);
	void readSersicFile(std::string filename, double z);
};

#endif
