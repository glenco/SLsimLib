#ifndef SKY_H_
#define SKY_H_

#include "source.h"
#include "utilities_slsim.h"
#include "InputParams.h"
#include "raw_data.h"

#include <vector>
#include <cstddef>

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
	
	/// returns field-of-view in deg^2 assuming region is square
	double getFov(); // TODO: make const once Source is sorted out
	
	/// Read data from all sources.
	void serialize(RawData& d) const;
	
	/// Write data to all sources.
	void unserialize(RawData& d);
	
	/// Randomize all sources.
	void randomize(double step, long* seed);
	
private:
	Utilities::MixedVector<Source> srcs;
	
	void checkType(const Source&) {}
	
private: /* input methods */
	void readGalaxyFile(std::string filename, Band band, double mag_limit);
};

#endif
