#ifndef SOURCE_TYPE_H_
#define SOURCE_TYPE_H_

#include <functional>

class Source;

/// Macro for easy source type definition.
#define SOURCE_TYPE(T_) \
	SourceType type() const { return SourceType(SourceType::T_); } \
	static SourceType source_type() { return SourceType(SourceType::T_); }

/// Class registry for sources of complete type.
class SourceType
{
public:
	enum Type
	{
		MultiSource,
		MultiSourceAnaGalaxy,
		SersicSource,
		OverzierSource,
		PixelledSource,
		SourceGaussian,
		SourceUniform,
		SourceBLRDisk,
		SourceBLRSph1,
		SourceBLRSph2
	};
	
	/// Get type of a source class.
	template<typename SourceT>
	static SourceType typeof()
	{
		return SourceT::source_type();
	}
	
	SourceType(const SourceType&);
	SourceType(SourceType::Type);
	SourceType(::Source*);
	
	bool operator==(const SourceType&);
	bool operator!=(const SourceType&);
	
private:
	Type t;
	
	/// Friend class template specialization of std::less
	friend class std::less<SourceType>;
};

/// Get the type of a source class. Must have been registered with the SOURCE_TYPE macro.
template<typename SourceT>
SourceType source_type_of()
{
	return SourceT::source_type();
}

/// Comparator specialization for SourceType.
template<>
class std::less<SourceType>
{
public:
	bool operator()(const SourceType& lhs, const SourceType& rhs) const
	{
		return lhs.t < rhs.t;
	}
};

#endif
