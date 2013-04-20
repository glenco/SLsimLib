#ifndef SOURCE_TYPE_H_
#define SOURCE_TYPE_H_

#include <functional>

class Source;

/// Macro for easy source type definition.
#define SOURCE_TYPE(T_) \
	static SourceType source_type() { return SourceType(SourceType::T_); } \
	static const char* source_name() { return #T_; } \
	SourceType type() const { return SourceType(SourceType::T_); } \
	const char* name() const { return #T_; } \
	T_* clone() const { return new T_(*this); }

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
	
	SourceType(SourceType::Type);
	
	bool operator==(const SourceType&);
	bool operator!=(const SourceType&);
	
private:
	Type t;
	
	/// Friend class template specialization of std::less
	friend class std::less<SourceType>;
};

/// Get the type of a source class. Must have been registered with the SOURCE_TYPE macro.
template<typename SourceT>
inline SourceType source_type_of()
{
	return SourceT::source_type();
}

/// Get the name of a source class. Must have been registered with the SOURCE_TYPE macro.
template<typename SourceT>
inline const char* source_name_of()
{
	return SourceT::source_name();
}

namespace std
{
/// Comparator specialization for SourceType.
template<>
class less<SourceType>
{
public:
	bool operator()(const SourceType& lhs, const SourceType& rhs) const
	{
		return lhs.t < rhs.t;
	}
};
}

#endif
