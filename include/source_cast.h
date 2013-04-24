#ifndef SOURCE_CAST_H_
#define SOURCE_CAST_H_

#include <stdexcept>

/// Exception thrown when a source is trying to be casted from a bad type.
template<typename SourceT>
class BadSourceCast : public std::runtime_error
{
public:
	BadSourceCast(const Source& s) : std::runtime_error(std::string() + "Bad source cast from " + s.name() + " to " + source_name_of<SourceT>() + ".") {}
};

template<typename SourceT>
class BadSourceCast<const SourceT> : public std::runtime_error
{
public:
	BadSourceCast(const Source& s) : std::runtime_error(std::string() + "Bad source cast from " + s.name() + " to " + source_name_of<SourceT>() + ".") {}
};

template<typename SourceT>
class BadSourceCast<SourceT&> : public std::runtime_error
{
public:
	BadSourceCast(const Source& s) : std::runtime_error(std::string() + "Bad source cast from " + s.name() + " to " + source_name_of<SourceT>() + ".") {}
};

template<typename SourceT>
class BadSourceCast<const SourceT&> : public std::runtime_error
{
public:
	BadSourceCast(const Source& s) : std::runtime_error(std::string() + "Bad source cast from " + s.name() + " to " + source_name_of<SourceT>() + ".") {}
};

/// Cast a source reference into a given type. Throws BadSourceCast on error.
template<typename RefSourceT>
RefSourceT source_cast(Source& s)
{
	if(s.type() != source_type_of<RefSourceT>())
		throw BadSourceCast<RefSourceT>(s);
	return (RefSourceT)s;
}

/// Cast a const source reference into a given type. Throws BadSourceCast on error.
template<typename ConstRefSourceT>
ConstRefSourceT source_cast(const Source& s)
{
	if(s.type() != source_type_of<ConstRefSourceT>())
		throw BadSourceCast<ConstRefSourceT>(s);
	return (ConstRefSourceT)s;
}

/// Cast a source pointer into a given type. Returns 0 on error.
template<typename PtrSourceT>
PtrSourceT source_cast(Source* s)
{
	if(s->type() != source_type_of<PtrSourceT>())
		return 0;
	return (PtrSourceT)s;
}

/// Cast a const source pointer into a given type. Returns 0 on error.
template<typename ConstPtrSourceT>
ConstPtrSourceT source_cast(const Source* s)
{
	if(s->type() != source_type_of<ConstPtrSourceT>())
		return 0;
	return (ConstPtrSourceT)s;
}

#endif
