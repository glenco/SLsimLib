#include "../include/source_type.h"
#include "../include/source.h"

SourceType::SourceType(const SourceType& other)
: t(other.t)
{
}

SourceType::SourceType(SourceType::Type type)
: t(type)
{
}

SourceType::SourceType(::Source* source)
: t(source->type().t)
{
}

bool SourceType::operator==(const SourceType& rhs)
{
	return (t == rhs.t);
}

bool SourceType::operator!=(const SourceType& rhs)
{
	return (t != rhs.t);
}
