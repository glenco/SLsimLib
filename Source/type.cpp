#include "../include/source_type.h"
#include "../include/source.h"

SourceType::SourceType(SourceType::Type type)
: t(type)
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
