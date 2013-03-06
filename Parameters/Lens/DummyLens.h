#ifndef PARAMETERS_DUMMYLENS_H_
#define PARAMETERS_DUMMYLENS_H_

#include "../../include/parameters.h"

class DummyLens;

template<>
class Parameters<DummyLens>
{
};

void operator<<(Parameters<DummyLens>&, const DummyLens&)
{
}

void operator>>(const Parameters<DummyLens>&, DummyLens&)
{
}

#endif
