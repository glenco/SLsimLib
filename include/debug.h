#ifndef _slsimlib_debug_declare_
#define _slsimlib_debug_declare_

#ifndef NDEBUG
#include <iostream>
#define SLSIMLIB_DEBUG(x) do { std::cerr << x << std::endl; } while(0)
#else
#define SLSIMLIB_DEBUG(x) do {} while(0)
#endif

#endif
