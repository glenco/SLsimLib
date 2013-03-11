#ifndef SLFITLIB_H
#define SLFITLIB_H

// basic types
#include "parameters.h"
#include "generator.h"
#include "mcmc.h"

// image-based mcmc
#include "image_processing.h"
#include "image_likelihood.h"

/****
 * parameter specializations
 */

// sources
#include "../Parameters/Source/OverGalaxy.h"
#include "../Parameters/Source/MultiSourceAnaGalaxy.h"

// lenses
#include "../Parameters/Lens/DummyLens.h"

#endif
