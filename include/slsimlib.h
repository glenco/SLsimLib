/** \file
 * \ingroup \HighLevel
 *
 * \brief Master header file for all routines in SLsimLib.  Should be the only header file that needs to be included.
 *
 * assert() statements can be removed here.
 *
 * slsimlib.h
 *
 *  Created on: Oct 25, 2011
 *      Author: bmetcalf
 */

#ifndef _SLSIMLIB_DECLARE_
#define _SLSIMLIB_DECLARE_

//#define NDEBUG  // Un-commenting this line will remove all assert() statements in the executable.

#include <iostream>
#include <cstdio>
#include <fstream>
#include <cfloat>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <pthread.h>

#include <math.h>
#include <assert.h>
#include <time.h>
#include <limits>
#include <cstring>

#include <nrutil.h>
#include <nr.h>
#include <nrD.h>

#include <cosmo.h>
#include <halo.h>

#include <point.h>
#include <List1.h>
#include <Tree.h>
#include <Kist.h>
#include <KistDriver.h>

#include <lens.h>
#include <source.h>
#include <sourceAnaGalaxy.h>
#include <divide_images.h>
#include <tree_maintenance.h>
#include <grid_maintenance.h>
#include <image_info.h>
#include <source_models.h>
#include <analytic_lens.h>

#include <model.h>
#include <TreeNB.h>
#include <forceTree.h>
#include <simpleTree.h>
#include <multiplane.h>
#include <singlelens.h>
#include <MOKAlens.h>
#include <MOKAfits.h>
#include <quadTree.h>

#include <peak_refinement.h>
#include <map_images.h>

#include <fitlens.h>
#include <galaxies.h>
#include <image_processing.h>
#include <nsie.h>
#include <InputParameters.h>

#ifndef pi
#define pi  3.1415926
#endif

#ifndef Grav
#define Grav  4.7788e-20  // G/c^2 in Mpc
#endif

#endif
