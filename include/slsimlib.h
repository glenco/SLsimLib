/** \file
 * \ingroup \HighLevel
 *
 * \brief Header file for all routines in SLsimLib.
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

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <limits.h>
#include <string.h>

#include <nrutil.h>
#include <nr.h>
#include <nrD.h>

#include <cosmo.h>
#include <point.h>
#include <List.h>
#include <Tree.h>
#include <Kist.h>
#include <KistDriver.h>
#include <peak_refinement.h>
#include <analytic_lens.h>
#include <tree_maintenance.h>
#include <divide_images.h>
#include <tree_maintenance.h>
#include <grid_maintenance.h>
#include <map_images.h>
#include <source_models.h>

#include <TreeNB.h>

#include <fitlens.h>
#include <galaxies.h>
#include <image_processing.h>
#include <nsie.h>

#ifndef pi
#define pi  3.1415926
#endif

#endif
