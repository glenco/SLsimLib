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

#include "utilities_slsim.h"

#include "lens_halos.h"

#include <KistDriver.h>

#include <divide_images.h>
#include <tree_maintenance.h>
#include <grid_maintenance.h>
#include <source_models.h>

#include <model.h>

#include <peak_refinement.h>
#include <map_images.h>

#include <fitlens.h>
#include <image_processing.h>
#include <nsie.h>

/**** halos ****/

#include "analytic_lens.h"
#include "uniform_lens.h"
#include "MOKAlens.h"

/**** sources ****/

#include "sourceAnaGalaxy.h"
#include "overzier_source.h"
#include "sersic_source.h"
#include "causticdata.h"

#endif
