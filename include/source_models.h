/*
 * source_models.h
 *
 *  Created on: Oct 8, 2011
 *      Author: bmetcalf
 */

#include <analytic_lens.h>

#ifndef source_declare
#define source_declare

// Separate files
double blr_surface_brightness_spherical(double x,ModelHndl model);
double blr_surface_brightness_spherical_random_motions(double x,ModelHndl model);
double blr_surface_brightness_spherical_circular_motions(double x,ModelHndl model);
double blr_surface_brightness_disk(double x[],ModelHndl model);

#endif
