/*
 * source_models.h
 *
 *  Created on: Oct 8, 2011
 *      Author: bmetcalf
 */

#ifndef source_declare
#define source_declare

#include <lens.h>

// Separate files
double blr_surface_brightness_spherical(double x,SourceBLR *source);
double blr_surface_brightness_spherical_random_motions(double x,SourceBLR *source);
double blr_surface_brightness_spherical_circular_motions(double x,SourceBLR *source);
double blr_surface_brightness_disk(double x[],SourceBLR *source);

#endif
