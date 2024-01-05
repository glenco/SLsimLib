/*
 * source_models.h
 *
 *  Created on: Oct 8, 2011
 *      Author: bmetcalf
 */

#ifndef source_model_declare
#define source_model_declare

#include "source.h"

// Separate files
double blr_surface_brightness_spherical(double x,const SourceBLR *source);
double blr_surface_brightness_spherical_random_motions(double x,const SourceBLR *source);
double blr_surface_brightness_spherical_circular_motions(double x,const SourceBLR *source);
double blr_surface_brightness_disk(double x[],const SourceBLR *source);

#endif
