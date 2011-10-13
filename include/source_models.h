/*
 * source_models.h
 *
 *  Created on: Oct 8, 2011
 *      Author: bmetcalf
 */

#include <analytic_lens.h>
#include <cosmo.h>

#ifndef source_declare
#define source_declare

// Separate files
double blr_surface_brightness_spherical(double x,AnaLens *lens);
double blr_surface_brightness_spherical_random_motions(double x,AnaLens *lens,COSMOLOGY *cosmo);
double blr_surface_brightness_spherical_circular_motions(double x,AnaLens *lens,COSMOLOGY *cosmo);
double blr_surface_brightness_disk(double x[],AnaLens *lens,COSMOLOGY *cosmo);

#endif
