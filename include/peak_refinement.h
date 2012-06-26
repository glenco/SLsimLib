/*
 * peak_refinement.h
 *
 *  Created on: Apr 14, 2011
 *      Author: bmetcalf
 */

#ifndef beamtypes_declare
#define beamtypes_declare

#include <lens.h>
#include <grid_maintenance.h>
#include <image_info.h>

/** \brief A simplified data structure for use in interface with other codes. */
typedef struct Beam{
	double source[2];
	double image[2];
	double kappa;
	double gamma[2];
} Beam;

typedef Beam * BeamHndl;

short find_peaks(LensHndl lens,GridHndl grid,double rEinsteinMin,double kappa_max,ImageInfo *imageinfo, int* Nimages);

#endif
