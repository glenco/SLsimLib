/*
 * peak_refinement.h
 *
 *  Created on: Apr 14, 2011
 *      Author: bmetcalf
 */

#include <Tree.h>

#ifndef beamtypes_declare
#define beamtypes_declare

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
