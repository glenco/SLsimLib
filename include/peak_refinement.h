/*
 * peak_refinement.h
 *
 *  Created on: Apr 14, 2011
 *      Author: bmetcalf
 */

#ifndef beamtypes_declare
#define beamtypes_declare

#include <grid_maintenance.h>

/** \brief A simplified data structure for use in interface with other codes. */
typedef struct Beam{
	double source[2];
	double image[2];
	double kappa;
	double gamma[2];
} Beam;

typedef Beam * BeamHndl;

namespace FindImages {
  short find_peaks(LensHndl lens,GridHndl grid,double rEinsteinMin,double kappa_max,std::vector<ImageInfo> &imageinfo, int* Nimages);
}

#endif
