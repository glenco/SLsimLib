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
struct Beam{
  std::vector<Point_2d> x;
  std::vector<KappaType> kappa;
  std::vector<Point_3d<> > gamma;
};

namespace FindImages {
  short find_peaks(LensHndl lens,GridHndl grid,double rEinsteinMin,double kappa_max,std::vector<ImageInfo> &imageinfo, int* Nimages);
}

#endif
