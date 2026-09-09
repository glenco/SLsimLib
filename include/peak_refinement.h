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
  Point_3d<> x;   // position of the end of beam
  Point_3d<> dx;  // tangent vector of beam at end
  Point_3d<> dxo; // tangent vector of beam at start
  Matrix2x2<double> A;  // beam magnification matrix

  void propogate(double ds){
    x+=dx*ds;
  }
};

namespace FindImages {
  short find_peaks(LensHndl lens,GridHndl grid,double rEinsteinMin,double kappa_max,std::vector<ImageInfo> &imageinfo, int* Nimages);
}

#endif
