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

BeamHndl find_peaks(double center[],double range,unsigned long Ngrid,double rEinsteinMin,double kappa_max
		,unsigned long *Nbeams,void (*rayshooter)(unsigned long N,Beam *beam));
void copyPointToBeam(ListHndl pointlist,Beam *beam);
void copyBeamToPoint(Beam *beam,Point *pointarr,unsigned long N);
void copyPointArrayToBeams(Point *pointarr,Beam *beams,unsigned long N);
void RayShooterRap(unsigned long N,Point *points,void (*rayshooter)(unsigned long N,Beam *beams));

#endif
