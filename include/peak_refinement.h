/*
 * peak_refinement.h
 *
 *  Created on: Apr 14, 2011
 *      Author: bmetcalf
 */

#include <Tree.h>

#ifndef beamtypes_declare
#define beamtypes_declare

typedef struct Beam{
	double source[2];
	double image[2];
	double kappa;
	double gamma[2];
} Beam;

typedef Beam * BeamHndl;

BeamHndl find_peaks(double center[],double range,unsigned long Ngrid,double rEinsteinMin,double kappa_max
		,unsigned long *Nbeams,void (*rayshooter)(unsigned long N,Point *i_point));
void copyPointToBeam(ListHndl pointlist,Beam *beam);
void copyBeamToPoint(Beam *beam,Point *pointarr,unsigned long N);

#endif
