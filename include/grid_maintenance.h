/*
 * grid_maintenance.h
 *
 *  Created on: Oct 12, 2011
 *      Author: bmetcalf
 */


#ifndef _grid_maintenance_declare_
#define _grid_maintenance_declare_

#include <model.h>

/** \ingroup ImageFinding
 * \brief Structure to contain both source and image trees.
 * It is not yet used, but may be useful.
 */
typedef struct Grid{

	Grid(LensHndl lens,int N1d,double center[2],double range);
	~Grid();

	void ReInitializeGrid(LensHndl lens);
	void TrimGrid(double highestres,bool useSB);
	void RefreshSurfaceBrightnesses(SourceHndl source);
	unsigned long NumberOfPoints();

	/// tree on image plane
	TreeHndl i_tree;
	/// tree on source plane
	TreeHndl s_tree;
	/// one dimensional size of grid
	int Ngrid;
	bool initialized;
};

typedef struct Grid *GridHndl;

//GridHndl NewGrid(LensHndl lens,int Ngrid,double center[2],double range);
//void FreeGrid(GridHndl grid);

#endif
