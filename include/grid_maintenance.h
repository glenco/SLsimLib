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
	unsigned long PruneTrees(double resolution,bool useSB,double fluxlimit);
	unsigned long PrunePointsOutside(double resolution,double *y,double r_in ,double r_out);

	double RefreshSurfaceBrightnesses(SourceHndl source);
	unsigned long NumberOfPoints();


	/// tree on image plane
	TreeHndl i_tree;
	/// tree on source plane
	TreeHndl s_tree;

	/// return initial number of grid points in each direction
	int getNgrid(){return Ngrid;}
	/// return number of cells in each dimension into which each cell is divided when a refinement is made
	int getNgrid_block(){return Ngrid_block;}
	/// return initial range of gridded region
	double getInitRange(){return i_tree->top->boundary_p2[0] - i_tree->top->boundary_p1[0];}

private:
	/// one dimensional size of grid
	int Ngrid;
	int Ngrid_block;
	bool initialized;
	KistHndl trashkist;
};

typedef struct Grid *GridHndl;

//GridHndl NewGrid(LensHndl lens,int Ngrid,double center[2],double range);
//void FreeGrid(GridHndl grid);

#endif
