/*
 * grid_maintenance.h
 *
 *  Created on: Oct 12, 2011
 *      Author: bmetcalf
 */

#ifndef _grid_maintenance_declare_
#define _grid_maintenance_declare_

#include <lens.h>
#include <source.h>
#include <Kist.h>
#include <Tree.h>
#include <image_info.h>

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

// in image_finder_kist.c

void find_images_kist(LensHndl lens,double *y_source,double r_source,GridHndl grid
		,int *Nimages,ImageInfo *imageinfo,const int NimageMax,unsigned long *Nimagepoints
		,double initial_size,bool splitimages,short edge_refinement
		,bool verbose,bool kappa_off);

short image_finder_kist(LensHndl lens, double *y_source,double r_source,GridHndl grid
		,int *Nimages,ImageInfo *imageinfo,const int NimageMax,unsigned long *Nimagepoints
		,short splitparities,short true_images);

int refine_grid_kist(LensHndl lens,GridHndl grid,ImageInfo *imageinfo
		,unsigned long Nimages,double res_target,short criterion,bool kappa_off,KistHndl newpointkist = NULL);

ImageInfo *find_crit(LensHndl lens,GridHndl grid,int *Ncrits,double resolution,bool *orderingsuccess
		,bool ordercurve,bool verbose);

void find_crit_kist(LensHndl lens,GridHndl grid,ImageInfo *critcurve,int maxNcrits,int *Ncrits
		,double resolution,bool *orderingsuccess,bool ordercurve,bool verbose);
#endif
