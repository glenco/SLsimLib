/*
 * grid_maintenance.h
 *
 *  Created on: Oct 12, 2011
 *      Author: bmetcalf
 */

#ifndef _grid_maintenance_declare_
#define _grid_maintenance_declare_

#include "lens.h"
#include "Tree.h"

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
	unsigned long getNumberOfPoints();


	/// tree on image plane
	TreeHndl i_tree;
	/// tree on source plane
	TreeHndl s_tree;

	/// return initial number of grid points in each direction
	int getInitNgrid(){return Ngrid_init;}
	/// return number of cells in each dimension into which each cell is divided when a refinement is made
	int getNgrid_block(){return Ngrid_block;}
	/// return initial range of gridded region
	double getInitRange(){return i_tree->top->boundary_p2[0] - i_tree->top->boundary_p1[0];}
	Point * RefineLeaf(LensHndl lens,Point *point,bool kappa_off);
	Point * RefineLeaves(LensHndl lens,std::vector<Point *>& points,bool kappa_off);
	void ClearAllMarks();

	void test_mag_matrix();

private:
	/// one dimensional size of initial grid
	int Ngrid_init;
	/// one dimensional number of cells a cell will be divided into on each refinement step
	int Ngrid_block;
	bool initialized;
	KistHndl trashkist;

	double maglimit;
	KistHndl neighbors;
	bool find_mag_matrix(double *a,Point *p0,Point *p1,Point *p2);

	bool uniform_mag_from_deflect(double *a,Point *point);
	bool uniform_mag_from_shooter(double *a,Point *point);

};

typedef struct Grid* GridHndl;

// in image_finder_kist.c

void find_images_kist(LensHndl lens,double *y_source,double r_source,GridHndl grid
		,int *Nimages,ImageInfo *imageinfo,const int NimageMax,unsigned long *Nimagepoints
		,double initial_size,bool splitimages,short edge_refinement
		,bool verbose,bool kappa_off);

short image_finder_kist(LensHndl lens, double *y_source,double r_source,GridHndl grid
		,int *Nimages,ImageInfo *imageinfo,const int NimageMax,unsigned long *Nimagepoints
		,short splitparities,short true_images);

int refine_grid_kist(LensHndl lens,GridHndl grid,ImageInfo *imageinfo
		,unsigned long Nimages,double res_target,short criterion,bool kappa_off,KistHndl newpointkist = NULL,bool batch=true);

ImageInfo *find_crit(LensHndl lens,GridHndl grid,int *Ncrits,double resolution,bool *orderingsuccess
		,bool ordercurve,bool verbose);

void find_crit_kist(LensHndl lens,GridHndl grid,ImageInfo *critcurve,int maxNcrits,int *Ncrits
		,double resolution,bool *orderingsuccess,bool ordercurve,bool verbose);

int refine_grid(LensHndl lens,GridHndl grid,OldImageInfo *imageinfo
		,unsigned long Nimages,double res_target,short criterion,bool kappa_off,bool batch=true);
long refine_edges(LensHndl lens,GridHndl grid,ImageInfo *imageinfo
		,unsigned long Nimages,double res_target,short criterion,bool kappa_off,bool batch=true);
long refine_edges2(LensHndl lens,double *y_source,double r_source,GridHndl grid
		,ImageInfo *imageinfo,bool *image_overlap,unsigned long Nimages,double res_target
		,short criterion,bool kappa_off,bool batch=true);
void sort_out_points(Point *i_points,ImageInfo *imageinfo,double r_source,double y_source[]);
void xygridpoints(Point *points,double range,double *center,long Ngrid
		,short remove_center);
//void initialize_grid(double center[],double range,long Ngrid,TreeHndl s_tree,TreeHndl i_tree);

#endif
