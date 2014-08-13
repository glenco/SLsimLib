/*
 * grid_maintenance.h
 *
 *  Created on: Oct 12, 2011
 *      Author: bmetcalf
 */

#ifndef _grid_maintenance_declare_
#define _grid_maintenance_declare_

#include "lens.h"
#include "point.h"
#include "Tree.h"
#include <mutex>

class LensHaloBaseNSIE;
class LensHaloMOKA;

/** \ingroup ImageFinding
 * \brief Structure to contain both source and image trees.
 * It is not yet used, but may be useful.
 */
struct Grid{

	Grid(LensHndl lens,unsigned long N1d,const double center[2],double range);
  Grid(LensHndl lens ,unsigned long Nx ,const PosType center[2] ,PosType rangeX ,PosType rangeY);
	~Grid();

	void ReInitializeGrid(LensHndl lens);
    void ReShoot(LensHndl lens);
	void zoom(LensHndl lens,double *center,double scale,Branch *top = NULL);

	unsigned long PruneTrees(double resolution,bool useSB,double fluxlimit);
	unsigned long PrunePointsOutside(double resolution,double *y,double r_in ,double r_out);

	double RefreshSurfaceBrightnesses(SourceHndl source);
    double ClearSurfaceBrightnesses();
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
	Point * RefineLeaf(LensHndl lens,Point *point);
	Point * RefineLeaves(LensHndl lens,std::vector<Point *>& points);
	void ClearAllMarks();

	void test_mag_matrix();
  void writeFits(const double center[],size_t Npixels,double resolution,LensingVariable lensvar,std::string filename);
  void writeFits(const double center[],size_t Nx,size_t Ny,double resolution,LensingVariable lensvar,std::string filename);
  void writeFitsVector(const double center[],size_t Npixels,double resolution,LensingVariable lensvar,std::string filename);
  PixelMap writePixelMap(const double center[],size_t Npixels,double resolution,LensingVariable lensvar);
  PixelMap writePixelMap(const double center[],size_t Nx,size_t Ny,double resolution,LensingVariable lensvar);
  PixelMap writePixelMapUniform(const PosType center[],size_t Nx,size_t Ny,LensingVariable lensvar);

  void xygridpoints(Point *points,double range,const double *center,long Ngrid
                          ,short remove_center);

private:
	/// one dimensional size of initial grid
	const int Ngrid_init;
  int Ngrid_init2;
  
	/// one dimensional number of cells a cell will be divided into on each refinement step
	const int Ngrid_block;
	bool initialized;
	Kist<Point> * trashkist;

	double maglimit;
	Kist<Point> * neighbors;
	bool find_mag_matrix(double *a,Point *p0,Point *p1,Point *p2);

	bool uniform_mag_from_deflect(double *a,Point *point);
	bool uniform_mag_from_shooter(double *a,Point *point);
  
  unsigned long pointID;
  PosType axisratio;
  void writePixelMapUniform_(PointList list,PixelMap *map,LensingVariable val);
  
  static std::mutex grid_mutex;
};

typedef struct Grid* GridHndl;

// in image_finder_kist.c
namespace ImageFinding{
    
void find_images_kist(LensHndl lens,double *y_source,double r_source,GridHndl grid
		,int *Nimages,ImageInfo *imageinfo,int NimageMax,unsigned long *Nimagepoints
		,double initial_size,bool splitimages,short edge_refinement
		,bool verbose = false);

void find_images_microlens(LensHndl lens,double *y_source,double r_source,GridHndl grid
		,int *Nimages,ImageInfo *imageinfo,const int NimageMax,unsigned long *Nimagepoints
		,double initial_size,double mu_min,bool splitimages,short edge_refinement
		,bool verbose);

void find_images_microlens_exper(LensHndl lens,PosType *y_source,PosType r_source
        ,GridHndl grid,int *Nimages,ImageInfo *imageinfo,const int NimageMax
        ,unsigned long *Nimagepoints,PosType initial_size ,PosType mu_min
        ,bool splitimages,short edge_refinement,bool verbose);

void image_finder_kist(LensHndl lens, double *y_source,double r_source,GridHndl grid
		,int *Nimages,ImageInfo *imageinfo,const int NimageMax,unsigned long *Nimagepoints
		,short splitparities,short true_images);

int refine_grid_kist(LensHndl lens,GridHndl grid,ImageInfo *imageinfo
		,int Nimages,double res_target,short criterion
		,Kist<Point> * newpointkist = NULL,bool batch=true);

void find_crit(LensHndl lens,GridHndl grid,ImageInfo *critcurve,int maxNcrits,int *Ncrits
		,double resolution,bool *orderingsuccess,bool ordercurve,bool dividecurves,double invmag_min = 0.0,bool verbose = false);
void find_crit2(LensHndl lens,GridHndl grid,ImageInfo *critcurve,int maxNcrits,int *Ncrits
		,double resolution,bool *orderingsuccess,bool ordercurve,bool dividecurves,double invmag_min = 0.0,bool verbose = false);
}
void refine_crit_in_image(LensHndl lens,GridHndl grid,double r_source,double x_source[],double resolution);

int refine_grid(LensHndl lens,GridHndl grid,OldImageInfo *imageinfo
		,unsigned long Nimages,double res_target,short criterion,bool batch=true);

long refine_edges(LensHndl lens,GridHndl grid,ImageInfo *imageinfo
		,int Nimages,double res_target,short criterion
		,Kist<Point> * newpointkist = NULL,bool batch=true);

long refine_edges2(LensHndl lens,double *y_source,double r_source,GridHndl grid
		,ImageInfo *imageinfo,bool *image_overlap,int Nimages,double res_target
		,short criterion,bool batch=true);

void sort_out_points(Point *i_points,ImageInfo *imageinfo,double r_source,double y_source[]);


void saveImage(LensHaloMOKA *mokahalo, GridHndl grid, bool saveprofile=true);

#endif
