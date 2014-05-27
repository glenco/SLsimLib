/*
 * map_images.h
 *
 *  Created on: Oct 6, 2010
 *      Author: bmetcalf
 */

#ifndef MAP_IMAGES_H_
#define MAP_IMAGES_H_

#include <lens.h>
#include <grid_maintenance.h>

namespace ImageFinding{

void map_images(LensHndl lens,Source *source,GridHndl grid,int *Nimages,ImageInfo *imageinfo
		,int Nimagesmax,double xmax,double xmin,double initial_size
		,ExitCriterion criterion,bool FindCenter,bool divide_images);
  
void map_images_fixedgrid(Source *source,GridHndl grid ,int *Nimages ,ImageInfo *imageinfo
                          ,int NimageMax ,double xmax ,bool divide_images,bool find_borders);

int refine_grid_on_image(Lens *lens,Source *source,GridHndl grid,double maxflux
                         ,ImageInfo *imageinfo,int *Nimages,ImageInfo *sourceinfo,int Nsources
                         ,int NimageMax,const double res_target,ExitCriterion criterion
                         ,bool divide_images,bool batch=true);
void check_sb_add(Source *source,ImageInfo *imageinfo,Point *i_points,double maxflux,unsigned long &Nold,int &number_of_refined);

bool RefinePoint(Point *point,TreeHndl i_tree,double image_area,double total_area,ExitCriterion criterion
		,double res_target,Kist<Point> * nearest);
bool RefinePoint2(Point *point,TreeHndl i_tree,double image_area,double total_area,double maxflux,ExitCriterion criterion
		,double res_target,Kist<Point> * nearest);
  
bool RefinePoint_sb(Point *point,TreeHndl i_tree,double image_area,double total_area
                      ,double sb_limit,PosType maxflux,Kist<Point> * nearest);

bool RefinePoint_smallsize(Point *point,TreeHndl i_tree,double image_area,double total_area
                           ,double smallsize,PosType maxflux,Kist<Point> * nearest);
bool IntegrateFluxInCell(Point *point,TreeHndl i_tree,Source *source);
void interpfrom2Points(Point const * p1,Point const * p2,PosType *x,PosType *y);
    
void UniformMagCheck(ImageInfo *imageinfo);
}
#endif
