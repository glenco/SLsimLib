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

void map_images(LensHndl lens,Source *source,GridHndl grid,int *Nimages,ImageInfo *imageinfo
		,int Nimagesmax,double xmax,double xmin,double initial_size
		,ExitCriterion criterion,bool kappa_off,bool FindCenter,bool divide_images);

int refine_grid_on_image(Lens *lens,Source *source,GridHndl grid,ImageInfo *imageinfo,int *Nimages
		,ImageInfo *sourceinfo,int Nsources,int NimageMax,const double res_target,ExitCriterion criterion
		,bool kappa_off,bool divide_images,bool batch=true);
void check_sb_add(Source *source,ImageInfo *imageinfo,Point *i_points,unsigned long &Nold,int &number_of_refined);

bool RefinePoint(Point *point,TreeHndl i_tree,double image_area,double total_area,ExitCriterion criterion
		,double res_target,Kist<Point> * nearest);
bool RefinePoint2(Point *point,TreeHndl i_tree,double image_area,double total_area,ExitCriterion criterion
		,double res_target,Kist<Point> * nearest);
void UniformMagCheck(ImageInfo *imageinfo);

#endif
