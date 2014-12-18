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

/**  \brief The ImageFinding namespace is for functions related to finding and mapping images.
 */
namespace ImageFinding{
  
  struct CriticalCurve{

    std::vector<PosType [2]> crit_curve;
    std::vector<PosType [2]> caustic_curve;
    
    PosType critical_center[2];      /// center of critical curve
    PosType caustic_center[2];   /// center of caustic curve
    PosType critical_area;        /// area of critical curve (radians^2)
    PosType caustic_area;        /// area of caustic curve (radians^2)
  };
  
  void map_images(LensHndl lens,Source *source,GridHndl grid,int *Nimages
                  ,std::vector<ImageInfo> &imageinfo
                  ,double xmax,double xmin,double initial_size
                  ,ExitCriterion criterion,bool FindCenter,bool divide_images);
  
  void map_images_fixedgrid(Source *source,GridHndl grid ,int *Nimages
                            ,std::vector<ImageInfo> &imageinfo
                            ,double xmax ,bool divide_images,bool find_borders);
  
  int refine_grid_on_image(Lens *lens,Source *source,GridHndl grid,double maxflux
                           ,std::vector<ImageInfo> &imageinfo,int *Nimages
                           ,std::vector<ImageInfo> &sourceinfo,int Nsources
                           ,const double res_target,ExitCriterion criterion
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
  void IntegrateFluxInCell(Point *point,Source &source,float tolerance,Boo &outcome);
  void IntegrateCellsParallel(Kist<Point>::iterator it1
                              ,Kist<Point>::iterator it2,Source *source,PosType *area,size_t *count);
  
  void interpfrom2Points(Point const * p1,Point const * p2,PosType *x,PosType *y);
  
  void UniformMagCheck(ImageInfo *imageinfo);
  
  
  void map_imagesISOP(LensHndl lens,Source *source,GridHndl grid,int *Nimages
                      ,std::vector<ImageInfo> &imageinfo,double rmax
                      ,double res_min,double initial_size,ExitCriterion criterion
                      ,bool divide_images,bool int_on = true,bool verbos=false);
  
  
  int refine_grid_on_imageISOP(Lens *lens,Source *source,GridHndl grid
                               ,std::vector<ImageInfo> &imageinfo,int *Nimages,int Nsources
                               ,double res_target,double res_min
                               ,double res_source_area,ExitCriterion criterion
                               ,bool divide_images);
  
  
}
#endif
