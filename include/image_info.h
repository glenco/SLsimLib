/*
 * image_info.h
 *
 *  Created on: Jun 26, 2012
 *      Author: mpetkova
 */

#ifndef IMAGE_INFO_H_
#define IMAGE_INFO_H_

#include "Kist.h"
#include "pointlist.h"
#include "Tree.h"
//#include <grid_maintenance.h>

#ifndef gate_declare
#define gate_declare
typedef enum{unchecked,yes,no} GATE;
#endif

/** \brief Structure for storing information about images or curves */
struct ImageInfo{

	ImageInfo();
  ImageInfo(const ImageInfo & rhs){
    imagekist = new Kist<Point>;
    innerborder = new Kist<Point>;
    outerborder = new Kist<Point>;

    copy(rhs,true);
  }
	~ImageInfo();
  ImageInfo & operator=(const ImageInfo & rhs){
    if(this == &rhs) return *this;
    copy(rhs,true);
    return *this;
  }
    /// Array of points in image,  SHOULD NOT BE USED IN FAVOR OF imagekist!  Still used by caustic finding routines.
  //Point *points;
  /// Number of points in image, SHOULD NOT BE USED IN FAVOR OF imagekist->Nunits().  Still used by caustic finding routines.
  //unsigned long Npoints;
  /// later addition, holds all points in image, will replace points eventually
  Kist<Point> * imagekist;
  /// gridrange[2] minimum grid size in image, gridrange[0] maximum grid size in outerborder, gridrange[1] maximum grid size in image
  PosType gridrange[3];
  /// Centroid of image
  PosType centroid[2];
  /// area of image or, when using map_images(), the total brightness of the image
  PosType area;
  /// error on the estimate of area
  PosType area_error;
  /// the points on the inner border of the image
  Kist<Point> * innerborder;
  /// the points on the outer border of the image, i.e. not in the image
  Kist<Point> * outerborder;
  //short Nencircled;

  /// Flag for stopping refinement but also used for other temporary purposes.
  short ShouldNotRefine;
  /// Flag for showing when the distortion of an image can be considered linear.
  GATE uniform_mag;

  /// returns number of points currently in the image
  unsigned long getNimagePoints(){return imagekist->Nunits();}
  /// Computes the time delay averaged over the image
  KappaType aveTimeDelay();
  /// Computes the inverse magnification averaged over the image
  KappaType aveInvMag();
  /// Print information about the image
  void PrintImageInfo();
  void copy(const ImageInfo & image,bool copykists = true);
  
  bool IsMergedImages();

  void ArcInfo(PosType *area,PosType *area_circ,PosType theta);
  void FindArc(PosType &radius,PosType *xc,PosType *arc_c,PosType &arclength,PosType &width
            ,PosType resolution,PosType threshold);
  bool constant(LensingVariable lenvar,PosType tol);
  
  /** \brief Area of the image using the concave hull of the points on the image plane.
   
   If useborder is true it will use the inner border to save time.  When a bourder has 
   not been calculated,as with a critical cureve, it should be set to false.
   */
  PosType ConcaveHullImageArea(bool useborder);
  /** \brief Area of the image using the concave hull of the points on the image plane.
   
   If useborder is true it will use the inner border to save time.  When a bourder has
   not been calculated,as with a critical cureve, it should be set to false.
   */
  PosType ConcaveHullSourceArea(bool useborder);
};

/** \brief This is an old version that should not be used anymore in favor of ImageInfo.
 *
 * This type is still used in the critical curve finding and handling routines.  The results of these
 * routines should be returned in the form of and array of ImageInfo's now.
 */
typedef struct OldImageInfo{

	OldImageInfo();
	~OldImageInfo();

    /// Array of points in image,  SHOULD NOT BE USED IN FAVOR OF imagekist!  Still used by caustic finding routines.
  Point *points;
  /// Number of points in image, SHOULD NOT BE USED IN FAVOR OF imagekist->Nunits().  Still used by caustic finding routines.
  unsigned long Npoints;
  /// later addition, holds all points in image, will replace points eventually
  //Kist<Point> * imagekist;
  /// gridrange[2] minimum grid size in image, gridrange[0] maximum grid size in outerborder, gridrange[1] maximum grid size in image
  PosType gridrange[3];
  /// Centroid of image
  PosType centroid[2];
  /// area of image or, when using map_images(), the total brightness of the image
  PosType area;
  /// error on the estimate of area
  PosType area_error;
  /// the points on the inner border of the image
  Kist<Point> * innerborder;
  /// the points on the outer border of the image, i.e. not in the image
  Kist<Point> * outerborder;
  short ShouldNotRefine;

} OldImageInfo;


#endif /* IMAGE_INFO_H_ */
