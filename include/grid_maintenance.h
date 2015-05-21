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
#include <utilities_slsim.h>

class LensHaloBaseNSIE;
class LensHaloMassMap;

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
  unsigned long getNumberOfPoints() const;
  
  
  /// tree on image plane
  TreeHndl i_tree;
  /// tree on source plane
  TreeHndl s_tree;
  
  /// return initial number of grid points in each direction
  int getInitNgrid(){return Ngrid_init;}
  /// return number of cells in each dimension into which each cell is divided when a refinement is made
  int getNgrid_block(){return Ngrid_block;}
  /// return initial range of gridded region
  double getInitRange(){return i_tree->getTop()->boundary_p2[0] - i_tree->getTop()->boundary_p1[0];}
  Point_2d getInitCenter();
  Point * RefineLeaf(LensHndl lens,Point *point);
  Point * RefineLeaves(LensHndl lens,std::vector<Point *>& points);
  void ClearAllMarks();
  
  void test_mag_matrix();
  void writeFits(const double center[],size_t Npixels,double resolution,LensingVariable lensvar,std::string filename);
  void writeFits(const double center[],size_t Nx,size_t Ny,double resolution,LensingVariable lensvar,std::string filename);
  void writePixeFits(size_t Nx           /// number of pixels in image in x dimension
                    ,LensingVariable lensvar  /// which quantity is to be displayed
                    ,std::string filename     /// file name for image -- .kappa.fits, .gamma1.fits, etc will be appended
                    );
  void writeFitsVector(const double center[],size_t Npixels,double resolution,LensingVariable lensvar,std::string filename);
  PixelMap writePixelMap(const double center[],size_t Npixels,double resolution,LensingVariable lensvar);
  PixelMap writePixelMap(const double center[],size_t Nx,size_t Ny,double resolution,LensingVariable lensvar);
  
  PixelMap writePixelMapUniform(const PosType center[],size_t Nx,size_t Ny,LensingVariable lensvar);
  void writePixelMapUniform(PixelMap &map,LensingVariable lensvar);
  void writeFitsUniform(const PosType center[],size_t Nx,size_t Ny,LensingVariable lensvar,std::string filename);
  
private:
  void xygridpoints(Point *points,double range,const double *center,long Ngrid
                    ,short remove_center);
  
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

/// enumerates the types of critical curves. ND is "not defined".
enum CritType {ND,radial,tangential,pseudo};
// in image_finder_kist.c
namespace ImageFinding{
  
  
  
  struct CriticalCurve{
    
    CriticalCurve(){
      critical_center[0] = critical_center[1] = 0.0;
      caustic_center[0] = caustic_center[1] = 0.0;
      critical_area = 0.0;
      caustic_area = 0.0;
      contour_ell = 0.0;
      ellipse_area = 0.0;
      z_source = 0.0;
      type = ND;
      caustic_intersections = -1;
    };
    CriticalCurve(const CriticalCurve &p){
      //critical_curve.resize(p.critical_curve.size());
      critical_curve = p.critical_curve;
      caustic_curve_outline = p.caustic_curve_outline;
      caustic_curve_intersecting = p.caustic_curve_intersecting;
      critical_center = p.critical_center;
      caustic_center = p.caustic_center;
      critical_area = p.critical_area;
      caustic_area = p.caustic_area;
      ellipse_curve = p.ellipse_curve;
      contour_ell = p.contour_ell;
      ellipse_area = p.ellipse_area;
      z_source = p.z_source;
      type = p.type;
      caustic_intersections = p.caustic_intersections;
   }

    CriticalCurve & operator=(const CriticalCurve &p){
      if(this == &p) return *this;
      
      critical_curve = p.critical_curve;
      caustic_curve_outline = p.caustic_curve_outline;
      caustic_curve_intersecting = p.caustic_curve_intersecting;
      critical_center = p.critical_center;
      caustic_center = p.caustic_center;
      critical_area = p.critical_area;
      caustic_area = p.caustic_area;
      ellipse_curve = p.ellipse_curve;
      contour_ell = p.contour_ell;
      ellipse_area = p.ellipse_area;
      z_source = p.z_source;
      type = p.type;
      caustic_intersections = p.caustic_intersections;
      return *this;
    }
    
    std::vector<Point_2d> critical_curve;
    std::vector<Point_2d> caustic_curve_outline;
    std::vector<Point_2d> caustic_curve_intersecting;
    std::vector<Point_2d> ellipse_curve;
    
    PosType z_source;
    CritType type;
    int caustic_intersections;  /// estimated number of intersections of the caustic, -1 if not set

    Point_2d critical_center;      /// center of critical curve
    Point_2d caustic_center;   /// center of caustic curve
    
    PosType critical_area;        /// area of critical curve (radians^2)
    PosType caustic_area;        /// area of caustic curve (radians^2)
    
    
    PosType contour_ell;  /// axis ratio of a contour defined by the ratio of the max to min distance between center (as given by hull alg) and contour
    PosType ellipse_area;  /// area of an ellipse with axis ratio contour_ell and major axis = max distance between center (as given by hull alg) and contour
    
    /// return true if x is inside or on the border of the caustic curve
    bool inCausticCurve(Point_2d x){
      return Utilities::incurve(x.x,caustic_curve_outline);
    }
    
    /// return true if x is inside or on the border of the critical curve
    bool inCriticalCurve(Point_2d x){
      return Utilities::incurve(x.x,critical_curve);
    }

    /// return true if x is strictly inside (entirely) the caustic curve
    bool EntirelyinCausticCurve(Point_2d x, PosType sourceRadius)
    {
      // Testing if the center of the source is within the caustic :
      bool IsInCurve = Utilities::incurve(x.x,caustic_curve_outline);
      
      // Testing now that it is not too close from it (i.e. farther than source radius) :
      int i=0; // index going over the different points of the caustic curve
      PosType DistSourceToCautic; // distance between the source center and the considered point of the caustic line.
      
      if(IsInCurve == true) // center inside the caustic
      {
        while(i<caustic_curve_outline.size()) // testing that the source size does not overlap with the caustic line.
        {
          DistSourceToCautic = sqrt((caustic_curve_outline[i].x[0] - x.x[0])*(caustic_curve_outline[i].x[0] - x.x[0]) + (caustic_curve_outline[i].x[1] - x.x[1])*(caustic_curve_outline[i].x[1] - x.x[1]));
          if (DistSourceToCautic < sourceRadius) return false ; // source too close, we return false and don't consider the point.
          i++;
        }
        return true ; // if not we return true (the source is valid)
      }
      else return false ; // center not inside the caustic
    }
    
    /// return true if x is strictly inside (entirely) the critical curve
    bool EntirelyinCriticalCurve(Point_2d x, PosType sourceRadius)
    {
      // Testing if the center of the source is within the critical curve :
      bool IsInCurve = Utilities::incurve(x.x,critical_curve);
      
      // Testing now that it is not too close from it (i.e. farther than source radius) :
      int i=0; // index going over the different points of the critical curve
      PosType DistSourceToCritCurve; // distance between the source center and the considered point of the critical line.
      
      if(IsInCurve == true) // center inside the critical line
      {
        while(i<critical_curve.size()) // testing that the source size does not overlap with the critical line.
        {
          DistSourceToCritCurve = sqrt((critical_curve[i].x[0] - x.x[0])*(critical_curve[i].x[0] - x.x[0]) + (critical_curve[i].x[1] - x.x[1])*(critical_curve[i].x[1] - x.x[1]));
          if (DistSourceToCritCurve < sourceRadius) return false ; // source too close, we return false and don't consider the point.
          i++;
        }
        return true ; // if not we return true (the source is valid)
      }
      else return false ; // center not inside the critical line
    }
    
    
    /** \brief Returns a vector of random point within the caustic.  It is more efficient to call this once for many point rather than repeatedly one at a time.
     */
    void RandomSourceWithinCaustic(
                                   int N                   /// number of points needed
                                   ,std::vector<Point_2d> &y  /// output vector of points
                                   ,Utilities::RandomNumbers_NR &rng  /// random number generator
                                   );
    /** \brief Returns a vector of random point strictly within the caustic (i.e. not touching the border).  It is more efficient to call this once for many point rather than repeatedly one at a time.
     */
    void RandomSourceStrictlyWithinCaustic(int N                              /// number of points needed
                                           ,std::vector<Point_2d> &y          /// output vector of points
                                           ,Utilities::RandomNumbers_NR &rng  /// random number generator
                                           ,PosType sourceRadius              /// radius of the source
                                           ,PosType distSourceToCaustic       /// distance wanted between the source and the caustic line (must be larger than sourceRadius)
                                           );
    /// find rectangular region enclosing critical curve
    void CritRange(Point_2d &p1,Point_2d &p2);
    /// find rectangular region enclosing caustic curve
    void CausticRange(Point_2d &p1,Point_2d &p2);
    
    /// find 3 measures of the critical curve radius
    void CriticalRadius(PosType &rmax,PosType &rmin,PosType &rave) const{
      if(critical_curve.size() < 2){
        rave = rmax = rmin = 0.0;
        return;
      }
      
      rave = 0;
      rmax = 0;
      rmin = std::numeric_limits<double>::max();
      PosType rad;
      
      for(size_t ii=0; ii< critical_curve.size(); ++ii){
        Point_2d tmp = critical_center - critical_curve[ii];
        
        rad = tmp.length();
        
        rave += rad;
        if(rad < rmin) rmin = rad;
        if(rad > rmax) rmax = rad;
      }
      
      rave /= critical_curve.size();
    }
    /// find 3 measures of the caustic curve radius
    void CausticRadius(PosType &rmax,PosType &rmin,PosType &rave) const{
      if(caustic_curve_outline.size() < 2){
        rave = rmax = rmin = 0.0;
        return;
      }
      rave = 0;
      rmax = 0;
      rmin = std::numeric_limits<double>::max();
      PosType rad;
      
      for(size_t ii=0; ii< caustic_curve_outline.size(); ++ii){
        Point_2d tmp = caustic_center - caustic_curve_outline[ii];
        
        rad = tmp.length();
        
        rave += rad;
        if(rad < rmin) rmin = rad;
        if(rad > rmax) rmax = rad;
      }
      
      rave /= caustic_curve_outline.size();
    }

    
  private:
    Point_2d p1,p2;
  };
  
  void find_images_kist(LensHndl lens,PosType *y_source,PosType r_source,GridHndl grid
                        ,int *Nimages,std::vector<ImageInfo> &imageinfo,unsigned long *Nimagepoints
                        ,PosType initial_size,bool splitimages,short edge_refinement
                        ,bool verbose = false);
  
  void find_image_simple(LensHndl lens,Point_2d y_source,PosType z_source,Point_2d &image_x
                         ,PosType xtol2,PosType &fret);
  
  void find_images_microlens(LensHndl lens,double *y_source,double r_source,GridHndl grid
                             ,int *Nimages,std::vector<ImageInfo> &imageinfo,unsigned long *Nimagepoints
                             ,double initial_size,double mu_min,bool splitimages,short edge_refinement
                             ,bool verbose);
  
  void find_images_microlens_exper(LensHndl lens,PosType *y_source,PosType r_source
                                   ,GridHndl grid,int *Nimages,std::vector<ImageInfo> &imageinfo,unsigned long *Nimagepoints,PosType initial_size ,PosType mu_min
                                   ,bool splitimages,short edge_refinement,bool verbose);
  
  void image_finder_kist(LensHndl lens, PosType *y_source,PosType r_source,GridHndl grid
                         ,int *Nimages,std::vector<ImageInfo> &imageinfo,unsigned long *Nimagepoints
                         ,short splitparities,short true_images);
  
  
  void find_crit(LensHndl lens,GridHndl grid,std::vector<CriticalCurve> &crtcurve,int *Ncrits
                 ,double resolution,double invmag_min = 0.0,bool verbose = false,bool test=false);
 
  //void find_crit2(LensHndl lens,GridHndl grid,std::vector<CriticalCurve> &critcurve,int *Ncrits
  //                ,double resolution,bool *orderingsuccess,bool ordercurve,bool dividecurves,double invmag_min = 0.0,bool verbose = false);
  
  CritType find_pseudo(ImageInfo &pseudocurve,ImageInfo &negimage
                                 ,PosType pseudolimit,LensHndl lens,GridHndl grid
                                 ,PosType resolution,Kist<Point> &paritypoints);
  
  void find_contour(LensHndl lens,GridHndl grid,std::vector<CriticalCurve> &contour,int *Ncrits,PosType resolution,bool *orderingsuccess,bool ordercurve, bool dividecurves, double contour_value,LensingVariable contour_type,bool verbose = false);
  
  namespace Temporary{
    //PosType *y;
    //Lens * lens;
    //Point *point;
    
    PosType mindyfunc(PosType *x);
  }

  namespace IF_routines{
    int refine_grid_kist(LensHndl lens,GridHndl grid,ImageInfo *imageinfo
                       ,int Nimages,double res_target,short criterion
                       ,Kist<Point> * newpointkist = NULL,bool batch=true);
    

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

  }
}


void saveImage(LensHaloMassMap *mokahalo, GridHndl grid, bool saveprofile=true);


#endif
