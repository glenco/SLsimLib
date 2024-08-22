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
#include "source.h"
#include <mutex>
#include <utilities_slsim.h>
#include "concave_hull.h"

class LensHaloBaseNSIE;
class LensHaloMassMap;

/** 
 * \brief Structure to contain both source and image trees.
 */
struct Grid{
  
  Grid(LensHndl lens,unsigned long N1d,const double center[2],double range);
  Grid(LensHndl lens ,unsigned long Nx ,const PosType center[2] ,PosType rangeX ,PosType rangeY);
  ~Grid();
  
  Grid ReInitialize(LensHndl lens);
  //void ReShoot(LensHndl lens);
  void zoom(LensHndl lens,double *center,double scale,Branch *top = NULL);
  
  //unsigned long PruneTrees(double resolution,bool useSB,double fluxlimit);
  //unsigned long PrunePointsOutside(double resolution,double *y,double r_in ,double r_out);
  
  double RefreshSurfaceBrightnesses(Source* source);
  
  double AddSurfaceBrightnesses(Source* source);
  
  
  double mark_closest_point_source_images(
                                Point_2d y_source    /// angular position of source,
                               ,PosType r_source_max  /// points outside this radius on the source plane will not be considered as in the image
                               ,PosType luminosity
                               ,bool verbose=false
                         );
  void find_point_source_images(
                               Point_2d y_source    /// angular position of source,
                              ,PosType r_source  /// points outside this radius on the source plane will not be considered as in the image
                              ,PosType z_source  /// redhsift of source
                              ,std::vector<RAY> &images /// returned image rays
                              ,bool verbose=false
                           );
  
  double ClearSurfaceBrightnesses();
  
  /*** Refine the gris based on the smoothness of the surface brightness.
     Return new total flux.
   
   May be slow.
   */
  double refine_on_surfacebrightness(Lens &lens,Source &source);
  
  unsigned long getNumberOfPoints() const;
  /// area of region with negative magnification
  PosType EinsteinArea() const;
  
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
  
  //void test_mag_matrix();
  template <typename T>
  void writeFits(const double center[],size_t Npixels,double resolution,LensingVariable lensvar,std::string filename);
  template <typename T>
  void writeFits(const double center[],size_t Nx,size_t Ny,double resolution,LensingVariable lensvar,std::string filename);
  
  /// make a fits image of whole grid region
  template <typename T>
  void writeFits(
                 double strech  /// resolution relative to the initial resolution
                 ,LensingVariable lensvar /// which quantity is to be displayed
                 ,std::string filename    /// file name for image -- .kappa.fits, .gamma1.fits, etc will be appended
                 );
  template <typename T>
  void writePixelFits(size_t Nx           /// number of pixels in image in x dimension
                    ,LensingVariable lensvar  /// which quantity is to be displayed
                    ,std::string filename     /// file name for image -- .kappa.fits, .gamma1.fits, etc will be appended
                    );
  template <typename T>
  void writeFitsVector(const double center[],size_t Npixels,double resolution,LensingVariable lensvar,std::string filename);
  template <typename T>
  PixelMap<T> writePixelMap(const double center[],size_t Npixels,double resolution,LensingVariable lensvar);
  template <typename T>
  PixelMap<T> writePixelMap(const double center[],size_t Nx,size_t Ny,double resolution,LensingVariable lensvar);
  
  /// With the initial boundaries and resolution, ie no refinement
  template <typename T>
  PixelMap<T> writePixelMap(LensingVariable lensvar);
  
  /// make image of surface brightness
  template<typename T>
  void MapSurfaceBrightness(PixelMap<T> &map){
    map.Clean();
    map.AddGridBrightness(*this);
  }
  /// make a map of the whole gridded area with given resolution
  template <typename T>
  PixelMap<T> MapSurfaceBrightness(double resolution);

  template <typename T>
  PixelMap<T> writePixelMapUniform(const PosType center[],size_t Nx,size_t Ny,LensingVariable lensvar);
  template<typename T>
  void writePixelMapUniform(PixelMap<T> &map,LensingVariable lensvar);
  
  template <typename T>
  void writeFitsUniform(const PosType center[],size_t Nx,size_t Ny,LensingVariable lensvar,std::string filename);
  
  void find_images(
                   PosType *y_source
                   ,PosType r_source
                   ,int &Nimages
                   ,std::vector<ImageInfo> &imageinfo
                   ,unsigned long &Nimagepoints
                   );
 
  void map_images(
                        Lens* lens
                        ,Source *source
                        ,int *Nimages
                        ,std::vector<ImageInfo> &imageinfo
                        ,PosType xmax
                        ,PosType xmin
                        ,PosType initial_size
                        ,ExitCriterion criterion
                        ,bool FindCenter
                        ,bool divide_images
                        );
 
  Grid(Grid &&grid){
    *this = std::move(grid);
  }
  
  Grid operator=(Grid &grid) = delete;
  Grid(Grid &grid) =  delete;
  
  Grid & operator=(Grid &&grid){
    assert(&grid != this);
    
    i_tree = grid.i_tree;
    grid.i_tree = nullptr;
    s_tree = grid.s_tree;
    grid.s_tree = nullptr;
    neighbors = grid.neighbors;
    grid.neighbors = nullptr;
    //trashkist = grid.trashkist;
    //grid.trashkist = nullptr;

    Ngrid_init = grid.Ngrid_init;
    Ngrid_init2 = grid.Ngrid_init2;
    Ngrid_block = grid.Ngrid_block;
    initialized = grid.initialized;
    maglimit = grid.maglimit;
    pointID = grid.pointID;
    axisratio = grid.axisratio;
    
    point_factory.clear();
    point_factory=std::move(grid.point_factory);
    
    return *this;
  }
  
  /// flux weighted local magnification that does not take multiple imaging into effect
  PosType magnification(double sblimit=-1.0e12) const;
  PosType UnlensedFlux(double sblimit=-1.0e12) const;
  PosType LensedFlux(double sblimit=-1.0e12) const;
  
  //PosType magnification2() const;
  //PosType magnification3() const;
 /// centroid of flux
  Point_2d centroid() const;
  
  private:
  void xygridpoints(Point *points,double range,const double *center,long Ngrid
                    ,short remove_center);
  
  /// one dimensional size of initial grid
  int Ngrid_init;
  int Ngrid_init2;
  
  /// one dimensional number of cells a cell will be divided into on each refinement step
  int Ngrid_block;
  bool initialized;
  //Kist<Point> * trashkist;
  
  double maglimit;
  Kist<Point> * neighbors;
  bool find_mag_matrix(double *a,Point *p0,Point *p1,Point *p2);
  
  bool uniform_mag_from_deflect(double *a,Point *point);
  //bool uniform_mag_from_shooter(double *a,Point *point);
  
  double mag_from_deflect(Point *point) const;
  
  unsigned long pointID;
  PosType axisratio;

  template <typename T>
  void writePixelMapUniform_(Point *head,size_t N
                             ,PixelMap<T> *map,LensingVariable val);

  // cluge to make compatible with old method of producing points
  Point * NewPointArray(size_t N){
    Point * p = point_factory(N);
    p[0].head = N;
    for(size_t i=1; i < N ; ++i) p[i].head = 0;
    return p;
  }
  MemmoryBank<Point> point_factory;
  static std::mutex grid_mutex;
};

typedef struct Grid* GridHndl;

/// converts CritType into descriptive string
std::string to_string(CritType crit);
 
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
      type = CritType::ND;
      caustic_intersections = -1;
      touches_edge = false;
    };
    
    CriticalCurve(const CriticalCurve &p){
      critcurve = p.critcurve;
      //critical_curve = p.critical_curve;
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
      touches_edge = p.touches_edge;
   }

    CriticalCurve & operator=(const CriticalCurve &p){
      if(this == &p) return *this;
      
      critcurve = p.critcurve;
      //critical_curve = p.critical_curve;
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
      touches_edge = p.touches_edge;
      return *this;
    }
    
    std::vector<RAY> critcurve;
    std::vector<Point_2d> caustic_curve_outline;
    std::vector<Point_2d> caustic_curve_intersecting;
    std::vector<Point_2d> ellipse_curve;
    
    PosType z_source;
    /// type of caustic, 0 -not defined, 1 -radial, 2 - tangential,3 - pseudo
    CritType type;
      /// estimated number of intersections of the caustic, -1 if not set
    int caustic_intersections;

    /// center of critical curve
    Point_2d critical_center;
    /// center of caustic curve
    Point_2d caustic_center;
    
    /// area of critical curve (radians^2)
    PosType critical_area;
    /// area of caustic curve (radians^2)
    PosType caustic_area;
    
      /// axis ratio of a contour defined by the ratio of the max to min distance between center (as given by hull alg) and contour
    PosType contour_ell;
      /// area of an ellipse with axis ratio contour_ell and major axis = max distance between center (as given by hull alg) and contour
    PosType ellipse_area;
    
    /// touches the edge of the gridded region
    bool touches_edge;
    
    /// return true if x is inside or on the border of the caustic curve

    bool inCausticCurve(Point_2d &x){
      return Utilities::inhull(x.x,caustic_curve_outline);
    }
    
    /// returns true if a circle of radius r around the point x intersects with the caustic curve
    bool intersectingCausticCurve(Point_2d &x,double r){
      return Utilities::circleIntersetsCurve(x, r, caustic_curve_outline);
    }
 
    
    /// return true if x is inside or on the border of the critical curve
    bool inCriticalCurve(Point_2d &x){
      return Utilities::inhull<RAY>(x.x,critcurve);
    }

    /// return true if x is strictly inside (entirely) the caustic curve
    bool EntirelyinCausticCurve(Point_2d &x, PosType sourceRadius)
    {
      // Testing if the center of the source is within the caustic :
      bool IsInCurve = Utilities::inCurve(x,caustic_curve_outline);
      
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
      bool IsInCurve = Utilities::inCurve(x,caustic_curve_outline);
      
      // Testing now that it is not too close from it (i.e. farther than source radius) :
      int i=0; // index going over the different points of the critical curve
      PosType DistSourceToCritCurve; // distance between the source center and the considered point of the critical line.
      
      if(IsInCurve == true) // center inside the critical line
      {
        while(i<critcurve.size()) // testing that the source size does not overlap with the critical line.
        {
          DistSourceToCritCurve = sqrt((critcurve[i].x[0] - x.x[0])*(critcurve[i].x[0] - x.x[0]) + (critcurve[i].x[1] - x.x[1])*(critcurve[i].x[1] - x.x[1]));
          if (DistSourceToCritCurve < sourceRadius) return false ; // source too close, we return false and don't consider the point.
          i++;
        }
        return true ; // if not we return true (the source is valid)
      }
      else return false ; // center not inside the critical line
    }
    
    
    /** \brief Returns a vector of random point within the caustic.  It is more efficient to call this once for many point rather than repeatedly one at a time.
     */
    void RandomSourcesWithinCaustic(
                                   int N                   /// number of points needed
                                   ,std::vector<Point_2d> &y  /// output vector of points
                                   ,Utilities::RandomNumbers_NR &rng  /// random number generator
                                   );
    
    Point_2d RandomSourceWithinCaustic(
                                   Utilities::RandomNumbers_NR &rng  /// random number generator
                                   );

    /** \brief Returns a vector of random point within distance R or filly within the caustic.
     */
    void RandomSourcesNearCaustic(double R
                                   ,int N                   /// number of points needed
                                   ,std::vector<Point_2d> &y  /// output vector of points
                                   ,Utilities::RandomNumbers_NR &rng  /// random number generator
                                   );
    /** \brief Returns a random point within distance R or filly within the caustic.
     */
    Point_2d RandomSourceNearCaustic(double R
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
      if(critcurve.size() < 2){
        rave = rmax = rmin = 0.0;
        return;
      }
      
      rave = rmin = rmax = (critical_center - critcurve[0].x).length();
      PosType rad;
      
      for(size_t ii=1; ii< critcurve.size(); ++ii){
        rad = (critical_center - critcurve[ii].x).length();
        
        rave += rad;
        if(rad < rmin) rmin = rad;
        if(rad > rmax) rmax = rad;
      }
      
      rave /= critcurve.size();
    }
    /// find 3 measures of the caustic curve radius
    void CausticRadius(PosType &rmax,PosType &rmin,PosType &rave) const{
      if(caustic_curve_outline.size() < 2){
        rave = rmax = rmin = 0.0;
        return;
      }
      
      rave = rmin = rmax = (caustic_center - caustic_curve_outline[0]).length();

      PosType rad;
      
      for(size_t ii=1; ii< caustic_curve_outline.size(); ++ii){
        rad = (caustic_center - caustic_curve_outline[ii]).length();
        
        rave += rad;
        if(rad < rmin) rmin = rad;
        if(rad > rmax) rmax = rad;
      }
      
      rave /= caustic_curve_outline.size();
    }

    /// returns an estimate of the area inside and within distance R of the caustic
    double AreaNearCaustic(double R /// distance in radians
                           );
      
  private:
    Point_2d p1,p2;
  };
  
  void find_images_kist(LensHndl lens,PosType *y_source,PosType r_source,GridHndl grid
                        ,int *Nimages,std::vector<ImageInfo> &imageinfo,unsigned long *Nimagepoints
                        ,PosType initial_size,bool splitimages,short edge_refinement
                        ,bool verbose = false);
  
  //void find_image_simple(LensHndl lens,Point_2d y_source,PosType z_source,Point_2d &image_x
  //                       ,PosType xtol2,PosType &fret);
  
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
  void find_crit(Lens &lens,GridMap &gridmap,std::vector<CriticalCurve> &crtcurves,bool verbose = false);
  
  // finds the contours of magnification and source plane curve
  void find_magnification_contour(
    Lens &lens
    ,GridMap &gridmap
    ,double invmag
    ,std::vector<std::vector<RAY> > &contour
    ,std::vector<bool> &hits_boundary
  );

  //void find_crit2(LensHndl lens,GridHndl grid,std::vector<CriticalCurve> &critcurve,int *Ncrits
  //                ,double resolution,bool *orderingsuccess,bool ordercurve,bool dividecurves,double invmag_min = 0.0,bool verbose = false);
  
  CritType find_pseudo(ImageInfo &pseudocurve,ImageInfo &negimage
                                 ,PosType pseudolimit,LensHndl lens,GridHndl grid
                                 ,PosType resolution,Kist<Point> &paritypoints,bool TEST=false);
  
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

  void printCriticalCurves(std::string filename
                           ,const std::vector<ImageFinding::CriticalCurve> &critcurves);
  
  /** \brief Makes an image of the critical curves.  The map will encompose all curves found.  The
   pixel values are the caustic type + 1 ( 2=radial,3=tangential,4=pseudo )
   */
  template<typename T>
  PixelMap<T> mapCriticalCurves(
                             /// list of critical curves
                             const std::vector<ImageFinding::CriticalCurve> &critcurves,
                             /// number of pixels to each size
                             int Nx
                             );
  
  /** \brief Makes an image of the caustic curves.  The map will encompose all curves found.  The
   pixel values are the caustic type + 1 ( 2=radial,3=tangential,4=pseudo )
   */
  template<typename T>
  PixelMap<T> mapCausticCurves(const std::vector<ImageFinding::CriticalCurve> &critcurves /// list of critical curves
                            ,int Nx /// number of pixels to each size
                            );
}


std::ostream &operator<<(std::ostream &os, const ImageFinding::CriticalCurve &p);

void saveImage(LensHaloMassMap *mokahalo, GridHndl grid, bool saveprofile=true);


/// Outputs a fits image of a lensing variable of choice
template <typename T>
void Grid::writeFits(
                     const PosType center[]     /// center of image
                     ,size_t Npixels           /// number of pixels in image in on dimension
                     ,PosType resolution        /// resolution of image in radians
                     ,LensingVariable lensvar  /// which quantity is to be displayed
                     ,std::string filename     /// file name for image -- .kappa.fits, .gamma1.fits, etc will be appended
                     ){
  writeFits<T>(center,Npixels,Npixels,resolution,lensvar,filename);
}
  /// Outputs a fits image of a lensing variable of choice
template <typename T>
void Grid::writeFits(
                     const PosType center[]     /// center of image
                     ,size_t Nx           /// number of pixels in image in x dimension
                     ,size_t Ny           /// number of pixels in image in y dimension
                     ,PosType resolution        /// resolution of image in radians
                     ,LensingVariable lensvar  /// which quantity is to be displayed
                     ,std::string filename     /// file name for image -- .kappa.fits, .gamma1.fits, etc will be appended
                     ){
  PixelMap<T> map(center, Nx,Ny, resolution);
  std::string tag;
  
  switch (lensvar) {
    case LensingVariable::DELAYT:
      tag = ".dt.fits";
      break;
    case LensingVariable::ALPHA1:
      tag = ".alpha1.fits";
      break;
    case LensingVariable::ALPHA2:
      tag = ".alpha2.fits";
      break;
    case LensingVariable::ALPHA:
      tag = ".alpha.fits";
      break;
    case LensingVariable::KAPPA:
      tag = ".kappa.fits";
      break;
    case LensingVariable::GAMMA1:
      tag = ".gamma1.fits";
      break;
    case LensingVariable::GAMMA2:
      tag = ".gamma2.fits";
      break;
    case LensingVariable::GAMMA3:
      tag = ".gamma3.fits";
      break;
    case LensingVariable::GAMMA:
      tag = ".gamma.fits";
      break;
    case LensingVariable::INVMAG:
      tag = ".invmag.fits";
      break;
    case LensingVariable::SurfBrightness:
      tag = ".surfbright.fits";
      break;
    default:
      break;
  }
  
  map.AddGrid(*this,lensvar);
  map.printFITS(filename + tag);

  return;
}

/// Outputs a PixelMap of the lensing quantities of a fixed grid
template <typename T>
PixelMap<T> Grid::writePixelMap(
                             const PosType center[]     /// center of image
                             ,size_t Npixels           /// number of pixels in image in on dimension
                             ,PosType resolution        /// resolution of image in radians
                             ,LensingVariable lensvar  /// which quantity is to be displayed
                             ){
  
  return writePixelMap<T>(center,Npixels,Npixels,resolution,lensvar);
}
/// Outputs a PixelMap of the lensing quantities of a fixed grid
template <typename T>
PixelMap<T> Grid::writePixelMap(
                             const PosType center[]     /// center of image
                             ,size_t Nx           /// number of pixels in image in on dimension
                             ,size_t Ny           /// number of pixels in image in on dimension
                             ,PosType resolution        /// resolution of image in radians
                             ,LensingVariable lensvar  /// which quantity is to be displayed
){
  PixelMap<T> map(center, Nx, Ny, resolution);
  map.AddGrid(*this,lensvar);
  
  return map;
}
/// Outputs a PixelMap of the lensing quantities of a fixed grid
template <typename T>
PixelMap<T> Grid::writePixelMap(
                             LensingVariable lensvar  /// which quantity is to be displayed
){
  
  Branch *branch = i_tree->getTop();
  double resolution = (branch->boundary_p2[0] - branch->boundary_p1[0])/Ngrid_init;
  PixelMap<T> map(branch->center, Ngrid_init, Ngrid_init2, resolution);
  map.AddGrid(*this,lensvar);
  
  return map;
}

template <typename T>
PixelMap<T>  Grid::MapSurfaceBrightness(double resolution){
  Branch *branch = i_tree->getTop();
  int Nx = (int)( (branch->boundary_p2[0] - branch->boundary_p1[0])/resolution );
  int Ny = (int)( (branch->boundary_p2[1] - branch->boundary_p1[1])/resolution );
  
  PixelMap<T> map(branch->center,Nx,Ny,resolution);
  map.AddGridBrightness(*this);

  return map;
}

/** \brief Make a fits map that is automatically centered on the grid and has approximately the same range as the grid.  Nx can be used to change the resolution.  Nx = grid.getInitNgrid() will give the initial grid resolution
 */
template <typename T>
void Grid::writePixelFits(
                         size_t Nx           /// number of pixels in image in x dimension
                         ,LensingVariable lensvar  /// which quantity is to be displayed
                         ,std::string filename     /// file name for image -- .kappa.fits, .gamma1.fits, etc will be appended
){
  
  Point_2d center = getInitCenter();
  PosType resolution =  getInitRange()/Nx;
  size_t Ny = (size_t)(Nx*axisratio);
  writeFits<T>(center.x,Nx,Ny,resolution, lensvar, filename);
  
  return;
}

/** \brief Output a fits map of the without distribution the pixels.
 *
 *  This will be faster than Grid::writePixelMap() and Grid::writeFits().
 *  But it puts each grid pixel in one pixelmap pixel and if there are two
 *  grid pixels in one pixelmap pixel it uses one at random.  This is meant
 *  for uniform maps to make equal sized PixelMaps.
 */
template <typename T>
void Grid::writeFitsUniform(
                                const PosType center[]  /// center of image
                                ,size_t Nx       /// number of pixels in image in on dimension
                                ,size_t Ny       /// number of pixels in image in on dimension
                                ,LensingVariable lensvar  /// which quantity is to be displayed
                                ,std::string filename     /// file name for image -- .kappa.fits, .gamma1.fits, etc will be appended
                                ){
  std::string tag;
  
  switch (lensvar) {
    case LensingVariable::DELAYT:
      tag = ".dt.fits";
      break;
    case LensingVariable::ALPHA1:
      tag = ".alpha1.fits";
      break;
    case LensingVariable::ALPHA2:
      tag = ".alpha2.fits";
      break;
    case LensingVariable::ALPHA:
      tag = ".alpha.fits";
      break;
    case LensingVariable::KAPPA:
      tag = ".kappa.fits";
      break;
    case LensingVariable::GAMMA1:
      tag = ".gamma1.fits";
      break;
    case LensingVariable::GAMMA2:
      tag = ".gamma2.fits";
      break;
    case LensingVariable::GAMMA3:
      tag = ".gamma3.fits";
      break;
    case LensingVariable::GAMMA:
      tag = ".gamma.fits";
      break;
    case LensingVariable::INVMAG:
      tag = ".invmag.fits";
      break;
    case LensingVariable::SurfBrightness:
      tag = ".surfbright.fits";
      break;
    default:
      break;
  }

  PixelMap<T> map = writePixelMapUniform<T>(center,Nx,Ny,lensvar);
  map.printFITS(filename + tag);
}

/** \brief Make a Pixel map of the without distribution the pixels.
 *
 *  This will be faster than Grid::writePixelMap() and Grid::writeFits().
 *  But it puts each grid pixel in one pixelmap pixel and if there are two
 *  grid pixels in one pixelmap pixel it uses one at random.  This is meant
 *  for uniform maps to make equal sized PixelMaps.
 */
template <typename T>
PixelMap<T> Grid::writePixelMapUniform(
                                    const PosType center[]  /// center of image
                                    ,size_t Nx       /// number of pixels in image in on dimension
                                    ,size_t Ny       /// number of pixels in image in on dimension
                                    ,LensingVariable lensvar  /// which quantity is to be displayed
                                    ){
  
  if(getNumberOfPoints() ==0 ) return PixelMap<T>();
  PixelMap<T> map(center, Nx, Ny,i_tree->pointlist.Top()->gridsize);
  map.Clean();
  
  //int Nblocks = Utilities::GetNThreads();
  int Nblocks = 16;
                                      
  //std::vector<PointList> lists(Nblocks);
                                      
  std::vector<Point *> heads(Nblocks);
  std::vector<size_t> sizes(Nblocks,0);
  
  bool allowDecent;
  TreeStruct::iterator i_tree_it(i_tree);
  int i = 0;
                                      
  do{
    if((*i_tree_it)->level == 4){
      
      heads[i] = (*i_tree_it)->points;
      sizes[i] = (*i_tree_it)->npoints;
      
      //lists[i].setTop( (*i_tree_it)->points );
      //lists[i].setN( (*i_tree_it)->npoints );
      ++i;
      allowDecent = false;
    }else{
      allowDecent = true;
    }
  }while(i_tree_it.TreeWalkStep(allowDecent) && i < Nblocks);
  
  std::vector<std::thread> thrs;
  for(int ii = 0; ii < i ;++ii){
  //writePixelMapUniform_(heads[ii],sizes[ii],&map,lensvar);
  //thrs.push_back(std::thread(&Grid::writePixelMapUniform_,this,lists[ii],&map,lensvar));

    thrs.push_back(std::thread(&Grid::writePixelMapUniform_<T>,this,heads[ii],sizes[ii],&map,lensvar));
  }
  for(int ii = 0; ii < i ;++ii) thrs[ii].join();
  
  return map;
}

template <typename T>
void Grid::writePixelMapUniform(
                                    PixelMap<T> &map
                                    ,LensingVariable lensvar  /// which quantity is to be displayed
                                    ){
  
  if(getNumberOfPoints() ==0 ) return;
  
  map.Clean();
  int Nblocks = 16;
  //std::vector<PointList> lists(Nblocks);
  TreeStruct::iterator i_tree_it(i_tree);

  std::vector<Point *> heads(Nblocks);
  std::vector<size_t> sizes(Nblocks,0);

  
  bool allowDecent;
  i_tree_it.movetop();
  int i = 0;
  do{
    if((*i_tree_it)->level == 4){
      assert(i < 16);
      //lists[i].setTop( (*i_tree_it)->points );
      //lists[i].setN( (*i_tree_it)->npoints );
      heads[i] = (*i_tree_it)->points;
      sizes[i] = (*i_tree_it)->npoints;
      
      ++i;
      allowDecent = false;
    }else{
      allowDecent = true;
    }
  }while(i_tree_it.TreeWalkStep(allowDecent) && i < Nblocks);
  
  std::vector<std::thread> thr;
  for(int ii = 0; ii < i ;++ii){
    thr.push_back(std::thread(&Grid::writePixelMapUniform_<T>,this,heads[ii],sizes[ii],&map,lensvar));
  }
  for(auto &t : thr) t.join();
}

template <typename T>
void Grid::writePixelMapUniform_(Point *head,size_t N,PixelMap<T> *map,LensingVariable val){
  double tmp;
  PosType tmp2[2];
  long index;
  
  Point *ppoint = head;
  
  for(size_t i = 0; i< N; ++i){
    
    switch (val) {
      case LensingVariable::ALPHA:
        tmp2[0] = ppoint->x[0] - ppoint->image->x[0];
        tmp2[1] = ppoint->x[1] - ppoint->image->x[1];
        tmp = sqrt(tmp2[0]*tmp2[0] + tmp2[1]*tmp2[1]);
        break;
      case LensingVariable::ALPHA1:
        tmp = (ppoint->x[0] - ppoint->image->x[0]);
        break;
      case LensingVariable::ALPHA2:
        tmp = (ppoint->x[1] - ppoint->image->x[1]);
        break;
      case LensingVariable::KAPPA:
        tmp = ppoint->kappa();
        break;
      case LensingVariable::GAMMA:
        tmp2[0] = ppoint->gamma1();
        tmp2[1] = ppoint->gamma2();
        tmp = sqrt(tmp2[0]*tmp2[0] + tmp2[1]*tmp2[1]);
        break;
      case LensingVariable::GAMMA1:
        tmp = ppoint->gamma1();
        break;
      case LensingVariable::GAMMA2:
        tmp = ppoint->gamma2();
        break;
      case LensingVariable::GAMMA3:
        tmp = ppoint->gamma3();
        break;
      case LensingVariable::INVMAG:
        tmp = ppoint->invmag();
        break;
      case LensingVariable::DELAYT:
        tmp = ppoint->dt;
        break;
      case LensingVariable::SurfBrightness:
        tmp = ppoint->surface_brightness;
        break;
      default:
        std::cerr << "PixelMap<T>::AddGrid() does not work for the input LensingVariable" << std::endl;
        throw std::runtime_error("PixelMap<T>::AddGrid() does not work for the input LensingVariable");
        break;
        // If this list is to be expanded to include ALPHA or GAMMA take care to add them as vectors
    }
    
    index = map->find_index(ppoint->x);
    if(index != -1)(*map)[index] = tmp;
    
    ppoint = ppoint->next;
  }
}

/// Outputs a fits file for making plots of vector fields
template <typename T>
void Grid::writeFitsVector(
                     const PosType center[]     /// center of image
                     ,size_t Npixels           /// number of pixels in image in on dimension
                     ,PosType resolution        /// resolution of image in radians
                     ,LensingVariable lensvar  /// which quantity is to be displayed
                     ,std::string filename     /// file name for image -- .kappa.fits, .gamma1.fits, etc will be appended
                     ){
  //throw std::runtime_error("Not done yet!");
  
  PosType range = Npixels*resolution,tmp_x[2];
  ImageInfo tmp_image,tmp_image_theta;
  size_t i;
  std::string tag;
  
  i_tree->PointsWithinKist(center,range/sqrt(2.),tmp_image.imagekist,0);
  i_tree->PointsWithinKist(center,range/sqrt(2.),tmp_image_theta.imagekist,0);
  
  std::vector<PosType> tmp_sb_vec(tmp_image.imagekist->Nunits());
  
  for(tmp_image.imagekist->MoveToTop(),i=0;i<tmp_sb_vec.size();++i,tmp_image.imagekist->Down()){
    tmp_sb_vec[i] = tmp_image.imagekist->getCurrent()->surface_brightness;
    switch (lensvar) {
      case LensingVariable::ALPHA1:
        tmp_x[0] = tmp_image.imagekist->getCurrent()->x[0]
            - tmp_image.imagekist->getCurrent()->image->x[0];

        tmp_x[1] = tmp_image.imagekist->getCurrent()->x[1]
            - tmp_image.imagekist->getCurrent()->image->x[1];
      
        tmp_image.imagekist->getCurrent()->surface_brightness = sqrt( tmp_x[0]*tmp_x[0] + tmp_x[1]*tmp_x[1]);
        tmp_image_theta.imagekist->getCurrent()->surface_brightness = atan2(tmp_x[1],tmp_x[0]);
            
        tag = ".alphaV.fits";
        break;
      case LensingVariable::GAMMA:
        
        tmp_x[0] = tmp_image.imagekist->getCurrent()->gamma1();
        tmp_x[1] = tmp_image.imagekist->getCurrent()->gamma2();

        tmp_image.imagekist->getCurrent()->surface_brightness = sqrt( tmp_x[0]*tmp_x[0] + tmp_x[1]*tmp_x[1]);
        tmp_image_theta.imagekist->getCurrent()->surface_brightness = atan2(tmp_x[1],tmp_x[0])/2;
            
        tag = ".gammaV.fits";
        break;
      default:
        std::cout << "Grid::writeFitsVector() does not support the LensVariable you are using." << std::endl;
        return;
    }
  }
  
  PixelMap<T> map_m(center, Npixels, resolution),map_t(center,Npixels,resolution);
  
  map_m.Clean();
  map_m.AddImages(&tmp_image,1,-1);
  map_m = PixelMap<T>(map_m,4);
  map_m = PixelMap<T>(map_m,1/4.);
    
  map_t.Clean();
  map_t.AddImages(&tmp_image_theta,1,-1);

  map_m.printFITS(filename + tag);
  
  for(tmp_image.imagekist->MoveToTop(),i=0;i<tmp_sb_vec.size();++i,tmp_image.imagekist->Down())
    tmp_image.imagekist->getCurrent()->surface_brightness = tmp_sb_vec[i];
}

template <typename T>
void Grid::writeFits(double strech,LensingVariable lensvar ,std::string filename){
  Point_2d center = getInitCenter();
  size_t N1 = (size_t)(strech*Ngrid_init);
  size_t N2 = (size_t)(strech*Ngrid_init2);
  
  writeFits<T>(center.x,N1,N2,getInitRange()/N1,lensvar,filename);
}



#endif
