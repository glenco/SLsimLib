//
//  gridmap.h
//  GLAMER
//
//  Created by bmetcalf on 8/19/14.
//
//

#ifndef __GLAMER__gridmap__
#define __GLAMER__gridmap__

#include <iostream>
#include "lens.h"
#include "point.h"
#include "concave_hull.h"
#include "Tree.h"
#include "source.h"
#include <mutex>

/** 
 * \brief A simplified version of the Grid structure for making non-adaptive maps of the lensing quantities (kappa, gamma, etc...)
 *
 *  GripMap is faster and uses less memory than Grid.  It does not construct the tree structures for the points 
 *  and thus cannot be used for adaptive mapping or image finding.
 *
 *  The distance between the left (lower) most and right (upper) most ray is range so the resolution is range/(N-1).
 *  The lower left pixel is at center[]-0.5*range and the upper right is at center[]+0.5*range
 */

struct GridMap{
  
  friend class Lens;
  
	GridMap(LensHndl lens,unsigned long N1d,const double center[2],double range);
  GridMap(LensHndl lens ,unsigned long Nx ,const PosType center[2] ,PosType rangeX ,PosType rangeY);
  /// this makes a dumy GridMap that has no lensing
  GridMap(unsigned long N1d,const double center[2],double range);
	~GridMap();
  
    /// reshoot the rays for example when the source plane has been changed
  GridMap ReInitialize(LensHndl lens);
  
  /// resets to state without lensing
  void deLens();
  /**
   * \brief Recalculate surface brightness at every point without changing the positions of the gridmap or any lens properties.
   *
   *  Recalculate the surface brightness at all points on the gridmap.
   * This is useful when changing the source model while preserving changes in the grid.
   * Both i_tree and s_tree are both changed although only s_tree shows up here.
   *
   * returns the sum of the surface brightnesses
   */
  double RefreshSurfaceBrightnesses(SourceHndl source);
  /**
   Oversample some pixels where the usrface brightness is not smooth and update surface brighnesses to be the average inside the pixel.
   
   May be slow.
   */
  double AdaptiveRefreshSurfaceBrightnesses(Lens &lens,Source &source);
  
 /**
   * \brief Recalculate surface brightness just like GridMap::RefreshSurfaceBrightness but
   * the new source is added to any sources that were already there.
   *
   * returns total flux from the new source
   */
  double AddSurfaceBrightnesses(SourceHndl source);

  /// get the image point for a index number
  Point_2d image_point(size_t index){return i_points[index];}
  /// get the image point for a index number
  Point_2d source_point(size_t index){return s_points[index];}
 
  void ClearSurfaceBrightnesses();
  void assertNAN(); // check for nan in surface prightness
	size_t getNumberOfPoints() const {return Ngrid_init*Ngrid_init2;}
  
	/// return initial number of grid points in each direction
	int getInitNgrid() const {return Ngrid_init;}
	/// return initial range of gridded region.  This is the distance from the first ray in a row to the last (unlike PixelMap)
	double getXRange() const {return x_range;}
	double getYRange() const {return x_range*axisratio;}
  /// resolution in radians, this is range / (N-1)
  double getResolution() const {return x_range/(Ngrid_init-1);}
  
   /// make pixel map of lensing quantities at the resolution of the GridMap
  PixelMap writePixelMap(LensingVariable lensvar);
   /// fits output of lensing quantities at the resolution of the GridMap
  void writeFits(LensingVariable lensvar,std::string filensame);

  void writePixelMapUniform(PixelMap &map,LensingVariable lensvar);
  void writeFitsUniform(const PosType center[],size_t Nx,size_t Ny,LensingVariable lensvar,std::string filename);
  PixelMap writePixelMapUniform(const PosType center[],size_t Nx,size_t Ny,LensingVariable lensvar);

  /// this will make a fits map of the grid as is.
  void writeFitsUniform(
                        LensingVariable lensvar    ///< quantity to be output
                        ,std::string filename     ///< name of output fits file
                        ){
    PixelMap map = writePixelMap(lensvar);
    map.printFITS(filename);
  }
  
  /// returns a PixelMap with the flux in pixels at a resolution of res times the original resolution
  PixelMap getPixelMapFlux(int res) const;
  /// update a PixelMap with the flux in pixels at a resolution of res times the original resolution.
  /// The map must have precisely the right size and center to match or an exception will be thrown.
  /// Constructing the map with PixelMap getPixelMapFlux(int res) will insure that it does.
  void getPixelMapFlux(PixelMap &map) const;
  
  /// returns the area (radians^2) of the region with negative magnification at resolution of fixed grid
  PosType EinsteinArea() const;

  /** flux weighted local magnification with current surface brightness averaged on the image plane,
   */
  //PosType magnification_invmag() const;
  //PosType magnification2() const;
  
  /// returns centroid of flux on the grid
  Point_2d centroid() const;
  
  Point_2d getCenter(){return center;}
  
  Point * operator[](size_t i){return i_points + i;};
  
  GridMap(GridMap &&grid){
    *this = std::move(grid);
  }
  
  GridMap & operator=(GridMap &&grid){
    Ngrid_init = grid.Ngrid_init;
    Ngrid_init2 = grid.Ngrid_init2;
    pointID = grid.pointID;
    axisratio = grid.axisratio;
    x_range = grid.x_range;
    
    i_points = grid.i_points;
    grid.i_points = nullptr;
    s_points = grid.s_points;
    grid.s_points = nullptr;
    
    center = grid.center;
    
    point_factory = std::move(grid.point_factory);

    return *this;
  }

  struct Triangle{
    Triangle(size_t i,size_t j,size_t k){
      index[0] = i;
      index[1] = j;
      index[2] = k;
    }
    size_t index[3];
    size_t & operator[](int i){return index[i];}
  };
  
    
  /*** \brief Returns a list of  RAYs from a set of source positions.
   
   The image positions are found in parallel by the triangle method.  The order of the
   output rays will be the same as the sources with multiple images consecutive.
   */
  std::list<RAY> find_images(std::vector<Point_2d> &ys
                             ,std::vector<int> &multiplicity
                             ) const;

  /** find all images by triangle method
   */
  void find_images(Point_2d y
                   ,std::vector<Point_2d> &image_points  /// positions of the images limited by resolution of the gridmap
                   ,std::vector<Triangle> &triangles     /// index's of the points that form the triangles that the images are in
  ) const ;
  
  /** \brief finds the boundary of the region on the source plane where there are more than one image
   
   This uses the triangle method to determin which points in a source plane grid of the same size and resolution as the image plane grid have multiple images.  This boundary will surround all caustics unlike for GridMap::find_crit.

   This should not be as susceptible to missing the radial caustic because of resolution in the image plane.
   The resolution of the curve might not be so good.
   */
  void find_boundaries_of_caustics(std::vector<std::vector<Point_2d> > &boundaries
                           ,std::vector<bool> &hits_edge
                           ){
    
    size_t N=Ngrid_init*Ngrid_init2;
    std::vector<Point_2d> source_points(N);
    for(size_t i=0 ; i<N ; ++i) source_points[i] = i_points[i];
    
    std::vector<int> multiplicities(N);
    find_images(source_points,multiplicities);
    
    std::vector<bool> bitmap(N,false);
    for(size_t i=0 ; i<N ; ++i){
      if(multiplicities[i] > 1) bitmap[i] = true;
    }
    
    find_boundaries(bitmap,boundaries,hits_edge);
  }
  
  /**
   Calculate the magnification of one source by adding up its flux for the lensed image and an image made on an unlensed regulare grid
   */
  PosType magnificationFlux(Source &source) const ;
  
  /**\brief calculate the LOCAL magnification by triangel method weighted by interpolated surface brightness
   
   This is done by finding the area of every half cell triangle on the source plane and multiplying by the surface bightness interpolated to the
   center of the triangle on the image plane.  This does not use the point-wise magnification calculated by the rayshooter beacuse this can be highly unstable.
   
   NOTE: This will not equal the ratio of the lensed flux to the unlensed flux except in the case of one image (assuming the source is well resolved).
   */
  double magnificationTr() const ;

  /** \brief Same as `magnificationTr()`  but for a limited number of cells.  Problematic when cell is intersected by critical curve.
   */
  double magnificationTr(std::vector<size_t> &pixels) const ;

  /// area of a cell  (pixel size region with its lower left at point k) on source plane - calculated by triangal method
  double AreaCellOnSourcePlane(size_t k) const;

  
// depricated version of PixalMap::magnificationTr()
//  double magnificationTr2() const {
//
//    size_t k;
//    size_t k1;
//    size_t k2;
//
//    double flux_source=0.0,flux_image=0.0;
//    for(size_t i=1 ; i< Ngrid_init-1 ; i++){
//      for(size_t j=1 ; j< Ngrid_init2-1 ; j++){
//
//        k = i + j*Ngrid_init;
//
//        assert( s_points[k].surface_brightness == i_points[k].surface_brightness);
//        double sb = s_points[k].surface_brightness;
//
//        if(sb !=0.0){
//          flux_source += AreaCellOnSourcePlane(k) *sb;
//          flux_source += AreaCellOnSourcePlane(k-1) *sb;
//          flux_source += AreaCellOnSourcePlane(k-1-Ngrid_init) *sb;
//          flux_source += AreaCellOnSourcePlane(k-Ngrid_init) *sb;
//        }
//
//        flux_image += sb;
//      }
//    }
//
//    return flux_image * getResolution() * getResolution() * 4 / flux_source;
//  }

  /** \brief add flux to the rays that are nearest to the source on the source plane for each image
  *
   * This uses GridMap::find_images to find the images.  It then finds the point that is closest to the source position.
   *  The flux is added to one point per image.  The total flux added is returned.  No further refinement of the grid is done
   *  so it is limited by the resolution of the GridMap.  Some spurious low magnification images can be found.
   */
  double AddPointSource(const Point_2d &y,double flux);
  
  void find_crit(std::vector<std::vector<Point_2d> > &points
                 ,std::vector<bool> &hits_boundary
                 ,std::vector<CritType> &crit_type
                 );
  
  /** finds ordered boundaries to regions where bitmap == true

   This can be used to find critical curves or contours.
   `bitmap` should be the same size as the `Gridmap`
   If the boundary curve  touches the edge of the `GridMap` it will be indicated in `hits_boundary` as
   `true`.
   
   Boundaries will never cross or lead off the grid.  On the edges they will leave the edge pixels out even if they should be in.  This is a technical comprimise.
  */
  void find_boundaries(std::vector<bool> &bitmap  // = true inside
                       ,std::vector<std::vector<Point_2d> > &points
                       ,std::vector<bool> &hits_edge
                       ,bool add_to_vector=false
                       );
  
private:
  
  // curve must be in pixel units
  bool  incurve(long k,std::vector<Point_2d> &curve) const;
  
  // cluge to make compatible with old method of producing points
  Point * NewPointArray(size_t N){
    Point * p = point_factory(N);
    p[0].head = N;
    for(size_t i=1; i < N ; ++i) p[i].head = 0;
    return p;
  }
  
  void xygridpoints(Point *points,double range,const double *center,long Ngrid
                    ,short remove_center);
  
	/// one dimensional size of initial grid
	int Ngrid_init;
  int Ngrid_init2;
  
  unsigned long pointID;
  PosType axisratio;
  PosType x_range;
  void writePixelMapUniform_(Point* points,size_t size,PixelMap *map,LensingVariable val);
  
  Point *i_points;
  Point *s_points;
  Point_2d center;
  
  bool to_refine(long i,long j,double total,double f) const ;
  static std::mutex grid_mutex;
  
  MemmoryBank<Point> point_factory;
  
  void _find_images_(Point_2d *ys,int *multiplicity,long Nys,std::list<RAY> &rays) const;

  // find if there are images of y in specific cells
  void limited_image_search(Point_2d &y
                   ,std::vector<size_t> &cell_numbers  /// positions of the images limited by resolution of the gridmap
                   ,std::vector<Triangle> &triangles     /// index's of the points that form the triangles that the images are in
  ) const;
};

#endif // defined(__GLAMER__gridmap__)
