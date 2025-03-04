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
#include <memory>
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
  double RefreshSurfaceBrightnesses(Source* source);
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
  double AddSurfaceBrightnesses(Source* source);

  /// get the image point for a index number
  Point_2d image_point(size_t index){return i_points[index];}
  /// get the image point for a index number
  Point_2d source_point(size_t index){return s_points[index];}
  /// get the image point for a index number
  RAY ray(size_t index){return i_points[index];}
 
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
  template<typename T>
  PixelMap<T> writePixelMap(LensingVariable lensvar);
   /// fits output of lensing quantities at the resolution of the GridMap
  template <typename T>
  void writeFits(LensingVariable lensvar,std::string filensame);

  template<typename T>
  void writePixelMapUniform(PixelMap<T> &map,LensingVariable lensvar);
  template <typename T>
  void writeFitsUniform(const PosType center[],size_t Nx,size_t Ny,LensingVariable lensvar,std::string filename);
  template<typename T>
  PixelMap<T> writePixelMapUniform(const PosType center[],size_t Nx,size_t Ny,LensingVariable lensvar);

  /// this will make a fits map of the grid as is.
  template<typename T>
  void writeFitsUniform(
                        LensingVariable lensvar    ///< quantity to be output
                        ,std::string filename     ///< name of output fits file
                        ){
    PixelMap<T> map = writePixelMap<T>(lensvar);
    map.printFITS(filename);
  }
  
  /// returns a PixelMap with the flux in pixels at a resolution of res times the original resolution
  template<typename T>
  PixelMap<T> getPixelMapFlux() const;
  
  /// update a PixelMap with the flux in pixels.
  /// The map must have precisely the right size and center to match or an exception will be thrown.
  /// Constructing the map with PixelMap getPixelMapFlux(i) will insure that it does.
  template<typename T>
  void getPixelMapFlux(PixelMap<T> &map) const;
  
  /// returns the area (radians^2) of the region with negative magnification at resolution of fixed grid
  PosType EinsteinArea() const;

  /** flux weighted local magnification with current surface brightness averaged on the image plane,
   */
  //PosType magnification_invmag() const;
  //PosType magnification2() const;
  
  /// returns centroid of flux on the grid
  Point_2d centroid() const;
  
  Point_2d getCenter(){return center;}
  
  Point * operator[](size_t i){return i_points.data() + i;};
  
  GridMap(GridMap &&grid){
    
    Ngrid_init = grid.Ngrid_init;
    Ngrid_init2 = grid.Ngrid_init2;
    pointID = grid.pointID;
    axisratio = grid.axisratio;
    x_range = grid.x_range;
    
    std::swap(i_points,grid.i_points);
    std::swap(s_points,grid.s_points);
    
    center = grid.center;
  }
  
  GridMap(GridMap &grid){
    
    Ngrid_init = grid.Ngrid_init;
    Ngrid_init2 = grid.Ngrid_init2;
    pointID = grid.pointID;
    axisratio = grid.axisratio;
    x_range = grid.x_range;
    
    i_points = grid.i_points;
    s_points = grid.s_points;
    
    center = grid.center;
  }
  
  GridMap & operator=(GridMap &&grid){
    Ngrid_init = grid.Ngrid_init;
    Ngrid_init2 = grid.Ngrid_init2;
    pointID = grid.pointID;
    axisratio = grid.axisratio;
    x_range = grid.x_range;
    
    std::swap(i_points,grid.i_points);
    std::swap(s_points,grid.s_points);
    
    center = grid.center;

    return *this;
  }

  struct Triangle{
    Triangle(size_t i,size_t j,size_t k){
      index[0] = i;
      index[1] = j;
      index[2] = k;
    }
    Triangle(){
      index[0] = 0;
      index[1] = 0;
      index[2] = 0;
    }
    Triangle(const Triangle &tri){
      index[0] = tri[0];
      index[1] = tri[1];
      index[2] = tri[2];
    }
    Triangle & operator=(const Triangle &tri){
      index[0] = tri[0];
      index[1] = tri[1];
      index[2] = tri[2];
      return *this;
    }
    ~Triangle(){};
    std::vector<long> index = {0,0,0};
    long & operator[](int i){return index[i];}
    long operator[](int i) const {return index[i];}
  };
  struct Rectangle{
    Rectangle(size_t i,size_t j){
      index[0] = i;
      index[1] = j;
    }
    Rectangle(){
      index[0] = 0;
      index[1] = 0;
    }
    Rectangle(const Rectangle &rec){
      index[0] = rec[0];
      index[1] = rec[1];
    }
    Rectangle & operator=(Rectangle &rec){
      index[0] = rec[0];
      index[1] = rec[1];
      return *this;
    }
    ~Rectangle(){};
    std::vector<long> index = {0,0};
    long & operator[](int i){return index[i];}
    long operator[](int i) const {return index[i];}
  };
  
  // determines if two rectangles touch
  inline bool touch(const Rectangle &tr1,const Rectangle &tr2) const;
  
  // returns a vector of rectangles that encompase all triangles with one rectangle for each touching group
  std::vector<Rectangle> merge_boxes(
                   std::list<Triangle> &triangles
                   ) const;
  std::vector<Rectangle> merge_boxes(
                   std::vector<Triangle> &triangles
                   ) const;
  /*** \brief Returns a list of  RAYs from a set of source positions.
   
   The image positions are found in parallel by the triangle method.  The order of the
   output rays will be the same as the sources with multiple images consecutive.  The number
      of images for each source position is given by the `multiplicity` array.
   
   No new rays are shot.  The image positions and magnification matrix are interpolated from the nearest
      image points already in the GridMap.
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
  
  /** parallel version of find_images()
   */
  void find_images2(Point_2d y
                   ,std::vector<Point_2d> &image_points  /// positions of the images limited by resolution of the gridmap
                   ,std::vector<Triangle> &triangles     /// index's of the points that form the triangles that the images are in
  ) const ;
  
  /** \brief finds the boundary of the region on the source plane where there are more than one image

   Warning : slow but perhaps more reliable than find_caustics() when no radial caustic is found.
   
   This uses the triangle method to determine which points in a source plane grid of the same size and resolution as the image plane grid have multiple images.  This boundary will surround all caustics unlike for GridMap::find_crit.

   This should not be as susceptible to missing the radial caustic because of resolution in the image plane.
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
    
    Utilities::find_boundaries<Point_2d>(bitmap,Ngrid_init,boundaries,hits_edge);
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
  
  /** \brief Find critical curves.  This is usually not used outside of ImageFinding::find_crit()
   
   This will find all the resolved tangential and radial critical curves.  If a radial critical curve is not found
   within a tangential one, curves around the maxima are used to estimate a radial or pseudo caustic.  These are labeled CritType::pseudo.  The out put is ordered so that the radia/pseudo curves within a tangent curve imediately follow it.
   */
  void find_crit(std::vector<std::vector<Point_2d> > &points
                 ,std::vector<bool> &hits_boundary
                 ,std::vector<CritType> &crit_type
                 );

  /** \brief Find image-plane contours of magnification.  
   * This is usually only used within ImageFinding:: functions where it will also find the contours on the source plane.
   */
  void find_magnification_contour(
      std::vector<std::vector<Point_2d> > &curves
      ,std::vector<bool> &hits_boundary
      ,double invmag
  );
  
//  void find_crit_boundary(std::vector<std::vector<Point_2d> > &points
//                          ,std::vector<bool> &hits_boundary
//                          ) const;
  int getNx(){return Ngrid_init;}
  int getNy(){return Ngrid_init2;}
private:
  GridMap & operator=(GridMap &grid);
  
  // curve must be in pixel units
  bool  incurve(long k,std::vector<Point_2d> &curve) const;
    
  // cluge to make compatible with old method of producing points
  std::vector<Point> NewPointArray(size_t N){
    std::vector<Point> p(N);
    p[0].head = N;
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
  template<typename T>
  void writePixelMapUniform_(Point* points,size_t size,PixelMap<T> *map,LensingVariable val);
  
  std::vector<Point> i_points;
  std::vector<Point> s_points;
  Point_2d center;
  
  bool to_refine(long i,long j,double total,double f) const ;
  static std::mutex grid_mutex;
  
  void _find_images_(Point_2d *ys,int *multiplicity,long Nys,std::list<RAY> &rays) const;
  void _find_images2_(size_t j1
                              ,size_t j2
                              ,std::list<Point_2d> &image_points
                              ,std::list<Triangle> &triangles
                              ,Point_2d y
                              ) const;

  // find if there are images of y in specific cells
  void limited_image_search(Point_2d &y
                   ,std::vector<size_t> &cell_numbers  /// positions of the images limited by resolution of the gridmap
                   ,std::vector<Triangle> &triangles     /// index's of the points that form the triangles that the images are in
  ) const;
};

/// Output a PixelMap of the surface brightness with same res as the GridMap
template<typename T>
PixelMap<T> GridMap::getPixelMapFlux() const{
  
  // The number of pixels on a side of the new map will be
  // N = (Ngrid_init-1)/resf + 1;
  // so that the resolution is resf x the GridMap resolution
  
  PixelMap<T> map(center.x
                  ,Ngrid_init
                  ,Ngrid_init2
                  ,x_range/(Ngrid_init-1));
  
  size_t index;
   size_t n = Ngrid_init*Ngrid_init2;
   for(size_t i=0 ; i<n ; ++i){
     index = map.find_index(i_points[i].x);
     map.data()[index] = i_points[i].surface_brightness;
   }
  
  //for(size_t i = 0 ; i < Ngrid_init ; ++i){
  //  for(size_t j = 0 ; j < Ngrid_init2 ; ++j){
  //    map.data()[i/resf + map.getNx() * (j / resf)] +=
  //    i_points[ i + Ngrid_init * j].surface_brightness/factor;
  //  }
  //}
  
  map.Renormalize(map.getResolution()*map.getResolution());
  
  return map;
}

/// Flux in pixels map
template<typename T>
void GridMap::getPixelMapFlux(PixelMap<T> &map) const{
  
  //int resf = (Ngrid_init-1)/(map.getNx()-1);
  
  if(map.getNx() != Ngrid_init) throw std::invalid_argument("PixelMap does not match GripMap! Use the other GridMap::getPixelMapFlux() to contruct a PixelMap.");
  if(map.getNy() != Ngrid_init2) throw std::invalid_argument("PixelMap does not match GripMap! Use the other GridMap::getPixelMapFlux() to contruct a PixelMap.");
  //if(map.getResolution() != x_range*resf/(Ngrid_init-1)) throw std::invalid_argument("PixelMap does not match GripMap resolution! Use the other GridMap::getPixelMapFlux() to contruct a PixelMap.");
  
  if(map.getCenter()[0] != center[0]) throw std::invalid_argument("PixelMap does not match GripMap!");
  if(map.getCenter()[1] != center[1]) throw std::invalid_argument("PixelMap does not match GripMap!");
  
  double res = getResolution();
  if((map.getResolution()-res) < 1.0e-6*res ) throw std::invalid_argument("PixelMap resolution does not match GripMap!");
  
  map.Clean();
  
  //long nx = map.getNx();
  //for(size_t i = 0 ; i < Ngrid_init ; ++i){
  //  for(size_t j = 0 ; j < Ngrid_init2 ; ++j){
  //    size_t k = i + nx * j;
  //    map.data()[k] += i_points[k].surface_brightness;
  //  }
  //}
  
   size_t index;
   size_t n = Ngrid_init*Ngrid_init2;
   for(size_t i=0 ; i<n ; ++i){
     index = map.find_index(i_points[i].x);
     map.data()[index] = i_points[i].surface_brightness;
   }
  
  
  map.Renormalize(map.getResolution()*map.getResolution());
}

/** \brief Make a Pixel map of the without distribution the pixels.
 *
 *  This will be faster than Grid::writePixelMap() and Grid::writeFits().
 *  But it puts each grid pixel in one pixelmap pixel and if there are two
 *  grid pixels in one pixelmap pixel it uses one at random.  This is meant
 *  for uniform maps to make equal sized PixelMaps.
 */
template<typename T>
PixelMap<T> GridMap::writePixelMapUniform(
                                       const PosType center[]  /// center of image
                                       ,size_t Nx       /// number of pixels in image in on dimension
                                       ,size_t Ny       /// number of pixels in image in on dimension
                                       ,LensingVariable lensvar  /// which quantity is to be displayed
){
  
  if(getNumberOfPoints() == 0 ) return PixelMap<T>();
  PixelMap<T> map(center, Nx, Ny,x_range/(Nx-1));
  
  map.Clean();
  
  writePixelMapUniform(map,lensvar);
  
  return map;
}

template<typename T>
void GridMap::writeFits(
                        LensingVariable lensvar /// which quantity is to be displayed
                        ,std::string filename  /// output files
                        ){
                          PixelMap<T> map = writePixelMap<T>(lensvar);
                          map.printFITS(filename);
}

template<typename T>
PixelMap<T> GridMap::writePixelMap(
              LensingVariable lensvar  /// which quantity is to be displayed
){
  size_t Nx =  Ngrid_init;
  size_t Ny = Ngrid_init2;
  
  PixelMap<T> map( center.x, Nx, Ny,x_range/(Nx-1) );
  
  size_t N = map.size();
  assert(N == Nx*Ny);
  
  double tmp2[2];
  switch (lensvar) {
    case LensingVariable::ALPHA:
      for(size_t i=0 ; i<N ; ++i){
        tmp2[0] = i_points[i].x[0] - i_points[i].image->x[0];
        tmp2[1] = i_points[i].x[1] - i_points[i].image->x[1];
        map[i] = sqrt(tmp2[0]*tmp2[0] + tmp2[1]*tmp2[1]);
      }
      break;
    case LensingVariable::ALPHA1:
      for(size_t i=0 ; i<N ; ++i)
        map[i] = (i_points[i].x[0] - i_points[i].image->x[0]);
      break;
    case LensingVariable::ALPHA2:
      for(size_t i=0 ; i<N ; ++i)
        map[i] = (i_points[i].x[1] - i_points[i].image->x[1]);
      break;
    case LensingVariable::KAPPA:
      for(size_t i=0 ; i<N ; ++i)
        map[i] = i_points[i].kappa();
      break;
    case LensingVariable::GAMMA:
       for(size_t i=0 ; i<N ; ++i){
         tmp2[0] = i_points[i].gamma1();
         tmp2[1] = i_points[i].gamma2();
         map[i] = sqrt(tmp2[0]*tmp2[0] + tmp2[1]*tmp2[1]);
      }
      break;
    case LensingVariable::GAMMA1:
      for(size_t i=0 ; i<N ; ++i)
        map[i] = i_points[i].gamma1();
      break;
    case LensingVariable::GAMMA2:
      for(size_t i=0 ; i<N ; ++i)
        map[i] = i_points[i].gamma2();
      break;
    case LensingVariable::GAMMA3:
      for(size_t i=0 ; i<N ; ++i)
        map[i] = i_points[i].gamma3();
      break;
    case LensingVariable::INVMAG:
      for(size_t i=0 ; i<N ; ++i)
        map[i] = i_points[i].invmag();
      break;
    case LensingVariable::DELAYT:
      for(size_t i=0 ; i<N ; ++i)
        map[i] = i_points[i].dt;
      break;
    case LensingVariable::SurfBrightness:
      for(size_t i=0 ; i<N ; ++i)
        map[i] = i_points[i].surface_brightness;
      break;
    default:
      std::cerr << "GridMap::writePixelMapUniform() does not work for the input LensingVariable" << std::endl;
      throw std::runtime_error("GridMap::writePixelMapUniform() does not work for the input LensingVariable");
      break;
      // If this list is to be expanded to include ALPHA or GAMMA take care to add them as vectors
  }
  return map;
}
template <typename T>
void GridMap::writePixelMapUniform(
                                   PixelMap<T> &map
                                   ,LensingVariable lensvar  /// which quantity is to be displayed
){
  
  if(getNumberOfPoints() ==0 ) return;
  
  map.Clean();
  
  //writePixelMapUniform_(i_points,getNumberOfPoints(),&map,lensvar);
  //return;
  
  std::vector<std::thread> thr;
  int nthreads = Utilities::GetNThreads();
  
  int chunk_size;
  do{
    chunk_size =  getNumberOfPoints()/nthreads;
    if(chunk_size == 0) nthreads /= 2;
  }while(chunk_size == 0);
  
  size_t size = chunk_size;
  for(int ii = 0; ii < nthreads ;++ii){
    if(ii == nthreads-1)
    size = getNumberOfPoints() - (nthreads-1)*chunk_size;
    thr.push_back(std::thread(&GridMap::writePixelMapUniform_<T>,this,&(i_points[ii*chunk_size]),size,&map,lensvar));
  }
  for(auto &t : thr) t.join();
}

template <typename T>
void GridMap::writePixelMapUniform_(Point* points,size_t size,PixelMap<T> *map,LensingVariable val){
  double tmp;
  PosType tmp2[2];
  long index;
  
  for(size_t i = 0; i< size; ++i){
    switch (val) {
      case LensingVariable::ALPHA:
        tmp2[0] = points[i].x[0] - points[i].image->x[0];
        tmp2[1] = points[i].x[1] - points[i].image->x[1];
        tmp = sqrt(tmp2[0]*tmp2[0] + tmp2[1]*tmp2[1]);
        break;
      case LensingVariable::ALPHA1:
        tmp = (points[i].x[0] - points[i].image->x[0]);
        break;
      case LensingVariable::ALPHA2:
        tmp = (points[i].x[1] - points[i].image->x[1]);
        break;
      case LensingVariable::KAPPA:
        tmp = points[i].kappa();
        break;
      case LensingVariable::GAMMA:
        tmp2[0] = points[i].gamma1();
        tmp2[1] = points[i].gamma2();
        tmp = sqrt(tmp2[0]*tmp2[0] + tmp2[1]*tmp2[1]);
        break;
      case LensingVariable::GAMMA1:
        tmp = points[i].gamma1();
        break;
      case LensingVariable::GAMMA2:
        tmp = points[i].gamma2();
        break;
      case LensingVariable::GAMMA3:
        tmp = points[i].gamma3();
        break;
      case LensingVariable::INVMAG:
        tmp = points[i].invmag();
        break;
      case LensingVariable::DELAYT:
        tmp = points[i].dt;
        break;
      case LensingVariable::SurfBrightness:
        tmp = points[i].surface_brightness;
        break;
      default:
        std::cerr << "PixelMap<T>::AddGrid() does not work for the input LensingVariable" << std::endl;
        throw std::runtime_error("PixelMap<T>::AddGrid() does not work for the input LensingVariable");
        break;
        // If this list is to be expanded to include ALPHA or GAMMA take care to add them as vectors
    }
    
    index = map->find_index(points[i].x);
    if(index != -1)(*map)[index] = tmp;
  }
}

template <typename T>
void GridMap::writeFitsUniform(
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

#endif // defined(__GLAMER__gridmap__)
