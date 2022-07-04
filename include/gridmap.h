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
#include "Tree.h"
#include "source.h"
#include <mutex>

/** 
 * \brief A simplified version of the Grid structure for making non-adaptive maps of the lensing quantities (kappa, gamma, etc...)
 *
 *  GripMap is faster and uses less memory than Grid.  It does not construct the tree structures for the points 
 *  and thus cannot be used for adaptive mapping or image finding.
 */

struct GridMap{
  
	GridMap(LensHndl lens,unsigned long N1d,const double center[2],double range);
  GridMap(LensHndl lens ,unsigned long Nx ,const PosType center[2] ,PosType rangeX ,PosType rangeY);
	~GridMap();
  
    /// reshoot the rays for example when the source plane has been changed
  void ReInitializeGrid(LensHndl lens);
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

  void ClearSurfaceBrightnesses();
	size_t getNumberOfPoints() const {return Ngrid_init*Ngrid_init2;}
  
	/// return initial number of grid points in each direction
	int getInitNgrid(){return Ngrid_init;}
	/// return initial range of gridded region.  This is the distance from the first ray in a row to the last (unlike PixelMap)
	double getXRange(){return x_range;}
	double getYRange(){return x_range*axisratio;}
  /// resolution in radians, this is range / (N-1)
  double getResolution(){return x_range/(Ngrid_init-1);}
  
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

  /// flux weighted magnification with current surface brightness averaged on the image plane
  PosType magnification() const;
  PosType magnification2() const;
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

    return *this;
  }

private:
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
};

#endif // defined(__GLAMER__gridmap__)
