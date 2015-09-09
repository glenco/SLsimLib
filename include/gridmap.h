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

/** \ingroup ImageFinding
 * \brief A simplified version of the Grid structure for making non-adaptive maps of the lensing quantities (kappa, gamma, etc...)
 *
 *  GripMap is faster and uses less memory than Grid.  It does not construct the tree structures for the points 
 *  and thus cannot be used for adaptive mapping or image finding.
 */

//class PixelMap;

struct GridMap{
  
	GridMap(LensHndl lens,unsigned long N1d,const double center[2],double range);
  GridMap(LensHndl lens ,unsigned long Nx ,const PosType center[2] ,PosType rangeX ,PosType rangeY);
	~GridMap();
  
  /// reshoot the rays for example when the source plane has been changed
  void ReInitializeGrid(LensHndl lens);
	double RefreshSurfaceBrightnesses(SourceHndl source);
  void ClearSurfaceBrightnesses();
	size_t getNumberOfPoints() const {return Ngrid_init*Ngrid_init2;}
  
	/// return initial number of grid points in each direction
	int getInitNgrid(){return Ngrid_init;}
	/// return initial range of gridded region
	double getXRange(){return x_range;}
	double getYRange(){return x_range*axisratio;}
  
  PixelMap writePixelMapUniform(const PosType center[],size_t Nx,size_t Ny,LensingVariable lensvar);
  void writePixelMapUniform(PixelMap &map,LensingVariable lensvar);
  void writeFitsUniform(const PosType center[],size_t Nx,size_t Ny,LensingVariable lensvar,std::string filename);
  
  /// returns a PixelMap with the flux in pixels at a resolution of res times the original resolution
  PixelMap getPixelMap(int res);
  /// update a PixelMap with the flux in pixels at a resolution of res times the original resolution.
  /// The map must have precisely the right size and center to match or an exception will be thrown.
  /// Constructing the map with PixelMap getPixelMap(int res) will insure that it does.
  void getPixelMap(PixelMap &map);
  
private:
  void xygridpoints(Point *points,double range,const double *center,long Ngrid
                    ,short remove_center);
  
	/// one dimensional size of initial grid
	const int Ngrid_init;
  int Ngrid_init2;
  
  unsigned long pointID;
  PosType axisratio;
  PosType x_range;
  void writePixelMapUniform_(Point* points,size_t size,PixelMap *map,LensingVariable val);
  
  Point *i_points;
  Point *s_points;
  
  static std::mutex grid_mutex;
};

#endif /* defined(__GLAMER__gridmap__) */
