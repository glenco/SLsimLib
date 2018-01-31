//
//  gridmap.cpp
//  GLAMER
//
//  Created by bmetcalf on 8/19/14.
//
//

#include "gridmap.h"
#include <mutex>
#include <thread>

std::mutex GridMap::grid_mutex;

/** \ingroup Constructor
 * \brief Constructor for initializing rectangular grid.
 *
 * Cells of grid will always be square with initial resolution rangeX/(Nx-1).
 * The Y range may not be exactly rangeY, but will be the nearest value that
 * is a whole number of cells.
 *
 * Note: Deflection solver must be specified before creating a GridMap.
 */
GridMap::GridMap(
                 LensHndl lens       /// lens model for initializing grid
                 ,unsigned long Nx   /// Initial number of grid points in X dimension.
                 ,const PosType my_center[2]  /// Center of grid.
                 ,PosType rangeX   /// Full width of grid in x direction in whatever units will be used.
                 ,PosType rangeY  /// Full width of grid in y direction in whatever units will be used.
): Ngrid_init(Nx),axisratio(rangeY/rangeX),x_range(rangeX){
  
  pointID = 0;
  
  center = my_center;

  assert(Nx > 0);
  assert(rangeX > 0 && rangeY >0);
  
  if(Nx <= 0){ERROR_MESSAGE();
    std::cout << "cannot make GridMap with no points" << std::endl;
    throw std::runtime_error("");
  }
  if(rangeX <= 0 || rangeY <= 0 ){ERROR_MESSAGE();
    throw std::runtime_error("");
    std::cout << "cannot make GridMap with no range" << std::endl;
  }
  
  //if(Ngrid_init % 2 == 1 ) ++Ngrid_init;
  
  Ngrid_init2 = (int)(Ngrid_init*axisratio);
  //if(Ngrid_init2 % 2 == 1) ++Ngrid_init2;
  
  i_points = NewPointArray(Ngrid_init*Ngrid_init2);
  
  int i;
  // set grid of positions
  for(int ii=0;ii<Ngrid_init;++ii){
    for(int jj=0;jj<Ngrid_init2;++jj){
      
      i = ii + jj*Ngrid_init;
      
      i_points[i].id=pointID;
      ++pointID;
      
      i_points[i].x[0] = center[0] + rangeX*ii/(Ngrid_init-1) - 0.5*rangeX;
      i_points[i].x[1] = center[1] + rangeX*jj/(Ngrid_init-1) - 0.5*rangeY;
      i_points[i].gridsize=rangeX/(Ngrid_init-1);
    }
  }
  
  assert(i == Ngrid_init*Ngrid_init2-1);
  
  s_points=LinkToSourcePoints(i_points,Ngrid_init*Ngrid_init2);
  {
    std::lock_guard<std::mutex> hold(grid_mutex);
    lens->rayshooterInternal(Ngrid_init*Ngrid_init2,i_points);
  }
}

/** \ingroup Constructor
 * \brief Constructor for initializing square grid.
 *
 * Note: Deflection solver must be specified before creating a GridMap.
 */
GridMap::GridMap(
                 LensHndl lens               /// lens model for initializing grid
                 ,unsigned long N1d          /// Initial number of grid points in each dimension.
                 ,const PosType my_center[2]    /// Center of grid.
                 ,PosType range              /// Full width of grid in whatever units will be used.
): Ngrid_init(N1d),Ngrid_init2(N1d),axisratio(1.0),x_range(range){
  
  pointID = 0;
  
  assert(N1d > 0);
  assert(range > 0);
  
  center = my_center;
  
  if(N1d <= 0){ERROR_MESSAGE(); std::cout << "cannot make GridMap with no points" << std::endl; exit(1);}
  if(range <= 0){ERROR_MESSAGE(); std::cout << "cannot make GridMap with no range" << std::endl; exit(1);}
  
  i_points = NewPointArray(Ngrid_init*Ngrid_init);
  xygridpoints(i_points,range,center.x,Ngrid_init,0);
  s_points=LinkToSourcePoints(i_points,Ngrid_init*Ngrid_init);
    
  {
    std::lock_guard<std::mutex> hold(grid_mutex);
    lens->rayshooterInternal(Ngrid_init*Ngrid_init,i_points);
  }
  
}

GridMap::~GridMap(){
  FreePointArray(i_points);
  FreePointArray(s_points);
}

void GridMap::ReInitializeGrid(LensHndl lens){
  
  {
    std::lock_guard<std::mutex> hold(grid_mutex);
    lens->rayshooterInternal(Ngrid_init*Ngrid_init,i_points);
  }
  ClearSurfaceBrightnesses();
}

/// Output a PixelMap of the surface brightness with same res as the GridMap
PixelMap GridMap::getPixelMap(int resf) const{
  
  if(resf <=0){
    ERROR_MESSAGE();
    throw std::invalid_argument("resf must be > 0");
  }
  
  // The number of pixels on a side of the new map will be
  // N = (Ngrid_init-1)/resf + 1;
  // so that the resolution is resf x the GridMap resolution
  
  PixelMap map(center.x,(Ngrid_init-1)/resf + 1 ,(Ngrid_init2-1)/resf + 1,resf*x_range/(Ngrid_init-1));
  
  int factor = resf*resf;
  for(size_t i = 0 ; i < Ngrid_init ; ++i){
    for(size_t j = 0 ; j < Ngrid_init2 ; ++j){
      map.data()[i/resf + map.getNx() * (j / resf)] +=
      i_points[ i + Ngrid_init * j].surface_brightness/factor;
    }
  }
  
  map.Renormalize(map.getResolution()*map.getResolution());
  
  return map;
}

/// surface brightness map
void GridMap::getPixelMap(PixelMap &map) const{
  
  int resf = (Ngrid_init-1)/(map.getNx()-1);
  
  if(resf*map.getNx() != Ngrid_init-1+resf) throw std::invalid_argument("PixelMap does not match GripMap! Use the other GridMap::getPixelMap() to contruct a PixelMap.");
  if(resf*map.getNy() != Ngrid_init2-1+resf) throw std::invalid_argument("PixelMap does not match GripMap! Use the other GridMap::getPixelMap() to contruct a PixelMap.");
  if(map.getResolution() != x_range*resf/(Ngrid_init-1)) throw std::invalid_argument("PixelMap does not match GripMap resolution! Use the other GridMap::getPixelMap() to contruct a PixelMap.");
  
  if(map.getCenter()[0] != center[0]) throw std::invalid_argument("PixelMap does not match GripMap!");
  if(map.getCenter()[1] != center[1]) throw std::invalid_argument("PixelMap does not match GripMap!");
  
  if(resf <=0){
    ERROR_MESSAGE();
    throw std::invalid_argument("resf must be > 0");
  }
  
  map.Clean();
  
  int factor = resf*resf;
  for(size_t i = 0 ; i < Ngrid_init ; ++i){
    for(size_t j = 0 ; j < Ngrid_init2 ; ++j){
      map.data()[i/resf + map.getNx() * (j / resf)] +=
      i_points[ i + Ngrid_init * j].surface_brightness/factor;
    }
  }
  
  map.Renormalize(map.getResolution()*map.getResolution());
}

double GridMap::RefreshSurfaceBrightnesses(SourceHndl source){
  PosType total=0,tmp;
  
  for(size_t i=0;i <s_points[0].head;++i){
    tmp = source->SurfaceBrightness(s_points[i].x);
    s_points[i].surface_brightness = s_points[i].image->surface_brightness
    = tmp;
    total += tmp;
    s_points[i].in_image = s_points[i].image->in_image = NO;
  }
  
  return total;
}

void GridMap::ClearSurfaceBrightnesses(){
  
  for(size_t i=0;i <s_points[0].head;++i){
    s_points[i].surface_brightness = s_points[i].image->surface_brightness
    = 0.0;
    s_points[i].in_image = s_points[i].image->in_image = NO;
  }
  
}

/** \brief Make a Pixel map of the without distribution the pixels.
 *
 *  This will be faster than Grid::writePixelMap() and Grid::writeFits().
 *  But it puts each grid pixel in one pixelmap pixel and if there are two
 *  grid pixels in one pixelmap pixel it uses one at random.  This is meant
 *  for uniform maps to make equal sized PixelMaps.
 */
PixelMap GridMap::writePixelMapUniform(
                                       const PosType center[]  /// center of image
                                       ,size_t Nx       /// number of pixels in image in on dimension
                                       ,size_t Ny       /// number of pixels in image in on dimension
                                       ,LensingVariable lensvar  /// which quantity is to be displayed
){
  
  if(getNumberOfPoints() == 0 ) return PixelMap();
  PixelMap map(center, Nx, Ny,x_range/(Nx-1));
  
  map.Clean();
  
  writePixelMapUniform(map,lensvar);
  
  return map;
}

PixelMap GridMap::writePixelMapUniform(
              LensingVariable lensvar  /// which quantity is to be displayed
){
  size_t Nx =  Ngrid_init;
  size_t Ny = Ngrid_init2;
  
  PixelMap map( center.x, Nx, Ny,x_range/(Nx-1) );
  map.Clean();
  
  writePixelMapUniform(map,lensvar);
  
  return map;
}
void GridMap::writePixelMapUniform(
                                   PixelMap &map
                                   ,LensingVariable lensvar  /// which quantity is to be displayed
){
  
  if(getNumberOfPoints() ==0 ) return;
  
  map.Clean();
  
  std::thread thr[20];
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
    thr[ii] = std::thread(&GridMap::writePixelMapUniform_,this,&(i_points[ii*chunk_size]),size,&map,lensvar);
  }
  for(int ii = 0; ii < nthreads ;++ii) thr[ii].join();
  
}

void GridMap::writePixelMapUniform_(Point* points,size_t size,PixelMap *map,LensingVariable val){
  double tmp;
  PosType tmp2[2];
  long index;
  
  for(size_t i = 0; i< size; ++i){
    switch (val) {
      case ALPHA:
        tmp2[0] = points[i].x[0] - points[i].image->x[0];
        tmp2[1] = points[i].x[1] - points[i].image->x[1];
        tmp = sqrt(tmp2[0]*tmp2[0] + tmp2[1]*tmp2[1]);
        break;
      case ALPHA1:
        tmp = (points[i].x[0] - points[i].image->x[0]);
        break;
      case ALPHA2:
        tmp = (points[i].x[1] - points[i].image->x[1]);
        break;
      case KAPPA:
        tmp = points[i].kappa;
        break;
      case GAMMA:
        tmp2[0] = points[i].gamma[0];
        tmp2[1] = points[i].gamma[1];
        tmp = sqrt(tmp2[0]*tmp2[0] + tmp2[1]*tmp2[1]);
        break;
      case GAMMA1:
        tmp = points[i].gamma[0];
        break;
      case GAMMA2:
        tmp = points[i].gamma[1];
        break;
      case GAMMA3:
        tmp = points[i].gamma[2];
        break;
      case INVMAG:
        tmp = points[i].invmag;
        break;
      case DT:
        tmp = points[i].dt;
        break;
      default:
        std::cerr << "PixelMap::AddGrid() does not work for the input LensingVariable" << std::endl;
        throw std::runtime_error("PixelMap::AddGrid() does not work for the input LensingVariable");
        break;
        // If this list is to be expanded to include ALPHA or GAMMA take care to add them as vectors
    }
    
    index = map->find_index(points[i].x);
    if(index != -1)(*map)[index] = tmp;
  }
}

void GridMap::writeFitsUniform(
                               const PosType center[]  /// center of image
                               ,size_t Nx       /// number of pixels in image in on dimension
                               ,size_t Ny       /// number of pixels in image in on dimension
                               ,LensingVariable lensvar  /// which quantity is to be displayed
                               ,std::string filename     /// file name for image -- .kappa.fits, .gamma1.fits, etc will be appended
){
  std::string tag;
  
  switch (lensvar) {
    case DT:
      tag = ".dt.fits";
      break;
    case ALPHA1:
      tag = ".alpha1.fits";
      break;
    case ALPHA2:
      tag = ".alpha2.fits";
      break;
    case ALPHA:
      tag = ".alpha.fits";
      break;
    case KAPPA:
      tag = ".kappa.fits";
      break;
    case GAMMA1:
      tag = ".gamma1.fits";
      break;
    case GAMMA2:
      tag = ".gamma2.fits";
      break;
    case GAMMA3:
      tag = ".gamma3.fits";
      break;
    case GAMMA:
      tag = ".gamma.fits";
      break;
    case INVMAG:
      tag = ".invmag.fits";
      break;
    default:
      break;
  }
  
  PixelMap map = this->writePixelMapUniform(center,Nx,Ny,lensvar);
  map.printFITS(filename + tag);
}

void GridMap::xygridpoints(Point *i_points,PosType range,const PosType *center,long Ngrid_1d,short remove_center){
  /* make a new rectolinear grid of points on the image plane **/
  /* and link them to points on the source plane **/
  /* remove_center = 0 include center point of grid */
  /*              != 1 leave out center point of grid if Ngrid_1d is odd, Ngrid_1d*Ngrid_1d-1 points outputted */
  /* warning: memory for i_points must be allocated before entering */
  long i,j;
  
  if(remove_center && (Ngrid_1d%2 == 1)){
    /*i_points=NewPointArray(Ngrid_1d*Ngrid_1d-1);*/
    for(i=0,j=0;i<Ngrid_1d*Ngrid_1d;++i){
      
      if( (2*(i/Ngrid_1d)/(Ngrid_1d-1) == 1) && (i%Ngrid_1d == Ngrid_1d/2+1) ) j=1;
      i_points[i-j].id=pointID;
      ++pointID;
      Utilities::PositionFromIndex(i,i_points[i-j].x,Ngrid_1d,range,center);
      i_points[i-j].gridsize=range/(Ngrid_1d-1);
    }
    
  }else{
    /*i_points=NewPointArray(Ngrid_1d*Ngrid_1d);*/
    for(i=0;i<Ngrid_1d*Ngrid_1d;++i){
      i_points[i].id=pointID;
      ++pointID;
      Utilities::PositionFromIndex(i,i_points[i].x,Ngrid_1d,range,center);
      i_points[i].gridsize=range/(Ngrid_1d-1);
    }
  }
  
  return;
}

PosType GridMap::EisnsteinArea() const{
  size_t count = 0;
  size_t N = Ngrid_init*Ngrid_init2;
  for(size_t i=0;i<N;++i){
    if(i_points[i].invmag < 0) ++count;
  }
  
  return count*x_range*x_range/Ngrid_init/Ngrid_init;
}

PosType GridMap::magnification() const{
  double mag = 0,flux = 0;
  size_t N = Ngrid_init*Ngrid_init2;
  for(size_t i=0;i<N;++i){
    mag += i_points[i].surface_brightness*fabs(i_points[i].invmag);
    flux += i_points[i].surface_brightness;
  }
  return flux/mag;
}

