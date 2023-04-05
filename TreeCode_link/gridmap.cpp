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

/** 
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
  //i_points = point_factory(Ngrid_init*Ngrid_init2);
 
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
  
  s_points = NewPointArray(Ngrid_init*Ngrid_init2);
  //s_points = point_factory(Ngrid_init*Ngrid_init2);
  LinkToSourcePoints(i_points,s_points,Ngrid_init*Ngrid_init2);
  {
    std::lock_guard<std::mutex> hold(grid_mutex);
    lens->rayshooterInternal(Ngrid_init*Ngrid_init2,i_points);
  }
}

/** 
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
  
  if(N1d < 1) throw std::runtime_error("GridMap size is < 1!");
  if(range <= 0) throw std::runtime_error("GridMap range is <= 0!");
  
  center = my_center;
  
  if(N1d <= 0){ERROR_MESSAGE(); std::cout << "cannot make GridMap with no points" << std::endl; exit(1);}
  if(range <= 0){ERROR_MESSAGE(); std::cout << "cannot make GridMap with no range" << std::endl; exit(1);}
  
  i_points = NewPointArray(Ngrid_init*Ngrid_init);
  xygridpoints(i_points,range,center.x,Ngrid_init,0);
  s_points = NewPointArray(Ngrid_init*Ngrid_init);
  LinkToSourcePoints(i_points,s_points,Ngrid_init*Ngrid_init);
    
  {
    std::lock_guard<std::mutex> hold(grid_mutex);
    lens->rayshooterInternal(Ngrid_init*Ngrid_init,i_points);
  }
}

GridMap::GridMap(
                 unsigned long N1d          /// Initial number of grid points in each dimension.
                 ,const PosType my_center[2]    /// Center of grid.
                 ,PosType range              /// Full width of grid in whatever units will be used.
): Ngrid_init(N1d),Ngrid_init2(N1d),axisratio(1.0),x_range(range){
  
  pointID = 0;
  
  if(N1d < 1) throw std::runtime_error("GridMap size is < 1!");
  if(range <= 0) throw std::runtime_error("GridMap range is <= 0!");
  
  center = my_center;
  
  if(N1d <= 0){ERROR_MESSAGE(); std::cout << "cannot make GridMap with no points" << std::endl; exit(1);}
  if(range <= 0){ERROR_MESSAGE(); std::cout << "cannot make GridMap with no range" << std::endl; exit(1);}
  
  i_points = NewPointArray(Ngrid_init*Ngrid_init);
  xygridpoints(i_points,range,center.x,Ngrid_init,0);
  s_points = NewPointArray(Ngrid_init*Ngrid_init);
  xygridpoints(s_points,range,center.x,Ngrid_init,0);
  LinkToSourcePoints(i_points,s_points,Ngrid_init*Ngrid_init);
}

GridMap::~GridMap(){}

GridMap GridMap::ReInitialize(LensHndl lens){
  
  GridMap newgrid(lens,Ngrid_init,center.x,x_range,getYRange());
  
  return newgrid;
  /*
  {
    std::lock_guard<std::mutex> hold(grid_mutex);
    lens->rayshooterInternal(Ngrid_init*Ngrid_init,i_points);
  }
  ClearSurfaceBrightnesses();
   */
}

/// Output a PixelMap of the surface brightness with same res as the GridMap
PixelMap GridMap::getPixelMapFlux(int resf) const{
  
  if(resf <=0){
    ERROR_MESSAGE();
    throw std::invalid_argument("resf must be > 0");
  }
  
  // The number of pixels on a side of the new map will be
  // N = (Ngrid_init-1)/resf + 1;
  // so that the resolution is resf x the GridMap resolution
  
  PixelMap map(center.x,(Ngrid_init-1)/resf + 1 ,(Ngrid_init2-1)/resf + 1,resf*x_range/(Ngrid_init-1));
  
  int factor = resf*resf;
  size_t index;
   size_t n = Ngrid_init*Ngrid_init2;
   for(size_t i=0 ; i<n ; ++i){
     index = map.find_index(i_points[i].x);
     map.data()[index] = i_points[i].surface_brightness/factor;
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

/// surface brightness map
void GridMap::getPixelMapFlux(PixelMap &map) const{
  
  int resf = (Ngrid_init-1)/(map.getNx()-1);
  
  if(resf*map.getNx() != Ngrid_init-1+resf) throw std::invalid_argument("PixelMap does not match GripMap! Use the other GridMap::getPixelMapFlux() to contruct a PixelMap.");
  if(resf*map.getNy() != Ngrid_init2-1+resf) throw std::invalid_argument("PixelMap does not match GripMap! Use the other GridMap::getPixelMapFlux() to contruct a PixelMap.");
  if(map.getResolution() != x_range*resf/(Ngrid_init-1)) throw std::invalid_argument("PixelMap does not match GripMap resolution! Use the other GridMap::getPixelMapFlux() to contruct a PixelMap.");
  
  if(map.getCenter()[0] != center[0]) throw std::invalid_argument("PixelMap does not match GripMap!");
  if(map.getCenter()[1] != center[1]) throw std::invalid_argument("PixelMap does not match GripMap!");
  
  if(resf <=0){
    ERROR_MESSAGE();
    throw std::invalid_argument("resf must be > 0");
  }
  
  map.Clean();
  
  int factor = resf*resf;
  size_t index;
  size_t n = Ngrid_init*Ngrid_init2;
  for(size_t i=0 ; i<n ; ++i){
    index = map.find_index(i_points[i].x);
    map.data()[index] = i_points[i].surface_brightness/factor;
  }
  
  //for(size_t i = 0 ; i < Ngrid_init ; ++i){
  //  for(size_t j = 0 ; j < Ngrid_init2 ; ++j){
  //    map.data()[i/resf + map.getNx() * (j / resf)] +=
  //    i_points[ i + Ngrid_init * j].surface_brightness/factor;
  //  }
  //}
  
  map.Renormalize(map.getResolution()*map.getResolution());
}

double GridMap::RefreshSurfaceBrightnesses(SourceHndl source){
  PosType total=0,tmp;
  
  double res2 = pow(getResolution(),2);
  for(size_t i=0;i <s_points[0].head;++i){
    tmp = source->SurfaceBrightness(s_points[i].x);
    s_points[i].surface_brightness = s_points[i].image->surface_brightness
    = tmp;
    total += tmp * res2;
    s_points[i].in_image = s_points[i].image->in_image = NO;
  }
  
  return total;
}

double GridMap::AdaptiveRefreshSurfaceBrightnesses(Lens &lens,Source &source){
  PosType f=1.0e-4;
  
  double resolution = getResolution();
  LinkedPoint point;

  double total_flux = RefreshSurfaceBrightnesses(&source)/resolution/resolution;
  double original_total_flux = total_flux; // ????
  if(total_flux == 0) return 0;
  
  for(long j=0 ; j < Ngrid_init2 ; ++j){
    size_t k = j*Ngrid_init;
    for(long i=0 ; i < Ngrid_init ; ++i,++k){

      if(to_refine(i,j,total_flux,f)){
 
        Point_2d ll = *(i_points[k].image);
        ll[0] -= resolution/2;
        ll[1] -= resolution/2;

        double new_flux = i_points[k].surface_brightness;
        double old_flux=0;
        int n = 1;

        total_flux -= new_flux;
        while ( fabs(old_flux-new_flux) > 1.0e-3*new_flux ){
          old_flux = new_flux;
          n *= 3;
          
          double res = resolution/n;
          new_flux = 0;

          for(int u=0 ; u<n ; ++u){
            point.x[0] = ll[0] + (0.5 + u)*res;
            for(int v=0 ; v<n ; ++v){
               point.x[1] = ll[1] + (0.5 + v)*res;

               lens.rayshooterInternal(1,&point);
               new_flux += source.SurfaceBrightness(point.image->x);
            }
          }
          
          new_flux /= n*n;
        }
        
        //std::cout << " change in pixel flux = " << (new_flux-original) / original << std::endl;
        i_points[k].surface_brightness = new_flux;
        total_flux += new_flux;
      }
      
    }
  }
  
  //std::cout << " change in total flux = " << (total_flux-original_total_flux) / original_total_flux << std::endl;

  assert( (total_flux-original_total_flux) <  0.1*original_total_flux);
  return total_flux * resolution * resolution;
}

bool GridMap::to_refine(long i,long j,double total,double f) const {
  double flux = i_points[i + j*Ngrid_init].surface_brightness;
  
  for(size_t ii=MAX(i-1,0) ; ii < MIN(i+2,Ngrid_init) ;++ii){
    for(size_t jj=MAX(j-1,0) ; jj < MIN(j+2,Ngrid_init2) ;++jj){
      size_t kk = ii + jj*Ngrid_init;
      if( fabs(flux - i_points[kk].surface_brightness ) / total > f ) return true;
    }
  }
  return false;
}


double GridMap::AddSurfaceBrightnesses(SourceHndl source){
  PosType total=0,tmp;
  
  double res2 = pow(getResolution(),2);
  for(size_t i=0;i <s_points[0].head;++i){
    tmp = source->SurfaceBrightness(s_points[i].x);
    s_points[i].surface_brightness += tmp;
    s_points[i].image->surface_brightness += tmp;
    total += tmp * res2;
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

void GridMap::assertNAN(){
  for(size_t i=0;i <s_points[0].head;++i){
    assert(!isnan(s_points[i].surface_brightness));
  }
  for(size_t i=0;i <i_points[0].head;++i){
     assert(!isnan(i_points[i].surface_brightness));
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

void GridMap::writeFits(
                        LensingVariable lensvar /// which quantity is to be displayed
                        ,std::string filename  /// output files
                        ){
                          PixelMap map = writePixelMap(lensvar);
                          map.printFITS(filename);
}

PixelMap GridMap::writePixelMap(
              LensingVariable lensvar  /// which quantity is to be displayed
){
  size_t Nx =  Ngrid_init;
  size_t Ny = Ngrid_init2;
  
  PixelMap map( center.x, Nx, Ny,x_range/(Nx-1) );
  
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
void GridMap::writePixelMapUniform(
                                   PixelMap &map
                                   ,LensingVariable lensvar  /// which quantity is to be displayed
){
  
  if(getNumberOfPoints() ==0 ) return;
  
  map.Clean();
  
  //writePixelMapUniform_(i_points,getNumberOfPoints(),&map,lensvar);
  //return;
  
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
  
  PixelMap map = writePixelMapUniform(center,Nx,Ny,lensvar);
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
    for(i=0,j=0;i<Ngrid_1d*Ngrid_1d;++i){
      
      if( (2*(i/Ngrid_1d)/(Ngrid_1d-1) == 1) && (i%Ngrid_1d == Ngrid_1d/2+1) ) j=1;
      i_points[i-j].id=pointID;
      ++pointID;
      Utilities::PositionFromIndex(i,i_points[i-j].x,Ngrid_1d,range,center);
      i_points[i-j].gridsize=range/(Ngrid_1d-1);
    }
    
  }else{
    for(i=0;i<Ngrid_1d*Ngrid_1d;++i){
      i_points[i].id=pointID;
      ++pointID;
      Utilities::PositionFromIndex(i,i_points[i].x,Ngrid_1d,range,center);
      i_points[i].gridsize=range/(Ngrid_1d-1);
    }
  }
  
  return;
}

PosType GridMap::EinsteinArea() const{
  size_t count = 0;
  size_t N = Ngrid_init*Ngrid_init2;
  for(size_t i=0;i<N;++i){
    if(i_points[i].invmag() < 0) ++count;
  }
  
  return count*x_range*x_range/Ngrid_init/Ngrid_init;
}

PosType GridMap::magnification() const{
  
  double mag = 0,flux = 0;
  size_t N = Ngrid_init*Ngrid_init2;
  for(size_t i=0;i<N;++i){
    mag += i_points[i].surface_brightness*fabs(i_points[i].invmag());
    flux += i_points[i].surface_brightness;
  }
  return flux/mag;
}


PosType GridMap::magnification2() const{
  double mag = 0,flux = 0;
  size_t N = Ngrid_init*Ngrid_init2;
  for(size_t i=0;i<N;++i){
    mag += i_points[i].surface_brightness;
    flux += i_points[i].surface_brightness/fabs(i_points[i].invmag());
  }
  return flux/mag;
}

Point_2d GridMap::centroid() const{
  double flux = 0;
  Point_2d centroid(0,0);
  
  size_t N = Ngrid_init*Ngrid_init2;
  for(size_t i=0;i<N;++i){
    centroid += i_points[i]*i_points[i].surface_brightness;
    flux += i_points[i].surface_brightness;
  }
  return centroid/flux;
}


