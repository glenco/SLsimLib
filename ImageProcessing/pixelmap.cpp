//
//  pixelmap.cpp
//  GLAMER
//
//  Created by Robert Benton Metcalf on 26/04/24.
//

#include <stdio.h>
#include <utility>

#include "pixelmap.h"
#include "grid_maintenance.h"
#include "gridmap.h"
#include "utilities_slsim.h"


template <typename T>
PixelMap<T>::PixelMap()
: map(), Nx(0), Ny(0), resolution(0), rangeX(0), rangeY(0),units(PixelMapUnits::ndef)
{
  center[0] = 0;
  center[1] = 0;
  
  map_boundary_p1[0] = 0;
  map_boundary_p1[1] = 0;
  map_boundary_p2[0] = 0;
  map_boundary_p2[1] = 0;
}
template <typename T>
PixelMap<T>::PixelMap(const PixelMap& other)
: map(other.map),
Nx(other.Nx), Ny(other.Ny), resolution(other.resolution), rangeX(other.rangeX), rangeY(other.rangeY),units(other.units)
{
  center[0] = other.center[0];
  center[1] = other.center[1];
  map_boundary_p1[0] = other.map_boundary_p1[0];
  map_boundary_p1[1] = other.map_boundary_p1[1];
  map_boundary_p2[0] = other.map_boundary_p2[0];
  map_boundary_p2[1] = other.map_boundary_p2[1];
}
/*
template<>
template<>
PixelMap<float>::PixelMap(const PixelMap<double>& other):
Nx(other.getNx()), Ny(other.getNy()), resolution(other.getResolution())
, rangeX(other.getRangeX()), rangeY(other.getRangeY()),units(other.getUnits())
{
  size_t n = Nx*Ny;
  map.resize(n);
  for(int i=0 ; i<n ; ++i) map[i] = other[i];
  Point_2d c=other.getCenter();
  center[0] = c[0];
  center[1] = c[1];
  
  map_boundary_p1[0] = center[0]-(Nx*resolution)/2.;
  map_boundary_p1[1] = center[1]-(Ny*resolution)/2.;
  map_boundary_p2[0] = center[0]+(Nx*resolution)/2.;
  map_boundary_p2[1] = center[1]+(Ny*resolution)/2.;
}
*/
// move constructor
template <typename T>
PixelMap<T>::PixelMap(PixelMap&& other)
:map(std::move(other.map)),Nx(0),Ny(0),resolution(0), rangeX(0), rangeY(0){
 
  Nx = other.Nx;
  Ny = other.Ny;
  resolution = other.resolution;
  rangeX = other.rangeX;
  rangeY = other.rangeY;
  std::copy(other.center, other.center + 2, center);
  std::copy(other.map_boundary_p1, other.map_boundary_p1 + 2, map_boundary_p1);
  std::copy(other.map_boundary_p2, other.map_boundary_p2 + 2, map_boundary_p2);
  
  other.Nx = 0;
  other.Ny = 0;
  other.resolution = 0;
  other.center[0] = 0;
  other.center[1] = 0;
  
  other.map_boundary_p1[0] = 0;
  other.map_boundary_p1[1] = 0;
  other.map_boundary_p2[0] = 0;
  other.map_boundary_p2[1] = 0;
}

/// make square PixelMap
template <typename T>
PixelMap<T>::PixelMap(
                   const PosType* center,  /// The location of the center of the map
                   std::size_t Npixels,  /// Number of pixels in one dimension of map.
                   PosType resolution        /// One dimensional range of map in whatever units the point positions are in
                   ,PixelMapUnits u
)
: map(0.0, Npixels*Npixels),
Nx(Npixels), Ny(Npixels), resolution(resolution),units(u)
{
  if(Npixels == 0 || resolution <=0) throw std::invalid_argument("invalid arguments");
  std::copy(center, center + 2, this->center);
  rangeX = resolution*Nx;
  rangeY = resolution*Ny;
  
  map_boundary_p1[0] = center[0]-(Npixels*resolution)/2.;
  map_boundary_p1[1] = center[1]-(Npixels*resolution)/2.;
  map_boundary_p2[0] = center[0]+(Npixels*resolution)/2.;
  map_boundary_p2[1] = center[1]+(Npixels*resolution)/2.;
}

/// make rectangular PixelMap with square pixels
template <typename T>
PixelMap<T>::PixelMap(
                   const PosType* center,  /// The location of the center of the map
                   std::size_t myNx,  /// Number of pixels in x dimension of map.
                   std::size_t myNy,  /// Number of pixels in y dimension of map.
                   PosType resolution  /// One dimensional range of map in whatever units the point positions are in
                   ,PixelMapUnits u
)
: map(0.0, myNx*myNy),
Nx(myNx), Ny(myNy), resolution(resolution),units(u)
{
  
  std::copy(center, center + 2, this->center);
  rangeX = resolution*Nx;
  rangeY = resolution*Ny;
  
  map_boundary_p1[0] = center[0]-(Nx*resolution)/2.;
  map_boundary_p1[1] = center[1]-(Ny*resolution)/2.;
  map_boundary_p2[0] = center[0]+(Nx*resolution)/2.;
  map_boundary_p2[1] = center[1]+(Ny*resolution)/2.;
}

/** \brief Constructs a PixelMap reading in a fits file
 * Infos about resolution, Npixels and center are read from the header.
 */
template <typename T>
PixelMap<T>::PixelMap(
                   std::string fitsfilename   /// file name of fits file to be read
                    ,double my_res         /// resolution (rad) of fits image if not given in fits file, use default or -1 otherwise
                   ,PixelMapUnits u
                   ,std::string extension  /// fits extension 1,2,,...
):units(u)
{
  
  if(fitsfilename.empty())
    throw std::invalid_argument("Please enter a valid filename for the FITS file input");

        
    if(!Utilities::IO::file_exists(fitsfilename)){
        std::cerr << "Problem with inputfile " << fitsfilename << std::endl;
        throw std::invalid_argument("bad file");
    }

  std::vector<long> cpsize;
  
  CPFITS_READ cpfits(fitsfilename,extension);
  
  //int bitpix;
  cpfits.imageDimensions(cpsize);
  
  Nx = cpsize[0];
  Ny = cpsize[1];
  
  int err = 0;
  
  err += cpfits.readKey("RA", center[0]);
  err += cpfits.readKey("DEC", center[1]);

  if(err){
    err = 0;
    err += cpfits.readKey("CRVAL1", center[0]);
    err += cpfits.readKey("CRVAL2", center[1]);
  }

  if(err)
  {
    center[0] = 0.0;
    center[1] = 0.0;
  }
  
  if(my_res == -1){
    // read the resolution
    err = 0;
    {
      double cdelt2;
      err += cpfits.readKey("CDELT1", my_res);
      err += cpfits.readKey("CDELT2", cdelt2);
      if(err == 0 && std::abs(my_res) - std::abs(cdelt2) > 1e-6)
        throw std::runtime_error("non-square pixels in FITS file " + fitsfilename);
    }
    if(err != 0)
    {
      err = 0;
      {
        err += cpfits.readKey("CD1_1", my_res);
        //err += cpfits.readKey("CD1_2", cd12);
        //err += cpfits.readKey("CD2_1", cd21);
        //err += cpfits.readKey("CD2_2", cd22);
        //if(err==0 && std::abs(my_res) - std::abs(cd22) > 1e-6)
        //  throw std::runtime_error("non-square pixels in FITS file " + fitsfilename);
        //if(cd12 || cd21)
        //  throw std::runtime_error("pixels not aligned with coordinates in FITS file " + fitsfilename);
      }
      if(err != 0){
        double ps;
        
        err = cpfits.readKey("PHYSICALSIZE",ps);
        
        if(err != 0){
          std::cerr << "PixelMap input fits field must have header keywords:" << std::endl
          << " PHYSICALSIZE - size of map in degrees" <<std::endl
          << " or CDELT1 and CDELT2 or CD1_1, DC1_2, CD2_1 and CD2_2" << std::endl;
          
          throw std::invalid_argument("bad header");
        }
        my_res = ps/Nx;
      }
    }
    resolution = fabs(my_res)*PI/180.;
  }else{
    resolution = my_res;
  }
  
  rangeX = resolution*Nx;
  rangeY = resolution*Ny;
  map_boundary_p1[0] = center[0] - (Nx*resolution)/2.;
  map_boundary_p1[1] = center[1] - (Ny*resolution)/2.;
  map_boundary_p2[0] = center[0] + (Nx*resolution)/2.;
  map_boundary_p2[1] = center[1] + (Ny*resolution)/2.;
  
  cpfits.read(map,cpsize);
  //std::cout << "map size : " << map[0] << std::endl;
  //std::cout << "map size : " << map.size() << std::endl;
}

/** \brief Creates a new PixelMap from a square region of a PixelMap.
 * If the region exceeds the boundaries of the original map, the new map is completed with zeros.
 */
template <typename T>
PixelMap<T>::PixelMap(const PixelMap& pmap,  /// Input PixelMap (from which the stamp is taken)
                   const PosType* my_center, /// center of the region to be duplicated (in rads)
                   std::size_t my_Npixels /// size of the region to be duplicated (in pixels)
)
: map(0.0, my_Npixels*my_Npixels),
Nx(my_Npixels), Ny(my_Npixels), resolution(pmap.resolution)
,units(pmap.units)
{
    std::copy(my_center, my_center + 2, this->center);
    rangeX = resolution*Nx;
    rangeY = resolution*Ny;
  
    map_boundary_p1[0] = my_center[0]-(Nx*resolution)/2.;
    map_boundary_p1[1] = my_center[1]-(Ny*resolution)/2.;
    map_boundary_p2[0] = my_center[0]+(Nx*resolution)/2.;
    map_boundary_p2[1] = my_center[1]+(Ny*resolution)/2.;
    
    units = pmap.getUnits();
    
    int edge[2];
    edge[0] = (my_center[0]-pmap.map_boundary_p1[0])/resolution - Nx/2;
    edge[1] = (my_center[1]-pmap.map_boundary_p1[1])/resolution - Ny/2;
    if (edge[0] > int(pmap.Nx) || edge[1] > int(pmap.Ny) || edge[0]+int(Nx) < 0 || edge[1]+int(Ny) < 0)
    {
      std::cout << "The region you selected is completely outside PixelMap!" << std::endl;
      throw std::runtime_error("Attempting to make Sub-PixelMap outside of parent PixelMap!");
    }
    for (unsigned long i=0; i < map.size(); ++i)
    {
      int ix = i%Nx;
      int iy = i/Nx;
      map[i] = 0;
      if (ix+edge[0] > 0 && ix+edge[0] < pmap.Nx && iy+edge[1] > 0 && iy+edge[1] < pmap.Ny)
        map[i] = pmap.map[ix+edge[0]+(iy+edge[1])*pmap.Nx];
    }
}

/// Produces a square cut-out of the input PixelMap
template <typename T>
PixelMap<T>::PixelMap(const PixelMap& in_map,  /// Input PixelMap (from which the stamp is taken)
                   long nx, /// lower left  pixels of in pmap
                   long ny, /// lower left  pixels of in  pmap
                   std::size_t my_Npixels /// size of the region to be duplicated (in pixels)
)
: map(0.0, my_Npixels*my_Npixels),
Nx(my_Npixels), Ny(my_Npixels), resolution(in_map.resolution)
,units(in_map.units)
{

  long iimin =  MAX<long>(0,nx);
  long iimax =  MIN<long>(iimin+Nx,in_map.Nx);
  long jjmin =  MAX<long>(0,ny);
  long jjmax =  MIN(jjmin+Nx,in_map.Ny);

  for(long ii=iimin ; ii < iimax ; ++ii){
    long i = ii - nx;
    for(long jj=jjmin; jj < jjmax ; ++jj){
      long j = jj - ny;
        
      map[i + Nx*j] = in_map(ii,jj);
    }
  }
  
  rangeX = resolution*Nx;
  rangeY = resolution*Ny;
  
  map_boundary_p1[0] = in_map.map_boundary_p1[0] + resolution*nx;
  map_boundary_p1[1] = in_map.map_boundary_p1[1] + resolution*ny;
  
  map_boundary_p2[0] = map_boundary_p1[0] + resolution * Nx;
  map_boundary_p2[1] = map_boundary_p1[1] + resolution * Ny;
  
  center[0] = (map_boundary_p2[0] + map_boundary_p1[0])/2;
  center[1] = (map_boundary_p2[1] + map_boundary_p1[1])/2;
}

/** \brief Creates a PixelMap at a different resolution.
 * The new counts are calculated integrating over the input pixels.
 * No interpolation or smoothing is performed.
 */
template <typename T>
PixelMap<T>::PixelMap(
                   const PixelMap& pmap
                   ,PosType res_ratio     /// resolution of map is res_ratio times the resolution of the input map
)
{
  resolution = res_ratio*pmap.resolution;
  Nx = pmap.Nx/res_ratio + .5;
  Ny = pmap.Ny/res_ratio + .5;
  rangeX = resolution*Nx;
  rangeY = resolution*Ny;
  center[0] = pmap.center[0];
  center[1] = pmap.center[1];
  map_boundary_p1[0] = center[0] - (Nx*resolution)/2.;
  map_boundary_p1[1] = center[1] - (Ny*resolution)/2.;
  map_boundary_p2[0] = center[0] + (Nx*resolution)/2.;
  map_boundary_p2[1] = center[1] + (Ny*resolution)/2.;
  units = pmap.units;
  
  map.resize(Nx*Ny);
  
  long old_Nx = pmap.Nx;
  long old_Ny = pmap.Ny;
  long ix, iy;
  PosType area;
  PosType old_p1[2];
  PosType old_p2[2];
  
  for(unsigned long i=0;i < map.size(); ++i)
  {
    ix = i%Nx;
    iy = i/Nx;
    map[ix+Nx*iy] = 0.;
    old_p1[0] = std::max(0,int(ix*res_ratio));
    old_p1[1] = std::max(0,int(iy*res_ratio));
    
    old_p2[0] = MIN(old_Nx-1,long((ix+1.)*res_ratio));
    old_p2[1] = MIN(old_Ny-1,long((iy+1.)*res_ratio));
    
    for (int old_iy = old_p1[1]; old_iy <= old_p2[1]; ++old_iy)
    {
      for (int old_ix = old_p1[0]; old_ix<= old_p2[0]; ++old_ix)
      {
        area = MIN(old_ix+0.5,(ix+1.)*res_ratio-0.5) - MAX(old_ix-0.5,ix*res_ratio-0.5);
        area *= MIN(old_iy+0.5,(iy+1.)*res_ratio-0.5) - MAX(old_iy-0.5,iy*res_ratio-0.5);
        map[ix+Nx*iy] += area*pmap.map[old_ix+old_Nx*old_iy];
      }
    }
  }
}

template <typename T>
PixelMap<T> PixelMap<T>::downsize(int n){

  size_t nx = Nx/n,ny = Ny/n;

  PixelMap new_map(center,nx,ny,resolution*n);

  new_map.units = units;

  for(size_t i=0 ; i < nx ; ++i){
    for(size_t j=0 ; j < ny ; ++j){

      T &m = new_map(i,j);
      m = 0;

      for(size_t jj = j*n ; jj < (j+1)*n ; ++jj){
        for(size_t ii = i*n ; ii < (i+1)*n ; ++ii){
          m += map[ii + Nx*jj];
        }
      }

    }
  }

  return new_map;
}

template <typename T>
PixelMap<T> PixelMap<T>::interpolate(int n){
  size_t nx = n*(Nx-1)+1,ny = n*(Ny-1)+1;
  //size_t nx = n*Nx,ny = n*Ny;
  PixelMap new_map(center,nx,ny,resolution/n);
  new_map.units = units;
  
  long nn=n*n;
  //long off = (1+nx)*(n-1)/2;
  
  for(long jj = 0 ; jj<Ny-1 ; ++jj){
    long index = Nx*jj;
    for(long ii = 0 ; ii<Nx-1 ; ++ii){

      for(long j=0 ; j<n ; ++j){
        long k = ii*n + (jj*n + j )*nx;
        for(long i=0 ; i<n ; ++i){
       
          new_map[k] = (n-i)*(n-j)*map[index]
            + i*(n-j)*map[index+1]
            + i*j*map[index+1+Nx]
            + (n-i)*j*map[index+Nx];
          
          new_map[k] /= nn;
          ++k;
        }
      }
      ++index;
    }
  }
  {
    long jj=Ny-1 ;
    long index = Nx*jj;
    for(long ii = 0 ; ii<Nx-1 ; ++ii){

      long k = ii*n + (jj*n)*nx;
      for(long i=0 ; i<n ; ++i){
        
        new_map[k] = (n-i)*map[index]
            + i*map[index+1];
          
        new_map[k] /= n;
        ++k;
      }
      ++index;
    }
  }
  {
    long ii = Nx-1;
    for(long jj = 0 ; jj<Ny-1 ; ++jj){
      long index = ii + Nx*jj;

      long k = ii*n + (jj*n)*nx;
      for(long j=0 ; j<n ; ++j){
        
        new_map[k] = (n-j)*map[index]
          + j*map[index+Nx];
          
        new_map[k] /= n;
        k += nx;
      }
    }
  }
  
  new_map *= 1.0/nn;
  
  return new_map;
}

template <typename T>
PixelMap<T>& PixelMap<T>::operator=(const PixelMap &other)
{
  if(this != &other){
    PixelMap<T> copy(other);
    PixelMap<T>::swap(*this, copy);
  }
  return *this;
}
template <typename T>
PixelMap<T>& PixelMap<T>::operator=(PixelMap &&other)
{
  if(this != &other){
    PixelMap<T>::swap(*this, other);
  }
  return *this;
}

template <typename T>
void PixelMap<T>::swap(PixelMap<T> &map1,PixelMap<T> &map2)
{

  std::swap(map1.map,map2.map);
  std::swap(map1.Nx,map2.Nx);
  std::swap(map1.Ny,map2.Ny);
  std::swap(map1.resolution,map2.resolution);
  std::swap(map1.rangeX,map2.rangeX);
  std::swap(map1.rangeY,map2.rangeY);
  std::swap(map1.center[0],map2.center[0]);
  std::swap(map1.center[1],map2.center[1]);

  std::swap(map1.map_boundary_p1[0],map2.map_boundary_p1[0]);
  std::swap(map1.map_boundary_p1[1],map2.map_boundary_p1[1]);

  std::swap(map1.map_boundary_p2[0],map2.map_boundary_p2[0]);
  std::swap(map1.map_boundary_p2[1],map2.map_boundary_p2[1]);
  std::swap(map1.units,map2.units);

  return;
}

/// Multiplies the whole map by a scalar factor
template <typename T>
void PixelMap<T>::Renormalize(T factor)
{
  map *= factor;
}

/// Adds a value to the i-th pixel
template <typename T>
void PixelMap<T>::AddValue(std::size_t i, T value)
{
  map[i] += value;
}

/// Assigns a value to the i-th pixel
template <typename T>
void PixelMap<T>::AssignValue(std::size_t i, T value)
{
  map[i] = value;
}

template <typename T>
bool PixelMap<T>::agrees(const PixelMap<T>& other) const
{
  return
  (Nx == other.Nx) &&
  (Ny == other.Ny) &&
  (resolution == other.resolution) &&
  (center[0] == other.center[0]) &&
  (center[1] == other.center[1]) &&
  (units == other.units);
}

/// Add the values of another PixelMap to this one.
template <typename T>
PixelMap<T>& PixelMap<T>::operator+=(const PixelMap& rhs)
{
  if(Nx != rhs.getNx() || Ny != rhs.getNy())
    throw std::runtime_error("Dimensions of maps are not compatible");
  if(units != rhs.units)
    throw std::runtime_error("Units of maps are not compatible");

  for(size_t i=0;i<map.size();++i) map[i] += rhs.map[i];
  return *this;
}

/// Add two PixelMaps.
template <typename T>
PixelMap<T> PixelMap<T>::operator+(const PixelMap& a) const
{
  if(a.units != units)
    throw std::runtime_error("Units of maps are not compatible");
  PixelMap sum(a);
  sum += *this;
  return sum;
}

/// Subtract the values of another PixelMap from this one.
template <typename T>
PixelMap<T>& PixelMap<T>::operator-=(const PixelMap& rhs)
{
  if(Nx != rhs.getNx() || Ny != rhs.getNy())
    throw std::runtime_error("Dimensions of maps are not compatible");
  if(units != rhs.units)
    throw std::runtime_error("Units of maps are not compatible");
  for(size_t i=0;i<map.size();++i) map[i] -= rhs.map[i];
  return *this;
}

/// Subtract two PixelMaps.
template <typename T>
PixelMap<T> PixelMap<T>::operator-(const PixelMap& a) const
{
  if(units != a.units)
    throw std::runtime_error("Units of maps are not compatible");
  PixelMap diff(*this);
  diff -= a;
  return diff;
}

/// Multiply the values of another PixelMap by this one.
template <typename T>
PixelMap<T>& PixelMap<T>::operator*=(const PixelMap& rhs)
{
  if(Nx != rhs.getNx() || Ny != rhs.getNy())
    throw std::runtime_error("Dimensions of maps are not compatible");
  for(size_t i=0;i<map.size();++i) map[i] *= rhs.map[i];
  //map *= rhs.map;
  return *this;
}

// multiply each pixel by a constant
template <typename T>
PixelMap<T>& PixelMap<T>::operator*=(PosType b)
{
  for(size_t i=0;i<map.size();++i) map[i] *= b;
  //map *= rhs.map;
  return *this;
}

// multiply each pixel by a constant
template <typename T>
PixelMap<T> PixelMap<T>::operator*(PosType b) const
{
  PixelMap map(*this);
  map *= b;
  return map;
}

/// Multiply two PixelMaps.
template <typename T>
PixelMap<T> PixelMap<T>::operator*(const PixelMap& a) const
{
  PixelMap diff(a);
  diff *= *this;
  return diff;
}

/// Multiply two PixelMaps.
template <typename T>
PixelMap<T> PixelMap<T>::operator/(const PixelMap& a) const
{
  PixelMap diff(*this);
  for(size_t i=0;i<map.size();++i) diff[i] /= a.map[i];
  return diff;
}

/** \brief Add an image to the map
 *
 *  If rescale==0 gives constant surface brightness, if < 0
 *  the surface brightness is not scaled by the pixel area as for the flux (default: 1).
 *  Negative values are good for mapping some quantity independant of the pixel size
 */
template <typename T>
void PixelMap<T>::AddImages(
                         ImageInfo *imageinfo   /// An array of ImageInfo-s.  There is no reason to separate images for this routine
                         ,int Nimages           /// Number of images on input.
                         ,float rescale         /// rescales the surface brightness while leaving the image unchanged,
///  see full notes

){
  
  if(units != PixelMapUnits::surfb) throw std::invalid_argument("wrong units");
  if(Nimages <= 0) return;
  if(imageinfo->imagekist->Nunits() == 0) return;
  
  PosType sb = 1;
  float area = 1;
  std::list <unsigned long> neighborlist;
  std::list<unsigned long>::iterator it;
  for(long ii=0;ii<Nimages;++ii){
    
    if(imageinfo->imagekist->Nunits() > 0){
      imageinfo[ii].imagekist->MoveToTop();
      do{
        if(rescale != 0.0) sb = fabs(rescale)*imageinfo[ii].imagekist->getCurrent()->surface_brightness;
        
        assert(imageinfo[ii].imagekist->getCurrent()->leaf);
        
        if ((inMapBox(imageinfo[ii].imagekist->getCurrent()->leaf)) == true){
          PointsWithinLeaf(imageinfo[ii].imagekist->getCurrent()->leaf,neighborlist);
          for(it = neighborlist.begin();it != neighborlist.end();it++){
            area = LeafPixelArea(*it,imageinfo[ii].imagekist->getCurrent()->leaf);
            map[*it] += sb*area;
          }
        }
      }while(imageinfo[ii].imagekist->Down());
    }
  }
  
  if(rescale < 0){
    for(size_t i=0; i< Nx*Ny ;++i) map[i] /= resolution*resolution;
  }
  
  return;
}

template <typename T>
void PixelMap<T>::AddGridBrightness(Grid &grid){
  
  //if(units != photon_flux) throw std::invalid_argument("wrong units");
  PointList *plist = &(grid.i_tree->pointlist);
  
  if(plist->size() == 0) return;
  
  PointList::iterator listit= plist->Top();// = plist->begin();
  
  PosType sb = 1;
  float area = 1;
  
  std::list <unsigned long> neighborlist;
  std::list<unsigned long>::iterator it;
  //for(long ii=0;ii<Nimages;++ii){

  do{
    //for(listit = plist->begin() ; listit != plist->end() ; ++listit ){
    sb = (*listit)->surface_brightness;
    
    if (sb != 0.0 && (inMapBox((*listit)->leaf)) == true){
      PointsWithinLeaf((*listit)->leaf,neighborlist);
      for(it = neighborlist.begin();it != neighborlist.end();it++){
        area = LeafPixelArea(*it,(*listit)->leaf);
        map[*it] += sb*area;
      }
    }
    
  }while(--listit);
  
  return;
}


template <typename T>
void PixelMap<T>::AddGridMapBrightness(const GridMap &grid){
  
  //if(units != PixelMapUnits::surfb) throw std::invalid_argument("wrong units");
  try {
    // if GridMap res is an integer multiple of PixelMap res and they are aligned this will go
    grid.getPixelMapFlux(*this);
  } catch (const std::invalid_argument& ia) {
    // dimensions and/or alignment do not match
    PixelMap<T> newmap = grid.getPixelMapFlux<T>(1);
    copy_in(newmap);
  }
  
  return;
}

template <typename T>
void PixelMap<T>::AddImages(
                         std::vector<ImageInfo> &imageinfo   /// An array of ImageInfo-s.  There is no reason to separate images for this routine
                          ,int Nimages           /// Number of images on input.
                         ,float rescale         /// rescales the surface brightness while leaving the image unchanged,
///  see full notes
){
  AddImages(imageinfo.data(),Nimages,rescale);
}

template <typename T>
void PixelMap<T>::AddPointSource(const Point_2d &x,T flux){
  long index = find_index(x.x);
  if(index > -1) map[index] += flux;
}

/*
template <typename T>
void PixelMap<T>::AddImages(const GridMap &map){
  Point_2d x;
  for(size_t i=0; i < map.getNumberOfPoints(); ++i){
    Utilities::PositionFromIndex(i,x.x,map.getInitNgrid(),map.getXRange,center);
  }
}
*/


/** \brief Add images with uniform surface brightness set by input parameter value.
 *
 *   This does not use the surface brightnesses stored in the image points.
 */
template <typename T>
void PixelMap<T>::AddUniformImages(
                      ImageInfo *imageinfo   /// An array of ImageInfo-s.  There is no reason to separate images for this routine
                      ,int Nimages,T value){
  
  if(units != PixelMapUnits::surfb) throw std::invalid_argument("wrong units");
  if(Nimages <= 0) return;
  if(imageinfo->imagekist->Nunits() == 0) return;
  
  float area = 1;
  std::list <unsigned long> neighborlist;
  std::list<unsigned long>::iterator it;
  for(long ii=0;ii<Nimages;++ii){
    
    if(imageinfo->imagekist->Nunits() > 0){
      imageinfo[ii].imagekist->MoveToTop();
      do{
        
        assert(imageinfo[ii].imagekist->getCurrent()->leaf);
        
        if ((inMapBox(imageinfo[ii].imagekist->getCurrent()->leaf)) == true){
          PointsWithinLeaf(imageinfo[ii].imagekist->getCurrent()->leaf,neighborlist);
          for(it = neighborlist.begin();it != neighborlist.end();it++){
            area = LeafPixelArea(*it,imageinfo[ii].imagekist->getCurrent()->leaf);
            map[*it] += value*area/resolution/resolution;
          }
        }
      }while(imageinfo[ii].imagekist->Down());
    }
  }
  
  return;
}
/// returns the grid points within the branch
template <typename T>
void PixelMap<T>::PointsWithinLeaf(Branch * branch1, std::list <unsigned long> &neighborlist){
  
  neighborlist.clear();
  
  long line_s,line_e,col_s,col_e;
  
  find_index(branch1->boundary_p1,line_s,col_s);
  find_index(branch1->boundary_p2,line_e,col_e);
  
  if(line_s < 0) line_s = 0;
  if(col_s < 0) col_s = 0;
  
  //line_s = std::max(0,int(Utilities::IndexFromPosition(branch1->boundary_p1[0],Nx,range,center[0])));
  //col_s = std::max(0,int(Utilities::IndexFromPosition(branch1->boundary_p1[1],Ny,range*Ny/Nx,center[1])));
  
  //line_e = Utilities::IndexFromPosition(branch1->boundary_p2[0],Nx,range,center[0]);
  //col_e = Utilities::IndexFromPosition(branch1->boundary_p2[1],Ny,range*Ny/Nx,center[1]);
  
  if (line_e < 0) line_e = Nx-1;
  if (col_e < 0) col_e = Ny-1;
  
  for (int iy = col_s; iy<= col_e; ++iy)
  {
    for (int ix = line_s; ix <= line_e; ++ix)
    {
      neighborlist.push_back(ix+Nx*iy);
    }
    }
}
/// checks if the branch is within map boundaries
template <typename T>
bool PixelMap<T>::inMapBox(Branch * branch1) const{
  if (branch1->boundary_p1[0] > map_boundary_p2[0] || branch1->boundary_p2[0] < map_boundary_p1[0]) return false;
  if (branch1->boundary_p1[1] > map_boundary_p2[1] || branch1->boundary_p2[1] < map_boundary_p1[1]) return false;
  return true;
}
/// checks if point is within map boundaries
template <typename T>
bool PixelMap<T>::inMapBox(PosType * x) const{
  if (x[0] > map_boundary_p2[0] || x[0] < map_boundary_p1[0]) return false;
  if (x[1] > map_boundary_p2[1] || x[1] < map_boundary_p1[1]) return false;
  return true;
}

template <typename T>
bool PixelMap<T>::pixels_are_neighbors(size_t i,size_t j) const{
  
  long x = i%Nx - j%Nx;
  if(std::abs(x) > 1) return false;
  x = i/Nx - j/Nx;
  if(std::abs(x) > 1) return false;
  return true;
}
template <typename T>
void PixelMap<T>::find_contour(T level
                            ,std::vector<std::vector<Point_2d> > &points
                            ,std::vector<bool> &hits_edge
                            ) const {
  
  std::vector<bool> bitmap( map.size() );
  
  for(size_t i = 0 ; i<map.size() ; ++i) bitmap[i] = (map[i] > level);
  
  Utilities::find_boundaries<Point_2d>(bitmap,Nx,points,hits_edge,false);
  
  // rescale the points to PixelMap coordinates
  Point_2d xo(map_boundary_p1[0] + resolution*0.5
             ,map_boundary_p1[1] + resolution*0.5);
  for(std::vector<Point_2d> &v : points){
    for(Point_2d &p : v) p = p * resolution + xo;
  }
}

template <typename T>
void PixelMap<T>::find_islands_holes(T level,
                            std::vector<std::vector<size_t> > &points
                            ) const {
  
  std::vector<bool> bitmap( map.size() );
  std::vector<size_t> points_in;
  
  // excludes boundaries that will be set to false in Utilities::find_boundaries
  for(size_t i = 1 ; i<Nx-1 ; ++i){
    for(size_t j = 1 ; j<Ny-1 ; ++j){
      size_t k = i + Nx*j;
      if (map[k] > level){
        bitmap[k] = true;
        points_in.push_back(k);
      }else{
        bitmap[k] = false;
      }
    }
  }
  
  std::vector<bool> hits_edge;
  std::vector<std::vector<Point_2d> > boundaries;
  Utilities::find_boundaries<Point_2d>(bitmap,Nx,boundaries,hits_edge,false);
  points.resize(boundaries.size());
  
  if(boundaries.size() == 1){
    std::swap(points[0],points_in);
    return;
  }
  
  for(auto &v : points) v.clear();
  
  size_t n=points_in.size();
  size_t m=0;
  for(size_t k=0 ; k<n ; ++k){
    for(int i=0 ; i<boundaries.size() ; ++i){
      if( incurve(points_in[k],boundaries[i]) ){
        points[i].push_back(points_in[k]);
        ++m;
        break;
      }
    }
  }
  
  // remove holes
//  int i=0,k=points.size();
//  while( i < k){
//    m += points[i].size();
//    if(points[i].size() == 0){
//      std::swap(points[i],points[k-1]);
//      --k;
//    }else{
//      ++i;
//    }
//  }
//
//  points.resize(k);
  
  assert(m == n && "In PixelMap<T>::find_islands_holes");
}

/// find the index of the pixels that are larger than all its neighbors
template <typename T>
std::vector<size_t> PixelMap<T>::maxima(T minlevel
                            ) const {
  std::vector<size_t> indexes;
  if(map.max() < minlevel) return indexes;
  
  for(size_t j=1 ; j < Ny-1 ; ++j){
    size_t k = Nx*j+1;
    for(size_t i=1 ; i < Nx-1 ; ++i,++k){
      if(map[k] > minlevel){
        T val = map[k];
        if(
           val > map[k-1] &&
           val > map[k+1] &&
           val > map[k+Nx] &&
           val > map[k+Nx-1] &&
           val > map[k+Nx+1] &&
           val > map[k-Nx] &&
           val > map[k-Nx-1] &&
           val > map[k-Nx+1]
           ){
             indexes.push_back(k);
           }
      }
    }
  }
  return indexes;
}

template <typename T>
bool  PixelMap<T>::incurve(long k,std::vector<Point_2d> &curve) const{
  int n=0;
  long i = k % Nx , j = k / Nx;
  for(Point_2d &p : curve){
    if( p[0] > i && fabs(j - p[1]) < 0.1 ) ++n;
  }
  
  return n%2 == 1;
}

template <typename T>
void PixelMap<T>::lens_definition(
                            T min_sn_per_image
                            ,T pixel_threshold
                            ,int &Nimages
                            ,T &total_sig_noise_source
                            ,std::vector<size_t> &maxima_indexes
                            ,std::vector<std::vector<size_t> > &image_points
                            ,bool &lens_TF
                            ,T &level
                            ,size_t &n_pix_in_source
                            ,bool verbose
                            ){
    
  
  // default
  total_sig_noise_source = 0;
  maxima_indexes.clear();
  image_points.clear();
  lens_TF = false;
  level = 0;
  n_pix_in_source = 0;
  Nimages = 0;
  
  T sn_max = map.max();
  if(sn_max < pixel_threshold) return;
  
  find_islands_holes(pixel_threshold,image_points);
 
  total_sig_noise_source = 0;
  if(verbose) std::cout << "Initial Number of islands : " << image_points.size() << std::endl;
  std::vector<T> sig_noise(image_points.size(),0);
  for(size_t i=0 ; i<image_points.size() ; ++i ){
    for(size_t k : image_points[i] ){
      sig_noise[i] += map[k];
    }
    if(verbose) std::cout << "    signal-to-noise : " << sig_noise[i] << "  " << image_points[i].size() << std::endl;
    if(sig_noise[i] >= min_sn_per_image) total_sig_noise_source += sig_noise[i];
  }
  
  
  // remove holes

  bool ring = false;
  { // remove holes
    int i=0,k=image_points.size();
    while( i < k){
      if(image_points[i].size() == 0){
        std::swap(sig_noise[i],sig_noise[k-1]);
        std::swap(image_points[i],image_points[k-1]);
        --k;
        ring = true;
        if(verbose) std::cout << "There is a hole ! " << std::endl;
      }else{
        ++i;
      }
    }
    image_points.resize(k);
    sig_noise.resize(k);
  }
  
  {
    // remove low s/n images
    long i=0,k=sig_noise.size();
    while( i < k){
      if(sig_noise[i] < min_sn_per_image){
        std::swap(sig_noise[i],sig_noise[k-1]);
        std::swap(image_points[i],image_points[k-1]);
        --k;
      }else{
        ++i;
      }
    }
    image_points.resize(k);
    sig_noise.resize(k);
  }
  
  for(auto &v : image_points){
    for(size_t k : v){
      T val = map[k];
       if(
          val > map[k-1] &&
          val > map[k+1] &&
          val > map[k+Nx] &&
          val > map[k+Nx-1] &&
          val > map[k+Nx+1] &&
          val > map[k-Nx] &&
          val > map[k-Nx-1] &&
          val > map[k-Nx+1]
          ){
            maxima_indexes.push_back(k);
          }
    }
  }
  if(verbose) std::cout << "Number of maxima : " << maxima_indexes.size() << std::endl;
  
  
  n_pix_in_source = 0;
  for(auto &v : image_points) n_pix_in_source += v.size();
  if(verbose) std::cout << "              total : " << total_sig_noise_source << "  " << n_pix_in_source << std::endl;

  level = pixel_threshold;
  Nimages = image_points.size();
  
  if( image_points.size() == 1 && !ring){
    std::vector<std::vector<size_t> > new_image_points=image_points;
    
    while( new_image_points.size() == 1
          && new_image_points[0].size() > 1
          && !ring
          && level < sn_max
          ){
      
      std::sort(new_image_points[0].begin(),new_image_points[0].end()
                ,[this](size_t i,size_t j){return map[i] < map[j];});
      level = (map[ new_image_points[0][1] ] +  map[ new_image_points[0][0] ])/2;
      
      find_islands_holes(level,new_image_points);
      if(verbose) std::cout << "    Number of islands : " << new_image_points.size() << "   level " << level << std::endl;
      
      // detect holes and delete holes
      for(size_t i=0 ; i<new_image_points.size() ; ++i ){
        if(new_image_points[i].size() == 0){
          ring = true;
          for(long k = i ; k<new_image_points.size()-1 ; ++k){
            std::swap(new_image_points[k],new_image_points[k+1]);
          }
          new_image_points.pop_back();
        }
      }
      
      sig_noise.resize(new_image_points.size());
      for(size_t i=0 ; i<new_image_points.size() ; ++i ){
        sig_noise[i] = 0;
        for(size_t k : new_image_points[i] ){
          sig_noise[i] += map[k];
        }
        if(verbose) std::cout << "    signal-to-noise : " << sig_noise[i] << "  " << new_image_points[i].size() << std::endl;
      }

      {
        // remove low s/n images
        int i=0,k=sig_noise.size();
        while( i < k){
          if(sig_noise[i] < 2 * pixel_threshold ){
            std::swap(sig_noise[i],sig_noise[k-1]);
            std::swap(new_image_points[i],new_image_points[k-1]);
            --k;
          }else{
            ++i;
          }
        }

        sig_noise.resize(k);
        new_image_points.resize(k);
      }
      
    }
    // returns to one image in cases where peaks do not have enough S/N
    Nimages = MAX(new_image_points.size(),image_points.size());
  }
  
  lens_TF = Nimages > 1 || ring;
  if(verbose && lens_TF) std::cout << " IS OBSERVABLE LENS" << std::endl;
}

template <typename T>
int PixelMap<T>::count_islands(std::vector<size_t> &pixel_index) const{
  
  if(pixel_index.size() == 0) return 0;
  if(pixel_index.size() == 1) return 1;
  
  size_t *end = pixel_index.data() + pixel_index.size();
  size_t *current = pixel_index.data();
  int Ngroups = 1;
  size_t *group_boundary = current + 1;
  
  while(group_boundary != end){
    long ic = *current%Nx;
    long jc = *current/Nx;
    
    int Neighbors = 0;
    for(size_t *test = group_boundary
        ; test != end && Neighbors < 8 && group_boundary != end
        ; ++test
        ){
      long it = *test%Nx;
      long jt = *test/Nx;

      if( abs(it - ic) <= 1 && abs(jt - jc) <= 1  ){
        ++Neighbors;
        // swap test for group boundary
        size_t tmp = *group_boundary;
        *group_boundary = *test;
        *test = tmp;
        ++group_boundary;
      }
    }
    ++current;
    if(current == group_boundary && group_boundary != end ){
      ++group_boundary;
      ++Ngroups;
    }
  }
  return Ngroups;
}

  /*
  template <typename T>
int PixelMap<T>::count_islands(std::list<size_t> &pixel_index,std::vector<std::list<size_t>::iterator> &heads) const{
  
  heads.clear();
  if(pixel_index.size() == 0) return 0;
  
  if(pixel_index.size() == 1){
    heads.push_back(pixel_index.begin());
    return 1;
  }
  
  
  int ngroups = 0;
  pixel_index.sort();
  
  if(pixel_index.back() > Nx*Ny ){
    throw std::invalid_argument("index out of range");
  }
  
  size_t current;
  std::list<size_t>::iterator group = pixel_index.begin();
  
  while(group != pixel_index.end()){
    heads.push_back(group);
    current = *group;
    ++group;
    _count_islands_(current, pixel_index, group);
    
    ++ngroups;
  }
  
  heads.push_back(pixel_index.end());
  
  assert(ngroups == heads.size()-1);
  
  return ngroups;
}
*/

template <typename T>
void PixelMap<T>::_count_islands_(size_t current,std::list<size_t> &reservoir
                     ,std::list<size_t>::iterator &group) const{
  
  std::list<size_t>::iterator it = group;
  
  size_t imax = current + Nx + 1;  // maximum value of an index that can be a neighbor to current
  
  while( it != reservoir.end() && *it <= imax){
    if(pixels_are_neighbors(current,*it)){

      size_t tmp = *it;

      if(group == it){
        ++group;
      }else{
        reservoir.erase(it);
        reservoir.insert(group,tmp);
      }
      
      _count_islands_(tmp,reservoir,group);
      it = group;
      while(*it <= tmp && it != reservoir.end() ) ++it;  // skip forward to the next on in the list that hasn't been tested
    }else{
      ++it;
    }
  }
  return;
}


//// Finds the area of the intersection between pixel i and branch1
template <typename T>
PosType PixelMap<T>::LeafPixelArea(IndexType i,Branch * branch1){
  PosType area=0;
  PosType p[2],p1[2],p2[2];
  
  //Utilities::PositionFromIndex(i,p,Nx,range,center);
  find_position(p,i);
  p1[0] = p[0] - .5*resolution;
  p1[1] = p[1] - .5*resolution;
  p2[0] = p[0] + .5*resolution;
  p2[1] = p[1] + .5*resolution;
  area = MIN(p2[0],branch1->boundary_p2[0])
  - MAX(p1[0],branch1->boundary_p1[0]);
  if(area < 0) return 0.0;
  
  area *= MIN(p2[1],branch1->boundary_p2[1])
  - MAX(p1[1],branch1->boundary_p1[1]);
  if(area < 0) return 0.0;
  
  return area;
}

/// Print an ASCII table of all the pixel values.
template <typename T>
void PixelMap<T>::printASCII() const
{
  std::cout << Nx << " " << Ny << "  " << rangeX << std::endl;
  for(std::size_t i=0;i < map.size(); ++i) std::cout << map[i] << std::endl;
  std::cout << Nx << " " << Ny << "  " << rangeX << std::endl;
  
  //map.resize(0);
  return;
}
/// Print an ASCII table of all the pixel values.
template <typename T>
void PixelMap<T>::printASCIItoFile(std::string filename) const
{
  std::ofstream file_map(filename.c_str());
  
  if(!file_map){
    std::cout << "unable to open file " << filename << std::endl;
    exit(0);
  }
  
  std::cout << Nx << " " << Ny << "  " << rangeX << std::endl;
  for(std::size_t i=0;i < map.size(); ++i) file_map << std::scientific << map[i] << std::endl;
  std::cout << Nx << " " << Ny << "  " << rangeX << std::endl;
  
  //map.resize(0);
  
  file_map.close();
  
  return;
}
/// Output the pixel map as a fits file.
template <typename T>
void PixelMap<T>::printFITS(std::string filename,bool flipX,bool ctype,bool verbose)
{

  if(filename.empty())
    throw std::invalid_argument("Please enter a valid filename for the FITS file output");
  
  CPFITS_WRITE cpfits(filename,false);
  
  std::vector<long> naxex(2);
  naxex[0] = Nx;
  naxex[1] = Ny;

  if(flipX){
    std::valarray<T> map_inv(map.size());
    size_t s = 0;
    for(size_t i=0 ; i<Ny ; ++i){
      for(size_t j=1 ; j<=Nx ; ++j){
        map_inv[Nx-j + i*Nx] = map[s++];
      }
    }
    cpfits.write_image(map_inv,naxex);  // write the map
  }else{
    cpfits.write_image(map,naxex);  // write the map
  }
  cpfits.writeKey("WCSAXES", 2, "number of World Coordinate System axes");
  cpfits.writeKey("CRPIX1", 0.5*(naxex[0]+1), "x-coordinate of reference pixel");
  cpfits.writeKey("CRPIX2", 0.5*(naxex[1]+1), "y-coordinate of reference pixel");

  //cpfits.writeKey("CDELT1", 180*resolution/PI, "partial of first axis coordinate w.r.t. x");
  //cpfits.writeKey("CDELT2", 180*resolution/PI, "partial of second axis coordinate w.r.t. y");
  
  cpfits.writeKey("CROTA2", 0.0, "");
  cpfits.writeKey("CD1_1", -180*resolution/PI, "partial of first axis coordinate w.r.t. x (deg)");
  cpfits.writeKey("CD1_2", 0.0, "partial of first axis coordinate w.r.t. y");
  cpfits.writeKey("CD2_1", 0.0, "partial of second axis coordinate w.r.t. x");
  cpfits.writeKey("CD2_2", 180*resolution/PI, "partial of second axis coordinate w.r.t. y (deg)");
  
  cpfits.writeKey("Nx", Nx, "");
  cpfits.writeKey("Ny", Ny, "");
  cpfits.writeKey("range x", map_boundary_p2[0]-map_boundary_p1[0], "radians");

  cpfits.writeKey("RA_global", RA, "radians, center");
  cpfits.writeKey("DEC_global",DEC, "radians, center");
  cpfits.writeKey("center_x", center[0], "radians, center");
  cpfits.writeKey("center_y", center[1], "radians, center");
  
  cpfits.writeKey("CRVAL1", center[0]/degreesTOradians, "RA, degrees");
  cpfits.writeKey("CRVAL2", center[1]/degreesTOradians, "DEC, degrees");
  
  cpfits.writeKey("UNITS",to_string(units),"");
  
  if(ctype){
    cpfits.writeKey("CTYPE1", "RA---TAN", "the coordinate type for the first axis");
    cpfits.writeKey("CTYPE2", "DEC--TAN", "the coordinate type for the second axis");
    cpfits.writeKey("RADESYS", "ICRS", "");
    cpfits.writeKey("CUNIT1", "deg     ", "the coordinate unit for the first axis");
    cpfits.writeKey("CUNIT2", "deg     ", "the coordinate unit for the second axis");
  }
  
  for(auto &h : headers_float){
    cpfits.writeKey(std::get<0>(h),std::get<1>(h),std::get<2>(h));
  }
  for(auto &h : headers_long){
    cpfits.writeKey(std::get<0>(h),std::get<1>(h),std::get<2>(h));
  }
  for(auto &h : headers_string){
    cpfits.writeKey(std::get<0>(h),std::get<1>(h),std::get<2>(h));
  }
  
}

template <typename T>
void PixelMap<T>::printFITS(std::string filename
                         ,std::vector<std::tuple<std::string,double,std::string> > &extra_header_info, bool verbose)
{

  if(filename.empty())
    throw std::invalid_argument("Please enter a valid filename for the FITS file output");
  
  CPFITS_WRITE cpfits(filename,false);
  

  std::vector<long> naxex(2);
  naxex[0] = Nx;
  naxex[1] = Ny;

  cpfits.write_image(map,naxex);

  cpfits.writeKey("WCSAXES", 2, "number of World Coordinate System axes");
  cpfits.writeKey("CRPIX1", 0.5*(naxex[0]+1), "x-coordinate of reference pixel");
  cpfits.writeKey("CRPIX2", 0.5*(naxex[1]+1), "y-coordinate of reference pixel");
  cpfits.writeKey("CRVAL1", center[0]/degreesTOradians, "RA, degrees");
  cpfits.writeKey("CRVAL2", center[1]/degreesTOradians, "DEC, degrees");
  //cpfits.writeKey("CTYPE1", "RA---TAN", "the coordinate type for the first axis");
  //cpfits.writeKey("CTYPE2", "DEC--TAN", "the coordinate type for the second axis");
  //cpfits.writeKey("CUNIT1", "deg     ", "the coordinate unit for the first axis");
  //cpfits.writeKey("CUNIT2", "deg     ", "the coordinate unit for the second axis");
  //cpfits.writeKey("CDELT1", 180*resolution/PI, "partial of first axis coordinate w.r.t. x");
  //cpfits.writeKey("CDELT2", 180*resolution/PI, "partial of second axis coordinate w.r.t. y");
  cpfits.writeKey("CROTA2", 0.0, "");

  cpfits.writeKey("CD1_1", -180*resolution/PI, "partial of first axis coordinate w.r.t. x (deg)");
  //cpfits.writeKey("CD1_2", 0.0, "partial of first axis coordinate w.r.t. y");
  //cpfits.writeKey("CD2_1", 0.0, "partial of second axis coordinate w.r.t. x");
  cpfits.writeKey("CD2_2", 180*resolution/PI, "partial of second axis coordinate w.r.t. y (deg)");
  
  cpfits.writeKey("Nx", Nx, "");
  cpfits.writeKey("Ny", Ny, "");
  cpfits.writeKey("range x", map_boundary_p2[0]-map_boundary_p1[0], "radians");
  
  cpfits.writeKey("RA_global", RA, "radians, center");
  cpfits.writeKey("DEC_global",DEC, "radians, center");
  cpfits.writeKey("center_x", center[0], "radians, center");
  cpfits.writeKey("center_y", center[1], "radians, center");

  cpfits.writeKey("CRVAL1", center[0]/degreesTOradians, "RA, degrees");
  cpfits.writeKey("CRVAL2", center[1]/degreesTOradians, "DEC, degrees");

  cpfits.writeKey("UNITS",to_string(units),"");
 
  for(auto hp : extra_header_info){
    cpfits.writeKey(std::get<0>(hp),std::get<1>(hp),std::get<2>(hp));
  }
}
template <typename T>
void PixelMap<T>::printFITS(std::string filename
                         ,std::vector<std::string> &headercards)
{

  if(filename.empty())
    throw std::invalid_argument("Please enter a valid filename for the FITS file output");
  
  CPFITS_WRITE cpfits(filename,false);
  
  std::vector<long> naxex(2);
  naxex[0] = Nx;
  naxex[1] = Ny;

  cpfits.write_image(map,naxex);
  cpfits.writeHeader(headercards);
}

/**
 *
 * \brief Smoothes a map with a Gaussian kernel of width sigma (in arcseconds)
 */
template <typename T>
void PixelMap<T>::smooth(PosType sigma){
  PosType sum=0,**mask;
  int ix,iy;
  int Nmask,Nmask_half;
  int j_cen, k_cen;
  
  sigma *= arcsecTOradians;
  Nmask=2*(int)(3*sigma/resolution + 1);
  //::cout << Nmask << std::endl;
  if(Nmask < 4 ) std::cout << "WARNING: pixels are large compare to psf Nmask=" << Nmask << std::endl;
  
  Nmask_half = int(Nmask/2);
  mask = new PosType*[Nmask];
  for (int j = 0; j <Nmask; j++)
    mask[j] = new PosType[Nmask];
  
  for(int j=0;j<Nmask;j++)
  {
    for(int k=0;k<Nmask;k++)
    {
      j_cen = j - Nmask_half;
      k_cen = k - Nmask_half;
      mask[j][k]= exp(-(pow(j_cen*resolution,2) + pow(k_cen*resolution,2))/2/pow(sigma,2) );
      sum+=mask[j][k];
    }
  }
  for(int j=0;j<Nmask;j++)
  {
    for(int k=0;k<Nmask;k++)
    {
      mask[j][k]/=sum;
    }
  }
  
  std::valarray<T> map_out(0.0, map.size());
  
  for(long i=0;i<map.size();i++){
    for(int j=0;j<Nmask;j++){
      ix=i%Nx + j-Nmask_half;
      if( (ix>-1) && (ix<Nx) ){
        for(int k=0;k<Nmask;k++){
          iy=i/Nx + k-Nmask_half;
          if( (iy>-1) && (iy<Ny) ){
            map_out[ix+Nx*iy] += mask[j][k]*map[i];
          }
        }
      }
    }
  }
  
  std::swap(map,map_out);
  
  for (int j = 0; j <Nmask; j++)
    delete[] mask[j];
  delete[] mask;
}

/**
 * \brief Draws a line between two points on the image by setting
 * the pixels equal to value.
 *
 * TODO: Could be improved by detecting if the line passes through the map
 * at all before starting.  Could also be improved by making the line fatter
 * by including neighbor points.
 */
template <typename T>
void PixelMap<T>::drawline(
                        PosType x1[]     /// one end point of line
                        ,PosType x2[]    /// other end point of line
                        ,T value   /// value that it is set to on the map
                        ,bool add        /// true : add value, false replace with value
){
  
  //PosType x[2],s1,s2,r;
  //long index;
  //PosType d = 0;
  
  long index0 = find_index(x1);
  long index1 = find_index(x2);
  
  DrawLineGS(index0 % Nx,index1 % Nx,index0 / Nx,index1 / Nx,value,add);
  
//  r = sqrt( (x2[0] - x1[0])*(x2[0] - x1[0]) + (x2[1] - x1[1])*(x2[1] - x1[1]) );
//
//  if(r==0.0){
//    if(inMapBox(x1)){
//      //index = Utilities::IndexFromPosition(x1,Nx,range,center);
//      index = find_index(x1);
//      map[index] = value;
//    }
//    return;
//  }
//
//  s1 = (x2[0] - x1[0])/r;
//  s2 = (x2[1] - x1[1])/r;
//
//  x[0] = x1[0];
//  x[1] = x1[1];
//  while(d <= r){
//    if(inMapBox(x)){
//      //index = Utilities::IndexFromPosition(x,Nx,range,center);
//      index = find_index(x);
//      if(index != -1) map[index] = value;
//    }
//    x[0] += s1*resolution;
//    x[1] += s2*resolution;
//    d += resolution;
//  }
//
  return;
}

template <typename T>
void PixelMap<T>::DrawLine(long x0,long x1,long y0,long y1,T value,bool add) {
  
  x0 = MAX<long>(x0,0);
  x0 = MIN<long>(x0,Nx-1);
  x1 = MAX<long>(x1,0);
  x1 = MIN<long>(x1,Nx-1);
  
  y0 = MAX<long>(y0,0);
  y0 = MIN(y0,Ny-1);
  y1 = MAX<long>(y1,0);
  y1 = MIN<long>(y1,Ny-1);

  long dx = abs(x1 - x0);
  int sx = x0 < x1 ? 1 : -1;
  long dy = -abs(y1 - y0);
  int sy = y0 < y1 ? 1 : -1;
  long error = dx + dy;
   
  while(true){
    if(add){
      (*this)(x0,y0) += value;
    }else{
      (*this)(x0,y0) = value;
    }
    if(x0 == x1 && y0 == y1) break;
    long e2 = 2 * error;
    if(e2 >= dy){
      if(x0 == x1) break;
      error = error + dy;
      x0 = x0 + sx;
    }
    if(e2 <= dx){
      if(y0 == y1) break;
      error = error + dx;
      y0 = y0 + sy;
    }
  }
}


template <typename T>
void PixelMap<T>::DrawLineGS(long x0,long x1,long y0,long y1,T value,bool add) {

    long x = x0;
    long y = y0;
    long dx = x1-x0;
    long dy = y1-y0;
    long d = 2 * dy-dx; // discriminator
    
    // Euclidean distance of point (x,y) from line (signed)
    double D = 0;
    
    // Euclidean distance between points (x1, y1) and (x2, y2)
    double length = sqrt(dx * dx + dy * dy);
    
    double s = dx / length;
    double c = dy / length;
    while (x <= x1) {
      (*this)(x,MAX<long>(y-1,0)) = value ;//* ( abs(D + c) < 1);
      (*this)(x,y) = value * (D < 1);
      (*this)(x,MIN(y+1,Ny-1)) = value ;//* ( abs(D - c) < 1);
      
      x = x + 1;
      if (d <= 0) {
        D = D + s;
        d = d + 2 * dy;
      } else {
        D = D + s-c;
        d = d + 2 * (dy-dx);
        y = y + 1;
        
      }
    }
}

/**
 * \brief Draws a circle
 */
template <typename T>
void PixelMap<T>::drawcircle(
                          PosType r_center[]    /// center of circle
                          ,PosType radius       /// radius of circle
                          ,PosType value        /// value that it is set to on the map
){
  
  PosType x1[2],x2[2];
  PosType dtheta = resolution/fabs(radius);
  
  for(float theta = 0; theta < 2*PI; theta += dtheta){
    x1[0] = r_center[0] + radius*cos(theta);
    x1[1] = r_center[1] + radius*sin(theta);
    x2[0] = r_center[0] + radius*cos(theta+dtheta);
    x2[1] = r_center[1] + radius*sin(theta+dtheta);
    drawline(x1,x2,value,false);
  }
  
  return;
}

/**
 * \brief Draws a disk
 */
template <typename T>
void PixelMap<T>::drawdisk(
                          PosType r_center[]    /// center of disk
                          ,PosType radius       /// radius of disk
                          ,PosType value        /// value that it is set to on the map
                          ,int Nstrip           /// number of lines we want
){
  
  PosType x1[2],x2[2];
  
  // To do the circle (easy) :
  // group---------------------
  drawcircle(r_center,radius,value);
  
  // To fill the circle :
  // ------------------==
  
/*  for(float theta = 0; theta < 2*PI; theta += pi/N){
    x1[0] = r_center[0] - radius*cos(theta);
    x2[0] = r_center[0] + radius*cos(theta);
    x1[1] = x2[1] = r_center[1] + radius*sin(theta);
    drawline(x1,x2,value);
  }
  */
  for(float y = -radius + resolution/2 ; y <= radius; y += resolution){
    x1[0] = sqrt(radius*radius - y*y) + r_center[0];
    x2[0] = -sqrt(radius*radius - y*y) + r_center[0];
    x1[1] = x2[1] = r_center[1] + y;
    drawline(x1,x2,value,false);
  }
  return;
}

/**
 * \brief Draws a grid
 */
template <typename T>
void PixelMap<T>::drawgrid(int N,PosType value){
  
  PosType x1[2],x2[2];
  x1[1] = map_boundary_p1[1];
  x2[1] = map_boundary_p2[1];
  for(int i=1;i<N;++i){
    x1[0] = x2[0] = map_boundary_p1[0] + i*rangeX/N;
    drawline(x1,x2,value,false);
  }
  
  x1[0] = map_boundary_p1[0];
  x2[0] = map_boundary_p2[0];
  for(int i=1;i<N;++i){
    x1[1] = x2[1] = map_boundary_p1[1] + i*rangeY/N;
    drawline(x1,x2,value,false);
  }
}
template <typename T>
void PixelMap<T>::drawPoints(std::vector<Point *> points,PosType size,PosType value){
  if(size < resolution*3){
    size_t index;
    for(int i=0;i<points.size();++i){
      if(inMapBox(points[i]->x)){
        //index = Utilities::IndexFromPosition(x1,Nx,range,center);
        index = find_index(points[i]->x);
        map[index] = value;
      }
    }
  }else
    for(int i=0;i<points.size();++i) drawcircle(points[i]->x,0.01*rangeX,value);
  
}
template <typename T>
void PixelMap<T>::drawPoints(std::vector<Point> points,PosType size,PosType value){
  if(size < resolution*3){
    size_t index;
    for(int i=0;i<points.size();++i){
      if(inMapBox(points[i].x)){
        //index = Utilities::IndexFromPosition(x1,Nx,range,center);
        index = find_index(points[i].x);
        map[index] = value;
      }
    }
  }else
    for(int i=0;i<points.size();++i) drawcircle(points[i].x,0.01*rangeX,value);
  
}
template <typename T>
void PixelMap<T>::drawPoints(std::vector<Point_2d> points,PosType size,PosType value){
  if(size < resolution*3){
    size_t index;
    for(int i=0;i<points.size();++i){
      if(inMapBox(points[i].x)){
        //index = Utilities::IndexFromPosition(x1,Nx,range,center);
        index = find_index(points[i].x);
        map[index] = value;
      }
    }
  }else
    for(int i=0;i<points.size();++i) drawcircle(points[i].x,0.01*rangeX,value);
  
}
/**
 * \brief Draws a square
 */
template <typename T>
void PixelMap<T>::drawSquare(PosType p1[],PosType p2[],PosType value){
  PosType x1[2],x2[2];
  
  x1[0] = p1[0];
  x1[1] = p1[1];
  x2[0] = p2[0];
  x2[1] = p1[1];
  drawline(x1,x2,value,false);
  
  x1[0] = p2[0];
  x1[1] = p1[1];
  x2[0] = p2[0];
  x2[1] = p2[1];
  drawline(x1,x2,value,false);
  
  x1[0] = p2[0];
  x1[1] = p2[1];
  x2[0] = p1[0];
  x2[1] = p2[1];
  drawline(x1,x2,value,false);
  
  x1[0] = p1[0];
  x1[1] = p2[1];
  x2[0] = p1[0];
  x2[1] = p1[1];
  drawline(x1,x2,value,false);
}

/**
 * \brief Draws a box (filling the inside with horizontal lines, starting from the top)
 */
template <typename T>
void PixelMap<T>::drawBox(PosType p1[],PosType p2[],PosType value,int Nstrip)
{
  PosType x1ini[2],x2ini[2];
  PosType x1[2],x2[2];
  double N = double(Nstrip);
  
  // To do the frame (easy) :
  // ------------------------
  drawSquare(p1,p2,value);
  
  // To fill the square :
  // ------------------==

  // Initiating :
  if(p2[1]-p1[1]<0)
  {
    x1ini[0] = p1[0]; x1ini[1] = p1[1];
    x2ini[0] = p2[0]; x2ini[1] = p1[1];
    N *= -1. ;
  }
  else if(p2[1]-p1[1]>0)
  {
    x1ini[0] = p1[0]; x1ini[1] = p2[1];
    x2ini[0] = p2[0]; x2ini[1] = p2[1];
  }
  else
  {
    ERROR_MESSAGE();
    std::cout << "Error with drawbox." << std::endl;
    exit(0);
  }

  // Filling :
  x1[0] = x1ini[0] ;
  x2[0] = x2ini[0] ;
  for(int i=1;i<Nstrip;i++)
  {
    x1[1] = x1ini[1]-i*(p2[1]-p1[1])/N;
    x2[1] = x2ini[1]-i*(p2[1]-p1[1])/N;
    drawline(x1,x2,value,false);
  }

  return ;
}

/**
 * \brief Draws a closed curve through the points in curve->imagekist
 *
 * This differs form PixelMap<T>::AddImage() in that it draws lines between the points
 * and takes no account of the cells that the points are in or the surface brightness.
 * The points must be ordered already.  Particularly useful for drawing the caustics
 * that may have irregular cell sizes.  The last point will be connected to the first point.
 */
template <typename T>
void PixelMap<T>::AddCurve(ImageInfo *curve,T value){
  AddCurve(curve->imagekist,value);
  return;
}

template <typename T>
void PixelMap<T>::AddCurve(Kist<Point> *imagekist,T value){
  PosType x[2];
  
  if(imagekist->Nunits() == 0 ) return;
  
  imagekist->MoveToTop();
  x[0] = imagekist->getCurrent()->x[0];
  x[1] = imagekist->getCurrent()->x[1];
  imagekist->Down();
  for(;!(imagekist->OffBottom());imagekist->Down()){
    drawline(x,imagekist->getCurrent()->x,value,false);
    x[0] = imagekist->getCurrent()->x[0];
    x[1] = imagekist->getCurrent()->x[1];
  }
  imagekist->MoveToTop();
  drawline(x,imagekist->getCurrent()->x,value,false);
  
  return;
}
template <typename T>
void PixelMap<T>::AddCurve(std::vector<Point_2d> &curve,T value){
  PosType x[2];
  
  if(curve.size() == 0 ) return;
  
  x[0] = curve[0][0];
  x[1] = curve[0][1];
  for(size_t ii=1;ii<curve.size();++ii){
    drawline(x,curve[ii].x,value,false);
    x[0] = curve[ii][0];
    x[1] = curve[ii][1];
  }
  drawline(x,curve[0].x,value,false);
  
  return;
  
}
template <typename T>
void PixelMap<T>::AddCurve(std::vector<RAY> &curve,T value){
  PosType x[2];
  
  if(curve.size() == 0 ) return;
  
  x[0] = curve[0].x[0];
  x[1] = curve[0].x[1];
  for(size_t ii=1;ii<curve.size();++ii){
    drawline(x,curve[ii].x.x,value,false);
    x[0] = curve[ii].x[0];
    x[1] = curve[ii].x[1];
  }
  drawline(x,curve[0].x.x,value,false);
  
  return;
}

/**
 *  \brief Fills in pixels where the image plane points in the grid are located with the value given.
 
 This is for lensing quantities and not surface brightness.  If you want surface brightness use PixelMap<T>::AddGridBrightness()
 */
template <typename T>
void PixelMap<T>::AddGrid(const Grid &grid,T value){
  if(grid.getNumberOfPoints() == 0) return;

  PointList* list = &(grid.i_tree->pointlist);
  size_t index;
  
  PointList::iterator list_current = list->Top();
  do{
    if(inMapBox((*list_current)->x)){
      index = find_index((*list_current)->x);
      map[index] = value;
    }
  }while(--list_current);
  
}
/**
 *  \brief Fills in pixels with the selected quantity from the grid points.
 *
 *  The grid and PixelMap do not need to be related in any way.
 *  Using this function multiple grids can be added to the same image.
 *
 * This is for lensing quantities and not surface brightness.  If you want surface brightness use PixelMap<T>::AddGridBrightness()

 *
 *  Warning: When adding a new grid it should not overlap with any of the previously added grids.
 */
template <typename T>
void PixelMap<T>::AddGrid(const Grid &grid,LensingVariable val){
  
  if(grid.getNumberOfPoints() == 0 ) return;
  
  AddGrid_(grid.i_tree->pointlist,val);
  
  return;
  //***********************************************************************
  
  int Nblocks = 16;
  std::vector<PointList> lists(Nblocks);
  PointList::iterator list_current;
  
  Kist<Point> kist;
  
  bool allowDecent;
  //grid.i_tree->moveTop();
  TreeStruct::iterator treeit(grid.i_tree);
  int i = 0;
  do{
    if((*treeit)->level == 4){
      assert(i < 16);
      
      lists[i].setTop((*treeit)->points);
      lists[i].setN((*treeit)->npoints);
      list_current.current = lists[i].Top();
      list_current.JumpDownList( (*treeit)->npoints -1);
      lists[i].setBottom(*list_current);

      ++i;
      allowDecent = false;
    }else{
      allowDecent = true;
    }
    //  }while(grid.i_tree->TreeWalkStep(allowDecent) && i < Nblocks);
  }while(treeit.TreeWalkStep(allowDecent) && i < Nblocks);
  
   std::vector<std::thread> thr;
  for(int i = 0; i< Nblocks ;++i){
    thr.push_back(std::thread(&PixelMap<T>::AddGrid_,this,lists[i],val));
  }
  for(auto &t : thr) t.join();
}
template <typename T>
void PixelMap<T>::AddGrid_(const PointList &list,LensingVariable val){
  double tmp,area;
  PosType tmp2[2];
  KappaType tmp3[3];
 
  PointList::iterator pl_it = list.Top();
  do{
  //for(PointList::iterator pl_it = list.begin() ; pl_it != list.end() ; ++pl_it){
  //for(size_t i = 0; i< list.size(); ++i){
    
    switch (val) {
      case LensingVariable::ALPHA:
        tmp2[0] = (*pl_it)->x[0] - (*pl_it)->image->x[0];
        tmp2[1] = (*pl_it)->x[1] - (*pl_it)->image->x[1];
        tmp = sqrt(tmp2[0]*tmp2[0] + tmp2[1]*tmp2[1])/resolution/resolution;
        break;
      case LensingVariable::ALPHA1:
        tmp = ((*pl_it)->x[0] - (*pl_it)->image->x[0])/resolution/resolution;
        break;
      case LensingVariable::ALPHA2:
        tmp = ((*pl_it)->x[1] - (*pl_it)->image->x[1])/resolution/resolution;
        break;
      case LensingVariable::KAPPA:
        tmp = (*pl_it)->kappa()/resolution/resolution;
        break;
      case LensingVariable::GAMMA:
        tmp2[0] = (*pl_it)->gamma1();
        tmp2[1] = (*pl_it)->gamma2();
        tmp = sqrt(tmp3[0]*tmp3[0] + tmp3[1]*tmp3[1])/resolution/resolution;
        break;
      case LensingVariable::GAMMA1:
        tmp = (*pl_it)->gamma1()/resolution/resolution;
        break;
      case LensingVariable::GAMMA2:
        tmp = (*pl_it)->gamma2()/resolution/resolution;
        break;
      case LensingVariable::GAMMA3:
        tmp = (*pl_it)->gamma3()/resolution/resolution;
        break;
      case LensingVariable::INVMAG:
        tmp = (*pl_it)->invmag()/resolution/resolution;
        break;
      case LensingVariable::DELAYT:
        tmp = (*pl_it)->dt/resolution/resolution;
        break;
      case LensingVariable::SurfBrightness:
        tmp = (*pl_it)->surface_brightness;
        break;
      default:
        std::cerr << "PixelMap<T>::AddGrid() does not work for the input LensingVariable" << std::endl;
        throw std::runtime_error("PixelMap<T>::AddGrid() does not work for the input LensingVariable");
        break;
        // If this list is to be expanded to include ALPHA or GAMMA take care to add them as vectors
    }
    
    std::list <unsigned long> neighborlist;
    std::list <unsigned long>::iterator it;

    if(tmp != 0.0){
      if( inMapBox((*pl_it)->leaf) == true){
        PointsWithinLeaf((*pl_it)->leaf,neighborlist);
        for(it = neighborlist.begin();it != neighborlist.end();it++){
          area = LeafPixelArea(*it,(*pl_it)->leaf);
          map[*it] += tmp*area;
        }
      }
    }
    
  }while(--pl_it);
}



/// Find arcs in image  WARNING: THIS IS UNDER CONSTRUCTION!
template <typename T>
void PixelMap<T>::FindArc(
                       PosType &radius
                       ,PosType *xc
                       ,PosType *arc_center
                       ,PosType &arclength
                       ,PosType &width
                       ,PosType threshold    // threshold in pixal value
){
  
  if(Nx != Ny){
    std::cout << "PixelMap<T>::FindArc() Doesn't work on nonsquare maps" << std::endl;
    throw std::runtime_error("nonsquare");
  }
  std::vector<size_t> mask(Nx*Nx);
  size_t j=0;
  long k=0;
  PosType const tmp_center[2] = {0,0};
  
  // mask pixels below threshhold
  T maxval = map[0],minval = map[0];
  for(size_t i=0;i<Nx*Nx;i++){
    if(map[i] > threshold){
      if(j==0) minval = map[i];
      else minval = MIN(minval,map[i]);
      mask[j++]=i;
    }
    maxval = MAX(maxval,map[i]);
  }
  mask.resize(j);
  minval *= 0.99;
  
  if(j == 0 || j == Nx*Nx){
    std::cout << "PixelMap<T>::FindArc() - No pixels above surface brighness limit" << std::endl;
    radius = arclength = width = 0.0;
    xc[0] = xc[1] = 0.0;
    return;
  }
  
  PosType Rmax,Rmin,r2;
  Rmax = Nx;
  Rmin = 2;
  
  size_t Nc,Nr;
  
  Nc = (size_t)(2*Rmax);
  Nr = (size_t)(Rmax-Rmin)/2;
  
  std::vector<PosType> x(Nc),y(Nc),R2(Nr);
  Utilities::D3Matrix<float> votes(Nc,Nc,Nr);
  for(size_t i = 0;i<Nc*Nc*Nr;++i) votes(i) = 0;
  
  for(size_t i = 0;i<Nc;++i) x[i] = i*2*Rmax/(Nc-1) - Rmax;
  for(size_t i = 0;i<Nc;++i) y[i] = i*2*Rmax/(Nc-1) - Rmax;
  for(size_t i = 0;i<Nr;++i) R2[i] = pow(Rmin + i*(Rmax-Rmin)/(Nr-1),2);
  
  const PosType range = 1.0*Nx;
  PosType RmaxSqr = Rmax*Rmax;
  PosType rminSqr = RmaxSqr,rmax2=0;
  for(size_t m=0;m<mask.size();++m){
    
    Utilities::PositionFromIndex(mask[m], xc, Nx, range, tmp_center);
    
    for(size_t ii=0;ii<Nc;++ii){
      for(size_t jj=0;jj<Nc;++jj){
        
        r2 = (xc[0]-x[ii])*(xc[0]-x[ii]) + (xc[1]-y[jj])*(xc[1]-y[jj]);
        
        rminSqr = MIN(rminSqr,r2);
        rmax2 = MAX(rmax2,r2);
        
        if(r2 < RmaxSqr){
          k = Utilities::locate<PosType>(R2,r2);
          if(k > -1 && k < Nr){
            //votes(ii,jj,k) += log(map[mask[m]]/minval);
            votes(ii,jj,k) += map[mask[m]];
            //std::cout << "vote = " << votes(ii,jj,k) << std::endl;
          }
        }
      }
    }
  }
  
  printFITS("!fit_test.fits");
  
  // find maximum votes
  size_t kmax=0,ksecond=0;
  PosType maxvotes,secondvotes;
  maxvotes = votes(0);
  for(size_t kk=0;kk < Nc*Nc*Nr;++kk){
    
    if(votes(kk) >= maxvotes ){
      
      secondvotes = maxvotes;
      ksecond = kmax;
      
      maxvotes = votes(kk);
      kmax = kk;
    }
  }
  
  xc[0] = x[votes.xindex(kmax)];
  xc[1] = y[votes.yindex(kmax)];
  radius = sqrt(R2[votes.zindex(kmax)]);
  
  PosType x_tmp[2],r_tmp,rmax,rmin;
  double xave[2] = {0,0};
  rmax = 0.0;
  rmin = radius;
  arclength = 0;
  
  // find arc length, width and center
  for(size_t m=0;m<mask.size();++m){
    Utilities::PositionFromIndex(mask[m], x_tmp, Nx, range, tmp_center);
    r_tmp = sqrt( (xc[0] - x_tmp[0])*(xc[0] - x_tmp[0]) + (xc[1] - x_tmp[1])*(xc[1] - x_tmp[1]) );
    
    if(fabs(r_tmp-radius) < 1.0){
      arclength += 1;
      xave[0] += x_tmp[0];
      xave[1] += x_tmp[1];
    }
    rmax = MAX(rmax,r_tmp);
    rmin = MIN(rmin,r_tmp);
  }
  xave[0] /= arclength;
  xave[1] /= arclength;
  
  double tmp = sqrt( (xave[0] - xc[0])*(xave[0] - xc[0]) + (xave[1] - xc[1])*(xave[1] - xc[1]) );
  arc_center[0] = radius*(xave[0] - xc[0])/tmp + xc[0];
  arc_center[1] = radius*(xave[1] - xc[1])/tmp + xc[1];
  
  // convert from pixel units to angular units
  width = MAX(rmax-rmin,1.0)*resolution;
  arclength *= resolution;
  radius *= resolution;
  
  xc[0] = xc[0]*resolution + center[0];
  xc[1] = xc[1]*resolution + center[1];
  
  arc_center[0] = arc_center[0]*resolution + center[0];
  arc_center[1] = arc_center[1]*resolution + center[1];
  
}

/// get the index for a position, returns -1 if out of map
template <typename T>
long PixelMap<T>::find_index(PosType const x[],long &ix,long &iy) const{
  
  ix = (long)((x[0] - map_boundary_p1[0])/resolution);
  iy = (long)((x[1] - map_boundary_p1[1])/resolution);
  
  if( ix < 0 || ix >= Nx){
    ix = iy = -1;
    return -1;
  }
  if( iy < 0 || iy >= Ny){
    ix = iy = -1;
    return -1;
  }
  
  return ix + Nx*iy;
}
/// get the index for a position, returns -1 if out of map
template <typename T>
long PixelMap<T>::find_index(PosType x,PosType y,long &ix,long &iy) const{
  
  //ix = (long)((x - map_boundary_p1[0])/resolution + 0.5);
  //iy = (long)((y - map_boundary_p1[1])/resolution + 0.5);
  
  ix = (long)((x - map_boundary_p1[0])/resolution );
  iy = (long)((y - map_boundary_p1[1])/resolution );
  
  if( ix < 0 || ix >= Nx){
    ix = iy = -1;
    return -1;
  }
  if( iy < 0 || iy >= Ny){
    ix = iy = -1;
    return -1;
  }
  
  return ix + Nx*iy;
}
/// get the index for a position, returns -1 if out of map
template <typename T>
long PixelMap<T>::find_index(PosType const x[]) const{
  long ix,iy;
  return find_index(x,ix,iy);
}
/// get the index for a position, returns -1 if out of map
template <typename T>
long PixelMap<T>::find_index(PosType x,PosType y) const{
  long ix,iy;
  return find_index(x,y,ix,iy);
}
/// get the index for a position, returns -1 if out of map
template <typename T>
void PixelMap<T>::find_position(PosType x[],std::size_t const index) const{
  if(Nx == 1){
    x[0] = center[0];
    x[1] = center[1];
    return;
  }
  x[0] = map_boundary_p1[0] + resolution*( index%Nx + 0.5);
  x[1] = map_boundary_p1[1] + resolution*( index/Nx + 0.5);
  return;
}
/// get the index for a position, returns -1 if out of map
template <typename T>
void PixelMap<T>::find_position(PosType x[],std::size_t const ix,std::size_t const iy) const{
  if(Nx == 1){
    x[0] = center[0];
    x[1] = center[1];
    return;
  }
  //x[0] = center[0] + rangeX*( 1.0*ix/(Nx-1) - 0.5);
  //x[1] = center[1] + rangeY*( iy*1.0/(Ny-1) - 0.5);
  
  x[0] = map_boundary_p1[0] + resolution*(ix + 0.5);
  x[1] = map_boundary_p1[1] + resolution*(iy + 0.5);
  return;
}

template <typename T>
PixelMap<T> PixelMap<T>::rotate(PosType theta,T scale){

  double s=-sin(theta);
  double c=cos(theta);
  
  long center[2] = {Nx/2,Ny/2};
  Point_2d map_center = getCenter();
  double f[2];
  
  PixelMap rot_map(map_center.data(),Nx,Ny,resolution,units);
  
  size_t N = Nx*Ny;
  for(long k=0 ; k<N ; ++k){
    long i = k%Nx-center[0];
    long j = k/Nx-center[1];
    
    double x = (i*c - j*s)/scale + center[0];
    if(x>=0 && x < Nx-1){
      double y = (i*s + j*c)/scale + center[1];
      if(y>=0 && y < Ny-1){
        
        long kk = (long)(x) + (long)(y)*Nx; // lower left
        
        f[0]= x - (long)(x);
        f[1]= y - (long)(y);
        
        rot_map[k] = (1-f[0])*(1-f[1])*map[kk] + f[0]*(1-f[1])*map[kk+1] + f[0]*f[1]*map[kk+1+Nx]
        + (1-f[0])*f[1]*map[kk+Nx];
      }
    }
  }
  
  return rot_map;
}
template <typename T>
T PixelMap<T>::linear_interpolate(PosType x[]){
  long ix,iy;
  PosType f[2];
  long index;
  
  //f[0] = ((x[0] - center[0])/rangeX + 0.5)*(Nx-1);
  //f[1] = ((x[1] - center[1])/rangeY + 0.5)*(Ny-1);
  
  /*f[0] = ((x[0] - map_boundary_p1[0])/resolution + 0.5);
   f[1] = ((x[1] - map_boundary_p1[1])/resolution + 0.5);
   //std::cout << "(  " << fx << " " << fy << "   ";
   
   if (f[0] < 0. || f[0] > Nx-1){return 0;}
   else ix = (unsigned long)(f[0]);
   
   if (f[1] < 0. || f[1] > Ny-1){return 0;}
   else iy = (unsigned long)(f[1]);
   */
  
  ix = (long)((x[0] - map_boundary_p1[0])/resolution - 0.5);
  iy = (long)((x[1] - map_boundary_p1[1])/resolution - 0.5);
  
  if(ix < 0 || iy < 0 || ix > Nx-1 || iy > Ny-1) return 0;
  
  if(ix == Nx-1) ix = Nx-2;
  if(iy == Ny-1) iy = Ny-2;
  
  // index of nearest grid point to the lower left
  index = ix + Nx*iy;
  
  find_position(f,index);
  
  /** bilinear interpolation */
  f[0]=(x[0] - f[0])/resolution;
  f[1]=(x[1] - f[1])/resolution;
  
  assert(f[0] > 0 || ix == 0);
  assert(f[1] > 0 || iy == 0);
  //assert(f[0] <= 1.0 && f[1] <= 1.0);
  
  return (1-f[0])*(1-f[1])*map[index] + f[0]*(1-f[1])*map[index+1] + f[0]*f[1]*map[index+1+Nx]
  + (1-f[0])*f[1]*map[index+Nx];
}

/*// get the index for a position, returns -1 if out of map
 long PixelMap<T>::find_index(PosType const x[],long &ix,long &iy){
 PosType fx, fy;
 
 fx = (x[0] - center[0])/rangeX;
 fy = (x[1] - center[1])/rangeY;
 if(fabs(fx) > 0.5 || fabs(fy) > 0.5){
 ix = iy = -1;
 return -1;
 }
 
 fx = (fx + 0.5)*(Nx-1) + 0.5;
 fy = (fy + 0.5)*(Ny-1) + 0.5;
 
 ix = (long)(fx);
 iy = (long)(fy);
 
 if( (ix<Nx) && (iy<Ny) ) return ix+Nx*iy;
 return -1;
 }
 */


template <typename T>
PosType PixelMap<T>::AddSource(Source &source){
  //if(units != PixelMapUnits::surfb) throw std::invalid_argument("wrong units");
  Point_2d s_center;
  source.getTheta(s_center);
  
  if( (s_center[0] + source.getRadius()) < map_boundary_p1[0] ) return 0.0;
  if( (s_center[0] - source.getRadius()) > map_boundary_p2[0] ) return 0.0;
  if( (s_center[1] + source.getRadius()) < map_boundary_p1[1] ) return 0.0;
  if( (s_center[1] - source.getRadius()) > map_boundary_p2[1] ) return 0.0;

  PosType y[2];
  PosType tmp = resolution*resolution;
  PosType total = 0;
  
  double sb;
  for(size_t index =0 ;index < map.size(); ++index){
    find_position(y,index);
    sb = source.SurfaceBrightness(y);
    map[index] += sb*tmp;
    total += sb*tmp;
  }
  
  return total;
}

template <typename T>
PosType PixelMap<T>::AddSource(Source &source,int oversample){
  if(units != PixelMapUnits::surfb) throw std::invalid_argument("wrong units");

  Point_2d s_center;
  source.getTheta(s_center);
  
  if( (s_center[0] + source.getRadius()) < map_boundary_p1[0] ) return 0.0;
  if( (s_center[0] - source.getRadius()) > map_boundary_p2[0] ) return 0.0;
  if( (s_center[1] + source.getRadius()) < map_boundary_p1[1] ) return 0.0;
  if( (s_center[1] - source.getRadius()) > map_boundary_p2[1] ) return 0.0;

  PosType y[2],x[2],bl;
  PosType tmp_res = resolution*1.0/oversample;
  PosType tmp = tmp_res*tmp_res;
  PosType total = 0.0;
  
  bl = resolution /2 - 0.5*tmp_res;
  
  for(size_t index =0 ;index < map.size(); ++index){
    find_position(y,index);
    y[0] -= bl;
    y[1] -= bl;
    for(int i = 0 ; i < oversample ; ++i){
      x[0] = y[0] + i*tmp_res;
      for(int j=0; j < oversample;++j){
        x[1] = y[1] + j*tmp_res;
        map[index] += source.SurfaceBrightness(x)*tmp;
        total += source.SurfaceBrightness(x)*tmp;
      }
    }
  }
  return total;
}
template <typename T>
void PixelMap<T>::duplicate(
                       const PixelMap& pmap
  ){
  
  if(!agrees(pmap)){
    throw std::invalid_argument("Maps not the same");
  }
  
  for(size_t i=0;i<map.size();++i) map[i] = pmap.map[i];
  return;
}
template <typename T>
void PixelMap<T>::copy_in(
                   const PixelMap& pmap
)
{
  
  if(agrees(pmap)){  // maps are the same dimensions and position
    for(size_t i=0;i<map.size();++i) map[i] += pmap.map[i];
    return;
  }
  double res_ratio = resolution / pmap.resolution;
  if(abs(res_ratio -1) < 1.0e-4) res_ratio = 1.0;
  double res_ratio2 = res_ratio*res_ratio;
  
  // check is maps overlap
  if(map_boundary_p1[0] > pmap.map_boundary_p2[0] ) return;
  if(map_boundary_p2[0] < pmap.map_boundary_p1[0] ) return;
  if(map_boundary_p1[1] > pmap.map_boundary_p2[1] ) return;
  if(map_boundary_p2[1] < pmap.map_boundary_p1[1] ) return;

  
  double halfpixel = res_ratio/2;
  PosType x[2];
  size_t NNx = pmap.getNx(),NNy = pmap.getNy();
  //size_t Npmap = NNx*NNy;
  
  for(size_t ii=0 ; ii < map.size() ; ++ii){
    
    find_position(x,ii);
      // find range if this pixel in pmap's pixel space
    double ix = (x[0] - pmap.map_boundary_p1[0])/pmap.resolution;
    double iy = (x[1] - pmap.map_boundary_p1[1])/pmap.resolution;
    
    double xmin = MAX(0.0,ix - halfpixel);
    double xmax = MIN<double>(NNx,ix + halfpixel);
    if(xmin >= xmax) continue;

    double ymin = MAX(0.0,iy - halfpixel);
    double ymax = MIN<double>(NNy,iy + halfpixel);
    if(ymin >= ymax) continue;
    
    long imin = MIN<long>((long)(xmin),NNx-1);
    long imax = MIN<long>((long)(xmax),NNx-1);
    long jmin = MIN<long>((long)(ymin),NNy-1);
    long jmax = MIN<long>((long)(ymax),NNy-1);
    long jj;
    
    for(long j = jmin ; j <= jmax ; ++j ){
      double area1 = MIN<double>(ymax,j+1) -  MAX<double>(ymin,j);
      if(area1 <= 0.0) continue;
      jj = NNx*j;
      for(long i = imin ; i <= imax ; ++i ){
        double area = (MIN<double>(xmax,i+1) -  MAX<double>(xmin,i)) * area1 ;
        map[ii] += pmap.map[i + jj]*area/res_ratio2;
      }
    }
  }
}
template <typename T>
void PixelMap<T>::paste(const PixelMap<T>& pmap){
  
  if(resolution < pmap.resolution){
    std::cerr << "PixelMap<T>::paste() resolution of image pasted in must of equal or higher resolution" << std::endl;
    std::cerr << resolution << " " << pmap.resolution << " dres/res " << (pmap.resolution-pmap.resolution)/pmap.resolution << std::endl;
    throw std::invalid_argument("low resolution");
  }
  
  // case where the maps do not overlap
  if( (map_boundary_p1[0] > pmap.map_boundary_p2[0]) || (pmap.map_boundary_p1[0] > map_boundary_p2[0])
     || (map_boundary_p1[1] > pmap.map_boundary_p2[1]) || (pmap.map_boundary_p1[1] > map_boundary_p2[1])
     ){
    return;
  }
  
  double x[2];
  
  if(map.size() < pmap.map.size() ){
    for(size_t i=0 ; i < map.size() ; ++i ){
      find_position(x,i);
      long j = pmap.find_index(x[0], x[1]);
      if(j >= 0 ){
        map[i] = pmap(j);
      }
    }
  }else{
    for(size_t j=0 ; j < pmap.map.size() ; ++j ){
      pmap.find_position(x,j);
      long i = find_index(x[0], x[1]);
      if(i >= 0 ){
        map[i] = pmap(j);
      }
    }
  }
}

template <typename T>
void PixelMap<T>::paste(const PixelMap& pmap
                     ,long nx_ll    // lower left pixel of this
                     ,long ny_ll
                     ){
  
  if(resolution != pmap.resolution){
    std::cerr << "PixelMap<T>::paste() resolution of image pasted in must of equal or higher resolution" << std::endl;
    std::cerr << resolution << " " << pmap.resolution << " dres/res " << (pmap.resolution-pmap.resolution)/pmap.resolution << std::endl;
    throw std::invalid_argument("resolution");
  }
  
  if((nx_ll > Nx-1) || (ny_ll > Ny-1) ) return;
  
  long nx_ur = nx_ll + pmap.getNx();
  long ny_ur = ny_ll + pmap.getNy();
  
  if((nx_ur < 0) || (ny_ur < 0) ) return;

  if(nx_ur > Nx) nx_ur = Nx;
  if(ny_ur > Ny) ny_ur = Ny;
 
  for(long j = MAX<long>(ny_ll,0) ; j<ny_ur ; ++j ){
    for(long i= MAX<long>(nx_ll,0) ; i<nx_ur ; ++i ){
      map[i+j*Nx] += pmap(i-nx_ll,j-ny_ll);
    }
  }
}

template <typename T>
PixelMap<T> PixelMap<T>::convolve(const PixelMap<T>& kernel){
  long nx = kernel.getNx();
  long ny = kernel.getNy();
  
  //long dx = nx/2+1;
  //long dy = ny/2+1;
  
  long dx = nx/2;
  long dy = ny/2;
  
  PixelMap<T> copy(center,Nx,Ny,resolution);
  
  for(long ii = 0 ; ii< Nx ; ++ii){
    long imin = (ii > dx) ? 0 : dx - ii;
    long imax = (Nx-ii > dx) ? nx-1 : Nx-ii-1+dx ;

    for(long jj = 0 ; jj<Ny ; ++jj){
      long jmin = (jj > dy) ? 0 : dy - jj ;
      long jmax = (Ny-jj > dx) ? ny-1 : Ny-jj-1+dy ;

      T &a = copy[ii + Nx*jj];
      a=0;
      
      for(long i=imin  ; i<=imax ; ++i){
        for(long j=jmin ; j<=jmax ; ++j){
          a += kernel[ i + nx*j ] * map[ ii + i - dx + Nx*(jj + j - dy) ];
        }
      }
    }
  }
  
  copy.units = units;
  
  return copy;
}

template <typename T>
PixelMap<T> PixelMap<T>::convolve2(const PixelMap<T>& kernel){
  long nx = kernel.getNx();
  long ny = kernel.getNy();
  
  long dx = nx/2;
  long dy = ny/2;
  
  PixelMap<T> copy(center,Nx,Ny,resolution);
  
  for(long ii = 0 ; ii< Nx ; ++ii){
    for(long jj = 0 ; jj<Ny ; ++jj){

      T &a = copy[ii + Nx*jj];
      a=0;
      
      for(long i=0  ; i<nx ; ++i){
        long k = ii + i - dx;
        if(k > -1 && k < Nx){
          for(long j=0 ; j<ny ; ++j){
            long kk = jj + j - dy;
            if(kk > -1 && kk < Ny ) a += kernel[ i + nx*j ] * map[ k + Nx*kk ];
          }
        }
      }
    }
  }
  
  copy.units = units;
  
  return copy;
}

template <typename T>
PixelMap<T> PixelMap<T>::cutout(long xmin,long xmax,long ymin,long ymax){
  long nx = xmax-xmin;
  long ny = ymax-ymin;
  
  PixelMap<T> copy(center,nx,ny,resolution);
  copy.units = units;
  
  for(long i=0  ; i<nx ; ++i){
    for(long j=0 ; j<ny ; ++j){
      copy[i + nx*j] = map[ (xmin+i) + Nx*(ymin+j) ];
    }
  }
  
  return copy;
}


template <typename T>
void PixelMap<T>::recenter(PosType c[2] /// new center
){
  double dc[2];
  dc[0] = c[0] - center[0];
  dc[1] = c[1] - center[1];
  
  map_boundary_p1[0] = map_boundary_p1[0] + dc[0];
  map_boundary_p1[1] = map_boundary_p1[1] + dc[1];
  map_boundary_p2[0] = map_boundary_p2[0] + dc[0];
  map_boundary_p2[1] = map_boundary_p2[1] + dc[1];
  
  center[0] = c[0];
  center[1] = c[1];
  
  return;
}
template <typename T>
void PixelMap<T>::recenter(Point_2d c /// new center
){
  recenter(c.x);
}

template <typename T>
void PixelMap<T>::PowerSpectrum(std::vector<PosType> &power_spectrum   /// output power spectrum
                     ,std::vector<PosType> &lvec            /// output l values of bands
                     ,bool overwrite                /// if false add power to existing power_spectrum (used for averaging over many fields
                     ){
                       
    power_spectrum.resize(lvec.size());
    
    if(overwrite){
      Utilities::powerspectrum2d(map,Nx,Ny,rangeX,rangeY,lvec,power_spectrum);
    }else{
      std::vector<PosType> tmp_power(power_spectrum.size());
      Utilities::powerspectrum2d(map,Nx,Ny,rangeX,rangeY,lvec,tmp_power);
      for(size_t ii=0;ii<power_spectrum.size();++ii) power_spectrum[ii] += tmp_power[ii];
    }
}

template <typename T>
void PixelMap<T>::PowerSpectrum(std::vector<PosType> &power_spectrum   /// output power spectrum
                     ,const std::vector<PosType> &lbins            /// input l values of bands
                     ,std::vector<PosType> &lave            /// output l values of bands
                     ,bool overwrite                /// if false add power to existing power_spectrum (used for averaging over many fields
                     ){
    
    if(overwrite){
      Utilities::powerspectrum2dprebin(map,Nx,Ny,rangeX,rangeY,lbins,power_spectrum,lave);
    }else{
      if(power_spectrum.size() != lbins.size()-1) throw std::invalid_argument("these must be the same size");
      std::vector<PosType> tmp_power(power_spectrum.size());
      Utilities::powerspectrum2dprebin(map,Nx,Ny,rangeX,rangeY,lbins,tmp_power,lave);
      for(size_t ii=0;ii<power_spectrum.size();++ii) power_spectrum[ii] += tmp_power[ii];
    }
  }

template <typename T>
void PixelMap<T>::AdaptiveSmooth(PosType value){
    std::valarray<T> tmp = Utilities::AdaptiveSmooth<T>(data(),Nx,Ny,value);
    map = tmp;
  }

/*
template <typename T>
void PixelMap<T>::AddSource(Source &source,int oversample){
  Point_2d s_center;
  source.getTheta(s_center);
  
  if( s_center[0] + source.getRadius() < map_boundary_p1[0] ) return;
  if( s_center[0] - source.getRadius() > map_boundary_p2[0] ) return;
  if( s_center[1] + source.getRadius() < map_boundary_p1[1] ) return;
  if( s_center[1] - source.getRadius() > map_boundary_p2[1] ) return;
  
  PosType tmp_res = resolution*1.0/oversample;
  PosType tmp = resolution;
  
  Point_2d dx,y;
  PosType r = source.getRadius();
  
  map[find_index(s_center[0],s_center[1])] += source.SurfaceBrightness(s_center.x)*tmp;
  
  dx[0] = 0;
  for(dx[1] = r/oversample ; dx[1] <= r; dx[1] += r/oversample){
      y = s_center + dx;
      map[find_index(y[0],y[1])] += source.SurfaceBrightness(y.x)*tmp;
      y = s_center - dx;
      map[find_index(y[0],y[1])] += source.SurfaceBrightness(y.x)*tmp;
  }
  
  for(dx[0] = r/oversample; dx[0] <= source.getRadius() ; dx[0] += r/oversample ){
    
    PosType range = sqrt(r*r - dx[0]*dx[0]);
    for(dx[1] = 0 ; dx[1] <= range; dx[1] += r/oversample){
      y = s_center + dx;
      map[find_index(y[0],y[1])] += source.SurfaceBrightness(y.x)*tmp;
      y = s_center - dx;
      map[find_index(y[0],y[1])] += source.SurfaceBrightness(y.x)*tmp;
    }
  }
  
  for(dx[0] = -r/oversample; dx[0] >= -source.getRadius() ; dx[0] -= r/oversample ){
    PosType range = sqrt(r*r - dx[0]*dx[0]);
    for(dx[1] = r/oversample ; dx[1] <= range; dx[1] += r/oversample){
      y = s_center + dx;
      map[find_index(y[0],y[1])] += source.SurfaceBrightness(y.x)*tmp;
      y = s_center - dx;
      map[find_index(y[0],y[1])] += source.SurfaceBrightness(y.x)*tmp;
    }
  }
}
*/

/* \brief convolve the image with a kernel.
 
 It is assumed that the size of the kernel is much smaller than the image and
 that the kernal has the same pixel size as the image.
 **
template <typename T>
void PixelMap<T>::convolve(PixelMap<T> &kernel,long center_x,long center_y){
  std::valarray<T> output(Nx*Ny);
  
  //std::cout << output.size() << " " << map.size() << std::endl;
  
  if( center_x == 0 ){
    center_x = kernel.getNx()/2;
    center_y = kernel.getNy()/2;
  }
  
  size_t Nxk = kernel.getNx();
  size_t Nyk = kernel.getNy();
  
  //std::cout << "sum of map : " << sum() << std::endl;
  //std::cout << "sum of kernel : " << kernel.sum() << std::endl;
  //double total=0;
  for(size_t k1 = 0 ; k1 < Nx ; ++k1){
    long bl1 = k1 - center_x;
    
    for(size_t k2 = 0 ; k2 < Ny ; ++k2){
      long bl2 = k2 - center_y;
      
      long k = k1 + Nx * k2;
      output[k] = 0;
      
      for(size_t j = 0 ; j < Nyk ; ++j){
        long kk2 = bl2 + j;
        if(kk2 < 0 || kk2 >= Ny) continue;
        
        for(size_t i = 0; i < Nxk ; ++i){
          long kk1 = bl1 + i;
          if(kk1 < 0 || kk1 >= Nx) continue;
          
          output[k] += map[ kk1 + Nx * kk2 ] * kernel[i + Nxk * j];
        }
      }
      //total += output[k];
    }
  }
  
  std::swap(map,output);
}*/

template class PixelMap<double>;
template class PixelMap<float>;

