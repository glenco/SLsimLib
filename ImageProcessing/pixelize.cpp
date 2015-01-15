/*
 * pixelize.c
 *
 *  Created on: Feb 27, 2010
 *      Author: R.B. Metcalf
 */

#include "slsimlib.h"

#ifdef ENABLE_FITS
#include <CCfits/CCfits>
//#include <CCfits>
#endif

#include <fstream>
#include <algorithm>
#include <utility>
#include <stdexcept>
#include <thread>
#include "image_processing.h"
#include "point.h"

#if __cplusplus < 201103L
template<typename T>
void swap(std::valarray<T>& x, std::valarray<T>& y)
{
  std::valarray<T> z(x);
  
  x.resize(y.size());
  x = y;
  
  y.resize(z.size());
  y = z;
}
#endif

void swap(PixelMap& x, PixelMap& y)
{
  using std::swap;
  
  swap(x.map,y.map);
  
  swap(x.Nx, y.Nx);
  swap(x.Ny, y.Ny);
  swap(x.resolution, y.resolution);
  swap(x.rangeX, y.rangeX);
  swap(x.rangeY, y.rangeY);
  
  swap(x.center[0], y.center[0]);
  swap(x.center[1], y.center[1]);
  
  swap(x.map_boundary_p1[0], y.map_boundary_p1[0]);
  swap(x.map_boundary_p1[1], y.map_boundary_p1[1]);
  swap(x.map_boundary_p2[0], y.map_boundary_p2[0]);
  swap(x.map_boundary_p2[1], y.map_boundary_p2[1]);
}

PixelMap::PixelMap()
: map(), Nx(0), Ny(0), resolution(0), rangeX(0), rangeY(0)
{
  center[0] = 0;
  center[1] = 0;
  
  map_boundary_p1[0] = 0;
  map_boundary_p1[1] = 0;
  map_boundary_p2[0] = 0;
  map_boundary_p2[1] = 0;
}

PixelMap::PixelMap(const PixelMap& other)
: map(other.map),
Nx(other.Nx), Ny(other.Ny), resolution(other.resolution), rangeX(other.rangeX), rangeY(other.rangeY)
{
  std::copy(other.center, other.center + 2, center);
  
  std::copy(other.map_boundary_p1, other.map_boundary_p1 + 2, map_boundary_p1);
  std::copy(other.map_boundary_p2, other.map_boundary_p2 + 2, map_boundary_p2);
}

/// make square PixelMap
PixelMap::PixelMap(
                   const PosType* center,  /// The location of the center of the map
                   std::size_t Npixels,  /// Number of pixels in one dimension of map.
                   PosType resolution        /// One dimensional range of map in whatever units the point positions are in
)
: map(0.0, Npixels*Npixels),
Nx(Npixels), Ny(Npixels), resolution(resolution)
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
PixelMap::PixelMap(
                   const PosType* center,  /// The location of the center of the map
                   std::size_t Nx,  /// Number of pixels in x dimension of map.
                   std::size_t Ny,  /// Number of pixels in y dimension of map.
                   PosType resolution        /// One dimensional range of map in whatever units the point positions are in
)
: map(0.0, Nx*Ny),
Nx(Nx), Ny(Ny), resolution(resolution)
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
PixelMap::PixelMap(
                   std::string fitsfilename   /// file name of fits file to be read
                   ,double my_res         /// resolution (rad) of fits image if not given in fits file, use default or -1 otherwise
)
{
#ifdef ENABLE_FITS
  if(fitsfilename.empty())
    throw std::invalid_argument("Please enter a valid filename for the FITS file input");
  
  std::auto_ptr<CCfits::FITS> fp(new CCfits::FITS(fitsfilename, CCfits::Read));
  
  CCfits::PHDU& h0 = fp->pHDU();
  
  //const CCfits::ExtMap *h1=&fp->extension();
  
  Nx = h0.axis(0);
  Ny = h0.axis(1);
  
  try
  {
    h0.readKey("CRVAL1", center[0]);
    h0.readKey("CRVAL2", center[1]);
  }
  catch(CCfits::HDU::NoSuchKeyword&)
  {
    center[0] = 0.0;
    center[1] = 0.0;
  }
  
  if(my_res == -1){
    // read the resolution
    try
    {
      double cdelt2;
      h0.readKey("CDELT1", my_res);
      h0.readKey("CDELT2", cdelt2);
      if(std::abs(my_res) - std::abs(cdelt2) > 1e-6)
        throw std::runtime_error("non-square pixels in FITS file " + fitsfilename);
    }
    catch(CCfits::HDU::NoSuchKeyword&)
    {
      try{
        double cd12, cd21, cd22;
        h0.readKey("CD1_1", my_res);
        h0.readKey("CD1_2", cd12);
        h0.readKey("CD2_1", cd21);
        h0.readKey("CD2_2", cd22);
        if(std::abs(my_res) - std::abs(cd22) > 1e-6)
          throw std::runtime_error("non-square pixels in FITS file " + fitsfilename);
        if(cd12 || cd21)
          throw std::runtime_error("pixels not aligned with coordingate in FITS file " + fitsfilename);
      }
      catch(CCfits::HDU::NoSuchKeyword&){
        double ps;
        h0.readKey("PHYSICALSIZE",ps);
        my_res = ps/Nx;
      }
    }
    resolution = fabs(my_res)*pi/180.;
  }
  
  rangeX = resolution*Nx;
  rangeY = resolution*Ny;
  map_boundary_p1[0] = center[0] - (Nx*resolution)/2.;
  map_boundary_p1[1] = center[1] - (Ny*resolution)/2.;
  map_boundary_p2[0] = center[0] + (Nx*resolution)/2.;
  map_boundary_p2[1] = center[1] + (Ny*resolution)/2.;
  
  h0.read(map);
#else
  std::cerr << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
  exit(1);
#endif
}

/** \brief Creates a new PixelMap from a square region of a PixelMap.
 * If the region exceeds the boundaries of the original map, the new map is completed with zeros.
 */
PixelMap::PixelMap(const PixelMap& pmap,  /// Input PixelMap (from which the stamp is taken)
                   const PosType* center, /// center of the region to be duplicated (in rads)
                   std::size_t my_Npixels /// size of the region to be duplicated (in pixels)
)
: map(0.0, my_Npixels*my_Npixels),
Nx(my_Npixels), Ny(my_Npixels), resolution(pmap.resolution)
{
		std::copy(center, center + 2, this->center);
		rangeX = resolution*Nx;
		rangeY = resolution*Ny;
  
		map_boundary_p1[0] = center[0]-(Nx*resolution)/2.;
		map_boundary_p1[1] = center[1]-(Ny*resolution)/2.;
		map_boundary_p2[0] = center[0]+(Nx*resolution)/2.;
		map_boundary_p2[1] = center[1]+(Ny*resolution)/2.;
  
		int edge[2];
		edge[0] = (center[0]-pmap.map_boundary_p1[0])/resolution - Nx/2;
		edge[1] = (center[1]-pmap.map_boundary_p1[1])/resolution - Ny/2;
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

/** \brief Creates a PixelMap at a different resolution.
 * The new counts are calculated integrating over the input pixels.
 * No interpolation or smoothing is performed.
 */
PixelMap::PixelMap(
                   const PixelMap& pmap
                   , PosType res_ratio     /// resolution of map is res_ratio times the resolution of the input map
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
  
  map.resize(Nx*Ny);
  
  int old_Nx = pmap.Nx;
  int old_Ny = pmap.Ny;
  int ix, iy;
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
    
    old_p2[0] = std::min(old_Nx-1,int((ix+1.)*res_ratio));
    old_p2[1] = std::min(old_Ny-1,int((iy+1.)*res_ratio));
    
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

PixelMap::~PixelMap()
{
  map.resize(0);
}

PixelMap& PixelMap::operator=(PixelMap other)
{
  swap(*this, other);
  return *this;
}

/// Zero the whole map
void PixelMap::Clean()
{
  map = 0;
}

/// Multiplies the whole map by a scalar factor
void PixelMap::Renormalize(PosType factor)
{
  map *= factor;
}

/// Adds a value to the i-th pixel
void PixelMap::AddValue(std::size_t i, PosType value)
{
  map[i] += value;
}

/// Assigns a value to the i-th pixel
void PixelMap::AssignValue(std::size_t i, PosType value)
{
  map[i] = value;
}

/// Check whether two PixelMaps agree in their physical dimensions.
bool PixelMap::agrees(const PixelMap& other) const
{
  return
  (Nx == other.Nx) &&
  (Ny == other.Ny) &&
  (resolution == other.resolution) &&
  (center[0] == other.center[0]) &&
  (center[1] == other.center[1]);
}

/// Add the values of another PixelMap to this one.
PixelMap& PixelMap::operator+=(const PixelMap& rhs)
{
  // TODO: maybe check if PixelMaps agree, but is slower
  if(Nx != rhs.getNx() || Ny != rhs.getNy())
    throw std::runtime_error("Dimensions of maps are not compatible");
  for(size_t i=0;i<map.size();++i) map[i] += rhs.map[i];
  return *this;
}

/// Add two PixelMaps.
PixelMap operator+(const PixelMap& a, const PixelMap& b)
{
  PixelMap sum(a);
  sum += b;
  return sum;
}

/// Subtract the values of another PixelMap from this one.
PixelMap& PixelMap::operator-=(const PixelMap& rhs)
{
  // TODO: maybe check if PixelMaps agree, but is slower
  if(Nx != rhs.getNx() || Ny != rhs.getNy())
    throw std::runtime_error("Dimensions of maps are not compatible");
  for(size_t i=0;i<map.size();++i) map[i] -= rhs.map[i];
  return *this;
}

/// Subtract two PixelMaps.
PixelMap operator-(const PixelMap& a, const PixelMap& b)
{
  PixelMap diff(a);
  diff -= b;
  return diff;
}

/// Multiply the values of another PixelMap by this one.
PixelMap& PixelMap::operator*=(const PixelMap& rhs)
{
  if(Nx != rhs.getNx() || Ny != rhs.getNy())
    throw std::runtime_error("Dimensions of maps are not compatible");
  for(size_t i=0;i<map.size();++i) map[i] *= rhs.map[i];
  //map *= rhs.map;
  return *this;
}

PixelMap& PixelMap::operator*=(PosType b)
{
  for(size_t i=0;i<map.size();++i) map[i] *= b;
  //map *= rhs.map;
  return *this;
}

/// Multiply two PixelMaps.
PixelMap operator*(const PixelMap& a, const PixelMap& b)
{
  PixelMap diff(a);
  diff *= b;
  return diff;
}

/// Multiply two PixelMaps.
PixelMap operator*(const PixelMap& a, PosType b)
{
  PixelMap diff(a);
  diff *= b;
  return diff;
}

PosType PixelMap::ave() const{
  PosType tmp=0;
  for(size_t i=0;i<map.size();++i){
    tmp += map[i];
  }
  return tmp/map.size();
}


/** \brief Add an image to the map
 *
 *  If rescale==0 gives constant surface brightness, if < 0
 *  the surface brightness is not scaled by the pixel area as for the flux (default: 1).
 *  Negative values are good for mapping some quantity independant of the pixel size
 */
void PixelMap::AddImages(
                         ImageInfo *imageinfo   /// An array of ImageInfo-s.  There is no reason to separate images for this routine
                         ,int Nimages           /// Number of images on input.
                         ,float rescale         /// rescales the surface brightness while leaving the image unchanged,
///  see full notes
){
  
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

void PixelMap::AddImages(
                         std::vector<ImageInfo> &imageinfo   /// An array of ImageInfo-s.  There is no reason to separate images for this routine
                         ,int Nimages           /// Number of images on input.
                         ,float rescale         /// rescales the surface brightness while leaving the image unchanged,
///  see full notes
){
  AddImages(imageinfo.data(),Nimages,rescale);
}

/// returns the grid points within the branch
void PixelMap::PointsWithinLeaf(Branch * branch1, std::list <unsigned long> &neighborlist){
  
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
bool PixelMap::inMapBox(Branch * branch1) const{
  if (branch1->boundary_p1[0] > map_boundary_p2[0] || branch1->boundary_p2[0] < map_boundary_p1[0]) return false;
  if (branch1->boundary_p1[1] > map_boundary_p2[1] || branch1->boundary_p2[1] < map_boundary_p1[1]) return false;
  return true;
}
/// checks if point is within map boundaries
bool PixelMap::inMapBox(PosType * x) const{
  if (x[0] > map_boundary_p2[0] || x[0] < map_boundary_p1[0]) return false;
  if (x[1] > map_boundary_p2[1] || x[1] < map_boundary_p1[1]) return false;
  return true;
}

//// Finds the area of the intersection between pixel i and branch1
PosType PixelMap::LeafPixelArea(IndexType i,Branch * branch1){
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
void PixelMap::printASCII() const
{
  std::cout << Nx << " " << Ny << "  " << rangeX << std::endl;
  for(std::size_t i=0;i < map.size(); ++i) std::cout << map[i] << std::endl;
  std::cout << Nx << " " << Ny << "  " << rangeX << std::endl;
  
  //map.resize(0);
  return;
}
/// Print an ASCII table of all the pixel values.
void PixelMap::printASCIItoFile(std::string filename) const
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
void PixelMap::printFITS(std::string filename, bool verbose) const
{
#ifdef ENABLE_FITS
  if(filename.empty())
    throw std::invalid_argument("Please enter a valid filename for the FITS file output");
  
  long naxis = 2;
  long naxes[2] = {(long)Nx, (long)Ny};
  
  // might throw CCfits::FITS::CantCreate
  std::auto_ptr<CCfits::FITS> fout(new CCfits::FITS(filename, FLOAT_IMG, naxis, naxes));
  
  std::vector<long> naxex(2);
  naxex[0] = Nx;
  naxex[1] = Ny;
  
  CCfits::PHDU& phout = fout->pHDU();
  
  phout.write(1, map.size(), map);
  
  phout.addKey("WCSAXES", 2, "number of World Coordinate System axes");
  phout.addKey("CRPIX1", 0.5*(naxex[0]+1), "x-coordinate of reference pixel");
  phout.addKey("CRPIX2", 0.5*(naxex[1]+1), "y-coordinate of reference pixel");
  phout.addKey("CRVAL1", 0.0, "first axis value at reference pixel");
  phout.addKey("CRVAL2", 0.0, "second axis value at reference pixel");
  phout.addKey("CTYPE1", "RA---TAN", "the coordinate type for the first axis");
  phout.addKey("CTYPE2", "DEC--TAN", "the coordinate type for the second axis");
  phout.addKey("CUNIT1", "deg     ", "the coordinate unit for the first axis");
  phout.addKey("CUNIT2", "deg     ", "the coordinate unit for the second axis");
  phout.addKey("CDELT1", -180*resolution/pi, "partial of first axis coordinate w.r.t. x");
  phout.addKey("CDELT2", 180*resolution/pi, "partial of second axis coordinate w.r.t. y");
  phout.addKey("CROTA2", 0.0, "");
  phout.addKey("CD1_1", -180*resolution/pi, "partial of first axis coordinate w.r.t. x");
  phout.addKey("CD1_2", 0.0, "partial of first axis coordinate w.r.t. y");
  phout.addKey("CD2_1", 0.0, "partial of second axis coordinate w.r.t. x");
  phout.addKey("CD2_2", 180*resolution/pi, "partial of second axis coordinate w.r.t. y");
  
  phout.addKey("Nx", Nx, "");
  phout.addKey("Ny", Ny, "");
  phout.addKey("range x", map_boundary_p2[0]-map_boundary_p1[0], "radians");
  phout.addKey("RA", center[0], "radians");
  phout.addKey("DEC", center[1], "radians");
  
  if(verbose)
    std::cout << phout << std::endl;
#else
  std::cerr << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
  exit(1);
#endif
}

/** \ingroup Image
 *
 * \brief Smoothes a map with a Gaussian kernel of width sigma (in arcseconds)
 */
void PixelMap::smooth(PosType sigma){
  PosType sum=0,**mask;
  int ix,iy;
  int Nmask,Nmask_half;
  int j_cen, k_cen;
  
  sigma /= 3600.*180/pi;
  Nmask=2*(int)(3*sigma/resolution + 1);
  std::cout << Nmask << std::endl;
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
  
  std::valarray<PosType> map_out(0.0, map.size());
  
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
  
  using std::swap;
  swap(map, map_out);
  
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
void PixelMap::drawline(
                        PosType x1[]     /// one end point of line
                        ,PosType x2[]    /// other end point of line
                        ,PosType value   /// value that it is set to on the map
){
  
  PosType x[2],s1,s2,r;
  size_t index;
  PosType d = 0;
  
  r = sqrt( (x2[0] - x1[0])*(x2[0] - x1[0]) + (x2[1] - x1[1])*(x2[1] - x1[1]) );
  
  if(r==0.0){
    if(inMapBox(x1)){
      //index = Utilities::IndexFromPosition(x1,Nx,range,center);
      index = find_index(x1);
      map[index] = value;
    }
    return;
  }
  
  s1 = (x2[0] - x1[0])/r;
  s2 = (x2[1] - x1[1])/r;
  
  x[0] = x1[0];
  x[1] = x1[1];
  while(d <= r){
    if(inMapBox(x)){
      //index = Utilities::IndexFromPosition(x,Nx,range,center);
      index = find_index(x);
      map[index] = value;
    }
    x[0] += s1*resolution;
    x[1] += s2*resolution;
    d += resolution;
  }
  
  return;
}

/**
 * \brief Draws a circle
 */
void PixelMap::drawcircle(
                          PosType r_center[]    /// center of circle
                          ,PosType radius       ///  radius of circle
                          ,PosType value        /// value that it is set to on the map
){
  
  PosType x1[2],x2[2];
  PosType dtheta = resolution/fabs(radius);
  
  for(float theta = 0; theta < 2*pi; theta += dtheta){
    x1[0] = r_center[0] + radius*cos(theta);
    x1[1] = r_center[1] + radius*sin(theta);
    x2[0] = r_center[0] + radius*cos(theta+dtheta);
    x2[1] = r_center[1] + radius*sin(theta+dtheta);
    drawline(x1,x2,value);
  }
  
  return;
}
void PixelMap::drawgrid(int N,PosType value){
  
  PosType x1[2],x2[2];
  x1[1] = map_boundary_p1[1];
  x2[1] = map_boundary_p2[1];
  for(int i=1;i<N;++i){
    x1[0] = x2[0] = map_boundary_p1[0] + i*rangeX/N;
    drawline(x1,x2,value);
  }
  
  x1[0] = map_boundary_p1[0];
  x2[0] = map_boundary_p2[0];
  for(int i=1;i<N;++i){
    x1[1] = x2[1] = map_boundary_p1[1] + i*rangeY/N;
    drawline(x1,x2,value);
  }
}
void PixelMap::drawBox(PosType p1[],PosType p2[],PosType value){
  PosType x1[2],x2[2];
  
  x1[0] = p1[0];
  x1[1] = p1[1];
  x2[0] = p2[0];
  x2[1] = p1[1];
  drawline(x1,x2,value);
  
  x1[0] = p2[0];
  x1[1] = p1[1];
  x2[0] = p2[0];
  x2[1] = p2[1];
  drawline(x1,x2,value);
  
  x1[0] = p2[0];
  x1[1] = p2[1];
  x2[0] = p1[0];
  x2[1] = p2[1];
  drawline(x1,x2,value);
  
  x1[0] = p1[0];
  x1[1] = p2[1];
  x2[0] = p1[0];
  x2[1] = p1[1];
  drawline(x1,x2,value);
  
  
}

/**
 * \brief Draws a closed curve through the points in curve->imagekist
 *
 * This differs form PixelMap::AddImage() in that it draws lines between the points
 * and takes no account of the cells that the points are in or the surface brightness.
 * The points must be ordered already.  Particularly useful for drawing the caustics
 * that may have irregular cell sizes.  The last point will be connected to the first point.
 */
void PixelMap::AddCurve(ImageInfo *curve,PosType value){
  PosType x[2];
  
  if(curve->imagekist->Nunits() == 0 ) return;
  
  curve->imagekist->MoveToTop();
  x[0] = curve->imagekist->getCurrent()->x[0];
  x[1] = curve->imagekist->getCurrent()->x[1];
  curve->imagekist->Down();
  for(;!(curve->imagekist->OffBottom());curve->imagekist->Down()){
    drawline(x,curve->imagekist->getCurrent()->x,value);
    x[0] = curve->imagekist->getCurrent()->x[0];
    x[1] = curve->imagekist->getCurrent()->x[1];
  }
  curve->imagekist->MoveToTop();
  drawline(x,curve->imagekist->getCurrent()->x,value);
  
  return;
}

void PixelMap::AddCurve(std::vector<Point_2d> &curve,double value){
  PosType x[2];
  
  if(curve.size() == 0 ) return;
  
  x[0] = curve[0][0];
  x[1] = curve[0][1];
  for(size_t ii=1;ii<curve.size();++ii){
    drawline(x,curve[ii].x,value);
    x[0] = curve[ii][0];
    x[1] = curve[ii][1];
  }
  drawline(x,curve[0].x,value);
  
  return;
  
}


/**
 *  \brief Fills in pixels where the image plane points in the grid are located with the value given
 */
void PixelMap::AddGrid(const Grid &grid,PosType value){
  PointList* list = grid.i_tree->pointlist;
  size_t index;
  
  PointList::iterator list_current(list->Top());
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
 *
 *  Warning: When adding a new grid it should not overlap with any of the previously added grids.
 */
void PixelMap::AddGrid(const Grid &grid,LensingVariable val){
  
  if(grid.getNumberOfPoints() == 0 ) return;
  
  AddGrid_(*(grid.i_tree->pointlist),val);
  
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
      list_current = lists[i].Top();
      list_current.JumpDownList( (*treeit)->npoints -1);
      lists[i].setBottom(*list_current);

      ++i;
      allowDecent = false;
    }else{
      allowDecent = true;
    }
    //  }while(grid.i_tree->TreeWalkStep(allowDecent) && i < Nblocks);
  }while(treeit.TreeWalkStep(allowDecent) && i < Nblocks);
  
  std::thread thr[16];
  
  for(int i = 0; i< Nblocks ;++i){
    thr[i] = std::thread(&PixelMap::AddGrid_,this,lists[i],val);
  }
  for(int ii=0;ii<Nblocks;++ii) thr[ii].join();
}

void PixelMap::AddGrid_(const PointList &list,LensingVariable val){
  double tmp,area;
  PosType tmp2[2];
  
  PointList::iterator pl_it(list.Top());
  for(size_t i = 0; i< list.size(); ++i){
    
    switch (val) {
      case ALPHA:
        tmp2[0] = (*pl_it)->x[0] - (*pl_it)->image->x[0];
        tmp2[1] = (*pl_it)->x[1] - (*pl_it)->image->x[1];
        tmp = sqrt(tmp2[0]*tmp2[0] + tmp2[1]*tmp2[1])/resolution/resolution;
        break;
      case ALPHA1:
        tmp = ((*pl_it)->x[0] - (*pl_it)->image->x[0])/resolution/resolution;
        break;
      case ALPHA2:
        tmp = ((*pl_it)->x[1] - (*pl_it)->image->x[1])/resolution/resolution;
        break;
      case KAPPA:
        tmp = (*pl_it)->kappa/resolution/resolution;
        break;
      case GAMMA:
        tmp2[0] = (*pl_it)->gamma[0];
        tmp2[1] = (*pl_it)->gamma[1];
        tmp = sqrt(tmp2[0]*tmp2[0] + tmp2[1]*tmp2[1])/resolution/resolution;
        break;
      case GAMMA1:
        tmp = (*pl_it)->gamma[0]/resolution/resolution;
        break;
      case GAMMA2:
        tmp = (*pl_it)->gamma[1]/resolution/resolution;
        break;
      case GAMMA3:
        tmp = (*pl_it)->gamma[2]/resolution/resolution;
        break;
      case INVMAG:
        tmp = (*pl_it)->invmag/resolution/resolution;
        break;
      case DT:
        tmp = (*pl_it)->dt/resolution/resolution;
        break;
      default:
        std::cerr << "PixelMap::AddGrid() does not work for the input LensingVariable" << std::endl;
        throw std::runtime_error("PixelMap::AddGrid() does not work for the input LensingVariable");
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
    --pl_it;
    
  }
}



/// Find arcs in image  WARNING: THIS IS UNDER CONSTRUCTION!
void PixelMap::FindArc(
                       PosType &radius
                       ,PosType *xc
                       ,PosType *arc_center
                       ,PosType &arclength
                       ,PosType &width
                       ,PosType threshold    // threshold in pixal value
){
  
  if(Nx != Ny){
    std::cout << "PixelMap::FindArc() Doesn't work on nonsquare maps" << std::endl;
    throw std::runtime_error("nonsquare");
  }
  std::vector<size_t> mask(Nx*Nx);
  size_t j=0;
  long k=0;
  PosType const tmp_center[2] = {0,0};
  
  // mask pixels below threshhold
  PosType maxval = map[0],minval = map[0];
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
    std::cout << "PixelMap::FindArc() - No pixels above surface brighness limit" << std::endl;
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

/** \brief Reads all the fits files in a directory into a vector of PixelMaps.
 *
 *  The input fits files must have .fits in their names in addition to the string filespec.
 */
void Utilities::LoadFitsImages(
                               std::string dir              /// path to directory containing fits files
                               ,const std::string& filespec /// string of charactors in fits file name that are matched
                               ,std::vector<PixelMap> & images  /// output vector of PixelMaps
                               ,int maxN       /// maximum number of images that will be read in
                               ,double resolution  /// resolution (rad) of fits image if not given in fits file, use default or -1 otherwise
                               ,bool verbose   /// lists files to stdout
){
  
  DIR *dp = opendir( dir.c_str() );
  struct dirent *dirp;
  struct stat filestat;
  std::string filepath,filename;
  size_t count = 0;
  
  if (dp == NULL)
  {
    throw std::runtime_error("error opening directory");
    return;
  }
  
  while ((dirp = readdir( dp )) && count < maxN)
  {
    filepath = dir + "/" + dirp->d_name;
    
    // If the file is a directory (or is in some way invalid) we'll skip it
    if (stat( filepath.c_str(), &filestat )) continue;
    if (S_ISDIR( filestat.st_mode ))         continue;
    
    filename = dirp->d_name;
    if(filename.find(".fits") !=  std::string::npos){
      if(filename.find(filespec) !=  std::string::npos){
        if(verbose) std::cout << "reading " << filepath << std::endl;
        //PixelMap map(filepath,resolution);
        //images.push_back(std::move(map));
        images.push_back(PixelMap(filepath,resolution));
        ++count;
      }
    }
  }
  
  closedir( dp );
  
  std::cout << count << " fits files read." << std::endl;
  return ;
}
/** \brief Reads all the fits files in a directory into a vector of PixelMaps.
 *
 *  The input fits files must have .fits in their names in addition to the string filespec.
 */
void Utilities::LoadFitsImages(
                               std::string dir              /// path to directory containing fits files
                               ,std::vector<std::string> filespecs /// string of charactors in fits file name that are matched
                               ,std::vector<std::string> file_non_specs /// string of charactors in fits file name cannot have
                               ,std::vector<PixelMap> & images  /// output vector of PixelMaps
                               ,std::vector<std::string> & names  /// file names
                               ,int maxN       /// maximum number of images that will be read in
                               ,double resolution  /// resolution (rad) of fits image if not given in fits file, use default or -1 otherwise
                               ,bool verbose   /// lists files to stdout
){
  
  DIR *dp = opendir( dir.c_str() );
  struct dirent *dirp;
  struct stat filestat;
  std::string filepath,filename;
  size_t count = 0;
  
  
  if (dp == NULL)
  {
    throw std::runtime_error("error opening directory");
    return;
  }
  
  while ((dirp = readdir( dp )) && count < maxN)
  {
    filepath = dir + "/" + dirp->d_name;
    
    // If the file is a directory (or is in some way invalid) we'll skip it
    if (stat( filepath.c_str(), &filestat )) continue;
    if (S_ISDIR( filestat.st_mode ))         continue;
    
    filename = dirp->d_name;
    if(filename.find(".fits") !=  std::string::npos){
      bool read =true;
      for(int i=0;i<filespecs.size();++i) if(filename.find(filespecs[i]) ==  std::string::npos) read = false;
      for(int i=0;i<file_non_specs.size();++i) if(filename.find(file_non_specs[i]) !=  std::string::npos) read = false;
      
      if(read){
        if(verbose) std::cout << "reading " << filepath << std::endl;
        //PixelMap map(filepath,resolution);
        //images.push_back(std::move(map));
        images.push_back(PixelMap(filepath,resolution));
        names.push_back(filepath);
        ++count;
      }
    }
  }
  
  closedir( dp );
  
  std::cout << count << " fits files read." << std::endl;
  return ;
}

/** \brief Reads the file names in a directory that contain a specific sub string.
 
 */
void Utilities::ReadFileNames(
                              std::string dir              /// path to directory containing fits files
                              ,const std::string filespec /// string of charactors in file name that are matched. It can be an empty string.
                              ,std::vector<std::string> & filenames  /// output vector of PixelMaps
                              ,bool verbose){
  
  DIR *dp = opendir( dir.c_str() );
  struct dirent *dirp;
  struct stat filestat;
  std::string filepath,filename;
  
  if (dp == NULL)
  {
    std::cerr << "Cannot find directory" << std::endl;
    throw std::runtime_error("error opening directory");
    return;
  }
  
  while ((dirp = readdir( dp )) )
  {
    filepath = dir + "/" + dirp->d_name;
    
    // If the file is a directory (or is in some way invalid) we'll skip it
    if (stat( filepath.c_str(), &filestat )) continue;
    if (S_ISDIR( filestat.st_mode ))         continue;
    
    filename = dirp->d_name;
    if(filename.find(filespec) !=  std::string::npos){
      if(verbose) std::cout << "adding " << filepath << std::endl;
      filenames.push_back(filename);
    }
  }
  
  closedir( dp );
  
  std::cout << filenames.size() << " file names." << std::endl;
  return ;
}

/*// get the index for a position, returns -1 if out of map
 long PixelMap::find_index(PosType const x[],long &ix,long &iy){
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

/// get the index for a position, returns -1 if out of map
long PixelMap::find_index(PosType const x[],long &ix,long &iy){
  
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
long PixelMap::find_index(PosType const x,PosType const y,long &ix,long &iy){
  
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
long PixelMap::find_index(PosType const x[]){
  long ix,iy;
  return find_index(x,ix,iy);
}
/// get the index for a position, returns -1 if out of map
long PixelMap::find_index(PosType const x,PosType const y){
  long ix,iy;
  return find_index(x,y,ix,iy);
}
/// get the index for a position, returns -1 if out of map
void PixelMap::find_position(PosType x[],std::size_t const index){
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
void PixelMap::find_position(PosType x[],std::size_t const ix,std::size_t const iy){
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

PosType PixelMap::linear_interpolate(PosType x[]){
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

MultiGridSmoother::MultiGridSmoother(
                                     double center[]    /// center of region to be gridded
                                     ,std::size_t Nx    /// number of pixels on x-axis in the highest resolution grid
                                     ,std::size_t Ny    /// number of pixels on y-axis in the highest resolution grid
                                     ,double resolution /// highest resolution to be used, usually the final desired resolution
)
{
  throw std::runtime_error("does not conserve mass yet");
  
  if( (Nx & (Nx-1)) != 0){
    ERROR_MESSAGE();
    std::printf("ERROR: MultiGridSmoother, Nx must be a power of 2\n");
    throw std::runtime_error("ERROR: MultiGridSmoother, Nx must be a power of 2\n");
  }
  if( (Ny & (Ny-1)) != 0){
    ERROR_MESSAGE();
    std::printf("ERROR: MultiGridSmoother, Ny must be a power of 2\n");
    throw std::runtime_error("ERROR: MultiGridSmoother, Ny must be a power of 2\n");
  }
  
  maps.push_back(PixelMap(center,Nx,Ny,resolution));
  while(Nx > 16 && Ny > 16){
    Nx /= 2;
    Ny /= 2;
    resolution *= 2;
    maps.push_back(PixelMap(center,Nx,Ny,resolution));
  }
}
MultiGridSmoother::MultiGridSmoother(double center[],std::size_t Nx,double resolution)
{
  //throw std::runtime_error("does not conserve mass yet");
  
  if( (Nx & (Nx-1)) != 0){
    ERROR_MESSAGE();
    std::printf("ERROR: MultiGridSmoother, Nx must be a power of 2\n");
    throw std::runtime_error("ERROR: MultiGridSmoother, Nx must be a power of 2\n");
  }
  PosType x[2] = {0,0};
  int k=0;
  maps.push_back(PixelMap(center,Nx,resolution));
  interpolators.push_back(Utilities::Interpolator<PixelMap>(x,maps[k].getNx(),maps[k].getRangeX(),maps[k].getNy(),maps[k].getRangeX(),center));
  while(Nx > 16 ){
    Nx /= 2;
    resolution *= 2;
    // resolution = maps[0].getRangeX()/(Nx-1);
    maps.push_back(PixelMap(center,Nx,resolution));
    interpolators.push_back(Utilities::Interpolator<PixelMap>(x,maps[k].getNx(),maps[k].getRangeX(),maps[k].getNy(),maps[k].getRangeX(),center));
  }
}

void MultiGridSmoother::add_particles(std::vector<PosType> x,std::vector<PosType> y){
  long index;
  
  for(int i=0;i<maps.size();++i){
    for(size_t j=0;j<x.size();++j){
      if(x[j] < 0 && y[j] < 0 && i == 8){
        std::cout << "should be in there" << std::endl;
      }
      index = maps[i].find_index(x[j],y[j]);
      assert(index != -1);
      //assert(index != 5);
      if(index != -1) maps[i][index] += 1;
    }
  }
  
  for(int i=0;i<maps[8].size();++i)
    std::cout << "i = " << i << " " << maps[8][i] << std::endl;
}
void MultiGridSmoother::output_map(PixelMap &map,int Nsmooth){
  
  int k;
  double res,x[2];
  long index,ix,iy;
  
  /// find if gids overlap at all
  
  for(size_t i = 0;i < map.size();++i){
    map.find_position(x,i);
    
    k =  maps.size()-1;
    index = maps[k].find_index(x);
    
    if(index != -1){
      
      while(maps[k][index] > Nsmooth && k > 0){
        --k;
        index = maps[k].find_index(x,ix,iy);
      }
      
      res = maps[k].getResolution();
      //PosType c[2] = {maps[k].getCenter()[0],maps[k].getCenter()[1]};
      //Utilities::Interpolator<PixelMap>
      //interp(x,maps[k].getNx(),maps[k].getRangeX(),maps[k].getNy(),maps[k].getRangeX(),c);
      //map[i] = interp.interpolate(x,maps[k])/res/res;
      map[i] = maps[k][index]/res/res;
    }
  }
}

void MultiGridSmoother::_smooth_(int k,size_t i,size_t j,int Nsmooth,PixelMap &map){
  
  if(k==0){  // highest level grid is reached
    map[i+j*map.getNx()] = maps[k][i + j*maps[k].getNx()]/maps[k].getResolution()/maps[k].getResolution();
    return;
  }
  
  if(maps[k][i + j*maps[k].getNx()] > Nsmooth){
    for(int ii = 0;ii < 2;++ii){
      for(int jj = 0;jj < 2;++jj){
        _smooth_(k-1,2*i + ii,2*j + jj,Nsmooth,map);
      }
    }
  }else{
    
    double tmp;
    size_t index;
    PosType x[2];
    // check neighbors
    
    // look through all
    size_t Nratio = map.getNx()/maps[k].getNx();
    //tmp = maps[k][ i + j*maps[k].getNx() ]/maps[k].getResolution()/maps[k].getResolution();
    tmp = maps[k].getResolution()*maps[k].getResolution();
    
    for(int ii = Nratio*i ; ii < Nratio*(i+1) ;++ii){
      for(int jj = Nratio*j ; jj < Nratio*(j+1) ; ++jj){
        // interpolate
        index = ii + jj*map.getNx();
        map.find_position(x, index);
        map[ index ] = maps[k].linear_interpolate(x)/tmp;
        //map[ index ] = tmp;
        
      }
    }
    
  }
  
  return;
}

void MultiGridSmoother::smooth(int Nsmooth,PixelMap &map){
  
  if(!map.agrees(maps[0])){
    throw std::runtime_error("MultiGridSmoother::smooth() map needs to be of the same size as smoother");
  }
  
  int k = maps.size()-1;
  for(size_t i = 0;i < maps[k].getNx();++i){
    for(size_t j = 0;j < maps[k].getNy();++j){
      _smooth_(k,i,j,Nsmooth,map);
    }
  }
  
}






