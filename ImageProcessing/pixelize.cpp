/*
 * pixelize.c
 *
 *  Created on: Feb 27, 2010
 *      Author: R.B. Metcalf
 */

//#include "slsimlib.h"
#include <fstream>
#include <algorithm>
#include <utility>
#include <stdexcept>
#include <thread>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "cpfits.h"
#include "image_processing.h"
#include "point.h"
#include "source.h"
#include "gridmap.h"
#include "grid_maintenance.h"

/*#if __cplusplus < 201103L
template<typename T>
void swap(std::valarray<T>& x, std::valarray<T>& y)
{
  std::valarray<T> z(x);
  
  x.resize(y.size());
  x = y;
  
  y.resize(z.size());
  y = z;
}
//#endif
*/

std::string to_string(PixelMapUnits unit){
  switch (unit) {
    case PixelMapUnits::ndef:
      return "not defined";
      break;
    case PixelMapUnits::surfb:
      return "surface brightness (ergs / s / cm**2) ";
      break;
    case PixelMapUnits::count_per_sec:
      return "counts per sec";
      break;
    case PixelMapUnits::mass:
      return "mass";
      break;
    case PixelMapUnits::mass_density:
      return "mass density";
      break;
    case PixelMapUnits::ADU:
      return "ADU";
      break;
    default:
      throw std::invalid_argument("No such unit");
      break;
  }
};


/** \brief Reads all the fits files in a directory into a vector of PixelMaps.
 *
 *  The input fits files must have .fits in their names in addition to the string filespec.
 *
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
}*/
/** \brief Reads all the fits files in a directory into a vector of PixelMaps.
 *
 *  The input fits files must have .fits in their names in addition to the string filespec.
 *
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
}*/

/*** \brief Reads the file names in a directory that contain a specific sub string.
 
 */
void Utilities::ReadFileNames(
                              std::string dir              /// path to directory containing fits files
                              ,const std::string filespec /// string of charactors in file name that are matched. It can be an empty string.
                              ,std::vector<std::string> & filenames  /// output vector of PixelMaps
                              ,const std::string file_non_spec /// string of charactors in file name that file must not have. 
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
    if(filename.find(filespec) !=  std::string::npos &&
       filename.find(file_non_spec) ==  std::string::npos ){
      if(verbose) std::cout << "adding " << filepath << std::endl;
      filenames.push_back(filename);
    }
  }
  
  closedir( dp );
  
  std::cout << filenames.size() << " file names." << std::endl;
  return ;
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
  
  maps.push_back(PixelMap<double>(center,Nx,Ny,resolution));
  while(Nx > 16 && Ny > 16){
    Nx /= 2;
    Ny /= 2;
    resolution *= 2;
    maps.push_back(PixelMap<double>(center,Nx,Ny,resolution));
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
  maps.push_back(PixelMap<double>(center,Nx,resolution));
  interpolators.push_back(Utilities::Interpolator<PixelMap<double> >(x,maps[k].getNx(),maps[k].getRangeX(),maps[k].getNy(),maps[k].getRangeX(),center));
  while(Nx > 16 ){
    Nx /= 2;
    resolution *= 2;
    // resolution = maps[0].getRangeX()/(Nx-1);
    maps.push_back(PixelMap<double>(center,Nx,resolution));
    interpolators.push_back(Utilities::Interpolator<PixelMap<double> >(x,maps[k].getNx(),maps[k].getRangeX(),maps[k].getNy(),maps[k].getRangeX(),center));
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
void MultiGridSmoother::output_map(PixelMap<double> &map,int Nsmooth){
  
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

void MultiGridSmoother::_smooth_(int k,size_t i,size_t j,int Nsmooth,PixelMap<double> &map){
  
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

void MultiGridSmoother::smooth(int Nsmooth,PixelMap<double> &map){
  
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
