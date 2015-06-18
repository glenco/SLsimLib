#include "slsimlib.h"
#include "particle_halo.h"
#include <fstream>

#ifdef ENABLE_FITS
#include <CCfits/CCfits>
using namespace CCfits;
#endif

LensHaloParticles::LensHaloParticles(
            const std::string& simulation_filename
            ,PosType redshift     /// redshift of origin
            ,int Nsmooth         /// number of neighbours for adaptive smoothing
            ,const COSMOLOGY& cosmo /// cosmology
            ,Point_2d theta_rotate /// rotation of particles around the origin
            ,bool recenter
): simfile(simulation_filename)
{
  
  LensHalo::setZlens(redshift);
  LensHalo::setCosmology(cosmo);
  LensHalo::set_Rmax(1.0e3);
  LensHalo::set_flag_elliptical(false);
  
  readPositionFileASCII(simulation_filename);
  sizefile = simfile + "." + std::to_string(Nsmooth) + "sizes";
  if(!readSizesFile(sizefile,Nsmooth)){
    // calculate sizes
    sizes = new float[Npoints];
    calculate_smoothing(Nsmooth);
  }

  LensHalo::set_mass(Npoints*mass);
  
  // convert from comoving to physical coordinates
  PosType scale_factor = 1/(1+redshift);
  for(size_t i=0;i<Npoints;++i){
    xp[i][0] *= scale_factor;
    xp[i][1] *= scale_factor;
    xp[i][2] *= scale_factor;
    
    mcenter[0] += xp[i][0];
    mcenter[1] += xp[i][1];
    mcenter[2] += xp[i][2];
    
  }
  
  mcenter /= Npoints;
  
  if(recenter){
    for(size_t i=0;i<Npoints;++i){
      xp[i][0] -= mcenter[0];
      xp[i][1] -= mcenter[1];
      xp[i][2] -= mcenter[2];
    }
  }
  
  // rotate positions
  rotate_particles(theta_rotate[0],theta_rotate[1]);

  qtree = new TreeQuad(xp,&mass,sizes,Npoints,false,true,0,20);
  
}

LensHaloParticles::~LensHaloParticles(){
  delete qtree;
  delete [] sizes;
  Utilities::free_PosTypeMatrix(xp,Npoints,3);
}

void LensHaloParticles::force_halo(double *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,double const *xcm
                ,bool subtract_point,PosType screening){
  qtree->force2D_recur(xcm,alpha,kappa,gamma,phi);
}

void LensHaloParticles::rotate(Point_2d theta){
  rotate_particles(theta[0],theta[1]);
  delete qtree;
  qtree =new TreeQuad(xp,&mass,sizes,Npoints,false,true,0,20);
}

/** \breaf Reads number of particle and particle positons into Npoint and xp from a ASCII file.
 *
 * Data file must have the lines "# nparticles ***" and "# mass ***" in the header.  All header
 * lines must begin with a "# "
 *
 * Coordinates of particles are in ???? units.
 */
void LensHaloParticles::readPositionFileASCII(const std::string &filename){
  
  std::ifstream myfile(filename);
  
  // find number of particles
  
  if (myfile.is_open()){
    
    std::string str,label;
    int count =0;
    while(std::getline(myfile, str)){
      std::stringstream ss(str);
      ss >> label;
      if(label == "#"){
        ss >> label;
        if(label == "nparticles"){
          ss >> Npoints;
          ++count;
        }
        if(label == "mass"){
          ss >> mass;
          ++count;
        }
      }else break;
      if(count == 2) break;
    }
    
    if(count != 2){
      std::cerr << "File " << filename << " must have the header lines: " << std::endl
      << "# nparticles   ****" << std::endl << "# mass   ****" << std::endl;
      throw std::runtime_error("file reading error");
    }
    
    xp = Utilities::PosTypeMatrix(Npoints,3);
    
    size_t row = 0;
    
    // read in particle positions
    while(std::getline(myfile, str)){
      if(str[0] == '#') continue; //for comments
      std::stringstream ss(str);
      
      ss >> xp[row][0];
      if(!(ss >> xp[row][1])) std::cerr << "3 columns are expected in line " << row
        << " of " << filename << std::endl;
      if(!(ss >> xp[row][2])) std::cerr << "3 columns are expected in line " << row
        << " of " << filename << std::endl;
      
      row++;
    }
    
    if(row != Npoints){
      std::cerr << "Number of data rows in " << filename << " does not match expected number of particles."
      << std::endl;
      throw std::runtime_error("file reading error");
    }
  }else{
    std::cerr << "Unable to open file " << filename << std::endl;
    throw std::runtime_error("file reading error");
  }
  
  std::cout << Npoints << " particle positions read from file " << filename << std::endl;
  
}

bool LensHaloParticles::readSizesFile(const std::string &filename,int Nsmooth){
  
  std::ifstream myfile(filename);
  
  // find number of particles
  
  if (myfile.is_open()){
    
    std::string str,label;
    int count =0;
    size_t Ntmp;
    int NStmp;
    while(std::getline(myfile, str)){
      std::stringstream ss(str);
      ss >> label;
      if(label == "#"){
        ss >> label;
        if(label == "nparticles"){
          ss >> Ntmp;
          if(Ntmp != Npoints){
            std::cerr << "Number of particles in " << filename << " does not match expected number" << std::endl;
            throw std::runtime_error("file reading error");
          }
          ++count;
        }
        if(label == "nsmooth"){
          ss >> NStmp;
          if(NStmp != Nsmooth) return false;
          ++count;
        }
        
       }else break;
      if(count == 2) break;
    }
    
    if(count != 2){
      std::cerr << "File " << filename << " must have the header lines: " << std::endl
      << "# nparticles   ****" << std::endl;
      throw std::runtime_error("file reading error");
    }
    
    sizes = new float[Npoints];
    
    size_t row = 0;
    
    std::cout << "reading in particle sizes from " << filename << "..." << std::endl;
    
    // read in particle sizes
    while(std::getline(myfile, str)){
      if(str[0] == '#') continue; //for comments
      std::stringstream ss(str);
      
      ss >> sizes[row];
      
      row++;
    }
    
    if(row != Npoints){
      std::cerr << "Number of data rows in " << filename << " does not match expected number of particles."
      << std::endl;
      throw std::runtime_error("file reading error");
    }
  }else{
    return false;
  }
  
  std::cout << Npoints << " particle sizes read from file " << filename << std::endl;
  return true;
}
void LensHaloParticles::rotate_particles(PosType theta_x,PosType theta_y){
  
  if(theta_x == 0.0 && theta_y == 0.0) return;
    
  PosType coord[3][3];
  PosType cx,cy,sx,sy;
  
  cx = cos(theta_x); sx = sin(theta_x);
  cy = cos(theta_y); sy = sin(theta_y);
  
  coord[0][0] = cy;  coord[1][0] = -sy*sx; coord[2][0] = cx;
  coord[0][1] = 0;   coord[1][1] = cx;     coord[2][1] = sx;
  coord[0][2] = -sy; coord[1][2] = -cy*sx; coord[2][2] = cy*cx;
  
  PosType tmp[3];
  int j;
  /* rotate particle positions */
  for(size_t i=0;i<Npoints;++i){
    for(j=0;j<3;++j) tmp[j]=0.0;
    for(j=0;j<3;++j){
      tmp[0]+=coord[0][j]*xp[i][j];
      tmp[1]+=coord[1][j]*xp[i][j];
      tmp[2]+=coord[2][j]*xp[i][j];
    }
    for(j=0;j<3;++j) xp[i][j]=tmp[j];
  }
}

void LensHaloParticles::calculate_smoothing(int Nsmooth){
  std::cout << "Calculating smoothing of particles ..." << std::endl
  << "  If there are a lot of particles this could take a while." << std::endl;
  
  // make 3d tree of particle postions
  TreeSimple tree3d(xp,Npoints,10,3,true);
  
  // find distance to nth neighbour for every particle
  IndexType neighbors[Nsmooth];
  for(size_t i=0;i<Npoints;++i){
    tree3d.NearestNeighbors(xp[i],Nsmooth,sizes + i,neighbors);
  }
  std::cout << "done" << std::endl;

  // save result to a file for future use
  writeSizes(sizefile,Nsmooth);
}

void LensHaloParticles::writeSizes(const std::string &filename,int Nsmooth){
  
  std::ofstream myfile(filename);
  
  // find number of particles
  
  if (myfile.is_open()){
    
    std::cout << "Writing particle size information to file " << filename << " ...." << std::endl;
    
    myfile << "# nparticles " << Npoints << std::endl;
    myfile << "# nsmooth " << Nsmooth << std::endl;
    for(size_t i=0;i<Npoints;++i){
      myfile << sizes[i] << std::endl;
      if(!myfile){
        std::cerr << "Unable to write to file " << filename << std::endl;
        throw std::runtime_error("file writing error");
      }
    }
    
    std::cout << "done" << std::endl;
    
  }else{
    std::cerr << "Unable to write to file " << filename << std::endl;
    throw std::runtime_error("file writing error");
  }
}
