#include "slsimlib.h"
#include "particle_halo.h"
#include <fstream>
#include <mutex>
#include <thread>

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
            ,bool my_multimass   /// Set to true is particles have different sizes
            ,PosType MinPSize    
        ):min_size(MinPSize),multimass(my_multimass),simfile(simulation_filename)
{
  
  LensHalo::setZlens(redshift);
  LensHalo::setCosmology(cosmo);
  LensHalo::set_Rsize(1.0e3);
  LensHalo::set_flag_elliptical(false);
  stars_N = 0;
  stars_implanted = false;
  
  Rsize = Rmax = 1.0e3;
  
  readPositionFileASCII(simulation_filename);
  
  sizefile = simfile + "." + std::to_string(Nsmooth) + "sizes";
  if(!readSizesFile(sizefile,Nsmooth,min_size)){
    // calculate sizes
    sizes.resize(Npoints);
    //calculate_smoothing(Nsmooth);

    std::cout << "Calculating smoothing of particles ..." << std::endl
    << Nsmooth << " neighbors.  If there are a lot of particles this could take a while." << std::endl;

    find_smoothing(xp,Npoints,sizes,Nsmooth);
    for(size_t i=0; i<Npoints ; ++i) if(sizes[i] < min_size) sizes[i] = min_size;
    
    // save result to a file for future use
    writeSizes(sizefile,Nsmooth);
  }
  
  // convert from comoving to physical coordinates
  PosType scale_factor = 1/(1+redshift);
  mass = 0.0;
  mcenter *= 0.0;
  PosType max_mass = 0.0,min_mass = HUGE_VALF;
  for(size_t i=0;i<Npoints;++i){
    xp[i][0] *= scale_factor;
    xp[i][1] *= scale_factor;
    xp[i][2] *= scale_factor;
    
    mcenter[0] += xp[i][0]*masses[multimass*i];
    mcenter[1] += xp[i][1]*masses[multimass*i];
    mcenter[2] += xp[i][2]*masses[multimass*i];
    
    mass += masses[multimass*i];

    max_mass = (masses[multimass*i] > max_mass) ? masses[multimass*i] : max_mass;
    min_mass = (masses[multimass*i] < min_mass) ? masses[multimass*i] : min_mass;
  }
  
  mcenter /= mass;
  
  std::cout << "   Particle mass range : " << min_mass << " to " << max_mass << "  ratio of : " << max_mass/min_mass << std::endl;
  

  if(recenter){
    PosType r2,r2max=0;
    for(size_t i=0;i<Npoints;++i){
      xp[i][0] -= mcenter[0];
      xp[i][1] -= mcenter[1];
      xp[i][2] -= mcenter[2];
      
      r2 = xp[i][0]*xp[i][0] + xp[i][1]*xp[i][1] + xp[i][2]*xp[i][2];
      if(r2 > r2max) r2max = r2;
    }
    
    Rsize = sqrt(r2max);
  }
  
  // rotate positions
  rotate_particles(theta_rotate[0],theta_rotate[1]);
  
  qtree = new TreeQuad(xp,masses.data(),sizes.data(),Npoints,multimass,true,0,20);
}

LensHaloParticles::LensHaloParticles(
                                     PosType **positions
                                     ,std::vector<float> &my_sizes
                                     ,std::vector<float> &my_masses
                                     ,PosType redshift
                                     ,const COSMOLOGY& cosmo
                                     ,bool my_multimass
                                     ,PosType sigma_back
                                     ):
xp(positions),min_size(0),multimass(my_multimass)
{
  
  LensHalo::setZlens(redshift);
  LensHalo::setCosmology(cosmo);
  LensHalo::set_Rsize(1.0e3);
  LensHalo::set_flag_elliptical(false);
  stars_N = 0;
  stars_implanted = false;
  
  Rsize = Rmax = 1.0e3;
  
  std::swap(sizes,my_sizes);
  std::swap(masses,my_masses);
  
  // convert from comoving to physical coordinates
  PosType scale_factor = 1/(1+redshift);
  mass = 0.0;
  mcenter *= 0.0;
  
  qtree = new TreeQuad(xp,masses.data(),sizes.data(),Npoints,multimass,true,sigma_back,20);
}

LensHaloParticles::~LensHaloParticles(){
  delete qtree;
  Utilities::free_PosTypeMatrix(xp,Npoints,3);
}

void LensHaloParticles::force_halo(double *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,double const *xcm
                ,bool subtract_point,PosType screening){
  qtree->force2D_recur(xcm,alpha,kappa,gamma,phi);
  
  alpha[0] *= -1;
  alpha[1] *= -1;
}

void LensHaloParticles::rotate(Point_2d theta){
  rotate_particles(theta[0],theta[1]);
  delete qtree;
  qtree =new TreeQuad(xp,masses.data(),sizes.data(),Npoints,multimass,true,0,20);
}

/** \brief Reads number of particle and particle positons into Npoint and xp from a ASCII file.
 *
 * Data file must have the lines "# nparticles ***" and "# mass ***" in the header.  All header
 * lines must begin with a "# "
 *
 * Coordinates of particles are in physical Mpc units.
 */
void LensHaloParticles::readPositionFileASCII(const std::string &filename){
  
  std::ifstream myfile(filename);
  
  // find number of particles
  
  if (myfile.is_open()){
    
    float tmp_mass = 0.0;
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
        if(!multimass){
          if(label == "mass"){
            ss >> tmp_mass;
            ++count;
          }
        }
      }else break;
      if(multimass && count == 1 ) break;
      if(!multimass && count == 2 ) break;
    }
    
    if(count == 0){
      if(multimass) std::cerr << "File " << filename << " must have the header lines: " << std::endl
        << "# nparticles   ****" << std::endl << "# mass   ****" << std::endl;
      if(!multimass) std::cerr << "File " << filename << " must have the header lines: " << std::endl
        << "# nparticles   ****" << std::endl;
      throw std::runtime_error("file reading error");
    }
    
    xp = Utilities::PosTypeMatrix(Npoints,3);
    if(multimass) masses.resize(Npoints);
    else masses.push_back(tmp_mass);
    
    size_t row = 0;
    
    // read in particle positions
    if(!multimass){
      while(std::getline(myfile, str) && row < Npoints){
        if(str[0] == '#') continue; //for comments
        std::stringstream ss(str);
      
        ss >> xp[row][0];
        if(!(ss >> xp[row][1])) std::cerr << "3 columns are expected in line " << row
          << " of " << filename << std::endl;
        if(!(ss >> xp[row][2])) std::cerr << "3 columns are expected in line " << row
          << " of " << filename << std::endl;
      
        row++;
      }
    }else{
      while(std::getline(myfile, str) && row < Npoints){
        if(str[0] == '#') continue; //for comments
        std::stringstream ss(str);
        
        ss >> xp[row][0];
        if(!(ss >> xp[row][1])) std::cerr << "4 columns are expected in line " << row
          << " of " << filename << std::endl;
        if(!(ss >> xp[row][2])) std::cerr << "4 columns are expected in line " << row
          << " of " << filename << std::endl;
        if(!(ss >> masses[row])) std::cerr << "4 columns are expected in line " << row
          << " of " << filename << std::endl;
        
        row++;
      }
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

bool LensHaloParticles::readSizesFile(const std::string &filename,int Nsmooth
                                      ,PosType min_size){
  
  std::ifstream myfile(filename);
  
  // find number of particles
  
  PosType min=HUGE_VALF,max=0.0;
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
    
    sizes.resize(Npoints);
    
    size_t row = 0;
    
    std::cout << "reading in particle sizes from " << filename << "..." << std::endl;
    
    // read in particle sizes
    while(std::getline(myfile, str)){
      if(str[0] == '#') continue; //for comments
      std::stringstream ss(str);
      
      ss >> sizes[row];
      if(min_size > sizes[row] ) sizes[row] = min_size;
      min = min < sizes[row] ? min :  sizes[row];
      max = max > sizes[row] ? max :  sizes[row];
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
  std::cout << "   maximun particle sizes " << max << " minimum " << min << " Mpc" << std::endl;
  
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

void LensHaloParticles::find_smoothing(PosType **xp,size_t N,std::vector<float> &size,int Nneighbors){
  // make 3d tree of particle postions
  TreeSimple tree3d(xp,N,10,3,true);
  
  // find distance to nth neighbour for every particle
  if(N < 1000){
    IndexType neighbors[Nneighbors];
    for(size_t i=0;i<N;++i){
      tree3d.NearestNeighbors(xp[i],Nneighbors,size.data() + i,neighbors);
    }
  }else{
    size_t chunksize = N/N_THREADS;
    std::thread thr[N_THREADS];
    
    size_t NN;
    for(int ii = 0; ii < N_THREADS ;++ii){
      if(ii == N_THREADS-1){
        NN = N - ii*chunksize;
      }else NN = chunksize;
      
      thr[ii] = std::thread(LensHaloParticles::find_smoothing_,&tree3d
                            ,&(xp[ii*chunksize]),&(size[ii*chunksize]),NN,Nneighbors);
    }
    for(int ii = 0; ii < N_THREADS ;++ii) thr[ii].join();
  }

}

void LensHaloParticles::find_smoothing_(TreeSimple *tree3d,PosType **xp,float *sizesp,size_t N,int Nsmooth){
  
  IndexType neighbors[Nsmooth];
  for(size_t i=0;i<N;++i){
    tree3d->NearestNeighbors(xp[i],Nsmooth,sizesp + i,neighbors);
  }
}

/*void LensHaloParticles::calculate_smoothing(int Nsmooth){
  std::cout << "Calculating smoothing of particles ..." << std::endl
  << Nsmooth << " neighbors.  If there are a lot of particles this could take a while." << std::endl;
  
  // make 3d tree of particle postions
  TreeSimple tree3d(xp,Npoints,10,3,true);
  
  // find distance to nth neighbour for every particle
  if(Npoints < 1000){
    IndexType neighbors[Nsmooth];
    for(size_t i=0;i<Npoints;++i){
      tree3d.NearestNeighbors(xp[i],Nsmooth,sizes.data() + i,neighbors);
    }
  }else{
    size_t chunksize = Npoints/N_THREADS;
    std::thread thr[N_THREADS];
    
    size_t N;
    for(int ii = 0; ii < N_THREADS ;++ii){
      if(ii == N_THREADS-1){
        N = Npoints - ii*chunksize;
      }else N = chunksize;
      
      thr[ii] = std::thread(&LensHaloParticles::smooth_,this,&tree3d
                            ,&(xp[ii*chunksize]),&(sizes[ii*chunksize]),N,Nsmooth);
    }
    for(int ii = 0; ii < N_THREADS ;++ii) thr[ii].join();
  }
  std::cout << "done" << std::endl;

  // save result to a file for future use
  writeSizes(sizefile,Nsmooth);
}

void LensHaloParticles::smooth_(TreeSimple *tree3d,PosType **xp,float *sizesp,size_t N,int Nsmooth){

  IndexType neighbors[Nsmooth];
  for(size_t i=0;i<N;++i){
    tree3d->NearestNeighbors(xp[i],Nsmooth,sizesp + i,neighbors);
  }
}*/

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
