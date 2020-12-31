//
//  particle_halo.h
//  GLAMER
//
//  Created by bmetcalf on 16/06/15.
//
//

#ifndef GLAMER_particle_halo_h
#define GLAMER_particle_halo_h

#include "geometry.h"
#include "quadTree.h"
#include "simpleTree.h"
#include "particle_types.h"
#include "utilities_slsim.h"
#include "lens_halos.h"
#include "quadTreeHalos.h"

/**
 *  \brief A class that represents the lensing by a collection of simulation particles.
 
   You can create a LensHaloParticles<> directly from a file, but it is recommended that you use the MakeParticleLenses class to create them and then move them to a Lens object.
 
   Smoothing is done according to the density of particles in 3D.  Smoothing sizes are
   either read in from a file (names simulation_filename + "." + Nsmooth + "sizes") or calculated
   if the file does not exist (in which case the file is created).  This can be
   time and memory consuming when there are a large number of particles.
 
   Input format:
     ASCII - a table of three floats for positions in comoving Mpc (no h factor).
             The lines "# nparticles ...." and "# mass ...." must be in
             header at the top of the file. # is otherwise a comment character.
             Only one type of particle in a single input file.
 
    More input formats will be added in the future.
*/

template<typename PType>
class LensHaloParticles : public LensHalo
{
public:
  
  LensHaloParticles(const std::string& simulation_filename /// name of data files
                    ,SimFileFormat format   /// format of data file
                    ,PosType redshift        /// redshift of origin
                    ,int Nsmooth             /// number of neighbours for adaptive smoothing
                    ,const COSMOLOGY& cosmo  /// cosmology
                    ,Point_2d theta_rotate   /// rotation of particles around the origin
                    ,bool recenter           /// center on center of mass
                    ,bool my_multimass       /// set to true is particles have different sizes
                    ,PosType MinPSize        /// minimum particle size
                    ,PosType rescale_mass = 1.0   /// rescale particle masses
                    ,bool verbose=false
  );
 
  LensHaloParticles(std::vector<PType> &pvector /// list of particles pdata[][i] should be the position in physical Mpc, the class takes possession of the data and leaves the vector empty
                    ,float redshift        /// redshift of origin
                    ,const COSMOLOGY& cosmo  /// cosmology
                    ,Point_2d theta_rotate   /// rotation of particles around the origin
                    ,bool recenter           /// center on center of mass
                    ,float MinPSize        /// minimum particle size
                    ,bool verbose=false
  ):LensHalo(redshift,cosmo), min_size(MinPSize),multimass(true)
  {
    std::swap(pvector,trash_collector);
    pp = trash_collector.data();
    Npoints = trash_collector.size();
    set_up(redshift,cosmo,theta_rotate,recenter,verbose);
  }
  ~LensHaloParticles();
  
  LensHaloParticles(LensHaloParticles &&h):LensHalo(std::move(h)){
    mcenter = h.mcenter;
    densest_point = h.densest_point;
    trash_collector =std::move(h.trash_collector);
    pp = trash_collector.data();
    h.pp = nullptr;
    
    min_size = h.min_size;
    multimass = h.multimass;
    center = h.center;
    Npoints = h.Npoints;
    simfile = h.simfile;
    sizefile = h.sizefile;
    
    qtree = h.qtree;
    h.qtree = nullptr;
  }
  LensHaloParticles<PType> & operator=(LensHaloParticles<PType> &&h);


  void force_halo(double *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,double const *xcm
                  ,bool subtract_point=false,PosType screening = 1.0);

  size_t getN() const { return Npoints; };
  
  /// rotate the simulation around the origin of the simulation coordinates, (radians)
  void rotate(Point_2d theta);
  
  /// get current center of mass in input coordinates

   Point_3d<> CenterOfMass(){return mcenter;}
  /// get the densistt point in input coordinates
   Point_3d<> DensestPoint(){return densest_point;}

  /** \brief This is a test class that makes a truncated SIE out of particles and puts it into a file in the right format for constructing a LensHaloParticles.
   
   This is useful for calculating the level of shot noise and finite source size.  The particles are distributed in 3D according to the SIE profile with only the perpendicular coordinates (1st and 2nd) distorted into an elliptical shape. If the halo is rotated from the original orientation it will not be a oblate spheroid.
   */
  static void makeSIE(
                      std::string new_filename  /// file name to store the particles
                      ,PosType redshift     /// redshift of particles
                      ,double particle_mass /// particle mass
                      ,double total_mass  /// total mass of SIE
                      ,double sigma       /// velocity dispersion in km/s
                      ,double q  /// axis ratio
                      ,Utilities::RandomNumbers_NR &ran
                      );
  
  static void calculate_smoothing(int Nsmooth
                                  ,PType *pp
                                  ,size_t Npoints
                                  ,bool verbose = false);
  
  void readPositionFileASCII(const std::string &filename);
  
  static void writeSizes(const std::string &filename
                         ,int Nsmooth
                         ,const PType *pp
                         ,size_t Npoints
                         );

  static bool readSizesFile(const std::string &filename
                            ,PType * pp
                            ,size_t Npoints
                            ,int Nsmooth
                            ,PosType min_size);

  friend class MakeParticleLenses;
  
protected:
  // constructure for derived classes
  LensHaloParticles(float redshift        /// redshift of origin
                    ,const COSMOLOGY& cosmo  /// cosmology
  ): LensHalo(redshift,cosmo){}

  // This constructor is really only for use by MakeParticleLenses. It does not take
  // possession of the data and so they will not be deleted on destruction
  LensHaloParticles(PType  *pdata          /// particle data (all physical distances)
                    ,size_t Nparticles
                    ,float redshift        /// redshift of origin
                    ,const COSMOLOGY& cosmo  /// cosmology
                    ,Point_2d theta_rotate   /// rotation of particles around the origin
                    ,bool recenter           /// center on center of mass
                    ,float MinPSize        /// minimum particle size
                    ,bool verbose
  ):LensHalo(redshift,cosmo),pp(pdata),min_size(MinPSize),multimass(true),Npoints(Nparticles)
  {
    set_up(redshift,cosmo,theta_rotate,recenter,verbose);
  }

  void rotate_particles(PosType theta_x,PosType theta_y);
  static void smooth_(TreeSimple<PType> *tree3d,PType *xp,size_t N,int Nsmooth);
  void assignParams(InputParams& params);

  Point_3d<> mcenter;
  Point_3d<> densest_point;
  
  PType *pp;
  std::vector<PType> trash_collector;
  
  PosType min_size;
  bool multimass;
  
  Utilities::Geometry::SphericalPoint<> center;
  
  size_t Npoints;
  
  std::string simfile;
  std::string sizefile;
  
  TreeQuadParticles<PType> * qtree;
  void set_up(float redshift,const COSMOLOGY& cosmo,Point_2d theta_rotate,bool recenter,bool verbose);
};

template<typename PType>
LensHaloParticles<PType>::LensHaloParticles(const std::string& simulation_filename
                                            ,SimFileFormat format
                                            ,PosType redshift
                                            ,int Nsmooth
                                            ,const COSMOLOGY& cosmo
                                            ,Point_2d theta_rotate
                                            ,bool recenter
                                            ,bool my_multimass
                                            ,PosType MinPSize
                                            ,PosType massscaling
                                            ,bool verbose
                                            )
:LensHalo(redshift,cosmo),min_size(MinPSize),multimass(my_multimass),simfile(simulation_filename)
{
  
  LensHalo::setZlens(redshift,cosmo);
  LensHalo::setCosmology(cosmo);
  LensHalo::set_flag_elliptical(false);
  
  //stars_N = 0;
  //stars_implanted = false;
  
  Rmax = 1.0e3;
  LensHalo::setRsize(Rmax);
  
  switch (format) {
    case ascii:
      readPositionFileASCII(simulation_filename);
      break;
    default:
      std::cerr << "LensHaloParticles does not accept format of particle data file." << std::endl;
      throw std::invalid_argument("bad format");
  }
  
  sizefile = simfile + "." + std::to_string(Nsmooth) + "sizes";
  
  if(!readSizesFile(sizefile,pp,Npoints,Nsmooth,min_size)){
    
    // calculate sizes
    calculate_smoothing(Nsmooth,pp,Npoints);
    
    // save result to a file for future use
    writeSizes(sizefile,Nsmooth,pp,Npoints);
    //for(size_t i=0; i<Npoints ; ++i) if(sizes[i] < min_size) sizes[i] = min_size;
    for(size_t i=0; i<Npoints ; ++i) if(pp[i].size() < min_size) pp[i].Size = min_size;
  }
  
  double min_s = pp[0].Size;
  size_t smallest_part = 0;
  for(size_t i = 1 ; i<Npoints ; ++i){
    if(pp[i].Size < min_s){
      min_s = pp[i].Size;
      smallest_part = i;
    }
    pp[i].Mass *= massscaling;
  }
  
  densest_point[0] = pp[smallest_part].x[0];
  densest_point[1] = pp[smallest_part].x[1];
  densest_point[2] = pp[smallest_part].x[2];

  
  // convert from comoving to physical coordinates
  PosType scale_factor = 1/(1+redshift);
  mcenter *= 0.0;
  PosType max_mass = 0.0,min_mass = HUGE_VALF,mass=0;
  for(size_t i=0;i<Npoints;++i){
    pp[i][0] *= scale_factor;
    pp[i][1] *= scale_factor;
    pp[i][2] *= scale_factor;
    
    mcenter[0] += pp[i][0]*pp[multimass*i].mass();
    mcenter[1] += pp[i][1]*pp[multimass*i].mass();
    mcenter[2] += pp[i][2]*pp[multimass*i].mass();
    
    mass += pp[multimass*i].mass();
    
    max_mass = (pp[multimass*i].mass() > max_mass) ? pp[multimass*i].mass() : max_mass;
    min_mass = (pp[multimass*i].mass() < min_mass) ? pp[multimass*i].mass() : min_mass;
  }
  LensHalo::setMass(mass);
  
  mcenter /= mass;
  
  if(verbose) std::cout << "   Particle mass range : " << min_mass << " to " << max_mass << "  ratio of : " << max_mass/min_mass << std::endl;
  
  
  if(recenter){
    PosType r2,r2max=0;
    for(size_t i=0;i<Npoints;++i){
      pp[i][0] -= mcenter[0];
      pp[i][1] -= mcenter[1];
      pp[i][2] -= mcenter[2];
      
      r2 = pp[i][0]*pp[i][0] + pp[i][1]*pp[i][1] + pp[i][2]*pp[i][2];
      if(r2 > r2max) r2max = r2;
    }
    
    LensHalo::setRsize( sqrt(r2max) );
    
    densest_point -= mcenter;
    
    mcenter *= 0;
  }
  
  // rotate positions
  rotate_particles(theta_rotate[0],theta_rotate[1]);
  
  qtree = new TreeQuadParticles<PType>(pp,Npoints,-1,-1,0,20);
}

template<typename PType>
void LensHaloParticles<PType>::set_up(
                                 float redshift        /// redshift of origin
                                 ,const COSMOLOGY& cosmo  /// cosmology
                                 ,Point_2d theta_rotate   /// rotation of particles around the origin
                                 ,bool recenter           /// center on center of mass
                                 ,bool verbose
){

  //LensHalo::setZlens(redshift);
  //LensHalo::setCosmology(cosmo);
  LensHalo::set_flag_elliptical(false);
  
  //stars_N = 0;
  //stars_implanted = false;
  
  Rmax = 1.0e3;
  LensHalo::setRsize(Rmax);
  
  double min_s = pp[0].Size;
  size_t smallest_part = 0;
  for(size_t i =0 ; i<Npoints ; ++i){
    if(pp[i].Size < min_s){
      min_s = pp[i].Size;
      smallest_part = i;
    }
  }
  
  densest_point[0] = pp[smallest_part].x[0];
  densest_point[1] = pp[smallest_part].x[1];
  densest_point[2] = pp[smallest_part].x[2];

  // convert from comoving to physical coordinates
  //PosType scale_factor = 1/(1+redshift);
  mcenter *= 0.0;
  PosType max_mass = 0.0,min_mass = HUGE_VALF,mass=0;
  for(size_t i=0;i<Npoints;++i){
  
    mcenter[0] += pp[i][0]*pp[i].mass();
    mcenter[1] += pp[i][1]*pp[i].mass();
    mcenter[2] += pp[i][2]*pp[i].mass();
    
    mass += pp[i].mass();
    
    max_mass = (pp[i].mass() > max_mass) ? pp[i].mass() : max_mass;
    min_mass = (pp[i].mass() < min_mass) ? pp[i].mass() : min_mass;
  }
  LensHalo::setMass(mass);
  
  mcenter /= mass;
  
  if(verbose) std::cout << "   Particle mass range : " << min_mass << " to " << max_mass << "  ratio of : " << max_mass/min_mass << std::endl;
  
  if(recenter){
    PosType r2,r2max=0;
    for(size_t i=0;i<Npoints;++i){
      pp[i][0] -= mcenter[0];
      pp[i][1] -= mcenter[1];
      pp[i][2] -= mcenter[2];
      
      r2 = pp[i][0]*pp[i][0] + pp[i][1]*pp[i][1] + pp[i][2]*pp[i][2];
      if(r2 > r2max) r2max = r2;
      
      densest_point -= mcenter;
      
      mcenter *= 0;
    }
    
    LensHalo::setRsize( sqrt(r2max) );
  }
  
  // rotate positions
  rotate_particles(theta_rotate[0],theta_rotate[1]);
  
  qtree = new TreeQuadParticles<PType>(pp,Npoints,-1,-1,0,20);
}



template<typename PType>
LensHaloParticles<PType>::~LensHaloParticles(){
  delete qtree;
}

template<typename PType>
void LensHaloParticles<PType>::force_halo(double *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi
                                          ,double const *xcm
                                          ,bool subtract_point,PosType screening){
  qtree->force2D_recur(xcm,alpha,kappa,gamma,phi);
  
  alpha[0] *= -1;
  alpha[1] *= -1;
}

/// rotate simulation
template<typename PType>
void LensHaloParticles<PType>::rotate(Point_2d theta){
  rotate_particles(theta[0],theta[1]);
  delete qtree;
  //qtree = new TreeQuadParticles<ParticleType<float> >(pp,Npoints,multimass,true,0,20);
  qtree = new TreeQuadParticles<ParticleType<float> >(pp,Npoints,-1,-1,0,20);
}

/** \brief Reads number of particle and particle positons into Npoint and xp from a ASCII file.
 *
 * Data file must have the lines "# nparticles ***" and "# mass ***" in the header.  All header
 * lines must begin with a "# "
 *
 * Coordinates of particles are in physical Mpc units.
 */
template<typename PType>
void LensHaloParticles<PType>::readPositionFileASCII(const std::string &filename
                                                     ){
  
  int ncoll = Utilities::IO::CountColumns(filename);
  if(!multimass && ncoll != 3 ){
    std::cerr << filename << " should have three columns!" << std::endl;
  }
  if(multimass && ncoll != 4 ){
    std::cerr << filename << " should have four columns!" << std::endl;
  }
  
  std::ifstream myfile(filename);
  
  //size_t Npoints = 0;
  
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
    
    trash_collector.resize(Npoints);  // this is here just to make sure this is deleted
    pp = trash_collector.data();
    //xp = Utilities::PosTypeMatrix(Npoints,3);
    //if(multimass) masses.resize(Npoints);
    //else masses.push_back(tmp_mass);
    
    size_t row = 0;
    
    // read in particle positions
    if(!multimass){
      while(std::getline(myfile, str) && row < Npoints){
        if(str[0] == '#') continue; //for comments
        std::stringstream ss(str);
        
        ss >> pp[row][0];
        if(!(ss >> pp[row][1])) std::cerr << "3 columns are expected in line " << row
          << " of " << filename << std::endl;
        if(!(ss >> pp[row][2])) std::cerr << "3 columns are expected in line " << row
          << " of " << filename << std::endl;
        
        row++;
      }
    }else{
      while(std::getline(myfile, str) && row < Npoints){
        if(str[0] == '#') continue; //for comments
        std::stringstream ss(str);
        
        ss >> pp[row][0];
        if(!(ss >> pp[row][1])) std::cerr << "4 columns are expected in line " << row
          << " of " << filename << std::endl;
        if(!(ss >> pp[row][2])) std::cerr << "4 columns are expected in line " << row
          << " of " << filename << std::endl;
        if(!(ss >> pp[row].Mass)) std::cerr << "4 columns are expected in line " << row
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

template<typename PType>
bool LensHaloParticles<PType>::readSizesFile(const std::string &filename
                                             ,PType * pp
                                             ,size_t Npoints
                                             ,int Nsmooth
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
    
    size_t row = 0;
    
    std::cout << "reading in particle sizes from " << filename << "..." << std::endl;
    
    // read in particle sizes
    while(std::getline(myfile, str)){
      if(str[0] == '#') continue; //for comments
      std::stringstream ss(str);
      
      ss >> pp[row].Size;
      if(min_size > pp[row].size() ) pp[row].Size = min_size;
      min = min < pp[row].size() ? min : pp[row].size();
      max = max > pp[row].size() ? max : pp[row].size();
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

template<typename PType>
void LensHaloParticles<PType>::rotate_particles(PosType theta_x,PosType theta_y){
  
  if(theta_x == 0.0 && theta_y == 0.0) return;
  
  PosType coord[3][3];
  PosType cx,cy,sx,sy;
  
  cx = cos(theta_x); sx = sin(theta_x);
  cy = cos(theta_y); sy = sin(theta_y);
  
  coord[0][0] = cy;  coord[1][0] = -sy*sx; coord[2][0] = cx;
  coord[0][1] = 0;   coord[1][1] = cx;     coord[2][1] = sx;
  coord[0][2] = -sy; coord[1][2] = -cy*sx; coord[2][2] = cy*cx;
  
  //PosType tmp[3];
  Point_3d<> tmp;
  int j;
  // rotate particle positions */
  for(size_t i=0;i<Npoints;++i){
    tmp *= 0;
    for(j=0;j<3;++j){
      tmp[0] += coord[0][j]*pp[i][j];
      tmp[1] += coord[1][j]*pp[i][j];
      tmp[2] += coord[2][j]*pp[i][j];
    }
    for(j=0;j<3;++j) pp[i][j]=tmp[j];
  }
  
  // rotate center of mass
  tmp *= 0.0;
  for(j=0;j<3;++j){
    tmp[0] += coord[0][j]*mcenter[j];
    tmp[1] += coord[1][j]*mcenter[j];
    tmp[2] += coord[2][j]*mcenter[j];
  }
  mcenter = tmp;

  // rotate center of mass
  tmp *= 0.0;
  for(j=0;j<3;++j){
    tmp[0] += coord[0][j]*densest_point[j];
    tmp[1] += coord[1][j]*densest_point[j];
    tmp[2] += coord[2][j]*densest_point[j];
  }
  densest_point = tmp;
}

template<typename PType>
void LensHaloParticles<PType>::calculate_smoothing(int Nsmooth,PType *pp
                                                   ,size_t Npoints
                                                   ,bool verbose
                                                   ){
  
  int nthreads = Utilities::GetNThreads();
  
  if(verbose) std::cout << "Calculating smoothing of particles ..." << std::endl
  << Nsmooth << " neighbors.  If there are a lot of particles this could take a while." << std::endl;
  
  time_t to,t;
  time(&to);
  
  // make 3d tree of particle postions
  TreeSimple<PType> tree3d(pp,Npoints,2*Nsmooth,3,true);
  // find distance to nth neighbour for every particle
  if(Npoints < 1000){
    //IndexType neighbors[Nsmooth];
    for(size_t i=0;i<Npoints;++i){
      pp[i].Size = tree3d.NNDistance(&pp[i][0],Nsmooth + 1);
    }
  }else{
    size_t chunksize = Npoints/nthreads;
    std::vector<std::thread> thr(nthreads);
    
    size_t N;
    for(int ii = 0; ii < nthreads ;++ii){
      if(ii == nthreads - 1){
        N = Npoints - ii*chunksize;
      }else N = chunksize;
      
      thr[ii] = std::thread(LensHaloParticles<PType>::smooth_,&tree3d
                            ,&(pp[ii*chunksize]),N,Nsmooth);
    }
    for(int ii = 0; ii < nthreads ;++ii) thr[ii].join();
  }
  time(&t);
  if(verbose) std::cout << "done in " << difftime(t,to) << " secs" << std::endl;
}

template<typename PType>
void LensHaloParticles<PType>::smooth_(TreeSimple<PType> *tree3d,PType *pp,size_t N,int Nsmooth){
  
  //IndexType neighbors[Nsmooth];
  for(size_t i=0;i<N;++i){
    pp[i].Size = tree3d->NNDistance(&(pp[i][0]),Nsmooth + 1);
  }
}


template<typename PType>
void LensHaloParticles<PType>::writeSizes(const std::string &filename,int Nsmooth
                                          ,const PType *pp,size_t Npoints
                                          ){
  
  std::ofstream myfile(filename);
  
  // find number of particles
  
  if (myfile.is_open()){
    
    std::cout << "Writing particle size information to file " << filename << " ...." << std::endl;
    
    myfile << "# nparticles " << Npoints << std::endl;
    myfile << "# nsmooth " << Nsmooth << std::endl;
    for(size_t i=0;i<Npoints;++i){
      myfile << pp[i].Size << std::endl;
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

template<typename PType>
void LensHaloParticles<PType>::makeSIE(
                                       std::string new_filename  /// file name
                                       ,PosType redshift     /// redshift of particles
                                       ,double particle_mass /// particle mass
                                       ,double total_mass  /// total mass of SIE
                                       ,double sigma       /// velocity dispersion in km/s
                                       ,double q  /// axis ratio
                                       ,Utilities::RandomNumbers_NR &ran
                                       ){
  
  size_t Npoints = total_mass/particle_mass;
  PosType Rmax = (1+redshift)*total_mass*Grav*lightspeed*lightspeed/sigma/sigma/2;
  Point_3d<> point;
  double qq = sqrt(q);
  
  std::ofstream datafile;
  datafile.open(new_filename);
  
  datafile << "# nparticles " << Npoints << std::endl;
  datafile << "# mass " << particle_mass << std::endl;
  // create particles
  for(size_t i=0; i< Npoints ;++i){
    point[0] = ran.gauss();
    point[1] = ran.gauss();
    point[2] = ran.gauss();
    
    point *= Rmax*ran()/point.length();
    
    point[0] *= qq;
    point[1] /= qq;
    
    datafile << point[0] << " " << point[1] << " "
    << point[2] << " " << std::endl;
    
  }
  
  datafile.close();
}

template<typename PType>
LensHaloParticles<PType> & LensHaloParticles<PType>::operator=(LensHaloParticles<PType> &&h){
  if(this == &h) return *this;
  LensHalo::operator=(std::move(h));
  mcenter = h.mcenter;
  densest_point = h.densest_point;
  trash_collector =std::move(h.trash_collector);
  pp = h.pp;  // note: this depends on std:move keeping the pointer valid if it is constructed with the public constructor
  h.pp = nullptr;
  
  min_size = h.min_size;
  multimass = h.multimass;
  center = h.center;
  Npoints = h.Npoints;
  simfile = h.simfile;
  sizefile = h.sizefile;
  
  qtree = h.qtree;
  h.qtree = nullptr;
  
  return *this;
}


/** \brief A class for constructing LensHalos from particles in a data file.
 
 <p>
 The particle data is stored in this structure so the LensHaloParticles should not be copied and then this object allowed to be destroyed.  The halos will be destroyed when this structure is destroyed.
 
 A separate LensHaloParticles is made for each type of particle that is present in the gadget file.
 The nearest N neighbour smoothing is done in 3D on construction separately for
 each type of particle.  The smoothing sizes are automatically saved to files and used again
 if the class is constructed again with the same file and smoothing number.
 
 On construction the LensHalos are not constructed.  You nead to run the `CreateHalos()` method
 the them to be created.
 
 The data file formats are :
 
 gadget2 - standard Gadget-2 output.  A different LensHalo is
 made for each type of particle.  The nearest neighbour
 smoothing scales are calculated within particle type.
 The sph density of gas particles are not used.
 
 csv3,csv4,csv5,csv6 - CSV ascii format without header.  The first three columns
 are the positions.  Next columns are used for the other formats being and
 interpreted as (column 4) masses are in Msun/h, (column 5) the paricle smoothing
 size in Mpc/h and (column 6) an integer for type of particle.  There can be more
 columns in the file than are uesed.  In the case of csv6, when there are more then one
 type of halo each type will be in a differeent LensHaloParticles with differnt smoothing.
 
 glmb - This is a binary format internal to GLAMER used to store
 the positions, masses and sizes of the particles.  If
 GLAMER has generated one, it should be all that is
 needed to recreate the LensHaloParticles.
 
 ascii2 - This is the original ascii GLAMER format.
 Three floats for positions in comoving Mpc (no h factor).
 The lines "# nparticles ...." and "# mass ...." must be
 in header at the top of the file. # is otherwise a
 comment character.  Only one type of particle in a single
 input file.
 
 example of use:

 COSMOLOGY cosmo(Planck);
 
 double zs = 2,zl = 0.5;
 double Dl = cosmo.coorDist(zl);
 
 
 std::string filename = "DataFiles/snap_058_centered.txt";
 
 MakeParticleLenses halomaker(filename,csv4,30,false);
 
 Point_3d Xmax,Xmin;

 halomaker.getBoundingBox(Xmin, Xmax);
 
 Point_3d c_mass = halomaker.getCenterOfMass();
 Point_2d center;

 center[0] = c_mass[0];
 center[1] = c_mass[1];

 // cut out a cylinder, could also do a ball
 halomaker.cylindricalCut(center,(Xmax[0]-Xmin[0]/2));
 
 long seed = 88277394;
 Lens lens(&seed,zs);
 
 double range = (Xmax[0]-Xmin[0])*1.05*cosmo.gethubble()/Dl; // angular range of simulation

 center *= cosmo.gethubble()/Dl; // convert to angular coordinates
 
 halomaker.CreateHalos(cosmo,zl);

 for(auto h : halomaker.halos){
 lens.insertMainHalo(h,zl, true);
 }
 
 GridMap gridmap(&lens, 2049,center.x,range);
 
 PixelMap pmap = gridmap.writePixelMapUniform(KAPPA);
 pmap.printFITS("!" + filename + ".kappa.fits");
 
 pmap = gridmap.writePixelMapUniform(ALPHA1);
 pmap.printFITS("!" + filename + ".alpha1.fits");
 
 pmap = gridmap.writePixelMapUniform(ALPHA2);
 pmap.printFITS("!" + filename + ".alpha2.fits");
 
 <\p>
*/
class MakeParticleLenses{
  
public:
  
  /// vector of LensHalos, one for each type of particle type
  std::vector<LensHaloParticles<ParticleType<float> > *> halos;
  
  /// returns number of particles of each type
  std::vector<size_t> getnp(){return nparticles;}
  
  /// returns mass of particles of each type, if 0 they can have different masses
  std::vector<float> getmp(){return masses;}
  
  /// returns original redshift of snapshot, redshifts of the halos can be changed
  double sim_redshift(){return z_original;}
  
  MakeParticleLenses(const std::string &filename  /// path / root name of gadget-2 snapshot
                     ,SimFileFormat format
                     ,int Nsmooth   /// number of nearest neighbors used for smoothing
                     ,bool recenter /// recenter so that the LenHalos are centered on the center of mass
                     ,bool ignore_type_in_smoothing = false /// used only when format == gadget2, nearest neighbour smoothing is done amongst particles by type if set to false
                     );
  
  MakeParticleLenses(const std::string &filename  /// path / name of glmb file
                     ,bool recenter /// recenter so that the LenHalos are centered on the center of mass
                     );
  
  ~MakeParticleLenses(){
    for(auto p : halos) delete p;
  }
  
  /// recenter the particles to 3d point in physical Mpc/h units  If the halos have already been created they will be destroyed.

  void Recenter(Point_3d<> x);

  void CreateHalos(const COSMOLOGY &cosmo,double redshift);
  
  /// remove particles that are beyond radius (Mpc/h) from center
  void radialCut(Point_3d<> center,double radius);
  /// remove particles that are beyond cylindrical radius (Mpc/h) of center
  void cylindricalCut(Point_2d center,double radius);

  /// returns the original center of mass of all the particles
  Point_3d<> getCenterOfMass() const{return cm;}
  
  /// returns the location of the densest particle in (Mpc/h)
  Point_3d<> densest_particle() const;

  /// return the maximum and minimum coordinates of the particles in each dimension in for the original simulation in Mpc/h
  void getBoundingBox(Point_3d<> &Xmin,Point_3d<> &Xmax) const{
    Xmin = bbox_ll;
    Xmax = bbox_ur;
  }
  
  double getZoriginal(){return z_original;}

  std::vector<ParticleType<float> > data;
private:
  const std::string filename;
  int Nsmooth;
  
  Point_3d<> bbox_ll;  // minumum coordinate values of particles
  Point_3d<> bbox_ur;  // maximim coordinate values of particles
  
  double z_original = -1;
  std::vector<size_t> nparticles;
  std::vector<float> masses;
  Point_3d<> cm;

  // write a glamB format file with all required particle data
  static void writeSizesB(const std::string &filename
                   ,std::vector<ParticleType<float> > &pv
                   ,int Nsmooth,std::vector<size_t> numbytype,double redshift){
    
    assert(numbytype.size() == 6);
    size_t ntot = pv.size();
    std::ofstream myfile(filename, std::ios::out | std::ios::binary);
    if(!myfile.write((char*)&Nsmooth,sizeof(int))){
      std::cerr << "Unable to write to file " << filename << std::endl;
      throw std::runtime_error("file writing error");
    }else{
      myfile.write((char*)&ntot,sizeof(size_t));
      myfile.write((char*)(numbytype.data()),sizeof(size_t)*6);
      myfile.write((char*)&redshift,sizeof(double));
      myfile.write((char*)(pv.data()),sizeof(ParticleType<float>)*ntot);
    }
  }
  
  // read a glamB format file
  static bool readSizesB(const std::string &filename
                   ,std::vector<ParticleType<float> > &pv
                  ,int &Nsmooth,std::vector<size_t> &numbytype,double &redshift){
    
    numbytype.resize(6);
    size_t ntot = pv.size();
    std::ifstream myfile(filename, std::ios::in | std::ios::binary);
    if(!myfile.read((char*)&Nsmooth,sizeof(int))){
      return false;
    }else{
      myfile.read((char*)&ntot,sizeof(size_t));
      myfile.read((char*)(numbytype.data()),sizeof(size_t)*6);
      myfile.read((char*)&redshift,sizeof(double));
      pv.resize(ntot);
      myfile.read((char*)(pv.data()),sizeof(ParticleType<float>)*ntot);
    }
    return true;
  }
  
  // Read particle from Gadget-2 format file
  bool readGadget2(bool ignore_type);
  
  // Reads particles from first 4 columns of csv file
  bool readCSV(int columns_used);

#ifdef ENABLE_HDF5
  bool readHDF5();
#endif

};

/** \brief
 
 */

template<typename HType>
class LensHaloHalos : public LensHalo
{
public:
  
  LensHaloHalos(std::vector<HType> &pvector /// list of particles pdata[][i] should be the position in physical Mpc, the class takes possession of the data and leaves the vector empty
                    ,float redshift        /// redshift of origin
                    ,const COSMOLOGY& cosmo  /// cosmology
                    ,bool verbose=false
  ):LensHalo(redshift,cosmo)
  {
    std::swap(pvector,trash_collector);
    Nhalos = trash_collector.size();
    vpp.resize(Nhalos);

    for(size_t ii = 0 ; ii < Nhalos ; ++ii) vpp[ii] = &trash_collector[ii];
    set_up(redshift,cosmo,verbose);
  }
  ~LensHaloHalos();
  
  LensHaloHalos(LensHaloHalos &&h):LensHalo(std::move(h)){
    mcenter = h.mcenter;
    trash_collector =std::move(h.trash_collector);
    vpp = std::move(h.vpp);
    
    center = h.center;
    Nhalos = h.Nhalos;
    
    qtree = h.qtree;
    h.qtree = nullptr;
  }
  LensHaloHalos & operator=(LensHaloHalos &&h);


  void force_halo(double *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi
                  ,double const *xcm
                  ,bool subtract_point=false,PosType screening = 1.0);

  size_t getN() const { return Nhalos; };
  
  /// get current center of mass in input coordinates
  Point_2d CenterOfMass(){return mcenter;}
  
protected:
  Point_2d mcenter;
  
  std::vector<HType *> vpp;
  std::vector<HType> trash_collector;
  
  Utilities::Geometry::SphericalPoint<> center;
  
  size_t Nhalos;
  
  TreeQuadHalos<HType> * qtree;
  
  //TreeQuadParticles<HType> * qtree;
  void set_up(float redshift,const COSMOLOGY& cosmo,bool verbose);
};

template<typename HType>
void LensHaloHalos<HType>::set_up(
                                 float redshift        /// redshift of origin
                                 ,const COSMOLOGY& cosmo  /// cosmology
                                  ,bool verbose
){

  LensHalo::set_flag_elliptical(false);
  
  Rmax = 1.0e3;  // ????
  LensHalo::setRsize(Rmax);

  // convert from comoving to physical coordinates
  //PosType scale_factor = 1/(1+redshift);
  mcenter *= 0.0;
  PosType max_mass = 0.0,min_mass = HUGE_VALF,mass=0;

  for(HType &h : trash_collector){
    
    mcenter[0] += h[0]*h.get_mass();
    mcenter[1] += h[1]*h.get_mass();
    
    mass += h.get_mass();
    
    max_mass = (h.get_mass() > max_mass) ? h.get_mass() : max_mass;
    min_mass = (h.get_mass() < min_mass) ? h.get_mass() : min_mass;
  }
  LensHalo::setMass(mass);
  
  mcenter /= mass;
  
  if(verbose) std::cout << "   Particle mass range : " << min_mass << " to " << max_mass << "  ratio of : " << max_mass/min_mass << std::endl;
  
  qtree = new TreeQuadHalos<HType>(vpp.data(),Nhalos);
  //qtree = new TreeQuadParticles<HType>(pp.data(),Nhalos,-1,-1,0,20);
}

template<typename HType>
LensHaloHalos<HType>::~LensHaloHalos(){
  delete qtree;
}

template<typename HType>
void LensHaloHalos<HType>::force_halo(double *alpha,KappaType *kappa
                               ,KappaType *gamma,KappaType *phi
                                          ,double const *xcm
                                          ,bool subtract_point,PosType screening){
  qtree->force2D_recur(xcm,alpha,kappa,gamma,phi);

//   ?????
//  alpha[0] *= -1;
//  alpha[1] *= -1;
}

template<typename HType>
LensHaloHalos<HType> & LensHaloHalos<HType>::operator=(LensHaloHalos<HType> &&h){
  if(this == &h) return *this;
  LensHalo::operator=(std::move(h));
  mcenter = h.mcenter;

  trash_collector =std::move(h.trash_collector);
  vpp = std::move(h.vpp);  // note: this depends on std:move keeping the pointer valid if it is constructed with the public constructor
  
  center = h.center;
  Nhalos = h.Nhalos;
  
  qtree = h.qtree;
  h.qtree = nullptr;
  
  return *this;
}

#endif
