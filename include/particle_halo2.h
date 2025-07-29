#ifndef GLAMER_particle_halo2_h
#define GLAMER_particle_halo2_h
#include <thread>
#include <mutex>

#include "geometry.h"
#include "quadTree.h"
#include "simpleTree.h"
#include "particle_types.h"
#include "utilities_slsim.h"
#include "lens_halos.h"
#include "oTreeNB.h"

typedef double PType;

class LensHaloParticlesO : public LensHalo
{
  public:
    LensHaloParticlesO(
      const std::string& simulation_filename /// name of data files
      ,SimFileFormat format   /// format of data file
      ,PosType redshift        /// redshift of origin
      ,double my_inv_area           /// inverse area for mass compensation
      ,PosType mass_particle  /// rescale particle masses
      ,const COSMOLOGY& cosmo  /// cosmology
      ,int my_Nsmooth = 5             /// number of neighbours for adaptive smoothing
      ,float Nbucket = 8  /// buckets size in tree
      ,float theta = 0.1      /// opening angle for tree
      ,bool recenter = false /// re-center on center of mass
      ,bool verbose = false
    ) : LensHalo(redshift,cosmo),particle_mass(mass_particle)
       ,inv_area(my_inv_area),Nsmooth(my_Nsmooth),theta2(theta*theta)
       ,Nbucket(Nbucket)
    {

      switch (format) {
        case  SimFileFormat::ascii:
          readPositionFileASCII(simulation_filename);
          break;
        default:
          std::cerr << "LensHaloParticles does not accept format of particle data file." << std::endl;
          throw std::invalid_argument("bad format");
      }
      if(verbose) std::cout << "   Read " << pp.size() << " particles from " << simulation_filename << std::endl;

      setUp(recenter,verbose);
    }

    /// create from a vector of particles
    LensHaloParticlesO(
      std::vector<Point_3d<PType> > &pvector /// list of particles pdata[][i] should be the position in physical Mpc, the class takes possession of the data and leaves the vector empty
      ,PosType redshift        /// redshift of origin
      ,double my_inv_area           /// inverse area for mass compensation
      ,PosType mass_particle  /// rescale particle masses
      ,const COSMOLOGY& cosmo  /// cosmology
      ,int my_Nsmooth = 5             /// number of neighbours for adaptive smoothing
      ,float Nbucket = 8  /// buckets size in tree
      ,float theta = 0.1      /// opening angle for tree
      ,bool recenter = false /// re-center on center of mass
      ,bool verbose = false
    ) : LensHalo(redshift,cosmo),particle_mass(mass_particle)
       ,inv_area(my_inv_area),Nsmooth(my_Nsmooth),theta2(theta*theta)
       ,Nbucket(Nbucket)
    {
      std::swap(pp, pvector); // take ownership of the data
      setUp(recenter,verbose);
    }

    // construct tree, particles positions must already by stored in comoving Mpc
    void setUp(
        bool recenter           /// center on center of mass
        ,bool verbose
    ){
        long Npoints = pp.size();
  
          // convert from comoving to physical coordinates
        PosType scale_factor = 1/(1 + getZlens());
        Point_3d<double> mcenter(0,0,0);

        for(auto &p : pp){
            p[0] *= scale_factor;
            p[1] *= scale_factor;
            p[2] *= scale_factor;

            mcenter[0] += p[0];
            mcenter[1] += p[1];
            mcenter[2] += p[2];
        }
  
        mcenter /= pp.size();
  
        if(verbose) std::cout << "   Center of mass : " << mcenter << std::endl;
  
        if(recenter){
            PosType r2,r2max=0;
            for(auto &p : pp){
                p[0] -= mcenter[0];
                p[1] -= mcenter[1];
                p[2] -= mcenter[2];
            }
            mcenter *= 0;
        }
  
        LensHalo::setTheta(mcenter[0]/Dist,mcenter[1]/Dist);
        LensHalo::mass = particle_mass * pp.size();
        LensHalo::set_flag_elliptical(false); // shouldn't do anything

        otree.reset( new OTreeNB<Point_3d<PType> >(pp.data(),pp.size()) );
        otree->build(Nbucket);
        //otree->calcMoments_point();
    }

    ~LensHaloParticlesO(){};
    
    /// does not zero lens quantities
    void force_halo(double *alpha
        ,KappaType *kappa
        ,KappaType *gamma
        ,KappaType *phi
        ,double const *xcm 
        ,bool subtract_point=false,PosType screening = 1     // here so that it overrides the LensHalo::force_halo                               
    ){

      *kappa = *phi = 0.0;
      gamma[0] = gamma[1] = 0.0;
      alpha[0] = alpha[1] = 0.0; // ?????
      otree->force2D(xcm,particle_mass,Nsmooth,theta2,inv_area
        ,alpha,kappa,gamma,phi);
  
      alpha[0] *= -1;
      alpha[1] *= -1;
    };

    size_t getN() const { return pp.size(); };

    /// get current center of mass in input coordinates
    Point_3d<> CenterOfMass(){return mcenter;}

    LensHaloParticlesO(LensHaloParticlesO &&h):LensHalo(std::move(h)){
      mcenter = h.mcenter;
      pp = std::move(h.pp);
      inv_area = h.inv_area;
      Nsmooth = h.Nsmooth;
      theta2 = h.theta2;
      particle_mass = h.particle_mass;
      Nbucket = h.Nbucket;
      otree = std::move(h.otree);
   }
    LensHaloParticlesO & operator=(LensHaloParticlesO &&h){
    if(this == &h) return *this;
      LensHalo::operator=(std::move(h));
      mcenter = h.mcenter;
      pp = std::move(h.pp);
      inv_area = h.inv_area;
      Nsmooth = h.Nsmooth;
      theta2 = h.theta2;
      particle_mass = h.particle_mass;
      Nbucket = h.Nbucket;
      otree = std::move(h.otree);
      return *this;
    }

    friend class MakeParticleLenses;
protected :
    Point_3d<double> mcenter;
    //PType *pp;
    std::vector<Point_3d<PType> > pp;
    int Nsmooth;  ///< number of neighbours for adaptive smoothing
    int Nbucket;  ///< number of buckets in tree
    double theta2;  ///< square of opening angle for tree
    double particle_mass;  ///< mass
  
    PosType inv_area;
    std::unique_ptr<OTreeNB<Point_3d<double> > > otree;

    void readPositionFileASCII(const std::string &filename);
};

void LensHaloParticlesO::readPositionFileASCII(const std::string &filename
                                                     ){
  
  int ncoll = Utilities::IO::CountColumns(filename);
  if(ncoll < 3 ){
    std::cerr << filename << " should have three columns!" << std::endl;
  }
  if(ncoll > 4 ){
    std::cerr << filename << " Warning! : LensHaloParticles is not using masses." << std::endl;
  }
  
  std::ifstream myfile(filename);
  long Npoints = 0;
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
        {
          if(label == "mass"){
            ss >> tmp_mass;
            ++count;
          }
        }
      }else break;
      if(count == 2 ) break;
    }
    
    if(count == 0){
      std::cerr << "File " << filename << " must have the header lines: " << std::endl
        << "# nparticles   ****" << std::endl;
      throw std::runtime_error("file reading error");
    }
    
    assert(Npoints > 0);
    pp.resize(Npoints);
    size_t row = 0;
    
    // read in particle positions
    do{
        if(str[0] == '#') continue; //for comments
        std::stringstream ss(str);
        
        ss >> pp[row][0];
        if(!(ss >> pp[row][1])) std::cerr << "3 columns are expected in line " << row
          << " of " << filename << std::endl;
        if(!(ss >> pp[row][2])) std::cerr << "3 columns are expected in line " << row
          << " of " << filename << std::endl;
        
        row++;
    }while(std::getline(myfile, str) && row < Npoints);
    
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
};

#endif
