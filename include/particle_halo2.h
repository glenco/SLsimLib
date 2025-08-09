/**
 * @file particle_halo2.h
 * @brief Defines classes for representing gravitational lensing by collections of simulation particles using tree codes.
 *
 * This header provides the LensHaloParticles and LensHaloParticlesO template classes, which model lensing effects from
 * particle distributions, such as those found in cosmological simulations. The classes support reading particle data from
 * files (ASCII or Gadget2 formats), adaptive smoothing, tree-based force calculations, and geometric transformations.
 *
 * Key Features:
 * - Construction from simulation files or direct particle vectors.
 * - Adaptive smoothing based on local particle density.
 * - Tree-based calculation of lensing quantities (deflection, convergence, shear, potential).
 * - Support for recentering, rotation, and translation of particle distributions.
 * - Optional exclusion of a central region ("hole") for specialized lensing scenarios.
 * - Static factory method for creating halos from simulation snapshots, supporting multiple particle types.
 *
 * Classes:
 * - LensHaloParticles<DType>: Main class for lensing by particle collections.
 * - LensHaloParticlesO<DType>: Extension of LensHaloParticles with support for a central hole.
 *
 * Dependencies:
 * - geometry.h, quadTree.h, simpleTree.h, particle_types.h, utilities_slsim.h, lens_halos.h, oTreeNB.h, gadget.hh
 *
 * Usage:
 * - Instantiate LensHaloParticles or LensHaloParticlesO with appropriate parameters and data sources.
 * - Use provided methods to compute lensing quantities, manipulate particle distributions, and query properties.
 *
 * @author [RB METCALF]
 * @date [Date]
 */
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
#include "gadget.hh"

/**  \brief A class that represents the lensing by a collection of simulation 
  particles by tree code.

   You can create a LensHaloParticles<> directly from a vector of paritcle positions,
   from a file, or use the LensHaloParticles::MakeLensHaloParticle() static function 
   to create a vector of LensHaloParticles from a simulation snapshot, one for each 
   type of particle.  The latter should be used to interface with Gardget2 format 
   files.

   Smoothing is done according to the density of particles in 3D.  All the particles 
   in a single LensHaloParticles object must have the same mass.

   Rays are always shot in the direction of the z-axis, the third coordinate, but 
   the particles can be rotated at the expence of having to rebuild the tree.

   Input format:
     SimFileFormat::ascii - a table of three floats for positions in comoving Mpc (no h factor).
             The lines "# nparticles ...." and "# mass ...." must be in
             header at the top of the file. # is otherwise a comment character.
             Only one type of particle in a single input file.
*/
template<typename DType = float >
class LensHaloParticles : public LensHalo
{
  public:
    LensHaloParticles(
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

      // convert from comoving to physical coordinates
      double scale_factor = 1.0/(1 + redshift);
      for(auto &p : pp){
            p[0] *= scale_factor;
            p[1] *= scale_factor;
            p[2] *= scale_factor;
      }
      setUp(recenter,verbose);
    }

    /// create from a vector of particles
    LensHaloParticles(
      std::vector<Point_3d<DType> > &pvector /// list of particles pdata[][i] should be the position in physical Mpc, the class takes possession of the data and leaves the vector empty
      ,PosType redshift        /// redshift of origin
      ,double my_inv_area      /// inverse area for mass compensation
      ,PosType mass_particle   /// rescale particle masses
      ,const COSMOLOGY& cosmo  /// cosmology
      ,int my_Nsmooth = 5      /// number of neighbours for adaptive smoothing
      ,float Nbucket = 8       /// buckets size in tree
      ,float theta = 0.1       /// opening angle for tree
      ,bool recenter = false   /// re-center on center of mass
      ,bool verbose = false
    ) : LensHalo(redshift,cosmo),particle_mass(mass_particle)
       ,inv_area(my_inv_area),Nsmooth(my_Nsmooth),theta2(theta*theta)
       ,Nbucket(Nbucket)
    {
      std::swap(pp, pvector); // take ownership of the data
      setUp(recenter,verbose);
    }

    ~LensHaloParticles(){};
    
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
    Point_3d<DType> CenterOfMass(){return mcenter;}

    LensHaloParticles(LensHaloParticles &&h):LensHalo(std::move(h)){
      mcenter = h.mcenter;
      pp = std::move(h.pp);
      inv_area = h.inv_area;
      Nsmooth = h.Nsmooth;
      theta2 = h.theta2;
      particle_mass = h.particle_mass;
      Nbucket = h.Nbucket;
      otree = std::move(h.otree);
   }
    LensHaloParticles & operator=(LensHaloParticles &&h){
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

    /// Rotate particles around the point xo by theta_x and then theta_y.
    /// Rays are always shot in the direction of the z-axis, the third coordinate.
    void rotate(
      PosType theta_x
      ,PosType theta_y
      ,Point_3d<DType> &xo);
    
    // reset the opening angle used in the tree force calculation
    void resetForceAngle(double theta){
      theta2 = theta*theta;
    }

    /// Translate particles by the vector xo.
    ///
    /// position (0,0) witll correspond to 
    /// the angular position that can be retreaved
    /// by getTheta() and set with setTheta()
    void translate_particles(Point_3d<DType> &xo);

    /// This is a static function that creats a vector of LensHaloParticles 
    /// from a simulations output file.  Each one needs to be moved into a 
    /// Lens object before it can be used.
    ///
    /// Each type of particle will be put into a separate LensHaloParticles object.
    /// The particles are assumed to be in comoving coordinates.
    /// In the case of Gadget2 files, the masses and reshift are read from the 
    /// file assuming stadard units (comoving kpc/h and 1e10 Msun/h).
    ///  
    /// The redshift and position of the lens can be changed retroactively
    /// by calling setZlens() and setTheta().  The relative PROPER (not comoving) 
    /// positions of the particles are held fixed when changing the redshift.
    ///
    /// The particles can be limited to a sphere of radius rmax centered on xo.
    /// These are comoving Mpc, i.e. simulation units times a.
    static std::vector< LensHaloParticles<DType> > make(
      std::string filename     /// name of data files
      ,SimFileFormat format    /// format of data file
      ,float umass             /// rescale particle masses, this times value in file is in solar masses / h
      ,float ulength           /// rescale length, this times value in file is in pc / h
      ,double inv_area         /// inverse area for mass compensation
      ,const COSMOLOGY& cosmo  /// cosmology
      ,int Nsmooth = 5         /// number of neighbours for adaptive smoothing
      ,float Nbucket = 8       /// buckets size in tree
      ,float theta = 0.1       /// opening angle for tree
      ,DType rmax = 0          /// radius of excised region in comoving Mpc/h, if <=0 no excision is done
      ,Point_3d<DType> xo = Point_3d<DType>(0,0,0) /// center of excised region in comoving Mpc/h
      ,bool verbose = false
    ){
      double r2max= rmax*rmax;
      
        switch (format) {
          case SimFileFormat::gadget2 : 
            {

              std::vector<Point_3d<DType> > tmpdata;
              GadgetFile<Point_3d<DType> > gadget_file(filename,tmpdata);
              gadget_file.openFile();
              gadget_file.readPOS();
              gadget_file.closeFile();

              if(rmax <= 0){
                r2max = 0; // no excision
                xo[0] = gadget_file.boxl/2; // center of box
                xo[1] = gadget_file.boxl/2;
                xo[2] = gadget_file.boxl/2;
              }
              // convert to physical Mpc
              const double a = ulength*1.0e-6/(1 + gadget_file.redshift)/gadget_file.hubble;
              // split particles into different types
              std::vector< std::vector<Point_3d<DType> > > pos(6);
              size_t n=0;
              int ntypes = 0;
              Point_3d<DType> dx;
              for(int i=0; i<6 ; ++i){
                pos[i].resize(gadget_file.npart[i]);
                if(gadget_file.npart[i] > 0) ++ntypes;
                if(r2max > 0){
                  size_t k = 0;
                  for(size_t j=0; j<gadget_file.npart[i]; ++j){
                      pos[i][k] = tmpdata[n++];
                      dx = pos[i][k] - xo;
                      if(dx.length_sqr() < r2max){ 
                        pos[i][k] -= xo;
                        pos[i][k] *= a; // convert to physical Mpc
                        ++k;
                      }
                  }
                  pos[i].resize(k); // resize to number of particles in sphere
                  gadget_file.npart[i] = k; // update number of particles
                }else{
                  for(size_t j=0; j<gadget_file.npart[i]; ++j){
                      pos[i][j] = (tmpdata[n++]-xo)*a;
                  }
                }
              }
              tmpdata.clear();

              std::vector< LensHaloParticles<DType> > halos;
              ntypes = 0;
              for(int i=0; i<6 ; ++i){
                if(gadget_file.npart[i] > 0){
                  halos.push_back( LensHaloParticles<DType>(
                    pos[i]
                    ,gadget_file.redshift
                    ,inv_area
                    ,gadget_file.masstab[i]*umass/gadget_file.hubble
                    ,cosmo
                    ,Nsmooth
                    ,Nbucket
                    ,theta
                    ,false
                    ,verbose)
                  );
                  ++ntypes;
                }
              }
              return halos;
            }
            break;
          default:
            std::cerr << "LensHaloParticles does not accept format of particle data file." << std::endl;
            throw std::invalid_argument("bad format");
        }
      }

      void getBoundingBox(
        Point_3d<DType> &p1,  ///< bottom, left, front corner of box
        Point_3d<DType> &p2   ///< top, right, back corner of box
      ) const {
        otree->getBoundingBox(p1,p2);
      }
protected :
    Point_3d<DType> mcenter;
    //DType *pp;
    std::vector<Point_3d<DType> > pp;
    int Nsmooth;  ///< number of neighbours for adaptive smoothing
    int Nbucket;  ///< number of buckets in tree
    double theta2;  ///< square of opening angle for tree
    double particle_mass;  ///< mass
  
    PosType inv_area;
    std::unique_ptr<OTreeNB<Point_3d<DType> > > otree;

    void readPositionFileASCII(const std::string &filename);
    // construct tree, particles positions must already by stored in comoving Mpc
    void setUp(
        bool recenter           /// center on center of mass
        ,bool verbose
    );
};

// construct tree, particles positions must already by stored in comoving Mpc
template<typename DType>
void LensHaloParticles<DType>::setUp(
        bool recenter           /// center on center of mass
        ,bool verbose
    ){
        long Npoints = pp.size();
  
        Point_3d<double> mcenter(0,0,0);

        for(auto &p : pp){
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

        if(verbose) std::cout << "   Creating particle tree ... " << std::endl;
        otree.reset( new OTreeNB<Point_3d<DType> >(pp.data(),pp.size()) );
        otree->build(Nbucket);
        if(verbose) std::cout << "done." << std::endl;
};

template<typename DType>
void LensHaloParticles<DType>::readPositionFileASCII(const std::string &filename
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

template<typename DType>
void LensHaloParticles<DType>::rotate(PosType theta_x,PosType theta_y
      ,Point_3d<DType> &xo){

  if(theta_x == 0.0 && theta_y == 0.0) return;
  
  PosType coord[3][3];
  PosType cx,cy,sx,sy;
  
  cx = cos(theta_x); sx = sin(theta_x);
  cy = cos(theta_y); sy = sin(theta_y);
  
  coord[0][0] = cy;  coord[1][0] = -sy*sx; coord[2][0] = cx;
  coord[0][1] = 0;   coord[1][1] = cx;     coord[2][1] = sx;
  coord[0][2] = -sy; coord[1][2] = -cy*sx; coord[2][2] = cy*cx;
  
  Point_3d<DType> tmp;
      // rotate particle positions 
  for(auto &p : pp){
    tmp *= 0;
    for(int j=0;j<3;++j){
      double tmp2 = p[j] - xo[j];
      tmp[0] += coord[0][j]*tmp2;
      tmp[1] += coord[1][j]*tmp2;
      tmp[2] += coord[2][j]*tmp2;
    }
    for(int j=0;j<3;++j) p[j]=tmp[j] + xo[j];
  }

      // rotate center of mass
  tmp *= 0.0;
  for(int j=0;j<3;++j){
    tmp[0] += coord[0][j]*(mcenter[j] - xo[j]);
    tmp[1] += coord[1][j]*(mcenter[j] - xo[j]);
    tmp[2] += coord[2][j]*(mcenter[j] - xo[j]);
  }
  mcenter = tmp + xo;

  otree.reset( new OTreeNB<Point_3d<DType> >(pp.data(),pp.size()) );
  otree->build(Nbucket);
};

template<typename DType>
void LensHaloParticles<DType>::translate_particles(Point_3d<DType> &xo){

  if(xo[0] == 0.0 && xo[1] == 0.0) return;
  
      // translate particle positions 
  for(auto &p : pp){
    for(int j=0;j<3;++j) p[j] += xo[j];
  }

  mcenter += xo;

  /// reconstruct the tree
  otree.reset( new OTreeNB<Point_3d<DType> >(pp.data(),pp.size()) );
  otree->build(Nbucket);
};

/** \brief LensHaloParticles but where the "particles" have different masses.
 * 
 */
template<typename DType = float >
class LensHaloParticlesM : public LensHaloParticles<DType>
{
public:

  LensHaloParticlesM(
      std::vector<Point_3d<DType> > &pvector /// list of particles pdata[][i] should be the position in physical Mpc, the class takes possession of the data and leaves the vector empty
      ,PosType redshift        /// redshift of origin
      ,double my_inv_area      /// inverse area for mass compensation
      ,std::vector<DType> Masses   /// rescale particle masses
      ,const COSMOLOGY& cosmo  /// cosmology
      ,int my_Nsmooth = 5      /// number of neighbours for adaptive smoothing
      ,float Nbucket = 8       /// buckets size in tree
      ,float theta = 0.1       /// opening angle for tree
      ,bool recenter = false   /// re-center on center of mass
      ,bool verbose = false
    ) : LensHaloParticles<DType>(pvector,redshift,my_inv_area,1,cosmo,my_Nsmooth
      ,Nbucket,theta,recenter,verbose)
    {
      std::swap(masses,Masses);
      this->otree->calcMoments(masses.data());
    };

  LensHaloParticlesM(LensHaloParticlesM &&h)
        : LensHaloParticles<DType>(std::move(h)),
          masses(std::move(h.masses)) 
  {
  }
  
  LensHaloParticlesM & operator=(LensHaloParticlesM &&h) {
    if (this == &h) return *this;
    LensHaloParticles<DType>::operator=(std::move(h));
    masses = std::move(h.masses);
    return *this;
  }
  /// does not zero lens quantities
  void force_halo(double *alpha
        ,KappaType *kappa
        ,KappaType *gamma
        ,KappaType *phi
        ,double const *xcm 
        ,bool subtract_point=false
        ,PosType screening = 1     // here so that it overrides the LensHalo::force_halo                               
  ){

    *kappa = *phi = 0.0;
    gamma[0] = gamma[1] = 0.0;
    alpha[0] = alpha[1] = 0.0; // ?????
    this->otree->force2D(xcm,masses.data(),this->Nsmooth,this->theta2,this->inv_area
        ,alpha,kappa,gamma,phi);
  
    alpha[0] *= -1;
    alpha[1] *= -1;
  };
 private:
  std::vector<DType> masses;
};


/** \brief LensHaloParticles only with a hole cut out of it.
 * 
 */
template<typename DType = float >
class LensHaloParticlesO : public LensHaloParticles<DType>
{
public:
  LensHaloParticlesO(
    const std::string& simulation_filename /// name of data files
    ,SimFileFormat format    /// format of data file
    ,PosType redshift        /// redshift of origin
    ,double my_inv_area      /// inverse area for mass compensation
    ,PosType mass_particle   /// rescale particle masses
    ,Point_2d x_hole         /// center of hole in comoving Mpc
    ,DType hole_radius       /// radius of hole physical Mpc
    ,const COSMOLOGY& cosmo  /// cosmology
    ,int my_Nsmooth = 5      /// number of neighbours for adaptive smoothing
    ,float Nbucket = 8  /// buckets size in tree
    ,float theta = 0.1      /// opening angle for tree
    ,bool recenter = false /// re-center on center of mass
    ,bool verbose = false
  ) : LensHaloParticles<DType>(simulation_filename,format,redshift,my_inv_area,mass_particle,cosmo,my_Nsmooth,Nbucket,theta,recenter,verbose)
  ,xc(x_hole),rhole(hole_radius)
  {};

  LensHaloParticlesO(
      std::vector<Point_3d<DType> > &pvector /// list of particles pdata[][i] should be the position in physical Mpc, the class takes possession of the data and leaves the vector empty
      ,PosType redshift        /// redshift of origin
      ,double my_inv_area      /// inverse area for mass compensation
      ,PosType mass_particle   /// rescale particle masses
      ,Point_2d x_hole         /// center of hole in comoving Mpc
      ,DType hole_radius       /// radius of hole physical Mpc
      ,const COSMOLOGY& cosmo  /// cosmology
      ,int my_Nsmooth = 5      /// number of neighbours for adaptive smoothing
      ,float Nbucket = 8       /// buckets size in tree
      ,float theta = 0.1       /// opening angle for tree
      ,bool recenter = false   /// re-center on center of mass
      ,bool verbose = false
    ) : LensHaloParticles<DType>(pvector,redshift,my_inv_area,mass_particle,cosmo,my_Nsmooth
      ,Nbucket,theta,recenter,verbose)
    {};

  LensHaloParticlesO(LensHaloParticlesO &&h)
        : LensHaloParticles<DType>(std::move(h)),
          xc(std::move(h.xc)),
          rhole(h.rhole)
  {}
  
  LensHaloParticlesO & operator=(LensHaloParticlesO &&h) {
    if (this == &h) return *this;
    LensHaloParticles<DType>::operator=(std::move(h));
    xc = std::move(h.xc);
    rhole = h.rhole;
    return *this;
  }
  /// does not zero lens quantities
  void force_halo(double *alpha
        ,KappaType *kappa
        ,KappaType *gamma
        ,KappaType *phi
        ,double const *xcm 
        ,bool subtract_point=false
        ,PosType screening = 1     // here so that it overrides the LensHalo::force_halo                               
  ){

    *kappa = *phi = 0.0;
    gamma[0] = gamma[1] = 0.0;
    alpha[0] = alpha[1] = 0.0; // ?????
    this->otree->force2D_hole(xcm,this->particle_mass,this->Nsmooth,this->theta2,this->inv_area
        ,rhole,xc.data()
        ,alpha,kappa,gamma,phi);
  
    alpha[0] *= -1;
    alpha[1] *= -1;
  };

  void setHole(
    Point_2d &x_hole,DType hole_radius
  ){
    xc[0] = x_hole[0];
    xc[1] = x_hole[1];
    rhole = hole_radius;
  };

 private:
  Point_2d xc;
  PosType rhole;
};

/** \brief LensHaloParticlesM only with a hole cut out of it.
 * 
 */
template<typename DType = float >
class LensHaloParticlesMO : public LensHaloParticles<DType>
{
public:

  LensHaloParticlesMO(
      std::vector<Point_3d<DType> > &pvector /// list of particles pdata[][i] should be the position in physical Mpc, the class takes possession of the data and leaves the vector empty
      ,PosType redshift        /// redshift of origin
      ,double my_inv_area      /// inverse area for mass compensation
      ,std::vector<DType> Masses   /// rescale particle masses
      ,Point_2d x_hole         /// center of hole in comoving Mpc
      ,DType hole_radius       /// radius of hole physical Mpc
      ,const COSMOLOGY& cosmo  /// cosmology
      ,int my_Nsmooth = 5      /// number of neighbours for adaptive smoothing
      ,float Nbucket = 8       /// buckets size in tree
      ,float theta = 0.1       /// opening angle for tree
      ,bool recenter = false   /// re-center on center of mass
      ,bool verbose = false
    ) : LensHaloParticles<DType>(pvector,redshift,my_inv_area,1,cosmo,my_Nsmooth
      ,Nbucket,theta,recenter,verbose)
    {
      std::swap(masses,Masses);
      this->otree->calcMoments(masses.data());
    };

  LensHaloParticlesMO(LensHaloParticlesMO &&h)
        : LensHaloParticles<DType>(std::move(h)),
          xc(std::move(h.xc)),
          rhole(h.rhole),
          masses(std::move(h.masses)) 
  {
    
  }
  
  LensHaloParticlesMO & operator=(LensHaloParticlesMO &&h) {
    if (this == &h) return *this;
    LensHaloParticles<DType>::operator=(std::move(h));
    xc = std::move(h.xc);
    rhole = h.rhole;
    masses = std::move(h.masses);
    return *this;
  }
  /// does not zero lens quantities
  void force_halo(double *alpha
        ,KappaType *kappa
        ,KappaType *gamma
        ,KappaType *phi
        ,double const *xcm 
        ,bool subtract_point=false
        ,PosType screening = 1     // here so that it overrides the LensHalo::force_halo                               
  ){

    *kappa = *phi = 0.0;
    gamma[0] = gamma[1] = 0.0;
    alpha[0] = alpha[1] = 0.0; // ?????
    this->otree->force2D_hole(xcm,masses.data(),this->Nsmooth,this->theta2,this->inv_area
        ,rhole,xc.data()
        ,alpha,kappa,gamma,phi);
  
    alpha[0] *= -1;
    alpha[1] *= -1;
  };

  void setHole(
    Point_2d &x_hole,DType hole_radius
  ){
    xc[0] = x_hole[0];
    xc[1] = x_hole[1];
    rhole = hole_radius;
  };

 private:
  Point_2d xc;
  PosType rhole;
  std::vector<DType> masses;
};

#endif
