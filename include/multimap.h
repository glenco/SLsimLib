/*
 * MOKAlens.h
 *
 */

#ifndef MultiLENS_H_
#define MultiLENS_H_

#include "standard.h"
#include "profile.h"
#include "InputParams.h"
#include "lens_halos.h"
#include "grid_maintenance.h"
#include "cpfits.h"

#include <stdexcept>

// this is just to make sure no plan is executed with the wrong dimensions
struct my_fftw_plan{
  fftw_plan plan_r2c;
  fftw_plan plan_c2r;
  size_t nx;
  size_t ny;
  
  void initialize(size_t _nx_,size_t _ny_){
    nx = _nx_;
    ny = _ny_;
    
    std::vector<double> rdummy(nx*ny);
    size_t Nkx = (nx/2+1);
    fftw_complex *cdummy = new fftw_complex[ny*Nkx];
    
    plan_r2c = fftw_plan_dft_r2c_2d(ny,nx,rdummy.data(),cdummy
                                          ,FFTW_ESTIMATE);
    plan_c2r = fftw_plan_dft_c2r_2d(ny,nx,cdummy,rdummy.data()
                                          ,FFTW_MEASURE);
    delete[] cdummy;
  }
  
  ~my_fftw_plan(){
    fftw_destroy_plan(plan_r2c);
    fftw_destroy_plan(plan_c2r);
  }
};

/**
 * \brief The MOKA map structure, containing all quantities that define it
 *
 * The MOKA map, that is read in from the fits file. Its components include the
 * lensing properties, as well as the cosmology, the size of the field of view,
 * the redshifts of the lens and source, the properties of the cluster, etc.
 *
 * Note: To use this class requires setting the ENABLE_FITS compiler flag and linking
 * the cfits library.
 */
struct LensMap{
  
  LensMap():nx(0),ny(0),boxlMpc(0),angular_pixel_size(0){};
  
  // move operators
  LensMap(LensMap &&m);
  LensMap& operator=(LensMap &&m);
  
	/// values for the map
	std::valarray<double> surface_density;  // Msun / Mpc^2
	std::valarray<float> alpha1_bar;  // Msun / Mpc
	std::valarray<float> alpha2_bar;  // Msun / Mpc
	std::valarray<float> gamma1_bar;  // Msun / Mpc^2
	std::valarray<float> gamma2_bar;  // Msun / Mpc^2
  std::valarray<float> phi_bar;     // Msun  
  int nx,ny;

  double boxlMpc;
  double angular_pixel_size;  // in radians
	Point_2d center;
  Point_2d lowerleft;  /// boundery with centred grid
  Point_2d upperright; ///
  
  double x_resolution(){return boxlMpc / nx ;}
  double y_resolution(){return (upperright[1]-lowerleft[1])/ny;}
  // # of pixels times resolution
  double x_range(){return boxlMpc;}
  // # of pixels times resolution
  double y_range(){return (upperright[1]-lowerleft[1]);}

  bool evaluate(const double *x,float &sigma,float *gamma,double *alpha);
  

  LensMap(std::string fits_input_file,double angDist){
    read(fits_input_file,angDist);
  }
  
  /// read an entire map
  void read(std::string input_fits,double angDist);//,float h,float z);
  /// read from a file that has been generated with LensMap::write()
  void Myread(std::string fits_input_file);

  
  /// read only header information
  //void read_header(std::string input_fits,float h,float z);
  void read_header(std::string input_fits,double angDist);

  /// read a subsection of the fits map
//  void read_sub(std::string input_fits
//                ,const std::vector<long> &first
//                ,const std::vector<long> &last
//                ,double Dist
//                );
  
    void read_sub(CPFITS_READ &cpfits
                ,std::vector<long> &first
                ,std::vector<long> &last
                ,double Dist
                );

  void write(std::string filename);
  /// meant to output directly in angulare units and lensing quantities
  void write(std::string filename,LensingVariable quant);

  // this calculates the other lensing quantities from the density map
  
  //std::mutex mutex_lensmap;

  // these are thread safe versions
  template <class T>
  void ProcessFFTs(float zerosize
                   ,T Wphi_of_k
                   ,my_fftw_plan &plan_padded // the plan must be initialized with thr right size arrays
                   ,bool do_alpha = true);
  
  template <class T>
  void ProcessFFTs(T Wphi_of_k
                   ,my_fftw_plan &plan   // the plan must be initialized with thr right size arrays
                   ,bool do_alpha = true);

  
  template <class T>
   void ProcessFFTs(float zerosize
                    ,T Wphi_of_k
                     ,bool do_alpha = true);
   
   template <class T>
   void ProcessFFTs(T Wphi_of_k
                     ,bool do_alpha = true);

  void make_fftw_plans(
                       my_fftw_plan &plans
                       ,double zerosize
                       ){
  
    size_t Nnx=int(zerosize*nx);
    size_t Nny=int(zerosize*ny);

    plans.nx = Nnx;
    plans.ny = Nny;
    
    std::vector<double> rdummy(Nnx*Nny);
    size_t Nkx = (Nnx/2+1);
    fftw_complex *cdummy = new fftw_complex[Nny*Nkx];

    plans.plan_r2c = fftw_plan_dft_r2c_2d(Nny,Nnx,rdummy.data(),cdummy
                                    ,FFTW_ESTIMATE);
    plans.plan_c2r = fftw_plan_dft_c2r_2d(Nny,Nnx,cdummy,rdummy.data()
                                    ,FFTW_MEASURE);
    
    delete[] cdummy;
  }
};

/** \brief A lens halo that calculates all lensing quantities on two grids - a low res long range grid
 *   and a high res short range grid.  This is done to reduce the required memory required.
 *
 * Note: To use this class requires setting the ENABLE_FITS compiler flag and linking
 * the cfits library.
 */
class LensHaloMultiMap : public LensHalo
{
public:

  LensHaloMultiMap(
                   std::string fitsfile       /// Original fits map of the density
                   ,std::string dir_data
                   ,double redshift
                   ,double mass_unit      /// should include h factors
                   ,double npane          /// number of submaps in full map
                   ,int number_of_subfields   /// number of subfields
                   ,COSMOLOGY &c
                   ,bool write_subfields = false   /// write subfields to be read if they already exist
                   ,std::string dir_scratch = ""   /// directory for saving long rang force if different than directory where fitsfile is
                   ,bool subtract_ave = true       /// subtract the average of the full field
                   ,double ffactor = 5   // coarse grid size in units of smoothing size
                   ,double gfactor = 5   // ratio of the border size to the smoothing scale
                   );

  ~LensHaloMultiMap(){
    --count;
    if(count == 0){
      plans_set_up = false;
    }
  };
	
  /// these plans will be shared between all instances of LensHaloMultiMap
  
  static my_fftw_plan plan_short_range;
  //static fftw_plan plan_c2r_short_range;
  static my_fftw_plan plan_long_range;
  //static fftw_plan plan_c2r_long_range;
  static std::mutex mutex_multimap;
  static bool plans_set_up;
  static size_t nx_sub;  // size of used highres field
  static size_t ny_sub;
  static size_t nx_long;  // size of used highres field
  static size_t ny_long;
  static size_t nx_sub_extended; // size of highres field with borders
  static size_t ny_sub_extended;
  static int count;

  //const double ffactor = 5,gfactor = 5;
  double ffactor;
  double gfactor;
  //const double ffactor = 10,gfactor = 10;

  // Add highres map be specifying the corners in pixel values
  //void push_back_submap(
  //            const std::vector<long> &lower_left
  //            ,const std::vector<long> &upper_right
  //);
  
  void resetsubmap(int i,
              const std::vector<long> &lower_left
              //,const std::vector<long> &upper_right
              );
  /// Sets the least highres smaller map in physical coordinates relative to the center of the original map, periodic boundary conditions apply
  //void push_back_submapPhys(Point_2d ll,Point_2d ur);
  void resetsubmapPhys(int i,Point_2d ll);//,Point_2d ur);
  
  //void push_back_submapAngular(Point_2d ll,Point_2d ur){
  //  double D = getDist();
  //  push_back_submapPhys(ll*D,ur*D);
  //}
  void resetsubmapAngular(int i,Point_2d ll){ //},Point_2d ur){
    if(i >= short_range_maps.size()) throw std::invalid_argument("");
    double D = getDist();
    resetsubmapPhys(i,ll*D);//,ur*D);
  }

	void force_halo(double *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,double const *xcm,bool subtract_point=false,PosType screening = 1.0);
  
	void writeImage(std::string fn);
  
  /// lower left of long range map in physical Mpc
  Point_2d getLowerLeft_lr() const { return long_range_map.lowerleft; }
  /// upper right of short range map in physical Mpc
  Point_2d getUpperRight_lr() const { return long_range_map.upperright; }
  /// center of short range map in physical Mpc
  Point_2d getCenter_lr() const { return long_range_map.center; }

  /// lower left of short range map in physical Mpc
  Point_2d getLowerLeft_sr(int i) const { return short_range_maps[i].lowerleft; }
  /// upper right of short range map in physical Mpc
  Point_2d getUpperRight_sr(int i) const { return short_range_maps[i].upperright; }
  /// center of short range map in physical Mpc
  Point_2d getCenter_sr(int i) const { return short_range_maps[i].center; }
  
  /// the number of short range grids
  int NumberOfShortRangeMaps(){return short_range_maps.size();}
  
  /// return range of long range map in physical Mpc
  double getRangeMpc_lr() const { return long_range_map.boxlMpc; }
  /// return range of long range map in physical Mpc
  double getRangeMpc_sr() const { return nx_sub*resolution_mpc; }

  /// return number of pixels on a x-axis side in original map
	size_t getNx_lr() const { return long_range_map.nx; }
	/// return number of pixels on a y-axis side in original map
	size_t getNy_lr() const { return long_range_map.ny; }
	
  /// return number of pixels on a x-axis side in short range map
  size_t getNx_sr() const { return nx_sub; }
  /// return number of pixels on a y-axis side in short range map
  size_t getNy_sr() const { return ny_sub; }

  /// return number of pixels on a x-axis side in original map
  size_t getNx() const { return Noriginal[0]; }
  /// return number of pixels on a y-axis side in original map
  size_t getNy() const { return Noriginal[1]; }
  
  double getMax() const {return max_pix;}
  double getMin() const {return min_pix;}
  double getResolutionMpc() const {return resolution_mpc;}
  double getResolutionAngular() const {return angular_resolution;}

  void operator =(LensHaloMultiMap &&m){
    LensHalo::operator=(std::move(m));
    cosmo = m.cosmo;
    long_range_map = std::move(m.long_range_map);
    short_range_maps = std::move(m.short_range_maps);
    //single_grid = m.single_grid;
    cpfits = std::move(m.cpfits);
    //m.cpfits = nullptr;
    max_pix = m.max_pix;
    min_pix = m.min_pix;
    mass_unit = m.mass_unit;
    Noriginal[0] = m.Noriginal[0];
    Noriginal[1] = m.Noriginal[1];
    resolution_mpc = m.resolution_mpc;
    angular_resolution = m.angular_resolution;
    border_width_pix = m.border_width_pix;
    fitsfilename = m.fitsfilename;
    rs2 = m.rs2;
    zerosize = m.zerosize;
    unit = m.unit;
    wsr = m.wsr;
    wlr = m.wlr;
    ave_ang_sd = m.ave_ang_sd;
    subfield_filename = m.subfield_filename;
    write_shorts = m.write_shorts;
    ffactor = m.ffactor;
    gfactor = m.gfactor;
  }
  
  LensHaloMultiMap(LensHaloMultiMap &&m):
  LensHalo(std::move(m)),cosmo(m.cosmo),cpfits(std::move(m.cpfits))
  {
    long_range_map = std::move(m.long_range_map);
    short_range_maps = std::move(m.short_range_maps);
    //single_grid = m.single_grid;
    //cpfits = std::move(m.cpfits);
    //cpfits = m.cpfits;
    //m.cpfits = nullptr;
    max_pix = m.max_pix;
    min_pix = m.min_pix;
    mass_unit = m.mass_unit;
    Noriginal[0] = m.Noriginal[0];
    Noriginal[1] = m.Noriginal[1];
    resolution_mpc = m.resolution_mpc;
    angular_resolution = m.angular_resolution;
    border_width_pix = m.border_width_pix;
    fitsfilename = m.fitsfilename;
    rs2 = m.rs2;
    zerosize = m.zerosize;
    unit = m.unit;
    wsr = m.wsr;
    wlr = m.wlr;
    ave_ang_sd = m.ave_ang_sd;
    subfield_filename = m.subfield_filename;
    write_shorts = m.write_shorts;
    ffactor = m.ffactor;
    gfactor = m.gfactor;
  }

public:
  LensMap long_range_map;
  std::vector<LensMap> short_range_maps;

private:

  bool write_shorts;
  //bool single_grid;
  COSMOLOGY &cosmo;
  CPFITS_READ cpfits;
  
  double ave_ang_sd;
  
  double max_pix = std::numeric_limits<double>::lowest();
  double min_pix = std::numeric_limits<double>::max();
  
  double mass_unit;
  
  size_t Noriginal[2]; // number of pixels in each dimension in original image
  double resolution_mpc;   // resolution of original image and short range image in Mpc
  double angular_resolution;  // angular resolution of original image
  long border_width_pix;   // width of short range maps padding
  std::string fitsfilename;
  std::string subfield_filename;
  
  double rs2;
  
	//const COSMOLOGY& cosmo;
  int zerosize;

  // setsup everything for the given short range map.
  void setsubmap(LensMap &short_range_map
          ,const std::vector<long> &lower_left
          //,const std::vector<long> &upper_right
                                   );
  struct UNIT{
    int operator()(float k2){return 1;}
  };
  struct WLR{
    float rs2;
    float operator()(float k2){return exp(-k2*rs2);}
  };
  struct WSR{
    float rs2;
    float operator()(float k2){return 1 - exp(-k2*rs2);}
  };
  
  UNIT unit;
  WSR wsr;
  WLR wlr;
};


/**
 * \brief pre-process surface mass density map computing deflection angles and shear in FFT,
 *  generalized to work with rectangular maps
 */
template <typename T>
void LensMap::ProcessFFTs(
                          float zerosize
                          ,T Wphi_of_k
                           ,bool do_alpha){
  my_fftw_plan plan_padded;
  
  make_fftw_plans(plan_padded,zerosize);

  ProcessFFTs<T>(zerosize,Wphi_of_k,plan_padded,do_alpha);
}

/// a thread safe version of ProcessFFTs
template <typename T>
void LensMap::ProcessFFTs(
                          float zerosize
                          ,T Wphi_of_k
                          ,my_fftw_plan &plan_padded // the plan must be initialized with thr right size arrays
                          ,bool do_alpha){
  
  assert(surface_density.size() == nx*ny);

  // size of the new map in x and y directions, factor by which each size is increased
  int Nnx=int(zerosize*nx);
  int Nny=int(zerosize*ny);
  size_t NN = Nnx*Nny;
  double boxlx = boxlMpc*zerosize;
  double boxly = ny*boxlMpc*zerosize/nx;
  
  assert( (plan_padded.nx == Nnx)*(plan_padded.ny == Nny) );

  int imin = (Nnx-nx)/2;
  int imax = (Nnx+nx)/2;
  int jmin = (Nny-ny)/2;
  int jmax = (Nny+ny)/2;
  
  size_t Nkx = (Nnx/2+1);
  
  double tmp = 2.*M_PI/boxlx;
  std::vector<double> kxs(Nkx);
  for( int i=0; i<Nkx; i++ ){
    kxs[i] = i*tmp;
  }
  tmp = 2.*M_PI/boxly;
  std::vector<double> kys(Nny);
  for( int j=0; j<Nny; j++ ){
    kys[j]=(j<Nny/2)?double(j):double(j-Nny);
    kys[j] *= tmp;
  }

  std::vector<double> extended_map( NN );
  
  // assume locate in a rectangular map and build up the new one
  for( int j=0; j<Nny; j++ ){
    for( int i=0; i<Nnx; i++ ){
      if(i>=imin && i<imax && j>=jmin && j<jmax){
        int ii = i-imin;
        int jj = j-jmin;
        
        if(ii>=nx || jj>=ny){
          std::cerr << " 1 error mapping " << ii << "  " << jj << std::endl;
          exit(1);
        }
        if(ii<0 || jj<0){
          std::cerr << " 2 error mapping " << ii << "  " << jj << std::endl;
          exit(1);
        }
        assert(ii+nx*jj < surface_density.size());
        assert(i+Nnx*j < extended_map.size());
        extended_map[i+Nnx*j] = surface_density[ii+nx*jj];
        //float tmp = extended_map[i+Nnx*j];
        assert(!isnan(extended_map[i+Nnx*j]));
        if(isinf(extended_map[i+Nnx*j])){
          extended_map[i+Nnx*j] = 0;
        }
    }else{
        extended_map[i+Nnx*j] = 0;
      }
  
    }
  }
  
  //std::vector<fftw_complex> fNmap(Nny*(Nnx/2+1));
  //std::vector<fftw_complex> fphi( Nny*(Nnx/2+1) );
  fftw_complex *fphi = new fftw_complex[Nny*Nkx];

  fftw_execute_dft_r2c(plan_padded.plan_r2c
                       ,extended_map.data(),fphi);

 
  // fourier space
  // std:: cout << " allocating fourier space maps " << std:: endl;
  
  // build modes for each pixel in the fourier space
  for( int i=0; i<Nkx; i++ ){
    for( int j=0; j<Nny; j++ ){
      
      double k2 = kxs[i]*kxs[i] + kys[j]*kys[j];
      size_t k = i+(Nkx)*j;

      // null for k2 = 0 no divergence
      if(k2 == 0){
        fphi[k][0] = 0.;
        fphi[k][1] = 0.;
      }else{
      
        //assert(k < Nny*Nkx);
        //assert(!isnan(fphi[k][0]));
      
        // fphi
        //fphi[i+(Nkx)*j][0]= -2.*fNmap[i+(Nkx)*j][0]/k2;
        //fphi[i+(Nkx)*j][1]= -2.*fNmap[i+(Nkx)*j][1]/k2;

        // fphi
        fphi[k][0] *= -2./k2;
        fphi[k][1] *= -2./k2;

        // apply window function
        double w = Wphi_of_k(k2);
        fphi[k][0] *= w;
        fphi[k][1] *= w;
      }
    }
  }
  
  fftw_complex *fft= new fftw_complex[Nny*(Nkx)];
  //double *realsp = new double[Nnx*Nny];
  
  //std::vector<fftw_complex> fft( Nny*(Nkx) );
  std::vector<double> realsp(Nnx*Nny);
  
  if(do_alpha){
  // alpha1
  {
    // build modes for each pixel in the fourier space
    for( int i=0; i<Nkx; i++ ){
      for( int j=0; j<Nny; j++ ){
        
        size_t k = i+(Nkx)*j;
        fft[k][0] = -kxs[i]*fphi[k][1];
        fft[k][1] =  kxs[i]*fphi[k][0];
        //assert(!isnan(fft[k][0]));
      }
    }
    
    fftw_execute_dft_c2r(plan_padded.plan_c2r,fft,realsp.data());
    
    alpha1_bar.resize(nx*ny,0);
    
    for( int j=jmin; j<jmax; j++ ){
      int jj = j-jmin;
      for( int i=imin; i<imax; i++ ){
        int ii = i-imin;

          alpha1_bar[ii+nx*jj] = -1*float(realsp[i+Nnx*j]/NN);
      }
    }
  }

  // alpha2
  {
    // build modes for each pixel in the fourier space
    for( int j=0; j<Nny; j++ ){
      for( int i=0; i<Nkx; i++ ){
        size_t k = i+(Nkx)*j;
        
        // alpha
        fft[k][0] = -kys[j]*fphi[k][1];
        fft[k][1] =  kys[j]*fphi[k][0];
        
      }
    }
    
    fftw_execute_dft_c2r(plan_padded.plan_c2r,fft,realsp.data());

    alpha2_bar.resize(nx*ny,0);
    
    for( int j=jmin; j<jmax; j++ ){
      int jj = j-jmin;
      for( int i=imin; i<imax; i++ ){
        int ii = i-imin;
        
        alpha2_bar[ii+nx*jj] = -1*float(realsp[i+Nnx*j]/NN);
      }
    }
  }
  }
  // gamma1
  {
    // build modes for each pixel in the fourier space
    for( int i=0; i<Nkx; i++ ){
       for( int j=0; j<Nny; j++ ){
        
        size_t k = i+(Nkx)*j;
        // gamma
         fft[k][0] = 0.5*(kxs[i]*kxs[i]-kys[j]*kys[j])*fphi[k][0];
         fft[k][1] = 0.5*(kxs[i]*kxs[i]-kys[j]*kys[j])*fphi[k][1];
      }
    }
    
    fftw_execute_dft_c2r(plan_padded.plan_c2r,fft,realsp.data());

    gamma1_bar.resize(nx*ny,0);
    
    for( int j=jmin; j<jmax; j++ ){
      int jj = j-jmin;
      for( int i=imin; i<imax; i++ ){
        int ii = i-imin;
        
        gamma1_bar[ii+nx*jj] = float( realsp[i+Nnx*j]/NN);
      }
    }
  }
  // gamma2
  {
    // build modes for each pixel in the fourier space
    for( int i=0; i<Nkx; i++ ){
       for( int j=0; j<Nny; j++ ){
        
        size_t k = i+(Nkx)*j;
        
        // gamma
        fft[k][0] = kxs[i]*kys[j]*fphi[k][0];
        fft[k][1] = kxs[i]*kys[j]*fphi[k][1];
        
      }
    }
    
    fftw_execute_dft_c2r(plan_padded.plan_c2r,fft,realsp.data());

    gamma2_bar.resize(nx*ny);
    
    for( int j=0; j<ny; j++ ){
      int jj = j+jmin;
      for( int i=0; i<nx; i++ ){
        int ii = i+imin;
        
        gamma2_bar[i+nx*j] = float(-realsp[ii+Nnx*jj]/NN);
      }
    }
  }

  // kappa - this is done over because of the window in Fourier space
  {
    
    // build modes for each pixel in the fourier space
    for( int i=0; i<Nkx; i++ ){
      for( int j=0; j<Nny; j++ ){
        
        double k2 = -(kxs[i]*kxs[i] + kys[j]*kys[j])/2;
        
        size_t k = i+(Nkx)*j;
        
        // surface density
        fft[k][0] = k2*fphi[k][0];
        fft[k][1] = k2*fphi[k][1];
        
      }
    }
    
    assert( (plan_padded.nx == Nnx)*(plan_padded.ny == Nny) );
    fftw_execute_dft_c2r(plan_padded.plan_c2r,fft,realsp.data());

    for( int j=jmin; j<jmax; j++ ){
      int jj = j-jmin;
      for( int i=imin; i<imax; i++ ){
        int ii = i-imin;
        
        surface_density[ii+nx*jj] = float(realsp[i+Nnx*j]/NN);
        
      }
    }
  }

  // phi
  {
    
    // build modes for each pixel in the fourier space
    for( int i=0; i<Nkx; i++ ){
      for( int j=0; j<Nny; j++ ){
        size_t k = i+(Nkx)*j;
        
        // surface density
        fft[k][0] = fphi[k][0];
        fft[k][1] = fphi[k][1];
        
      }
    }
    
    assert( (plan_padded.nx == Nnx)*(plan_padded.ny == Nny) );
    fftw_execute_dft_c2r(plan_padded.plan_c2r,fft,realsp.data());
    
    phi_bar.resize(nx*ny);
    for( int j=jmin; j<jmax; j++ ){
      int jj = j-jmin;
      for( int i=imin; i<imax; i++ ){
        int ii = i-imin;
        
        phi_bar[ii+nx*jj] = float(realsp[i+Nnx*j]/NN);
        
      }
    }
  }

  // std:: cout << " remapping the map in the original size " << std:: endl;
  delete[] fft;
  delete[] fphi;
}

/// no padding
template <class T>
 void LensMap::ProcessFFTs(T Wphi_of_k
                  ,bool do_alpha){
   my_fftw_plan plan;
  
   make_fftw_plans(plan,1.0);
   ProcessFFTs<T>(Wphi_of_k,plan,do_alpha);
 }

template <class T>
void LensMap::ProcessFFTs(T Wphi_of_k
                    ,my_fftw_plan &plan   // the plan must be initialized with thr right size arrays
                    ,bool do_alpha){

  
  assert(surface_density.size() == nx*ny);
  assert( (plan.nx == nx)*(plan.ny == ny) );

  // size of the new map in x and y directions, factor by which each size is increased
  
  size_t Nkx = (nx/2+1);
  size_t NN = nx*ny;
  
  double tmp = 2.*M_PI/boxlMpc;
  
  std::vector<double> kxs(Nkx);
  for( int i=0; i<Nkx; i++ ){
    kxs[i] = i*tmp;
  }
  tmp = 2.*M_PI * nx /boxlMpc / ny;
  std::vector<double> kys(ny);
  for( int j=0; j<ny; j++ ){
    kys[j]=(j<ny/2)?double(j):double(j-ny);
    kys[j] *= tmp;
  }
  
  //std::vector<double> extended_map( Nnx*Nny );
  fftw_complex *fphi   = new fftw_complex[ny*Nkx];
  fftw_complex *fft= new fftw_complex[ny*(Nkx)];
  std::vector<double> realsp(NN);
 
  /*
  fftw_plan plan_r2c,plan_c2r;
  {
    std::lock_guard<std::mutex> hold(mu);
    plan_r2c = fftw_plan_dft_r2c_2d(ny,nx,&(surface_density[0])
                                     ,fphi,FFTW_ESTIMATE);
    plan_c2r = fftw_plan_dft_c2r_2d(ny,nx,fft,realsp.data(),FFTW_MEASURE);
  }
*/
  
  fftw_execute_dft_r2c(plan.plan_r2c,&(surface_density[0]),fphi);
  
  // fourier space
  // std:: cout << " allocating fourier space maps " << std:: endl;
  
  // build modes for each pixel in the fourier space
  for( int i=0; i<Nkx; i++ ){
    for( int j=0; j<ny; j++ ){
      
      double k2 = kxs[i]*kxs[i] + kys[j]*kys[j];
      size_t k = i+(Nkx)*j;
      
      // null for k2 = 0 no divergence
      if(k2 == 0){
        fphi[k][0] = 0.;
        fphi[k][1] = 0.;
      }else{
        
        //assert(k < Nny*Nkx);
        //assert(!isnan(fphi[k][0]));
        
        // fphi
        //fphi[i+(Nkx)*j][0]= -2.*fNmap[i+(Nkx)*j][0]/k2;
        //fphi[i+(Nkx)*j][1]= -2.*fNmap[i+(Nkx)*j][1]/k2;
        
        // fphi
        fphi[k][0] *= -2./k2;
        fphi[k][1] *= -2./k2;
        
        // apply window function
        double w = Wphi_of_k(k2);
        fphi[k][0] *= w;
        fphi[k][1] *= w;
      }
    }
  }
  
  if(do_alpha){
  // alpha1
  {
    
    // build modes for each pixel in the fourier space
    for( int i=0; i<Nkx; i++ ){
      for( int j=0; j<ny; j++ ){
        
        size_t k = i + Nkx * j;
        fft[k][0] = -kxs[i]*fphi[k][1];
        fft[k][1] =  kxs[i]*fphi[k][0];
        //assert(!isnan(fft[k][0]));
      }
    }
    
    fftw_execute_dft_c2r(plan.plan_c2r,fft,realsp.data());

    alpha1_bar.resize(NN);
    for( size_t i=0; i<NN; i++ ) alpha1_bar[i] = -1*float(realsp[i]/NN);
  }
  
  // alpha2
  {
    
    // build modes for each pixel in the fourier space
    for( int j=0; j<ny; j++ ){
      for( int i=0; i<Nkx; i++ ){
        size_t k = i+(Nkx)*j;
        
        // alpha
        fft[k][0] = -kys[j]*fphi[k][1];
        fft[k][1] =  kys[j]*fphi[k][0];
        
      }
    }
    
    fftw_execute_dft_c2r(plan.plan_c2r,fft,realsp.data());

    alpha2_bar.resize(NN);
    for( size_t i=0; i<NN; i++ ) alpha2_bar[i] = -1*float(realsp[i]/NN);
  }
    
    // phi 
    {
      
      // build modes for each pixel in the fourier space
      for( int i=0; i<Nkx; i++ ){
        for( int j=0; j<ny; j++ ){
          size_t k = i+(Nkx)*j;
          fft[k][0] = fphi[k][0];
          fft[k][1] = fphi[k][1];
        }
      }
      
      fftw_execute_dft_c2r(plan.plan_c2r,fft,realsp.data());

      phi_bar.resize(NN);
      for( size_t i=0; i<NN; i++ ) phi_bar[i] = float(realsp[i]/NN);
    }
  }
  // gamma1
  {
    
    // build modes for each pixel in the fourier space
    for( int i=0; i<Nkx; i++ ){
      for( int j=0; j<ny; j++ ){
        
        size_t k = i+(Nkx)*j;
        double tmp = 0.5*(kxs[i]*kxs[i] - kys[j]*kys[j]);
        // gamma
        fft[k][0] = tmp * fphi[k][0];
        fft[k][1] = tmp * fphi[k][1];
      }
    }
    
    fftw_execute_dft_c2r(plan.plan_c2r,fft,realsp.data());

    gamma1_bar.resize(NN);
    for(size_t i=0; i<NN; i++ ) gamma1_bar[i] = float( realsp[i]/NN);
  }
  // gamma2
  {
    
    // build modes for each pixel in the fourier space
    for( int i=0; i<Nkx; i++ ){
      for( int j=0; j<ny; j++ ){
        
        size_t k = i+(Nkx)*j;
        double tmp = kxs[i]*kys[j];
        // gamma
        fft[k][0] = tmp * fphi[k][0];
        fft[k][1] = tmp * fphi[k][1];
        
      }
    }
    
    fftw_execute_dft_c2r(plan.plan_c2r,fft,realsp.data());

    gamma2_bar.resize(NN);
    for( size_t i=0; i<NN; i++ ) gamma2_bar[i] = float(-realsp[i]/NN);
  }
  
  // kappa - this is done over because of the window in Fourier space
  {
    
    // build modes for each pixel in the fourier space
    for( int i=0; i<Nkx; i++ ){
      for( int j=0; j<ny; j++ ){
        
        double k2 = -(kxs[i]*kxs[i] + kys[j]*kys[j])/2;
        
        size_t k = i+(Nkx)*j;
        
        // surface density
        fft[k][0] = k2*fphi[k][0];
        fft[k][1] = k2*fphi[k][1];
      }
    }
    
    fftw_execute_dft_c2r(plan.plan_c2r,fft,realsp.data());

    for( size_t i=0; i<NN; i++ ) surface_density[i] = (realsp[i]/NN);
  }
  
  // std:: cout << " remapping the map in the original size " << std:: endl;
  delete[] fft;
  delete[] fphi;
}

#endif
/* MultiLENS_H_ */



