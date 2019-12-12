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
  
#ifdef ENABLE_FITS

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

#endif

#ifdef ENABLE_FFTW
  // this calculates the other lensing quantities from the density map
  
  template <class T>
  void PreProcessFFTWMap(float zerosize,T Wphi_of_k,bool do_alpha = true);
  template <class T>
  void PreProcessFFTWMap(T Wphi_of_k,bool do_alpha = true);
 #endif
  
};

/** \brief A lens halo that calculates all lensing qunatities on two grids - a low res long range grid
 *   and a high res short range grid.  This is done to reduce the required memory required.
 *
 * Note: To use this class requires setting the ENABLE_FITS compiler flag and linking
 * the cfits library.
 */
class LensHaloMultiMap : public LensHalo
{
public:

  LensHaloMultiMap(
                   std::string fitsfile            /// Original fits map of the density
                   ,std::string dir_data
                   ,double redshift
                   ,double mass_unit               /// should include h factors
                   ,COSMOLOGY &c
                   ,bool write_subfields = false   /// write subfields to be read if they already exist
                   ,std::string dir_scratch = ""   /// directory for saving long rang force if different than directory where fitsfile is
                   ,bool subtract_ave = true       /// subtract the average of the full field
                   ,bool single_grid_mode = false  /// don't do the short & long range decomposition
                    );

  ~LensHaloMultiMap(){
  };
	
  const double ffactor = 5,gfactor = 5;
  //const double ffactor = 10,gfactor = 10;

  /// Set highres map be specifying the corners in pixel values
  void submap(
              const std::vector<long> &lower_left
              ,const std::vector<long> &upper_right
              );
  /// Sets the highres smaller map in physical coordinates relative to the center of the original map, periodic boundary conditions apply
  void submapPhys(Point_2d ll,Point_2d ur);
  
  void submapAngular(Point_2d ll,Point_2d ur){
    double D = getDist();
    submapPhys(ll*D,ur*D);
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
  Point_2d getLowerLeft_sr() const { return short_range_map.lowerleft; }
  /// upper right of short range map in physical Mpc
  Point_2d getUpperRight_sr() const { return short_range_map.upperright; }
  /// center of short range map in physical Mpc
  Point_2d getCenter_sr() const { return short_range_map.center; }
  
  /// return range of long range map in physical Mpc
  double getRangeMpc_lr() const { return long_range_map.boxlMpc; }
  /// return range of long range map in physical Mpc
  double getRangeMpc_sr() const { return short_range_map.boxlMpc; }

  /// return number of pixels on a x-axis side in original map
	size_t getNx_lr() const { return long_range_map.nx; }
	/// return number of pixels on a y-axis side in original map
	size_t getNy_lr() const { return long_range_map.ny; }
	
  /// return number of pixels on a x-axis side in short range map
  size_t getNx_sr() const { return short_range_map.nx; }
  /// return number of pixels on a y-axis side in short range map
  size_t getNy_sr() const { return short_range_map.ny; }

  /// return number of pixels on a x-axis side in original map
  size_t getNx() const { return Noriginal[0]; }
  /// return number of pixels on a y-axis side in original map
  size_t getNy() const { return Noriginal[1]; }
  
  double getMax() const {return max_pix;}
  double getMin() const {return min_pix;}
  double getResolutionMpc() const {return resolution;}
  double getResolutionAngular() const {return angular_resolution;}

  void operator =(LensHaloMultiMap &&m){
    LensHalo::operator=(std::move(m));
    cosmo = m.cosmo;
    long_range_map = std::move(m.long_range_map);
    short_range_map = std::move(m.short_range_map);
    single_grid = m.single_grid;
    cpfits = std::move(m.cpfits);
    //m.cpfits = nullptr;
    max_pix = m.max_pix;
    min_pix = m.min_pix;
    mass_unit = m.mass_unit;
    Noriginal[0] = m.Noriginal[0];
    Noriginal[1] = m.Noriginal[1];
    resolution = m.resolution;
    angular_resolution = m.angular_resolution;
    border_width = m.border_width;
    fitsfilename = m.fitsfilename;
    rs2 = m.rs2;
    zerosize = m.zerosize;
    unit = m.unit;
    wsr = m.wsr;
    wlr = m.wlr;
    ave_ang_sd = m.ave_ang_sd;
    subfield_filename = m.subfield_filename;
    write_shorts = m.write_shorts;
  }
  
  LensHaloMultiMap(LensHaloMultiMap &&m):
  LensHalo(std::move(m)),cosmo(m.cosmo),cpfits(std::move(m.cpfits))
  {
    long_range_map = std::move(m.long_range_map);
    short_range_map = std::move(m.short_range_map);
    single_grid = m.single_grid;
    //cpfits = std::move(m.cpfits);
    //cpfits = m.cpfits;
    //m.cpfits = nullptr;
    max_pix = m.max_pix;
    min_pix = m.min_pix;
    mass_unit = m.mass_unit;
    Noriginal[0] = m.Noriginal[0];
    Noriginal[1] = m.Noriginal[1];
    resolution = m.resolution;
    angular_resolution = m.angular_resolution;
    border_width = m.border_width;
    fitsfilename = m.fitsfilename;
    rs2 = m.rs2;
    zerosize = m.zerosize;
    unit = m.unit;
    wsr = m.wsr;
    wlr = m.wlr;
    ave_ang_sd = m.ave_ang_sd;
    subfield_filename = m.subfield_filename;
    write_shorts = m.write_shorts;
  }

public:
  LensMap long_range_map;
  LensMap short_range_map;

private:
  
  bool write_shorts;
  bool single_grid;
  COSMOLOGY &cosmo;
  CPFITS_READ cpfits;
  
  double ave_ang_sd;
  
  double max_pix = std::numeric_limits<double>::lowest();
  double min_pix = std::numeric_limits<double>::max();
  
  double mass_unit;
  
  size_t Noriginal[2]; // number of pixels in each dimension in original image
  double resolution;   // resolution of original image and short range image in Mpc
  double angular_resolution;  // angular resolution of original image
  long border_width;   // width of short range maps padding
  std::string fitsfilename;
  std::string subfield_filename;
  
  double rs2;
  
	//const COSMOLOGY& cosmo;
  int zerosize;

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
void LensMap::PreProcessFFTWMap(float zerosize,T Wphi_of_k,bool do_alpha){
  
  assert(surface_density.size() == nx*ny);

  // size of the new map in x and y directions, factor by which each size is increased
  int Nnx=int(zerosize*nx);
  int Nny=int(zerosize*ny);
  size_t NN = Nnx*Nny;
  double boxlx = boxlMpc*zerosize;
  double boxly = ny*boxlMpc*zerosize/nx;
  
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
          std::cout << " 1 error mapping " << ii << "  " << jj << std::endl;
          exit(1);
        }
        if(ii<0 || jj<0){
          std::cout << " 2 error mapping " << ii << "  " << jj << std::endl;
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
  fftw_complex *fphi   = new fftw_complex[Nny*Nkx];

  fftw_plan p = fftw_plan_dft_r2c_2d(Nny,Nnx,extended_map.data(),fphi,FFTW_ESTIMATE);
  
  fftw_execute( p );
 
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

   fftw_plan pp = fftw_plan_dft_c2r_2d(Nny,Nnx,fft,realsp.data(),FFTW_MEASURE);
  
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
    
    fftw_execute( pp );
    
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
    
    fftw_execute( pp );
    
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
    
    fftw_execute( pp );
    
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
    
    fftw_execute( pp );
    
    gamma2_bar.resize(nx*ny);
    
    for( int j=0; j<ny; j++ ){
      int jj = j+jmin;
      for( int i=0; i<nx; i++ ){
        int ii = i+imin;
        
        gamma2_bar[i+nx*j] = float(-realsp[ii+Nnx*jj]/NN);
      }
    }
    //for(auto &a : gamma2_bar) assert(!isnan(a)); // ???
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
    
    fftw_execute( pp );

    for( int j=jmin; j<jmax; j++ ){
      int jj = j-jmin;
      for( int i=imin; i<imax; i++ ){
        int ii = i-imin;
        
        surface_density[ii+nx*jj] = float(realsp[i+Nnx*j]/NN);
        
      }
    }
  }

  // std:: cout << " remapping the map in the original size " << std:: endl;
  delete[] fft;
  delete[] fphi;
  
  //phi_bar.resize(nx*ny,0);  // ??? this needs to be calculated in the future
}

/// no padding
template <typename T>
void LensMap::PreProcessFFTWMap(T Wphi_of_k,bool do_alpha){
  
  assert(surface_density.size() == nx*ny);
  
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
  
  //float *fp = &(surface_density[0]);
  fftw_plan p = fftw_plan_dft_r2c_2d(ny,nx,&(surface_density[0])
                                     ,fphi,FFTW_ESTIMATE);
  
  fftw_execute( p );
  
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
  
  fftw_complex *fft= new fftw_complex[ny*(Nkx)];
  std::vector<double> realsp(NN);
  
  fftw_plan pp = fftw_plan_dft_c2r_2d(ny,nx,fft,realsp.data(),FFTW_MEASURE);
  
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
    
    fftw_execute( pp );
    
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
    
    fftw_execute( pp );
    
    alpha2_bar.resize(NN);
    for( size_t i=0; i<NN; i++ ) alpha2_bar[i] = -1*float(realsp[i]/NN);
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
    
    fftw_execute( pp );

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
    
    fftw_execute( pp );
    
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
    
    fftw_execute( pp );
    
    for( size_t i=0; i<NN; i++ ) surface_density[i] = (realsp[i]/NN);
  }
  
  // std:: cout << " remapping the map in the original size " << std:: endl;
  delete[] fft;
  delete[] fphi;
  
  //phi_bar.resize(NN,0);  // ??? this needs to be calculated in the future
}

#endif
/* MultiLENS_H_ */



