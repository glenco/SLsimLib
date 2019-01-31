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

#include <stdexcept>

#ifdef ENABLE_FITS
#include <CCfits/CCfits>
#endif

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
  
  LensMap():nx(0),ny(0),boxlMpc(0),z(0){};
  
	/// values for the map
	std::valarray<float> surface_density;  // Msun / Mpc^2
	std::valarray<float> alpha1_bar;
	std::valarray<float> alpha2_bar;
	std::valarray<float> gamma1_bar;
	std::valarray<float> gamma2_bar;
	//std::valarray<double> gamma3_bar;
  std::valarray<float> phi_bar;
	//std::vector<double> x;
  int nx,ny;
  // boxlMpc is Mpc/h for MOKA
	/// lens and source properties
  //double zlens,m,zsource,Dlens,DLS,DS,c,cS,fsub,mstar,minsubmass;
  /// range in x direction, pixels are square
  //double boxlarcsec,boxlrad;
  double boxlMpc;
  /// cosmology
  //double omegam,omegal,h,wq;
	//double inarcsec;
	Point_2d center;
  Point_2d lowerleft;
  Point_2d upperright;
  
  double z;
  
#ifdef ENABLE_FITS

  LensMap(std::string fits_input_file,float h){
    read(fits_input_file,h);
  }
  
  /// read an entire map
  void read(std::string input_fits,float h);
  
  /// read only header information
  void read_header(std::string input_fits,float h);

  /// read a subsection of the fits map
  void read_sub(std::string input_fits
                         ,const std::vector<long> &first
                         ,const std::vector<long> &last
                         ,float h
                         );
  
  
  void read_header(std::unique_ptr<CCfits::FITS> ff,float h);

  void read_sub(CCfits::FITS *ff
                ,const std::vector<long> &first
                ,const std::vector<long> &last
                ,float h
                );

  void write(std::string filename);

#endif

#ifdef ENABLE_FFTW
  // this calculates the other lensing quantities from the density map
  //void PreProcessFFTWMap(float zerosize);

  //static double identity(double x){return 1;}
  
  template <class T>
  void PreProcessFFTWMap(float zerosize,T Wphi_of_k);
  //void PreProcessFFTWMap(float zerosize,std::function<double(double)> Wphi_of_k = identity);
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
                   std::string fitsfile  /// Original fits map of the density
                   ,double mass_unit     /// shoudl include h factors
                   ,const COSMOLOGY &c
                   ,bool single_grid_mode = false
                   );

  //double Wlr(double k2){return exp(-k2*rs2);}
    
  ~LensHaloMultiMap(){
    delete ff;
  };
	
  const double f = 5,g = 6;
  
  void submap(
              const std::vector<long> &lower_left
              ,const std::vector<long> &upper_right
              );
  /// in physical coordinates relative to the center of the original map
  void submap(Point_2d ll,Point_2d ur);

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

private:
  
  bool single_grid;
  const COSMOLOGY &cosmo;
  CCfits::FITS *ff;
  
public:  // ?????
  LensMap long_range_map;
  LensMap short_range_map;
private:  // ?????
  double mass_unit;
  
  size_t Noriginal[2]; // number of pixels in each dimension in original image
  double resolution;          // resolution of original image and short rnage image in Mpc
  long border_width;   // width of short range maps padding
  std::string fitsfilename;

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

//void LensMap::PreProcessFFTWMap(float zerosize,std::function<double(double)> Wphi_of_k){
template <typename T>
void LensMap::PreProcessFFTWMap(float zerosize,T Wphi_of_k){
  
  assert(surface_density.size() == nx*ny);

  // size of the new map in x and y directions, factor by which each size is increased
  int Nnx=int(zerosize*nx);
  int Nny=int(zerosize*ny);
  double boxlx = boxlMpc*zerosize;
  double boxly = ny*boxlMpc*zerosize/nx;
  
  int imin = (Nnx-nx)/2;
  int imax = (Nnx+nx)/2;
  int jmin = (Nny-ny)/2;
  int jmax = (Nny+ny)/2;
  
  size_t Nkx = (Nnx/2+1);
  
  std::vector<double> kxs(Nkx);
  for( int i=0; i<Nkx; i++ ){
    kxs[i] = i*2.*M_PI/boxlx;
  }
  std::vector<double> kys(Nny);
  for( int j=0; j<Nny; j++ ){
    kys[j]=(j<Nny/2)?double(j):double(j-Nny);
    kys[j] *= 2.*M_PI/boxly;
  }


  std::vector<double> extended_map( Nnx*Nny );
  
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
        //assert(ii+nx*jj < surface_density.size());
        //assert(i+Nnx*j < extended_map.size());
        extended_map[i+Nnx*j] = surface_density[ii+nx*jj];
        //float tmp = extended_map[i+Nnx*j];
        //assert(!isnan(extended_map[i+Nnx*j]));
        if(isinf(extended_map[i+Nnx*j])){
          extended_map[i+Nnx*j] = 0; /// ????
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
    
    alpha1_bar.resize(nx*ny);
    
    for( int j=jmin; j<jmax; j++ ){
      int jj = j-jmin;
      for( int i=imin; i<imax; i++ ){
        int ii = i-imin;

          alpha1_bar[ii+nx*jj] = -1*float(realsp[i+Nnx*j]/Nnx/Nny);
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
    
    alpha2_bar.resize(nx*ny);
    
    for( int j=jmin; j<jmax; j++ ){
      int jj = j-jmin;
      for( int i=imin; i<imax; i++ ){
        int ii = i-imin;
        
        alpha2_bar[ii+nx*jj] = -1*float(realsp[i+Nnx*j]/Nnx/Nny);
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
    
    gamma1_bar.resize(nx*ny);
    
    for( int j=jmin; j<jmax; j++ ){
      int jj = j-jmin;
      for( int i=imin; i<imax; i++ ){
        int ii = i-imin;
        
        gamma1_bar[ii+nx*jj] = float( realsp[i+Nnx*j]/Nnx/Nny);
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
    
    for( int j=jmin; j<jmax; j++ ){
      int jj = j-jmin;
      for( int i=imin; i<imax; i++ ){
        int ii = i-imin;
        
        gamma2_bar[ii+nx*jj] = float(-realsp[i+Nnx*j]/Nnx/Nny);
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
    
    fftw_execute( pp );

    for( int j=jmin; j<jmax; j++ ){
      int jj = j-jmin;
      for( int i=imin; i<imax; i++ ){
        int ii = i-imin;
        
        surface_density[ii+nx*jj] = float(-realsp[i+Nnx*j]/Nnx/Nny);
        
      }
    }
  }

  // std:: cout << " remapping the map in the original size " << std:: endl;
  delete[] fft;
  delete[] fphi;
  
  phi_bar.resize(nx*ny);  // ??? this needs to be calculated in the future
}

#endif
/* MultiLENS_H_ */



