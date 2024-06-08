/*
 * image_processing.h
 *
 *  Created on: Feb 28, 2010
 *      Author: R.B. Metcalf
 */

#ifndef IMAGE_PROCESSING_H_
#define IMAGE_PROCESSING_H_

#include <complex>
#include <vector>
#include <tuple>

#include "point.h"
#include "image_info.h"
#include "pixelmap.h"
#include "fftw3.h"
//#include "Tree.h"
//#include "utilities_slsim.h"
//#include "utilities_slsim.h"
//#include "image_processing.h"
//#include "source.h"

class Source;
#include "fftw3.h"

// forward declaration
//struct Grid;
//struct GridMap;
//class Source;

enum class Telescope {Euclid_VIS,Euclid_Y,Euclid_J,Euclid_H,KiDS_u,KiDS_g,KiDS_r,KiDS_i,HST_ACS_I,CFHT_u,CFHT_g,CFHT_r,CFHT_i,CFHT_z};

enum class UnitType {counts_x_sec, flux} ;


class Obs{
public:
  
  Obs(size_t Npix_xx,size_t Npix_yy  /// number of pixels in observation
      ,double pix_size               /// pixel size (in rad)
      ,int oversample          /// oversampling for input image
      ,float seeing = 0 // seeing in arcsec
  );
  
  virtual ~Obs(){};
  
  size_t getNxInput() const { return Npix_x_input;}
  size_t getNyInput() const { return Npix_y_input;}

  size_t getNxOutput() const { return Npix_x_output;}
  size_t getNyOutput() const { return Npix_y_output;}

  std::valarray<double> getPSF(){return map_psf;}
  //void setPSF(std::string psf_file);
  void setPSF(std::string psf_file,double resolution=0);
  template <typename T>
  void setPSF(PixelMap<T> &psf_map);
  /// rotate and scale the psf from the original
  void rotatePSF(double theta   /// counter-clockwise rotation (radians)
                 ,double scale_x=1  /// scale <1 shrinks it
                 ,double scale_y=1  /// scale <1 shrinks it
  );
  
  /// add two PSFs to simulate stacking
  void coaddPSF(double f         /// relative weight of the PSFs, 1 being equal weight
                ,double theta1   /// ratation of first PSF
                ,double theta2   /// rotation of second PSF
                ,double scale_x  /// scale <1 shrinks it
                ,double scale_y  /// scale <1 shrinks it
                );
 
  template <typename T>
  void ApplyPSF(PixelMap<T> &map_in,PixelMap<T> &map_out);
  
  float getPixelSize() const {return pix_size;}
  void setNoiseCorrelation(std::string nc_file);
  
  // virtual methods
//  template <typename T>
//  virtual void AddNoise(PixelMap<T> &pmap
//                        ,PixelMap<T> &error_map
//                        ,Utilities::RandomNumbers_NR &ran,bool cosmics) = 0;
//  template <typename T>
//  virtual void Convert(PixelMap<T> &map_in
//                           ,PixelMap<T> &map_out
//                           ,PixelMap<T> &error_map
//                           ,bool psf
//                           ,bool noise
//                           ,Utilities::RandomNumbers_NR &ran,bool cosmics) = 0;
  
  virtual float getBackgroundNoise() const = 0;

  /// convert using stan
  virtual double mag_to_counts(double m) const = 0;
  virtual double counts_to_mag(double flux) const = 0;
  virtual double zeropoint() const = 0;
  virtual void setZeropoint(double zpoint) = 0;


protected:

  double pix_size; // pixel size (in rad)
  template <typename T>
  void CorrelateNoise(PixelMap<T> &pmap);
  float seeing;  // full-width at half maximum of the gaussian smoothing
  
  // the number of pixels in the real image
  size_t Npix_x_output,Npix_y_output;
  // the number of pixels in the oversamples image
  size_t Npix_x_input,Npix_y_input;
 
  float psf_oversample; // psf oversampling factor
  template <typename T>
  void downsample(PixelMap<T> &map_in,PixelMap<T> &map_out) const;  // downsize from Npix_input to Npix_output
  
private:
  double input_psf_pixel_size;
  size_t side_ncorr; // pixels on a side of input noise correlation function

  void fftpsf();  // FFT the psf for later use
  std::valarray<double> map_psf;  // array of the point spread function
  std::valarray<double> map_psfo;  // initial array of the point spread function
  
  std::vector<std::complex<double> > fft_psf;
  std::vector<std::complex<double> > fft_padded;
  std::vector<double> image_padded;
  std::vector<double> sqrt_noise_power;  // stores sqrt root of power noise spectrum

  // size of borders for psf convolution
  size_t nborder_x = 0;
  size_t nborder_y = 0;
  
  // size of padded images
  size_t n_x = 0;
  size_t n_y = 0;

  fftw_plan image_to_fft;
  fftw_plan fft_to_image;

  //PixelMap noise_correlation;
  std::vector<std::complex<double> > noise_fft_image;
  std::vector<double> noise_in_zeropad;
  fftw_plan p_noise_r2c;
  std::vector<double> noise_image_out;
  fftw_plan p_noise_c2r;
};

/**
 * \brief It creates a realistic image from the output of a ray-tracing simulation.
 *
 * It translates pixel values in observed units (counts/sec), applies PSF and noise.
 * Input must be in ergs/(s*cm^2*Hz).
 *
 *  see https://www.ucolick.org/~bolte/AY257/s_n.pdf
 */

class ObsVIS : public Obs{
private:
  
  // standard from magnitude to e- per sec
  double zero_point = 24.4;
  //double sigma_back_per_qsrttime = 0.00267 * sqrt(5.085000000000E+03);
  
  //double gain = 11160; // e-/ADU (Analog Digital Units)
  //double exp_num = 4;
  //double exp_time = 2260.;  // seconds
  //double l = 7103.43;
  //double dl = 3318.28;
  //double seeing = 0.18;
  
  // derived parameters;
  double sigma_background2;  // background variance x time
  //double sb_to_e;  // approximate convertion between ergs / cm^2 / s and e-

  // adds random cosmic rays to the noise map
  template <typename T>
  void cosmics(PixelMap<T> &error_map
                ,double inv_sigma2 // for one dither
                ,int nc // number of cosmics to be added
                ,Utilities::RandomNumbers_NR &ran) const ;
  
public:
  
  // exposure times are set to wide survey expectations
  ObsVIS(size_t Npix_x,size_t Npix_y
         ,int oversample
         ,double resolution = 0.1*arcsecTOradians
         //,double t = 5.085000000000E+03  // observation time in seconds. default is for SC8
  );
  
  ObsVIS(size_t Npix_x
         ,size_t Npix_y
         ,const std::vector<double> &exposure_times  // in seconds
         ,int oversample
         );
  
  ObsVIS(size_t Npix_x
         ,size_t Npix_y
         ,const std::vector<double> &exposure_times  // in seconds
         ,int oversample
         ,double resolution
         ,double background_sigma
         ,double calibration_exposure_time
         );
  ~ObsVIS(){};
  
  /// add poisson noise to an image that is in units of electrons
  template <typename T>
  void AddPoisson(PixelMap<T> &pmap
                         ,Utilities::RandomNumbers_NR &ran
                     );
  
  /// Applies  noise (read-out + Poisson) on an image, returns noise map
  template <typename T>
  void AddNoise(PixelMap<T> &pmap
                 ,PixelMap<T> &error_map
                 ,Utilities::RandomNumbers_NR &ran
                ,bool cosmic=true);

  template <typename T>
  void Convert(PixelMap<T> &map_in
               ,PixelMap<T> &map_out
               ,PixelMap<T> &error_map  // this is sigma
               ,bool psf
               ,bool noise
               ,Utilities::RandomNumbers_NR &ran
               ,bool cosmic=true);
  
 
  double mag_to_counts(double m) const{
    if(m == 100) return 0;
    return pow(10,-0.4*(m + zero_point));
  }
  double counts_to_mag(double flux) const{
    if(flux <=0) return 100;
    return -2.5 * log10(flux) - zero_point;
  }

  double zeropoint() const {return zero_point;}
  void setZeropoint(double zpoint){zero_point=zpoint;}
 
  /// returns std of pixels in e-
  float getBackgroundNoise() const {
    double dt = Utilities::vec_sum(t_exp);// 3 * t1 + t2;
    
    return sqrt( sigma_background2 / dt );
  }
private:
  //double t1;
  //double t2;
  std::vector<double> t_exp;  // exposure times

};

/** 
 * \brief It creates a realistic image from the output of a ray-tracing simulation.
 *
 * It translates pixel values in observed units (counts/sec), applies PSF and noise.
 * Input must be in ergs/(s*cm^2*Hz).
 *
 *  see https://www.ucolick.org/~bolte/AY257/s_n.pdf
 */

class Observation : public Obs
{
public:
	Observation(Telescope tel_name,double exposure_time,int exposure_num
              ,size_t Npix_x,size_t Npix_y, float oversample);
	Observation(float diameter, float transmission, float exp_time, int exp_num, float back_mag, float read_out_noise
              ,size_t Npix_x,size_t Npix_y,double pix_size,float seeing = 0.);
	Observation(float diameter, float transmission, float exp_time, int exp_num, float back_mag, float read_out_noise ,std::string psf_file,size_t Npix_x,size_t Npix_y,double pix_size, float oversample = 1.);
  
  Observation(float zeropoint_mag, float exp_time, int exp_num, float back_mag, float read_out_noise, size_t Npix_x,size_t Npix_y,double pix_size,float seeing=0);
  Observation(float zeropoint_mag, float exp_time, int exp_num, float back_mag, float read_out_noise ,std::string psf_file,size_t Npix_x,size_t Npix_y,double pix_size, float oversample = 1.);
  
  ~Observation(){};
  
  float getExpTime() const {return exp_time;}
	int getExpNum() const {return exp_num;}
	float getBackMag() const {return back_mag;}
   /// read-out noise in electrons/pixel
	float getRon() const {return read_out_noise;}
  /// seeing in arcsecs
	float getSeeing() const {return seeing;}
	float getZeropoint() const {return mag_zeropoint;}
  void setZeropoint(double zpoint){mag_zeropoint=zpoint;}
    /// pixel size in radians
  float getBackgroundNoise(float resolution, UnitType unit = UnitType::counts_x_sec) const;
  float getBackgroundNoise() const {return 0;};

  template <typename T>
  void AddNoise(PixelMap<T> &pmap,PixelMap<T> &error_map
  ,Utilities::RandomNumbers_NR &ran,bool dummy);

  template <typename T>
  void Convert(PixelMap<T> &map_in
               ,PixelMap<T> &map_out
               ,PixelMap<T> &error_map
               ,bool psf
               ,bool noise
               ,Utilities::RandomNumbers_NR &ran
               ,bool cosmic=false
               );
 
  /// returns factor by which code image units need to be multiplied by to get flux units
  //double flux_convertion_factor(){ return pow(10,-0.4*mag_zeropoint); }

  void setExpTime(float time){exp_time = time;}
  void setPixelSize(float pixel_size){pix_size=pixel_size;}
 
  double mag_to_counts(double m) const {
    if(m == 100) return 0;
    return pow(10,-0.4*(m - mag_zeropoint));
  }
  double counts_to_mag(double flux) const{
    if(flux <=0) return 100;
    return -2.5 * log10(flux) + mag_zeropoint;
  }
  double zeropoint() const{
     return mag_zeropoint;
   }
private:
  
	//float diameter;  // diameter of telescope (in cm)
	//float transmission;  // total transmission of the instrument
	float mag_zeropoint;  // magnitude of a source that produces one count/sec in the image
	float exp_time;  // total exposure time (in sec)
	int exp_num;  // number of exposures
	float back_mag;  // sky (or background) magnitude in mag/arcsec^2
	float read_out_noise;  // read-out noise in electrons/pixel
  float gain;
  
	bool telescope; // was the observation created from a default telescope?
  float e_per_s_to_ergs_s_cm2;  // e- / s   for zero magnitudes
  float background_flux;  // e- / s / arcsec

  void set_up();
  
  template <typename T>
  void ToCounts(PixelMap<T> &pmap);
  template <typename T>
  void ToSurfaceBrightness(PixelMap<T> &pmap);
  template <typename T>
  void ToADU(PixelMap<T> &pmap);
};

void pixelize(double *map,long Npixels,double range,double *center
		,ImageInfo *imageinfo,int Nimages,bool constant_sb,bool cleanmap
		,bool write_for_skymaker = false, std::string filename="");
void _SplitFluxIntoPixels(TreeHndl ptree,Branch *leaf,double *leaf_sb);
//void smoothmap(double *map_out,double *map_in,long Npixels,double range,double sigma);

namespace Utilities{
    //void LoadFitsImages(std::string dir,const std::string& filespec,std::vector<PixelMap> & images,int maxN,double resolution = -1,bool verbose = false);
    //void LoadFitsImages(std::string dir,std::vector<std::string> filespecs,std::vector<std::string> file_non_specs                                  ,std::vector<PixelMap> & images,std::vector<std::string> & names,int maxN,double resolution = -1,bool verbose = false);
    void ReadFileNames(std::string dir,const std::string filespec
                       ,std::vector<std::string> & filenames
                       ,const std::string file_non_spec = " "
                       ,bool verbose = false);
}

/** \brief Warning: Not tested yet. Class for doing adaptive smoothing using multiply resolution grids.
 */
class MultiGridSmoother{
public:
  MultiGridSmoother(double center[],std::size_t Nx,std::size_t Ny,double resolution);
  MultiGridSmoother(double center[],std::size_t Nx,double resolution);
  ~MultiGridSmoother(void){
    maps.clear();
  }
  
  /// resolution of finest grid from which interpolation is done
  PosType getHighestRes(){return maps[0].getResolution();}
  /// resolution of coarsest grid from which interpolation is done
  PosType getLowestRes(){return maps.back().getResolution();}
  
  /// Add particles to the map.  These do not need to be kept in memory after they are added.
  void add_particles(std::vector<PosType> x,std::vector<PosType> y);
  /// Output a map at the resolution of the map smoothed so that no superpixel as less than Nsmooth particles
  void output_map(PixelMap<double> &map,int Nsmooth);
  void smooth(int Nsmooth,PixelMap<double> &map);
  
private:
  void _smooth_(int k,size_t i,size_t j,int Nsmooth,PixelMap<double> &map);
  std::vector<PixelMap<double> > maps;
  std::vector<Utilities::Interpolator<PixelMap<double> >> interpolators;
};


template <typename T>
void ObsVIS::AddPoisson(PixelMap<T> &pmap
                        ,Utilities::RandomNumbers_NR &ran
                        ){
   
  double dt =  Utilities::vec_sum(t_exp);
  if(ran() < 0.2){
    // select a frame
    int missing_frame = (int)(ran()*t_exp.size());
    
    dt -= t_exp[missing_frame];
  }
  
  for(auto &p : pmap.data() ){
    if(p>0) p = ran.poisson(p*dt)/dt;
  }
}

template <typename T>
void ObsVIS::AddNoise(PixelMap<T> &pmap
                      ,PixelMap<T> &error_map
                      ,Utilities::RandomNumbers_NR &ran,bool cosmic
                      ){
  if(pmap.getNx() != pmap.getNy()){
    std::cerr << "Observation::AddNoise() Doesn't work on nonsquare maps" << std::endl;
    throw std::runtime_error("nonsquare");
  }
  
  double dt = Utilities::vec_sum(t_exp);
  int missing_frame = -1;
  //int drop=0;
  if(ran() < 0.2){
    // select a frame
    missing_frame = (int)(ran()*t_exp.size());
    dt -= t_exp[missing_frame];
    //drop=1;
  }

  double inv_sigma2 = (dt)/sigma_background2;
  size_t N = pmap.size();
  for (unsigned long i = 0; i < N ; i++){
    error_map[i] = inv_sigma2;
  }
  
  //double p = (t_exp.size()-drop)*t1/dt;
  if(cosmic){
    
    double t = (int)(ran()*dt);
    for(int j=0 ; j<100 ; ++j){
      // chooses an exposure at random weighted by exposure time
      int i=0;
      double tt=t_exp[0];
      while(t>tt){
        ++i;
        tt += t_exp[i];
      }
    
      if(i != missing_frame) cosmics(error_map,t_exp[i]/sigma_background2,1,ran);
    }
    
//    if(ran() < p ){
//      cosmics(error_map,t1/sigma2,ran.poisson(100),ran);
//    }else{
//      cosmics(error_map,t2/sigma2,ran.poisson(100*(dt-t1)/dt),ran);
//    }
  }
  for (unsigned long i = 0; i < N ; i++){
    pmap[i] = ran.poisson(MAX<float>(pmap[i] * dt,0))/dt;
    pmap[i] += ran.gauss() / sqrt(error_map[i]);
    
    error_map[i] = sqrt( 1.0 / error_map[i] + MAX<float>(pmap[i] / dt,0) ) ;
  }

  return;
}

template <typename T>
void  ObsVIS::cosmics(
                      PixelMap<T> &error_map  // inverse error
                      ,double inv_sigma2 // for one dither
                      ,int nc // number of cosmics to be added
                      ,Utilities::RandomNumbers_NR &ran
                      ) const {
  size_t N = error_map.size();
  size_t Nx = error_map.getNx();
  size_t Ny = error_map.getNy();
  
  int length = 2;
  size_t ic=0;
  
  while (ic < nc) {
  
    size_t no = (size_t)(N*ran());
    double theta = 2*PI*ran();
    double c=cos(theta),s=sin(theta);
  
    size_t io = no % Nx;
    size_t jo = no / Nx;
 
    double L = MAX<double>(length/tan(PI*ran()/2),1);
    //error_map.DrawLine(io, io + L * c, jo, jo + L * s, -inv_sigma2, true); // ???
    //error_map.DrawLine(io + 0.5*s, io + L * c + 0.5*s, jo - 0.5*c, jo + L * s - 0.5*c, -inv_sigma2, true);

    double t = tan(theta);
    long x,y;
    if(abs(t) < 1){
      long x1 = (long)(io + L * c );
      int sgn = sign(c);
      x1 = MAX<long>(0,x1);
      x1 = MIN<long>(Nx-1,x1);
      for(x = io ; x !=  x1 ; x = x + sgn){
        y = MAX<long>((long)(t*(x-io) + jo),0);
        y= MIN<long>(y,Ny-2);
        error_map(x,y) += -inv_sigma2;
        error_map(x,y+1) += -inv_sigma2;
      }
    }else{
      long y1 = (long)(jo + L * s );
      int sgn = sign(s);
      y1 = MAX<long>(0,y1);
      y1 = MIN<long>(Ny-1,y1);
      for(y = jo ; y !=  y1 ; y = y + sgn){
        x = MAX<long>((long)((y-jo)/t + io),0);
        x = MIN<long>(x,Nx-2);
        error_map(x,y) += -inv_sigma2;
        error_map(x+1,y) += -inv_sigma2;
      }
    }
    ++ic;
  }
}

template <typename T>
void ObsVIS::Convert(
              PixelMap<T> &map_in
             ,PixelMap<T> &map_out
             ,PixelMap<T> &error_map
             ,bool psf
             ,bool noise
             ,Utilities::RandomNumbers_NR &ran
             ,bool cosmic
                     ){
  assert(map_in.getNx() == Npix_x_input);
  assert(map_in.getNy() == Npix_y_input);
  
  if (fabs(map_in.getResolution()*psf_oversample - pix_size) > pix_size*1.0e-5)
  {
    std::cout << "The resolution of the input map is different from the one of the simulated instrument in Observation::Convert!" << std::endl;
    throw std::runtime_error("The resolution of the input map is different from the one of the simulated instrument!");
  }

  map_out.Clean();
  
  if (psf == true){
    PixelMap<T> map_scratch(Point_2d(0,0).x
                            , Npix_x_input
                            , Npix_y_input, pix_size);
    ApplyPSF<T>(map_in,map_scratch);
    downsample<T>(map_scratch,map_out);
  }else{
    downsample<T>(map_in,map_out);
  }
 
  if (noise == true) AddNoise<T>(map_out,error_map,ran,cosmic);
  
  return;
}


/// Reads in and sets the PSF from a fits file. If the pixel size of the fits is different (smaller) than the one of the telescope, it must be specified.
template <typename T>
void Obs::setPSF(PixelMap<T> &psf_map/// name of fits file with psf
                 ){
 
  input_psf_pixel_size = psf_map.getResolution();
  
  if( (input_psf_pixel_size - pix_size/psf_oversample)/input_psf_pixel_size > 1.0e-3){
    std::cout << "Obs::setPSF() - psf is not resolved." << std::endl;
    throw std::runtime_error("");
  }
  
  map_psf.resize(psf_map.size());
  for(size_t i=0 ; i<psf_map.size() ; ++i){
    map_psf[i] = psf_map[i];
  }
  
  std::vector<size_t> size = {psf_map.getNx(),psf_map.getNy()};
  
  long max_x,max_y;
  {
    long i=0,imax=0;
    double amax = map_psf[0];
    for(auto &a : map_psf){
      if(a > amax){imax=i;amax=a;}
      ++i;
    }
    max_y = imax % size[0];
    max_x = imax / size[0];
  }
  long Lx;
  long Ly = Lx = 2*MIN<double>( MIN<double>(max_y,size[1]-max_y)
                               ,MIN<double>(max_x,size[0]-max_x));
  
  std::valarray<double> tmp(Lx*Ly);
  long ii=0;
  for(long i=max_x - Lx/2 ; i<max_x + Lx/2 ;++i){
    assert(i<size[0]);
    long jj=0;
    for(long j=max_y - Ly/2 ; j<max_y + Ly/2 ;++j){
      assert(j<size[1]);
      tmp[ii + Lx*jj] = map_psf[j + size[1] * i];
      ++jj;
    }
    ++ii;
  }
  
  std::swap(tmp,map_psf);
  map_psfo=map_psf;
  
  fftpsf();
}

template <typename T>
void Obs::downsample(PixelMap<T> &map_in,PixelMap<T> &map_out) const{
  
  assert(map_in.getNx() == Npix_x_input);
  assert(map_in.getNy() == Npix_y_input);
  assert(map_out.getNx() == Npix_x_output);
  assert(map_out.getNy() == Npix_y_output);
 
  if(psf_oversample == 1){
    size_t n=map_in.size();
    for(size_t i=0; i<n ; ++i ) map_out[i] = map_in[i]; // keep bounding box information
    return;
  }
 
  map_out.Clean();
//  Point_2d x;
//  long n = map_in.size();
//  for(long i=0 ; i<n ; ++i){
//    map_in.find_position(x.x, i);
//    long k = map_out.find_index(x.x);
//    if(k > -1) map_out[k] += map_in[i];
//  }
  
  for(size_t i=0 ; i<Npix_x_input ; ++i){
    size_t ii = MIN<long>(i / psf_oversample + 0.5, Npix_x_output - 1 ) ;
    for(size_t j=0 ; j<Npix_y_input ; ++j){
      size_t jj = MIN<long>(j / psf_oversample + 0.5, Npix_y_output - 1 ) ;

      map_out(ii,jj) += map_in(i,j);
    }
  }

  return;
}

/** * \brief Smooths the image with a PSF map.
*
*/
template <typename T>
void Obs::ApplyPSF(PixelMap<T> &map_in,PixelMap<T> &map_out)
{
  if(map_in.getNx() != map_in.getNy()){
    std::cout << "Obs::ApplyPSF() Doesn't work on nonsquare maps" << std::endl;
    throw std::runtime_error("nonsquare");
  }
  if(map_in.getNx() != Npix_y_input){
    std::cout << "Obs::ApplyPSF() Input map it the wrong size." << std::endl;
    throw std::runtime_error("nonsquare");
  }
  if(map_in.getNx() != map_out.getNx()){
    std::cout << "Obs::ApplyPSF() Input map output are not the same size." << std::endl;
    throw std::runtime_error("nonsquare");
  }

  
  if (map_psf.size() == 0){
    
    map_out = map_in;
    
    if (seeing > 0.){
      //PixelMap outmap(pmap);
      map_out.smooth(seeing/2.355);
      return;
    }else{
      return;
    }
    
  }else{
    map_out.ChangeUnits(map_in.getUnits());
    
    // paste image into image with padding
    for(double &a :image_padded) a=0;
    for(int i=0 ; i<Npix_x_input ; ++i){
      for(int j=0 ; j<Npix_y_input ; ++j){
        image_padded[ (i+nborder_x)*n_x + j + nborder_y] = map_in(i,j);
      }
    }
 
    fftw_execute(image_to_fft);
    
    // multiply by DFT of PSF
    assert(fft_psf.size() == fft_padded.size());
    size_t i=0;
    for(std::complex<double> &a : fft_padded) a *= fft_psf[i++];
    
    // inverse DFT
    fftw_execute(fft_to_image);
    
    // copy region within padding
    size_t N = n_x*n_y;
    for(int i=0 ; i< Npix_x_input ; ++i){
      for(int j=0 ; j< Npix_y_input ; ++j){
        map_out(i,j) = image_padded[ (i+nborder_x)*n_x + j + nborder_y]/N;
      }
    }

    return;

  }
}

/** * \brief Correlate a noise map
 *
 */
template <typename T>
void Obs::CorrelateNoise(PixelMap<T> &pmap)
{
  
  if(sqrt_noise_power.size()==0) return;
  
  if(pmap.getNx() != pmap.getNy()){
    std::cout << "Observation::CorrelateNoise() Doesn't work on nonsquare maps" << std::endl;
    throw std::runtime_error("nonsquare");
  }
  if(pmap.getNx() != Npix_x_output){
    std::cout << "Observation::CorrelateNoise() Map must have the same dimensions as the observation." << std::endl;
    throw std::runtime_error("nonsquare");
  }
  
    // creates plane for fft of map, sets properly input and output data, then performs fft
    assert(Npix_x_output == Npix_y_output);
    size_t Npix = Npix_x_output;
    size_t ncorr_big_zeropad_Npixels = Npix + side_ncorr;
    long Npix_zeropad = Npix + side_ncorr;
  
    size_t fftsize = sqrt_noise_power.size();
    
    // rows and columns between first_p and last_p are copied in the zero-padded version
    long first_p = side_ncorr/2;
    long last_p = first_p + (Npix-1);
    
    // add zero-padding
     for (int i = 0; i < Npix_zeropad*Npix_zeropad; i++)
    {
      long ix = i/Npix_zeropad;
      long iy = i%Npix_zeropad;
      if (ix >= first_p && ix <= last_p && iy >= first_p && iy <= last_p)
        noise_in_zeropad[i] = pmap[(ix-side_ncorr/2)*Npix+(iy-side_ncorr/2)];
      else
        noise_in_zeropad[i] = 0.;
    }
    
     fftw_execute(p_noise_r2c);
  
    // performs convolution in Fourier space , and transforms back to real space
    for (unsigned long i = 0; i < fftsize ; i++)
    {
      size_t ix = i/(Npix_zeropad/2+1);
      size_t iy = i%(Npix_zeropad/2+1);
      if (ix>Npix_zeropad/2)
        noise_fft_image[i] *= sqrt_noise_power[(ncorr_big_zeropad_Npixels-(Npix_zeropad-ix))*(ncorr_big_zeropad_Npixels/2+1)+iy];
      else
        noise_fft_image[i] *= sqrt_noise_power[ix*(ncorr_big_zeropad_Npixels/2+1)+iy];
    }
  
  // creats plane for perform backward fft after convolution, sets output data
  //double* image_out = new double[Npix_zeropad*Npix_zeropad];
  
    fftw_execute(p_noise_c2r);
    
    // translates array of data in (normalised) counts map
    for (unsigned long i = 0; i < Npix_zeropad*Npix_zeropad; i++)
    {
      size_t ix = i/Npix_zeropad;
      size_t iy = i%Npix_zeropad;
      if (ix >= first_p && ix <= last_p && iy >= first_p && iy <= last_p)
      {
        int ii = (ix-side_ncorr/2)*Npix + (iy-side_ncorr/2);
        pmap[ii] = noise_image_out[i]/(Npix_zeropad*Npix_zeropad);
      }
    }
    return;
}


/**  \brief Converts the input map to a realistic image
 *
 * \param map Input map in photons/(cm^2*Hz)
 * \param psf Decides if the psf smoothing is applied
 * \param noise Decides if noise is added
 * \param unit Decides units of output (if flux, output is in 10**(-0.4*mag))
 */
template <typename T>
void Observation::Convert(PixelMap<T> &map_in
                          ,PixelMap<T> &map_out
                          ,PixelMap<T> &error_map
                          ,bool psf
                          ,bool noise
                          ,Utilities::RandomNumbers_NR &ran
,bool cosmic)
{
  
  assert(map_in.getNx() == Npix_x_input);
  assert(map_in.getNy() == Npix_y_input);

  if(cosmic) throw std::runtime_error("cosmics not implemented");
  if (fabs(map_in.getResolution()*psf_oversample - pix_size) > pix_size*1.0e-5)
  {
    std::cout << "The resolution of the input map is different from the one of the simulated instrument in Observation::Convert!" << std::endl;
    throw std::runtime_error("The resolution of the input map is different from the one of the simulated instrument!");
  }
  
  map_out.Clean();
  if (psf == true){
    PixelMap<T> map_scratch(Point_2d(0,0).x
                            , Npix_x_input
                            , Npix_y_input, pix_size);
    ApplyPSF<T>(map_in,map_scratch);
    downsample<T>(map_scratch,map_out);
  }else{
    downsample<T>(map_in,map_out);
  }
  ToCounts(map_out);
 
  if (noise == true) AddNoise<T>(map_out,error_map,ran,true);
  //ToSurfaceBrightness(map_out);
  
  return;
}


/// Applies realistic noise (read-out + Poisson) on an image, returns noise map
template <typename T>
void Observation::AddNoise(
                           PixelMap<T> &pmap
                           ,PixelMap<T> &error_map
                           ,Utilities::RandomNumbers_NR &ran,bool dummy)
{
  if(pmap.getNx() != pmap.getNy()){
    std::cerr << "Observation::AddNoise() Doesn't work on nonsquare maps" << std::endl;
    throw std::runtime_error("nonsquare");
  }
  if(pmap.getUnits() != PixelMapUnits::count_per_sec){
    std::cerr << "Units need to be in counts per second in Observation::AddNoise." << std::endl;
    throw std::runtime_error("wrong units.");
  }
  
  double res_in_arcsec = pmap.getResolution() / arcsecTOradians;
  double back_mean = background_flux * res_in_arcsec*res_in_arcsec * exp_time ;

  double var_readout = exp_num*read_out_noise*read_out_noise;
  double norm_map;
  
  PixelMap<T> noise_map(pmap);
  //double sum=0,sum2=0;
  size_t N = pmap.size();
  for (unsigned long i = 0; i < N ; i++)
  {
    norm_map = pmap[i]*exp_time;  // in counts
    if (norm_map+back_mean > 500.)
    {
      double sdv = sqrt(var_readout + norm_map + back_mean);
      noise_map[i] = ran.gauss()*sdv/exp_time;  // back to counts per sec
      error_map[i] = sdv*sdv/exp_time/exp_time;
    }
    else
    {
      // photons from mean
      int k = ran.poisson(norm_map + back_mean);
      double noise = ran.gauss()*sqrt(var_readout);
      
      noise_map[i] = (k + noise - back_mean - norm_map)/exp_time;
      error_map[i] = (back_mean + norm_map + var_readout)/exp_time;
    }
  }

  CorrelateNoise(noise_map);
  pmap += noise_map;
  
  return;
}

/// Translates photon flux (in 1/(s*cm^2*Hz*hplanck)) into telescope counts per second
template <typename T>
void Observation::ToCounts(PixelMap<T> &pmap)
{
  
  //zero_point_flux = pow(10,-0.4*mag_zeropoint);  // erg/s/Hz/cm**2
  //background_flux = pow(10,-0.4*(back_mag-mag_zeropoint ));

  double Q;
  PixelMapUnits units = pmap.getUnits();
  if(units == PixelMapUnits::count_per_sec) return;

  if(units == PixelMapUnits::surfb){
    Q = e_per_s_to_ergs_s_cm2;
  }else if(pmap.getUnits() == PixelMapUnits::ADU){
    Q = 1.0/gain;
  }else{
    std::cerr << "Map needs to be in photon flux units." << std::endl;
    throw std::runtime_error("wrong units");
  }
  
  pmap.Renormalize(Q);
  pmap.ChangeUnits(PixelMapUnits::count_per_sec);
  return;
}

/// Translates photon flux (in 1/(s*cm^2*Hz*hplanck)) into telescope counts per second
template <typename T>
void Observation::ToSurfaceBrightness(PixelMap<T> &pmap)
{

  double Q;
  if(pmap.getUnits() == PixelMapUnits::count_per_sec){
    Q = 1.0/e_per_s_to_ergs_s_cm2 ;
  }else if(pmap.getUnits() == PixelMapUnits::ADU){
    Q = 1.0/gain/e_per_s_to_ergs_s_cm2;
  }else{
    std::cerr << "Map needs to be in photon flux units." << std::endl;
    throw std::runtime_error("wrong units");
  }
  
  pmap.Renormalize(Q);
  pmap.ChangeUnits(PixelMapUnits::surfb);
  return;
}

/// Translates photon flux (in 1/(s*cm^2*Hz*hplanck)) into telescope counts per second
template <typename T>
void Observation::ToADU(PixelMap<T> &pmap)
{
  
  ToCounts(pmap);
  pmap *= gain;
  pmap.ChangeUnits(PixelMapUnits::ADU);
  
  return;
}

#endif
