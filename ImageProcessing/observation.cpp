/*
 * observation.cpp
 *
 */

#include <complex>
//#include "slsimlib.h"
#include "image_processing.h"

#include "cpfits.h"
#include "fftw3.h"

#include <fstream>
#include <limits>


ObsVIS::ObsVIS(size_t Npix_x,size_t Npix_y,int oversample,double resolution)
:Obs(Npix_x,Npix_y,resolution,oversample,0)
{
  //sigma_background = 0.00365150 * sqrt(2366.);

  //sb_to_e = (119.*119.*PI/4.) * t * dl / l / hplanck;
  
  //t1 = 565;
  //t2 = 106;

  // new values 16/1/24
  //t1 = 560;
  //t2 = 89.5;
  t_exp = {560,560,560,560,89.5,89.5};
  sigma_background2 = 0.002 * 0.002 * Utilities::vec_sum(t_exp) ;
}

ObsVIS::ObsVIS(size_t Npix_x,size_t Npix_y
               ,const std::vector<double> &exposure_times  // in seconds
               ,int oversample
               )
:Obs(Npix_x,Npix_y,0.1*arcsecTOradians,oversample,0),t_exp(exposure_times)
{
  sigma_background2 = 0.0015 * 0.0015 * ( Utilities::vec_sum(t_exp) );
}

ObsVIS::ObsVIS(size_t Npix_x,size_t Npix_y
               ,const std::vector<double> &exposure_times  // in seconds
               ,int oversample
               ,double resolution
               ,double my_background_sigma
               ,double calibration_exposure_time)
:Obs(Npix_x,Npix_y,resolution,oversample,0),t_exp(exposure_times)
{
  sigma_background2 = my_background_sigma * my_background_sigma * calibration_exposure_time  ;
}


Obs::Obs(size_t Npix_xx,size_t Npix_yy  /// number of pixels in observation
    ,double pix_size               /// pixel size (in rad)
    ,int oversample          /// oversampling for input image
    ,float my_seeing // seeing in arcsec
    ):
  pix_size(pix_size)
,seeing(my_seeing)
,Npix_x_output(Npix_xx)
,Npix_y_output(Npix_yy)
,psf_oversample(oversample)
{
  Npix_x_input = oversample * Npix_x_output;
  Npix_y_input = oversample * Npix_y_output;
}

/// Reads in and sets the PSF from a fits file. If the pixel size of the fits is different (smaller) than the one of the telescope, it must be specified.
void Obs::setPSF(std::string psf_file  /// name of fits file with psf
                  ,double resolution  /// resolution in degrees if it isn't in fits header
                 )
{
 
  CPFITS_READ cpfits(psf_file);

  int exits = cpfits.readKey("CD1_1",input_psf_pixel_size);  // this is in degrees
  
  if(exits != 0){
    if(resolution > 0){
      std::cout << "No resolution found in PSF.  Using input resolution" << std::endl;
      input_psf_pixel_size = resolution;
    }else{
    std::cerr << " PSF " << psf_file << " does not have key word CD1_1" << std::endl;
    throw std::runtime_error("bad file");
    }
  }
  std::cout << "Obs::setPSF() - intrinsic psf resolution : " <<
  input_psf_pixel_size*60*60 << " arcsec" << std::endl;
  input_psf_pixel_size *= degreesTOradians;
  
  if( (input_psf_pixel_size - pix_size/psf_oversample)/input_psf_pixel_size > 1.0e-3){
    std::cout << "Obs::setPSF() - psf is not resolved." << std::endl;
    std::cout << (input_psf_pixel_size - pix_size/psf_oversample)/input_psf_pixel_size << std::endl;
    throw std::runtime_error("");
  }
                           
  std::vector<long> size;
  cpfits.read(map_psf, size);
  
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
  long Ly = Lx = 2*MIN( MIN(max_y,size[1]-max_y),MIN(max_x,size[0]-max_x));
  
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
  
  map_psfo = map_psf;
  
  fftpsf();
}

void Obs::rotatePSF(double theta,double scale_x,double scale_y){
  
  double s=-sin(theta);
  double c=cos(theta);
  
  int N_psf = map_psf.size();
  int n_side_psf = sqrt(N_psf);

  long center[2] = {n_side_psf/2,n_side_psf/2};
  double f[2];
  
  std::valarray<double> map_psf(N_psf);
  
  for(long k=0 ; k<N_psf ; ++k){
    long i = k%n_side_psf-center[0];
    long j = k/n_side_psf-center[1];
    
    double x = (i*c - j*s)/scale_x + center[0];
    if(x>=0 && x < n_side_psf-1){
      double y = (i*s + j*c)/scale_y + center[1];
      if(y>=0 && y < n_side_psf-1){
        
        long kk = (long)(x) + (long)(y)*n_side_psf; // lower left
        
        f[0]= x - (long)(x);
        f[1]= y - (long)(y);
        
        map_psf[k] = (1-f[0])*(1-f[1])*map_psfo[kk] + f[0]*(1-f[1])*map_psfo[kk+1]
        + f[0]*f[1]*map_psfo[kk+1+n_side_psf]
        + (1-f[0])*f[1]*map_psfo[kk+n_side_psf];
      }else{
        map_psf[k] = 0;
      }
    }else{
      map_psf[k] = 0;
    }
  }
  
  fftpsf();
}

void Obs::coaddPSF(double f,double theta1,double theta2,double scale_x,double scale_y){
  
  double a = f/(1+f);
  double b = 1/(1+f);
  double s1=-sin(theta1);
  double c1=cos(theta1);
  
  double s2=-sin(theta2);
  double c2=cos(theta2);

  int N_psf = map_psf.size();
  int n_side_psf = sqrt(N_psf);

  long center[2] = {n_side_psf/2,n_side_psf/2};
  double ff[2];
  
  std::valarray<double> map_psf(N_psf);
  
  for(long k=0 ; k<N_psf ; ++k){
    long i = k%n_side_psf-center[0];
    long j = k/n_side_psf-center[1];
    
    double x = (i*c1 - j*s1)/scale_x + center[0];
    if(x>=0 && x < n_side_psf-1){
      double y = (i*s1 + j*c1)/scale_y + center[1];
      if(y>=0 && y < n_side_psf-1){
        
        long kk = (long)(x) + (long)(y)*n_side_psf; // lower left
        
        ff[0]= x - (long)(x);
        ff[1]= y - (long)(y);
        
        map_psf[k] = a*( (1-ff[0])*(1-ff[1])*map_psfo[kk] + ff[0]*(1-ff[1])*map_psfo[kk+1]
        + ff[0]*ff[1]*map_psfo[kk+1+n_side_psf]
        + (1-ff[0])*ff[1]*map_psfo[kk+n_side_psf] );
      }else{
        map_psf[k] = 0;
      }
    }else{
      map_psf[k] = 0;
    }
    
    x = (i*c2 - j*s2)/scale_x + center[0];
    if(x>=0 && x < n_side_psf-1){
      double y = (i*s2 + j*c2)/scale_y + center[1];
      if(y>=0 && y < n_side_psf-1){
        
        long kk = (long)(x) + (long)(y)*n_side_psf; // lower left
        
        ff[0]= x - (long)(x);
        ff[1]= y - (long)(y);
        
        map_psf[k] += b*( (1-ff[0])*(1-ff[1])*map_psfo[kk] + ff[0]*(1-ff[1])*map_psfo[kk+1]
        + ff[0]*ff[1]*map_psfo[kk+1+n_side_psf]
        + (1-ff[0])*ff[1]*map_psfo[kk+n_side_psf] );
      }
    }
  }
  fftpsf();
}


/// Read in and set the noise correlation function.
void Obs::setNoiseCorrelation(std::string nc_file  /// name of fits file with noise correlation function in pixel units
)
{
  std::cout << nc_file << std::endl;
  PixelMap<double> noise_corr(nc_file,pix_size);
  size_t N = noise_corr.size();
  double sum=0,corr_max=0;
  // take the square root of the correlation function
  for(size_t i = 0 ; i < N ; ++i){
    corr_max = MAX(corr_max,noise_corr[i]);
    sum += noise_corr[i];
    //noise_corr[i] = sqrt(fabs(noise_corr[i]));
  }
  sum = sqrt(sum);
  for(size_t i = 0 ; i < N ; ++i) noise_corr[i] /= sum;
  
  std::cout << "zero lag noise correlation is : " << corr_max
  << " size : " << noise_corr.size() << std::endl;
  
  //PixelMap<T>::swap(noise_corr,nc_map);
  // now find the power spectrum of the noise
  
  // calculates normalisation of ncor
  int N_ncorr = noise_corr.size();
  side_ncorr = sqrt(N_ncorr);
  double map_norm = 0.;
  for (int i = 0; i < N_ncorr; i++)
  {
    map_norm += noise_corr[i];
  }
  fftw_plan p_ncorr;
  
  // arrange ncorr data for fft, creates plane, then performs fft
  // ncorr data are moved into the four corners of ncorr_big_zeropad
  
  assert(Npix_x_output == Npix_y_output);
  size_t Npix = Npix_x_output;
  size_t ncorr_big_zeropad_Npixels = Npix + side_ncorr;
  size_t ncorr_big_zeropad_Npixels2 = ncorr_big_zeropad_Npixels*ncorr_big_zeropad_Npixels;
  
  std::vector<double> ncorr_big_zeropad(ncorr_big_zeropad_Npixels2);
  std::vector<std::complex<double> > out_ncorr(ncorr_big_zeropad_Npixels*(ncorr_big_zeropad_Npixels/2+1));
  
  size_t fftsize = ncorr_big_zeropad_Npixels*(ncorr_big_zeropad_Npixels/2+1);
  sqrt_noise_power.resize(fftsize);
  
  p_ncorr = fftw_plan_dft_r2c_2d(ncorr_big_zeropad_Npixels,ncorr_big_zeropad_Npixels,ncorr_big_zeropad.data(), reinterpret_cast<fftw_complex*>(out_ncorr.data()), FFTW_ESTIMATE);
  long ix, iy;
  for (int i = 0; i < ncorr_big_zeropad_Npixels*ncorr_big_zeropad_Npixels; i++)
  {
    ix = i/ncorr_big_zeropad_Npixels;
    iy = i%ncorr_big_zeropad_Npixels;
    if(ix<side_ncorr/2 && iy<side_ncorr/2)
      ncorr_big_zeropad[i] = noise_corr[(ix+side_ncorr/2)*side_ncorr+(iy+side_ncorr/2)]/map_norm;
    else if(ix<side_ncorr/2 && iy>=ncorr_big_zeropad_Npixels-side_ncorr/2)
      ncorr_big_zeropad[i] = noise_corr[(ix+side_ncorr/2)*side_ncorr+(iy-(ncorr_big_zeropad_Npixels-side_ncorr/2))]/map_norm;
    else if(ix>=ncorr_big_zeropad_Npixels-side_ncorr/2 && iy<side_ncorr/2)
      ncorr_big_zeropad[i] = noise_corr[(ix-(ncorr_big_zeropad_Npixels-side_ncorr/2))*side_ncorr+(iy+side_ncorr/2)]/map_norm;
    else if(ix>=ncorr_big_zeropad_Npixels-side_ncorr/2 && iy>=ncorr_big_zeropad_Npixels-side_ncorr/2)
      ncorr_big_zeropad[i] = noise_corr[(ix-(ncorr_big_zeropad_Npixels-side_ncorr/2))*side_ncorr+(iy-(ncorr_big_zeropad_Npixels-side_ncorr/2))]/map_norm;
    else
      ncorr_big_zeropad[i] = 0.;
  }
  fftw_execute(p_ncorr);

  fftw_destroy_plan(p_ncorr);
  
  // normalized so that later the variance will be conserved
  sum = 0.0;
  for(size_t i = 0 ; i < sqrt_noise_power.size() ; ++i)
    sum += out_ncorr[i].real();
  sum /= sqrt_noise_power.size();
  for(size_t i = 0 ; i < sqrt_noise_power.size() ; ++i)
    sqrt_noise_power[i] = sqrt( out_ncorr[i].real()/sum );
  
  
  noise_fft_image.resize(fftsize);
  noise_in_zeropad.resize(ncorr_big_zeropad_Npixels2);
  p_noise_r2c = fftw_plan_dft_r2c_2d(ncorr_big_zeropad_Npixels,ncorr_big_zeropad_Npixels
                                     ,noise_in_zeropad.data()
                                     , reinterpret_cast<fftw_complex*>(noise_fft_image.data())
                                     , FFTW_ESTIMATE);

  noise_image_out.resize(ncorr_big_zeropad_Npixels2);
  p_noise_c2r = fftw_plan_dft_c2r_2d(ncorr_big_zeropad_Npixels,ncorr_big_zeropad_Npixels
                                     ,reinterpret_cast<fftw_complex*>(noise_fft_image.data())
                                     , noise_image_out.data()
                                     , FFTW_ESTIMATE);
}


void Obs::fftpsf(){
  
  // calculates normalisation of psf
  int N_psf = map_psf.size();
  int n_side_psf = sqrt(N_psf);

  nborder_x = Npix_x_input/2;
  nborder_y = Npix_y_input/2;

  n_x = Npix_x_input + 2*nborder_x;
  n_y = Npix_y_input + 2*nborder_y;
  
  double oversample_factor = pix_size / input_psf_pixel_size / psf_oversample;
  
  // make extended map of psf
  std::vector<double> psf_padded(n_x * n_y,0);
  std::vector<int> psf_count(n_x * n_y,0);

  // find maximum of psf
  long half_psf_x,half_psf_y;
  {
    long i=0,imax=0;
    double amax =map_psf[0];
    for(auto &a : map_psf){
      if(a > amax){imax=i;amax=a;}
      ++i;
    }
    half_psf_x = imax % n_side_psf;
    half_psf_y = imax / n_side_psf;
  }
  
  
  
  long n_side_psf_x = MIN<long>(2*half_psf_x , n_side_psf);
  long n_side_psf_y = MIN<long>(2*half_psf_y , n_side_psf);
 
  // shift center of psf to bottom left with a rap and down sample
  //long half_psf = n_side_psf/2;
  for(long i=0 ; i< n_side_psf_x ; ++i){
    size_t ii = (i >= half_psf_x) ? (i - half_psf_x)/oversample_factor + 0.5 :
                                   n_x + (i - half_psf_x)/oversample_factor + 0.5;
    if(ii >=n_x) ii=n_x-1;
    for(long j=0 ; j< n_side_psf_y ; ++j){
      size_t jj = (j >= half_psf_y) ? (j - half_psf_y)/oversample_factor + 0.5 :
                                   n_y + (j - half_psf_y)/oversample_factor + 0.5;
      
      if(jj >= n_y) jj=n_y-1;
      psf_padded[ii*n_x + jj] += map_psf[i*n_side_psf + j];
      psf_count[ii*n_x + jj] += 1;
    }
  }
  
  double psf_norm = 0.;
  for(long i=0;i<psf_count.size();++i){
    if(psf_count[i] >0){
      psf_padded[i] /= psf_count[i];
      psf_norm += psf_padded[i];
    }else{
      psf_padded[i] = 0;
    }
  }
  
  for(double &p : psf_padded) p /= psf_norm;
  
  fft_psf.resize(n_x*(n_y/2+1));
  fft_padded.resize(n_x*(n_y/2+1));
  image_padded.resize(n_x*n_y);
  for(double &a :image_padded) a=0;
  
  fftw_plan p_psf = fftw_plan_dft_r2c_2d(n_x ,n_y ,psf_padded.data()
                               ,reinterpret_cast<fftw_complex*>(fft_psf.data()), FFTW_ESTIMATE);
  fftw_execute(p_psf);
  
  //std::complex<double> phase = std::polar(1, i*half_psf_x / n_x + j*half_psf_y / n_y ); *****
  
  image_to_fft = fftw_plan_dft_r2c_2d(n_x,n_y ,image_padded.data()
                                      ,reinterpret_cast<fftw_complex*>(fft_padded.data()), FFTW_ESTIMATE);

  fft_to_image = fftw_plan_dft_c2r_2d(n_x,n_y,reinterpret_cast<fftw_complex*>(fft_padded.data())
                                      , image_padded.data(), FFTW_ESTIMATE);
}


/** * \brief Creates an observation setup that mimics a known instrument
 *
 */
Observation::Observation(Telescope tel_name
                         ,double exposure_time,int exposure_num
                         ,size_t Npix_x,size_t Npix_y, float oversample):
  Obs(Npix_x,Npix_y,0,oversample),exp_time(exposure_time),exp_num(exposure_num)
{
 
  float diameter,transmission;
  switch (tel_name) {
//    case Telescope::Euclid_VIS:
//      // from Eric
//      // Equivalent gain is 11160 e-/ADU (Analog Digital Units)
//      // Saturation is 71 ADU
//      // Equivalent exposure time with 4 frames is 2260s
//      // Magnitude Zeropoint is 23.9
//
//      gain = 11160;
//      //exp_time = 1800.;
//      //exp_time = 2260.;
//      //exp_num = 4;
//      mag_zeropoint = 23.9;
//
//      back_mag = 22.8;  //back_mag = 25.0; // ?????
//      read_out_noise = 5.; //read_out_noise = 1.0; //read_out_noise = 0.01; // ?????
//      seeing = 0.18;
//      pix_size = 0.1*arcsecTOradians;
//
//      break;
//    case Telescope::Euclid_Y:
//      diameter = 119.;
//      transmission = 0.0961;
//
//      gain = 0;
//      //exp_time = 264.;
//      //exp_num = 3;
//      back_mag = 22.57;
//      read_out_noise = 5.;
//      seeing = 0.3;
//      pix_size = .3*arcsecTOradians;
//
//      //mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
//      //mag_zeropoint = mag_to_counts(1.0/(diameter*diameter*transmission*PI/4.));
//
//      mag_zeropoint = 25.04;
//
//      break;
//    case Telescope::Euclid_J:
//      diameter = 119.;
//      transmission = 0.0814;
//
//      gain = 0;
//      //exp_time = 270.;
//      //exp_num = 3;
//      back_mag = 22.53;
//      read_out_noise = 5.;
//      seeing = 0.3;
//      pix_size = .3*arcsecTOradians;
//
//      //mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
//      //mag_zeropoint = mag_to_counts(1.0/(diameter*diameter*transmission*PI/4.));
//
//      mag_zeropoint = 25.26;
//      break;
//    case Telescope::Euclid_H:
//      diameter = 119.;
//      transmission = 0.1692;
//      gain = 0;
//      //exp_time = 162.;
//      //exp_num = 3;
//      back_mag = 22.59;
//      read_out_noise = 5.;
//      seeing = 0.3;
//      pix_size = .3/60./60./180.*PI;
//
//      //mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;   // convert from flux to magnitudes
//      //mag_zeropoint = mag_to_counts(1.0/(diameter*diameter*transmission*PI/4.));
//
//      mag_zeropoint = 25.21;
//      break;
    case Telescope::KiDS_u:
      diameter = 265.;
      transmission = 0.032;
      gain = 0;
      //exp_time = 1000.;
      //exp_num = 5;
      back_mag = 22.93;
      read_out_noise = 5.;
      seeing = 1.0;
      pix_size = .2/60./60./180.*PI;
            
      //mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      mag_zeropoint = mag_to_counts(1.0/(diameter*diameter*transmission*PI/4.));
 
      break;
    case Telescope::KiDS_g:
      diameter = 265.;
      transmission = 0.1220;
      gain = 0;
      //exp_time = 900.;
      //exp_num = 5;
      back_mag = 22.29;
      read_out_noise = 5.;
      seeing = 0.8;
      pix_size = .2/60./60./180.*PI;
      
      //mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      mag_zeropoint = mag_to_counts(1.0/(diameter*diameter*transmission*PI/4.));
      break;
    case Telescope::KiDS_r:
      diameter = 265.;
      transmission = 0.089;
      gain = 0;
      //exp_time = 1800.;
      //exp_num = 5;
      back_mag = 21.40;
      read_out_noise = 5.;
      seeing = 0.7;
      pix_size = .2/60./60./180.*PI;
      
      //mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      mag_zeropoint = mag_to_counts(1.0/(diameter*diameter*transmission*PI/4.));
      break;
    case Telescope::KiDS_i:
      diameter = 265.;
      transmission = 0.062;
      gain = 0;
      //exp_time = 1200.;
      //exp_num = 5;
      back_mag = 20.64;
      read_out_noise = 5.;
      seeing = 1.1;
      pix_size = .2/60./60./180.*PI;
      
      //mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      mag_zeropoint = mag_to_counts(1.0/(diameter*diameter*transmission*PI/4.));
      break;
    case Telescope::HST_ACS_I:
      diameter = 250.;
      transmission = 0.095;
      //exp_time = 420.;
      //exp_num = 1;
      back_mag = 22.8;
      read_out_noise = 3.;
      seeing = 0.1;
      pix_size = .05/60./60./180.*PI;
      
      //mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      mag_zeropoint = mag_to_counts(1.0/(diameter*diameter*transmission*PI/4.));
      break;
    case Telescope::CFHT_u:
      diameter = 358.;
      transmission = 0.0644;
      //exp_time = 3000.;
      //exp_num = 5;
      back_mag = 22.7;
      read_out_noise = 5.;
      seeing = 0.85;
      pix_size = .187/60./60./180.*PI;
      
      //mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      mag_zeropoint = mag_to_counts(1.0/(diameter*diameter*transmission*PI/4.));
      break;
    case Telescope::CFHT_g:
      diameter = 358.;
      transmission = 0.1736;
      //exp_time = 2500.;
      //exp_num = 5;
      back_mag = 22.0;
      read_out_noise = 5.;
      seeing = 0.78;
      pix_size = .187/60./60./180.*PI;
      
      //mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      mag_zeropoint = mag_to_counts(1.0/(diameter*diameter*transmission*PI/4.));
      break;
    case Telescope::CFHT_r:
      diameter = 358.;
      transmission = 0.0971;
      //exp_time = 2000.;
      //exp_num = 4;
      back_mag = 21.3;
      read_out_noise = 5.;
      seeing = 0.71;
      pix_size = .187/60./60./180.*PI;
      
      //mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      mag_zeropoint = mag_to_counts(1.0/(diameter*diameter*transmission*PI/4.));
      break;
    case Telescope::CFHT_i:
      diameter = 358.;
      transmission = 0.0861;
      //exp_time = 4300.;
      //exp_num = 5;
      back_mag = 20.3;
      read_out_noise = 5.;
      seeing = 0.64;
      pix_size = .187/60./60./180.*PI;
      
      //mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      mag_zeropoint = mag_to_counts(1.0/(diameter*diameter*transmission*PI/4.));
      break;
    case Telescope::CFHT_z:
      diameter = 358.;
      transmission = 0.0312;
      //exp_time = 3600.;
      //exp_num = 6;
      back_mag = 19.4;
      read_out_noise = 5.;
      seeing = 0.68;
      pix_size = .187/60./60./180.*PI;
      
      //mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      mag_zeropoint = mag_to_counts(1.0/(diameter*diameter*transmission*PI/4.));
      break;
    default:
      throw std::runtime_error("The Telescope selected is not available.");
      break;
	}

 	telescope = true;
  set_up();
}

/** *  Creates a custom observation setup with parameters decided by the user.
 *
 * \param diameter Diameter of telescope (in cm) (Collecting area is pi/4.*diameter^2)
 * \param transmission Total transmission of the telescope (/int T(\lambda)/\lambda d\lambda)
 * \param exp_time Total exposure time
 * \param exp_num Number of exposures
 * \param back_mag Flux due to background and/or sky in mag/arcsec^2
 * \param read_out_noise Read-out noise in electrons/pixel
 * \param seeing FWHM in arcsecs of the image
 */
Observation::Observation(float diameter, float transmission, float my_exp_time, int my_exp_num
, float my_back_mag, float my_read_out_noise
, size_t Npix_x,size_t Npix_y,double my_pix_size,float my_seeing):
    Obs(Npix_x,Npix_y,my_pix_size,1,my_seeing),exp_time(my_exp_time), exp_num(my_exp_num)
    , back_mag(my_back_mag),read_out_noise(my_read_out_noise)
{
  mag_zeropoint = mag_to_counts(1.0/(diameter*diameter*transmission*PI/4.));
  telescope = false;
 
  set_up();
}

/**  Creates a custom observation setup with parameters decided by the user. Allows for the use of a psf fits image.
 *
 * \param diameter Diameter of telescope (in cm) (Collecting area is pi/4.*diameter^2)
 * \param transmission Total transmission of the telescope (/int T(\lambda)/\lambda d\lambda)
 * \param exp_time Total exposure time
 * \param exp_num Number of exposures
 * \param back_mag Flux due to background and/or sky in mag/arcsec^2
 * \param read_out_noise Read-out noise in electrons/pixel
 * \param psf_file Input PSF image
 * \param oversample Oversampling rate of the PSF image
 */
Observation::Observation(float diameter, float transmission, float my_exp_time, int my_exp_num
, float my_back_mag, float ron, std::string psf_file,size_t Npix_x,size_t Npix_y
,double my_pix_size, float oversample):
    Obs(Npix_x,Npix_y,my_pix_size,oversample), exp_time(my_exp_time), exp_num(my_exp_num)
    , back_mag(my_back_mag) , read_out_noise(0)
		{
      //mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      mag_zeropoint = mag_to_counts(1.0/(diameter*diameter*transmission*PI/4.));
 
      setPSF(psf_file);

      telescope = false;
      
      set_up();
}

/**  Creates a custom observation setup with parameters decided by the user.
 *
 * \param zeropoint_mag Magnitude that produces one count per second on the detector
 * \param exp_time Total exposure time
 * \param exp_num Number of exposures
 * \param back_mag Flux due to background and/or sky in mag/arcsec^2
 * \param read_out_noise Read-out noise in electrons/pixel
 * \param seeing FWHM in arcsecs of the image
 */
Observation::Observation(float zeropoint_mag, float exp_time, int exp_num, float back_mag, float read_out_noise, size_t Npix_x,size_t Npix_y,double pix_size,float seeing):Obs(Npix_x,Npix_y,pix_size,1,seeing)
    ,mag_zeropoint(zeropoint_mag),exp_time(exp_time), exp_num(exp_num), back_mag(back_mag)
    ,read_out_noise(read_out_noise)
{
  telescope = false;
 
  set_up();
}

/**  Creates a custom observation setup with parameters decided by the user. Allows for the use of a psf fits image.
 *
 * \param zeropoint_mag Magnitude that produces one count per second on the detector
 * \param exp_time Total exposure time
 * \param exp_num Number of exposures
 * \param back_mag Flux due to background and/or sky in mag/arcsec^2
 * \param read_out_noise Read-out noise in electrons/pixel
 * \param psf_file Input PSF image
 * \param oversample Oversampling rate of the PSF image
 */
Observation::Observation(float zeropoint_mag, float exp_time, int exp_num, float back_mag, float myread_out_noise, std::string psf_file,size_t Npix_x,size_t Npix_y,double pix_size, float oversample):Obs(Npix_x,Npix_y,pix_size,oversample), mag_zeropoint(zeropoint_mag), exp_time(exp_time), exp_num(exp_num), back_mag(back_mag), read_out_noise(myread_out_noise)
    {
      setPSF(psf_file);
      telescope = false;
      
      set_up();
}



/// Outputs rms of noise counts due to background and instrument
/// in the unit decided by the user
float Observation::getBackgroundNoise(float resolution, UnitType unit)
const {
    if (telescope==true && fabs(resolution-pix_size) > pix_size*1.0e-5)
    {
        std::cout << "The resolution is different from the one of the simulated instrument in Observation::getBackgroundNoise!" << std::endl;
      throw std::runtime_error("The resolution is different from the one of the simulated instrument!");
    }

    double res_in_arcsec = resolution / arcsecTOradians ;
    double back_mean = background_flux * res_in_arcsec*res_in_arcsec * exp_time ;

    // in e-
    double rms = sqrt(exp_num*read_out_noise*read_out_noise + back_mean) / exp_time;
    if (unit==UnitType::flux) rms /= e_per_s_to_ergs_s_cm2;
  
    return rms;
}

void Observation::set_up(){
  e_per_s_to_ergs_s_cm2 = 1.0/(mag_to_counts(mag_zeropoint));
  background_flux = pow(10,-0.4*(back_mag-mag_zeropoint ));
 }
