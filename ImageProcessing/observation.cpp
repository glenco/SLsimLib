/*
 * observation.cpp
 *
 */

#include <complex>
//#include "slsimlib.h"
#include "image_processing.h"

#include "cpfits.h"

#ifdef ENABLE_FFTW
#include "fftw3.h"
#endif

#include <fstream>
#include <limits>


ObsVIS::ObsVIS(size_t Npix_x,size_t Npix_y,int oversample,double resolution)
:Obs(Npix_x,Npix_y,resolution,oversample,1)
{
  //sigma_background = 0.00365150 * sqrt(2366.);

  //sb_to_e = (119.*119.*PI/4.) * t * dl / l / hplanck;
  
  //t1 = 565;
  //t2 = 106;

  // new values 16/1/24
  //t1 = 560;
  //t2 = 89.5;
  t_exp = {560,560,560,560,89.5,89.5};
  sigma_background = 0.002 * sqrt( Utilities::vec_sum(t_exp) );

  //sigma_background = 0.002 * sqrt(4*t1 + 2*t2);
}

ObsVIS::ObsVIS(size_t Npix_x,size_t Npix_y
               ,const std::vector<double> &exposure_times  // in seconds
               ,int oversample
               )
:Obs(Npix_x,Npix_y,0.1*arcsecTOradians,oversample,1),t_exp(exposure_times)
{
  sigma_background = 0.0015 * sqrt( Utilities::vec_sum(t_exp) );
}

ObsVIS::ObsVIS(size_t Npix_x,size_t Npix_y
               ,const std::vector<double> &exposure_times  // in seconds
               ,int oversample
               ,double resolution
               ,double my_background_sigma)
:Obs(Npix_x,Npix_y,resolution,oversample,1),t_exp(exposure_times)
{
  sigma_background = my_background_sigma * sqrt( Utilities::vec_sum(t_exp) );
}


void ObsVIS::AddPoisson(PixelMap &pmap
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

void ObsVIS::AddNoise(PixelMap &pmap
                      ,PixelMap &error_map
                      ,Utilities::RandomNumbers_NR &ran,bool cosmic
                      ){
  if(pmap.getNx() != pmap.getNy()){
    std::cerr << "Observation::AddNoise() Doesn't work on nonsquare maps" << std::endl;
    throw std::runtime_error("nonsquare");
  }
  
  double sigma2 = sigma_background*sigma_background;
  //  + pmap[i] / sb_to_e );
  
  double dt = Utilities::vec_sum(t_exp);
  int missing_frame = -1;
  //int drop=0;
  if(ran() < 0.2){
    // select a frame
    missing_frame = (int)(ran()*t_exp.size());
    dt -= t_exp[missing_frame];
    //drop=1;
  }

  double inv_sigma2 = (dt)/sigma2;
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
    
      if(i != missing_frame) cosmics(error_map,t_exp[i]/sigma2,1,ran);
    }
    
//    if(ran() < p ){
//      cosmics(error_map,t1/sigma2,ran.poisson(100),ran);
//    }else{
//      cosmics(error_map,t2/sigma2,ran.poisson(100*(dt-t1)/dt),ran);
//    }
  }
  for (unsigned long i = 0; i < N ; i++){
    error_map[i] = sqrt( 1.0 / error_map[i] + MAX<float>(pmap[i] / dt,0) ) ;
    pmap[i] += ran.gauss() * error_map[i];
  }

  return;
}

void  ObsVIS::cosmics(
                      PixelMap &error_map  // inverse error
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

void ObsVIS::Convert(
              PixelMap &map_in
             ,PixelMap &map_out
             ,PixelMap &error_map
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
    ApplyPSF(map_in,map_scratch);
    downsample(map_scratch,map_out);
  }else{
    downsample(map_in,map_out);
  }
 
  if (noise == true) AddNoise(map_out,error_map,ran,cosmic);
  
  return;
}

Obs::Obs(size_t Npix_xx,size_t Npix_yy  /// number of pixels in observation
    ,double pix_size               /// pixel size (in rad)
    ,int oversample          /// oversampling for input image
    ,float seeing // seeing in arcsec
    ):
  pix_size(pix_size)
,seeing(seeing)
,Npix_x_output(Npix_xx)
,Npix_y_output(Npix_yy)
,psf_oversample(oversample)
,map_scratch(Point_2d(0,0).x, oversample * Npix_xx,  oversample * Npix_yy, pix_size){
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


/// Reads in and sets the PSF from a fits file. If the pixel size of the fits is different (smaller) than the one of the telescope, it must be specified.
void Obs::setPSF(PixelMap &psf_map/// name of fits file with psf
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

/// Read in and set the noise correlation function.
void Obs::setNoiseCorrelation(std::string nc_file  /// name of fits file with noise correlation function in pixel units
)
{
  std::cout << nc_file << std::endl;
  PixelMap noise_corr(nc_file,pix_size);
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
  
  //PixelMap::swap(noise_corr,nc_map);
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

void Obs::downsample(PixelMap &map_in,PixelMap &map_out) const{
  
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

/** * \brief Smooths the image with a PSF map.
*
*/
void Obs::ApplyPSF(PixelMap &map_in,PixelMap &map_out)
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
void Obs::CorrelateNoise(PixelMap &pmap)
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

#ifdef ENABLE_FFTW
  
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
#else
    std::cerr << "Please enable the preprocessor flag ENABLE_FFTW !" << std::endl;
    exit(1);
#endif
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
    case Telescope::Euclid_VIS:
      // from Eric
      // Equivalent gain is 11160 e-/ADU (Analog Digital Units)
      // Saturation is 71 ADU
      // Equivalent exposure time with 4 frames is 2260s
      // Magnitude Zeropoint is 23.9
      
      gain = 11160;
      //exp_time = 1800.;
      //exp_time = 2260.;
      //exp_num = 4;
      mag_zeropoint = 23.9;
      
      back_mag = 22.8;  //back_mag = 25.0; // ?????
      read_out_noise = 5.; //read_out_noise = 1.0; //read_out_noise = 0.01; // ?????
      seeing = 0.18;
      pix_size = 0.1*arcsecTOradians;
      
      break;
    case Telescope::Euclid_Y:
      diameter = 119.;
      transmission = 0.0961;
      
      gain = 0;
      //exp_time = 264.;
      //exp_num = 3;
      back_mag = 22.57;
      read_out_noise = 5.;
      seeing = 0.3;
      pix_size = .3*arcsecTOradians;
      
      //mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      //mag_zeropoint = mag_to_counts(1.0/(diameter*diameter*transmission*PI/4.));
      
      mag_zeropoint = 25.04;
      
      break;
    case Telescope::Euclid_J:
      diameter = 119.;
      transmission = 0.0814;
      
      gain = 0;
      //exp_time = 270.;
      //exp_num = 3;
      back_mag = 22.53;
      read_out_noise = 5.;
      seeing = 0.3;
      pix_size = .3*arcsecTOradians;
      
      //mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      //mag_zeropoint = mag_to_counts(1.0/(diameter*diameter*transmission*PI/4.));
      
      mag_zeropoint = 25.26;
      break;
    case Telescope::Euclid_H:
      diameter = 119.;
      transmission = 0.1692;
      gain = 0;
      //exp_time = 162.;
      //exp_num = 3;
      back_mag = 22.59;
      read_out_noise = 5.;
      seeing = 0.3;
      pix_size = .3/60./60./180.*PI;
      
      //mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;   // convert from flux to magnitudes
      //mag_zeropoint = mag_to_counts(1.0/(diameter*diameter*transmission*PI/4.));
       
      mag_zeropoint = 25.21;
      break;
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
Observation::Observation(float diameter, float transmission, float exp_time, int exp_num, float back_mag, float read_out_noise, size_t Npix_x,size_t Npix_y,double pix_size,float seeing):
    Obs(Npix_x,Npix_y,pix_size,1,seeing),exp_time(exp_time), exp_num(exp_num), back_mag(back_mag),read_out_noise(read_out_noise)
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
Observation::Observation(float diameter, float transmission, float exp_time, int exp_num, float back_mag, float ron, std::string psf_file,size_t Npix_x,size_t Npix_y,double pix_size, float oversample):
    Obs(Npix_x,Npix_y,pix_size,oversample), exp_time(exp_time), exp_num(exp_num), back_mag(back_mag) , read_out_noise(read_out_noise)
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


/**  \brief Converts the input map to a realistic image
 *
 * \param map Input map in photons/(cm^2*Hz)
 * \param psf Decides if the psf smoothing is applied
 * \param noise Decides if noise is added
 * \param unit Decides units of output (if flux, output is in 10**(-0.4*mag)) 
 */
void Observation::Convert(PixelMap &map_in
                          ,PixelMap &map_out
                          ,PixelMap &error_map
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
    ApplyPSF(map_in,map_scratch);
    downsample(map_scratch,map_out);
  }else{
    downsample(map_in,map_out);
  }
  ToCounts(map_out);
 
  if (noise == true) AddNoise(map_out,error_map,ran,true);
  //ToSurfaceBrightness(map_out);
  
  return;
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
  //zero_point_flux = pow(10,-0.4*mag_zeropoint);  // erg/s/Hz/cm**2
  //e_per_s_to_ergs_s_cm2 = pow(10,0.4*(mag_zeropoint-AB_zeropoint));
  
  e_per_s_to_ergs_s_cm2 = 1.0/(mag_to_counts(mag_zeropoint));
  background_flux = pow(10,-0.4*(back_mag-mag_zeropoint ));
 }

/// Applies realistic noise (read-out + Poisson) on an image, returns noise map
void Observation::AddNoise(
                           PixelMap &pmap
                           ,PixelMap &error_map
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
  
  PixelMap noise_map(pmap);
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
void Observation::ToCounts(PixelMap &pmap)
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
void Observation::ToSurfaceBrightness(PixelMap &pmap)
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
void Observation::ToADU(PixelMap &pmap)
{
  
  ToCounts(pmap);
  pmap *= gain;
  pmap.ChangeUnits(PixelMapUnits::ADU);
  
  return;
}
