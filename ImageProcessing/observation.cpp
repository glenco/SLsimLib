/*
 * observation.cpp
 *
 *  Created on: Apr 29, 2013
 *      Author: F. Bellagamba
 */

#include <complex>
#include "slsimlib.h"

#include "cpfits.h"

#ifdef ENABLE_FFTW
#include "fftw3.h"
#endif

#include <fstream>
#include <limits>


/** * \brief Creates an observation setup that mimics a known instrument
 *
 */
Observation::Observation(Telescope tel_name,size_t Npix_x,size_t Npix_y):
Npix_x(Npix_x),Npix_y(Npix_y)
{
 
  float diameter,transmission;
  switch (tel_name) {
    case Euclid_VIS:
      // from Eric
      // Equivalent gain is 11160 e-/ADU (Analog Digital Units)
      // Saturation is 71 ADU
      // Equivalent exposure time with 4 frames is 2260s
      // Magnitude Zeropoint is 23.9
     
      //exp_time = 1800.;
      exp_time = 2260.;
      exp_num = 4;
      mag_zeropoint = 23.9;
      
      back_mag = 22.8;  //back_mag = 25.0; // ?????
      ron = 5.; //ron = 0.01; // ?????
      seeing = 0.18;
      pix_size = .1*arcsecTOradians;
      break;
   case Euclid_Y:
      diameter = 119.;
      transmission = 0.0961;
      
      
      exp_time = 264.;
      exp_num = 3;
      back_mag = 22.57;
      ron = 5.;
      seeing = 0.3;
      pix_size = .3*arcsecTOradians;
      
      mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      break;
    case Euclid_J:
      diameter = 119.;
      transmission = 0.0814;
      
      exp_time = 270.;
      exp_num = 3;
      back_mag = 22.53;
      ron = 5.;
      seeing = 0.3;
      pix_size = .3*arcsecTOradians;
      
      mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      break;
    case Euclid_H:
      diameter = 119.;
      transmission = 0.1692;
      exp_time = 162.;
      exp_num = 3;
      back_mag = 22.59;
      ron = 5.;
      seeing = 0.3;
      pix_size = .3/60./60./180.*PI;
      
      mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;   // convert from flux to magnitudes

      break;
    case KiDS_u:
      diameter = 265.;
      transmission = 0.032;
      exp_time = 1000.;
      exp_num = 5;
      back_mag = 22.93;
      ron = 5.;
      seeing = 1.0;
      pix_size = .2/60./60./180.*PI;
            
      mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      break;
    case KiDS_g:
      diameter = 265.;
      transmission = 0.1220;
      exp_time = 900.;
      exp_num = 5;
      back_mag = 22.29;
      ron = 5.;
      seeing = 0.8;
      pix_size = .2/60./60./180.*PI;
      
      mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      break;
    case KiDS_r:
      diameter = 265.;
      transmission = 0.089;
      exp_time = 1800.;
      exp_num = 5;
      back_mag = 21.40;
      ron = 5.;
      seeing = 0.7;
      pix_size = .2/60./60./180.*PI;
      
      mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      break;
    case KiDS_i:
      diameter = 265.;
      transmission = 0.062;
      exp_time = 1200.;
      exp_num = 5;
      back_mag = 20.64;
      ron = 5.;
      seeing = 1.1;
      pix_size = .2/60./60./180.*PI;
      
      mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      break;
    case HST_ACS_I:
      diameter = 250.;
      transmission = 0.095;
      exp_time = 420.;
      exp_num = 1;
      back_mag = 22.8;
      ron = 3.;
      seeing = 0.1;
      pix_size = .05/60./60./180.*PI;
      
      mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      break;
    case CFHT_u:
      diameter = 358.;
      transmission = 0.0644;
      exp_time = 3000.;
      exp_num = 5;
      back_mag = 22.7;
      ron = 5.;
      seeing = 0.85;
      pix_size = .187/60./60./180.*PI;
      
      mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      break;
    case CFHT_g:
      diameter = 358.;
      transmission = 0.1736;
      exp_time = 2500.;
      exp_num = 5;
      back_mag = 22.0;
      ron = 5.;
      seeing = 0.78;
      pix_size = .187/60./60./180.*PI;
      
      mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      break;
    case CFHT_r:
      diameter = 358.;
      transmission = 0.0971;
      exp_time = 2000.;
      exp_num = 4;
      back_mag = 21.3;
      ron = 5.;
      seeing = 0.71;
      pix_size = .187/60./60./180.*PI;
      
      mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
      break;
    case CFHT_i:
      diameter = 358.;
      transmission = 0.0861;
      exp_time = 4300.;
      exp_num = 5;
      back_mag = 20.3;
      ron = 5.;
      seeing = 0.64;
      pix_size = .187/60./60./180.*PI;
      
      mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
    break;
    case CFHT_z:
      diameter = 358.;
      transmission = 0.0312;
      exp_time = 3600.;
      exp_num = 6;
      back_mag = 19.4;
      ron = 5.;
      seeing = 0.68;
      pix_size = .187/60./60./180.*PI;
      
      mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
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
 * \param ron Read-out noise in electrons/pixel
 * \param seeing FWHM in arcsecs of the image
 */
Observation::Observation(float diameter, float transmission, float exp_time, int exp_num, float back_mag, float ron, size_t Npix_x,size_t Npix_y,float seeing):
Npix_x(Npix_x),Npix_y(Npix_y),exp_time(exp_time), exp_num(exp_num), back_mag(back_mag),ron(ron),seeing(seeing)
{
  mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;
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
 * \param ron Read-out noise in electrons/pixel
 * \param psf_file Input PSF image
 * \param oversample Oversampling rate of the PSF image
 */
Observation::Observation(float diameter, float transmission, float exp_time, int exp_num, float back_mag, float ron, std::string psf_file,size_t Npix_x,size_t Npix_y, float oversample):
Npix_x(Npix_x),Npix_y(Npix_y), exp_time(exp_time), exp_num(exp_num), back_mag(back_mag) , ron(ron), oversample(oversample)
		{
      mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) + AB_zeropoint;

      CPFITS_READ cpfits(psf_file);

      std::vector<long> sizes;
      cpfits.read(map_psf,sizes);
      //int side_psf = h0->axis(0);
      //int side_psf = sizes[0];
      //int N_psf = side_psf*side_psf;
      //map_psf.resize(N_psf);
      //h0->read(map_psf);
      telescope = false;
      
      set_up();
}

/**  Creates a custom observation setup with parameters decided by the user.
 *
 * \param zeropoint_mag Magnitude that produces one count per second on the detector
 * \param exp_time Total exposure time
 * \param exp_num Number of exposures
 * \param back_mag Flux due to background and/or sky in mag/arcsec^2
 * \param ron Read-out noise in electrons/pixel
 * \param seeing FWHM in arcsecs of the image
 */
Observation::Observation(float zeropoint_mag, float exp_time, int exp_num, float back_mag, float ron, size_t Npix_x,size_t Npix_y,float seeing):
Npix_x(Npix_x),Npix_y(Npix_y),mag_zeropoint(zeropoint_mag),exp_time(exp_time), exp_num(exp_num), back_mag(back_mag)
,ron(ron),seeing(seeing)
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
 * \param ron Read-out noise in electrons/pixel
 * \param psf_file Input PSF image
 * \param oversample Oversampling rate of the PSF image
 */
Observation::Observation(float zeropoint_mag, float exp_time, int exp_num, float back_mag, float ron, std::string psf_file,size_t Npix_x,size_t Npix_y, float oversample):
Npix_x(Npix_x),Npix_y(Npix_y), mag_zeropoint(zeropoint_mag), exp_time(exp_time), exp_num(exp_num), back_mag(back_mag)
, ron(ron), oversample(oversample)
    {
      CPFITS_READ cpfits(psf_file);

      std::vector<long> sizes;
      cpfits.read(map_psf,sizes);
      //int side_psf = h0->axis(0);
      //int side_psf = sizes[0];
      //int N_psf = side_psf*side_psf;
      //map_psf.resize(N_psf);
      //h0->read(map_psf);
      telescope = false;
      
      set_up();
}

/// Reads in and sets the PSF from a fits file. If the pixel size of the fits is different (smaller) than the one of the telescope, it must be specified.
void Observation::setPSF(std::string psf_file  /// name of fits file with psf
                         , float os  /// over sampling factor
                         )
{
 
  CPFITS_READ cpfits(psf_file);

  std::vector<long> size;
  cpfits.read(map_psf, size);

  //int side_psf = h0->axis(0);
  //int side_psf = size[0];
	//int N_psf = side_psf*side_psf;
	//map_psf.resize(N_psf);
	//h0->read(map_psf);

  oversample = os;
  
  fftpsf();
}

/// Read in and set the noise correlation function.
void Observation::setNoiseCorrelation(std::string nc_file  /// name of fits file with noise correlation function in pixel units
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
  
  assert(Npix_x == Npix_y);
  size_t Npix = Npix_x;
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

/**  \brief Converts the input map to a realistic image
 *
 * \param map Input map in photons/(cm^2*Hz)
 * \param psf Decides if the psf smoothing is applied
 * \param noise Decides if noise is added
 * \param unit Decides units of output (if flux, output is in 10**(-0.4*mag)) 
 */
PixelMap Observation::Convert(PixelMap &map, bool psf, bool noise,Utilities::RandomNumbers_NR &ran)//long *seed)
{
  
  if (telescope == true && fabs(map.getResolution()-pix_size) > pix_size*1.0e-5)
  {
    std::cout << "The resolution of the input map is different from the one of the simulated instrument in Observation::Convert!" << std::endl;
    throw std::runtime_error("The resolution of the input map is different from the one of the simulated instrument!");
  }
  
  PixelMap error_map;
  
  ToCounts(map);
  if (psf == true)  ApplyPSF(map);
  if (noise == true) error_map = AddNoise(map,ran);
  ToSurfaceBrightness(map);
  
  //if (unit == flux && map.getUnits() != photon_flux )
  //{
    //double counts_to_flux = pow(10,-0.4*mag_zeropoint);
    //map.Renormalize(zero_point_flux);
    //map.ChangeUnits(photon_flux);
  //}
  return error_map;
}

void Observation::fftpsf(){
  
  // calculates normalisation of psf
  int N_psf = map_psf.size();
  int n_side_psf = sqrt(N_psf);
  double map_norm = 0.;
  for (int i = 0; i < N_psf; i++)
  {
    map_norm += map_psf[i];
  }

  nborder_x = Npix_x/2;
  nborder_y = Npix_y/2;

  n_x = Npix_x + 2*nborder_x;
  n_y = Npix_y + 2*nborder_y;
  
  // make extended map of psf
  std::vector<double> psf_padded(n_x * n_y);
  for(double &a : psf_padded) a=0;
  
  // shift center of psf to bottom left whish a rap
  long half_psf = n_side_psf/2;
  for(long i=0 ; i< n_side_psf ; ++i){
    size_t ii = (i >= half_psf) ? (i - half_psf)/oversample : n_x + (i - half_psf)/oversample;
    for(long j=0 ; j< n_side_psf ; ++j){
      size_t jj = (j >= half_psf) ? (j - half_psf)/oversample : n_y + (j - half_psf)/oversample;
      psf_padded[ii*n_x + jj] += map_psf[i*n_side_psf + j]/map_norm;
    }
  }
  
  fft_psf.resize(n_x*(n_y/2+1));
  fft_padded.resize(n_x*(n_y/2+1));
  image_padded.resize(n_x*n_y);
  for(double &a :image_padded) a=0;
  
  fftw_plan p_psf;
  p_psf = fftw_plan_dft_r2c_2d(n_x ,n_y ,psf_padded.data()
                               ,reinterpret_cast<fftw_complex*>(fft_psf.data()), FFTW_ESTIMATE);
  fftw_execute(p_psf);
  
  image_to_fft = fftw_plan_dft_r2c_2d(n_x,n_y ,image_padded.data()
                                      ,reinterpret_cast<fftw_complex*>(fft_padded.data()), FFTW_ESTIMATE);

  fft_to_image = fftw_plan_dft_c2r_2d(n_x,n_y,reinterpret_cast<fftw_complex*>(fft_padded.data())
                                      , image_padded.data(), FFTW_ESTIMATE);
}

/** * \brief Smooths the image with a PSF map.
*
*/
void Observation::ApplyPSF(PixelMap &pmap)
{
  if(pmap.getNx() != pmap.getNy()){
    std::cout << "Observation::AddNoise() Doesn't work on nonsquare maps" << std::endl;
    throw std::runtime_error("nonsquare");
  }
  
	if (map_psf.size() == 0)
	{
		if (seeing > 0.)
		{
			//PixelMap outmap(pmap);
			pmap.smooth(seeing/2.355);
			return;
		}
		else
		{
			return;
		}
	}
	else
	{

#ifdef ENABLE_FFTW
    
    
    // paste image into image with padding
    for(double &a :image_padded) a=0;
    for(int i=0 ; i<Npix_x ; ++i){
      for(int j=0 ; j<Npix_y ; ++j){
        image_padded[(i+nborder_x)*n_x + j+nborder_y] = pmap(i,j);
      }
    }
 
    fftw_execute(image_to_fft);
    
    assert(fft_psf.size() == fft_padded.size());
    size_t i=0;
    for(std::complex<double> &a : fft_padded) a *= fft_psf[i++];
    
    fftw_execute(fft_to_image);
    
    size_t N = n_x*n_y;
    for(int i=0 ; i<Npix_x ; ++i){
      for(int j=0 ; j<Npix_y ; ++j){
        pmap(i,j) = image_padded[ (i+nborder_x)*n_x + j+nborder_y ]/N;
      }
    }

    return;
#else
		std::cerr << "Please enable the preprocessor flag ENABLE_FFTW !" << std::endl;
		exit(1);
#endif

	}
}

/** * \brief Smooths the image with a PSF map.
 *
 */
void Observation::CorrelateNoise(PixelMap &pmap)
{
  
  if(sqrt_noise_power.size()==0)return;
  
  if(pmap.getNx() != pmap.getNy()){
    std::cout << "Observation::AddNoise() Doesn't work on nonsquare maps" << std::endl;
    throw std::runtime_error("nonsquare");
  }
  if(pmap.getNx() != Npix_x){
    std::cout << "Observation::AddNoise() Map must have the same dimensions as Observation was constructed with." << std::endl;
    throw std::runtime_error("nonsquare");
  }

#ifdef ENABLE_FFTW
  
  
    // creates plane for fft of map, sets properly input and output data, then performs fft
    assert(Npix_x == Npix_y);
    size_t Npix = Npix_x;
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
/// Outputs rms of noise counts due to background and instrument
/// in the unit decided by the user
float Observation::getBackgroundNoise(float resolution, unitType unit)
{
    if (telescope==true && fabs(resolution-pix_size) > pix_size*1.0e-5)
    {
        std::cout << "The resolution is different from the one of the simulated instrument in Observation::getBackgroundNoise!" << std::endl;
      throw std::runtime_error("The resolution is different from the one of the simulated instrument!");
    }

    double res_in_arcsec = resolution / arcsecTOradians ;
    double back_mean = background_flux * res_in_arcsec*res_in_arcsec * exp_time ;

    double rms = sqrt(exp_num*ron*ron + back_mean) / exp_time;
    if (unit==flux) rms /= e_per_s_to_ergs_s_cm2 * hplanck;
  
    return rms;
}

void Observation::set_up(){
  //zero_point_flux = pow(10,-0.4*mag_zeropoint);  // erg/s/Hz/cm**2
  e_per_s_to_ergs_s_cm2 = pow(10,0.4*(mag_zeropoint-AB_zeropoint));
  background_flux = pow(10,-0.4*(back_mag-mag_zeropoint ));
 }

/// Applies realistic noise (read-out + Poisson) on an image, returns noise map
PixelMap Observation::AddNoise(PixelMap &pmap,Utilities::RandomNumbers_NR &ran)//,long *seed)
{
  if(pmap.getNx() != pmap.getNy()){
    std::cerr << "Observation::AddNoise() Doesn't work on nonsquare maps" << std::endl;
    throw std::runtime_error("nonsquare");
  }
  if(pmap.getUnits() != count_per_sec){
    std::cerr << "Units need to be in counts per second in Observation::AddNoise." << std::endl;
    throw std::runtime_error("wrong units.");
  }
  
  
  PixelMap error_map(pmap);
  
  double res_in_arcsec = pmap.getResolution() / arcsecTOradians;
  double back_mean = background_flux * res_in_arcsec*res_in_arcsec * exp_time ;

  double var_readout = exp_num*ron*ron;
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
//      int k = 0;
//      double p = 1.;
//      double L = exp(-(norm_map + back_mean));
//      do{
//        k++;
//        p *= ran();
//      } while (p > L);
//
      // photons from mean
      int k = ran.poisson(norm_map + back_mean);
      double noise = ran.gauss()*sqrt(var_readout);
      
      noise_map[i] = (k + noise - back_mean - norm_map)/exp_time;
      error_map[i] = (back_mean + norm_map + var_readout)/exp_time;
    }
  }

  CorrelateNoise(noise_map);
  pmap += noise_map;
  
  return error_map;
}

/// Translates photon flux (in 1/(s*cm^2*Hz*hplanck)) into telescope counts per second
void Observation::ToCounts(PixelMap &pmap)
{
  
  //zero_point_flux = pow(10,-0.4*mag_zeropoint);  // erg/s/Hz/cm**2
  //background_flux = pow(10,-0.4*(back_mag-mag_zeropoint ));

  double Q;
  PixelMapUnits units = pmap.getUnits();
  if(units == count_per_sec) return;

  if(units == surfb){
    Q = e_per_s_to_ergs_s_cm2 * hplanck;
    //Q = hplanck/zero_point_flux;
    //Q = pow(10,0.4*(mag_zeropoint - AB_zeropoint))*hplanck;
    //  }else if(pmap.getUnits() == photon_flux ){
    //    Q = hplanck/;
  //}else if(units==photon_flux){
  //  Q = 1/zero_point_flux;
  }else{
    std::cerr << "Map needs to be in photon flux units." << std::endl;
    throw std::runtime_error("wrong units");
  }
  
	pmap.Renormalize(Q);
  pmap.ChangeUnits(count_per_sec);
	return;
}

/// Translates photon flux (in 1/(s*cm^2*Hz*hplanck)) into telescope counts per second
void Observation::ToSurfaceBrightness(PixelMap &pmap)
{

  double Q;
  if(pmap.getUnits() == count_per_sec){
    Q = 1.0/e_per_s_to_ergs_s_cm2 / hplanck;
  }else{
    std::cerr << "Map needs to be in photon flux units." << std::endl;
    throw std::runtime_error("wrong units");
  }
  
  pmap.Renormalize(Q);
  pmap.ChangeUnits(surfb);
  return;
}


