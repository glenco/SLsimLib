/*
 * observation.cpp
 *
 *  Created on: Apr 29, 2013
 *      Author: F. Bellagamba
 */

#include <complex>
#include "slsimlib.h"

#ifdef ENABLE_FITS
#include <CCfits/CCfits>
//#include <CCfits>
#endif

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
  
 
  switch (tel_name) {
    case Euclid_VIS:
      // from Eric
      // Equivalent gain is 11160 e-/ADU
      // Saturation is 71 ADU
      // Equivalent exposure time with 4 frames is 2260s
      // Magnitude Zeropoint is 23.9
      diameter = 119.;
      transmission = 0.30;
      exp_time = 1800.;
      exp_num = 3;
      back_mag = 22.8;
      ron = 5.;
      seeing = 0.18;
      pix_size = .1/60./60./180.*PI;
      break;
   case Euclid_Y:
      diameter = 119.;
      transmission = 0.0961;
      exp_time = 264.;
      exp_num = 3;
      back_mag = 22.57;
      ron = 5.;
      seeing = 0.3;
      pix_size = .3/60./60./180.*PI;
      break;
    case Euclid_J:
      diameter = 119.;
      transmission = 0.0814;
      exp_time = 270.;
      exp_num = 3;
      back_mag = 22.53;
      ron = 5.;
      seeing = 0.3;
      pix_size = .3/60./60./180.*PI;
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
      break;
      
    default:
      throw std::runtime_error("The Telescope selected is not available.");
      break;
	}

	mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) - 48.6;
	telescope = true;
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
diameter(diameter), transmission(transmission), exp_time(exp_time), exp_num(exp_num), back_mag(back_mag), ron(ron), seeing(seeing),Npix_x(Npix_x),Npix_y(Npix_y)
		{
			mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) - 48.6;
			telescope = false;
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
diameter(diameter), transmission(transmission), exp_time(exp_time), exp_num(exp_num), back_mag(back_mag), ron(ron), oversample(oversample),Npix_x(Npix_x),Npix_y(Npix_y)
		{
	mag_zeropoint = 2.5*log10(diameter*diameter*transmission*PI/4./hplanck) - 48.6;

#ifdef ENABLE_FITS

	//std::auto_ptr<CCfits::FITS> fp (new CCfits::FITS (psf_file.c_str(), CCfits::Read));
      
      std::auto_ptr<CCfits::FITS> fp(0);
      try
      {
        fp.reset( new CCfits::FITS (psf_file.c_str(), CCfits::Read) );
      }
      catch (CCfits::FITS::CantOpen)
      {
        std::cerr << "Cannot open " << psf_file << std::endl;
        exit(1);
      }

      
	CCfits::PHDU *h0=&fp->pHDU();
	int side_psf = h0->axis(0);
	int N_psf = side_psf*side_psf;
	map_psf.resize(N_psf);
	h0->read(map_psf);

#else
		std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
		exit(1);
#endif
		telescope = false;
		}

/// Reads in and sets the PSF from a fits file. If the pixel size of the fits is different (smaller) than the one of the telescope, it must be specified.
void Observation::setPSF(std::string psf_file  /// name of fits file with psf
                         , float os  /// over sampling factor
                         )
{
#ifdef ENABLE_FITS
   // std::auto_ptr<CCfits::FITS> fp (new CCfits::FITS (psf_file.c_str(), CCfits::Read));
 
  std::auto_ptr<CCfits::FITS> fp(0);
  try
  {
    fp.reset( new CCfits::FITS (psf_file.c_str(), CCfits::Read) );
  }
  catch (CCfits::FITS::CantOpen)
  {
    std::cerr << "Cannot open " << psf_file << std::endl;
    exit(1);
  }

	CCfits::PHDU *h0=&fp->pHDU();
	int side_psf = h0->axis(0);
	int N_psf = side_psf*side_psf;
	map_psf.resize(N_psf);
	h0->read(map_psf);
    
#else
    std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
    exit(1);
#endif
    
  oversample = os;
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
  
  // calculates normalisation of ncorr
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
  
  sqrt_noise_power.resize(ncorr_big_zeropad_Npixels*(ncorr_big_zeropad_Npixels/2+1));
  
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
}

/**  \brief Converts the input map to a realistic image
 *
 * \param map Input map in photons/(cm^2*Hz)
 * \param psf Decides if the psf smoothing is applied
 * \param noise Decides if noise is added
 * \param unit Decides units of output (if flux, output is in 10**(-0.4*mag)) 
 */
void Observation::Convert(PixelMap &map, bool psf, bool noise, long *seed, unitType unit)
{
  if (telescope == true && fabs(map.getResolution()-pix_size) > pix_size*1.0e-5)
  {
    std::cout << "The resolution of the input map is different from the one of the simulated instrument in Observation::Convert!" << std::endl;
    throw std::runtime_error("The resolution of the input map is different from the one of the simulated instrument!");
  }
  //PixelMap outmap =
  PhotonToCounts(map);
  if (psf == true)  ApplyPSF(map);
  if (noise == true) AddNoise(map,seed);
  
  if (unit == flux)
  {
    double counts_to_flux = pow(10,-0.4*mag_zeropoint);
    map.Renormalize(counts_to_flux);
  }
  return;
}

/// Converts an observed image to the units of the lensing simulation
void Observation::Convert_back (PixelMap &map)
{
	//PixelMap outmap(map);
	double Q = pow(10,0.4*(mag_zeropoint+48.6))*hplanck;
	map.Renormalize(1./Q);
	return;
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
//#ifdef ENABLE_FITS
#ifdef ENABLE_FFTW
    
    // calculates normalisation of psf
    int N_psf = map_psf.size();
    int side_psf = sqrt(N_psf);
    double map_norm = 0.;
    for (int i = 0; i < N_psf; i++)
    {
      map_norm += map_psf[i];
    }
    fftw_plan p_psf;
    
    // creates plane for fft of map, sets properly input and output data, then performs fft
    fftw_plan p;
    long Npix = pmap.getNx();
    long Npix_zeropad = Npix + side_psf;
        
    // rows and columns between first_p and last_p are copied in the zero-padded version
    long first_p = side_psf/2;
    long last_p = first_p + (Npix-1);
    std::vector<std::complex<double> > out(Npix_zeropad*(Npix_zeropad/2+1));
    
    // add zero-padding
    std::vector<double> in_zeropad(Npix_zeropad*Npix_zeropad);
    for (int i = 0; i < Npix_zeropad*Npix_zeropad; i++)
    {
      long ix = i/Npix_zeropad;
      long iy = i%Npix_zeropad;
      if (ix >= first_p && ix <= last_p && iy >= first_p && iy <= last_p)
        in_zeropad[i] = pmap[(ix-side_psf/2)*Npix+(iy-side_psf/2)];
      else
        in_zeropad[i] = 0.;
    }
    
    p = fftw_plan_dft_r2c_2d(Npix_zeropad,Npix_zeropad,in_zeropad.data()
                             , reinterpret_cast<fftw_complex*>(out.data()), FFTW_ESTIMATE);
    fftw_execute(p);
    
    // creates plane for perform backward fft after convolution, sets output data 
    fftw_plan p2;
    //double* out2 = new double[Npix_zeropad*Npix_zeropad];
    std::vector<double> out2(Npix_zeropad*Npix_zeropad);
    p2 = fftw_plan_dft_c2r_2d(Npix_zeropad,Npix_zeropad,reinterpret_cast<fftw_complex*>(out.data()), out2.data(), FFTW_ESTIMATE);
    
    
    // arrange psf data for fft, creates plane, then performs fft
    // psf data are moved into the four corners of psf_big_zeropad
    long psf_big_zeropad_Npixels = static_cast<int>((Npix+side_psf)*oversample);
    std::vector<double> psf_big_zeropad(psf_big_zeropad_Npixels*psf_big_zeropad_Npixels);
    std::vector<std::complex<double> > out_psf(psf_big_zeropad_Npixels*(psf_big_zeropad_Npixels/2+1));
    
    p_psf = fftw_plan_dft_r2c_2d(psf_big_zeropad_Npixels,psf_big_zeropad_Npixels,psf_big_zeropad.data(), reinterpret_cast<fftw_complex*>(out_psf.data()), FFTW_ESTIMATE);
    long ix, iy;
    for (int i = 0; i < psf_big_zeropad_Npixels*psf_big_zeropad_Npixels; i++)
    {
      ix = i/psf_big_zeropad_Npixels;
      iy = i%psf_big_zeropad_Npixels;
      if(ix<side_psf/2 && iy<side_psf/2)
        psf_big_zeropad[i] = map_psf[(ix+side_psf/2)*side_psf+(iy+side_psf/2)]/map_norm;
      else if(ix<side_psf/2 && iy>=psf_big_zeropad_Npixels-side_psf/2)
        psf_big_zeropad[i] = map_psf[(ix+side_psf/2)*side_psf+(iy-(psf_big_zeropad_Npixels-side_psf/2))]/map_norm;
      else if(ix>=psf_big_zeropad_Npixels-side_psf/2 && iy<side_psf/2)
        psf_big_zeropad[i] = map_psf[(ix-(psf_big_zeropad_Npixels-side_psf/2))*side_psf+(iy+side_psf/2)]/map_norm;
      else if(ix>=psf_big_zeropad_Npixels-side_psf/2 && iy>=psf_big_zeropad_Npixels-side_psf/2)
        psf_big_zeropad[i] = map_psf[(ix-(psf_big_zeropad_Npixels-side_psf/2))*side_psf+(iy-(psf_big_zeropad_Npixels-side_psf/2))]/map_norm;
      else
        psf_big_zeropad[i] = 0.;
    }
    fftw_execute(p_psf);
    
    // performs convolution in Fourier space , and transforms back to real space
    for (unsigned long i = 0; i < Npix_zeropad*(Npix_zeropad/2+1); i++)
    {
      ix = i/(Npix_zeropad/2+1);
      iy = i%(Npix_zeropad/2+1);
      if (ix>Npix_zeropad/2)
        out[i] *= out_psf[(psf_big_zeropad_Npixels-(Npix_zeropad-ix))*(psf_big_zeropad_Npixels/2+1)+iy];
      else
        out[i] *= out_psf[ix*(psf_big_zeropad_Npixels/2+1)+iy];
    }
    fftw_execute(p2);
    
    // translates array of data in (normalised) counts map
    for (unsigned long i = 0; i < Npix_zeropad*Npix_zeropad; i++)
    {
      ix = i/Npix_zeropad;
      iy = i%Npix_zeropad;
      if (ix >= first_p && ix <= last_p && iy >= first_p && iy <= last_p)
      {
        int ii = (ix-side_psf/2)*Npix+(iy-side_psf/2);
        pmap.AssignValue(ii,out2[i]/double(Npix_zeropad*Npix_zeropad));
      }
    }
    return;
#else
		std::cout << "Please enable the preprocessor flag ENABLE_FFTW !" << std::endl;
		exit(1);
#endif
    
//#else
//		std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
//		exit(1);
//#endif
	}
}

/** * \brief Smooths the image with a PSF map.
 *
 */
void Observation::CorrelateNoise(PixelMap &pmap)
{
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
    std::vector<std::complex<double> > out(fftsize);
    
    // add zero-padding
    std::vector<double> in_zeropad(Npix_zeropad*Npix_zeropad);
    for (int i = 0; i < Npix_zeropad*Npix_zeropad; i++)
    {
      long ix = i/Npix_zeropad;
      long iy = i%Npix_zeropad;
      if (ix >= first_p && ix <= last_p && iy >= first_p && iy <= last_p)
        in_zeropad[i] = pmap[(ix-side_ncorr/2)*Npix+(iy-side_ncorr/2)];
      else
        in_zeropad[i] = 0.;
    }
    
    fftw_plan p = fftw_plan_dft_r2c_2d(Npix_zeropad,Npix_zeropad,in_zeropad.data()
                             , reinterpret_cast<fftw_complex*>(out.data()), FFTW_ESTIMATE);
    fftw_execute(p);
    
    // creates plane for perform backward fft after convolution, sets output data
    fftw_plan p2;
    //double* out2 = new double[Npix_zeropad*Npix_zeropad];
    std::vector<double> out2(Npix_zeropad*Npix_zeropad);
    p2 = fftw_plan_dft_c2r_2d(Npix_zeropad,Npix_zeropad,reinterpret_cast<fftw_complex*>(out.data()), out2.data(), FFTW_ESTIMATE);
    
  
    // performs convolution in Fourier space , and transforms back to real space
    for (unsigned long i = 0; i < fftsize ; i++)
    {
      size_t ix = i/(Npix_zeropad/2+1);
      size_t iy = i%(Npix_zeropad/2+1);
      if (ix>Npix_zeropad/2)
        out[i] *= sqrt_noise_power[(ncorr_big_zeropad_Npixels-(Npix_zeropad-ix))*(ncorr_big_zeropad_Npixels/2+1)+iy];
      else
        out[i] *= sqrt_noise_power[ix*(ncorr_big_zeropad_Npixels/2+1)+iy];
    }
    fftw_execute(p2);
    
    // translates array of data in (normalised) counts map
    for (unsigned long i = 0; i < Npix_zeropad*Npix_zeropad; i++)
    {
      size_t ix = i/Npix_zeropad;
      size_t iy = i%Npix_zeropad;
      if (ix >= first_p && ix <= last_p && iy >= first_p && iy <= last_p)
      {
        int ii = (ix-side_ncorr/2)*Npix+(iy-side_ncorr/2);
        pmap[ii] = out2[i]/(Npix_zeropad*Npix_zeropad);
      }
    }
    return;
#else
    std::cout << "Please enable the preprocessor flag ENABLE_FFTW !" << std::endl;
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

    double Q = pow(10,0.4*(mag_zeropoint+48.6));
    double res_in_arcsec = resolution*180.*60.*60/PI;
    double back_mean = pow(10,-0.4*(48.6+back_mag))*res_in_arcsec*res_in_arcsec*Q*exp_time;
    
    double rms = sqrt(exp_num*ron*ron+back_mean)/exp_time;
    
    if (unit==flux) rms *= pow(10,-0.4*mag_zeropoint);
    
    return rms;
}

/// Applies realistic noise (read-out + Poisson) on an image
void Observation::AddNoise(PixelMap &pmap,long *seed)
{
  if(pmap.getNx() != pmap.getNy()){
    std::cout << "Observation::AddNoise() Doesn't work on nonsquare maps" << std::endl;
    throw std::runtime_error("nonsquare");
  }
  
  //PixelMap outmap(pmap);
  double Q = pow(10,0.4*(mag_zeropoint+48.6));
  double res_in_arcsec = pmap.getResolution()*180.*60.*60/PI;
  double back_mean = pow(10,-0.4*(48.6+back_mag))*res_in_arcsec*res_in_arcsec*Q*exp_time;
  double rms, noise;
  double rms2 = sqrt(exp_num*ron*ron);
  double norm_map;
  
  PixelMap noise_map(pmap);
  //double sum=0,sum2=0;
  size_t N = pmap.size();
  for (unsigned long i = 0; i < N ; i++)
  {
    norm_map = pmap[i]*exp_time;
    if (norm_map+back_mean > 500.)
    {
      rms = sqrt(exp_num*ron*ron + norm_map + back_mean);
      noise = gasdev(seed)*rms;
      noise_map[i] = noise/exp_time;
    }
    else
    {
      int k = 0;
      double p = 1.;
      double L = exp(-(norm_map + back_mean));
      while (p > L)
      {
        k++;
        p *= ran2(seed);
      }
      noise = gasdev(seed)*rms2;
      noise_map[i] = (k - 1 + noise - back_mean - norm_map)/exp_time;
    }
    //sum += noise_map[i];
    //sum2 += noise_map[i]*noise_map[i];
  }
  
  //sum /= N;
  //sum2 = (sum2 - N*sum*sum)/(N-1);
  //std::cout << "Noise Variance in flux units is : " << sum2 << std::endl;
  
  CorrelateNoise(noise_map);
  
  pmap += noise_map;
  
  return;
}

/// Translates photon flux (in 1/(s*cm^2*Hz*hplanck)) into telescope pixel counts
void Observation::PhotonToCounts(PixelMap &pmap)
{
	//PixelMap outmap(pmap);
	double Q = pow(10,0.4*(mag_zeropoint+48.6))*hplanck;
	pmap.Renormalize(Q);
	return;
}

