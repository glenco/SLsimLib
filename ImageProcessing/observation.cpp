/*
 * observation.cpp
 *
 *  Created on: Apr 29, 2013
 *      Author: F. Bellagamba
 */

#include "slsimlib.h"

#ifdef ENABLE_FITS
#include <CCfits/CCfits>
//#include <CCfits>
#endif

#ifdef ENABLE_FFTW
#include "fftw3.h"
#endif

#include <fstream>

/*\brief Creates an observation setup that mimics a known instrument
 *
 */
Observation::Observation(Telescope tel_name)
{
	if (tel_name == Euclid_VIS)
	{
		diameter = 119.;
		transmission = 0.30;
		exp_time = 1800.;
		exp_num = 3;
		back_mag = 22.8;
		ron = 5.;
		seeing = 0.18;
	}
	mag_zeropoint = 2.5*log10(diameter*diameter*transmission*pi/4./hplanck) - 48.6;
}

/* Creates a custom observation setup with parameters decided by the user.
 *
 * \param diameter Diameter of telescope (in cm) (Collecting area is pi/4.*diameter^2)
 * \param transmission Total transmission of the telescope (/int T(\lambda)/\lambda d\lambda)
 * \param exp_time Total exposure time
 * \param exp_num Number of exposures
 * \param back_mag Flux due to background and/or sky in mag/arcsec^2
 * \param ron Read-out noise in electrons/pixel
 * \param seeing FWHM in arcsecs of the image
 */
Observation::Observation(float diameter, float transmission, float exp_time, int exp_num, float back_mag, float ron, float seeing):
		exp_time(exp_time), exp_num(exp_num), back_mag(back_mag), diameter(diameter), transmission(transmission), ron(ron), seeing(seeing)
		{
			mag_zeropoint = 2.5*log10(diameter*diameter*transmission*pi/4./hplanck) - 48.6;
		}

/* Creates a custom observation setup with parameters decided by the user. Allows for the use of a psf fits image.
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
Observation::Observation(float diameter, float transmission, float exp_time, int exp_num, float back_mag, float ron, std::string psf_file, float oversample):
		exp_time(exp_time), exp_num(exp_num), back_mag(back_mag), diameter(diameter), transmission(transmission), ron(ron), oversample(oversample)
		{
	mag_zeropoint = 2.5*log10(diameter*diameter*transmission*pi/4./hplanck) - 48.6;

#ifdef ENABLE_FITS

	std::auto_ptr<CCfits::FITS> fp (new CCfits::FITS (psf_file.c_str(), CCfits::Read));
	CCfits::PHDU *h0=&fp->pHDU();
	int side_psf = h0->axis(0);
	int N_psf = side_psf*side_psf;
	map_psf.resize(N_psf);
	h0->read(map_psf);

#else
		std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
		exit(1);
#endif
		}

/* \brief Converts the input map to a realistic image
 *
 * \param map Input map in photons/(cm^2*Hz)
 * \param psf Decides if the psf smoothing is applied
 * \param noise Decides if noise is added
 */
PixelMap Observation::Convert (PixelMap &map, bool psf, bool noise)
{
	PixelMap outmap = PhotonToCounts(map);
	if (psf == true)  outmap = ApplyPSF(outmap);
	if (noise == true) outmap = AddNoise(outmap);
	return outmap;
}

/// Converts an observed image to the units of the lensing simulation
PixelMap Observation::Convert_back (PixelMap &map)
{
	PixelMap outmap(map);
	double Q = pow(10,0.4*(mag_zeropoint+48.6))*hplanck;
	outmap.Renormalize(1./Q);
	return outmap;
}


/** \brief Smooths the image with a PSF map.
*
*/
PixelMap Observation::ApplyPSF(PixelMap &pmap)
{
	if (map_psf.size() == 0)
	{
		if (seeing > 0.)
		{
			PixelMap outmap(pmap);
			outmap.smooth(seeing/2.355);
			return outmap;
		}
		else
		{
			return pmap;
		}
	}
	else
	{
#ifdef ENABLE_FITS
#ifdef ENABLE_FFTW

	PixelMap outmap(pmap);
	// creates plane for fft of map, sets properly input and output data, then performs fft
	fftw_plan p;
	long Npix = outmap.getNpixels();
	double* in = new double[Npix*Npix];
	std::complex<double>* out=new std::complex<double> [Npix*(Npix/2+1)];
	for (unsigned long i = 0; i < Npix*Npix; i++)
	{
		in[i] = outmap[i];
	}
	p = fftw_plan_dft_r2c_2d(Npix,Npix,in, reinterpret_cast<fftw_complex*>(out), FFTW_ESTIMATE);
	fftw_execute(p);

	// creates plane for perform backward fft after convolution, sets output data
	fftw_plan p2;
	double* out2 = new double[Npix*Npix];
	p2 = fftw_plan_dft_c2r_2d(Npix,Npix,reinterpret_cast<fftw_complex*>(out), out2, FFTW_ESTIMATE);

	// calculates normalisation of psf
	int N_psf = map_psf.size();
	int side_psf = sqrt(N_psf);
	double map_norm = 0.;
	for (int i = 0; i < N_psf; i++)
	{
		map_norm += map_psf[i];
	}
	fftw_plan p_psf;

	// arrange psf data for fft, creates plane, then performs fft
	int psf_big_Npixels = static_cast<int>(Npix*oversample);
	double* psf_big = new double[psf_big_Npixels*psf_big_Npixels];
	std::complex<double>* out_psf=new std::complex<double> [psf_big_Npixels*(psf_big_Npixels/2+1)];
	p_psf = fftw_plan_dft_r2c_2d(psf_big_Npixels,psf_big_Npixels,psf_big, reinterpret_cast<fftw_complex*>(out_psf), FFTW_ESTIMATE);
	long ix, iy;
	for (int i = 0; i < psf_big_Npixels*psf_big_Npixels; i++)
	{
		ix = i/psf_big_Npixels;
		iy = i%psf_big_Npixels;
		if(ix<side_psf/2 && iy<side_psf/2)
			psf_big[i] = map_psf[(ix+side_psf/2)*side_psf+(iy+side_psf/2)]/map_norm;
		else if(ix<side_psf/2 && iy>=psf_big_Npixels-side_psf/2)
			psf_big[i] = map_psf[(ix+side_psf/2)*side_psf+(iy-(psf_big_Npixels-side_psf/2))]/map_norm;
		else if(ix>=psf_big_Npixels-side_psf/2 && iy<side_psf/2)
			psf_big[i] = map_psf[(ix-(psf_big_Npixels-side_psf/2))*side_psf+(iy+side_psf/2)]/map_norm;
		else if(ix>=psf_big_Npixels-side_psf/2 && iy>=psf_big_Npixels-side_psf/2)
			psf_big[i] = map_psf[(ix-(psf_big_Npixels-side_psf/2))*side_psf+(iy-(psf_big_Npixels-side_psf/2))]/map_norm;
		else
			psf_big[i] = 0.;
	}
	fftw_execute(p_psf);

	// performs convolution in Fourier space, and transforms back to real space
	for (unsigned long i = 0; i < Npix*(Npix/2+1); i++)
	{
		ix = i/(Npix/2+1);
		iy = i%(Npix/2+1);
		if (ix>Npix/2)
			out[i] *= out_psf[(psf_big_Npixels-(Npix-ix))*(psf_big_Npixels/2+1)+iy];
		else
			out[i] *= out_psf[ix*(psf_big_Npixels/2+1)+iy];
	}
	fftw_execute(p2);

	// translates array of data in (normalised) counts map
	for (unsigned long i = 0; i < Npix*Npix; i++)
	{
		ix = i/Npix;
		iy = i%Npix;
		outmap.AssignValue(i,out2[i]/double(Npix*Npix));
	}
	return outmap;
#else
		std::cout << "Please enable the preprocessor flag ENABLE_FFTW !" << std::endl;
		exit(1);
#endif

#else
		std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
		exit(1);
#endif
	}
}

/// Applies realistic noise (read-out + Poisson) on an image
PixelMap Observation::AddNoise(PixelMap &pmap)
{
	PixelMap outmap(pmap);
	double Q = pow(10,0.4*(mag_zeropoint+48.6));
	double res_in_arcsec = outmap.getResolution()*180.*60.*60/pi;
	double back_mean = pow(10,-0.4*(48.6+back_mag))*res_in_arcsec*res_in_arcsec*Q*exp_time;
	double rms, noise;
	long seed = 24;
	double norm_map;
	for (unsigned long i = 0; i < outmap.getNpixels()*outmap.getNpixels(); i++)
	{
		norm_map = outmap[i]*exp_time;
		if (norm_map+back_mean > 500.)
		{
			rms = sqrt(pow(exp_num*ron,2)+norm_map+back_mean);
			noise = gasdev(&seed)*rms;
			outmap.AssignValue(i,double(norm_map+noise)/exp_time);
		}
		else
		{
			int k = 0;
			double p = 1.;
			double L = exp(-(norm_map+back_mean));
			while	(p > L)
			{
				k++;
				p *= ran2(&seed);
			}
			outmap.AssignValue(i,double(k-1-back_mean)/exp_time);
		}
	}
	return outmap;
}

/// Translates photon flux (in photons/(cm^2*Hz)) into telescope pixel counts
PixelMap Observation::PhotonToCounts(PixelMap &pmap)
{
	PixelMap outmap(pmap);
	double Q = pow(10,0.4*(mag_zeropoint+48.6))*hplanck;
	outmap.Renormalize(Q);
	return outmap;
}

