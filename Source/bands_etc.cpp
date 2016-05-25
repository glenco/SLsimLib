//
//  bands_etc.cpp
//  GLAMER
//
//  Created by bmetcalf on 13/04/16.
//
//

#include <stdio.h>
#include "InputParams.h"
#include "bands_etc.h"

double BandInfo::lambda_central(Band band){   /// central wavelength in Angstroms
  
  int i=0;
  while(bands[i] != band && i < bands.size() ) ++i;
  
  return lambda[i];
}


/// absolute magnitude of the sun
double SunInfo::AbsMag(Band band){
  
  double wavelength = bandinfo.lambda_central(band);
  
  auto it = std::lower_bound(wavelengths.begin(),wavelengths.end(),wavelength);
  
  return -2.5 * std::log10(sed[it - wavelengths.begin()]) - 48.6;
}

/// k-corrected apparent magnitude, very crudely done without integration of spectrum
double SunInfo::mag(Band band,double z,const COSMOLOGY &cosmo){
  
  double wavelength = bandinfo.lambda_central(band);
  double lambda = wavelength/(1+z);
  
  auto it = std::lower_bound(wavelengths.begin(),wavelengths.end(),lambda);
  
  // remember that cosmo.lumDist(z) is the bolometric distance
  return -2.5 * std::log10(sed[it - wavelengths.begin()]*(1+z)) - 48.6 + 5 * (std::log10(cosmo.lumDist(z)*1.0e6) - 1);
}

/// k-corrected flux in ergs/s/cm^2/Hz
double SunInfo::flux(Band band,double z,const COSMOLOGY &cosmo){
  
  double wavelength = bandinfo.lambda_central(band);
  double lambda = wavelength/(1+z);
  
  auto it = std::lower_bound(wavelengths.begin(),wavelengths.end(),lambda);
  
  // remember that cosmo.lumDist(z) is the bolometric distance
  return (1+z)*sed[it - wavelengths.begin()]*pow(1.0e5/cosmo.lumDist(z),2);
}
