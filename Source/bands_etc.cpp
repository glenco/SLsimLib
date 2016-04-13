//
//  bands_etc.cpp
//  GLAMER
//
//  Created by bmetcalf on 13/04/16.
//
//

#include <stdio.h>
#include "bands_etc.h"

double BandInfo::lambdal(Band band){   /// central wavelength in Angstroms
  
  int i=0;
  while(bands[i] != band && i < bands.size() ) ++i;
  
  return lambda[i];
}


/// absolute magnitude of the sun
double SunInfo::AbsMag(Band band){
  
  double wavelength = bandinfo.lambda(band);
  
  int i = std::lower_bound(wavelengths.begin(),wavelengths.end(),wavelength);
  
  return -2.5 * std::log10(sed[i]) - 48.6;
}

/// k-corrected apparent magnitude
double SunInfo::mag(Band band,double z,const COSMOLOGY &como){
  
  double wavelength = bandinfo.central(band);
  double lambda = wavelength/(1+z);
  
  int i = std::lower_bound(wavelengths.begin(),wavelengths.end(),wavelength);
  
  return -2.5 * std::log10(sed[i]) - 48.6 + 5 * (std::log10(cosmo.lumDist(z)*1.0e6) - 1);
}
