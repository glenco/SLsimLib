/*
 * profile.cpp
 *
 *  Created on: May 23, 2012
 *      Author: cgiocoli
 */

#include "profile.h"

// create profile of the maps for each lensing component - spherical simmetry is assumed
double * estprof(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, double dr0, double xmax){
  int nbin = int(xmax/dr0);
  //std:: cout << " nbins (in estprof) = " << nbin << std:: endl;
  double *kr = new double[nbin];
  for (int k=0;k<nbin;k++){
    int contapx=0;
    kr[k] = 0;
    // for each bin in r estimate the mean value
    for( int i=0; i<nx; i++ )
      for( int j=0; j<ny; j++ ){
        if(r[i+ny*j]>dr0*double(k) && r[i+ny*j]<=dr0*double(k+1)){
          contapx = contapx + 1;
          kr[k] = kr[k] + q[i+ny*j];
        }
      }
    kr[k] = kr[k]/double(contapx);
    if(contapx==0) kr[k]=0.;
  }
  return kr; // return the pointer
}

// variance of the profile
double * estsigmaprof(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, double dr0, double xmax, \
 double* qm){
  int nbin = int(xmax/dr0);
  double *kr = new double[nbin];
  for (int k=0;k<nbin;k++){
    int contapx=0;
    kr[k] = 0;
    // for each bin in r estimate the mean value
    for( int i=0; i<nx; i++ )
      for( int j=0; j<ny; j++ ){
        if(r[i+ny*j]>dr0*double(k) && r[i+ny*j]<=dr0*double(k+1)){
          contapx = contapx + 1;
          kr[k] = (q[i+ny*j]-qm[k])*(q[i+ny*j]-qm[k]) + kr[k];
          if(kr[k]<0) {
            std:: cout << "negative " << kr[k] << std:: endl;
            exit(1);
          }
        }
      }
    kr[k] = sqrt(kr[k]/double(contapx));
    if(contapx==0) kr[k]=0.;
  }
  return kr; // return the pointer
}

// create cumulative profile of the maps for each lensing component - spherical simmetry is assumed
double * estcprof(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, double dr0, double xmax){
  int nbin = int(xmax/dr0);
  //std:: cout << " nbins (in estprof) = " << nbin << std:: endl;
  double *kr = new double[nbin];
  for (int k=0;k<nbin;k++){
    double ranulus = dr0*(double(k)+double(k+1))*0.5;
    int contapx=0;
    kr[k] = 0;
    // for each bin in r estimate the mean value
    for( int i=0; i<nx; i++ )
      for( int j=0; j<ny; j++ ){
        if(r[i+ny*j]<=ranulus){
          contapx = contapx + 1;
          kr[k] = kr[k] + q[i+ny*j];
        }
      }
    kr[k] = kr[k]/double(contapx);
    if(contapx==0) kr[k]=0.;
  }
  return kr; // return the pointer
}

// variance of the profile
double * estsigmacprof(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, double dr0, double xmax, \
 double* qm){
  int nbin = int(xmax/dr0);
  double *kr = new double[nbin];
  for (int k=0;k<nbin;k++){
    double ranulus = dr0*(double(k)+double(k+1))*0.5;
    int contapx=0;
    kr[k] = 0;
    // for each bin in r estimate the mean value
    for( int i=0; i<nx; i++ )
      for( int j=0; j<ny; j++ ){
        if(r[i+ny*j]<=ranulus){
          contapx = contapx + 1;
          kr[k] = (q[i+ny*j]-qm[k])*(q[i+ny*j]-qm[k]) + kr[k];
          if(kr[k]<0) {
            std:: cout << "negative " << kr[k] << std:: endl;
            exit(1);
          }
        }
      }
    kr[k] = sqrt(kr[k]/double(contapx));
    if(contapx==0) kr[k]=0.;
  }
  return kr; // return the pointer
}