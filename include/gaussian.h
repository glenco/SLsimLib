//
//  guassian.hpp
//  NR
//
//  Created by Ben Metcalf on 17/04/2018.
//

#ifndef guassian_hpp
#define guassian_hpp

#include <stdio.h>
#include <complex>
#include <fstream>
#include <vector>
#include "slsimlib.h"
#include <functional>


#ifdef ENABLE_FFTW
#include "fftw3.h"
#endif

class GaussianField{
public:
  GaussianField(size_t Nx,size_t Ny):
  Nx(Nx),Ny(Ny)
  {
    kappak.resize(Nx*(Ny/2+1));
    std::vector<double> image(Nx*Ny);
    
    pbackward = fftw_plan_dft_c2r_2d(Nx,Ny
                              ,reinterpret_cast<fftw_complex*>(kappak.data())
                              ,image.data(), FFTW_ESTIMATE);
  }

  ~GaussianField(){
    fftw_destroy_plan(pbackward);
  }
  void GaussianField2D(double *image,std::function<double(double,double)> PofK,
                       double boxlrad,Utilities::RandomNumbers_NR &ran);
  
  void GaussianField2D(std::vector<double> &image,std::function<double(double,double)> PofK,
                       double boxlrad,Utilities::RandomNumbers_NR &ran);
  
  void GaussianField2D(std::valarray<double> &image,std::function<double(double,double)> PofK,
                       double boxlrad,Utilities::RandomNumbers_NR &ran);

  
private:
  std::vector<std::complex<double> > kappak;
  fftw_plan pbackward;
  size_t Nx;
  size_t Ny;
};

/** \brief Discete Fourier k from 1d index using FFTW convention
 asd
 k is in units (box length)^-1
 */
void index_to_k_2d(size_t k  /// id index
                   ,double &kx  /// output kx in 1/Length_x units
                   ,double &ky  /// output ky in 1/Length_y units
                   ,size_t Nx   /// number of pixels on x-side
                   ,size_t Ny   /// number of pixels on y-side
                   );

void GaussianField2D(double *image,std::function<double(double,double)> PofK,
                     double boxlrad,int Nx,int Ny,Utilities::RandomNumbers_NR &ran);


#endif /* guassian_hpp */
