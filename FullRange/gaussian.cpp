//
//  guassian.cpp
//  NR
//
//  Created by Ben Metcalf on 17/04/2018.
//

#include "gaussian.h"
#include "slsimlib.h"
#include <complex>
#include <fstream>


#ifdef ENABLE_FFTW
#include "fftw3.h"
#endif

/** \brief Index of the discrete Fourier mode using FFTW convention
 */
double k_to_index_2d(double k){
  
}
/** \brief Discete Fourier k from 1d index using FFTW convention
 
 k is in units 
 */
double index_to_k_2d(size_t i){
  
}

template <typename T>
void gaussianField(std::vector<T> &im,std::vector<T> &power
        ,size_t Nx,size_t Ny,Utilities::RandomNumbers_NR &ran){
  
  size_t N = Nx*Ny;
  if(N != im.size()){
    std::cerr << "incompatable input in gaussianField" << std::endl;
    throw std::invalid_argument("");
  }
  std::vector<std::complex<T> > fftmap(Nx*(Ny/2+1));
 
  for(size_t i = 0 ; i < N ; ++i){
    fftmap[i].real(sqrt(power[i])*ran.gauss());
    fftmap[i].imag(sqrt(power[i])*ran.gauss());
  }
  
  fftw_plan p = fftw_plan_dft_c2r_2d(Nx,Ny,im.data()
                           , reinterpret_cast<fftw_complex*>(fftmap.data()), FFTW_ESTIMATE);
  fftw_execute(p);
  
  fftw_destroy_plan(p);
}
