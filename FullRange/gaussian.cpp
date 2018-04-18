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
#include <functional>



#ifdef ENABLE_FFTW
#include "fftw3.h"
#endif

void GaussianField::GaussianField2D(double *image
                                    ,std::function<double(double,double)> PofK
                                    ,double boxlrad
                                    ,Utilities::RandomNumbers_NR &ran
                                    ){
  
  size_t i,j,ii;
  double kx,ky;
  double p;
  
  double K = 2. * M_PI  / boxlrad;
  
  for (i=0; i < Nx; i++){
    for (j=0; j < Ny/2+1; j++){
      
      ii = i*(Ny/2+1) + j;
      
      kx= (i<Nx/2) ? double(i)*K : double(i-Nx)*K;
      ky= (j<Ny/2) ? double(j)*K : double(j-Ny)*K;

      if((kx == 0 && ky == 0) || ii == 0){
        kappak[ii].real(0);
        kappak[ii].imag(0);
      }else{

        p = sqrt(PofK(kx,ky)/2);
        
        kappak[ii].imag(ran.gauss() * p * K);
        kappak[ii].real(ran.gauss() * p * K);
      }
      //std::cout << kappak[ii] << std::endl;
    }
  }
  
  fftw_execute_dft_c2r(pbackward,reinterpret_cast<fftw_complex*>(kappak.data()),image);
  return;
}

void GaussianField::GaussianField2D(std::vector<double> &image
                                    ,std::function<double(double,double)> PofK
                                    ,double boxlrad
                                    ,Utilities::RandomNumbers_NR &ran
                                    ){
  image.resize(Nx*Ny);
  GaussianField2D(image.data(),PofK,boxlrad,ran);
}

void GaussianField::GaussianField2D(std::valarray<double> &image
                                    ,std::function<double(double,double)> PofK
                                    ,double boxlrad
                                    ,Utilities::RandomNumbers_NR &ran
                                    ){
  image.resize(Nx*Ny);
  GaussianField2D(&(image[0]),PofK,boxlrad,ran);
}
