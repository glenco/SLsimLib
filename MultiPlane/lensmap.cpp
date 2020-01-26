//
//  lensmap.c
//  GLAMER
//
//  Created by bmetcalf on 18.12.19.
//

#include <stdio.h>
#include "fftw3.h"
#include "lensmap.h"

bool LensMap::evaluate(const double *x,float &sigma,float *gamma,double *alpha) {
  
  double fx = (x[0] - lowerleft[0]) / x_resolution();
  double fy = (x[1] - lowerleft[1]) / y_resolution();
  
  if( (fx>=0)*(fx<nx)*(fy>=0)*(fy<ny) ){
    
    size_t ix = fx;
    if(ix == nx-1) ix = nx-2;
    size_t iy = fy;
    if(iy == ny-1) iy = ny-2;
    
    size_t index = ix + nx * iy;
    
    fx = fx - ix;
    fy = fy - iy;
    
    double a = (1-fx)*(1-fy);
    double b = fx*(1-fy);
    double c = fx*fy;
    double d = (1-fx)*fy;
    
    // bilinear interpolation
    sigma = a * surface_density[index] + b * surface_density[index+1]
    + c * surface_density[index+1+nx] + d * surface_density[index+nx];
    
    alpha[0] = a * alpha1_bar[index] + b * alpha1_bar[index+1]
    + c * alpha1_bar[index+1+nx] + d * alpha1_bar[index+nx];
    alpha[1] = a * alpha2_bar[index] + b * alpha2_bar[index+1]
    + c * alpha2_bar[index+1+nx] + d * alpha2_bar[index+nx];
    
    gamma[0] = a * gamma1_bar[index] + b * gamma1_bar[index+1]
    + c * gamma1_bar[index+1+nx] + d * gamma1_bar[index+nx];
    gamma[1] = a * gamma2_bar[index] + b * gamma2_bar[index+1]
    + c * gamma2_bar[index+1+nx] + d * gamma2_bar[index+nx];
    gamma[2] = 0.0;
    
    /*
     if(isnan(gamma[1])){
     std::cerr << index+1+nx << " < " << gamma2_bar.size() << std::endl;
     std::cerr << alpha[0] << " " << alpha[1] << std::endl;
     std::cerr << gamma[0] << " " << gamma[1] << " " << gamma[2] << std::endl;
     
     std::cerr << gamma2_bar[index] << " " << gamma2_bar[index+1]
     << " " << gamma2_bar[index+1+nx] << " " << gamma2_bar[index+nx] << std::endl;
     
     
     assert(!isnan(gamma[1]));
     }*/
    return false;
  }
  
  sigma = 0;
  alpha[0] = alpha[1] = 0.0;
  gamma[0] = gamma[1] = gamma[2] = 0.0;
  
  return true;
}


LensMap::LensMap(LensMap &&m){
  surface_density=std::move(m.surface_density);
  alpha1_bar =std::move(m.alpha1_bar);
  alpha2_bar =std::move(m.alpha2_bar);
  gamma1_bar =std::move(m.gamma1_bar);
  gamma2_bar =std::move(m.gamma2_bar);
  phi_bar = std::move(m.phi_bar);
  
  nx = m.nx;
  ny = m.ny;
  boxlMpc = m.boxlMpc;
  center = m.center;
  lowerleft = m.lowerleft;
  upperright = m.upperright;
  angular_pixel_size = m.angular_pixel_size;
}

LensMap& LensMap::operator=(LensMap &&m){
  if(&m==this) return *this;
  
  surface_density=std::move(m.surface_density);
  alpha1_bar =std::move(m.alpha1_bar);
  alpha2_bar =std::move(m.alpha2_bar);
  gamma1_bar =std::move(m.gamma1_bar);
  gamma2_bar =std::move(m.gamma2_bar);
  phi_bar = std::move(m.phi_bar);
  
  nx = m.nx;
  ny = m.ny;
  boxlMpc = m.boxlMpc;
  center = m.center;
  lowerleft = m.lowerleft;
  upperright = m.upperright;
  angular_pixel_size = m.angular_pixel_size;
  
  return *this;
}
LensMap& LensMap::operator=(const LensMap &m){
  if(&m==this) return *this;
  
  surface_density=m.surface_density;
  alpha1_bar =m.alpha1_bar;
  alpha2_bar =m.alpha2_bar;
  gamma1_bar =m.gamma1_bar;
  gamma2_bar =m.gamma2_bar;
  phi_bar = m.phi_bar;
  
  nx = m.nx;
  ny = m.ny;
  boxlMpc = m.boxlMpc;
  center = m.center;
  lowerleft = m.lowerleft;
  upperright = m.upperright;
  angular_pixel_size = m.angular_pixel_size;
  
  return *this;
}

void LensMap::read_header(std::string fits_input_file
                          ,double angDist){
  
  CPFITS_READ cpfits(fits_input_file);
  
  std::vector<long> size;
  //int bitpix;
  cpfits.imageDimensions(size);
  
  assert(size.size() ==2);
  
  nx = size[0];
  ny = size[1];
  
  double phys_res;
  
  if(cpfits.readKey("CD1_1",angular_pixel_size)){
    std::cerr << "LensMap fits map must have header keywords:" << std::endl
    << " CD1_1 - angular resolution must exit" << std::endl;
    
    exit(1);
  }
  angular_pixel_size *= degreesTOradians;
  phys_res = angular_pixel_size*angDist;
  boxlMpc = phys_res*nx;
  
  center *= 0;
  lowerleft = center;
  lowerleft[0] -= boxlMpc/2;
  lowerleft[1] -= phys_res*ny/2;
  
  upperright = center;
  upperright[0] += boxlMpc/2;
  upperright[1] += phys_res*ny/2.;
}

void LensMap::read(std::string fits_input_file,double angDist){
  
  std:: cout << " reading lens density map file: " << fits_input_file << std:: endl;
  
  CPFITS_READ cpfits(fits_input_file);
  
  int n_images = cpfits.get_num_hdus();
  
  std::vector<long> dims;
  cpfits.read(surface_density, dims);
  assert(dims.size() == 2);
  
  nx = dims[0];
  ny = dims[1];
  
  if(n_images > 1){  // file contains other lensing quantities
    
    assert(n_images == 6);
    
    cpfits.change_hdu(2);
    cpfits.read(alpha1_bar,dims);
    
    cpfits.change_hdu(3);
    cpfits.read(alpha2_bar,dims);
    
    cpfits.change_hdu(4);
    cpfits.read(gamma1_bar,dims);
    
    cpfits.change_hdu(5);
    cpfits.read(gamma2_bar,dims);
    
    cpfits.change_hdu(6);
    cpfits.read(phi_bar,dims);
  }
  
  if( cpfits.readKey ("CD1_1",angular_pixel_size) ) {// angular resolution degrees
    
    /* these are always present in each*/
    float res;
    angular_pixel_size *= degreesTOradians;
    res = angular_pixel_size*angDist;
    boxlMpc = res * nx;
  }else if(!cpfits.readKey("PHYSICALSIZE",boxlMpc)){
      angular_pixel_size = boxlMpc / angDist;
  }else{
      std::cerr << "LensMap fits map must have header keywords:" << std::endl
      << " CD1_1 - angular resolution" << std::endl
      << " or " << std::endl
      << " PHYSICALSIZE - physical size of map in Mpc = resolution * nx " << std::endl;
    exit(1);
  }
  
  double phys_res = boxlMpc/ nx;
  center *= 0;
  
  lowerleft = center;
  lowerleft[0] -= boxlMpc/2;
  lowerleft[1] -= phys_res * ny /2;
  
  upperright = center;
  upperright[0] += boxlMpc/2;
  upperright[1] += phys_res * ny /2.;
}

void LensMap::Myread(std::string fits_input_file){
  
  std:: cout << " reading lens density map file: "
  << fits_input_file << std:: endl;
  CPFITS_READ cpfits(fits_input_file);
  
  assert(cpfits.get_num_hdus() == 5);
  double yrange;
  cpfits.readKey("SIDEL1",boxlMpc);
  cpfits.readKey("SIDEL2",yrange);
  cpfits.readKey("CENTER_X",center[0]);
  cpfits.readKey("CENTER_Y",center[1]);
  
  std::vector<long> dims;
  cpfits.read(surface_density,dims);
  assert(dims.size() == 2 );
  nx = dims[0];
  ny = dims[1];
  
  cpfits.change_hdu(2);
  cpfits.read(alpha1_bar,dims);
  
  cpfits.change_hdu(3);
  cpfits.read(alpha2_bar,dims);
  
  cpfits.change_hdu(4);
  cpfits.read(gamma1_bar,dims);
  
  cpfits.change_hdu(5);
  cpfits.read(gamma2_bar,dims);
  
  lowerleft = center;
  lowerleft[0] -= boxlMpc/2;
  lowerleft[1] -= yrange/2;
  
  upperright = center;
  upperright[0] += boxlMpc/2;
  upperright[1] += yrange/2.;
}


/*
 void LensMap::read_sub(std::string fits_input_file
 ,const std::vector<long> &first   /// 2d vector for pixel of lower left, (1,1) offset
 ,const std::vector<long> &last    /// 2d vector for pixel of upper right, (1,1) offset
 ,double angDist
 ){
 
 //std:: cout << " reading lens density map file: " << fits_input_file << std:: endl;
 std::unique_ptr<CCfits::FITS> ff(new CCfits::FITS (fits_input_file, CCfits::Read));
 CCfits::PHDU &h0 = ff->pHDU();
 
 h0.readAllKeys();
 
 // these are always present in each
 float res;
 h0.readKey ("CD1_1",res);  // resolution in
 //h0.readKey ("REDSHIFT",z);
 //h0.readKey ("WLOW",wlow);
 //h0.readKey ("WUP",wup);
 
 //  double D = 3*(pow(wup,4) - pow(wlow,4))/(pow(wup,3) - pow(wlow,3))/4;
 //  if(wlow == wup) D = wup;
 //  D /= (1+z)*h;
 res *= degreesTOradians*angDist;
 boxlMpc = res*nx;
 
 assert(h0.axes() >= 2);
 std::vector<long> stride = {1,1};
 
 // principal HDU is read
 //std::cout << surface_density.size() << std::endl;
 h0.read(surface_density,first,last,stride);
 //std::cout << surface_density.size() << std::endl;
 
 // set the full field lower left
 lowerleft = {-boxlMpc/2,-res*ny/2};
 upperright = lowerleft;
 
 
 nx = h0.axis(0);
 ny = h0.axis(1);
 
 ff.release();
 
 lowerleft[0] += (first[0]-1)*res ;
 lowerleft[1] += (first[1]-1)*res ;
 
 //lowerleft[0] =   boxlMpc*(2*first[0] - nx)/2 - boxlMpc/2;
 //lowerleft[1] = boxlMpc*(2*first[1] - ny)/2 - res*ny/2;
 
 upperright[0] += last[0]*res;
 upperright[1] += last[1]*res;
 
 center = (lowerleft + upperright)/2;
 
 nx = last[0] - first[0] + 1;
 ny = last[1] - first[1] + 1;
 
 boxlMpc = res*nx;
 }
 */

void LensMap::read_sub(CPFITS_READ &cpfits
                       ,std::vector<long> &first   // 1 to n
                       ,std::vector<long> &last    // 1 to n
                       ,double angDist
                       ){
  
  //int bitpix;
  std::vector<long> sizes(2);
  cpfits.imageDimensions(sizes);
  
  //long nx_orig = h0.axis(0);
  //long ny_orig = h0.axis(1);
  
  // these are always present in each
  //float wlow,wup,res;
  float res;
  cpfits.readKey("CD1_1",res);
  
  res *= degreesTOradians * angDist;
  //double boxlMpc_orig = res * nx_orig;
  double boxlMpc_orig = res * sizes[0];
  
  assert(sizes.size() >= 2);
  std::vector<long> stride = {1,1};
  
  surface_density.resize( (last[0]-first[0]+1) * (last[1]-first[1]+1) );
  // principal HDU is read
  //std::cout << surface_density.size() << std::endl;
  //h0.read(surface_density,first,last,stride);
  cpfits.read_subset(&(surface_density[0])
                     ,first.data(),last.data());
  
  //std::cout << surface_density.size() << std::endl;
  
  // set the full field lower left
  //lowerleft = upperright = {-boxlMpc_orig/2,-res*ny_orig/2};
  lowerleft = upperright = {-boxlMpc_orig/2,-res*sizes[1]/2};
  
  lowerleft[0] += (first[0]-1)*res ;
  lowerleft[1] += (first[1]-1)*res ;
  
  upperright[0] += last[0]*res;
  upperright[1] += last[1]*res;
  
  center = (lowerleft + upperright)/2;
  
  nx = last[0] - first[0] + 1;
  ny = last[1] - first[1] + 1;
  
  boxlMpc = res*nx;
  
  //std::cout << "Subs map resolution : " << x_resolution() << " " << y_resolution() << std::endl;
}

/**
 * \brief write the fits file of the maps of all the lensing quantities.
 *
 * Unlike the read operations this will not have the h factors in everything so
 * when reading from a file created by is you should set h to 1
 */
void LensMap::write(std::string filename
                    ){
  CPFITS_WRITE cpfits(filename,false);
  
  std::vector<long> naxex(2);
  naxex[0]=nx;
  naxex[1]=ny;
  cpfits.write_image(surface_density,naxex);
  
  //PHDU *phout=&fout->pHDU();
  cpfits.writeKey("SIDEL1",boxlMpc,"x range in physical Mpc");
  cpfits.writeKey("SIDEL2",y_range(),"y range in physical Mpc");
  cpfits.writeKey("CENTER_X",center[0],"center of field x");
  cpfits.writeKey("CENTER_Y",center[1],"center of field y");
  cpfits.writeKey("QUANTITY","KAPPA","lensing quantity");
  cpfits.writeKey("UNITS"," Mpc^-2","units of values");
  
  cpfits.write_image(alpha1_bar,naxex);
  cpfits.writeKey("QUANTITY","ALPHA1","lensing quantity");
  cpfits.writeKey("UNITS"," Mpc^-1","units of values");
  cpfits.write_image(alpha2_bar,naxex);
  cpfits.writeKey("QUANTITY","ALPHA2","lensing quantity");
  cpfits.writeKey("UNITS"," Mpc^-1","units of values");
  cpfits.write_image(gamma1_bar,naxex);
  cpfits.writeKey("QUANTITY","GAMMA1","lensing quantity");
  cpfits.writeKey("UNITS"," Mpc^-2","units of values");
  cpfits.write_image(gamma2_bar,naxex);
  cpfits.writeKey("QUANTITY","GAMMA2","lensing quantity");
  cpfits.writeKey("UNITS"," Mpc^-2","units of values");
}

void LensMap::write(std::string filename
                    ,LensingVariable quant){
  
  if( boxlMpc != angular_pixel_size*nx){
    ERROR_MESSAGE();
    std::cerr << " This function is meant for and angular map and not a physical unit map " << std::endl;
    throw std::invalid_argument("");
  }
  
  CPFITS_WRITE cpfits(filename,false);
  
  std::vector<long> naxex(2);
  naxex[0]=nx;
  naxex[1]=ny;
  
  switch (quant) {
    case KAPPA:
      cpfits.write_image(surface_density,naxex);
      break;
    case GAMMA1:
      cpfits.write_image(gamma1_bar,naxex);
      break;
    case GAMMA2:
      cpfits.write_image(gamma2_bar,naxex);
      break;
    case GAMMA:
    {
      std::valarray<float>  gamma =  sqrt( gamma1_bar*gamma1_bar + gamma2_bar*gamma2_bar );
      cpfits.write_image(gamma,naxex);
    }
      break;
    case ALPHA1:
      cpfits.write_image(alpha1_bar,naxex);
      cpfits.writeKey("UNITS","radians","units of values");
      break;
    case ALPHA2:
      cpfits.write_image(alpha2_bar,naxex);
      cpfits.writeKey("UNITS","radians","units of values");
      break;
    default:
      break;
  }
  
  cpfits.writeKey("CD1_1",angular_pixel_size /degreesTOradians,"pixel size in degrees");
  cpfits.writeKey("CD1_1",angular_pixel_size /degreesTOradians,"pixel size in degrees");
  cpfits.writeKey("SIDEL1",boxlMpc,"x range in radians");
  cpfits.writeKey("SIDEL2",y_range(),"y range in radians");
  cpfits.writeKey("CENTER_X",center[0],"center of field x");
  cpfits.writeKey("CENTER_Y",center[1],"center of field y");
  
}


void MOKAmap::read(std::string MOKA_input_file,bool zeromean,const COSMOLOGY &cosmo){
  
  std:: cout << " reading lens density map file: " << MOKA_input_file << std:: endl;
  CPFITS_READ cpfits(MOKA_input_file);
  std::vector<long> size;
  cpfits.read(surface_density,size);
  
  assert(size.size() == 2);
  
  nx = size[0];
  ny = size[1];
  
  //size_t size = nx*ny;
  //surface_density.resize(size);
  
  // try to read MVIR, if it exists is a MOKA map
  bool moka;
  
  if(cpfits.readKey ("MVIR",m)){
    moka=true;
  }else{
    moka=false;
  }
  int n_images = cpfits.get_num_hdus();
  
  if(moka){
    // check if other quantities exist if not it has to compute them
    
    if(n_images >= 5){  // file contains other lensing quantities
      
      cpfits.change_hdu(2);
      cpfits.read(alpha1_bar,size);
      
      cpfits.change_hdu(3);
      cpfits.read(alpha2_bar,size);
      
      cpfits.change_hdu(4);
      cpfits.read(gamma1_bar,size);
      
      cpfits.change_hdu(5);
      cpfits.read(gamma2_bar,size);
      
      /*
       alpha1_bar.resize(size);
       alpha2_bar.resize(size);
       gamma1_bar.resize(size);
       gamma2_bar.resize(size);
       //gamma3_bar.resize(size);
       phi_bar.resize(size);
       
       ExtHDU &h1=ff->extension(1);
       h1.read(alpha1_bar);
       
       ExtHDU &h2=ff->extension(2);
       h2.read(alpha2_bar);
       
       ExtHDU &h3=ff->extension(3);
       h3.read(gamma1_bar);
       ExtHDU &h4=ff->extension(4);
       h4.read(gamma2_bar);
       std::cout << *h0 << h1 << h2 << h3  << h4 << std::endl;
       */
      
      if(n_images == 6){
        cpfits.change_hdu(6);
        cpfits.read(phi_bar,size);
      }
    }
    int error = 0;
    
    // these are always present in each fits file created by MOKA
    error += cpfits.readKey ("SIDEL",boxlarcsec);
    error += cpfits.readKey ("SIDEL2",boxlMpc);  // recall you that MOKA Mpc/h
    error += cpfits.readKey ("ZLENS",zlens);
    error += cpfits.readKey ("ZSOURCE",zsource);
    error += cpfits.readKey ("OMEGA",omegam);
    error += cpfits.readKey ("LAMBDA",omegal);
    error += cpfits.readKey ("H",h);
    error += cpfits.readKey ("W",wq);
    error += cpfits.readKey ("MSTAR",mstar);
    // error += cpfits.readKey ("MVIR",m);
    error += cpfits.readKey ("CONCENTRATION",c);
    error += cpfits.readKey ("DL",Dlens);
    error += cpfits.readKey ("DLS",DLS);
    error += cpfits.readKey ("DS",DS);
    
    if(error != 0){
      std::cerr << "MOKA fits map must have header keywords:" << std::endl
      << " SIDEL - length on a side in Mpc/h" << std::endl
      << " SIDEL2 - length on other side in Mpc/h" << std::endl
      << " ZLENS - redshift of lens" << std::endl
      << " ZSOURCE - redshift of source" << std::endl
      << " OMEGA - Omega matter" << std::endl
      << " LAMBDA - Omega lambda" << std::endl
      << " H - hubble constant" << std::endl
      << " W - " << std::endl
      << " MSTAR - MSTAR" << std::endl
      << " CONCENTRATION - " << std::endl
      << " DL - " << std::endl
      << " DLS - " << std::endl
      << " DS - " << std::endl
      << " W - " << std::endl;
      exit(1);
    }
    
  }else{  // Pixelized mass map
    
    int npixels = nx;
    
    // keep it like it is, even if it is a rectangle
    if(ny!=nx){
      std:: cout << " " << std:: endl;
      std:: cout << " the plane maps are rectangles! " << std:: endl;
      std:: cout << " " << std:: endl;
    }
    // FITS file must contain the following keywords:
    /*
     ZLENS or REDSHIFT                                 double
     PHYSICALSIZE                                      double
     */
    
    {
      double d1, d2;
      
      // for(int i=0;i<ni;i++) std:: cout << zi[i] << "  " << dli[i] << std:: endl;
      // std:: cout << dli[ni-1] << std:: endl;
      // exit(1);
      int error = cpfits.readKey("WLOW",d1);
      d1=d1/cosmo.gethubble();
      
      if(error){
        if(cpfits.readKey("DLLOW",d1)){
          d1 = d2 = 0;
        }
        
        error = cpfits.readKey("WUP",d2);
        d2=d2/cosmo.gethubble();
        if(error){
          if(cpfits.readKey("DLUP",d2)){
            d1 = d2 = 0;
          }
        }
        if(d2 != 0 ){
          zlens = cosmo.invComovingDist(( d1 + d2 )*0.5);
          
        }else{
          // if angular size distances are not set use ZLENS or REDSHIFT
          if(cpfits.readKey("ZLENS",zlens)){
            if(cpfits.readKey("REDSHIFT",zlens)){
              std::cout << "unable to read fits mass map header keywords" << std::endl <<  "  either DLUP and DLLOW need to be set or ZLENS or REDSHIFT" << std::endl;
              exit(1);
            }
          }
          
        }
        
      }
    }
    if(zlens <= 0.0){
      std::cerr << "Pixel map lens planes cannot have zero or negative redshifts!" << std::endl;
      throw std::runtime_error("Invalid header");
    }
    
    Dlens = cosmo.angDist(0.,zlens);  // physical
    double inarcsec  = 180./M_PI/Dlens*60.*60.;
    double pixLMpc,pixelunit;
    
    // physical size in degrees
    int error = cpfits.readKey("PHYSICALSIZE",boxlarcsec);
    boxlarcsec=boxlarcsec*60.*60.;
    pixLMpc = boxlarcsec/npixels/inarcsec;
    boxlMpc = pixLMpc*npixels;
    if(error){
      // physical size in degrees
      error = cpfits.readKey("PHYSICAL",boxlarcsec);
      boxlarcsec=boxlarcsec*60.*60.;
      pixLMpc = boxlarcsec/npixels/inarcsec;
      boxlMpc = pixLMpc*npixels;
      
      if(error){
        std::cerr << "fits mass map must have header keywords:" << std::endl
        << " REDSHIFT - redshift of map plane" << std::endl
        //<< " DLOW - ?? Mpc/h" << std::endl
        //<< " DLUP - ?? Mpc/h" << std::endl
        << " PHYSICALSIZE - size of map in the x-direction (degrees)" << std::endl
        << " PIXELUNIT - pixel units in solar masses" << std::endl;
        
        std::cout << " unable to read map PIXELUNITS" << std::endl;
        exit(1);
      }
    }
    
    error = cpfits.readKey("PIXELUNIT",pixelunit);
    pixelunit=pixelunit/pixLMpc/pixLMpc;
    
    if(error) {
      
      error = cpfits.readKey("PIXELUNI",pixelunit);
      pixelunit=pixelunit/pixLMpc/pixLMpc;
      
      if(error) {
        std::cerr << "fits mass map must have header keywords:" << std::endl
        << " REDSHIFT - redshift of map plane" << std::endl
        //<< " DLOW - ?? Mpc/h" << std::endl
        //<< " DLUP - ?? Mpc/h" << std::endl
        << " PHYSICALSIZE - size of map in the x-direction (degrees)" << std::endl
        << " PIXELUNIT - pixel units in solar masses" << std::endl;
        
        std::cout << " unable to read map PIXELUNITS" << std::endl;
        exit(1);
      }
    }
    
    for(int i=0;i<nx;i++) for(int j=0;j<ny;j++){
      surface_density[i+nx*j] *= pixelunit;
    }
    
    if(zeromean){
      double ave_sigma = 0;
      
      for(int i=0;i<nx;i++) for(int j=0;j<ny;j++){
        ave_sigma += surface_density[i+nx*j];
      }
      
      ave_sigma /= (nx*ny);
      
      for(int i=0;i<nx;i++) for(int j=0;j<ny;j++){
        surface_density[i+nx*j] -= ave_sigma;
      }
    }
    
    // kappa is not divided by the critical surface density
    // they don't need to be preprocessed by fact
    // create alpha and gamma arrays by FFT
    // valid only to force the map to be square nx = ny = npixels;
    //std:: cout << "  preProcessing Map " << std:: endl;
    // reducing the size of the maps
    // carlo test begin
    /*
     int nl = 2048;
     std:: cout << " resizing the map to " << nl << std:: endl;
     std:: vector<double> xh,xl;
     fill_linear(xh,npixels,0.,1.);
     fill_linear(xl,nl,0.,1.);
     std:: valarray<float> kl(nl*nl);
     for(int i=0;i<nl;i++) for(int j=0;j<nl;j++){
     double xi = xl[i];
     double yi = xl[j];
     
     double ki = getMapVal(convergence,npixels,xi,yi,xh,xh);
     kl[i+nl*j] = ki;
     if(ki!=ki){
     std:: cout << xi << "  " << yi << "  " << ki << std:: endl;
     exit(1);
     }
     }
     std:: cout << "  done " << std:: endl;
     nx = ny = nl;
     convergence.resize(nl*nl);
     
     std:: cout << "  remapping " << std:: endl;
     for(int i=0;i<nl;i++) for(int j=0;j<nl;j++){
     
     convergence[i+nl*j] = kl[i+nl*j];
     }
     */
    // carlo test end
    //PreProcessFFTWMap(zerosize);
  }
  // std:: cout << boxlMpc << "  " << boxlarcsec << std:: endl;
  
  for(size_t i=0;i<surface_density.size();++i){
    assert(surface_density[i] == surface_density[i]);
  }
}
void MOKAmap::write(std::string filename){
  
  CPFITS_WRITE cpfits(filename,false);
  
  std::vector<long> naxex(2);
  naxex[0]=nx;
  naxex[1]=ny;
  
  cpfits.write_image(surface_density, naxex);
  //PHDU *phout=&fout->pHDU();
  //phout->write( 1,nx*ny,surface_density );
  
  cpfits.writeKey ("SIDEL",boxlarcsec,"arcsec");
  cpfits.writeKey ("SIDEL2",boxlMpc,"Mpc/h");
  cpfits.writeKey ("ZLENS",zlens,"lens redshift");
  cpfits.writeKey ("ZSOURCE",zsource, "source redshift");
  cpfits.writeKey ("OMEGA",omegam,"omega matter");
  cpfits.writeKey ("LAMBDA",omegal,"omega lamda");
  cpfits.writeKey ("H",h,"hubble/100");
  cpfits.writeKey ("W",wq,"dark energy equation of state parameter");
  cpfits.writeKey ("MSTAR",mstar,"stellar mass of the BCG in Msun/h");
  cpfits.writeKey ("MVIR",m,"virial mass of the halo in Msun/h");
  cpfits.writeKey ("CONCENTRATION",c,"NFW concentration");
  cpfits.writeKey ("DL",Dlens,"Mpc/h");
  cpfits.writeKey ("DLS",DLS,"Mpc/h");
  cpfits.writeKey ("DS",DS,"Mpc/h");
  
  cpfits.write_image(gamma1_bar, naxex);
  cpfits.write_image(gamma2_bar, naxex);
  
  //ExtHDU *eh1=fout->addImage("gamma1", FLOAT_IMG, naxex);
  //eh1->write(1,nx*ny,gamma1_bar);
  //ExtHDU *eh2=fout->addImage("gamma2", FLOAT_IMG, naxex);
  //eh2->write(1,nx*ny,gamma2_bar);
  //ExtHDU *eh3=fout->addImage("gamma3", FLOAT_IMG, naxex);
  //eh3->write(1,nx*ny,gamma3_bar);
  //std::cout << *phout << std::endl;
}

/**
 * \brief pre-process surface mass density map computing deflection angles and shear in FFT,
 *  generalized to work with rectangular maps
 */
void MOKAmap::PreProcessFFTWMap(float zerosize){
  // initialize the quantities
  //int npix_filter = 0;   // filter the map if you want on a given number of pixels: CHECK IT NOT IMPLEMENTED YET
  
  // size of the new map in x and y directions, factor by which each size is increased
  
  int Nnx=int(zerosize*nx);
  int Nny=int(zerosize*ny);
  double Nboxlx = boxlMpc*zerosize;
  double Nboxly = boxlMpc*zerosize/nx*ny;
  
  std:: valarray<float> Nmap;
  try{
    Nmap.resize( Nnx*Nny );
  }catch(std::exception &e){
    std::cerr << "exception thrown in MOKAmap::PreProcessFFTWMap(): " << e.what() << std::endl;
  }
  // assume locate in a rectangular map and build up the new one
  for( int j=0; j<Nny; j++ ){
    for( int i=0; i<Nnx; i++ ){
      Nmap[i+Nnx*j]=0;
      if(i>=int(Nnx/2-nx/2) && i<int(Nnx/2+nx/2) && j>=int(Nny/2-ny/2) && j<int(Nny/2+ny/2)){
        int ii = i-int(Nnx/2-nx/2);
        int jj = j-int(Nny/2-ny/2);
        
        if(ii>=nx || jj>=ny){
          std::cout << " 1 error mapping " << ii << "  " << jj << std::endl;
          exit(1);
        }
        if(ii<0 || jj<0){
          std::cout << " 2 error mapping " << ii << "  " << jj << std::endl;
          exit(1);
        }
        Nmap[i+Nnx*j] = surface_density[ii+nx*jj];
      }
    }
  }
  
  double *dNmap=new double[Nnx*Nny];
  double *input=new double[Nnx*Nny];
  fftw_complex *fNmap=new fftw_complex[Nny*(Nnx/2+1)];
  fftw_complex *output=new fftw_complex[Nny*(Nnx/2+1)];
  for(int k=0;k<Nnx*Nny;k++) dNmap[k] = double(Nmap[k]);
  fftw_plan p;
  p=fftw_plan_dft_r2c_2d(Nny,Nnx,input,output,FFTW_ESTIMATE);
  
  for(int i=0;i<Nnx*Nny;i++) input[i] = dNmap[i];
  fftw_execute( p );
  for(int i=0; i<Nny*(Nnx/2+1);i++){
    fNmap[i][0] = output[i][0];
    fNmap[i][1] = output[i][1];
  }
  delete[] input;
  delete[] output;
  delete[] dNmap;
  fftw_destroy_plan(p);
  
  // fourier space
  // std:: cout << " allocating fourier space maps " << std:: endl;
  
  
  fftw_complex *fphi   = new fftw_complex[Nny*(Nnx/2+1)];
  // build modes for each pixel in the fourier space
  for( int i=0; i<Nnx/2+1; i++ ){
    double kx=double(i);
    kx=kx*2.*M_PI/Nboxlx;
    for( int j=0; j<Nny; j++ ){
      double ky=(j<Nny/2)?double(j):double(j-Nny);
      ky=ky*2.*M_PI/Nboxly;
      double k2 = kx*kx+ky*ky;
      //smooth if you want to IMPLEMENT
      // if(npix_filter>0){
      // fNmap[j+(Nnpixels/2+1)*i][0] = fNmap[j+(Nnpixels/2+1)*i][0]*exp(-k2/sigmag/sigmag/2.);
      // fNmap[j+(Nnpixels/2+1)*i][1] = fNmap[j+(Nnpixels/2+1)*i][1]*exp(-k2/sigmag/sigmag/2.);
      // }
      
      // fphi
      fphi[i+(Nnx/2+1)*j][0]= -2.*fNmap[i+(Nnx/2+1)*j][0]/k2;
      fphi[i+(Nnx/2+1)*j][1]= -2.*fNmap[i+(Nnx/2+1)*j][1]/k2;
      // null for k2 = 0 no divergence
      if(k2 == 0){
        fphi[i+(Nnx/2+1)*j][0] = 0.;
        fphi[i+(Nnx/2+1)*j][1] = 0.;
      }
    }
  }
  
  delete[] fNmap;
  
  
  fftw_complex *fft= new fftw_complex[Nny*(Nnx/2+1)];
  double *realsp = new double[Nnx*Nny];
  //fftw_plan pp = fftw_plan_dft_c2r_2d(Nny,Nnx,fft,realsp,FFTW_ESTIMATE);
  fftw_plan pp = fftw_plan_dft_c2r_2d(Nny,Nnx,fft,realsp,FFTW_MEASURE);
  
  // alpha1
  {
    
    // build modes for each pixel in the fourier space
    for( int i=0; i<Nnx/2+1; i++ ){
      double kx=double(i);
      kx=kx*2.*M_PI/Nboxlx;
      for( int j=0; j<Nny; j++ ){
        double ky=(j<Nny/2)?double(j):double(j-Nny);
        ky=ky*2.*M_PI/Nboxly;
        
        fft[i+(Nnx/2+1)*j][0] = -kx*fphi[i+(Nnx/2+1)*j][1];
        fft[i+(Nnx/2+1)*j][1] =  kx*fphi[i+(Nnx/2+1)*j][0];
      }
    }
    
    //pp=fftw_plan_dft_c2r_2d(Nny,Nnx,fft,realsp,FFTW_ESTIMATE);
    fftw_execute( pp );
    //fftw_destroy_plan(pp);
    
    alpha1_bar.resize(nx*ny);
    
    for( int j=Nny/2-ny/2; j<Nny/2+ny/2; j++ ){
      for( int i=Nnx/2-nx/2; i<Nnx/2+nx/2; i++ ){
        int ii = i-int(Nnx/2-nx/2);
        int jj = j-int(Nny/2-ny/2);
        alpha1_bar[ii+nx*jj] = -1*float(realsp[i+Nnx*j]/Nnx/Nny);
      }
    }
  }
  
  // alpha2
  {
    
    // build modes for each pixel in the fourier space
    for( int i=0; i<Nnx/2+1; i++ ){
      double kx=double(i);
      kx=kx*2.*M_PI/Nboxlx;
      for( int j=0; j<Nny; j++ ){
        double ky=(j<Nny/2)?double(j):double(j-Nny);
        ky=ky*2.*M_PI/Nboxly;
        
        // alpha
        fft[i+(Nnx/2+1)*j][0] = -ky*fphi[i+(Nnx/2+1)*j][1];
        fft[i+(Nnx/2+1)*j][1] =  ky*fphi[i+(Nnx/2+1)*j][0];
        
      }
    }
    
    //pp=fftw_plan_dft_c2r_2d(Nny,Nnx,fft,realsp,FFTW_ESTIMATE);
    fftw_execute( pp );
    //fftw_destroy_plan(pp);
    
    alpha2_bar.resize(nx*ny);
    
    for( int j=Nny/2-ny/2; j<Nny/2+ny/2; j++ ){
      for( int i=Nnx/2-nx/2; i<Nnx/2+nx/2; i++ ){
        int ii = i-int(Nnx/2-nx/2);
        int jj = j-int(Nny/2-ny/2);
        
        alpha2_bar[ii+nx*jj] = -1*float(realsp[i+Nnx*j]/Nnx/Nny);
      }
    }
  }
  // gamma1
  {
    
    // build modes for each pixel in the fourier space
    for( int i=0; i<Nnx/2+1; i++ ){
      double kx=double(i);
      kx=kx*2.*M_PI/Nboxlx;
      for( int j=0; j<Nny; j++ ){
        double ky=(j<Nny/2)?double(j):double(j-Nny);
        ky=ky*2.*M_PI/Nboxly;
        
        // gamma
        fft[i+(Nnx/2+1)*j][0] = 0.5*(kx*kx-ky*ky)*fphi[i+(Nnx/2+1)*j][0];
        fft[i+(Nnx/2+1)*j][1] = 0.5*(kx*kx-ky*ky)*fphi[i+(Nnx/2+1)*j][1];
        
      }
    }
    
    //pp=fftw_plan_dft_c2r_2d(Nny,Nnx,fft,realsp,FFTW_ESTIMATE);
    fftw_execute( pp );
    //fftw_destroy_plan(pp);
    
    gamma1_bar.resize(nx*ny);
    
    for( int j=Nny/2-ny/2; j<Nny/2+ny/2; j++ ){
      for( int i=Nnx/2-nx/2; i<Nnx/2+nx/2; i++ ){
        int ii = i-int(Nnx/2-nx/2);
        int jj = j-int(Nny/2-ny/2);
        
        gamma1_bar[ii+nx*jj] = float( realsp[i+Nnx*j]/Nnx/Nny);
        
      }
    }
  }
  // gamma2
  {
    
    // build modes for each pixel in the fourier space
    for( int i=0; i<Nnx/2+1; i++ ){
      double kx=double(i);
      kx=kx*2.*M_PI/Nboxlx;
      for( int j=0; j<Nny; j++ ){
        double ky=(j<Nny/2)?double(j):double(j-Nny);
        ky=ky*2.*M_PI/Nboxly;
        
        // gamma
        fft[i+(Nnx/2+1)*j][0] = kx*ky*fphi[i+(Nnx/2+1)*j][0];
        fft[i+(Nnx/2+1)*j][1] = kx*ky*fphi[i+(Nnx/2+1)*j][1];
        
      }
    }
    
    //pp=fftw_plan_dft_c2r_2d(Nny,Nnx,fft,realsp,FFTW_ESTIMATE);
    fftw_execute( pp );
    //fftw_destroy_plan(pp);
    
    gamma2_bar.resize(nx*ny);
    
    for( int j=Nny/2-ny/2; j<Nny/2+ny/2; j++ ){
      for( int i=Nnx/2-nx/2; i<Nnx/2+nx/2; i++ ){
        int ii = i-int(Nnx/2-nx/2);
        int jj = j-int(Nny/2-ny/2);
        
        gamma2_bar[ii+nx*jj] = float(-realsp[i+Nnx*j]/Nnx/Nny);
      }
    }
  }
  // phi - this is done over because of the window in Fourier space
  {
    
    for( int i=0; i<Nnx/2+1; i++ ){
      for( int j=0; j<Nny; j++ ){
        
        fft[i+(Nnx/2+1)*j][0] = fphi[i+(Nnx/2+1)*j][0];
        fft[i+(Nnx/2+1)*j][1] = fphi[i+(Nnx/2+1)*j][1];
        
      }
    }
    
    fftw_execute( pp );
    
    phi_bar.resize(nx*ny);
    for( int j=Nny/2-ny/2; j<Nny/2+ny/2; j++ ){
      for( int i=Nnx/2-nx/2; i<Nnx/2+nx/2; i++ ){
        int ii = i-int(Nnx/2-nx/2);
        int jj = j-int(Nny/2-ny/2);
        
        phi_bar[ii+nx*jj] = float(-realsp[i+Nnx*j]/Nnx/Nny);
      }
    }
    
  }
  
  
  // std:: cout << " remapping the map in the original size " << std:: endl;
  delete[] fft;
  delete[] realsp;
  delete[] fphi;
}

