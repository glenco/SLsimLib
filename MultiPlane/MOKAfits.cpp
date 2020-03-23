/*
 * fits.cpp
 *
 *  Created on: Jun 19, 2012
 *      Author: mpetkova
 */

#include "MOKAlens.h"
#include <fstream>

#include "cpfits.h"

#ifdef ENABLE_FFTW
#include "fftw3.h"
#endif

// carlo test begin
/*
 int getiMIN(int a, int b){
 return (a) < (b) ? (a):(b);
 }
 
 int getiMAX(int a, int b){
 return (a) > (b) ? (a):(b);
 }
 
 double getMapVal(std::valarray<double> &map,int npix,double x,double y, std:: vector<double> &xi,std:: vector<double> &yi){
 size_t i=locate (xi,x);
 i=getiMIN( getiMAX( i, 0 ), npix-2 );
 size_t j=locate (yi,y);
 j=getiMIN( getiMAX( j, 0 ), npix-2 );
 size_t tnpix = npix;
 double fQ11 = map[i+tnpix*j];
 double fQ21;
 if(i<npix-1) fQ21 = map[(i+1)+tnpix*j];
 else fQ21 = 0;
 double fQ12;
 if(j<npix-1) fQ12 = map[i+tnpix*(j+1)];
 else fQ12 = 0;
 double fQ22;
 if(i<npix-1 && j<npix-1) fQ22 = map[(i+1)+tnpix*(j+1)];
 else fQ22 = 0;
 double fxy;
 //if(i<npix-1 && j<npix-1){
 fxy = fQ11*(xi[i+1]-x)*(yi[j+1]-y)/((xi[i+1]-xi[i])*(yi[j+1]-yi[j])) +
 fQ21*(x-xi[i])*(yi[j+1]-y)/((xi[i+1]-xi[i])*(yi[j+1]-yi[j])) +
 fQ12*(xi[i+1]-x)*(y-yi[j])/((xi[i+1]-xi[i])*(yi[j+1]-yi[j])) +
 fQ22*(x-xi[i])*(y-yi[j])/((xi[i+1]-xi[i])*(yi[j+1]-yi[j]));
 //}
 
 
 if(fxy!=fxy){
 std:: cout << x << "  " << y  << std:: endl;
 std:: cout << i << "  " << j << "  " << npix << std:: endl;
 std:: cout << "x = " << xi[i] << "  " << xi[i+1] << std:: endl;
 std:: cout << "y = " << yi[j] << "  " << yi[j+1] << std:: endl;
 std:: cout << fQ11 << "  " << fQ21 << "  " << fQ12 << "  " << fQ22 << std:: endl;
 }
 return fxy;
 }
 */
// carlo test end

void LensHaloMassMap::getDims(){
  
  CPFITS_READ cpfits(MOKA_input_file);
  std::vector<long> size;
  //int bitpix;
  cpfits.imageDimensions(size);
 
  map.nx = size[0];
  map.ny = size[1];
  std:: cout << "nx           ny " << std:: endl;
  std:: cout << map.nx << "   " << map.ny << std:: endl;
}

/**
 * \brief reads in the fits file for the MOKA or mass map and saves it in the structure map
 */

void LensHaloMassMap::readMap(){
  map.read(MOKA_input_file,zeromean,cosmo);
  if(map.alpha1_bar.size() == 0){
    map.PreProcessFFTWMap(zerosize);
  }
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
      throw std::runtime_error("bad file");
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


/**
 * \brief write the fits file of the new MOKA map from the structure map
 */
void LensHaloMassMap::writeImage(std::string filename){
  map.write(filename);
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
 * routine used by fof to link nearby grid cell points
 */
void make_friendship(int ii,int ji,int np,std:: vector<int> &friends, std:: vector<double> &pointdist){
  for(int jj=0;jj<np;jj++){
    if(friends[ji+np*jj]!=0){
      if(friends[ji+np*jj]<0){
        friends[ii+np*jj]=-(ii+1);
      }
      else{
        friends[ii+np*jj]=(ii+1);
      }
      friends[ji+np*jj]=0;
    }
  }
  friends[ii+np*ji]=-(ii+1);
}

/*
 * given a a set of grid points xci and yci and an interpixeld distance l return the id of the
 * fof group nearest to the centre
 */

int fof(double l,std:: vector<double> xci, std:: vector<double> yci, std:: vector<int> &groupid){
  int np = xci.size();
  std:: vector<int> friends(np*np);
  std:: vector<double> pointdist(np*np);
  for(int ii = 0;ii<np; ii++) for(int ji = 0;ji<np; ji++){
    pointdist[ii+np*ji] = sqrt( pow(xci[ii] - xci[ji],2) + pow(yci[ii] - yci[ji],2));
    groupid[ii] = 0;
    friends[ii+np*ji]=0;
  }
  for(int ii=0;ii<np;ii++) for(int ji = 0;ji<np; ji++){
    if(pointdist[ii+np*ji]<=1.5*l) friends[ii+np*ji] = ii+1;
  }
  for(int ii=0;ii<np;ii++){
    int r = 0;
    while(r==0){
      r=1;
      for(int ji=0;ji<np;ji++){
        if(friends[ii+np*ji]>0){
          if(ii!=ji){
            make_friendship(ii,ji,np,friends,pointdist);
            r=0;
          }
        }
      }
    }
  }
  for(int ii=0;ii<np;ii++){
    int p=0;
    for(int ji=0;ji<np;ji++){
      if(friends[ji+np*ii]!=0) p++;
      if(p==2){
        std:: cout << ji << "  " << ii << ":  "  << friends[ji+np*ii] << "  " << friends[ii+np*ji] << std:: endl;
        exit(1);
      }
    }
  }
  // count the particles in each group
  int kt = 0;
  int ng= 0;
  std:: vector<double> distcentre;
  std:: vector<int> idgroup;
  for(int ii=0;ii<np;ii++){
    int k = 0;
    double xcm=0;
    double ycm=0;
    for(int ji=0;ji<np;ji++){
      if(friends[ii+np*ji]!=0){
        k++;
        groupid[ji]=ii+1;
        xcm += xci[ji];
        ycm += yci[ji];
      }
    }
    if(k>4){
      ng++;
      xcm/=k;
      ycm/=k;
      distcentre.push_back(sqrt(xcm*xcm+ycm*ycm));
      idgroup.push_back(ii+1);
      // std:: cout << "  " << ii+1 << "  " << k << "   " << sqrt(xcm*xcm+ycm*ycm) << std:: endl;
    }
    kt = kt + k;
  }
  if(kt != np){
    std:: cout << " number of screaned particles : " << kt << std:: endl;
    std:: cout << " differes from the number of particles : " << np << std:: endl;
    std:: cout << " number of group found : " << ng << std:: endl;
    std:: cout << "     " << std:: endl;
    std:: cout << " I will STOP here!!! " << std:: endl;
    exit(1);
  }
  if(idgroup.size()>0){
    std:: vector<double>::iterator it = min_element(distcentre.begin(), distcentre.end());
    int minpos = idgroup[distance(distcentre.begin(), it)];
    // std:: cout << " nearest to the centre " << minpos << std:: endl;
    return minpos;
  }
  else return 0;
  /* Make a histogram of the data */
  /*
   std::vector< int > histogram(np,0);
   std::vector< int >::iterator it = groupid.begin();
   while(it != groupid.end()) histogram[*it++]++;
   // Print out the frequencies of the values in v
   std::copy(histogram.begin(),histogram.end(),std::ostream_iterator< int >(std::cout, " "));
   std::cout << std::endl;
   // Find the mode
   int mode = std::max_element(histogram.begin(),histogram.end()) - histogram.begin();
   return mode;
   */
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


