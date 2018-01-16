/*
 * fits.cpp
 *
 *  Created on: Jun 19, 2012
 *      Author: mpetkova
 */

#include "MOKAlens.h"
#include <fstream>

#ifdef ENABLE_FITS
#include <CCfits/CCfits>
using namespace CCfits;

#endif

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
#ifdef ENABLE_FITS
  try{
    std::auto_ptr<FITS> ff(new FITS (MOKA_input_file, Read));
    
    PHDU *h0=&ff->pHDU();
    
    map->nx = h0->axis(0);
    map->ny = h0->axis(1);
    std:: cout << "nx           ny " << std:: endl;
    std:: cout << map->nx << "   " << map->ny << std:: endl;
  }
  catch(FITS::CantOpen){
    std::cout << "can not open " << MOKA_input_file << std::endl;
    exit(1);
  }
#else
  std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
  exit(1);
#endif
}

/**
 * \brief reads in the fits file for the MOKA or mass map and saves it in the structure map
 */
void LensHaloMassMap::readMap(){
#ifdef ENABLE_FITS
  std:: cout << " reading lens density map file: " << MOKA_input_file << std:: endl;
  std::auto_ptr<FITS> ff(new FITS (MOKA_input_file, Read));
  
  PHDU *h0=&ff->pHDU();
  
  h0->readAllKeys();
  
  // try to read MVIR, if it exists is a MOKA map
  bool moka;
  try {
    h0->readKey ("MVIR",map->m);
    moka=true;
  }
  catch(CCfits::HDU::NoSuchKeyword) {
    moka=false;
  }
  
  // principal HDU is read
  h0->read(map->convergence);
  
  if(moka){
    // check if they exist if not it has to compute them
    int nhdu = h0->axis(2);
    
    if(nhdu==1){
      std:: cout << "  preProcessing Map MOKA fits file has only KAPPA map" << std:: endl;
      PreProcessFFTWMap();
    }
    else{
      ExtHDU &h1=ff->extension(1);
      h1.read(map->alpha1);
      ExtHDU &h2=ff->extension(2);
      h2.read(map->alpha2);
      ExtHDU &h3=ff->extension(3);
      h3.read(map->gamma1);
      ExtHDU &h4=ff->extension(4);
      h4.read(map->gamma2);
      std::cout << *h0 << h1 << h2 << h3  << h4 << std::endl;
    }
    try{
      /* these are always present in each fits file created by MOKA */
      h0->readKey ("SIDEL",map->boxlarcsec);
      h0->readKey ("SIDEL2",map->boxlMpc);  // recall you that MOKA Mpc/h
      h0->readKey ("ZLENS",map->zlens);
      h0->readKey ("ZSOURCE",map->zsource);
      h0->readKey ("OMEGA",map->omegam);
      h0->readKey ("LAMBDA",map->omegal);
      h0->readKey ("H",map->h);
      h0->readKey ("W",map->wq);
      h0->readKey ("MSTAR",map->mstar);
      // h0->readKey ("MVIR",map->m);
      h0->readKey ("CONCENTRATION",map->c);
      h0->readKey ("DL",map->Dlens);
      h0->readKey ("DLS",map->DLS);
      h0->readKey ("DS",map->DS);
    }
    catch(CCfits::HDU::NoSuchKeyword){
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
    
    int npixels = map->nx;
    // cut a square!
    // if(map->ny<map->nx) npixels = map->ny;
    // std:: valarray<double> mapbut(map->nx*map->ny);
    // for(int i=0;i<map->nx;i++) for(int j=0;j<map->ny;j++){
    //    mapbut[i+map->nx*j] = map->convergence[i+map->nx*j];
    //  }
    // map->convergence.resize(npixels*npixels);
    // keep it like it is, even if it is a rectangle
    if(map->ny!=map->nx){
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
      try{
        h0->readKey("WLOW",d1);
        d1=d1/cosmo.gethubble();
      }
      catch(CCfits::HDU::NoSuchKeyword){
        try{
          h0->readKey("DLLOW",d1);
        }
        catch(CCfits::HDU::NoSuchKeyword){
          d1 = d2 = 0;
        }
        
        try{
          h0->readKey("WUP",d2);
          d2=d2/cosmo.gethubble();
        }
        catch(CCfits::HDU::NoSuchKeyword){
          try{
            h0->readKey("DLUP",d2);
          }
          catch(CCfits::HDU::NoSuchKeyword){
            d1 = d2 = 0;
          }
        }
        if(d2 != 0 ){
          /*double dll = ( d1 + d2 )*0.5; // comoving dists
          
          std:: vector<double> zi;
          int ni=2048;
          std:: vector<double> dli(ni);
          Utilities::fill_linear(zi,ni,0.,5.); // max redshift should be around 2.5!
          for(int i=0;i<ni;i++) dli[i] = cosmo.angDist(zi[i])*(1+zi[i]);
          
          if(dli[ni-1] < dll){
            std::cerr << "ERROR: redshift table in LensHaloMassMap::readMap() does not extend to high enough redshift" << std::endl;
            throw std::runtime_error("small redshift table");
          }
          
          // set the redshift of the plane half distance between
          // d1 and d2
          map->zlens = Utilities::InterpolateYvec(dli,zi,dll);
          */
          
          map->zlens = cosmo.invComovingDist(( d1 + d2 )*0.5);
          
        }else{
          
          // if angular size distances are not set use ZLENS or REDSHIFT
          
          try {
            h0->readKey("ZLENS",map->zlens);
          }
          catch(CCfits::HDU::NoSuchKeyword) {
            try {
              h0->readKey("REDSHIFT",map->zlens);
            }
            catch(CCfits::HDU::NoSuchKeyword){
              std::cout << "unable to read fits mass map header keywords" << std::endl <<  "  either DLUP and DLLOW need to be set or ZLENS or REDSHIFT" << std::endl;
              exit(1);
            }
            
            
          }
          
        }
        
      }
    }
    if(map->zlens <= 0.0){
      std::cerr << "Pixel map lens planes cannot have zero or negative redshifts!" << std::endl;
      throw std::runtime_error("Invalid header");
    }
    map->Dlens = cosmo.angDist(0.,map->zlens);  // physical
    double inarcsec  = 180./M_PI/map->Dlens*60.*60.;
    double pixLMpc,pixelunit;
    try{
      // physical size in degrees
      h0->readKey("PHYSICALSIZE",map->boxlarcsec);
      map->boxlarcsec=map->boxlarcsec*60.*60.;
      pixLMpc = map->boxlarcsec/npixels/inarcsec;
      map->boxlMpc = pixLMpc*npixels;
    }
    catch(CCfits::HDU::NoSuchKeyword) {
      try{
        // physical size in degrees
        h0->readKey("PHYSICAL",map->boxlarcsec);
        map->boxlarcsec=map->boxlarcsec*60.*60.;
        pixLMpc = map->boxlarcsec/npixels/inarcsec;
        map->boxlMpc = pixLMpc*npixels;
      }
      catch(CCfits::HDU::NoSuchKeyword) {
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
    
    
    try {
      h0->readKey("PIXELUNIT",pixelunit);
      pixelunit=pixelunit/pixLMpc/pixLMpc;
    }
    catch(CCfits::HDU::NoSuchKeyword) {
      
      try {
        h0->readKey("PIXELUNI",pixelunit);
        pixelunit=pixelunit/pixLMpc/pixLMpc;
      }
      catch(CCfits::HDU::NoSuchKeyword) {
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
    /*catch(CCfits::HDU::NoSuchKeyword) {
     
     std::cerr << "fits mass map should have header keywords:" << std::endl
     << " REDSHIFT - redshift of map plane" << std::endl
     //<< " DLOW - ?? Mpc/h" << std::endl
     //<< " DLUP - ?? Mpc/h" << std::endl
     << " PHYSICALSIZE - size of map in the x-direction (degrees)" << std::endl
     << " PIXELUNIT - pixel units in solar masses" << std::endl;
     
     std::cout << "unable to read map PHYSICALSIZE" << std::endl;
     std::cout << "assuming is the MultiDark file" << std::endl;
     map->boxlarcsec = 8.7*60.*60.;    // W1 x-field of view
     // map->boxlarcsec = 5.5*60.*60.;    // W4 x-field of view
     pixLMpc = map->boxlarcsec/npixels/inarcsec;
     map->boxlMpc = pixLMpc*npixels;
     pixelunit = 1.e+10/cosmo.gethubble()/pixLMpc/pixLMpc; // by hand
     }*/
    
    
    
    // made square // need to be
    // 1. take the part located in the left side
    //for(int i=0;i<npixels;i++) for(int j=0;j<npixels;j++){
    //    map->convergence[i+npixels*j] = mapbut[i+map->nx*j]*pixelunit;
    //    avkappa += map->convergence[i+npixels*j];
    //  }
    // 2. take the part located in the right side
    // for(int i=0;i<npixels;i++) for(int j=0;j<npixels;j++){
    //    map->convergence[i+npixels*j] = mapbut[(map->nx-npixels+i)+map->nx*j]*pixelunit;
    //    avkappa += map->convergence[i+npixels*j];
    //  }
    // avkappa /= (npixels*npixels);
    
    for(int i=0;i<map->nx;i++) for(int j=0;j<map->ny;j++){
      map->convergence[i+map->nx*j] *= pixelunit;
    }
    
    if(zeromean){
      double avkappa = 0;
      
      for(int i=0;i<map->nx;i++) for(int j=0;j<map->ny;j++){
        avkappa += map->convergence[i+map->nx*j];
      }
      avkappa /= (map->nx*map->ny);
      
      for(int i=0;i<map->nx;i++) for(int j=0;j<map->ny;j++){
        map->convergence[i+map->nx*j] = (map->convergence[i+map->nx*j] - avkappa);
      }
    }
    
    // kappa is not divided by the critical surface density
    // they don't need to be preprocessed by fact
    // create alpha and gamma arrays by FFT
    // valid only to force the map to be square map->nx = map->ny = npixels;
#ifdef ENABLE_FFTW
    std:: cout << "  preProcessing Map " << std:: endl;
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
     double ki = getMapVal(map->convergence,npixels,xi,yi,xh,xh);
     kl[i+nl*j] = ki;
     if(ki!=ki){
     std:: cout << xi << "  " << yi << "  " << ki << std:: endl;
     exit(1);
     }
     }
     std:: cout << "  done " << std:: endl;
     map->nx = map->ny = nl;
     map->convergence.resize(nl*nl);
     std:: cout << "  remapping " << std:: endl;
     for(int i=0;i<nl;i++) for(int j=0;j<nl;j++){
     map->convergence[i+nl*j] = kl[i+nl*j];
     }
     */
    // carlo test end
    PreProcessFFTWMap();
#else
    std::cout << "Please enable the preprocessor flag ENABLE_FFTW !" << std::endl;
    exit(1);
#endif
  }
  
  // std:: cout << map->boxlMpc << "  " << map->boxlarcsec << std:: endl;
  
  for(size_t i=0;i<map->convergence.size();++i){
    assert(map->convergence[i] == map->convergence[i]);
  }
#else
  
  std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
  exit(1);
#endif
}


/**
 * \brief write the fits file of the new MOKA map from the structure map
 */
void LensHaloMassMap::writeImage(std::string filename){
#ifdef ENABLE_FITS
  long naxis=2;
  long naxes[2]={map->nx,map->ny};
  
  std::auto_ptr<FITS> fout(0);
  
  try{
    fout.reset(new FITS(filename,FLOAT_IMG,naxis,naxes));
  }
  catch(FITS::CantCreate){
    std::cout << "Unable to open fits file " << filename << std::endl;
    ERROR_MESSAGE();
    exit(1);
  }
  
  std::vector<long> naxex(2);
  naxex[0]=map->nx;
  naxex[1]=map->ny;
  
  PHDU *phout=&fout->pHDU();
  
  phout->write( 1,map->nx*map->ny,map->convergence );
  
  phout->addKey ("SIDEL",map->boxlarcsec,"arcsec");
  phout->addKey ("SIDEL2",map->boxlMpc,"Mpc/h");
  phout->addKey ("ZLENS",map->zlens,"lens redshift");
  phout->addKey ("ZSOURCE",map->zsource, "source redshift");
  phout->addKey ("OMEGA",map->omegam,"omega matter");
  phout->addKey ("LAMBDA",map->omegal,"omega lamda");
  phout->addKey ("H",map->h,"hubble/100");
  phout->addKey ("W",map->wq,"dark energy equation of state parameter");
  phout->addKey ("MSTAR",map->mstar,"stellar mass of the BCG in Msun/h");
  phout->addKey ("MVIR",map->m,"virial mass of the halo in Msun/h");
  phout->addKey ("CONCENTRATION",map->c,"NFW concentration");
  phout->addKey ("DL",map->Dlens,"Mpc/h");
  phout->addKey ("DLS",map->DLS,"Mpc/h");
  phout->addKey ("DS",map->DS,"Mpc/h");
  
  ExtHDU *eh1=fout->addImage("gamma1", FLOAT_IMG, naxex);
  eh1->write(1,map->nx*map->ny,map->gamma1);
  ExtHDU *eh2=fout->addImage("gamma2", FLOAT_IMG, naxex);
  eh2->write(1,map->nx*map->ny,map->gamma2);
  ExtHDU *eh3=fout->addImage("gamma3", FLOAT_IMG, naxex);
  eh3->write(1,map->nx*map->ny,map->gamma3);
  
  std::cout << *phout << std::endl;
#else
  std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
  exit(1);
#endif
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
 * \brief pre-process sourface mass density map computing deflection angles and shear in FFT,
 *  generalized to work with rectangular maps
 */

void LensHaloMassMap::PreProcessFFTWMap(){
#ifdef ENABLE_FFTW
  // initialize the quantities
  //int npix_filter = 0;   // filter the map if you want on a given number of pixels: CHECK IT NOT IMPLEMENTED YET
  
  // size of the new map in x and y directions, factor by which each size is increased
  int Nnx=int(zerosize*map->nx);
  int Nny=int(zerosize*map->ny);
  double Nboxlx = map->boxlMpc*zerosize;
  double Nboxly = map->boxlMpc*zerosize/map->nx*map->ny;
  
  std:: valarray<float> Nmap;
  try{
    Nmap.resize( Nnx*Nny );
  }catch(std::exception &e){
    std::cerr << "exception thrown in LensHaloMassMap::PreProcessFFTWMap(): " << e.what() << std::endl;
  }
  // assume locate in a rectangular map and build up the new one
  for( int j=0; j<Nny; j++ ){
    for( int i=0; i<Nnx; i++ ){
      Nmap[i+Nnx*j]=0;
      if(i>=int(Nnx/2-map->nx/2) && i<int(Nnx/2+map->nx/2) && j>=int(Nny/2-map->ny/2) && j<int(Nny/2+map->ny/2)){
        int ii = i-int(Nnx/2-map->nx/2);
        int jj = j-int(Nny/2-map->ny/2);
        
        if(ii>=map->nx || jj>=map->ny){
          std::cout << " 1 error mapping " << ii << "  " << jj << std::endl;
          exit(1);
        }
        if(ii<0 || jj<0){
          std::cout << " 2 error mapping " << ii << "  " << jj << std::endl;
          exit(1);
        }
        Nmap[i+Nnx*j]=map->convergence[ii+map->nx*jj];
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
    
    map->alpha1.resize(map->nx*map->ny);
    
    for( int j=Nny/2-map->ny/2; j<Nny/2+map->ny/2; j++ ){
      for( int i=Nnx/2-map->nx/2; i<Nnx/2+map->nx/2; i++ ){
        int ii = i-int(Nnx/2-map->nx/2);
        int jj = j-int(Nny/2-map->ny/2);
        
        map->alpha1[ii+map->nx*jj] = -1*float(realsp[i+Nnx*j]/Nnx/Nny);
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
    
    map->alpha2.resize(map->nx*map->ny);
    
    for( int j=Nny/2-map->ny/2; j<Nny/2+map->ny/2; j++ ){
      for( int i=Nnx/2-map->nx/2; i<Nnx/2+map->nx/2; i++ ){
        int ii = i-int(Nnx/2-map->nx/2);
        int jj = j-int(Nny/2-map->ny/2);
        
        map->alpha2[ii+map->nx*jj] = -1*float(realsp[i+Nnx*j]/Nnx/Nny);
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
    
    map->gamma1.resize(map->nx*map->ny);
    
    for( int j=Nny/2-map->ny/2; j<Nny/2+map->ny/2; j++ ){
      for( int i=Nnx/2-map->nx/2; i<Nnx/2+map->nx/2; i++ ){
        int ii = i-int(Nnx/2-map->nx/2);
        int jj = j-int(Nny/2-map->ny/2);
        
        map->gamma1[ii+map->nx*jj] = float( realsp[i+Nnx*j]/Nnx/Nny);
        
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
    
    map->gamma2.resize(map->nx*map->ny);
    
    for( int j=Nny/2-map->ny/2; j<Nny/2+map->ny/2; j++ ){
      for( int i=Nnx/2-map->nx/2; i<Nnx/2+map->nx/2; i++ ){
        int ii = i-int(Nnx/2-map->nx/2);
        int jj = j-int(Nny/2-map->ny/2);
        
        map->gamma2[ii+map->nx*jj] = float(-realsp[i+Nnx*j]/Nnx/Nny);
        
      }
    }
  }
  
  /*
   double *phi    = new double[Nnx*Nny];
   fftw_plan pp;
   pp=fftw_plan_dft_c2r_2d(Nny,Nnx,fphi,phi,FFTW_ESTIMATE);
   fftw_execute( pp );
   fftw_destroy_plan(pp);
   delete[] phi;
   */
  
  // std:: cout << " remapping the map in the original size " << std:: endl;
  delete[] fft;
  delete[] realsp;
  delete[] fphi;
#endif
}


