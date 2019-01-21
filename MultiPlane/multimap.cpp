/*
 * multimap.cpp
 *
 */

#include "slsimlib.h"
#include "multimap.h"

using namespace std;
using namespace CCfits;

LensHaloMultiMap::LensHaloMultiMap(
                 std::string fitsfile  /// Original fits map of the density
                 ,double z
                 ,double mass_unit
                 ,const COSMOLOGY &c
                 ):LensHalo(),cosmo(c),mass_unit(mass_unit),fitsfilename(fitsfile)
{
  
  ff = new CCfits::FITS (fitsfilename, CCfits::Read);
  
  zerosize = 1;
  
  setZlens(z);
  
  rscale = 1.0;
  
  LensHalo::setTheta(0,0);
  setZlensDist(z,cosmo);
  
  if(std::numeric_limits<float>::has_infinity)
    Rmax = std::numeric_limits<float>::infinity();
  else
    Rmax = std::numeric_limits<float>::max();
  
  LensMap submap;
  submap.read_header(fitsfilename,cosmo.gethubble());
  
  long_range_map.boxlMpc = submap.boxlMpc;
  
  //std::size_t size = bigmap.nx*bigmap.ny;
  rs2 = submap.boxlMpc*submap.boxlMpc/( 2*PI*submap.nx );
  res = submap.boxlMpc/submap.nx;
  Noriginal[0] = submap.nx;
  Noriginal[1] = submap.ny;
  border_width = 4.5*sqrt(rs2)/res + 1;
  
  int desample = 2. * PI * sqrt(rs2) / res / 3. ;
  
  desample = sqrt(rs2) / res; // ???? test
  
  size_t nx = long_range_map.nx = submap.nx / desample ;
  size_t ny = long_range_map.ny = nx * submap.ny * 1./ submap.nx ;
  
  {
    Point_2d dmap = {submap.boxlMpc,ny*submap.boxlMpc/nx};
    dmap /= 2;
    long_range_map.center = {0,0};
    long_range_map.upperright = long_range_map.center + dmap;
    long_range_map.lowerleft = long_range_map.center - dmap;
  }
  
  long_range_map.surface_density.resize(nx*ny,0);
  
  long chunk = MIN(nx*ny/Noriginal[0],Noriginal[1]);  // same size as long range grid
  if( chunk == 0) chunk = MIN(100*desample,Noriginal[1]);
  
  std::vector<long> first(2);
  std::vector<long> last(2);
  
  first[0] = 1;
  last[0] = Noriginal[0];
  
  size_t kj = 0;
  for(size_t j = 0 ; j < Noriginal[1] ; ++j ){
    size_t jj = MIN(j/desample,ny-1); // ???
    //assert(jj < ny);
    //size_t kjj = nx*(j/desample);
    size_t kjj = nx*jj;
    
    if( j%chunk == 0 ){
      first[1] = j+1;
      last[1] = MIN(j + chunk,Noriginal[1]);
      //submap.read_sub(fitsfilename,first,last,cosmo.gethubble());
      submap.read_sub(ff,first,last,cosmo.gethubble());
      kj = 0;
    }else{
      kj += submap.nx;
      //assert(kj < submap.nx*submap.ny);
    }
    
    for(size_t i = 0 ; i < Noriginal[0] ; ++i ){
      size_t ii = MIN(i/desample,nx-1);
      //assert(ii + kjj < long_range_map.surface_density.size() );
      //assert(i + kj < submap.surface_density.size() );
      long_range_map.surface_density[ ii + kjj ] += submap.surface_density[ i + kj ];
      //assert(!isnan(submap.surface_density[ i + kj ]));
    }
  }
  
  double area = res*res/mass_unit; //*** units  ???
  // convert to
  for(auto &p : long_range_map.surface_density){
    p /= area;
    //assert(!isnan(p));
  }
  
  wlr.rs2 = wsr.rs2 = rs2;
  long_range_map.PreProcessFFTWMap<WLR>(1.0,wlr);
  
  // ????
//  UNIT w;
//  long_range_map.PreProcessFFTWMap<UNIT>(1.0,w);

  long_range_map.write("!" + fitsfile + "_lr.fits");
  
  //for(auto p : long_range_map.surface_density){ // ??????
  //  assert(!isnan(p));
  //}

  // ??? don
  //bigmap.boxlMpc = long_range_map.boxlMpc;
  //bigmap.PreProcessFFTWMap(1.0,wsr);
  //bigmap.write("!" + fitsfile + "_sr.fits");
};

void LensHaloMultiMap::submap(Point_2d ll,Point_2d ur){
  std::vector<long> lower_left(2);
  std::vector<long> upper_right(2);

  ll = (ll - long_range_map.lowerleft)/res;
  ur = (ur - long_range_map.lowerleft)/res;
  lower_left[0] = lround(ll[0]);
  lower_left[1] = lround(ll[1]);

  upper_right[0] = lround(ur[0]);
  upper_right[1] = lround(ur[1]);

  submap(lower_left,upper_right);
}

void LensHaloMultiMap::submap(
                              const std::vector<long> &lower_left
                              ,const std::vector<long> &upper_right
                              ){
  
  // check range
  
  if( (upper_right[0] < 0) || (upper_right[1] < 0) || (lower_left[0] >= Noriginal[0] )
      || (lower_left[1] >= Noriginal[1])  ){
    std::cerr << "LensHaloMap : sub map is out of bounds" << std::endl;
    throw std::invalid_argument("out of bounds");
  }
  if( (upper_right[0] - lower_left[0]) > Noriginal[0]
     || (upper_right[1] - lower_left[1]) > Noriginal[1] ){
    std::cerr << "LensHaloMap : sub map is too large" << std::endl;
    throw std::invalid_argument("out of bounds");
  }
  
  std::vector<long> first(2);
  std::vector<long> last(2);
  
  first[0] = lower_left[0] - border_width;
  first[1] = lower_left[1] - border_width;

  last[0] = upper_right[0] + border_width;
  last[1] = upper_right[1] + border_width;
  
  LensMap map;

  if( (first[0] > 0)*(first[1] > 0)*(last[0] < Noriginal[0])*(last[1] < Noriginal[1]) ){
    map.read_sub(fitsfilename,first,last,cosmo.gethubble());
  }else{
    
    size_t nx = map.nx = last[0] - first[0] + 1;
    size_t ny = map.ny = last[1] - first[1] + 1;
    
    // case where subfield overlaps edge
    std::vector<long> first_sub(2);
    std::vector<long> last_sub(2);
    
    first_sub[0] = 1;
    last_sub[0] = Noriginal[0];
    
    map.surface_density.resize(nx*ny);
    
    LensMap partmap;
    
    const size_t chunk = nx*ny/Noriginal[0];
    
    size_t kk,k = chunk + 1;
    long jj = first[1];
    for(size_t j = 0 ; j < ny ; ++j , ++jj ){
      size_t kj = nx*j;
      
      if( k >= chunk  || jj < 0 || jj >= Noriginal[1]  ){
        
        if(jj < 0) jj += Noriginal[1];
        if(jj >= Noriginal[1] ) jj -= Noriginal[1];

        first_sub[1] = jj + 1;  last_sub[1] = MIN(jj + chunk + 1,Noriginal[1]);
        //partmap.read_sub(fitsfilename,first_sub,last_sub,cosmo.gethubble());
        partmap.read_sub(ff,first_sub,last_sub,cosmo.gethubble());
        k = 0;
      }else{
        ++k;
      }
      
      kk = k*Noriginal[0];
      long ii = first[0];
      for(size_t i = 0 ; i < nx ; ++i,++ii ){
        if(ii < 0) ii += Noriginal[0];
        if(ii >= Noriginal[0] ) ii -= Noriginal[0];

        //assert( i + kj < map.surface_density.size() );
        //assert( ii + kk < partmap.surface_density.size() );
        map.surface_density[ i + kj ] = partmap.surface_density[ ii + kk ];
      }
    }

  }

  // need to do overlap region
  map.boxlMpc = map.nx*res;

  double area = res*res/mass_unit; //*** units  ???
  // convert to
  for(auto &p : map.surface_density){
    p /= area;
  }
  
  map.PreProcessFFTWMap(1.0,wsr);
  
  // cut off bourders
  
  size_t nx = short_range_map.nx = upper_right[0] - lower_left[0] + 1;
  size_t ny = short_range_map.ny = upper_right[1] - lower_left[1] + 1;
  
  short_range_map.surface_density.resize(nx*ny);
  short_range_map.alpha1_bar.resize(nx*ny);
  short_range_map.alpha2_bar.resize(nx*ny);
  short_range_map.gamma1_bar.resize(nx*ny);
  short_range_map.gamma2_bar.resize(nx*ny);

  for(long j=0 ; j < ny ; ++j){
    long kjj = (j + border_width )*map.nx;
    long kj = nx*j;
    long ii = border_width;
    for(long i=0 ; i < nx ; ++i,++ii){
      short_range_map.surface_density[i + kj] = map.surface_density[ii + kjj];
      short_range_map.alpha1_bar[i + kj] = map.alpha1_bar[ii + kjj];
      short_range_map.alpha2_bar[i + kj] = map.alpha2_bar[ii + kjj];
      short_range_map.gamma1_bar[i + kj] = map.gamma1_bar[ii + kjj];
      short_range_map.gamma2_bar[i + kj] = map.gamma2_bar[ii + kjj];
    }
  }

  short_range_map.lowerleft = short_range_map.upperright = long_range_map.lowerleft;
  short_range_map.lowerleft[0] += lower_left[0]*res;
  short_range_map.lowerleft[1] += lower_left[1]*res;
  
  short_range_map.upperright[0] += upper_right[0]*res;
  short_range_map.upperright[1] += upper_right[1]*res;
  
  short_range_map.center = (short_range_map.lowerleft + short_range_map.upperright)/2;
  short_range_map.boxlMpc = nx*res;
 
  short_range_map.write("!" + fitsfilename + "_sr.fits");
}


void LensHaloMultiMap::force_halo(double *alpha
                                 ,KappaType *kappa
                                 ,KappaType *gamma
                                 ,KappaType *phi
                                 ,double const *xx
                                 ,bool subtract_point
                                 ,PosType screening
                                 )
{
  
  // interpolate from the maps
  
  Utilities::Interpolator<valarray<float> > interp(xx
                                                    ,long_range_map.nx,long_range_map.boxlMpc
                                                    ,long_range_map.ny
                                                    ,long_range_map.ny*long_range_map.boxlMpc/long_range_map.nx
                                                    ,long_range_map.center.x);
 
  //assert(long_range_map.nx == long_range_map.ny);
  
  alpha[0] = interp.interpolate(long_range_map.alpha1_bar);
  alpha[1] = interp.interpolate(long_range_map.alpha2_bar);
  gamma[0] = interp.interpolate(long_range_map.gamma1_bar);
  gamma[1] = interp.interpolate(long_range_map.gamma2_bar);
  gamma[2] = 0.0;
  *kappa = interp.interpolate(long_range_map.surface_density);
  //*phi = interp.interpolate(long_range_map.phi_bar);
 
  //assert(alpha[0] == alpha[0] && alpha[1] == alpha[1]);
  //assert(gamma[0] == gamma[0] && gamma[1] == gamma[1]);
  //assert(kappa == kappa);

  if( (xx[0] > short_range_map.lowerleft[0])*(xx[0] < short_range_map.upperright[0])
     *(xx[1] > short_range_map.lowerleft[1])*(xx[1] < short_range_map.upperright[1])
     ){

    Utilities::Interpolator<valarray<float> > short_interp(xx
                                                         ,short_range_map.nx,short_range_map.boxlMpc
                                                         ,short_range_map.ny
                                                         ,short_range_map.ny*short_range_map.boxlMpc/short_range_map.nx
                                                         ,short_range_map.center.x);


    alpha[0] += short_interp.interpolate(short_range_map.alpha1_bar);
    alpha[1] += short_interp.interpolate(short_range_map.alpha2_bar);
    gamma[0] += short_interp.interpolate(short_range_map.gamma1_bar);
    gamma[1] += short_interp.interpolate(short_range_map.gamma2_bar);
    *kappa += short_interp.interpolate(short_range_map.surface_density);
  //*phi = interp.interpolate(short_range_map.phi_bar);
  
    //assert(alpha[0] == alpha[0] && alpha[1] == alpha[1]);
    //assert(gamma[0] == gamma[0] && gamma[1] == gamma[1]);
    //assert(kappa == kappa);
  }
  
  return;
}

void LensMap::read_header(std::string fits_input_file,float h){
  
   std::auto_ptr<CCfits::FITS> ff(new CCfits::FITS (fits_input_file, CCfits::Read));
  
  //CCfits::PHDU *h0=&ff->pHDU();
  CCfits::PHDU &h0 = ff->pHDU();

  h0.readAllKeys();
  
  assert(h0.axes() >= 2);
  
  nx = h0.axis(0);
  ny = h0.axis(1);
  
  try{
    /* these are always present in ea*/
    float wlow,wup,res;
    h0.readKey ("CD1_1",res);  // resolution in degrees
    h0.readKey ("REDSHIFT",z);
    h0.readKey ("WLOW",wlow);
    h0.readKey ("WUP",wup);
    
    double D = 3*(pow(wup,4) - pow(wlow,4))/(pow(wup,3) - pow(wlow,3))/4;
    if(wlow == wup) D = wlow;
    D /= (1+z)*h;
    //D /= (1+z);

    res *= degreesTOradians*D;
    boxlMpc = res*nx;
  }
  catch(CCfits::HDU::NoSuchKeyword){
    
    std::cerr << "LensMap fits map must have header keywords:" << std::endl
    << " CD1_1 - length on other side in Mpc/h" << std::endl
    << " REDSHIFT - redshift of lens" << std::endl
    << " WLOW - closest radial distance in cMpc/h" << std::endl
    << " WUP - furthest radial distance in cMpc/h" << std::endl;

    exit(1);
  }

  
  lowerleft = center;
  lowerleft[0] -= boxlMpc/2;
  lowerleft[1] -= boxlMpc*ny/nx/2;
  
  upperright = center;
  upperright[0] += boxlMpc/2;
  upperright[1] += boxlMpc*ny/nx/2;
}

void LensMap::read(std::string fits_input_file,float h){
  
  std:: cout << " reading lens density map file: " << fits_input_file << std:: endl;
  std::auto_ptr<CCfits::FITS> ff(new CCfits::FITS (fits_input_file, CCfits::Read));
  
  CCfits::PHDU &h0 = ff->pHDU();
  
  h0.readAllKeys();
  
  assert(h0.axes() >= 2);
  
  nx = h0.axis(0);
  ny = h0.axis(1);
  
  size_t size = nx*ny;
  
  // principal HDU is read
  h0.read(surface_density);
  int nhdu = h0.axes();

  if(nhdu > 2){  // file contains other lensing quantities
    alpha1_bar.resize(size);
    alpha2_bar.resize(size);
    gamma1_bar.resize(size);
    gamma2_bar.resize(size);
    //gamma3_bar.resize(size);
    phi_bar.resize(size);
      
    CCfits::ExtHDU &h1=ff->extension(1);
    h1.read(alpha1_bar);
    CCfits::ExtHDU &h2=ff->extension(2);
    h2.read(alpha2_bar);
    CCfits::ExtHDU &h3=ff->extension(3);
    h3.read(gamma1_bar);
    CCfits::ExtHDU &h4=ff->extension(4);
    h4.read(gamma2_bar);
    std::cout << h0 << h1 << h2 << h3  << h4 << std::endl;
  }
  try{
    /* these are always present in ea*/
    float wlow,wup,res;
    h0.readKey ("CD1_1",res);  // recall you that MOKA Mpc/h
    h0.readKey ("REDSHIFT",z);
    h0.readKey ("WLOW",wlow);
    h0.readKey ("WUP",wup);
    
    double D = 3*(pow(wup,4) - pow(wlow,4))/(pow(wup,3) - pow(wlow,3))/4;
    D /= (1+z)*h;

    res *= degreesTOradians*D;
    boxlMpc = res*nx;
  }
  catch(CCfits::HDU::NoSuchKeyword){
    
    std::cerr << "LensMap fits map must have header keywords:" << std::endl
    << " CD1_1 - length on other side in Mpc/h" << std::endl
    << " REDSHIFT - redshift of lens" << std::endl
    << " WLOW - closest radial distance in cMpc/h" << std::endl
    << " WUP - furthest radial distance in cMpc/h" << std::endl;
    exit(1);
  }
  
  /// ??? set center
  
  lowerleft = center;
  lowerleft[0] -= boxlMpc/2;
  lowerleft[1] -= boxlMpc*ny/nx/2;

  upperright = center;
  upperright[0] += boxlMpc/2;
  upperright[1] += boxlMpc*ny/nx/2;
}

void LensMap::read_sub(std::string fits_input_file
                       ,const std::vector<long> &first
                       ,const std::vector<long> &last
                       ,float h
){
  
  //std:: cout << " reading lens density map file: " << fits_input_file << std:: endl;
  std::unique_ptr<CCfits::FITS> ff(new CCfits::FITS (fits_input_file, CCfits::Read));
  /*std::unique_ptr<CCfits::FITS> ff(nullptr);
  try{
    //std::unique_ptr<CCfits::FITS> ff(new CCfits::FITS (fits_input_file, CCfits::Read));
    
    ff.reset(new CCfits::FITS (fits_input_file, CCfits::Read));
  }
  catch(CCfits::FitsError){
    std::cerr << " error in opening fits file " << fits_input_file << std::endl;
  }*/
  CCfits::PHDU &h0 = ff->pHDU();
  
  h0.readAllKeys();
  
  /* these are always present in ea*/
  float wlow,wup,res;
  h0.readKey ("CD1_1",res);  // resolution in
  h0.readKey ("REDSHIFT",z);
  h0.readKey ("WLOW",wlow);
  h0.readKey ("WUP",wup);
  
  double D = 3*(pow(wup,4) - pow(wlow,4))/(pow(wup,3) - pow(wlow,3))/4;
  if(wlow == wup) D = wup;
  D /= (1+z)*h;
  res *= degreesTOradians*D;
  boxlMpc = res*nx;
  
  assert(h0.axes() >= 2);
  std::vector<long> stride = {1,1};

  // principal HDU is read
  //std::cout << surface_density.size() << std::endl;
  h0.read(surface_density,first,last,stride);
  //std::cout << surface_density.size() << std::endl;

  nx = h0.axis(0);
  ny = h0.axis(1);
  
  ff.release();
  
  lowerleft[0] = boxlMpc*(2*first[0] - nx)/2;
  lowerleft[1] = boxlMpc*(2*first[1] - ny)/2;

  upperright[0] = boxlMpc*(2*last[0] - nx)/2;
  upperright[1] = boxlMpc*(2*last[1] - ny)/2;

  center = (lowerleft + upperright)/2;
  
  nx = last[0] - first[0] + 1;
  ny = last[1] - first[1] + 1;
  
  boxlMpc *= nx*1./h0.axis(0);
}


void LensMap::read_sub(CCfits::FITS *ff
                       ,const std::vector<long> &first
                       ,const std::vector<long> &last
                       ,float h
                       ){
  CCfits::PHDU &h0 = ff->pHDU();
  
  h0.readAllKeys();
  
  // these are always present in each
  float wlow,wup,res;
  h0.readKey ("CD1_1",res);  // resolution in
  h0.readKey ("REDSHIFT",z);
  h0.readKey ("WLOW",wlow);
  h0.readKey ("WUP",wup);
  
  double D = 3*(pow(wup,4) - pow(wlow,4))/(pow(wup,3) - pow(wlow,3))/4;
  if(wlow == wup) D = wup;
  D /= (1+z)*h;
  res *= degreesTOradians*D;
  boxlMpc = res*nx;
  
  assert(h0.axes() >= 2);
  std::vector<long> stride = {1,1};
  
  // principal HDU is read
  //std::cout << surface_density.size() << std::endl;
  h0.read(surface_density,first,last,stride);
  //std::cout << surface_density.size() << std::endl;
  
  nx = h0.axis(0);
  ny = h0.axis(1);
  
  lowerleft[0] = boxlMpc*(2*first[0] - nx)/2;
  lowerleft[1] = boxlMpc*(2*first[1] - ny)/2;
  
  upperright[0] = boxlMpc*(2*last[0] - nx)/2;
  upperright[1] = boxlMpc*(2*last[1] - ny)/2;
  
  center = (lowerleft + upperright)/2;
  
  nx = last[0] - first[0] + 1;
  ny = last[1] - first[1] + 1;
  
  boxlMpc *= nx*1./h0.axis(0);
}

/**
 * \brief write the fits file of the maps of all the lensing quantities
 */
void LensMap::write(std::string filename){
#ifdef ENABLE_FITS
  long naxis=2;
  long naxes[2]={nx,ny};
  
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
  naxex[0]=nx;
  naxex[1]=ny;
  
  PHDU *phout=&fout->pHDU();

  phout->addKey("SIDEL2",boxlMpc,"Mpc/h");
  
  phout->write( 1,nx*ny,surface_density );
  
  ExtHDU *eh1=fout->addImage("alpah1", FLOAT_IMG, naxex);
  eh1->write(1,nx*ny,alpha1_bar);
  ExtHDU *eh2=fout->addImage("alpha2", FLOAT_IMG, naxex);
  eh2->write(1,nx*ny,alpha2_bar);

  ExtHDU *eh3=fout->addImage("gamma1", FLOAT_IMG, naxex);
  eh3->write(1,nx*ny,gamma1_bar);
  ExtHDU *eh4=fout->addImage("gamma2", FLOAT_IMG, naxex);
  eh4->write(1,nx*ny,gamma2_bar);
  
  std::cout << *phout << std::endl;
#else
  std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
  exit(1);
#endif
}

