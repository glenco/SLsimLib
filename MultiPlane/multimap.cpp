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
                 ,double redshift
                 ,double mass_unit
                 ,COSMOLOGY &c
                 ,bool subtract_ave
                 ,bool single_grid_mode
                 ):LensHalo(redshift,c),single_grid(single_grid_mode),cosmo(c),mass_unit(mass_unit),fitsfilename(fitsfile)
{
  
  try{
    ff = new CCfits::FITS (fitsfilename, CCfits::Read);
  }catch(...){
    std::cerr << "CCfits throw an exception: probably file " << fitsfile << " does not exist." << std::endl;
  }
  zerosize = 1;
  rscale = 1.0;
  
  LensHalo::setTheta(0,0);
  
  if(std::numeric_limits<float>::has_infinity)
    Rmax = std::numeric_limits<float>::infinity();
  else
    Rmax = std::numeric_limits<float>::max();
 
  setZlens(redshift);
  setZlensDist(redshift,cosmo);

  double D = LensHalo::getDist();
  
  LensMap submap;
  submap.read_header(fitsfilename,D);
  
  //long_range_map.boxlMpc = submap.boxlMpc;
  
  //std::size_t size = bigmap.nx*bigmap.ny;
  //rs2 = submap.boxlMpc*submap.boxlMpc/( 2*PI*submap.nx );
  //rs2 = 2*4*submap.boxlMpc*submap.boxlMpc/( 2*PI*submap.nx );
  rs2 = submap.boxlMpc*submap.boxlMpc/( 2*submap.nx )*gfactor/ffactor;

  wlr.rs2 = wsr.rs2 = rs2;
  
  //border_width = 4.5*sqrt(rs2)/res + 1;
  border_width = ffactor * sqrt(rs2) / resolution + 1;
 
  int desample = sqrt(0.5*submap.nx) / sqrt(ffactor * gfactor);

  Noriginal[0] = submap.nx;
  Noriginal[1] = submap.ny;

  bool long_range_file_exists = Utilities::IO::file_exists(fitsfile + "_lr.fits");
  
  long_range_file_exists = false; // ?????
  
  if( long_range_file_exists && !single_grid ){
  
    std::cout << " reading file " << fitsfile + "_lr.fits .. " << std::endl;
    long_range_map.Myread(fitsfile + "_lr.fits");

  }else if(single_grid){
    std::vector<long> first = {1,1};
    std::vector<long> last = {(long)Noriginal[0],(long)Noriginal[1]};
    
    long_range_map.read_sub(ff,first,last,getDist());
    
    //double res = submap.boxlMpc/submap.nx;
    //Point_2d dmap = {submap.boxlMpc,submap.ny*resolution};
    //dmap /= 2;
    //long_range_map.center = {0,0};
    //long_range_map.upperright = long_range_map.center + dmap;
    //long_range_map.lowerleft = long_range_map.center - dmap;
    //long_range_map.boxlMpc = submap.boxlMpc;
    
    double area = long_range_map.x_resolution()
    *long_range_map.y_resolution()/mass_unit; //*** units  ???
    
    //double area = 1.0/mass_unit;
    
    // convert to
    double ave = 0;
    for(auto &p : long_range_map.surface_density){
      p /= area;
      ave += p;
    }
    if(subtract_ave){
      ave /= long_range_map.surface_density.size();
      for(float &p : long_range_map.surface_density){
        p -= ave;
      }
    }

    long_range_map.PreProcessFFTWMap<UNIT>(1.0,unit);

  }else{

    long_range_map.lowerleft = submap.lowerleft;
    long_range_map.upperright = submap.upperright;
    long_range_map.boxlMpc = submap.boxlMpc;
    
    long_range_map.center = (long_range_map.upperright + long_range_map.lowerleft)/2;

    Point_2d range = long_range_map.upperright - long_range_map.lowerleft;
    
    long_range_map.nx = submap.nx/desample;
    double res_x = long_range_map.boxlMpc / long_range_map.nx;

    // get the ny that makes the pixels clossest to square
    double Ly = submap.y_range();
    long_range_map.ny = (int)( Ly/res_x );
    //long_range_map.ny = (int)( Ly/res_x ) - 1;    // ????
    if(Ly - long_range_map.ny*res_x > res_x/2) long_range_map.ny += 1;

    size_t nx = long_range_map.nx;
    size_t ny = long_range_map.ny;
    
    double res_y = submap.y_range()/ny;
    
    std::cout << "ratio of low resolution pixel dimensions "
    << res_x/res_y << std::endl;
    
    size_t  N = long_range_map.ny*long_range_map.nx;
    long_range_map.surface_density.resize(N,0);
  
    long chunk = MIN(N/Noriginal[0],Noriginal[1]);  // same size as long range grid
    if( chunk == 0) chunk = MIN(100*desample,Noriginal[1]);
  
    std::vector<long> first(2);
    std::vector<long> last(2);
  
    first[0] = 1;
    last[0] = Noriginal[0];
  
    size_t kj = 0;
    size_t jj;
    for(size_t j = 0 ; j < Noriginal[1] ; ++j ){
      jj = j*resolution/res_y;
      if( jj >= ny) break;

      //assert(jj < ny);
      //size_t kjj = nx*(j/desample);
      size_t kjj = nx*jj;
    
      if( j%chunk == 0 ){
        first[1] = j+1;
        last[1] = MIN(j + chunk,Noriginal[1]);
        //submap.read_sub(fitsfilename,first,last,cosmo.gethubble());
        submap.read_sub(ff,first,last,getDist());

        kj = 0;
      }else{
        kj += submap.nx;
        //assert(kj < submap.nx*submap.ny);
      }
    
      double tmp;
      for(size_t i = 0 ; i < Noriginal[0] ; ++i ){
        size_t ii = i*resolution/res_x;
        if( ii >= nx) break;
        tmp = long_range_map.surface_density[ ii + kjj ] += submap.surface_density[ i + kj ];

        max_pix = MAX(max_pix,tmp);
        min_pix = MIN(min_pix,tmp);
      }
    }
    assert(jj == ny-1);
 
    //  double area = pow(long_range_map.boxlMpc/long_range_map.nx,2)/mass_unit; //*** units  ???
    //double area = long_range_map.x_resolution()
    //*long_range_map.y_resolution()/mass_unit/resolution/resolution; //*** units  ???
    double area = long_range_map.x_resolution()
    *long_range_map.y_resolution()/mass_unit; //*** units  ???
    //double area = 1.0/mass_unit;

    // convert to
    double ave = 0;
    for(auto &p : long_range_map.surface_density){
      p /= area;
      ave += p;
    }
    if(subtract_ave){
      ave /= long_range_map.surface_density.size();
      for(float &p : long_range_map.surface_density){
        p -= ave;
      }
    }
    
    long_range_map.PreProcessFFTWMap<WLR>(1.0,wlr);
    long_range_map.write("!" + fitsfile + "_lr.fits");
  }
};

void LensHaloMultiMap::submapPhys(Point_2d ll,Point_2d ur){
  if(single_grid) return;
  
  std::vector<long> lower_left(2);
  std::vector<long> upper_right(2);

  //Point_2d range = long_range_map.upperright - long_range_map.lowerleft;
  
  //std::cout << "ang res = " << resolution/getDist()/arcsecTOradians
  //<< " arcsec" << std::endl;
  //std::cout << "lower left angle " << ll/getDist() << std::endl;
  //std::cout << "lower left angle " << long_range_map.lowerleft/getDist() << std::endl;
  
  
  ll = (ll - long_range_map.lowerleft)/resolution;
  //std::cout << "lower left pixels = " << ll << std::endl;
  assert(ll[0] > -0.1 );
  assert(ll[1] > -0.1 );
  if(ll[0] < 0) ll[0] = 0;
  if(ll[1] < 0) ll[1] = 0;
  
  ur = (ur - long_range_map.lowerleft)/resolution;
  
  std::cout << "ur = " << ur << std::endl;
  std::cout << "dim. original " << Noriginal[0]
  << " " << Noriginal[1] << std::endl;
  assert(ur[0] < long_range_map.x_range()/resolution + 0.1);
  assert(ur[1] < long_range_map.y_range()/resolution + 0.1);
  
  lower_left[0] = floor(ll[0]);
  lower_left[1] = floor(ll[1]);

  upper_right[0] = floor(ur[0]) - 1 ;
  upper_right[1] = floor(ur[1]) - 1 ;

  submap(lower_left,upper_right);
}

void LensHaloMultiMap::submap(
                              const std::vector<long> &lower_left
                              ,const std::vector<long> &upper_right
                              ){
  
  // check range
  
  if(single_grid) return;

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
  
  first[0] = lower_left[0] - border_width + 1;
  first[1] = lower_left[1] - border_width + 1;

  last[0] = upper_right[0] + border_width + 1;
  last[1] = upper_right[1] + border_width + 1;
  
  LensMap map;

  if( (first[0] > 0)*(first[1] > 0)*(last[0] <= Noriginal[0])*(last[1] <= Noriginal[1]) ){
    //map.read_sub(fitsfilename,first,last,getDist());

    map.read_sub(ff,first,last,getDist());

  }else{
    
    size_t nx_big = map.nx = last[0] - first[0] + 1;
    size_t ny_big = map.ny = last[1] - first[1] + 1;
    
    // case where subfield overlaps edge
    std::vector<long> first_sub(2);
    std::vector<long> last_sub(2);
    
    first_sub[0] = 1;
    last_sub[0] = Noriginal[0];
    
    map.surface_density.resize(nx_big*ny_big);
    
    LensMap partmap;
    
    const size_t chunk = nx_big*ny_big/Noriginal[0];
    
    size_t kk,k = chunk + 1;
    long jj = first[1]-1;
    for(size_t j = 0 ; j < ny_big ; ++j , ++jj ){
      size_t kj = nx_big*j;
      
      if( k >= chunk  || jj < 0 || jj >= Noriginal[1]  ){
        
        if(jj < 0) jj += Noriginal[1];
        if(jj >= Noriginal[1] ) jj -= Noriginal[1];

        first_sub[1] = jj + 1;  last_sub[1] = MIN(jj + chunk + 1,Noriginal[1]);
        //partmap.read_sub(fitsfilename,first_sub,last_sub,cosmo.gethubble());
        partmap.read_sub(ff,first_sub,last_sub,getDist());
 
        k = 0;
      }else{
        ++k;
      }
      
      kk = k*Noriginal[0];
      long ii = first[0]-1;
      double tmp;
      for(size_t i = 0 ; i < nx_big ; ++i,++ii ){
        if(ii < 0) ii += Noriginal[0];
        if(ii >= Noriginal[0] ) ii -= Noriginal[0];

        //assert( i + kj < map.surface_density.size() );
        //assert( ii + kk < partmap.surface_density.size() );
        tmp = map.surface_density[ i + kj ] = partmap.surface_density[ ii + kk ];
        
        max_pix = MAX(max_pix,tmp);
        min_pix = MIN(min_pix,tmp);
      }
    }
    // need to do overlap region
    map.boxlMpc = map.nx*resolution;
  }


  double area = resolution*resolution/mass_unit; //*** units  ???
  // convert to
  //double area = 1.0/mass_unit; //*** units  ???
  // convert to
  for(auto &p : map.surface_density){
    p /= area;
  }
  
  map.PreProcessFFTWMap(1.0,wsr);
  //map.PreProcessFFTWMap(1.0,unit);
  
  // cut off bourders
  
  size_t nx = short_range_map.nx = upper_right[0] - lower_left[0] + 1;
  size_t ny = short_range_map.ny = upper_right[1] - lower_left[1] + 1;
  
  short_range_map.surface_density.resize(nx*ny);
  short_range_map.alpha1_bar.resize(nx*ny);
  short_range_map.alpha2_bar.resize(nx*ny);
  short_range_map.gamma1_bar.resize(nx*ny);
  short_range_map.gamma2_bar.resize(nx*ny);

  long jj = border_width;
  for(long j=0 ; j < ny ; ++j,++jj){
    long kjj = jj*map.nx;
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
  short_range_map.lowerleft[0] += lower_left[0]*resolution;
  short_range_map.lowerleft[1] += lower_left[1]*resolution;
  
  short_range_map.upperright[0] += (upper_right[0] + 1)*resolution;
  short_range_map.upperright[1] += (upper_right[1] + 1)*resolution;
  
  short_range_map.center = (short_range_map.lowerleft + short_range_map.upperright)/2;
  short_range_map.boxlMpc = short_range_map.nx*resolution;
 
  if(!single_grid) short_range_map.write("!" + fitsfilename + "_sr.fits");
}


bool LensMap::evaluate(const double *x,float &sigma,float *gamma,double *alpha) {
  
  //double fx = ((x[0] - center[0])/range + 0.5)*(nx-1);
  //double fy = ((x[1] - center[1])/range_y + 0.5)*(ny-1);
  //std::cout << "(  " << fx << " " << fy << "   ";
  
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
    + c * alpha1_bar[index+1+nx] + d * alpha1_bar[index+nx];
    gamma[1] = a * gamma2_bar[index] + b * gamma2_bar[index+1]
    + c * gamma2_bar[index+1+nx] + d * gamma2_bar[index+nx];
    gamma[2] = 0.0;
    
    return false;
  }
  
  sigma = 0;
  alpha[0] = alpha[1] = 0.0;
  gamma[0] = gamma[1] = gamma[2] = 0.0;
  
  return true;
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

  if(single_grid){
    
    long_range_map.evaluate(xx,*kappa,gamma,alpha);
    return;
  }
  
  if( (xx[0] >= short_range_map.lowerleft[0])*(xx[0] <= short_range_map.upperright[0])
     *(xx[1] >= short_range_map.lowerleft[1])*(xx[1] <= short_range_map.upperright[1])
     ){
    
    long_range_map.evaluate(xx,*kappa,gamma,alpha);
    
    float t_kappa,t_gamma[3];
    double t_alpha[2];
    short_range_map.evaluate(xx,t_kappa,t_gamma,t_alpha);

    alpha[0] += t_alpha[0];
    alpha[1] += t_alpha[1];
    gamma[0] += t_gamma[0];
    gamma[1] += t_gamma[1];
    *kappa += t_kappa;
    
  }else{
    
    alpha[0] = 0;
    alpha[1] = 0;
    gamma[0] = 0;
    gamma[1] = 0;
    *kappa = 0;
  }
  return;
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
  
  return *this;
}
void LensMap::read_header(std::string fits_input_file
                          ,double angDist){
  
   std::auto_ptr<CCfits::FITS> ff(new CCfits::FITS (fits_input_file, CCfits::Read));
  
  //CCfits::PHDU *h0=&ff->pHDU();
  CCfits::PHDU &h0 = ff->pHDU();

  h0.readAllKeys();
  
  assert(h0.axes() >= 2);
  
  nx = h0.axis(0);
  ny = h0.axis(1);
  double phys_res;
  
  try{
    /* these are always present in ea*/
    //float wlow,wup,
   h0.readKey ("CD1_1",phys_res);  // resolution in degrees
    //h0.readKey ("REDSHIFT",z);
    //h0.readKey ("WLOW",wlow);
    //h0.readKey ("WUP",wup);
    
    //double D = 3*(pow(wup,4) - pow(wlow,4))/(pow(wup,3) - pow(wlow,3))/4;
    //if(wlow == wup) D = wlow;
    //D /= (1+z)*h;
    //D /= (1+z);

    phys_res *= degreesTOradians*angDist;
    boxlMpc = phys_res*nx;
  }
  catch(CCfits::HDU::NoSuchKeyword){
    
    std::cerr << "LensMap fits map must have header keywords:" << std::endl
    << " CD1_1 - length on other side in Mpc/h" << std::endl;
    //<< " WLOW - closest radial distance in cMpc/h" << std::endl
    //<< " WUP - furthest radial distance in cMpc/h" << std::endl;

    exit(1);
  }

  //double phys_res = boxlMpc/nx;

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
    //std::cout << h0 << h1 << h2 << h3  << h4 << std::endl;
  }
  try{
    /* these are always present in ea*/
    float res;
    h0.readKey ("CD1_1",res);  // angular resolution degrees
    //h0.readKey ("REDSHIFT",z);
    //h0.readKey ("WLOW",wlow);
    //h0.readKey ("WUP",wup);
    
    //double D = 3*(pow(wup,4) - pow(wlow,4))/(pow(wup,3) - pow(wlow,3))/4;
    //D /= (1+z)*h;

    res *= degreesTOradians*angDist;
    boxlMpc = res * nx;
  }
  catch(CCfits::HDU::NoSuchKeyword){
    
    std::cerr << "LensMap fits map must have header keywords:" << std::endl
    << " CD1_1 - length on other side in Mpc/h" << std::endl;
    //<< " REDSHIFT - redshift of lens" << std::endl
    //<< " WLOW - closest radial distance in cMpc/h" << std::endl
    //<< " WUP - furthest radial distance in cMpc/h" << std::endl;
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
  //int nhdu = h0.axes();
  
  // file contains other lensing quantities
  alpha1_bar.resize(size);
  alpha2_bar.resize(size);
  gamma1_bar.resize(size);
  gamma2_bar.resize(size);
  //phi_bar.resize(size);
    
  CCfits::ExtHDU &h1=ff->extension(1);
  h1.read(alpha1_bar);
  CCfits::ExtHDU &h2=ff->extension(2);
  h2.read(alpha2_bar);
  CCfits::ExtHDU &h3=ff->extension(3);
  h3.read(gamma1_bar);
  CCfits::ExtHDU &h4=ff->extension(4);
  h4.read(gamma2_bar);
  //std::cout << h0 << h1 << h2 << h3  << h4 << std::endl;
  
  // these are always present in each
  //h0.readKey ("REDSHIFT",z);
  double yrange;
  h0.readKey("SIDEL1",boxlMpc);
  h0.readKey("SIDEL2",yrange);

  //double res = boxlMpc/nx;
  center *= 0;
  
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

void LensMap::read_sub(CCfits::FITS *ff
                       ,const std::vector<long> &first
                       ,const std::vector<long> &last
                       ,double angDist
                       ){

  CCfits::PHDU &h0 = ff->pHDU();
  
  long nx_orig = h0.axis(0);
  long ny_orig = h0.axis(1);

  h0.readAllKeys();
  
  // these are always present in each
  //float wlow,wup,res;
  float res;
  h0.readKey ("CD1_1",res);  // resolution in
  //h0.readKey ("REDSHIFT",z);
  //h0.readKey ("WLOW",wlow);
  //h0.readKey ("WUP",wup);
  
  //double D = 3*(pow(wup,4) - pow(wlow,4))/(pow(wup,3) - pow(wlow,3))/4;
  //if(wlow == wup) D = wup;
  //D /= (1+z)*h;
  res *= degreesTOradians * angDist;
  double boxlMpc_orig = res * nx_orig;
  
  assert(h0.axes() >= 2);
  std::vector<long> stride = {1,1};
  
  // principal HDU is read
  //std::cout << surface_density.size() << std::endl;
  h0.read(surface_density,first,last,stride);
  //std::cout << surface_density.size() << std::endl;
  
  // set the full field lower left
  lowerleft = upperright = {-boxlMpc_orig/2,-res*ny_orig/2};
  
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

  phout->addKey("SIDEL1",boxlMpc,"x range in physical Mpc");
  phout->addKey("SIDEL2",y_range(),"y range in physical Mpc");
  
  phout->write( 1,nx*ny,surface_density );
  
  ExtHDU *eh1=fout->addImage("alpah1", FLOAT_IMG, naxex);
  eh1->write(1,nx*ny,alpha1_bar);
  ExtHDU *eh2=fout->addImage("alpha2", FLOAT_IMG, naxex);
  eh2->write(1,nx*ny,alpha2_bar);

  ExtHDU *eh3=fout->addImage("gamma1", FLOAT_IMG, naxex);
  eh3->write(1,nx*ny,gamma1_bar);
  ExtHDU *eh4=fout->addImage("gamma2", FLOAT_IMG, naxex);
  eh4->write(1,nx*ny,gamma2_bar);
  
  //std::cout << *phout << std::endl;
#else
  std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
  exit(1);
#endif
}

