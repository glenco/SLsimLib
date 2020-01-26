/*
 * multimap.cpp
 *
 */

#include "slsimlib.h"
#include "multimap.h"
#include "Tree.h"

using namespace std;

LensHaloMultiMap::LensHaloMultiMap(
                 std::string fitsfile  /// Original fits map of the density
                 ,std::string dir_data
                 ,double redshift
                 ,double mass_unit
                 ,COSMOLOGY &c
                 ,bool write_subfields
                 ,std::string dir_scratch
                 ,bool subtract_ave
                 ,bool single_grid_mode
                 ):
LensHalo(redshift,c),write_shorts(write_subfields),single_grid(single_grid_mode),cosmo(c),cpfits(dir_data + fitsfile)
,ave_ang_sd(0),mass_unit(mass_unit),fitsfilename(dir_data + fitsfile)
{

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
  
  resolution = submap.x_resolution();
  angular_resolution = submap.angular_pixel_size;
  
  //border_width = 4.5*sqrt(rs2)/res + 1;
  border_width = ffactor * sqrt(rs2) / resolution + 1;
 
  int desample = sqrt(0.5*submap.nx) / sqrt(ffactor * gfactor);

  Noriginal[0] = submap.nx;
  Noriginal[1] = submap.ny;

  string lr_file = fitsfile + "_lr.fits";
  if(dir_scratch == ""){
    lr_file = dir_data + lr_file;
  }else{
    lr_file = dir_scratch + lr_file;
  }
  bool long_range_file_exists = Utilities::IO::file_exists(lr_file);
  
  ///long_range_file_exists = false;
  if( long_range_file_exists && !single_grid ){
  
    std::cout << " reading file " << lr_file << " .. " << std::endl;
    long_range_map.Myread(lr_file);
    CPFITS_READ cpfits(lr_file);
    if(cpfits.readKey("ave_ang_sd", ave_ang_sd)){
      std::cerr << "need ave_ang_sd in " << lr_file << std::endl;
      throw std::invalid_argument("need ave_ang_sd in ");
    }

  }else if(single_grid){
    std::vector<long> first = {1,1};
    std::vector<long> last = {(long)Noriginal[0],(long)Noriginal[1]};
    
    //long_range_map.read_sub(ff,first,last,getDist());
    long_range_map.read_sub(cpfits,first,last,getDist());

    //double res = submap.boxlMpc/submap.nx;
    //Point_2d dmap = {submap.boxlMpc,submap.ny*resolution};
    //dmap /= 2;
    //long_range_map.center = {0,0};
    //long_range_map.upperright = long_range_map.center + dmap;
    //long_range_map.lowerleft = long_range_map.center - dmap;
    //long_range_map.boxlMpc = submap.boxlMpc;
    
    double area = long_range_map.x_resolution()
    *long_range_map.y_resolution()/mass_unit;
    
    //double area = 1.0/mass_unit;
    
    // convert to
    ave_ang_sd = 0;
    for(auto &p : long_range_map.surface_density){
      p /= area;
      ave_ang_sd += p;
    }
    ave_ang_sd /= long_range_map.surface_density.size();
   if(subtract_ave){
       for(auto &p : long_range_map.surface_density){
        p -= ave_ang_sd;
      }
    }

    long_range_map.PreProcessFFTWMap<LensMap::UNIT>(1.0,unit);

  }else{

    long_range_map.lowerleft = submap.lowerleft;
    long_range_map.upperright = submap.upperright;
    long_range_map.boxlMpc = submap.boxlMpc;
    
    long_range_map.center = (long_range_map.upperright + long_range_map.lowerleft)/2;

    Point_2d range = long_range_map.upperright - long_range_map.lowerleft;
    
    long_range_map.nx = submap.nx/desample;
    double lr_res_x = long_range_map.boxlMpc / long_range_map.nx;

    // get the ny that makes the pixels closest to square
    double Ly = submap.y_range();
    long_range_map.ny = (int)( Ly/lr_res_x );
  
    if(Ly - long_range_map.ny*lr_res_x > lr_res_x/2) long_range_map.ny += 1;

    size_t nx = long_range_map.nx;
    size_t ny = long_range_map.ny;
    
    double lr_res_y = submap.y_range()/ny;
    
    std::cout << "ratio of low resolution pixel dimensions "
    << lr_res_x/lr_res_y << std::endl;
    
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
    
    double res_tmp = resolution / lr_res_y;
    //double res_tmp = resolution / res_y;
    for(size_t j = 0 ; j < Noriginal[1] ; ++j ){

      //jj = j / desample;
      jj = j * res_tmp;
      if( jj >= ny) break;

      //assert(jj < ny);
      //size_t kjj = nx*(j/desample);
      size_t kjj = nx*jj;
    
      if( j%chunk == 0 ){
        first[1] = j+1;
        last[1] = MIN(j + chunk,Noriginal[1]);
        //submap.read_sub(ff,first,last,getDist());
        submap.read_sub(cpfits,first,last,getDist());

        kj = 0;
      }else{
        kj += submap.nx;
        //assert(kj < submap.nx*submap.ny);
      }
    
      double tmp,res_tmp_x = resolution/lr_res_x;
      for(size_t i = 0 ; i < Noriginal[0] ; ++i ){
        size_t ii = i * res_tmp_x;
        //size_t ii = i / desample;
        
        if( ii >= nx) break;
        
        tmp = long_range_map.surface_density[ ii + kjj ] += submap.surface_density[ i + kj ];

        max_pix = MAX(max_pix,tmp);
        min_pix = MIN(min_pix,tmp);
      }
    }
    assert(jj == ny-1);
 
    double area = long_range_map.x_resolution()
    *long_range_map.y_resolution()/mass_unit;
    
    // convert to surface density
    ave_ang_sd = 0;
    for(auto &p : long_range_map.surface_density){
      p /= area;
      ave_ang_sd += p;
    }

    ave_ang_sd /= long_range_map.surface_density.size();

    if(subtract_ave){
      for(auto &p : long_range_map.surface_density){
        p -= ave_ang_sd;
      }
    }
    
    double padd_lr = 1 + 2 * border_width * resolution / long_range_map.x_range();
    padd_lr = MAX(padd_lr,2);
    
    long_range_map.PreProcessFFTWMap<LensMap::WLR>(padd_lr,wlr);
    long_range_map.write("!" + lr_file);
 
    CPFITS_WRITE tmp_cpfits(lr_file,true);
    tmp_cpfits.writeKey("ave_ang_sd", ave_ang_sd,"average angulare density");
  }

  /// set prefix to subfield maps
  size_t lastindex = fitsfile.find_last_of(".");
  
  subfield_filename = fitsfile.substr(0, lastindex) + "_sr";
  if(dir_scratch == ""){
    subfield_filename = dir_data + subfield_filename;
  }else{
    subfield_filename = dir_scratch + subfield_filename;
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
  
  //std::cout << "ur = " << ur << std::endl;
  //std::cout << "dim. original " << Noriginal[0]
  //<< " " << Noriginal[1] << std::endl;
  assert(ur[0] < long_range_map.x_range()/resolution + 0.1);
  assert(ur[1] < long_range_map.y_range()/resolution + 0.1);
  
  lower_left[0] = floor(ll[0]);
  lower_left[1] = floor(ll[1]);

  upper_right[0] = floor(ur[0]) + 1 ;
  upper_right[1] = floor(ur[1]) + 1 ;

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
  
  /// construct file name and have it printed
  std::string sr_file = subfield_filename + to_string(lower_left[0]) + "-" + to_string(lower_left[1]) + "-" +
  to_string(upper_right[0]) + "-" + to_string(upper_right[1]) + "_sr.fits";
  bool short_range_file_exists = Utilities::IO::file_exists(sr_file);

  if(short_range_file_exists && write_shorts){
    short_range_map.Myread( sr_file );
    
    short_range_map.lowerleft = short_range_map.upperright = long_range_map.lowerleft;
    short_range_map.lowerleft[0] += lower_left[0]*resolution;
    short_range_map.lowerleft[1] += lower_left[1]*resolution;
    
    short_range_map.upperright[0] += (upper_right[0] + 1)*resolution;
    short_range_map.upperright[1] += (upper_right[1] + 1)*resolution;
    
    short_range_map.center = (short_range_map.lowerleft + short_range_map.upperright)/2;
    short_range_map.boxlMpc = short_range_map.nx*resolution;
    
  }else{
  
    std::vector<long> first(2);
    std::vector<long> last(2);
  
    first[0] = lower_left[0] - border_width + 1;
    first[1] = lower_left[1] - border_width + 1;

    last[0] = upper_right[0] + border_width + 1;
    last[1] = upper_right[1] + border_width + 1;
  
    double area = resolution*resolution/mass_unit;
    LensMap map;

    if( (first[0] > 1)*(first[1] > 1)*(last[0] <= Noriginal[0])*(last[1] <= Noriginal[1]) ){
 
      map.nx = last[0] - first[0] + 1;
      map.ny = last[1] - first[1] + 1;

      map.surface_density.resize( map.nx * map.ny );
      cpfits.read_subset(&map.surface_density[0],first.data(),last.data() );
      
      // convert to relative surface angular surface density
      for(auto &p : map.surface_density){
        p /= area;
        p -= ave_ang_sd;
      }
      
      map.boxlMpc = map.nx*resolution;

    }else{
    
      // case where subfield overlaps edge
      size_t nx_big = map.nx = last[0] - first[0] + 1;
      size_t ny_big = map.ny = last[1] - first[1] + 1;
    
      std::vector<long> first_sub(2);
      std::vector<long> last_sub(2);

      map.surface_density.resize(nx_big * ny_big,0);
    
      std::valarray<float> v;
      first_sub[0] = MAX(first[0],1);
      first_sub[1] = MAX(first[1],1);
      last_sub[0] = MIN(last[0],Noriginal[0]);
      last_sub[1] = MIN(last[1],Noriginal[1]);
      
      Point_2d offset(first_sub[0] - first[0] , first_sub[1] - first[1] );
 
      long nx = last_sub[0] - first_sub[0] + 1;
      long ny = last_sub[1] - first_sub[1] + 1;
      
      v.resize(nx*ny);
      cpfits.read_subset(&v[0], first_sub.data(), last_sub.data() );
      
      long jj =  offset[1];
      for(long j = 0; j < ny ; ++j,++jj){
        long ii =  offset[0];
        for(long i = 0; i < nx ; ++i,++ii){
          map.surface_density[ii + jj*nx_big] = v[i + j*nx] / area - ave_ang_sd;;
        }
      }
    
      // need to do overlap region
      map.boxlMpc = map.nx*resolution;
    }
 
  map.PreProcessFFTWMap(wsr);

    // cut off bounders
  
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
 
  if(write_shorts) short_range_map.write("!" + sr_file);
}

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

