/*
 * multimap.cpp
 *
 */

#include "slsimlib.h"
#include "multimap.h"
#include "Tree.h"

using namespace std;

my_fftw_plan LensHaloMultiMap::plan_short_range;
//fftw_plan LensHaloMultiMap::plan_c2r_short_range;
my_fftw_plan LensHaloMultiMap::plan_long_range;
//fftw_plan LensHaloMultiMap::plan_c2r_long_range;
std::mutex LensHaloMultiMap::mutex_multimap;
bool LensHaloMultiMap::plans_set_up = false;
size_t LensHaloMultiMap::nx_sub = 0;
size_t LensHaloMultiMap::ny_sub = 0;
size_t LensHaloMultiMap::nx_long = 0;
size_t LensHaloMultiMap::ny_long = 0;
size_t LensHaloMultiMap::nx_sub_extended = 0;
size_t LensHaloMultiMap::ny_sub_extended = 0;
int LensHaloMultiMap::count;

LensHaloMultiMap::LensHaloMultiMap(
                 std::string fitsfile  /// Original fits map of the density
                 ,std::string dir_data
                 ,double redshift
                 ,double mass_unit
                 ,double npanes
                 ,int number_of_subfields
                 ,COSMOLOGY &c
                 ,bool write_subfields
                 ,std::string dir_scratch
                 ,bool subtract_ave
                 ,double ffactor
                 ,double gfactor
                 ,double rsmooth2
                 ):
LensHalo(redshift,c)
,ffactor(ffactor),gfactor(gfactor)
,write_shorts(write_subfields)
,cosmo(c),cpfits(dir_data + fitsfile)
,ave_ang_sd(0)
,rsmooth2(rsmooth2)
,mass_unit(mass_unit)
,fitsfilename(dir_data + fitsfile)
{
  
  ++count;
  zerosize = 1;
  rscale = 1.0;
  
  LensHalo::setTheta(0,0);
  
  if(std::numeric_limits<float>::has_infinity)
    Rmax = std::numeric_limits<float>::infinity();
  else
    Rmax = std::numeric_limits<float>::max();
 
  setZlens(redshift,c);
  setZlensDist(redshift,cosmo);

  double D = LensHalo::getDist();
  
  LensMap tmp_map;
  tmp_map.read_header(fitsfilename,D);
   
  resolution_mpc = tmp_map.x_resolution();
  angular_resolution = tmp_map.angular_pixel_size;
  
  //std::size_t size = bigmap.nx*bigmap.ny;
  //rs2 = submap.boxlMpc*submap.boxlMpc/( 2*PI*submap.nx );
  //rs2 = 2*4*submap.boxlMpc*submap.boxlMpc/( 2*PI*submap.nx );
  rs2 = tmp_map.boxlMpc*tmp_map.boxlMpc/( 2*tmp_map.nx )*gfactor/ffactor;
  //wlr.rs2 = wsr.rs2 = rs2;
  wlr.rs2 = rs2 + rsmooth2;
  wsr.rs2 = rs2;
  wsr_smooth.rs2 = rs2;
  wsr_smooth.r_sm2 = rsmooth2;

  //border_width = 4.5*sqrt(rs2)/res + 1;
  border_width_pix = ffactor * sqrt(rs2) / resolution_mpc + 1;
 
  int desample = sqrt(0.5*tmp_map.nx) / sqrt(ffactor * gfactor);

  Noriginal[0] = tmp_map.nx;
  Noriginal[1] = tmp_map.ny;

  string lr_file = fitsfile + "_lr.fits";
  if(dir_scratch == ""){
    lr_file = dir_data + lr_file;
  }else{
    lr_file = dir_scratch + lr_file;
  }
  bool long_range_file_exists = Utilities::IO::file_exists(lr_file);
  double padd_lr = 0;
  ///long_range_file_exists = false;
  if( long_range_file_exists ){
  
    std::cout << " reading file " << lr_file << " .. " << std::endl;
    long_range_map.Myread(lr_file);
    CPFITS_READ cpfits(lr_file);
    if(cpfits.readKey("ave_ang_sd", ave_ang_sd)){
      std::cerr << "need ave_ang_sd in " << lr_file << std::endl;
      throw std::invalid_argument("need ave_ang_sd in ");
    }
    // *****************************************************

    // check that the header information is what is expected
    double tmp_double;
    int tmp_int;
    long tmp_long;
    
    cpfits.readKey("redshift",tmp_double);
    if(fabs((tmp_double-redshift)/redshift) > 1.0e-5 ){
      std::cout << "missmatch z " << tmp_double << " " << redshift << std::endl;
      long_range_file_exists = false;
    }
    cpfits.readKey("mass_unit",tmp_double);
    if(tmp_double != mass_unit){
      std::cout << "missmatch mass_units " << tmp_double << " " << mass_unit << std::endl;
      long_range_file_exists = false;
    }
    cpfits.readKey("cosmology",tmp_int);
    if(tmp_int != int(c.ParamSet())){
      std::cout << "missmatch cosmology " << tmp_int << " " << c.ParamSet() << std::endl;
      long_range_file_exists = false;
    }
    cpfits.readKey("rs2",tmp_double);
    if( fabs((tmp_double-rs2)/rs2) > 1.0e-5 ){
      std::cout << "missmatch rs2 " << tmp_double << " " << rs2 << std::endl;
      long_range_file_exists = false;
    }
    cpfits.readKey("border_width",tmp_long);
    if(fabs((tmp_long-border_width_pix)*1.0/border_width_pix) > 1.0e-5 ){
      std::cout << "missmatch border_width " << tmp_long << " " << border_width_pix << std::endl;
      long_range_file_exists = false;
    }
    if(!long_range_file_exists) std::cout << "Long range file not compatible. Recalculating .. " << std::endl;
    // *****************************************************

 /* }else if(single_grid){
    std::vector<long> first = {1,1};
    std::vector<long> last = {(long)Noriginal[0],(long)Noriginal[1]};
    
    //long_range_map.read_sub(ff,first,last,getDist());
    long_range_map.read_sub(cpfits,first,last,getDist());

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

    long_range_map.ProcessFFTs<UNIT>(unit,plan_r2c_long_range
                                         ,plan_c2r_long_range);
    //long_range_map.PreProcessFFTWMap<UNIT>(unit,mutex_multimap);
*/
  }
  if( !long_range_file_exists ){

    long_range_map.lowerleft = tmp_map.lowerleft;
    long_range_map.upperright = tmp_map.upperright;
    long_range_map.boxlMpc = tmp_map.boxlMpc;
    
    long_range_map.center = (long_range_map.upperright + long_range_map.lowerleft)/2;

    Point_2d range = long_range_map.upperright - long_range_map.lowerleft;
    
    long_range_map.nx = tmp_map.nx/desample;
    double lr_res_x = long_range_map.boxlMpc / long_range_map.nx;

    // get the ny that makes the pixels closest to square
    double Ly = tmp_map.y_range();
    long_range_map.ny = (int)( Ly/lr_res_x );
  
    if(Ly - long_range_map.ny*lr_res_x > lr_res_x/2) long_range_map.ny += 1;

    size_t nx = long_range_map.nx;
    size_t ny = long_range_map.ny;
    
    double lr_res_y = tmp_map.y_range()/ny;
    
    std::cout << "ratio of low resolution pixel dimensions "
    << lr_res_x/lr_res_y << std::endl;
    
    size_t  N = long_range_map.ny*long_range_map.nx;
    long_range_map.surface_density.resize(N,0);
  
    long chunk = MIN(N/Noriginal[0],Noriginal[1]);  // same size as long range grid
    if( chunk == 0) chunk = MIN<long>(100*desample,Noriginal[1]);
  
    std::vector<long> first(2);
    std::vector<long> last(2);
  
    first[0] = 1;
    last[0] = Noriginal[0];
  
    size_t kj = 0;
    size_t jj;
    
    double res_tmp = resolution_mpc / lr_res_y;
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
        tmp_map.read_sub(cpfits,first,last,getDist());

        kj = 0;
      }else{
        kj += tmp_map.nx;
        //assert(kj < submap.nx*submap.ny);
      }
    
      double tmp,res_tmp_x = resolution_mpc / lr_res_x;
      for(size_t i = 0 ; i < Noriginal[0] ; ++i ){
        size_t ii = i * res_tmp_x;
        //size_t ii = i / desample;
        
        if( ii >= nx) break;
        
        tmp = long_range_map.surface_density[ ii + kjj ] += tmp_map.surface_density[ i + kj ];

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
    
    padd_lr = 1 + 2 * border_width_pix * resolution_mpc / long_range_map.x_range();
    padd_lr = MAX<double>(padd_lr,2);
  }  // long rang map deosn't exist

  /// set prefix to subfield maps
  size_t lastindex = fitsfile.find_last_of(".");
  
  subfield_filename = fitsfile.substr(0, lastindex) + "_sr";
  if(dir_scratch == ""){
    subfield_filename = dir_data + subfield_filename;
  }else{
    subfield_filename = dir_scratch + subfield_filename;
  }
  
  short_range_maps.resize(number_of_subfields);

  size_t nx_submap = (size_t)(Noriginal[0] * 1.0 / npanes + 0.5);
  size_t ny_submap = nx_submap;

  if( nx_submap > Noriginal[0]
     || ny_submap > Noriginal[1] ){
    std::cerr << "LensHaloMap : sub map is too large" << std::endl;
    throw std::invalid_argument("out of bounds");
  }

   for(LensMap &map : short_range_maps){
     map.nx = nx_submap;
     map.ny = ny_submap;
   }
   
   { /// set up fftw plans for short range maps
     std::lock_guard<std::mutex> lock(mutex_multimap);
     if(!plans_set_up){
       nx_sub = nx_submap;
       ny_sub = ny_submap;
       nx_sub_extended = nx_sub  + 2*border_width_pix;
       ny_sub_extended = ny_sub  + 2*border_width_pix;

       nx_long = long_range_map.nx;
       ny_long = long_range_map.ny;

       plan_short_range.initialize(nx_sub_extended,nx_sub_extended);
       long_range_map.make_fftw_plans(plan_long_range,padd_lr);

       plans_set_up = true;
     }else{
       // check to make sure the maps' dimension are the same as the static variables
       if( !(nx_sub == nx_submap)*(ny_sub == ny_submap)*
          (nx_sub_extended == short_range_maps[0].nx)*
          (ny_sub_extended == short_range_maps[0].ny)*
          (nx_long == long_range_map.nx)*
          (ny_long == long_range_map.ny)
          ){
         std::cerr << "All instances of LensHaloMultiMap must of the same dimensions" << std::endl;
         throw std::invalid_argument("wrong size");
       }
     }
   }
  
 
  if(!long_range_file_exists){
    long_range_map.ProcessFFTs<WLR>(padd_lr,wlr,plan_long_range);
    
    std::lock_guard<std::mutex> lock(mutex_multimap);
    long_range_map.write("!" + lr_file);
      
    CPFITS_WRITE tmp_cpfits(lr_file,true);
 
    // log parameters so that they can be checked later
    // *****************************************************
    tmp_cpfits.writeKey("ave_ang_sd", ave_ang_sd,"average angulare density");
    tmp_cpfits.writeKey("redshift",redshift,"");
    tmp_cpfits.writeKey("mass_unit",mass_unit,"");
    tmp_cpfits.writeKey("cosmology",to_string(c.ParamSet()),"cosmology");
    tmp_cpfits.writeKey("rs2",rs2,"");
    tmp_cpfits.writeKey("border_width",border_width_pix,"");
    // *****************************************************
                        
  }
};


/*
void LensHaloMultiMap::push_back_submapPhys(Point_2d ll,Point_2d ur){
  if(single_grid) return;
  
  std::vector<long> lower_left(2);
  std::vector<long> upper_right(2);
  
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

  push_back_submap(lower_left,upper_right);
}
*/
void LensHaloMultiMap::resetsubmapPhys(int i
                                       ,Point_2d ll
                                       ){
  
  std::vector<long> lower_left_pix(2);
  
  ll = (ll - long_range_map.lowerleft) / resolution_mpc;
  
  lower_left_pix[0] = floor(ll[0]);
  lower_left_pix[1] = floor(ll[1]);
  
  resetsubmap(i,lower_left_pix);
}

void LensHaloMultiMap::resetsubmap(
                              int i
                              ,const std::vector<long> &lower_left_pix
                                        ){
  if(i >= short_range_maps.size()){
    std::cerr << "Short range map has not been created yet." << std::endl;
    throw std::invalid_argument("out of range");
  }
                                          
  setsubmap(short_range_maps[i],lower_left_pix);
}

void LensHaloMultiMap::setsubmap(LensMap &short_range_map
                                   ,const std::vector<long> &lower_left
                                 ){

  std::vector<long> upper_right =  {lower_left[0] + (long)(nx_sub)-1
      ,lower_left[1] + (long)(ny_sub)-1};
    
  if( (upper_right[0] < 0) || (upper_right[1] < 0) || (lower_left[0] >= (long)(Noriginal[0]) )
      || (lower_left[1] >= (long)(Noriginal[1]))  ){

    std::cerr << "LensHaloMap : sub map is out of bounds" << std::endl;
    std::cerr << "  lower left " << lower_left[0] << " " << lower_left[1]
              << " upper_right " << upper_right[0] << " " << upper_right[1] << std::endl;
    std::cerr << " N " << Noriginal[0] << " " << Noriginal[1] << std::endl;
  }
  
  /// construct file name and have it printed
  std::string sr_file = subfield_filename + to_string(lower_left[0]) + "-"
    + to_string(lower_left[1]) + "-" + to_string(upper_right[0]) + "-"
    + to_string(upper_right[1]) + "_sr.fits";
  bool short_range_file_exists = Utilities::IO::file_exists(sr_file);

  if(short_range_file_exists && write_shorts){
    short_range_map.Myread( sr_file );
    
    short_range_map.lowerleft = short_range_map.upperright = long_range_map.lowerleft;
    short_range_map.lowerleft[0] += lower_left[0]*resolution_mpc;
    short_range_map.lowerleft[1] += lower_left[1]*resolution_mpc;
    
    short_range_map.upperright[0] += (upper_right[0] + 1)*resolution_mpc;
    short_range_map.upperright[1] += (upper_right[1] + 1)*resolution_mpc;
    
    short_range_map.center = (short_range_map.lowerleft + short_range_map.upperright)/2;
    short_range_map.boxlMpc = short_range_map.nx*resolution_mpc;
    
    if( !(nx_sub == short_range_maps[0].nx)*(ny_sub == short_range_maps[0].ny)
       ){
      std::cerr << "All instances of LensHaloMultiMap must of the same dimensions" << std::endl;
      throw std::invalid_argument("wrong size");
    }

  }else{
  
    std::vector<long> first(2);
    std::vector<long> last(2);
  
    first[0] = lower_left[0] - border_width_pix + 1;
    first[1] = lower_left[1] - border_width_pix + 1;

    //last[0] = upper_right[0] + border_width + 1;
    //last[1] = upper_right[1] + border_width + 1;

    last[0] = first[0] + nx_sub_extended - 1;
    last[1] = first[1] + ny_sub_extended - 1;

    double area = resolution_mpc * resolution_mpc / mass_unit;
    LensMap map;
    map.nx = nx_sub_extended;
    map.ny = ny_sub_extended;

    if( (first[0] > 1)*(first[1] > 1)*(last[0] <= Noriginal[0])*(last[1] <= Noriginal[1]) ){
  
      map.surface_density.resize( map.nx * map.ny );
      cpfits.read_subset(&map.surface_density[0]
                         ,first.data(),last.data() );
      
      // convert to relative surface angular surface density
      for(auto &p : map.surface_density){
        p /= area;
        p -= ave_ang_sd;
      }
      
      map.boxlMpc = map.nx * resolution_mpc;

    }else{
    
      // case where subfield overlaps edge
      size_t nx_big = map.nx;
      size_t ny_big = map.ny;
    
      std::vector<long> first_sub(2);
      std::vector<long> last_sub(2);

      map.surface_density.resize(nx_big * ny_big,0);
    
      std::valarray<float> v;
      first_sub[0] = MAX<long>(first[0],1);
      first_sub[1] = MAX<long>(first[1],1);
      last_sub[0] = MIN<long>(last[0],Noriginal[0]);
      last_sub[1] = MIN<long>(last[1],Noriginal[1]);
      
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
      map.boxlMpc = map.nx * resolution_mpc;
    }

    if(wsr_smooth.r_sm2 <= 0){
      map.ProcessFFTs<WSR>(wsr,plan_short_range);
    }else{
      map.ProcessFFTs<WSRSMOOTH>(wsr_smooth,plan_short_range);
    }
 
    //map.PreProcessFFTWMap(wsr,mutex_multimap);

    // cut off bounders
  
    assert(short_range_map.nx == nx_sub);
    assert(short_range_map.ny == ny_sub);
    size_t nx = nx_sub;
    size_t ny = ny_sub;
  
    short_range_map.surface_density.resize(nx*ny);
    short_range_map.alpha1_bar.resize(nx*ny);
    short_range_map.alpha2_bar.resize(nx*ny);
    short_range_map.gamma1_bar.resize(nx*ny);
    short_range_map.gamma2_bar.resize(nx*ny);

    long jj = border_width_pix;
    for(long j=0 ; j < ny ; ++j,++jj){
      long kjj = jj*map.nx;
      long kj = nx*j;
      long ii = border_width_pix;
      for(long i=0 ; i < nx ; ++i,++ii){
        short_range_map.surface_density[i + kj] = map.surface_density[ii + kjj];
        short_range_map.alpha1_bar[i + kj] = map.alpha1_bar[ii + kjj];
        short_range_map.alpha2_bar[i + kj] = map.alpha2_bar[ii + kjj];
        short_range_map.gamma1_bar[i + kj] = map.gamma1_bar[ii + kjj];
        short_range_map.gamma2_bar[i + kj] = map.gamma2_bar[ii + kjj];
      }
    }

    short_range_map.lowerleft = short_range_map.upperright = long_range_map.lowerleft;
    short_range_map.lowerleft[0] += lower_left[0]*resolution_mpc;
    short_range_map.lowerleft[1] += lower_left[1]*resolution_mpc;
  
    short_range_map.upperright[0] += (upper_right[0] + 1)*resolution_mpc;
    short_range_map.upperright[1] += (upper_right[1] + 1)*resolution_mpc;
  
    short_range_map.center = (short_range_map.lowerleft + short_range_map.upperright)/2;
    short_range_map.boxlMpc = short_range_map.nx*resolution_mpc;
 
    if(write_shorts) short_range_map.write("!" + sr_file);
  }
}

bool LensMap::evaluate(const PosType *x,KappaType &sigma,KappaType *gamma,PosType *alpha) {

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
  
  long_range_map.evaluate(xx,*kappa,gamma,alpha);
  
  for(auto &smap : short_range_maps){

    if(( xx[0] >= smap.lowerleft[0] )*( xx[0] <= smap.upperright[0] )
      *( xx[1] >= smap.lowerleft[1] )*( xx[1] <= smap.upperright[1] )
       ){

      KappaType t_kappa,t_gamma[3];
      PosType t_alpha[2];

      smap.evaluate(xx,t_kappa,t_gamma,t_alpha);
    
      alpha[0] += t_alpha[0];
      alpha[1] += t_alpha[1];
      gamma[0] += t_gamma[0];
      gamma[1] += t_gamma[1];
      *kappa += t_kappa;
    }else{
      std::cerr << "ray is off pane!" << std::endl;
      std::cerr << "    rays : " << xx[0] << "  " << xx[1] << std::endl;
      std::cerr << "   range : " << smap.lowerleft[0] << "  " << smap.upperright[0] << std::endl;
      std::cerr << "           " << smap.lowerleft[1] << "  " << smap.upperright[1] << std::endl;

      throw std::runtime_error("out of bounds");
    }
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
    << " CD1_1 - angular resolution in degrees must exit" << std::endl;
    
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
  
  int err = 0;
  {
    /* these are always present in each*/
    float res;
    err = cpfits.readKey ("CD1_1",angular_pixel_size);  // angular resolution degrees
  
    angular_pixel_size *= degreesTOradians;
    res = angular_pixel_size*angDist;
    boxlMpc = res * nx;
  }
  if(err){
    
    std::cerr << "LensMap fits map must have header keywords:" << std::endl
    << " CD1_1 - angular resolution in degrees" << std::endl;
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
    case LensingVariable::KAPPA:
      cpfits.write_image(surface_density,naxex);
      break;
    case LensingVariable::GAMMA1:
      cpfits.write_image(gamma1_bar,naxex);
      break;
    case LensingVariable::GAMMA2:
      cpfits.write_image(gamma2_bar,naxex);
      break;
    case LensingVariable::GAMMA:
    {
      std::valarray<float>  gamma =  sqrt( gamma1_bar*gamma1_bar + gamma2_bar*gamma2_bar );
      cpfits.write_image(gamma,naxex);
    }
      break;
    case LensingVariable::ALPHA1:
      cpfits.write_image(alpha1_bar,naxex);
      cpfits.writeKey("UNITS","radians","units of values");
     break;
    case LensingVariable::ALPHA2:
      cpfits.write_image(alpha2_bar,naxex);
      cpfits.writeKey("UNITS","radians","units of values");
      break;
    default:
      break;
  }
  
  cpfits.writeKey("CD1_1",-angular_pixel_size /degreesTOradians,"pixel size in degrees");
  cpfits.writeKey("CD2_2",angular_pixel_size /degreesTOradians,"pixel size in degrees");
  cpfits.writeKey("SIDEL1",boxlMpc,"x range in radians");
  cpfits.writeKey("SIDEL2",y_range(),"y range in radians");
  cpfits.writeKey("CENTER_X",center[0],"center of field x");
  cpfits.writeKey("CENTER_Y",center[1],"center of field y");

}

