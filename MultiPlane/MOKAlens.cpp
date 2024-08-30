/*
 * MOKAlens.cpp
 *
 *  Created on: Jun 8, 2012
 *      Author: mpetkova
 */


#include "slsimlib.h"

using namespace std;

/*
 * \brief center of mass of a map using the moving centre
 *
void cmass(int n, std::valarray<double> map, std:: vector<double> x, double &xcm, double &ycm){
  double shrink = 1.02;
  // 0.05% of the total number of pixels
  int nstop = int(0.05*double(n)*double(n)/100.);
  if(nstop<36) nstop = 36;
  std:: cout << "nstop = " << nstop << std:: endl;
  xcm = 0.;
  ycm = 0.;
  double tot = 0;
  //
  double rsel = x[n-1] - x[0];
  double xmin = x[0];
  double drpix = rsel/n;
  rsel/=2.;
  // I will set in the centre with t width equal to 128 pixels
  rsel = drpix*64.;
  //
  int nsel = n*n;
  int i,j;
  double dist;
  while(nsel>nstop){
    double nxcm = 0;
    double nycm = 0;
    double xi,yi;
    tot = 0;
    nsel = 0;
    for(i=0;i<n;i++) for(j=0;j<n;j++){
      xi = (xmin+(drpix*0.5)+i*drpix);
      yi = (xmin+(drpix*0.5)+j*drpix);
      dist = sqrt(pow(xi-xcm,2)+ pow(yi-ycm,2));
      if(dist<=rsel){
        nxcm+=map[i+n*j]*xi;
        nycm+=map[i+n*j]*yi;
        tot+=map[i+n*j];
        nsel++;
      }
    }
    nxcm/=tot;
    nycm/=tot;
    rsel = rsel / shrink;
    // center of mass
    xcm = nxcm;
    ycm = nycm;
  }
}*/

/**
 * \brief loads a mass map from a given filename
 */

LensHaloMassMap::LensHaloMassMap(const std::string& filename, PixelMapType my_maptype,int pixel_map_zeropad
                                 ,bool my_zeromean, const COSMOLOGY& lenscosmo)
: LensHalo(),
MOKA_input_file(filename), flag_MOKA_analyze(0), flag_background_field(0),
maptype(my_maptype), cosmo(lenscosmo),zerosize(pixel_map_zeropad),zeromean(my_zeromean)
{
  initMap();
  if(zerosize < 1.0 ){
    std::cerr << "pixel_map_zeropad in LensHaloMap cosntructor must be >= 1" << std::endl;
  }
  // set redshift to value from map
  setZlens(map.zlens,lenscosmo);
}

/** \brief Create a LensHalo from a PixelMap representing the mass.
 
 This is especially useful for representing the visible stars from an image in
 the lens model.
 
 constructor for making lens halo directly from a mass map
 
 The pixel map should be defined with physical Mpc dimensions.
 
 */

LensHaloMassMap::LensHaloMassMap(
                                 const PixelMap<double> &MassMap   /// mass map
                                 ,double massconvertion    /// convertion factor from pixel units to solar masses
                                 ,double redshift          /// redshift of lens
                                 ,int pixel_map_zeropad    /// factor by which to zero pad in FFTs, ex. 1 is no padding, 2 FFT grid is twice as big as original map
                                 ,bool my_zeromean         /// if true, subtracts average density
                                 ,const COSMOLOGY& lenscosmo  /// cosmology
)
:LensHalo()
, flag_MOKA_analyze(0),flag_background_field(0),maptype(PixelMapType::pix_map),cosmo(lenscosmo),zerosize(pixel_map_zeropad),zeromean(my_zeromean)
{
  rscale = 1.0;

  setMap<double>(MassMap,massconvertion,redshift);
  
  LensHalo::setTheta(MassMap.getCenter()[0],MassMap.getCenter()[1]);
  
  setZlensDist(map.zlens,cosmo);
}

LensHaloMassMap::LensHaloMassMap(
                                 const PixelMap<float> &MassMap   /// mass map
                                 ,double massconvertion    /// convertion factor from pixel units to solar masses
                                 ,double redshift          /// redshift of lens
                                 ,int pixel_map_zeropad    /// factor by which to zero pad in FFTs, ex. 1 is no padding, 2 FFT grid is twice as big as original map
                                 ,bool my_zeromean         /// if true, subtracts average density
                                 ,const COSMOLOGY& lenscosmo  /// cosmology
)
:LensHalo()
, flag_MOKA_analyze(0),flag_background_field(0),maptype(PixelMapType::pix_map),cosmo(lenscosmo),zerosize(pixel_map_zeropad),zeromean(my_zeromean)
{
  rscale = 1.0;

  setMap<float>(MassMap,massconvertion,redshift);
  
  LensHalo::setTheta(MassMap.getCenter()[0],MassMap.getCenter()[1]);
  
  setZlensDist(map.zlens,cosmo);
}

LensHaloMassMap::LensHaloMassMap(
                double mass          /// total mass on sheet
                ,Point_2d center
                ,Point_2d range      /// range in radians
                ,double resolution   /// resolution in radians
                ,int zeropadding   
                ,double redshift
                ,const COSMOLOGY &cosmo
                )
:LensHalo(),flag_MOKA_analyze(0),flag_background_field(0),maptype(PixelMapType::pix_map),cosmo(cosmo)
,zerosize(zeropadding),zeromean(false)
{
  rscale = 1.0;

  double d = cosmo.angDist(redshift);
  
  size_t Nx = abs(range[0])/resolution;
  size_t Ny = abs(range[1])/resolution;

  PixelMap<double> mass_map(center.x,Nx,Ny,resolution * d);  // Is this right ????

  double density = mass / Nx / Ny;
  
  for(double &a : mass_map.data()){
    a = density;
  }
  
  setMap<double>(mass_map,1,redshift);
  LensHalo::setTheta(mass_map.getCenter()[0],mass_map.getCenter()[1]);
  
  setZlensDist(map.zlens,cosmo);
}
  

/*
LensHaloMassMap::LensHaloMassMap(
                                 PixelMap &my_map        /// map of mass
                                 ,double massconvertion  /// factor that converts the units of my_map to solar masses
                                 ,double zlens           /// redshift of lens
                                 ,double zsource         /// redshit of source
                                 ,int pixel_map_zeropad  /// factor by which grid is expanded with padding when doing FFTs
                                 ,const COSMOLOGY& lenscosmo  /// cosmology
):
LensHalo(),MOKA_input_file(""),maptype(pix_map),cosmo(lenscosmo),zerosize(pixel_map_zeropad),zeromean(false)
{
  
  rscale = 1.0;
  Dist = lenscosmo.angDist(zlens);

  map = new MOKAmap();
  
  if(std::numeric_limits<float>::has_infinity)
    Rmax = std::numeric_limits<float>::infinity();
  else
    Rmax = std::numeric_limits<float>::max();
  
  map.nx = my_map.getNx();
  map.ny = my_map.getNy();
  std::cout << "nx           ny " << std::endl;
  std::cout << map.nx << "   " << map.ny << std::endl;
  
  std::size_t size = map.nx*map.ny;
 
  map.surface_density.resize(size);
  map.alpha1.resize(size,0.0);
  map.alpha2.resize(size,0.0);
  map.gamma1.resize(size,0.0);
  map.gamma2.resize(size,0.0);
  map.gamma3.resize(size,0.0);
  map.phi.resize(size,0.0);
  
  map.boxlarcsec = my_map.getRangeX()/arcsecTOradians;
  map.boxlrad = my_map.getRangeX();
  map.Dlens = Dist;
  map.boxlMpc = map.boxlrad/Dist; // physical
  
  // convertion to solar masses per Mpc^2
  double convert = massconvertion/(my_map.getResolution()*my_map.getResolution()*map.Dlens*map.Dlens);
  
  size_t Npixels = map.nx*map.ny;
  LensHalo::mass = 0.0;
  for(size_t i=0 ; i<Npixels ;++i){
    map.surface_density[i] = my_map(i)*convert;
    LensHalo::mass += my_map(i)*massconvertion;
  }
  
  map.zsource = zsource;
  map.zlens = zlens;
  
  map.center[0] = map.center[1] = 0.0;
  map.boxlrad = map.boxlarcsec*PI/180/3600.;

  PreProcessFFTWMap();
  
  assert(map.nx*map.ny == map.surface_density.size());
 
  LensHalo::setTheta(my_map.getCenter()[0],my_map.getCenter()[1]);
  
  setZlens(zlens);
}
*/
/*
 * \brief allocates and reads the MOKA map in
 *
 *  In the future this could be used to read in individual PixelDMaps or other types of maps if the type were specified in the paramfile.
 */
//LensHaloMassMap::LensHaloMassMap(InputParams& params, COSMOLOGY& lenscosmo)
//: LensHalo(), maptype(moka), cosmo(lenscosmo)
//{
//  // read in parameters
//  assignParams(params);
//  
//  // initialize MOKA map
//  initMap();
//  
//  // set redshift if necessary
//  if(LensHalo::getZlens() == -1)
//    setZlens(map.zlens,lenscosmo);
//}

LensHaloMassMap::~LensHaloMassMap()
{
}

void LensHaloMassMap::initMap()
{
  
  if(!(maptype == PixelMapType::pix_map || maptype == PixelMapType::moka)){
    ERROR_MESSAGE();
    throw runtime_error("Does not recognize input lens map type");
  }
  
#ifndef ENABLE_FITS
  std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
  exit(1);
#endif
  
  if(std::numeric_limits<float>::has_infinity)
    Rmax = std::numeric_limits<float>::infinity();
  else
    Rmax = std::numeric_limits<float>::max();
  
  getDims();
  
  readMap();
  
  if(flag_background_field == 1)
  {
    map.surface_density = 0;
    map.alpha1_bar = 0;
    map.alpha2_bar = 0;
    map.gamma1_bar = 0;
    map.gamma2_bar = 0;
    map.phi_bar = 0;
  }
  
  map.center[0] = map.center[1] = 0.0;
  map.boxlrad = map.boxlarcsec*PI/180/3600.;
  
  if(maptype == PixelMapType::moka){
    
    /// converts to the code units
    std::cout << "converting the units of the MOKA map" << std::endl;
    
    double fac = map.DS/map.DLS/map.Dlens*map.h/(4*PI*Grav);
    
    map.surface_density *= fac;
    map.gamma1_bar *= fac;
    map.gamma2_bar *= fac;
    map.phi_bar *= fac;

    fac = 1/(4*PI*Grav);
    
    map.alpha1_bar *= fac;
    map.alpha2_bar *= fac;
   
    checkCosmology();
  }else{
    // simulation case
    // not needed for now does all in MOKAfits.cpp
    // convertmap(map,maptype);
  }
}

/**
 * \brief reads in the fits file for the MOKA or mass map and saves it in the structure map
 */
template<typename T>
void LensHaloMassMap::setMap(
                                 const PixelMap<T> &inputmap  // mass map
                                 ,double massconvertion    // convertion factor from pixel units to solar masses
                                 ,double z
                                 ){
  
  // must be a square map
  // ???? assert(inputmap.getNx() == inputmap.getNy());
  
  if(std::numeric_limits<float>::has_infinity)
    Rmax = std::numeric_limits<float>::infinity();
  else
    Rmax = std::numeric_limits<float>::max();
  
  map.nx = inputmap.getNx();
  map.ny = inputmap.getNy();
  map.center[0] = map.center[1] = 0.0;
  
  std::size_t size = map.nx*map.ny;
  
  map.surface_density.resize(size);
  
  map.zlens = z;
  
  assert(map.nx !=0);
  // keep it like it is, even if is a rectangle
  
  map.Dlens = cosmo.angDist(0.,map.zlens);  // physical
  map.boxlrad = inputmap.getRangeX();
  map.boxlarcsec = inputmap.getRangeX()/arcsecTOradians;
  map.boxlMpc = inputmap.getRangeX()*map.Dlens;
  
  double pixelarea = inputmap.getResolution()*map.Dlens;
  pixelarea *= pixelarea;
  
  for(size_t i=0;i<size;++i){
    //assert(!isnan(inputmap(i)));
    map.surface_density[i] = massconvertion*inputmap(i)/pixelarea;
  }
  
  if(zeromean){
    double avkappa = 0;
    
    for(size_t i=0;i<size;i++){
      avkappa += map.surface_density[i];
    }
    avkappa /= size;
    
    for(size_t i=0;i<size;i++){
      map.surface_density[i] -= avkappa;
    }
  }
  
  // kappa is not divided by the critical surface density
  // they don't need to be preprocessed by fact
  // create alpha and gamma arrays by FFT
  // valid only to force the map to be square map.nx = map.ny = npixels;
  
  //std:: cout << "  preProcessing Map " << std:: endl;
  map.PreProcessFFTWMap(zerosize);
}

/** \brief checks the cosmology against the MOKA map parameters
 */
/// checks that cosmology in the header of the input fits map is the same as the one set
void LensHaloMassMap::checkCosmology()
{
  if(cosmo.getOmega_matter() == map.omegam)
    std::cerr << "LensHaloMassMap: Omega_matter " << cosmo.getOmega_matter() << " (cosmology) != " << map.omegam << " (MOKA)" << std::endl;
  if(cosmo.getOmega_lambda() == map.omegal)
    std::cerr << "LensHaloMassMap: Omega_lambda " << cosmo.getOmega_lambda() << " (cosmology) != " << map.omegal << " (MOKA)" << std::endl;
  if(cosmo.gethubble() == map.h)
    std::cerr << "LensHaloMassMap: hubble " << cosmo.gethubble() << " (cosmology) != " << map.h << " (MOKA)" << std::endl;
}

/**
 * Sets many parameters within the MOKA lens model
 */

void LensHaloMassMap::assignParams(InputParams& params)
{
  PosType tmp;
  if(!params.get("z_lens", tmp)){
    LensHalo::setZlens(-1,cosmo); // set to -1 so that it will be set to the MOKA map value
  }else{
    LensHalo::setZlens(tmp,cosmo);
  }
  if(!params.get("MOKA_input_file", MOKA_input_file))
  {
    ERROR_MESSAGE();
    std::cout << "parameter MOKA_input_file needs to be set in parameter file " << params.filename() << std::endl;
    exit(0);
  }
  
  if(!params.get("MOKA_background_field", flag_background_field))
    flag_background_field = 0;
  
  if(!params.get("MOKA_analyze",flag_MOKA_analyze))
    flag_MOKA_analyze = 0;
  
  zeromean = true;
  zerosize = 4;  /// not really used
  maptype = PixelMapType::moka;
}

/** 
 *
 * \brief Routine for obtaining the deflection and other lensing quantities for
 * a LensHaloMassMap for just one ray!!
 *
 */
void LensHaloMassMap::force_halo(double *alpha
                                 ,KappaType *kappa
                                 ,KappaType *gamma
                                 ,KappaType *phi
                                 ,double const *xx      /// position in physical Mpc
                                 ,bool subtract_point
                                 ,PosType screening
                                 )
{
  // interpolate from the maps
  
  Utilities::Interpolator<valarray<double> > interp(xx,map.nx,map.boxlMpc,map.ny
          ,map.ny*map.boxlMpc/map.nx,map.center.x);

  // ??? assert(map.nx == map.ny);
  
  alpha[0] += interp.interpolate(map.alpha1_bar);
  alpha[1] += interp.interpolate(map.alpha2_bar);
  gamma[0] += interp.interpolate(map.gamma1_bar);
  gamma[1] += interp.interpolate(map.gamma2_bar);
  gamma[2] += 0.0;

  *kappa += interp.interpolate(map.surface_density);
  *phi += interp.interpolate(map.phi_bar);
  
  return;
}

