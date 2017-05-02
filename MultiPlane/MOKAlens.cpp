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
 */
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
}

/**
 * \brief loads a mass map from a given filename
 */

LensHaloMassMap::LensHaloMassMap(const std::string& filename, PixelMapType my_maptype,int pixel_map_zeropad,bool my_zeromean, const COSMOLOGY& lenscosmo)
: LensHalo(),
MOKA_input_file(filename), flag_MOKA_analyze(0), flag_background_field(0),
maptype(my_maptype), cosmo(lenscosmo),zerosize(pixel_map_zeropad),zeromean(my_zeromean)
{
  initMap();
 	
  // set redshift to value from map
  setZlens(map->zlens);
}

/** \brief Create a LensHalo from a PixelMap representing the mass.
 
 This is especially useful for representing the visible stars from an image in
 the lens model.
 */
/// constructor for making lens halo directly from a mass map
LensHaloMassMap::LensHaloMassMap(
                                 const PixelMap &MassMap   /// mass map in solar mass units
                                 ,double massconvertion    /// convertion factor from pixel units to solar masses
                                 ,double redshift          /// redshift of lens
                                 ,int pixel_map_zeropad    /// factor by which to zero pad in FFTs, ex. 4
                                 ,bool my_zeromean         /// if true, subtracts average density
                                 ,const COSMOLOGY& lenscosmo  /// cosmology
)
:LensHalo()
, flag_MOKA_analyze(0), flag_background_field(0),maptype(pix_map),cosmo(lenscosmo),zerosize(pixel_map_zeropad),zeromean(my_zeromean)
{
  rscale = 1.0;

  setMap(MassMap,massconvertion,redshift);
  
  LensHalo::setTheta(MassMap.getCenter()[0],MassMap.getCenter()[1]);
  
  setZlensDist(map->zlens,cosmo);
  //setZlens(redshift);
  // set redshift to value from map
  //setZlens(map->zlens);
  
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
  
  map->nx = my_map.getNx();
  map->ny = my_map.getNy();
  std::cout << "nx           ny " << std::endl;
  std::cout << map->nx << "   " << map->ny << std::endl;
  
  std::size_t size = map->nx*map->ny;
  
  map->convergence.resize(size);
  map->alpha1.resize(size,0.0);
  map->alpha2.resize(size,0.0);
  map->gamma1.resize(size,0.0);
  map->gamma2.resize(size,0.0);
  map->gamma3.resize(size,0.0);
  map->phi.resize(size,0.0);
  
  
  map->boxlarcsec = my_map.getRangeX()/arcsecTOradians;
  map->boxlrad = my_map.getRangeX();
  map->Dlens = Dist;
  map->boxlMpc = map->boxlrad/Dist; // physical
  
  // convertion to solar masses per Mpc^2
  double convert = massconvertion/(my_map.getResolution()*my_map.getResolution()*map->Dlens*map->Dlens);
  
  size_t Npixels = map->nx*map->ny;
  LensHalo::mass = 0.0;
  for(size_t i=0 ; i<Npixels ;++i){
    map->convergence[i] = my_map(i)*convert;
    LensHalo::mass += my_map(i)*massconvertion;
  }
  
  map->zsource = zsource;
  map->zlens = zlens;
  
  map->center[0] = map->center[1] = 0.0;
  map->boxlrad = map->boxlarcsec*pi/180/3600.;

  PreProcessFFTWMap();
  
  assert(map->nx*map->ny == map->convergence.size());
  
  LensHalo::setTheta(my_map.getCenter()[0],my_map.getCenter()[1]);
  
  setZlens(zlens);
}
*/
/**
 * \brief allocates and reads the MOKA map in
 *
 *  In the future this could be used to read in individual PixelDMaps or other types of maps if the type were specified in the paramfile.
 */
LensHaloMassMap::LensHaloMassMap(InputParams& params, const COSMOLOGY& lenscosmo)
: LensHalo(), maptype(moka), cosmo(lenscosmo)
{
  // read in parameters
  assignParams(params);
  
  // initialize MOKA map
  initMap();
  
  // set redshift if necessary
  if(LensHalo::getZlens() == -1)
    setZlens(map->zlens);
}

LensHaloMassMap::~LensHaloMassMap()
{
  delete map;
}

void LensHaloMassMap::initMap()
{
  
  if(!(maptype == pix_map || maptype == moka)){
    ERROR_MESSAGE();
    throw runtime_error("Does not recognize input lens map type");
  }
  
#ifndef ENABLE_FITS
  std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
  exit(1);
#endif
  
  map = new MOKAmap();
  
  if(std::numeric_limits<float>::has_infinity)
    Rmax = std::numeric_limits<float>::infinity();
  else
    Rmax = std::numeric_limits<float>::max();
  
  getDims();
  
  std::size_t size = map->nx*map->ny;
  
  map->convergence.resize(size);
  map->alpha1.resize(size);
  map->alpha2.resize(size);
  map->gamma1.resize(size);
  map->gamma2.resize(size);
  map->gamma3.resize(size);
  map->phi.resize(size);
  
  readMap();
  
  if(flag_background_field == 1)
  {
    map->convergence = 0;
    map->alpha1 = 0;
    map->alpha2 = 0;
    map->gamma1 = 0;
    map->gamma2 = 0;
    map->phi = 0;
  }
  
  map->center[0] = map->center[1] = 0.0;
  map->boxlrad = map->boxlarcsec*pi/180/3600.;
  
  if(maptype == moka){
    
    Utilities::fill_linear(map->x, map->nx, -0.5*map->boxlMpc, 0.5*map->boxlMpc); // physical Mpc/h
    
    /// converts to the code units
    std::cout << "converting the units of the MOKA map" << std::endl;
    
    double fac = map->DS/map->DLS/map->Dlens*map->h/(4*pi*Grav);
    
    map->convergence *= fac;
    map->gamma1 *= fac;
    map->gamma2 *= fac;
    map->phi *= fac;
    
    fac = 1/(4*pi*Grav);
    
    map->alpha1 *= fac;
    map->alpha2 *= fac;
    
    checkCosmology();
  }else{
    // simulation case
    // not needed for now does all in MOKAfits.cpp
    // convertmap(map,maptype);
  }
}

void LensHaloMassMap::convertmap(MOKAmap *map,PixelMapType maptype){
  
  // TODO: convert units
  throw std::runtime_error("needs to be finished");
  
  map->center[0] *= map->Dlens;//(1+map->zlens);
  map->center[1] *= map->Dlens;//(1+map->zlens);
  
  //float pixLMpc = map->boxlMpc/map->nx;   // TODO: What if it isn't square?
  //float fac = 1.e+10/pixLMpc/pixLMpc/cosmo->gethubble();
  
  //map->convergence *= fac;
  //map->gamma1 *= fac;
  //map->gamma2 *= fac;
  
  // TODO: Need to check this
  //map->alpha1 *= fac*pixLMpc;
  //map->alpha2 *= fac*pixLMpc;
  
}

/**
 * \brief reads in the fits file for the MOKA or mass map and saves it in the structure map
 */
void LensHaloMassMap::setMap(
                                 const PixelMap &inputmap  // mass map
                                 ,double massconvertion    // convertion factor from pixel units to solar masses
                                 ,double z
                                 ){
  
  // must be a square map
  assert(inputmap.getNx() == inputmap.getNy());
  
  map = new MOKAmap();
  
  if(std::numeric_limits<float>::has_infinity)
    Rmax = std::numeric_limits<float>::infinity();
  else
    Rmax = std::numeric_limits<float>::max();
  
  map->nx = map->ny = inputmap.getNx();
  map->center[0] = map->center[1] = 0.0;
  
  std::size_t size = map->nx*map->ny;
  
  map->convergence.resize(size);
  map->alpha1.resize(size);
  map->alpha2.resize(size);
  map->gamma1.resize(size);
  map->gamma2.resize(size);
  map->gamma3.resize(size);
  map->phi.resize(size);
  
  map->zlens = z;
  
  assert(map->nx !=0);
  // keep it like it is, even if is a rectangle
  
  map->Dlens = cosmo.angDist(0.,map->zlens);  // physical
  map->boxlrad = inputmap.getRangeX();
  map->boxlarcsec = inputmap.getRangeX()/arcsecTOradians;
  map->boxlMpc = inputmap.getRangeX()/map->Dlens;
  
  double pixelarea = inputmap.getResolution()*map->Dlens;
  pixelarea *= pixelarea;
  
  for(size_t i=0;i<size;++i){
    map->convergence[i] = massconvertion*inputmap(i)/pixelarea;
  }
  
  if(zeromean){
    double avkappa = 0;
    
    for(size_t i=0;i<size;i++){
      avkappa += map->convergence[i];
    }
    avkappa /= size;
    
    for(size_t i=0;i<size;i++){
      map->convergence[i] -= avkappa;
    }
  }
  
  // kappa is not divided by the critical surface density
  // they don't need to be preprocessed by fact
  // create alpha and gamma arrays by FFT
  // valid only to force the map to be square map->nx = map->ny = npixels;
#ifdef ENABLE_FFTW
  std:: cout << "  preProcessing Map " << std:: endl;
  PreProcessFFTWMap();
#else
  std::cout << "Please enable the preprocessor flag ENABLE_FFTW !" << std::endl;
  exit(1);
#endif
  
}

/** \brief checks the cosmology against the MOKA map parameters
 */
/// checks that cosmology in the header of the input fits map is the same as the one set
void LensHaloMassMap::checkCosmology()
{
  if(cosmo.getOmega_matter() == map->omegam)
    std::cerr << "LensHaloMassMap: Omega_matter " << cosmo.getOmega_matter() << " (cosmology) != " << map->omegam << " (MOKA)" << std::endl;
  if(cosmo.getOmega_lambda() == map->omegal)
    std::cerr << "LensHaloMassMap: Omega_lambda " << cosmo.getOmega_lambda() << " (cosmology) != " << map->omegal << " (MOKA)" << std::endl;
  if(cosmo.gethubble() == map->h)
    std::cerr << "LensHaloMassMap: hubble " << cosmo.gethubble() << " (cosmology) != " << map->h << " (MOKA)" << std::endl;
}

/**
 * Sets many parameters within the MOKA lens model
 */

void LensHaloMassMap::assignParams(InputParams& params)
{
  PosType tmp;
  if(!params.get("z_lens", tmp)){
    LensHalo::setZlens(-1); // set to -1 so that it will be set to the MOKA map value
  }else{
    LensHalo::setZlens(tmp);
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
  maptype = moka;
}

/**
 * saves the image, by reading off the values from the image tree
 * and then saving to a fits file and computing the radial profile
 * of the convergence
 */
void LensHaloMassMap::saveImage(GridHndl grid,bool saveprofiles){
  std::stringstream f;
  std::string filename;
  
  if(flag_background_field==1) f << MOKA_input_file << "_only_noise.fits";
  else{
    if(flag_MOKA_analyze == 0) f << MOKA_input_file << "_noisy.fits";
    else f << MOKA_input_file << "_no_noise.fits";
  }
  filename = f.str();
  
  PointList::iterator i_tree_pointlist_current(grid->i_tree->pointlist->Top());
  
  do{
    long index = Utilities::IndexFromPosition((*i_tree_pointlist_current)->x,map->nx,map->boxlrad,map->center);
    if(index > -1){
      map->convergence[index] = (*i_tree_pointlist_current)->kappa;
      map->gamma1[index] = (*i_tree_pointlist_current)->gamma[0];
      map->gamma2[index] = (*i_tree_pointlist_current)->gamma[1];
      map->gamma3[index] = (*i_tree_pointlist_current)->gamma[2];
    }
  }while( (--i_tree_pointlist_current)==true );
  
  writeImage(filename);
  
  if(saveprofiles == true){
    
    std:: cout << " saving profile " << std:: endl;
    double RE3,xxc,yyc;
    saveProfiles(RE3,xxc,yyc);
    estSignLambdas();
    double RE1,RE2;
    EinsteinRadii(RE1,RE2,xxc,yyc);
    std::ostringstream fEinr;
    if(flag_background_field==1) fEinr << MOKA_input_file << "_only_noise_Einstein.radii.dat";
    else{
      if(flag_MOKA_analyze == 0) fEinr << MOKA_input_file << "_noisy_Einstein.radii.dat";
      else fEinr << MOKA_input_file << "_no_noise_Einstein.radii.dat";
    }
    std:: ofstream filoutEinr;
    std:: string filenameEinr = fEinr.str();
    filoutEinr.open(filenameEinr.c_str());
    filoutEinr << "# effective        median      from_profiles" << std:: endl;
    filoutEinr << RE1 << "   " << RE2 << "    " << RE3 << std:: endl;
    filoutEinr.close();
  }
}


/**
 * computing and saving the radial profile of the convergence, reduced tangential and parallel shear and of the shear
 *  */
void LensHaloMassMap::saveProfiles(double &RE3,double &xxc,double &yyc){
  /* measuring the differential and cumulative profile*/
  double xmin = -map->boxlMpc*0.5;
  double xmax =  map->boxlMpc*0.5;
  double drpix = map->boxlMpc/map->nx;
  
  int galaxiesPerBin = 64;
  
  std::valarray<double> pxdist(map->nx*map->ny);
  std::valarray<double> red_sgE(map->nx*map->ny),red_sgB(map->nx*map->ny),sgm(map->nx*map->ny);
  int i, j;
  /*
	  measure the center of mass
	  double xcm=0,ycm=0,tot=0;
	  for(i=0; i<map->nx; i++ ) for(j=0; j<map->ny; j++ ){
	  xcm+=map->convergence[i+map->ny*j]*(xmin+(drpix*0.5)+i*drpix);
	  ycm+=map->convergence[i+map->ny*j]*(xmin+(drpix*0.5)+j*drpix);
	  tot+=map->convergence[i+map->ny*j];
	  }
	  xcm/=tot;
	  ycm/=tot;
	  xxc = xcm;
	  yyc = ycm;
   */
  // moving center
  double xcm,ycm;
  if(flag_background_field==1){
    xcm = 0.;
    ycm = 0.;
  }
  else{
    cmass(map->ny,map->convergence,map->x,xcm,ycm);
  }
  xxc = xcm;
  yyc = ycm;
  int ai = Utilities::locate(map->x,xxc);
  ai = ((ai > 0) ? ai:0);
  ai = ((ai < map->nx-1) ? ai:map->nx-1);
  int bi = Utilities::locate(map->x,yyc);
  bi = ((bi > 0) ? bi:0);
  bi = ((bi < map->nx-1) ? bi:map->nx-1);
  std:: cout << "  ------------ center ------------- " << std:: endl;
  std:: cout << "    " << xxc << "  " << yyc << std:: endl;
  std:: cout << "    " << ai << "  " << bi << std:: endl;
  std:: cout << "  ------------------r ------------- " << std:: endl;
  for(i=0; i<map->nx; i++ ) for(j=0; j<map->ny; j++ ){
    pxdist[i+map->ny*j]= sqrt(pow((xmin+(drpix*0.5)+i*drpix-xcm),2) +
                              pow((xmin+(drpix*0.5)+j*drpix-ycm),2));
    // reduced shear E a B
    double dx=map->x[i];
    double dy=map->x[j];
    double p=atan2( dy, dx ); // check gamma 1 and gamma 2 definition
    red_sgE[i+map->ny*j] = (-map->gamma1[i+map->ny*j]*cos(2*p)-map->gamma2[i+map->ny*j]*sin(2*p))/(1.-map->convergence[i+map->ny*j]);
    red_sgB[i+map->ny*j] = (map->gamma1[i+map->ny*j]*sin(2*p)-map->gamma2[i+map->ny*j]*cos(2*p))/(1.-map->convergence[i+map->ny*j]);
    sgm[i+map->ny*j] = sqrt(pow(map->gamma1[i+map->ny*j],2) + pow(map->gamma2[i+map->ny*j],2));
  }
  
  double dr0 = 8.*(0.5*map->boxlMpc)/(map->nx/2.);
  int nbin = int(xmax/dr0);
  int nbggal=0;
  int ih;
  // check if backgroundgalXarcmin2.d file exist
  std:: string filebgXarcmin2 = "backgroundgalXarcmin2.d";
  std:: ifstream infilebgXarcmin2;
  infilebgXarcmin2.open(filebgXarcmin2.c_str());
  // ... if so it read it!
  if(infilebgXarcmin2.is_open()){
    std:: cout << " I will read  backgroundgalXarcmin2.d file" << std:: endl;
    infilebgXarcmin2 >> nbggal;
    infilebgXarcmin2 >> ih;
    std:: cout << "  and consider " << nbggal <<  " per arcmin2 as background points " << std:: endl;
    std:: cout << " in building the cluster profile" << std:: endl;
    infilebgXarcmin2.close();
  }
  // number of galaxies per arcminsquare
  double dntbggal=double(nbggal)*(map->boxlarcsec*map->boxlarcsec)/3600.;
  int ntbggal=int(dntbggal+0.5);
  std:: vector<int> runi,runj;
  int ibut;
  int lrun;
  if(ntbggal>0){
    std:: cout << " ~ ~ ~ " << ntbggal << " number of bg gal in the field" << std:: endl;
    // fill the vector with random coordinates of background galaxies
    srand(ih+94108);
    for(lrun=0;lrun<ntbggal;lrun++){
      ibut = rand()%map->nx;
      if(ibut>map->nx || ibut<0) std:: cout << "1. ibut out of npix range = " << ibut << std:: endl;
      runi.push_back(ibut);
      ibut = rand()%map->nx;
      if(ibut>map->nx || ibut<0) std:: cout << "2. ibut out of npix range = " << ibut << std:: endl;
      runj.push_back(ibut);
    }
    if(runi.size()>ntbggal || runj.size()>ntbggal){
      std:: cout << " random points in the field out of range " << std:: endl;
      std:: cout << " runi.size() = " << runi.size() << std:: endl;
      std:: cout << " runj.size() = " << runj.size() << std:: endl;
    }
    
    dr0 = map->boxlMpc*sqrt(galaxiesPerBin/ntbggal);
    
    nbin = int(xmax/dr0);
  }
  //
  std:: cout << "   " << std:: endl;
  std:: cout << " nbins = " << nbin << "  dr0 = " << dr0 << std:: endl;
  std:: cout << " ______________________________________________________ " << std:: endl;
  std:: cout << " computing profiles assuming spherical symmetry";
  // - - - - - - - - - - - - - - - - -
  
  // TODO: Carlo:  These are all memory leaks!  They are never deleted!
  double *kprofr = estprof(map->convergence,map->nx,map->ny,pxdist,dr0,xmax,runi,runj,ntbggal);
  double *sigmakprof = estsigmaprof(map->convergence,map->nx,map->ny,pxdist,dr0,xmax,runi,runj,ntbggal,kprofr);
  double *ckprofr = estcprof(map->convergence,map->nx,map->ny,pxdist,dr0,xmax,runi,runj,ntbggal);
  double *sigmackprof = estsigmacprof(map->convergence,map->nx,map->ny,pxdist,dr0,xmax,runi,runj,ntbggal,kprofr);
  double *gamma1profr = estprof(red_sgE,map->nx,map->ny,pxdist,dr0,xmax,runi,runj,ntbggal); // reduced shear
  double *sigmagamma1prof = estsigmaprof(red_sgE,map->nx,map->ny,pxdist,dr0,xmax,runi,runj,ntbggal,gamma1profr);
  double *gamma0profr = estprof(red_sgB,map->nx,map->ny,pxdist,dr0,xmax,runi,runj,ntbggal);
  double *sigmagamma0prof = estsigmaprof(red_sgB,map->nx,map->ny,pxdist,dr0,xmax,runi,runj,ntbggal,gamma0profr);
  double *gamma2profr = estprof(sgm,map->nx,map->ny,pxdist,dr0,xmax,runi,runj,ntbggal);
  double *sigmagamma2prof = estsigmaprof(sgm,map->nx,map->ny,pxdist,dr0,xmax,runi,runj,ntbggal,gamma2profr);
  std::ostringstream fprof;
  
  if(flag_background_field==1) fprof << MOKA_input_file << "_only_noise_MAP_radial_prof.dat";
  else{
    if(flag_MOKA_analyze == 0) fprof << MOKA_input_file << "_noisy_MAP_radial_prof.dat";
    else fprof << MOKA_input_file << "_no_noise_MAP_radial_prof.dat";
  }
  std:: ofstream filoutprof;
  std:: string filenameprof = fprof.str();
  filoutprof.open(filenameprof.c_str());
  filoutprof <<"# r      kappa     sig_k     ckappa     sig_ck    redgE   sig_redgE   redgB   sig_redgE   g   sig_g   theta   Anulus_area" << std:: endl;
  int l;
  std:: vector<double> lrOFr(nbin),lprofFORre(nbin);
  for(l=0;l<nbin;l++){
    double Aanulus = M_PI*((dr0*l+dr0)*(dr0*l+dr0)-(dr0*l)*(dr0*l));
    filoutprof << dr0*l + dr0/2. << "  "
    << kprofr[l] << "  " << sigmakprof[l] << "  "
    << ckprofr[l] << "  " << sigmackprof[l] << "   "
    << gamma1profr[l] << "  " << sigmagamma1prof[l] << "  "
    << gamma0profr[l] << "  " << sigmagamma0prof[l] << "  "
    << gamma2profr[l] << "  " << sigmagamma2prof[l] << "   "
    << (dr0*l + dr0/2.)*map->inarcsec << "   "  << Aanulus*map->inarcsec*map->inarcsec << "   " <<
    std:: endl;
    lprofFORre[l] = log10(kprofr[l] + gamma2profr[l]);
    lrOFr[l] = log10(dr0*l + dr0/2.);
  }
  filoutprof.close();
  RE3 = Utilities::InterpolateYvec(lprofFORre,lrOFr,0.);
  if((RE3-lrOFr[0])<1e-4 || (RE3-lrOFr[nbin-1])<1e-4) RE3=0;
  else RE3 = pow(10.,RE3)*map->inarcsec;
}

/** \ingroup DeflectionL2
 *
 * \brief Routine for obtaining the deflection and other lensing quantities for
 * a MOKA map (MOKALensHalo), for just one ray!!
 *
 */
void LensHaloMassMap::force_halo(double *alpha
                                 ,KappaType *kappa
                                 ,KappaType *gamma
                                 ,KappaType *phi
                                 ,double const *xx
                                 ,bool subtract_point
                                 ,PosType screening
                                 )
{
  
  /*
   long index = Utilities::IndexFromPosition(xx,map->nx,map->boxlMpc/map->h,map->center);
   
   if(index > -1){
   alpha[0] = map->alpha1[index];
   alpha[1] = map->alpha2[index];
   gamma[0] = map->gamma1[index];
   gamma[1] = map->gamma2[index];
   gamma[2] = 0.0;
   *kappa = map->convergence[index];
   }
   else{
   alpha[0] = alpha[1] = 0.0;
   gamma[0] = gamma[1] = gamma[2] = 0.0;
   *kappa = 0.0;
   }
   */
  
  // interpolate from the maps
  
  Utilities::Interpolator<valarray<double> > interp(xx,map->nx,map->boxlMpc,map->ny
                                                    ,map->ny*map->boxlMpc/map->nx,map->center);

  assert(map->nx == map->ny);
  
  size_t N = map->nx * map->ny;
  
  map->convergence.size();
  
  assert(map->alpha1.size() == N);
  assert(map->alpha2.size() == N);
  assert(map->gamma1.size() == N);
  assert(map->gamma2.size() == N);
  assert(map->convergence.size() == N);
  assert(map->phi.size() == N);

  alpha[0] = interp.interpolate(map->alpha1);
  alpha[1] = interp.interpolate(map->alpha2);
  gamma[0] = interp.interpolate(map->gamma1);
  gamma[1] = interp.interpolate(map->gamma2);
  gamma[2] = 0.0;
  *kappa = interp.interpolate(map->convergence);
  *phi = interp.interpolate(map->phi);
  
  return;
}

/**
 * compute the signal of \lambda_r and \lambda_t
 */
void LensHaloMassMap::estSignLambdas(){
  map->Signlambdar.resize(map->nx*map->ny);
  map->Signlambdat.resize(map->nx*map->ny);
  double gamma,lambdar,lambdat;
  int i, j;
  for(i=0;i<map->nx;i++)
    for(j=0;j<map->ny;j++){
      gamma = sqrt(pow(map->gamma1[i+map->ny*j],2) +
                   pow(map->gamma2[i+map->ny*j],2));
      lambdat=1-map->convergence[i+map->ny*j]-gamma;
      lambdar=1-map->convergence[i+map->ny*j]+gamma;
      
      if(lambdar>=0) map->Signlambdar[i+map->ny*j]=1;
      else map->Signlambdar[i+map->ny*j]=-1;
      
      if(lambdat>=0) map->Signlambdat[i+map->ny*j]=1;
      else map->Signlambdat[i+map->ny*j]=-1;
    }
}

/**
 * measure the effective and the median Einstein radii of the connected critical
 * points present at the halo center
 */
void LensHaloMassMap::EinsteinRadii(double &RE1, double &RE2, double &xxc, double &yyc){
  double signV;
  //  std:: vector<double> xci1,yci1;
  std:: vector<double> xci2,yci2;
  // open file readable by ds9
  std::ostringstream fcrit;
  if(flag_background_field==1) fcrit << MOKA_input_file << "_only_noise_Criticals.reg";
  else{
    if(flag_MOKA_analyze == 0) fcrit << MOKA_input_file << "_noisy_Criticals.reg";
    else fcrit << MOKA_input_file << "_no_noise_Criticals.reg";
  }
  std:: ofstream filoutcrit;
  std:: string filenamecrit = fcrit.str();
  filoutcrit.open(filenamecrit.c_str());
  // define the critical points in the map
  int i, j;
  for(i=1;i<map->nx-1;i++)
    for(j=1;j<map->ny-1;j++){
      signV=map->Signlambdar[i-1+map->ny*j]+map->Signlambdar[i+map->ny*(j-1)]+
      map->Signlambdar[i+1+map->ny*j]+map->Signlambdar[i+map->ny*(j+1)];
      if(fabs(signV)<4.){
        // xci1.push_back(map->x[i]);
        // yci1.push_back(map->x[j]);
        filoutcrit << "circle(" << i << "," << j << ",0.5)" << std:: endl;
      }
      signV=map->Signlambdat[i-1+map->ny*j]+map->Signlambdat[i+map->ny*(j-1)]+
      map->Signlambdat[i+1+map->ny*j]+map->Signlambdat[i+map->ny*(j+1)];
      if(fabs(signV)<4.){
        xci2.push_back(map->x[i]);
        yci2.push_back(map->x[j]);
        filoutcrit << "circle(" << i << "," << j << ",0.5)" << std:: endl;
      }
    }
  filoutcrit.close();
  double pixDinL = map->boxlMpc/double(map->nx);
  /* measure the Einstein radius */
  std:: vector<double> xci,yci;
  //for(int ii=0;ii<xci1.size();ii++){
  //  xci.push_back(xci1[ii]);
  //  yci.push_back(yci1[ii]);
  //}
  for(int ii=0;ii<xci2.size();ii++){
    xci.push_back(xci2[ii]);
    yci.push_back(yci2[ii]);
  }
  // xci1.clear();
  // yci1.clear();
  xci2.clear();
  yci2.clear();
  int nc = xci.size();
  std:: vector<int> groupid(nc);
  int nearest0groupid = fof(pixDinL,xci,yci,groupid);
  if(nearest0groupid>0){
    std:: vector<double> xcpoints,ycpoints;
    double xercm=0;
    double yercm=0;
    for(int ii=0;ii<nc;ii++){
      if(groupid[ii] == nearest0groupid){
        xcpoints.push_back(xci[ii]);
        ycpoints.push_back(yci[ii]);
        xercm+=xci[ii];
        yercm+=yci[ii];
      }
    }
    nc = xcpoints.size();
    xercm=xercm/double(nc);
    yercm=yercm/double(nc);
    double distcentre=sqrt((xercm-xxc)*(xercm-xxc)+(yercm-yyc)*(yercm-yyc))*map->inarcsec;
    std:: vector<double>::iterator maxit, minit;
    // find the min and max elements in the vector
    maxit = max_element(xcpoints.begin(), xcpoints.end());
    minit = min_element(xcpoints.begin(), xcpoints.end());
    double xmincpoints,xmaxcpoints;
    xmaxcpoints = *maxit;
    xmincpoints = *minit;
    int imin = Utilities::locate(map->x,xmincpoints);
    imin = ((imin > 0) ? imin:0);
    imin = ((imin < map->nx-1) ? imin:map->nx-1);
    // imin=GSL_MIN( GSL_MAX( imin, 0 ), map->nx-1 );
    int imax = Utilities::locate(map->x,xmaxcpoints);
    imax = ((imax > 0) ? imax:0);
    imax = ((imax < map->nx-1) ? imax:map->nx-1);
    // imax=GSL_MIN( GSL_MAX( imax, 0 ), map->nx-1 );
    std:: vector<double> ysup,yinf,xsup,xinf;
    for(int ii=imin;ii<=imax;ii++){
      std:: vector<double>ybut;
      int condition=0;
      for(int ji=0;ji<nc;ji++){
        if(fabs(xcpoints[ji]-map->x[ii])<pixDinL/2){
          if(condition==0){
            xsup.push_back(xcpoints[ji]);
            xinf.push_back(xcpoints[ji]);
            condition=1;
          }
          ybut.push_back(ycpoints[ji]);
        }
      }
      if(ybut.size()>0){
        std:: vector<double>::iterator ymax, ymin;
        // Find the min and max elements in the vector
        ymax = max_element(ybut.begin(), ybut.end());
        ymin = min_element(ybut.begin(), ybut.end());
        double ymincpoints,ymaxcpoints;
        ymaxcpoints = *ymax;
        ymincpoints = *ymin;
        ysup.push_back(ymaxcpoints);
        yinf.push_back(ymincpoints);
      }
      if(ybut.size()==1){
        double ymincpoints,ymaxcpoints;
        ymaxcpoints = ybut[0];
        ymincpoints = ybut[0];
        ysup.push_back(ymaxcpoints);
        yinf.push_back(ymincpoints);
      }
    }
    nc = yinf.size();
    int npixIN=0;
    std:: vector<double> RE;
    for(int ii=0;ii<nc;ii++){
      RE.push_back(sqrt(pow(xinf[ii]-xercm,2.) + pow(yinf[ii]-yercm,2)));
      std:: vector<double> ycounts;
      for(int ji=0;ji<map->nx;ji++){
        if(map->x[ji]>=yinf[ii] && map->x[ji]<=ysup[ii]){
          ycounts.push_back(map->x[ji]);
        }
      }
      int ncounts = ycounts.size();
      npixIN=npixIN+ncounts;
    }
    for(int ii=nc-1;ii>=0;ii--){
      RE.push_back(sqrt(pow(xsup[ii]-xercm,2) + pow(ysup[ii]-yercm,2)));
    }
    RE1=map->inarcsec*sqrt(pixDinL*pixDinL*npixIN/M_PI);
    RE2=map->inarcsec*Utilities::median(RE);
    // if is not in the centre
    std:: cout << "distance " <<  distcentre << std:: endl;
    if(distcentre>1.5*RE2){
      RE1=-RE1;
      RE2=-RE2;
    }
    if(RE2!=RE2) RE2=0.;
  }
  else{
    RE1=0.;
    RE2=0.;
  }
}
/**
 * saves MAP properties, computing the radial profile
 * of the convergence and shear
 */
void LensHaloMassMap::saveImage(bool saveprofiles){
  std::stringstream f;
  std::string filename;
  if(flag_background_field==1) f << MOKA_input_file << "_only_noise.fits";
  else{
    if(flag_MOKA_analyze == 0) f << MOKA_input_file << "_noisy.fits";
    else f << MOKA_input_file << "_no_noise.fits";
  }
  filename = f.str();
  
  writeImage(filename);
  
  if(saveprofiles == true){
    std:: cout << " saving profile " << std:: endl;
    double RE3,xxc,yyc;
    saveProfiles(RE3,xxc,yyc);
    estSignLambdas();
    double RE1,RE2;
    EinsteinRadii(RE1,RE2,xxc,yyc);
    std::ostringstream fEinr;
    if(flag_background_field==1) fEinr << MOKA_input_file << "_only_noise_Einstein.radii.dat";
    else{
      if(flag_MOKA_analyze == 0) fEinr << MOKA_input_file << "_noisy_Einstein.radii.dat";
      else fEinr << MOKA_input_file << "_no_noise_Einstein.radii.dat";
    }
    std:: ofstream filoutEinr;
    std:: string filenameEinr = fEinr.str();
    filoutEinr.open(filenameEinr.c_str());
    filoutEinr << "# effective        median      from_profles" << std:: endl;
    filoutEinr << RE1 << "   " << RE2 << "    " << RE3 << std:: endl;
    filoutEinr.close();
  }
}

