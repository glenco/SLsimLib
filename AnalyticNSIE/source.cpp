/*
 * source.cpp
 *
 *  Created on: Feb 8, 2012
 *      Author: mpetkova
 */

#include <typeinfo>
#include "cpfits.h"
#include "source.h"
#include "source_models.h"
#include "sourceAnaGalaxy.h"
#include "image_processing.h"

using namespace std;

double flux_to_mag_AB(double flux){
  if(flux <=0) return 100;
  return -2.5 * log10(flux) - 48.6;
}

double mag_to_flux_AB(double m){
  if(m == 100) return 0;
  return pow(10,-0.4*(m+48.6));
}

double jansky_to_mag_AB(double flux){
  if(flux <=0) return 100;
  return -2.5 * log10(flux) + 8.9;
}

double mag_to_jansky_AB(double m){
  if(m == 100) return 0;
  return pow(10,-0.4*(m - 8.9));
}

double flux_to_mag(double flux,double zeropoint){
  if(flux <=0) return 100;
  return -2.5 * log10(flux) + zeropoint;
}

double mag_to_flux(double m,double zeropoint){
  if(m == 100) return 0;
  return pow(10,-0.4*(m - zeropoint));
}

//SourceUniform::SourceUniform(InputParams& params) : Source(){
//  assignParams(params);
//}

SourceUniform::SourceUniform(Point_2d position,PosType z,PosType radius_in_radians):
Source(radius_in_radians,position,z,-1,-48.6)
{
}


//SourceGaussian::SourceGaussian(InputParams& params) : Source(){
//  assignParams(params);
//}

//SourceBLR::SourceBLR(InputParams& params) : Source(){
//  assignParams(params);
//}
//
//SourceBLRDisk::SourceBLRDisk(InputParams& params) : SourceBLR(params){
//
//}
//
//SourceBLRSph1::SourceBLRSph1(InputParams& params) : SourceBLR(params){
//
//}
//
//SourceBLRSph2::SourceBLRSph2(InputParams& params) : SourceBLR(params){
//
//}

Source::~Source(){
}

SourceUniform::~SourceUniform(){
}

SourceGaussian::~SourceGaussian(){
}

SourceBLR::~SourceBLR(){
}

SourceBLRDisk::~SourceBLRDisk(){
}

SourceBLRSph1::~SourceBLRSph1(){
}

SourceBLRSph2::~SourceBLRSph2(){
}

void SourceUniform::assignParams(InputParams& params){
  
  if(!params.get("z_source",zsource)){
		  ERROR_MESSAGE();
		  cout << "parameter z_source needs to be set in parameter file " << params.filename() << endl;
		  exit(0);
  }
  
  printSource();
}

void SourceGaussian::assignParams(InputParams& params){
  
  
  if(!params.get("z_source",zsource)){
		  ERROR_MESSAGE();
		  cout << "parameter z_source needs to be set in parameter file " << params.filename() << endl;
		  exit(0);
  }
  
  if(!params.get("gauss_r2",source_gauss_r2)){
		  ERROR_MESSAGE();
		  cout << "parameter gauss_r2 needs to be set in parameter file " << params.filename() << endl;
		  exit(0);
  }
  
  source_r = 5*sqrt(source_gauss_r2);
  
  printSource();
}

void SourceBLR::assignParams(InputParams& params){
  
  bool fail = false;
  if(!params.get("source_z",zsource)){
		  ERROR_MESSAGE();
		  cout << "parameter source_z_source needs to be set in parameter file " << params.filename() << endl;
		  fail = true;
  }
  if(!params.get("source_BHmass",source_BHmass)){
		  ERROR_MESSAGE();
		  cout << "parameter source_BHmass needs to be set in parameter file " << params.filename() << endl;
		  fail = true;
  }
  
  if(!params.get("source_gamma",source_gamma)){
		  ERROR_MESSAGE();
		  cout << "parameter source_gamma needs to be set in parameter file " << params.filename() << endl;
		  fail = true;
  }
  if(!params.get("source_inclin",source_inclination)){
		  ERROR_MESSAGE();
		  cout << "parameter source_inclin needs to be set in parameter file " << params.filename() << endl;
		  fail = true;
  }
  if(!params.get("source_opening_ang",source_opening_angle)){
		  ERROR_MESSAGE();
		  cout << "parameter source_opening_ang needs to be set in parameter file " << params.filename() << endl;
		  fail = true;
  }
  if(!params.get("source_r_in",source_r_in)){
		  ERROR_MESSAGE();
		  cout << "parameter source_r_in needs to be set in parameter file " << params.filename() << endl;
		  fail = true;
  }
  if(!params.get("source_r_out",source_r_out)){
		  ERROR_MESSAGE();
		  cout << "parameter source_r_out needs to be set in parameter file " << params.filename() << endl;
		  fail = true;
  }
  if(!params.get("source_nuo",source_nuo)){
		  ERROR_MESSAGE();
		  cout << "parameter source_nuo needs to be set in parameter file " << params.filename() << endl;
		  fail = true;
  }
  if(!params.get("source_fK",source_fK)){
		  ERROR_MESSAGE();
		  cout << "parameter source_fK needs to be set in parameter file " << params.filename() << endl;
		  fail = true;
  }
  
  if(fail) throw std::runtime_error("In SourceBLR");
  
  source_inclination *= PI/180;
  source_opening_angle *= PI/180;
  source_monocrome = false;
  
  printSource();
}

void SourceUniform::printSource(){
  cout << endl << "**Source model**" << endl;
  
  cout << "z_source " << zsource << endl << endl;
}

void SourceGaussian::printSource(){
  cout << endl << "**Source model**" << endl;
  
  cout << "z_source " << zsource << endl;
  cout << "gauss_r2 " << source_gauss_r2 << endl << endl;
}

void SourceBLR::printSource(){
  cout << endl << "**Source model**" << endl;
  
  cout << "z_source " << zsource << endl;
  cout << "BHmass " << source_BHmass << endl;
  cout << "gamma " << source_gamma << endl;
  cout << "incl " << source_inclination*180/PI << endl;
  cout << "opening angl " << source_opening_angle*180/PI << endl;
  cout << "r_in " << source_r_in << endl;
  cout << "r_out " << source_r_out << endl;
  cout << "nuo " << source_nuo << endl;
  cout << "source_fK " << source_fK << endl << endl;
}

PosType SourceUniform::SurfaceBrightness(const PosType *y) const{
  return (PosType)( (pow(y[0]-getTheta()[0],2) + pow(y[1]-getTheta()[1],2)) < source_r*source_r );
}

PosType SourceGaussian::SurfaceBrightness(const PosType *y) const{
  return exp( -(pow(y[0]-getTheta()[0],2) + pow(y[1]-getTheta()[1],2))/source_gauss_r2 );
}
// surface brightness for models of the Broad Line Region
PosType SourceBLRDisk::SurfaceBrightness(const PosType *y) const{
  PosType x[2] = {y[0]-getTheta()[0],y[1]-getTheta()[1]};
  return blr_surface_brightness_disk(x,this);
}

PosType SourceBLRSph1::SurfaceBrightness(const PosType *y) const{
  return blr_surface_brightness_spherical_circular_motions(sqrt((pow(y[0]-getTheta()[0],2) + pow(y[1]-getTheta()[1],2))),this);
}
PosType SourceBLRSph2::SurfaceBrightness(const PosType *y) const{
  return blr_surface_brightness_spherical_random_motions(sqrt((pow(y[0]-getTheta()[0],2) + pow(y[1]-getTheta()[1],2))),this);
}

//void in_source(PosType *y_source,ListHndl sourcelist){
//  return;
//}
//SourcePixelled::SourcePixelled(InputParams& params)
//{}

SourcePixelled::SourcePixelled(
                               PosType my_z            /// redshift of the source
                               , PosType* my_center  /// center (in rad)
                               , int my_Npixels           /// number of pixels per side
                               , PosType my_resolution  /// resolution (in rad)
                               , PosType* arr_val          /// array of pixel values (must be of size = Npixels*Npixels)
                               , PosType zero_point        /// magnitude zero point
)
:Source(0,my_center,my_z,-1,zero_point), resolution(my_resolution), Npixels (my_Npixels){
  zsource = my_z;
  
  range = resolution*(Npixels-1);
  
  values.resize(Npixels*Npixels);
  for (int i = 0; i < Npixels*Npixels; i++)
    values[i] = arr_val[i];
  source_r =  range/sqrt(2.);
  calcTotalFlux();
  calcCentroid();
  calcEll();
  calcSize();
}

/** \brief Creates a SourcePixelled from a PixelMap image.
 *  The idea is to use stamps of observed galaxies as input sources for simuations.
 *  Surface brightness of the source is conserved, taking into account the input pixel size.
 *  Factor allows for rescaling of the flux, in case one wants to simulate a different observation.
 */
SourcePixelled::SourcePixelled(
                               const PixelMap& gal_map  /// Input image and information
                               , PosType my_z                 /// redshift of the source
                               , PosType factor                /// optional rescaling factor for the flux
                               , PosType zero_point
)
:Source(0,Point_2d(0,0),0,-1,zero_point){
  if(gal_map.getNx() != gal_map.getNy()){
    std::cout << "SourcePixelled::SourcePixelled() Doesn't work on nonsquare maps" << std::endl;
    throw std::runtime_error("nonsquare");
  }
  
  zsource = my_z;
  resolution = gal_map.getResolution();
  Npixels = gal_map.getNx();
  range = resolution*(Npixels-1);
  source_x[0] = gal_map.getCenter()[0];
  source_x[1] = gal_map.getCenter()[1];
  source_r =  range/sqrt(2.);
  values.resize(Npixels*Npixels);
  
  double convertion = 1.0/resolution/resolution*factor;
  for (int i = 0; i < Npixels*Npixels; i++)
    values[i] = gal_map(i)*convertion;
  
  calcTotalFlux();
  calcCentroid();
  calcEll();
  calcSize();
}

SourcePixelled::~SourcePixelled(){
}

void SourcePixelled::calcCentroid(){
  PosType x_sum = 0;
  PosType y_sum = 0;
  PosType sum = 0;
  PosType x[2];
  for (unsigned long i = 0; i < Npixels*Npixels; i++)
  {
    Utilities::PositionFromIndex(i,x,Npixels,range,source_x.x);
    x_sum += x[0]*values[i];
    y_sum += x[1]*values[i];
    sum += values[i];
  }
  centroid[0] = x_sum/sum;
  centroid[1] = y_sum/sum;
}

void SourcePixelled::calcEll(){
  for (int j = 0; j < 2; j++)
  {
    for (int k = 0; k < 2; k++)
    {
      quad[j][k] = 0.;
    }
  }
  
  PosType x[2];
  for (unsigned long i = 0; i < Npixels*Npixels; i++)
  {
    Utilities::PositionFromIndex(i,x,Npixels,range,source_x.x);
    for (int j = 0; j < 2; j++)
    {
      for (int k = 0; k < 2; k++)
      {
        quad[j][k] += values[i]*(x[j]-centroid[j])*(x[k]-centroid[k]);
      }
    }
  }
  
  ell[0] = (quad[0][0] - quad[1][1])/(quad[0][0] + quad[1][1]);
  ell[1] = 2*quad[0][1]/(quad[0][0] + quad[1][1]);
}

void SourcePixelled::calcSize(){
  PosType r_sum = 0.;
  PosType sum = 0.;
  PosType rad;
  PosType x[2];
  
  for (unsigned long i = 0; i < Npixels*Npixels; i++)
  {
    Utilities::PositionFromIndex(i,x,Npixels,range,source_x.x);
    rad = sqrt((x[0]-centroid[0])*(x[0]-centroid[0])+(x[1]-centroid[1])*(x[1]-centroid[1]));
    r_sum += rad*values[i];
    sum += values[i];
  }
  size = r_sum/sum;
}

PosType SourcePixelled::SurfaceBrightness(const PosType *y) const{
  long ix = Utilities::IndexFromPosition(y[0],Npixels,range,source_x[0]);
  long iy = Utilities::IndexFromPosition(y[1],Npixels,range,source_x[1]);
  if (ix>-1 && iy>-1)
  {
    return values[iy*Npixels+ix];
  }
  else
    return 0.;
}

void SourcePixelled::calcTotalFlux(){
  PosType val_tot = 0.;
  for (int i = 0; i < Npixels*Npixels; i++)
    val_tot += values[i];
  flux = val_tot*resolution*resolution;
}

void SourcePixelled::printSource(){}
void SourcePixelled::assignParams(InputParams& params){}

/**  \brief Calculates the difference in magnitude when changing the observing filter
 *
 */
PosType Source::changeFilter(
                             std::string filter_in  /// file with the old observing filter
                             , std::string filter_out	/// file with the new observing filter
                             , std::string sed			/// file with the galaxy spectral energy distribution
)
{
  
  // reads in the input filter
  std::ifstream fin(filter_in.c_str());
  std::vector<PosType> wavel_in, ampl_in;
  fin.seekg(0,fin.beg);
  PosType x,y;
  for (int i = 0; ; i++)
  {
    if( fin.eof() ) break;
    if (fin.peek() == '#')
    {
      fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      continue;
    }
    fin >> x;
    wavel_in.push_back(x);
    fin >> y;
    ampl_in.push_back(y);
  }
  
  // reads in the output filter
  std::ifstream fout(filter_out.c_str());
  std::vector<PosType> wavel_out, ampl_out;
  fout.seekg(0,fout.beg);
  for (int i = 0; ; i++)
  {
    if (fout.eof() ) break;
    if (fout.peek() == '#')
    {
      fout.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      continue;
    }
    fout >> x;
    wavel_out.push_back(x);
    fout >> y;
    ampl_out.push_back(y);
  }
  
  // reads in the source sed
  std::ifstream sed_input(sed.c_str());
  std::vector<PosType> wavel_sed, ampl_sed;
  sed_input.seekg(0,sed_input.beg);
  for (int i = 0; ; i++)
  {
    if(sed_input.eof() ) break;
    if (sed_input.peek() == '#')
    {
      sed_input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      continue;
    }
    sed_input >> x;
    wavel_sed.push_back(x);
    sed_input >> y;
    ampl_sed.push_back(y);
    wavel_sed[i] *= 1+zsource;
  }
  
  // Applies Delta m = int (sed*fout) / int (sed*fin) * int fin / int fout
 	PosType fin_int = integrateFilter(wavel_in,ampl_in);
  PosType fout_int = integrateFilter(wavel_out,ampl_out);
  PosType sed_in = integrateFilterSED(wavel_in, ampl_in, wavel_sed, ampl_sed);
  PosType sed_out = integrateFilterSED(wavel_out, ampl_out, wavel_sed, ampl_sed);
  
  if (sed_in < std::numeric_limits<PosType>::epsilon() || sed_out < std::numeric_limits<PosType>::epsilon())
  {
    return 99.;
  }
  else
  {
    PosType delta_mag = -2.5 * log10(sed_out/sed_in*fin_int/fout_int);
    return std::min(delta_mag,99.);
  }
}

/**  \brief Calculates the integral of the filter curve given as an array of (x,y) values.
 *
 */
PosType Source::integrateFilter(std::vector<PosType> wavel_fil, std::vector<PosType> fil)
{
  PosType integr = 0.;
  for (int i = 0; i < wavel_fil.size()-1; i++)
    integr += (wavel_fil[i+1] - wavel_fil[i])*(fil[i+1]+fil[i]);
  integr /= 2.;
  return integr;
}

/**  \brief Calculates the integral of the sed multiplied by the filter curve.
 *
 */
PosType Source::integrateFilterSED(std::vector<PosType> wavel_fil, std::vector<PosType> fil, std::vector<PosType> wavel_sed, std::vector<PosType> sed)
{
  int wavel_new_size = 10000;
  vector<PosType> wavel_new(wavel_new_size), fil_val(wavel_new_size), sed_val(wavel_new_size);
  
  PosType integr = 0.;
  for (int i = 0; i < wavel_new_size; i++)
  {
    wavel_new[i] = MAX<PosType>(wavel_fil[0],wavel_sed[0]) + PosType(i)/PosType(wavel_new_size-1)*(MIN<PosType>(wavel_fil[wavel_fil.size()-1],wavel_sed[wavel_sed.size()-1])-MAX<PosType>(wavel_fil[0],wavel_sed[0]) );
    vector<PosType>::iterator it_f = lower_bound(wavel_fil.begin(), wavel_fil.end(), wavel_new[i]);
    long p = it_f-wavel_fil.begin()-1;
    vector<PosType>::iterator it_s = lower_bound(wavel_sed.begin(), wavel_sed.end(), wavel_new[i]);
    long q = it_s-wavel_sed.begin()-1;
    if (p >= 0&& p < wavel_fil.size()-1 && q >= 0 && q < wavel_sed.size())
    {
      fil_val[i] = fil[p] + (fil[p+1]-fil[p])*(wavel_new[i]-wavel_fil[p])/(wavel_fil[p+1]-wavel_fil[p]);
      sed_val[i] = sed[q] + (sed[q+1]-sed[q])*(wavel_new[i]-wavel_sed[q])/(wavel_sed[q+1]-wavel_sed[q]);
    }
    else
    {
      fil_val[i] = 0.;
      sed_val[i] = 0.;
    }
  }
  for (int i = 0; i < wavel_new_size-1; i++)
  {
    integr += (wavel_new[i+1] - wavel_new[i])*(fil_val[i+1]*sed_val[i+1]*wavel_new[i+1]*wavel_new[i+1]+fil_val[i]*sed_val[i]*wavel_new[i]*wavel_new[i]);
  }
  integr /= 2.;
  return integr;
}

/*
 SourceShapelets::SourceShapelets()
 : 	mag(0),ang(0),n1(0),n2(0)
 {
 setX(0, 0);
 zsource=0;
 }
 */

size_t SourceShapelets::count = 0;;

void SourceColored::setActiveBand(Band band)
{
  
  if(mag_map.size() == 0) return;
  
  mag = mag_map.at(band);
  if (mag < 0.)
      flux_total = std::numeric_limits<PosType>::epsilon();
    else
      flux_total = mag_to_flux(mag);
  
  assert(flux_total > 0);
  current_band = band;
  return;
}

SourceShapelets::SourceShapelets(
                                 PosType my_z                  /// redshift of the source
                                 , PosType my_mag							/// magnitude
                                 , PosType my_scale						/// scale of the shapelets decomposition
                                 , std::valarray<PosType> my_coeff  	/// coefficients of the shapelets decomposition
                                 , PosType* my_center           			/// center (in rad)
                                 , PosType my_ang					/// rotation angle (in rad)
                                 , PosType zeropoint       /// magnitude zero point
)
:SourceColored(my_mag,my_scale,Point_2d(my_center[0],my_center[1]),zsource,-1,zeropoint)
{

  assert(my_center != NULL);

  cos_sin[0] = cos(my_ang);
  cos_sin[1] = sin(my_ang);

  n1 = sqrt(my_coeff.size());
  n2 = n1;
  coeff = my_coeff;
  
  assert(flux_total > 0.0);
  
  NormalizeFlux();
  
  ++count;
}

SourceShapelets::SourceShapelets(
                                 PosType my_z							/// redshift of the source
                                 , PosType my_mag						/// magnitude
                                 , std::string shap_file				/// fits file with coefficients in a square array
                                 , PosType* my_center           			/// center (in rad)
                                 , PosType my_ang			/// rotation angle (in rad)
                                 , PosType zeropoint       /// magnitude zero point
)
:SourceColored(my_mag,0,Point_2d(my_center[0],my_center[1]),my_z,-1,zeropoint)
{
  
  assert(my_center != NULL);
//  if(my_center != NULL)
//    setTheta(my_center[0], my_center[1]);
//  else
//    setTheta(0, 0);
  
  cos_sin[0] = cos(my_ang);
  cos_sin[1] = sin(my_ang);
  
  if(shap_file.empty())
    throw std::invalid_argument("Please enter a valid filename for the FITS file input");
  
  CPFITS_READ cpfits(shap_file.c_str());
  
  cpfits.readKey("BETA", source_r);
  source_r *= 0.03*arcsecTOradians;
  cpfits.readKey("DIM", n1);
  cpfits.readKey("ID", id);
  n2 = n1;
  vector<long> size;
  cpfits.read(coeff,size);

  assert(flux_total > 0);
  
  NormalizeFlux();
  
  current_band = Band::NoBand;
  ++count;
}

SourceShapelets::SourceShapelets(
                                 std::string shap_file				/// fits file with coefficients in a square array. Mag and redshift are read from the header.
                                 , PosType my_ang         /// rotation angle (in rad)
                                 , PosType zeropoint       /// magnitude zero point
 )
:SourceColored(0,0,Point_2d(0,0),0,-1,zeropoint)
{
  
  cos_sin[0] = cos(my_ang);
  cos_sin[1] = sin(my_ang);

  if(shap_file.empty())
    throw std::invalid_argument("Please enter a valid filename for the FITS file input");
 
  CPFITS_READ cpfits(shap_file.c_str());
  
  cpfits.readKey("BETA", source_r);
  source_r *= 0.03*arcsecTOradians;

  cpfits.readKey("SED_TYPE",sed_type);
  
  cpfits.readKey("MAG_B",mag_map[Band::F435W]); // ACS F435W band magnitude
  cpfits.readKey("MAG_V",mag_map[Band::F606W]); // ACS F606W band magnitude
  cpfits.readKey("MAG_I",mag_map[Band::F775W]); // ACS F775W band magnitude
  cpfits.readKey("MAG_Z",mag_map[Band::F850LP]);// ACS F850LP band magnitude
  cpfits.readKey("MAG_J",mag_map[Band::F110W]); // NIC3 F110W band magnitude
  cpfits.readKey("MAG_H",mag_map[Band::F160W]);  // NIC3 F160W band magnitude
  cpfits.readKey("MAG_u_KIDS",mag_map[Band::KiDS_U]); // u band obtained from SED fitting
  cpfits.readKey("MAG_g_KIDS",mag_map[Band::KiDS_G]); // g band obtained from SED fitting
  cpfits.readKey("MAG_r_KIDS",mag_map[Band::KiDS_R]); // r band obtained from SED fitting
  cpfits.readKey("MAG_i_KIDS",mag_map[Band::KiDS_I]); // i band obtained from SED fitting

  
  // by default, the magnitude is the one in the i band,
  // whose image has been used for shapelets decomposition
  setActiveBand(Band::KiDS_I);
  
  cpfits.readKey("REDSHIFT", zsource);
  //cpfits.readKey("ID", id); // I'm not sure what this is.
  cpfits.readKey("DIM", n1);
  
  n2 = n1;
  std::vector<long> size;
  cpfits.read(coeff,size);

  
  // ??? kluge
  mag_map[Band::EUC_VIS] = mag_map.at(Band::KiDS_I);
  //mag_map[Band::EUC_Y] = mag_map.at(Band::F110W);
  mag_map[Band::EUC_J] = mag_map.at(Band::F110W);
  mag_map[Band::EUC_H] = mag_map.at(Band::F160W);
  
  NormalizeFlux();
  ++count;
}

/// Returns surface brightness in erg/cm2/sec/Hz, normalized by hplanck.
/// Given the units of hplanck, the final units are 1/sec/cm2.
PosType SourceShapelets::SurfaceBrightness(const PosType *y)
const{
  PosType sb = 0.;
  PosType y_norm[2],tmp;
  y_norm[0] = ((y[0]-source_x[0])*cos_sin[0]-(y[1]-source_x[1])*cos_sin[1])/source_r;
  y_norm[1] = ((y[0]-source_x[0])*cos_sin[1]+(y[1]-source_x[1])*cos_sin[0])/source_r;
  //PosType dist = sqrt(y_norm[0]*y_norm[0]+y_norm[1]*y_norm[1]);
  
  
  double r_norm2 = (y_norm[0]*y_norm[0]+y_norm[1]*y_norm[1])/2;
  
  if(r_norm2 > 10) return 0;
  
  std::vector<PosType> Her1,Her2;
  
  Hermite(Her1,n1,y_norm[0]);
  Hermite(Her2,n2,y_norm[1]);
  
  size_t coei=1,coej=1;
  
  for (int i = 0; i < n1; i++,coei *= 2,coej=1 )
  {
    tmp = factrl(i)*PI;
    for (int j = 0; j < n2; j++,coej *= 2 )
    {
      
      PosType norm = 1./sqrt(coei*coej*tmp*factrl(j));
      sb += norm*coeff[j*n1+i]*Her1[i]*Her2[j];
    }
  }
  sb *= exp(- r_norm2 )/source_r;
  sb *= flux_total/coeff_flux;

  assert(flux_total > 0);
  
  
  return MAX<double>(sb,0);
  //return max(sb,std::numeric_limits<PosType>::epsilon());
}

/// Returns the value of the Hermite polynomials from degree 0 to n at position x
void SourceShapelets::Hermite(std::vector<PosType> &hg,int N, PosType x) const
{
  hg.resize(N);
  hg[0] = 1.;
  for (int i = 1; i < N; i++)
  {
    if (i==1)
      hg[1] = 2.*x;
    else
      hg[i] = 2.*x*hg[i-1]-2.*(i-1)*hg[i-2];
  }
  return;
}

void SourceShapelets::printSource()
{
  cout << endl << "**Source model**" << endl;
  
  cout << "z_source " << zsource << endl;
  cout << "mag " << mag << endl;
  
};

void SourceShapelets::assignParams(InputParams& params){};

/// Rescales the coefficients to make the source bright as we want.
void SourceShapelets::NormalizeFlux()
{
  coeff_flux = 0.;
  for (int i = 0; i < n1; i=i+2)
  {
    for (int j = 0; j < n2; j=j+2)
    {
      coeff_flux += pow(2,0.5*(2-i-j))*sqrt(factrl(i))/factrl(i/2.)*sqrt(factrl(j))/factrl(j/2.)*coeff[j*n1+i];
    }
  }
  coeff_flux *= sqrt(PI)*source_r;
}

// Default constructor. Reads in sources from the default catalog. No magnitude limit.
//SourceMultiShapelets::SourceMultiShapelets(InputParams& params)
//: Source(),index(0)
//{
//  assignParams(params);
//  readCatalog();
//}

SourceMultiShapelets::SourceMultiShapelets(
                                           const std::string &my_shapelets_folder,Band my_band
                                           ,double my_max_mag_limit
                                           ,double my_min_mag_limit
                                           ,double my_z_max
                                           ,double my_sb_limit
                                           ,double maximum_radius
                                           ,double zeropoint
                                           )
: Source(0,Point_2d(0,0),0,my_sb_limit,zeropoint),index(0),max_mag_limit(my_max_mag_limit),min_mag_limit(my_min_mag_limit)
,band(my_band),radius_max(maximum_radius),shapelets_folder(my_shapelets_folder),z_max(my_z_max)
{
  
//  if(sb_limit == -1)
//    setSBlimit_magarcsec(30.);
//  else
//    sb_limit = flux_to_mag(sb_limit)*pow(180*60*60/PI,2);
  
  readCatalog();
}

void SourceMultiShapelets::input(const std::string &my_shapelets_folder,Band my_band
                                 ,double my_max_mag_limit,double my_min_mag_limit,double my_z_max
                                 ,double my_sb_limit,double maximum_radius,double zero_point
)
{
  
  index=0;
  max_mag_limit = my_max_mag_limit;
  min_mag_limit = my_min_mag_limit;
  band = my_band;
  radius_max = maximum_radius;
  shapelets_folder = my_shapelets_folder;
  setMagZeroPoint(zero_point);
  z_max = my_z_max;          /// maximum redshift

  
 // if(sb_limit == -1)
 //   setSBlimit_magarcsec(30.);
 // else
 //   sb_limit = flux_to_mag(sb_limit)*pow(180*60*60/PI,2);
  
  readCatalog();
}


SourceMultiShapelets::~SourceMultiShapelets()
{
}

/// Reads in the default catalog
void SourceMultiShapelets::readCatalog()
{
  
  // read in shaplet catalogs
  
  std::vector<std::vector<float> > viz_cat;
  std::vector<std::string> col_names;
  
  Utilities::IO::ReadCSVnumerical2<float>(
                                          shapelets_folder + "/euclid_cats/euclid_riz.cat"
                                          ,viz_cat
                                          ,col_names
                                          ,1000000
                                          ,'#'
                                          ,' '
                                          ,false);
  
  
  //sort by id number
  std::sort(viz_cat.begin(),viz_cat.end()
            ,[](const std::vector<float> &a,const std::vector<float> &b){return a[1] < b[1]; } );
  int ii=0;
  for(auto &a : viz_cat){
    if(a[1] >= 0) break;
    ++ii;
  }
  
  std::vector<std::vector<float> > h_cat;
  Utilities::IO::ReadCSVnumerical2<float>(
                                          shapelets_folder + "/euclid_cats/euclid_H.cat"
                                          ,h_cat
                                          ,col_names
                                          ,1000000
                                          ,'#'
                                          ,' '
                                          ,false);
  std::sort(h_cat.begin(),h_cat.end()
            ,[](const std::vector<float> &a,const std::vector<float> &b){return a[1] < b[1]; } );

  std::vector<std::vector<float> > y_cat;
  Utilities::IO::ReadCSVnumerical2<float>(
                                          shapelets_folder + "/euclid_cats/euclid_Y.cat"
                                          ,y_cat
                                          ,col_names
                                          ,1000000
                                          ,'#'
                                          ,' '
                                          ,false);
  std::sort(y_cat.begin(),y_cat.end()
            ,[](const std::vector<float> &a,const std::vector<float> &b){return a[1] < b[1]; } );

  std::vector<std::vector<float> > j_cat;
  Utilities::IO::ReadCSVnumerical2<float>(
                                          shapelets_folder + "/euclid_cats/euclid_J.cat"
                                          ,j_cat
                                          ,col_names
                                          ,1000000
                                          ,'#'
                                          ,' '
                                          ,false);
  std::sort(j_cat.begin(),j_cat.end()
            ,[](const std::vector<float> &a,const std::vector<float> &b){return a[1] < b[1]; } );
  
  
  int j=ii;
  int max_num = 37012;
  for (int i = 0; i < max_num+1; i++)
  {
    std::string shap_file = shapelets_folder+"/obj_"+ std::to_string(i)+"_mag_z.sif";
    std::ifstream shap_input(shap_file.c_str());
    if (shap_input)
    {
      SourceShapelets s(shap_file.c_str(),0,getMagZeroPoint());
      
      s.setID(i);

      s.sed_type = viz_cat[j][4];
      
      assert(viz_cat.size() > j );
      s.setBand(Band::EUC_VIS,viz_cat[j][2]);
      assert(y_cat.size() > j );
      s.setBand(Band::EUC_Y,y_cat[j][2]);
      assert(j_cat.size() > j );
      s.setBand(Band::EUC_J,j_cat[j][2]);
      assert(h_cat.size() > j );
      s.setBand(Band::EUC_H,h_cat[j++][2]);
      
      //s.setActiveBand(band);
      if (s.getMag() > 0.
          && s.getMag(Band::EUC_VIS) < max_mag_limit
          && s.getMag(Band::EUC_VIS) > min_mag_limit
          && s.getMag(Band::EUC_J) > 0
          && s.getMag(Band::EUC_H) > 0
          && s.getRadius() < radius_max
          && s.getZ() < z_max){
        galaxies.push_back(s);
      }
      shap_input.close();
      /*else if (i == 1)
      {
        std::cout << "Can't open file " << shap_file << std::endl;
        ERROR_MESSAGE();
        throw std::runtime_error(" Cannot open file.");
        exit(1);
      }*/
    }
  }
  std::cout << galaxies.size() << " shapelet sources out of "
  << max_num  << " passed selection." << std::endl;
  band = galaxies[0].getBand();
}


void SourceMultiShapelets::assignParams(InputParams& params){
  if(!params.get("source_mag_limit",max_mag_limit)){
    std::cerr << "ERROR: Must assign source_mag_limit in parameter file " << params.filename() << std::endl;
    exit(1);
  }
  if(!params.get("source_min_mag_limit",min_mag_limit)){
    min_mag_limit = -1;
  }

  
  if(!params.get("source_sb_limit",sb_limit))
    sb_limit = -1;
  else
    sb_limit = SBlimit_magarcsec(sb_limit);
  if(!params.get("shapelets_folder",shapelets_folder)){
    std::cerr << "ERROR: shapelets_folder not found in parameter file " << params.filename() << std::endl;
    exit(1);
  }
  
  if(!params.get("shapelets_band",band)){
    std::cerr << "ERROR: Must assign shapelets_band in parameter file " << params.filename() << std::endl;
    exit(1);
  }
  
  
}

/// Print info on current source parameters
void SourceMultiShapelets::printSource(){
  std::cout << "Shapelets" << std::endl;
  galaxies[index].printSource();
}
/// Sort the sources by redshift in assending order
void SourceMultiShapelets::sortInRedshift(){
  std::sort(galaxies.begin(),galaxies.end(),redshiftcompare_shap);
}
// used in MultiSourceShapelets::sortInRedshift()
bool redshiftcompare_shap(SourceShapelets s1,SourceShapelets s2){
  return (s1.getZ() < s2.getZ());
}
/// Sort the sources by magnitude in assending order
void SourceMultiShapelets::sortInMag(){
  std::sort(galaxies.begin(),galaxies.end(),magcompare_shap);
}
// used in MultiSourceShapelets::sortInRedshift()
bool magcompare_shap(SourceShapelets s1,SourceShapelets s2){
  return (s1.getMag() < s2.getMag());
}



