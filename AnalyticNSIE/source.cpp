/*
 * source.cpp
 *
 *  Created on: Feb 8, 2012
 *      Author: mpetkova
 */

#include "slsimlib.h"
#include <typeinfo>

#ifdef ENABLE_FITS
#include <CCfits/CCfits>
//#include <CCfits>
#endif

using namespace std;

SourceUniform::SourceUniform(InputParams& params) : Source(){
  assignParams(params);
}

SourceGaussian::SourceGaussian(InputParams& params) : Source(){
  assignParams(params);
}

SourceBLR::SourceBLR(InputParams& params) : Source(){
  assignParams(params);
}

SourceBLRDisk::SourceBLRDisk(InputParams& params) : SourceBLR(params){
  
}

SourceBLRSph1::SourceBLRSph1(InputParams& params) : SourceBLR(params){
  
}

SourceBLRSph2::SourceBLRSph2(InputParams& params) : SourceBLR(params){
  
}

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
  
  if(fail) exit(1);
  
  source_inclination *= pi/180;
  source_opening_angle *= pi/180;
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
  cout << "incl " << source_inclination*180/pi << endl;
  cout << "opening angl " << source_opening_angle*180/pi << endl;
  cout << "r_in " << source_r_in << endl;
  cout << "r_out " << source_r_out << endl;
  cout << "nuo " << source_nuo << endl;
  cout << "source_fK " << source_fK << endl << endl;
}

PosType SourceUniform::SurfaceBrightness(PosType *y){
  return (PosType)( (pow(y[0]-getX()[0],2) + pow(y[1]-getX()[1],2)) < source_r*source_r );
}

PosType SourceGaussian::SurfaceBrightness(PosType *y){
  return exp( -(pow(y[0]-getX()[0],2) + pow(y[1]-getX()[1],2))/source_gauss_r2 );
}
// surface brightness for models of the Broad Line Region
PosType SourceBLRDisk::SurfaceBrightness(PosType *y){
  PosType x[2] = {y[0]-getX()[0],y[1]-getX()[1]};
  return blr_surface_brightness_disk(x,this);
}

PosType SourceBLRSph1::SurfaceBrightness(PosType *y){
  return blr_surface_brightness_spherical_circular_motions(sqrt((pow(y[0]-getX()[0],2) + pow(y[1]-getX()[1],2))),this);
}
PosType SourceBLRSph2::SurfaceBrightness(PosType *y){
  return blr_surface_brightness_spherical_random_motions(sqrt((pow(y[0]-getX()[0],2) + pow(y[1]-getX()[1],2))),this);
}

//void in_source(PosType *y_source,ListHndl sourcelist){
//  return;
//}
SourcePixelled::SourcePixelled(InputParams& params)
{}

SourcePixelled::SourcePixelled(
                               PosType my_z            /// redshift of the source
                               , PosType* my_center  /// center (in rad)
                               , int my_Npixels           /// number of pixels per side
                               , PosType my_resolution  /// resolution (in rad)
                               , PosType* arr_val          /// array of pixel values (must be of size = Npixels*Npixels)
)
:Source(), resolution(my_resolution), Npixels (my_Npixels){
  zsource = my_z;
  
  range = resolution*(Npixels-1);
  
  values.resize(Npixels*Npixels);
  for (int i = 0; i < Npixels*Npixels; i++)
    values[i] = arr_val[i];
  source_r =  range/sqrt(2.);
  source_x[0] = my_center[0];
  source_x[1] = my_center[1];
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
)
:Source(){
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
  for (int i = 0; i < Npixels*Npixels; i++)
    values[i] = gal_map(i)/resolution/resolution*factor;
  
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

PosType SourcePixelled::SurfaceBrightness(PosType *y){
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
                             std::string filter_in  				/// file with the old observing filter
                             , std::string filter_out			/// file with the new observing filter
                             , std::string sed					/// file with the galaxy spectral energy distribution
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
  PosType sed_in = integrateFilterSed(wavel_in, ampl_in, wavel_sed, ampl_sed);
  PosType sed_out = integrateFilterSed(wavel_out, ampl_out, wavel_sed, ampl_sed);
  
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
PosType Source::integrateFilterSed(std::vector<PosType> wavel_fil, std::vector<PosType> fil, std::vector<PosType> wavel_sed, std::vector<PosType> sed)
{
  int wavel_new_size = 10000;
  vector<PosType> wavel_new(wavel_new_size), fil_val(wavel_new_size), sed_val(wavel_new_size);
  
  PosType integr = 0.;
  for (int i = 0; i < wavel_new_size; i++)
  {
    wavel_new[i] = max(wavel_fil[0],wavel_sed[0]) + PosType(i)/PosType(wavel_new_size-1)*(min(wavel_fil[wavel_fil.size()-1],wavel_sed[wavel_sed.size()-1])-max(wavel_fil[0],wavel_sed[0]) );
    vector<PosType>::iterator it_f = lower_bound(wavel_fil.begin(), wavel_fil.end(), wavel_new[i]);
    int p = it_f-wavel_fil.begin()-1;
    vector<PosType>::iterator it_s = lower_bound(wavel_sed.begin(), wavel_sed.end(), wavel_new[i]);
    int q = it_s-wavel_sed.begin()-1;
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
Band SourceShapelets::shape_band[10] = {F435W,F606W,F775W,F850LP,F110W,F160W,SDSS_U,SDSS_G,SDSS_R,SDSS_I};

void SourceShapelets::setActiveBand(Band band)
{
  for (int i = 0; i < 10; i++)
  {
    if (shape_band[i] == band)
    {
      mag = mags[i];
      flux = fluxes[i];
      return;
    }
  }
  std::cout << "The band is not available! Available bands are F435W,F606W,F775W,F850LP,F110W,F160W,SDSS_U,SDSS_G,SDSS_R,SDSS_I" << std::endl;
  exit(1);
}

SourceShapelets::SourceShapelets(
                                 PosType my_z                              /// redshift of the source
                                 , PosType my_mag							/// magnitude
                                 , PosType my_scale						/// scale of the shapelets decomposition
                                 , std::valarray<PosType> my_coeff  	/// coefficients of the shapelets decomposition
                                 , PosType* my_center           			/// center (in rad)
                                 , PosType my_ang					/// rotation angle (in rad)
)
:Source()
{
  zsource = my_z;
  mag = my_mag;
  if(my_center != NULL)
    setX(my_center[0], my_center[1]);
  else
    setX(0, 0);
  ang = my_ang;
  n1 = sqrt(my_coeff.size());
  n2 = n1;
  coeff = my_coeff;
  source_r = my_scale;
  
  flux = pow(10,-0.4*(mag+48.6))*inv_hplanck;
  
  NormalizeFlux();
}

SourceShapelets::SourceShapelets(
                                 PosType my_z							/// redshift of the source
                                 , PosType my_mag						/// magnitude
                                 , std::string shap_file				/// fits file with coefficients in a square array
                                 , PosType* my_center           			/// center (in rad)
                                 , PosType my_ang			/// rotation angle (in rad)
)
:Source()
{
  zsource = my_z;
  mag = my_mag;
  if(my_center != NULL)
    setX(my_center[0], my_center[1]);
  else
    setX(0, 0);
  ang = my_ang;
  
#ifdef ENABLE_FITS
  if(shap_file.empty())
    throw std::invalid_argument("Please enter a valid filename for the FITS file input");
  
  std::auto_ptr<CCfits::FITS> fp(new CCfits::FITS(shap_file.c_str(), CCfits::Read));
  
  CCfits::PHDU& h0 = fp->pHDU();
  
  h0.readKey("BETA", source_r);
  source_r *= 0.03/180./60./60.*pi;
  h0.readKey("DIM", n1);
  h0.readKey("ID", id);
  n2 = n1;
  h0.read(coeff);
  
#else
  std::cerr << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
  exit(1);
#endif
  
  flux = pow(10,-0.4*(mag+48.6))*inv_hplanck;
  
  NormalizeFlux();
}

SourceShapelets::SourceShapelets(
                                 std::string shap_file				/// fits file with coefficients in a square array. Mag and redshift are read from the header.
                                 , PosType* my_center  					/// center (in rad)
                                 , PosType my_ang				 /// rotation angle (in rad)
)
:Source()
{
  if(my_center != NULL)
    setX(my_center[0], my_center[1]);
  else
    setX(0, 0);
  ang = my_ang;
  
#ifdef ENABLE_FITS
  if(shap_file.empty())
    throw std::invalid_argument("Please enter a valid filename for the FITS file input");
  
  std::auto_ptr<CCfits::FITS> fp(new CCfits::FITS(shap_file.c_str(), CCfits::Read));
  
  CCfits::PHDU& h0 = fp->pHDU();
  
  h0.readKey("BETA", source_r);
  source_r *= 0.03/180./60./60.*pi;

  h0.readKey("MAG_B", mags[0]);
  h0.readKey("MAG_V", mags[1]);
  h0.readKey("MAG_I", mags[2]);
  h0.readKey("MAG_Z", mags[3]);
  h0.readKey("MAG_J", mags[4]);
  h0.readKey("MAG_H", mags[5]);
  h0.readKey("MAG_u_KIDS", mags[6]);
  h0.readKey("MAG_g_KIDS", mags[7]);
  h0.readKey("MAG_r_KIDS", mags[8]);
  h0.readKey("MAG_i_KIDS", mags[9]);
  
  for (int i = 0; i < 10; i++)
  {
    if (mags[i] < 0.)
      fluxes[i] = std::numeric_limits<PosType>::epsilon();
    else
      fluxes[i] = pow(10,-0.4*(mags[i]+48.6))*inv_hplanck;
  }
  
  // by default, the magnitude is the one in the i band,
  // whose image has been used for shapelets decomposition
  mag = mags[2];
  flux = fluxes[2];
  
  h0.readKey("REDSHIFT", zsource);
  h0.readKey("ID", id);
  h0.readKey("DIM", n1);
  n2 = n1;
  h0.read(coeff);
#else
  std::cerr << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
  exit(1);
#endif
  
  NormalizeFlux();
}


PosType SourceShapelets::SurfaceBrightness(PosType *y)
{
  PosType sb = 0.;
  PosType y_norm[2],tmp;
  y_norm[0] = ((y[0]-source_x[0])*cos(ang)+(y[1]-source_x[1])*sin(ang))/source_r;
  y_norm[1] = ((y[0]-source_x[0])*sin(ang)-(y[1]-source_x[1])*cos(ang))/source_r;
  //PosType dist = sqrt(y_norm[0]*y_norm[0]+y_norm[1]*y_norm[1]);
  std::vector<PosType> Her1,Her2;
  
  Hermite(Her1,n1,y_norm[0]);
  Hermite(Her2,n2,y_norm[1]);
  
  size_t coei=1,coej=1;
  
  for (int i = 0; i < n1; i++,coei *= 2,coej=1 )
  {
    tmp = factrl(i)*pi;
    for (int j = 0; j < n2; j++,coej *= 2 )
    {
      
      PosType norm = 1./sqrt(coei*coej*tmp*factrl(j));
      sb += norm*coeff[j*n1+i]*Her1[i]*Her2[j];
    }
  }
  sb *= exp(-0.5*(y_norm[0]*y_norm[0]+y_norm[1]*y_norm[1]) )/source_r;
  sb *= flux/coeff_flux;
  
  return max(sb,std::numeric_limits<PosType>::epsilon());
}

/// Returns the value of the Hermite polynomials from degree 0 to n at position x
void SourceShapelets::Hermite(std::vector<PosType> &hg,int N, PosType x)
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
  coeff_flux *= sqrt(pi)*source_r;
  
}

/// Default constructor. Reads in sources from the default catalog. No magnitude limit.
SourceMultiShapelets::SourceMultiShapelets(InputParams& params)
: Source(),index(0)
{
  assignParams(params);
  readCatalog();
}

SourceMultiShapelets::~SourceMultiShapelets()
{
}

/// Reads in the default catalog
void SourceMultiShapelets::readCatalog()
{
  int max_num = 37012;
  for (int i = 0; i < max_num+1; i++)
  {
    std::string shap_file = shapelets_folder+"/obj_"+to_string(i)+"_mag_z.sif";
    std::ifstream shap_input(shap_file.c_str());
    if (shap_input)
    {
      SourceShapelets s(shap_file.c_str());
      s.setActiveBand(band);
      if (s.getMag() > 0. && s.getMag() < mag_limit)
        galaxies.push_back(s);
      shap_input.close();
    }
    else if (i == 1)
    {
      std::cout << "Can't open file " << shap_file << std::endl;
      ERROR_MESSAGE();
      throw std::runtime_error(" Cannot open file.");
      exit(1);
    }
    
  }
  
}


void SourceMultiShapelets::assignParams(InputParams& params){
  if(!params.get("source_mag_limit",mag_limit)){
    std::cout << "ERROR: Must assign source_mag_limit in parameter file " << params.filename() << std::endl;
    exit(1);
  }
  
  
  if(!params.get("source_sb_limit",sb_limit))
    setSBlimit_magarcsec(30.);
  else
    sb_limit = pow(10,-0.4*(48.6+sb_limit))*pow(180*60*60/pi,2)/hplanck;
  
  if(!params.get("shapelets_folder",shapelets_folder)){
    std::cout << "ERROR: shapelets_folder not found in parameter file " << params.filename() << std::endl;
    exit(1);
  }
  
  if(!params.get("shapelets_band",band)){
    std::cout << "ERROR: Must assign shapelets_band in parameter file " << params.filename() << std::endl;
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



