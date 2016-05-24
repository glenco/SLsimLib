/*
 * MOKAlens.h
 *
 *  Created on: Jun 8, 2012
 *      Author: mpetkova
 */


#ifndef MOKALENS_H_
#define MOKALENS_H_

#include "standard.h"
#include "profile.h"
#include "InputParams.h"
#include "lens_halos.h"
#include "grid_maintenance.h"

#include <stdexcept>

/**
 * \brief The MOKA map structure, containing all quantities that define it
 *
 * The MOKA map, that is read in from the fits file. Its components include the
 * lensing properties, as well as the cosmology, the size of the field of view,
 * the redshifts of the lens and source, the properties of the cluster, etc.
 *
 * Note: To use this class requires setting the ENABLE_FITS compiler flag and linking
 * the cfits library.
 */
struct MOKAmap{
	/// values for the map
	std::valarray<double> convergence;
	std::valarray<double> alpha1;
	std::valarray<double> alpha2;
	std::valarray<double> gamma1;
	std::valarray<double> gamma2;
	std::valarray<double> gamma3;
  std::valarray<double> phi;
	std::valarray<double> Signlambdar;
	std::valarray<double> Signlambdat;
	std:: vector<double> x;
  int nx,ny;
  // boxlMpc is Mpc/h for MOKA
	/// lens and source properties
  double zlens,m,zsource,Dlens,DLS,DS,c,cS,fsub,mstar,minsubmass;
  /// range in x direction, pixels are square
  double boxlarcsec,boxlMpc,boxlrad;
  /// cosmology
  double omegam,omegal,h,wq;
	double inarcsec;
	double center[2];
};

/**
 *  \brief A class that includes the MOKA lens map
 *
 * A class, where the lens is represented by a MOKA cluster in the form of a
 * MOKAmap object (see the description of MOKAmap). It can either be used with the
 * Lens model or by itself.
 *
 * Note: To use this class requires setting the ENABLE_FITS compiler flag and linking
 * the cfits library.
 */
class LensHaloMassMap : public LensHalo
{
public:
	LensHaloMassMap(const std::string& filename
                  , PixelMapType maptype
                  ,int pixel_map_zeropad
                  ,bool my_zeromean
                  , const COSMOLOGY& lenscosmo
                  );
  
  LensHaloMassMap(InputParams& params, const COSMOLOGY& lenscosmo);

  LensHaloMassMap(
                  const PixelMap &MassMap   /// mass map in solar mass units
                  ,double redshift          /// redshift of lens
                  ,int pixel_map_zeropad    /// factor by which to zero pad in FFTs, ex. 4
                  ,bool my_zeromean         /// if true, subtracts average density
                  ,const COSMOLOGY& lenscosmo  /// cosmology
  );


	~LensHaloMassMap();
	
	std::string MOKA_input_file;
	/// if >=1 (true), do analyzis only; if = 0 (false) change units to internal GLAMER units and prepare for ray-shooting
	int flag_MOKA_analyze;
	int flag_background_field;
	
	void assignParams(InputParams& params);

	void checkCosmology();
	
	void saveImage(bool saveprofile=true);
	void saveKappaProfile();
	void saveGammaProfile();
	void saveProfiles(double &RE3, double &xxc, double &yyc);
    
	void force_halo(double *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,double const *xcm,bool subtract_point=false,PosType screening = 1.0);
    
	void saveImage(GridHndl grid,bool saveprofiles);
	
	void estSignLambdas();
	void EinsteinRadii(double &RE1, double &RE2, double &xxc, double &yyc);
	
	void getDims();
	void readMap();
  // set up map from a PixelMap of the surface density
  void convertMap(const PixelMap &inputmap,double z);

	void writeImage(std::string fn);
  
	/// return center in physical Mpc
	const double* getCenter() const { return map->center; }
	/// return range of input map in rad
	double getRangeRad() const { return map->boxlrad; }
	/// return range of input map in physical Mpc
	double getRangeMpc() const { return map->boxlMpc; }
	/// /// return number of pixels on a side in original map
	size_t getN() const
	{
		if(map->nx != map->ny)
			throw std::runtime_error("mass map is not square");
		return map->nx;
	}
	/// return number of pixels on a x-axis side in original map
	size_t getNx() const { return map->nx; }
	/// return number of pixels on a y-axis side in original map
	size_t getNy() const { return map->ny; }
	
private:
	PixelMapType maptype;
	void initMap();
	void convertmap(MOKAmap *map,PixelMapType maptype);
	MOKAmap* map;
	const COSMOLOGY& cosmo;
	void PreProcessFFTWMap();
  int zerosize;
  bool zeromean;
};
  


void make_friendship(int ii,int ji,int np,std:: vector<int> &friends, std:: vector<double> &pointdist);

int fof(double l,std:: vector<double> xci, std:: vector<double> yci, std:: vector<int> &groupid);
#endif /* MOKALENS_H_ */

/// constructor for making lens halo directly from a mass map
LensHaloMassMap::LensHaloMassMap(
                                 const PixelMap &MassMap   /// mass map in solar mass units
                                 ,double redshift          /// redshift of lens
                                 ,int pixel_map_zeropad    /// factor by which to zero pad in FFTs, ex. 4
                                 ,bool my_zeromean         /// if true, subtracts average density
                                 ,const COSMOLOGY& lenscosmo  /// cosmology
)
:LensHalo()
, flag_MOKA_analyze(0), flag_background_field(0),maptype(pix_map),cosmo(lenscosmo),zerosize(pixel_map_zeropad),zeromean(my_zeromean)
{
  convertMap(MassMap, redshift);
  
  posHalo[0] = MassMap.getCenter()[0];
  posHalo[1] = MassMap.getCenter()[1];
 	
  // set redshift to value from map
  setZlens(map->zlens);
}

/**
 * \brief reads in the fits file for the MOKA or mass map and saves it in the structure map
 */
void LensHaloMassMap::convertMap(
                                 const PixelMap &inputmap  // in units of solar masses
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
    map->convergence[i] = inputmap(i)/pixelarea;
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

