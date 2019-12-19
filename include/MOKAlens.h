/*
 * MOKAlens.h
 *
 */


#ifndef MOKALENS_H_
#define MOKALENS_H_

#include "standard.h"
#include "profile.h"
#include "InputParams.h"
#include "lens_halos.h"
#include "grid_maintenance.h"
#include "multimap.h"
#include "lensmap.h"

#include <stdexcept>


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
                  ,double redshift
                  ,double mass_unit       /// factor the pixels should be multiplied by to get solar masses
                  ,int pixel_map_zeropad
                  ,bool my_zeromean
                  ,COSMOLOGY& lenscosmo
                  );
  
  //LensHaloMassMap(PixelMap &map,double massconvertion,double zlens,double zsource,int pixel_map_zeropad,const COSMOLOGY& lenscosmo);
	LensHaloMassMap(InputParams& params, COSMOLOGY& lenscosmo);
	
  LensHaloMassMap(
                  const PixelMap &MassMap   /// mass map in solar mass units
                  ,double massconvertion    /// convertion factor from pixel units to solar masses
                  ,double redshift          /// redshift of lens
                  ,int pixel_map_zeropad    /// factor by which to zero pad in FFTs, ex. 4
                  ,bool my_zeromean         /// if true, subtracts average density
                  ,COSMOLOGY& lenscosmo  /// cosmology
  );
  
  LensHaloMassMap(const LensHaloMassMap &h):LensHalo(h),cosmo(h.cosmo){
//    LensHaloMassMap(const LensHaloMassMap &h){
//    maptype = h.maptype;
    map = h.map;
    zerosize = h.zerosize;
    zeromean = h.zeromean;
  }
  LensHaloMassMap(LensHaloMassMap &&h):LensHalo(std::move(h)),cosmo(h.cosmo){
    map = std::move(h.map);
    zerosize = h.zerosize;
    zeromean = h.zeromean;
  }
  
  LensHaloMassMap & operator=(LensHaloMassMap &h){
    if(&h != this){
      LensHalo::operator=(h);
      //cosmo = h.cosmo;
      //maptype = h.maptype;
      map = h.map;
      zerosize = h.zerosize;
      zeromean = h.zeromean;
    }
    return *this;
  }
  LensHaloMassMap & operator=(LensHaloMassMap &&h){
    if(&h != this){
      LensHalo::operator=(h);
      //cosmo = h.cosmo;
      //maptype = h.maptype;
      map = std::move(h.map);
      zerosize = h.zerosize;
      zeromean = h.zeromean;
    }
    return *this;
  }

	~LensHaloMassMap();
	
	std::string MOKA_input_file;
  
	void assignParams(InputParams& params);

	//void checkCosmology();
	
	void saveImage(bool saveprofile=true);
	//void saveKappaProfile();
	//void saveGammaProfile();
	//void saveProfiles(double &RE3, double &xxc, double &yyc);
    
	void force_halo(double *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,double const *xcm,bool subtract_point=false,PosType screening = 1.0);
    
	//void saveImage(GridHndl grid,bool saveprofiles);
	
	//void estSignLambdas();
	//void EinsteinRadii(double &RE1, double &RE2, double &xxc, double &yyc);
	
  std::vector<long> getDims();
  // set up map from a PixelMap of the surface density
  void setMap(const PixelMap &inputmap,double convertionfactor,double z);

	void writeImage(std::string fn);
  
	/// return center in physical Mpc
	Point_2d getCenter() const { return map.center; }

	/// return range of input map in rad
	double getRangeRad() const { return map.boxlMpc / cosmo.angDist( getZlens() ); }
	/// return range of input map in physical Mpc
	double getRangeMpc() const { return map.boxlMpc; }
	/// /// return number of pixels on a side in original map
	size_t getN() const
	{
		if(map.nx != map.ny)
			throw std::runtime_error("mass map is not square");
		return map.nx;
	}
	/// return number of pixels on a x-axis side in original map
	size_t getNx() const { return map.nx; }
	/// return number of pixels on a y-axis side in original map
	size_t getNy() const { return map.ny; }
  
private:
	//PixelMapType maptype;
	void initMap();
  double mass_unit;

	void convertmap(MOKAmap &map,PixelMapType maptype);
  //MOKAmap map;
  LensMap map;
	const COSMOLOGY& cosmo;
  int zerosize;
  bool zeromean;
};

void make_friendship(int ii,int ji,int np,std:: vector<int> &friends, std:: vector<double> &pointdist);

int fof(double l,std:: vector<double> xci, std:: vector<double> yci, std:: vector<int> &groupid);
#endif
/* MOKALENS_H_ */



