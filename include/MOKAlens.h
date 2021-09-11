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
	std::valarray<double> surface_density;  // Msun / Mpc^2
	std::valarray<double> alpha1_bar;
	std::valarray<double> alpha2_bar;
	std::valarray<double> gamma1_bar;
	std::valarray<double> gamma2_bar;
	//std::valarray<double> gamma3_bar;
  std::valarray<double> phi_bar;
	//std::vector<double> x;
  int nx,ny;
  // boxlMpc is Mpc/h for MOKA
	/// lens and source properties
  double zlens,m,zsource,Dlens,DLS,DS,c,cS,fsub,mstar,minsubmass;
  /// range in x direction, pixels are square
  double boxlarcsec,boxlMpc,boxlrad;
  /// cosmology
  double omegam,omegal,h,wq;
	double inarcsec;
  
  MOKAmap(std::string MOKA_input_file,bool zeromean,const COSMOLOGY &cosmo){
    read(MOKA_input_file,zeromean,cosmo);
  }
  
  void read(std::string MOKA_input_file,bool zeromean,const COSMOLOGY &cosmo);
  void write(std::string filename);

	Point_2d center;
  
  
  MOKAmap(){}
  MOKAmap(const MOKAmap &map){
      surface_density = map.surface_density;
      alpha1_bar = map.alpha1_bar;
      alpha2_bar = map.alpha2_bar;
      gamma1_bar = map.gamma1_bar;
      gamma2_bar = map.gamma2_bar;
      //gamma3_bar = map.gamma3;
      phi_bar = map.phi_bar;
      //Signlambdar = map.Signlambdar;
      //Signlambdat = map.Signlambdat;
      //x = map.x;
      nx = map.nx;
      ny = map.ny;
      zlens = map.zlens;
      m = map.m;
      zsource = map.zsource;
      Dlens = map.Dlens;
      DLS = map.DLS;
      DS = map.DS;
      c = map.c;
      cS = map.cS;
      fsub = map.fsub;
      mstar = map.mstar;
      minsubmass = map.minsubmass;
      boxlarcsec = map.boxlarcsec;
      boxlMpc = map.boxlMpc;
      boxlrad = map.boxlrad;
      omegam = map.omegam;
      omegal = map.omegal;
      h = map.h;
      wq = map.wq;
      inarcsec = map.inarcsec;
      center = map.center;
  }

  MOKAmap(MOKAmap &&map){
    *this = std::move(map);
  }
  
  MOKAmap & operator=(MOKAmap &&map){
    if(this != &map){
      surface_density = std::move(map.surface_density);
      alpha1_bar = std::move(map.alpha1_bar);
      alpha2_bar = std::move(map.alpha2_bar);
      gamma1_bar = std::move(map.gamma1_bar);
      gamma2_bar = std::move(map.gamma2_bar);
      //gamma3_bar = std::move(map.gamma3_bar);
      phi_bar = std::move(map.phi_bar);
      //Signlambdar = std::move(map.Signlambdar);
      //Signlambdat = std::move(map.Signlambdat);
      //x = std::move(map.x);
      nx = map.nx;
      ny = map.ny;
      zlens = map.zlens;
      m = map.m;
      zsource = map.zsource;
      Dlens = map.Dlens;
      DLS = map.DLS;
      DS = map.DS;
      c = map.c;
      cS = map.cS;
      fsub = map.fsub;
      mstar = map.mstar;
      minsubmass = map.minsubmass;
      boxlarcsec = map.boxlarcsec;
      boxlMpc = map.boxlMpc;
      boxlrad = map.boxlrad;
      omegam = map.omegam;
      omegal = map.omegal;
      h = map.h;
      wq = map.wq;
      inarcsec = map.inarcsec;
      center = map.center;
    }
    
    return *this;
  }
  MOKAmap & operator=(const MOKAmap &map){
    if(this != &map){
      surface_density = map.surface_density;
      alpha1_bar = map.alpha1_bar;
      alpha2_bar = map.alpha2_bar;
      gamma1_bar = map.gamma1_bar;
      gamma2_bar = map.gamma2_bar;
      //gamma3_bar = map.gamma3_bar;
      phi_bar = map.phi_bar;
      //Signlambdar = map.Signlambdar;
      //Signlambdat = map.Signlambdat;
      //x = map.x;
      nx = map.nx;
      ny = map.ny;
      zlens = map.zlens;
      m = map.m;
      zsource = map.zsource;
      Dlens = map.Dlens;
      DLS = map.DLS;
      DS = map.DS;
      c = map.c;
      cS = map.cS;
      fsub = map.fsub;
      mstar = map.mstar;
      minsubmass = map.minsubmass;
      boxlarcsec = map.boxlarcsec;
      boxlMpc = map.boxlMpc;
      boxlrad = map.boxlrad;
      omegam = map.omegam;
      omegal = map.omegal;
      h = map.h;
      wq = map.wq;
      inarcsec = map.inarcsec;
      center = map.center;
    }
    
    return *this;
  }
  
  void PreProcessFFTWMap(float zerosize);

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
                  ,PixelMapType maptype
                  ,int pixel_map_zeropad
                  ,bool my_zeromean
                  ,const COSMOLOGY& lenscosmo
                  );
  
  //LensHaloMassMap(PixelMap &map,double massconvertion,double zlens,double zsource,int pixel_map_zeropad,const COSMOLOGY& lenscosmo);
  
	//LensHaloMassMap(InputParams& params, COSMOLOGY& lenscosmo);
	
  LensHaloMassMap(
                  const PixelMap &MassMap   /// mass map in solar mass units
                  ,double massconvertion    /// convertion factor from pixel units to solar masses
                  ,double redshift          /// redshift of lens
                  ,int pixel_map_zeropad    /// factor by which to zero pad in FFTs, ex. 4
                  ,bool my_zeromean         /// if true, subtracts average density
                  ,const COSMOLOGY& lenscosmo  /// cosmology
  );
  
  /// This makes a uniform rectangular mass sheat
  LensHaloMassMap(double mass           /// total mass in rectangle
                  ,Point_2d center      /// center of rectangle
                  ,Point_2d range       /// width and hight in radians
                  ,double resolution    /// resolution in radians
                  ,int zeropadding      /// factor by which to zero pad in FFTs, ex. 1 is no padding, 2 FFT grid is twice as big as original map
                  ,double redshift      /// redshift of plane
                  ,const COSMOLOGY &cosmo
                  );
  
  LensHaloMassMap(const LensHaloMassMap &h):LensHalo(h),cosmo(h.cosmo){
    maptype = h.maptype;
    map = h.map;
    zerosize = h.zerosize;
    zeromean = h.zeromean;
  }
 LensHaloMassMap(LensHaloMassMap &&h):cosmo(h.cosmo){
    *this = std::move(h);
  }
  
  LensHaloMassMap & operator=(LensHaloMassMap &h){
    if(&h != this){
      LensHalo::operator=(h);
      //cosmo = h.cosmo;
      maptype = h.maptype;
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
      maptype = h.maptype;
      map = std::move(h.map);
      zerosize = h.zerosize;
      zeromean = h.zeromean;
    }
    return *this;
  }

	~LensHaloMassMap();
	
	std::string MOKA_input_file;
	/// if >=1 (true), do analyzis only; if = 0 (false) change units to internal GLAMER units and prepare for ray-shooting
	int flag_MOKA_analyze;
	int flag_background_field;
	
	void assignParams(InputParams& params);

	void checkCosmology();
	
	void saveImage(bool saveprofile=true);
	//void saveKappaProfile();
	//void saveGammaProfile();
	//void saveProfiles(double &RE3, double &xxc, double &yyc);
    
	void force_halo(double *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,double const *xcm,bool subtract_point=false,PosType screening = 1.0);
    
	//void saveImage(GridHndl grid,bool saveprofiles);
	
	//void estSignLambdas();
	//void EinsteinRadii(double &RE1, double &RE2, double &xxc, double &yyc);
	
	void getDims();
	void readMap();
  // set up map from a PixelMap of the surface density
  void setMap(const PixelMap &inputmap,double convertionfactor,double z);

	void writeImage(std::string fn);
  
	/// return center in physical Mpc
	Point_2d getCenter() const { return map.center; }

	/// return range of input map in rad
	double getRangeRad() const { return map.boxlrad; }
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
	PixelMapType maptype;
	void initMap();

	void convertmap(MOKAmap &map,PixelMapType maptype);
	MOKAmap map;
	const COSMOLOGY& cosmo;
  int zerosize;
  bool zeromean;
  
protected:
  LensHaloMassMap(COSMOLOGY &c):cosmo(c){}
};

void make_friendship(int ii,int ji,int np,std:: vector<int> &friends, std:: vector<double> &pointdist);

int fof(double l,std:: vector<double> xci, std:: vector<double> yci, std:: vector<int> &groupid);
#endif
/* MOKALENS_H_ */



