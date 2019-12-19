//
//  lensmap.h
//  GLAMER
//
//  Created by bmetcalf on 18.12.19.
//

#ifndef lensmap_h
#define lensmap_h

#include <valarray>
#include "point.h"
#include "cpfits.h"

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
struct LensMap{
  
  LensMap():nx(0),ny(0),boxlMpc(0),angular_pixel_size(0){};
  
  // move operators
  LensMap(LensMap &&m);
  LensMap& operator=(LensMap &&m);
  LensMap& operator=(const LensMap &m);

  /// values for the map
  std::valarray<double> surface_density;  // Msun / Mpc^2
  std::valarray<float> alpha1_bar;  // Msun / Mpc
  std::valarray<float> alpha2_bar;  // Msun / Mpc
  std::valarray<float> gamma1_bar;  // Msun / Mpc^2
  std::valarray<float> gamma2_bar;  // Msun / Mpc^2
  std::valarray<float> phi_bar;     // Msun
  int nx,ny;
  
  double boxlMpc;
  double angular_pixel_size;  // in radians
  Point_2d center;
  Point_2d lowerleft;  /// boundery with centred grid
  Point_2d upperright; ///
  
  double x_resolution(){return boxlMpc / nx ;}
  double y_resolution(){return (upperright[1]-lowerleft[1])/ny;}
  // # of pixels times resolution
  double x_range(){return boxlMpc;}
  // # of pixels times resolution
  double y_range(){return (upperright[1]-lowerleft[1]);}
  
  bool evaluate(const double *x,float &sigma,float *gamma,double *alpha);
  
  LensMap(std::string fits_input_file,double angDist){
    read(fits_input_file,angDist);
  }
  
  /// read an entire map
  void read(std::string input_fits,double angDist);//,float h,float z);
  /// read from a file that has been generated with LensMap::write()
  void Myread(std::string fits_input_file);
  
  
  /// read only header information
  //void read_header(std::string input_fits,float h,float z);
  void read_header(std::string input_fits,double angDist);
  
  /// read a subsection of the fits map
  //  void read_sub(std::string input_fits
  //                ,const std::vector<long> &first
  //                ,const std::vector<long> &last
  //                ,double Dist
  //                );
  
  void read_sub(CPFITS_READ &cpfits
                ,std::vector<long> &first
                ,std::vector<long> &last
                ,double Dist
                );
  
  void write(std::string filename);
  /// meant to output directly in angulare units and lensing quantities
  void write(std::string filename,LensingVariable quant);
  
  // this calculates the other lensing quantities from the density map
  
  template <class T>
  void PreProcessFFTWMap(float zerosize,T Wphi_of_k,bool do_alpha = true);
  template <class T>
  void PreProcessFFTWMap(T Wphi_of_k,bool do_alpha = true);
  
  
  struct UNIT{
    int operator()(float k2){return 1;}
  };
  struct WLR{
    float rs2;
    float operator()(float k2){return exp(-k2*rs2);}
  };
  struct WSR{
    float rs2;
    float operator()(float k2){return 1 - exp(-k2*rs2);}
  };
};

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

#endif /* lensmap_h */
