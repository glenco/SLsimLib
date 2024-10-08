/*
 * source.h
 *
 */

#ifndef SOURCE_H_
#define SOURCE_H_

#include <vector>
//#include "standard.h"
#include "utilities_slsim.h"
#include "InputParams.h"
#include "image_processing.h"
#include "point.h"

double mag_to_jansky_AB(double);
double jansky_to_mag_AB(double flux);

/// ergs / s /cm^2 / Hz
double mag_to_flux_AB(double);
/// ergs / s /cm^2 / Hz
double flux_to_mag_AB(double flux);

template <typename T>
T flux_to_mag(T flux,T zeropoint){
  if(flux <=0) return 100;
  return -2.5 * log10(flux) + zeropoint;
}

template <typename T>
T mag_to_flux(T m,T zeropoint){
  if(m == 100) return 0;
  return pow(10,-0.4*(m - zeropoint));
}

//double mag_to_flux(double m,double zeropoint);
//double flux_to_mag(double flux,double zeropoint);


/** \brief Base class for all sources.
 *
 */
class Source
{

public:
  /// shell constructor
	Source(PosType r
         ,Point_2d x
         ,PosType z
         ,PosType SBlimit  // surface brightness limit in mags per arcsecond
         ,PosType zero_point // the magnitude zero point
         ):source_r(r),source_x(x),zsource(z),sb_limit(SBlimit),mag_zero_point(zero_point)
  {
    //setSBlimit_magarcsec(SBlimit);
    id=-1;
  };
  
  virtual ~Source();

  Source(const Source &s):
    source_r(s.source_r),
    source_x(s.source_x),
    zsource(s.zsource),
    sb_limit(s.sb_limit),
    id(s.id),
    mag_zero_point(s.mag_zero_point){ }
  
  Source & operator=(const Source &s){
    if(this == &s) return *this;
    source_r = s.source_r;
    source_x = s.source_x;
    zsource = s.zsource;
    sb_limit = s.sb_limit;
    id = s.id;
    mag_zero_point = s.mag_zero_point;
    
    return *this;
  }

	// in lens.cpp
	/** Surface brightness of source in grid coordinates not source centered coordinates.
   
   * The units shuld be ergs / s / Hz / cm^2 
   */
  virtual PosType SurfaceBrightness(const PosType *y) const = 0;
  virtual PosType getTotalFlux() const = 0;
  virtual void printSource() = 0;

  /// convert mag/arcsec^2 to flux units
  double SBlimit_magarcsec(double limit) {return mag_to_flux(limit,mag_zero_point)*pow(180*60*60/PI,2);}
 
	/// Gets sb_limit in erg/cm^2/sec/rad^2/Hz
	PosType getSBlimit(){return sb_limit;}
	/// Gets sb_limit in mag/arcsec^2
  PosType getSBlimit_magarcsec(){return SBlimit_magarcsec(sb_limit);}

	// accessor functions that will sometimes be over ridden in class derivatives
	/// Redshift of source
	PosType getZ() const {return zsource;}
	void setZ(PosType my_z){zsource = my_z;}
	/// Radius of source in radians
	inline PosType getRadius() const {return source_r;}
  /// Reset the radius of the source in radians
	void setRadius(PosType my_radius){source_r = my_radius;}
  /// position of source in radians
  inline Point_2d getTheta() const {return source_x;}
  /// position of source in radians
  inline void getTheta(PosType *x) const {x[0] = source_x.x[0]; x[1] = source_x.x[0];}
  /// position of source in radians
  inline void getTheta(Point_2d &x) const {x = source_x;}
  
  /// Reset the position of the source in radians
	void setTheta(PosType *xx){source_x[0] = xx[0]; source_x[1] = xx[1];}
  void setTheta(PosType my_x,PosType my_y){source_x[0] = my_x; source_x[1] = my_y;}
  void setTheta(const Point_2d &p){source_x = p;}

	/// Sets sb_limit in erg/cm^2/sec/rad^2/Hz
	void setSBlimit(float limit) {sb_limit = limit;}
 
  void setMagZeroPoint(float zeropoint){mag_zero_point=zeropoint;}
  double getMagZeroPoint(){return mag_zero_point;}

	PosType changeFilter(std::string filter_in, std::string filter_out, std::string sed);
	PosType integrateFilter(std::vector<PosType> wavel_fil, std::vector<PosType> fil);
	PosType integrateFilterSED(std::vector<PosType> wavel_fil, std::vector<PosType> fil, std::vector<PosType> wavel_sed, std::vector<PosType> sed);

  static PosType *getx(Source &source){return source.source_x.x;}
  
  /// test if flux in pixels matches total flux
  double TEST_surface_brightness(double res,int N){
    
    Point_2d x = source_x;
    x[0] -= res*N/2.;
    x[1] -= res*N/2.;
    double total_flux = 0;
    
    for(size_t i=0 ; i<N ; ++i){
      for(size_t j=0 ; j<N ; ++j){
        total_flux += SurfaceBrightness(x.x);
        x[1] += res;
      }
      x[0] += res;
      x[1] = source_x[1] - res*N/2.;
    }
    total_flux *= res*res;
    
    std::cout << "------- check flux ----------" << std::endl;
    std::cout << " total flux in pixels : " << total_flux << std::endl;
    std::cout << " getTotalFlux() : " << getTotalFlux() << std::endl;
    std::cout << " fractional difference : " << (total_flux - getTotalFlux()) / getTotalFlux() << std::endl;
    std::cout << " magnitude from flux in pixels : " << flux_to_mag(total_flux) << std::endl;
    std::cout << " mag from getTotalFlux() : " <<  flux_to_mag(getTotalFlux()) << std::endl;
    std::cout << "-----------------------------" << std::endl;
     
    if(abs( (total_flux - getTotalFlux()) / getTotalFlux() ) > 0.1 ){
      //assert( abs( (total_flux - getTotalFlux()) / getTotalFlux() ) < 0.1 );
    }
    
    return total_flux;
  }
  
  long getID() const {return id;}
  void setID(long i){id=i;}
  
protected:
  virtual void assignParams(InputParams& params){};
	
	// source parameters
	/// charactoristic source size
	PosType source_r;
	/// center of source
	Point_2d source_x;
	
	/// redshift of source
	PosType zsource;
	PosType sb_limit;
  
  double flux_to_mag(double flux) const {
    if(flux <=0) return 100;
    return -2.5 * log10(flux) + mag_zero_point;
  }
  
//  double mag_to_flux(double mag) const {
//    if(mag == 100) return 0;
//    return pow(10,-0.4*(mag - mag_zero_point));
//  }
  
  long id;
  double mag_zero_point;
};

class SourceColored : public Source{
public:
  SourceColored(PosType magnitude,PosType r,Point_2d x,PosType z
                ,PosType SBlimit,Band band)
  :Source(r,x,z,SBlimit,zeropoints.at(band))
  {
    if(zeropoints.find(band) == zeropoints.end()){
      std::cerr << "SourceColored band " << to_string(band) << " zeropoint needs to be set." << std::endl;
      throw std::runtime_error("unset zeropoint");
    }
    //zeropoints[Band::NoBand] = zero_point;
    current_band = band;
    
    mag_map[current_band]=magnitude;
    mag = magnitude;
    flux_total = mag_to_flux(mag,zeropoints.at(band));
    setID(-1);
  }
  
  SourceColored(PosType magnitude,PosType r,Point_2d x,PosType z
                ,PosType SBlimit,double zero_point,Band band)
  :Source(r,x,z,SBlimit,zero_point)
  {
    zeropoints[band] = zero_point;
    current_band = band;
    
    mag_map[current_band]=magnitude;
    mag = magnitude;
    flux_total = mag_to_flux(mag,zero_point);
    setID(-1);
  }
  
  SourceColored(const SourceColored &s):Source(s){
    mag = s.mag;
    mag_map = s.mag_map;
    current_band = s.current_band;
    sed_type = s.sed_type;
    flux_total = s.flux_total;
  }

  SourceColored & operator= (const SourceColored &s){
    if(this == &s) return *this;
    
    Source::operator=(s);
    mag = s.mag;
    mag_map = s.mag_map;
    current_band = s.current_band;
    sed_type = s.sed_type;
    flux_total = s.flux_total;
    
    return *this;
  }
  
  /// this copies only the Source and SourceColored parts
  SourceColored & copy_color(const SourceColored &s){
    if(this == &s) return *this;
    
    Source::operator=(s);
    mag = s.mag;
    mag_map = s.mag_map;
    current_band = s.current_band;
    sed_type = s.sed_type;
    flux_total = s.flux_total;
    
    return *this;
  }
  
  virtual ~SourceColored(){};
  
  PosType getMag() const { assert(current_band != Band::NoBand); return mag;}
  PosType getMag(Band band) const {return (mag_map.size() > 0) ? mag_map.at(band) : mag;}
  Band getBand() const{return current_band;}
  float getSEDtype() const {return sed_type;}
  double getMagZeroPoint(Band band) const {return zeropoints.at(band);}
  
  void setSEDtype(float sed) {sed_type = sed;}
  static void setMagZeroPoint(Band band,double zeropoint){zeropoints[band]=zeropoint;}
  static void setMagZeroPoints(std::map<Band,PosType> &zero_points){
    for(auto &p : zero_points){
      zeropoints[p.first] = p.second;
    }
  }
 
  void setActiveBand(Band band);
  
  inline PosType getTotalFlux() const {return flux_total;}
  void setMag(float magnitude,Band band,double zeropoint){
    mag_map[band]=magnitude;
    zeropoints[band]=zeropoint;
    if(band==current_band) flux_total = mag_to_flux(mag,zeropoint);
  }
  
  /// this sets the magnitude in a band without changing the current band
  void setMag(float magnitude,Band band){
    mag_map[band]=magnitude;
    if(band==current_band) flux_total = mag_to_flux(mag,zeropoints.at(band));
  }
  
  // rotate on the sky
  virtual void rotate(PosType theta) = 0;

  /// shift all magnitudes my delta_mag and update flux_total, keeps color
  void shiftmag(float delta_mag){
    mag += delta_mag;
    flux_total = mag_to_flux(mag,zeropoints[current_band]);
    for(auto &b : mag_map){
      b.second += delta_mag;
    }
  }
  
  
protected:
  static std::map<Band,PosType> zeropoints;
  
  Band current_band;
  float sed_type;
  double mag;
  double flux_total;
private:
  std::map<Band,PosType> mag_map;
};



/** \brief Class for sources described by an array of pixels
 *
 *  The sources are created as a square array of pixels of dimension Npixels x Npixels and pixel size equal to resolution.
 *
 */
class SourcePixelled: public Source{
public:
	SourcePixelled(PosType my_z, PosType* center, int Npixels, PosType resolution, PosType* arr_val,PosType zero_point);
  template <typename T>
	SourcePixelled(const PixelMap<T>& gal_map, PosType z, PosType factor, PosType zero_point);
	//SourcePixelled(InputParams& params);
  
	~SourcePixelled();
	PosType SurfaceBrightness(const PosType *y) const;
	void printSource();
	inline PosType getTotalFlux() const {return flux;}
	inline PosType getRadius() const {return source_r;}
	inline PosType* getEll(){return ell;};
	inline PosType getQuad(int i, int j){return quad[i][j];};
	inline PosType getSize(){return size;};
	inline PosType* getCentroid(){return centroid;};
	inline PosType getMag(){return flux_to_mag(flux);};
private:
	void assignParams(InputParams& params);
	void calcEll();
	void calcSize();
	void calcCentroid();
	void calcTotalFlux();
	PosType resolution;
	PosType range;
	long Npixels;
	PosType flux;
	PosType quad[2][2];
	PosType ell[2];
	PosType size;
	PosType centroid[2];
	std::valarray<PosType> values;
};

class SourceMultiShapelets;

/** \brief Class for sources described by shapelets.
 *
 *  The sources are created from magnitude, scale radius, and the coefficients of their decomposition into the shapelets basis functions (Refregier et al., 2001).
 *  The coefficients can be read from a fits square array.
 *  In case the magnitude is not given as input, the constructor will read an array of values from the shapelet file header. One can then choose the desired magnitude with setActiveMag.
 *
 */
class SourceShapelets: public SourceColored{
public:
  
  friend SourceMultiShapelets;
    //SourceShapelets();
	SourceShapelets(PosType my_z, PosType my_mag, PosType my_scale, std::valarray<PosType> my_coeff, PosType* my_center, PosType my_ang,PosType zero_point,Band band);
	SourceShapelets(PosType my_z, PosType my_mag, std::string shap_file, PosType *my_center, PosType my_ang, PosType zero_point,Band band);
  SourceShapelets(std::string shap_file, PosType my_ang, PosType zero_point,Band band);
  
  ~SourceShapelets(){--count;}
  
  SourceShapelets(const SourceShapelets &s):SourceColored(s){
    coeff = s.coeff;
    n1 = s.n1;
    n2 = s.n2;
   //ang = s.ang;
    cos_sin[0] = s.cos_sin[0];
    cos_sin[1] = s.cos_sin[1];
    coeff_flux = s.coeff_flux;
    ++count;
  }

  SourceShapelets & operator= (const SourceShapelets &s){
    if(this == &s) return *this;
    
    SourceColored::operator=(s);
    coeff = s.coeff;
    n1 = s.n1;
    n2 = s.n2;
    //ang = s.ang;
    cos_sin[0] = s.cos_sin[0];
    cos_sin[1] = s.cos_sin[1];
    coeff_flux = s.coeff_flux;
    
    return *this;
  }
  
  /// rotate source
  void rotate(PosType ang  /// angle in radians
              ){
    cos_sin[0] = cos(ang);
    cos_sin[1] = sin(ang);
  }

	PosType SurfaceBrightness(const PosType *y) const;
	void printSource();
  /// maximum size
	inline PosType getRadius() const {return source_r*10.;}
  
  void pepper(int n,double s,Utilities::RandomNumbers_NR &ran){
    
    double f = SurfaceBrightness(source_x.x);
    
    double t=0.9,t2=0.3;
    
    for(int i=0 ; i < n ; ++i){
      double r = source_r*sqrt(ran());
      double theta = 2*PI*ran();
      double ff = f*( t - t2 * ran() );
      
      Point_2d x( r*cos(theta) , r*sin(theta) );
  
      grains.emplace_back(ff,s,x);
    }
  }
  
  double sub_structure_sb(double *x){
    
    double sum = 0;
    for(auto &g : grains){
      sum += g(x);
    }
    
    return sum;
  }

  class PepperCorn{
  public:
    PepperCorn(double f,double sig,Point_2d &x):
    x(x),sigma(sig),fo(f){};
    
    double operator()(double *y){
      return fo * exp(-( (x[0]-y[0])*(x[0]-y[0]) - (x[1]-y[1])*(x[1]-y[1]) ) /sigma/sigma);
    }
    
  private:
    Point_2d x;
    double sigma;
    double fo;
  };
  
private:
  
 
	void assignParams(InputParams& params);
  void Hermite(std::vector<PosType> &hg,int N, PosType x) const;

	void NormalizeFlux();
	std::valarray<PosType> coeff;
	int n1,n2;
 
	//PosType ang;
  PosType cos_sin[2];  // [0] is cos(ang), [1] is sin(ang)
  PosType coeff_flux;
  
  std::vector<SourceShapelets::PepperCorn> grains;
  
  static size_t count;
};

/// A uniform surface brightness circular source.
class SourceUniform : public Source{
public:
  
  SourceUniform(Point_2d position   /// postion on the sky in radians
                ,PosType z          /// redshift of source
                ,PosType radius_in_radians  /// radius of source in radians
                );
	~SourceUniform();

	PosType SurfaceBrightness(const PosType *y) const;
	void assignParams(InputParams& params);
	void printSource();
	PosType getTotalFlux() const {return PI*source_r*source_r;}
};


/// A uniform surface brightness circular source.
class SourcePoint : public SourceColored{
public:
  SourcePoint(
              Point_2d position           /// postion on the sky in radians
              ,PosType z                  /// redshift of source
              ,PosType magnitude         /// unlensed brightness
              ,PosType radius_in_radians  /// radius of source in radians
              ,double SBlimit      /// minimum surface brightness limit
              ,Band band  /// radius of source in radians
):SourceColored(magnitude,radius_in_radians,position,z,SBlimit,band)
{
    sed_type = 1;
  };

  ~SourcePoint(){};
  
  virtual PosType SurfaceBrightness(const PosType *y)const {return 0.0;};
  void printSource(){};
  void assignParams(InputParams& params){}; // do nothing
  void rotate(PosType t){};                 // do nothing
};

/*** \brief A source with a Gaussian surface brightness profile
 
 This is normalized so that it is equal to 1 at the center
 */
class SourceGaussian : public Source{
public:
	//SourceGaussian(InputParams& params);
  SourceGaussian(
                 Point_2d position  /// postion of source (radians)
                 ,double r_size  /// angular scale size (radians)
                 ,double z
                 ,double SBlimit
                 ,double zero_point
                 ):    /// redshift
  Source(0,position,z,SBlimit,zero_point),source_gauss_r2(r_size*r_size)
  {
    zsource = z;
    source_r = 5*sqrt(source_gauss_r2);
    //setSBlimit_magarcsec(100.);
  }

	~SourceGaussian();
	
	/// internal scale parameter
	PosType source_gauss_r2;
	
	PosType SurfaceBrightness(const PosType *y) const;
	void assignParams(InputParams& params);
	void printSource();
	PosType getTotalFlux() const {return 2*PI*source_gauss_r2;/*std::cout << "No total flux in SourceGaussian yet" << std::endl; exit(1);*/}
};

/// Base class for all sources representing the Broad Line Region (BLR) of a AGN/QSO
class SourceBLR : public Source{
public:
	//SourceBLR(InputParams& params);
	~SourceBLR();
	
	void printSource();
	PosType getTotalFlux() const {std::cout << "No total flux in SourceBLR yet" << std::endl; exit(1);}
	
	virtual inline PosType getRadius() const {return source_r_out;}
	
	/// lag time
	PosType source_tau;
	/// frequency
	PosType source_nu;
	float source_nuo;
	/// inner radius of BLR
	float source_r_in;
	/// outer radius of BLR
	float source_r_out;
	///inclination of BLR in radians, face on is
	float source_inclination;
	float source_opening_angle;
	float source_gamma;
	float source_BHmass;
	/// fraction of Keplerian velocity in random motions
	float source_fK;
	/// set to true to integrate over frequency
	bool source_monocrome;
  
  inline PosType getDlDs() const{return DlDs;}

private:
  PosType DlDs;

	void assignParams(InputParams& params);
};

/// A source representing a BLR with a Keplarian disk
class SourceBLRDisk : public SourceBLR{
public:
	PosType SurfaceBrightness(const PosType *y) const;
	PosType getTotalFlux() const {std::cout << "No total flux in SourceBLRDisk yet" << std::endl; exit(1);}
	
	//SourceBLRDisk(InputParams&);
	~SourceBLRDisk();
};

/// A source representing a BLR with a spherical symmetry and circular orbits
class SourceBLRSph1 : public SourceBLR{
public:
	PosType SurfaceBrightness(const PosType *y) const;
	PosType getTotalFlux() const {std::cout << "No total flux in SourceBLRSph1 yet" << std::endl; exit(1);}
	
	//SourceBLRSph1(InputParams&);
	~SourceBLRSph1();
};

/// A source representing a BLR with a spherical symmetry and random velocity dispersion
class SourceBLRSph2 : public SourceBLR{
public:
	PosType SurfaceBrightness(const PosType *y) const;
	PosType getTotalFlux() const {std::cout << "No total flux in SourceBLRSph2 yet" << std::endl; exit(1);}

	//SourceBLRSph2(InputParams&);
	~SourceBLRSph2();
};

/// pointer to surface brightness function
//PosType (Source::*SurfaceBrightness)(PosType *y);


/**  \brief Class to handle redshift-dependent quasar luminosity functions.
 *
 *   At the moment, only i band is available
 *   QLF from Ross et al. 2013
 *   k-correction and mean QSO colors from Richards et al. 2006 (the k-corr is very close to the one used by Ross, small differencies only for z > 3)
 */
class QuasarLF{
	public:
    //QuasarLF(PosType red, PosType mag_limit, InputParams &params);
    ~QuasarLF();
    // returns the integral of the luminosity function at redshift red
    PosType getNorm() {return pow(10,log_phi)*norm;}; // in Mpc^(-3)
    PosType getRandomMag(Utilities::RandomNumbers_NR &rand);
    PosType getRandomFlux(Band band,Utilities::RandomNumbers_NR &rand);
    PosType getColor(Band band);
    PosType getFluxRatio(Band band);

	private:
    PosType kcorr;
    //PosType red;
    //PosType mag_limit;
    //PosType mstar;
    PosType log_phi;
    //PosType alpha;
    //PosType beta;
    PosType norm;
    int arr_nbin;
    PosType* mag_arr;
    PosType* lf_arr;
    PosType dl;
    PosType colors[4];
    PosType ave_colors[4];
    PosType color_dev[4];

    std::string kcorr_file, colors_file;

	void assignParams(InputParams& params);
    
    //typedef PosType (QuasarLF::*pt2MemFunc)(PosType) const;
    //PosType nintegrateQLF(pt2MemFunc func, PosType a,PosType b,PosType tols) const;
    //PosType trapzQLFlocal(pt2MemFunc func, PosType a, PosType b, int n, PosType *s2) const;
    //PosType lf_kernel (PosType mag) const;
    
    struct LF_kernel
    {
        LF_kernel(PosType alpha,PosType beta,PosType mstar)
        : alpha(alpha),beta(beta),mstar(mstar){};
        
        PosType alpha;
        PosType beta;
        PosType mstar;
        
        double operator () (double mag) { 
            return 1.0/( pow(10,0.4*(alpha+1)*(mag-mstar)) + pow(10,0.4*(beta+1)*(mag-mstar)) );
        }
    };
        
};

/// Functor to turn sources into binary functions
struct SourceFunctor
{
	SourceFunctor(Source& source) : source(source) {}
	
	double operator()(double x, double y)
	{
		// TODO: make const double[2] as soon as possible
		double z[2] = {x, y};
		return source.SurfaceBrightness(z);
	}
	
	Source& source;
};


/** \brief Creates a SourcePixelled from a PixelMap image.
 *  The idea is to use stamps of observed galaxies as input sources for simuations.
 *  Surface brightness of the source is conserved, taking into account the input pixel size.
 *  Factor allows for rescaling of the flux, in case one wants to simulate a different observation.
 */
template <typename T>
SourcePixelled::SourcePixelled(
                               const PixelMap<T>& gal_map  /// Input image and information
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



#endif /* SOURCE_H_ */
