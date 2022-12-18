/*
 * source.h
 *
 *  Created on: Feb 6, 2012
 *      Author: mpetkova
 */

#ifndef SOURCE_H_
#define SOURCE_H_

#include "standard.h"
#include "InputParams.h"
#include "image_processing.h"
#include "utilities_slsim.h"

double flux_to_mag(double flux);
double mag_to_flux(double m);

class PixelMap;

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
         ,PosType SBlimit=100  // surface brightness limit in mags per arcsecond
         ):source_r(r),source_x(x),zsource(z)
  {
    setSBlimit_magarcsec(SBlimit);
  };
  
  virtual ~Source();

  Source(const Source &s):
    source_r(s.source_r),
    source_x(s.source_x),
    zsource(s.zsource),
    sb_limit(s.sb_limit){}
  
  Source & operator=(const Source &s){
    if(this == &s) return *this;
    source_r = s.source_r;
    source_x = s.source_x;
    zsource = s.zsource;
    sb_limit = s.sb_limit;
    
    return *this;
  }

	// in lens.cpp
	/** Surface brightness of source in grid coordinates not source centered coordinates.
   
   * The units shuld be ergs / s / Hz / cm^2 / hplanck
   */
  virtual PosType SurfaceBrightness(PosType *y) = 0;
  virtual PosType getTotalFlux() const = 0;
  virtual void printSource() = 0;

  
	/// Gets sb_limit in erg/cm^2/sec/rad^2/Hz
	PosType getSBlimit(){return sb_limit;}
	/// Gets sb_limit in mag/arcsec^2
  PosType getSBlimit_magarcsec(){return flux_to_mag(sb_limit/pow((180*60*60/PI),2));}

	// accessor functions that will sometimes be over ridden in class derivatives
	/// Redshift of source
	PosType getZ() const {return zsource;}
	void setZ(PosType my_z){zsource = my_z;}
	/// Radius of source in radians
	inline PosType getRadius() const {return source_r;}
  /// Reset the radius of the source in radians
	void setRadius(PosType my_radius){source_r = my_radius;}
  /// position of source in radians
  inline Point_2d getTheta(){return source_x;}
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
	/// Sets sb_limit in mag/arcsec^2
  void setSBlimit_magarcsec(float limit) {sb_limit = mag_to_flux(limit)*pow(180*60*60/PI,2);}

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
  
};

typedef Source *SourceHndl;

class SourceColored : public Source{
public:
  SourceColored(PosType magnitude,PosType r,Point_2d x,PosType z
                ,PosType SBlimit=40):Source(r,x,z,SBlimit)
  {
    current_band = Band::NoBand;
    setMag(magnitude);
    id = -1;
  }
  
  SourceColored(const SourceColored &s):Source(s){
    mag = s.mag;
    mag_map = s.mag_map;
    current_band = s.current_band;
    sed_type = s.sed_type;
    id = s.id;
    flux_total = s.flux_total;
  }

  SourceColored & operator= (const SourceColored &s){
    if(this == &s) return *this;
    
    Source::operator=(s);
    mag = s.mag;
    mag_map = s.mag_map;
    current_band = s.current_band;
    sed_type = s.sed_type;
    id = s.id;
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
    id = s.id;
    flux_total = s.flux_total;
    
    return *this;
  }
  
  virtual ~SourceColored(){};
  
  PosType getMag() const { assert(current_band != Band::NoBand); return mag;}
  PosType getMag(Band band) const {return (mag_map.size() > 0) ? mag_map.at(band) : mag;}
  Band getBand() const{return current_band;}
  long getID() const {return id;}
  float getSEDtype() const {return sed_type;}
  
  void setSEDtype(float sed) {sed_type = sed;}
  void setActiveBand(Band band);
  void setBand(Band band,float m){mag_map[band]=m;};
  inline PosType getTotalFlux() const {return flux_total;}
  
  void setMag(float magnitude){
    if(current_band == Band::NoBand){
      mag=magnitude;
      flux_total = mag_to_flux(mag);
    }else{
      std::cerr << "Cannot change mag of source that has multiple bands." << std::endl;
      throw std::runtime_error("source has ");
    }
  }
  
  // rotate on the sky
  virtual void rotate(PosType theta) = 0;

  /// shift all magnitudes my delta_mag and update flux_total, keeps color
  void shiftmag(float delta_mag){
    mag += delta_mag;
    flux_total = mag_to_flux(mag);
    for(auto &b : mag_map) b.second += delta_mag;
  }
protected:
  std::map<Band,PosType> mag_map;
  Band current_band;
  float sed_type;
  PosType mag;
  int id;
  PosType flux_total;
};



/** \brief Class for sources described by an array of pixels
 *
 *  The sources are created as a square array of pixels of dimension Npixels x Npixels and pixel size equal to resolution.
 *
 */
class SourcePixelled: public Source{
public:
	SourcePixelled(PosType my_z, PosType* center, int Npixels, PosType resolution, PosType* arr_val);
	SourcePixelled(const PixelMap& gal_map, PosType z, PosType factor = 1.);
	//SourcePixelled(InputParams& params);
  
	~SourcePixelled();
	PosType SurfaceBrightness(PosType *y);
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
	SourceShapelets(PosType my_z, PosType my_mag, PosType my_scale, std::valarray<PosType> my_coeff, PosType* my_center = 0, PosType my_ang = 0.);
	SourceShapelets(PosType my_z, PosType my_mag, std::string shap_file, PosType *my_center = 0, PosType my_ang = 0.);
  SourceShapelets(std::string shap_file, PosType my_ang = 0.);//, PosType* my_center = 0);
  
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

	PosType SurfaceBrightness(PosType *y);
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
  void Hermite(std::vector<PosType> &hg,int N, PosType x);

	void NormalizeFlux();
	std::valarray<PosType> coeff;
	int n1,n2;
 
	//PosType ang;
  PosType cos_sin[2];  // [0] is cos(ang), [1] is sin(ang)
  PosType coeff_flux;
  
  std::vector<PepperCorn> grains;
  
  static size_t count;
};

/// A uniform surface brightness circular source.
class SourceUniform : public Source{
public:
  //SourceUniform(InputParams& params);
  SourceUniform(Point_2d position   /// postion on the sky in radians
                ,PosType z          /// redshift of source
                ,PosType radius_in_radians  /// radius of source in radians
                );
	~SourceUniform();

	PosType SurfaceBrightness(PosType *y);
	void assignParams(InputParams& params);
	void printSource();
	PosType getTotalFlux() const {return PI*source_r*source_r;}
};


/// A uniform surface brightness circular source.
class SourcePoint : public SourceColored{
public:
//  SourcePoint(Point_2d position         /// postion on the sky in radians
//                ,PosType z              /// redshift of source
//                ,PosType brightness     /// unlensed brightness
//                ,PosType radius_in_radians  /// radius of source in radians
//                );
  SourcePoint(
                           Point_2d position           /// postion on the sky in radians
                           ,PosType z                  /// redshift of source
                           ,PosType magnitude         /// unlensed brightness
                           ,PosType radius_in_radians  /// radius of source in radians
                  ):SourceColored(magnitude,radius_in_radians,position,z)
  {
    sed_type = 1;
  };

  ~SourcePoint(){};
  
  PosType SurfaceBrightness(PosType *y){return 0.0;};
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
                 ,double z):    /// redshift
  Source(0,position,z),source_gauss_r2(r_size*r_size)
  {
    zsource = z;
    source_r = 5*sqrt(source_gauss_r2);
    setSBlimit_magarcsec(100.);
  }

	~SourceGaussian();
	
	/// internal scale parameter
	PosType source_gauss_r2;
	
	PosType SurfaceBrightness(PosType *y);
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
  
  inline PosType getDlDs(){return DlDs;}

private:
  PosType DlDs;

	void assignParams(InputParams& params);
};

/// A source representing a BLR with a Keplarian disk
class SourceBLRDisk : public SourceBLR{
public:
	PosType SurfaceBrightness(PosType *y);
	PosType getTotalFlux() const {std::cout << "No total flux in SourceBLRDisk yet" << std::endl; exit(1);}
	
	//SourceBLRDisk(InputParams&);
	~SourceBLRDisk();
};

/// A source representing a BLR with a spherical symmetry and circular orbits
class SourceBLRSph1 : public SourceBLR{
public:
	PosType SurfaceBrightness(PosType *y);
	PosType getTotalFlux() const {std::cout << "No total flux in SourceBLRSph1 yet" << std::endl; exit(1);}
	
	//SourceBLRSph1(InputParams&);
	~SourceBLRSph1();
};

/// A source representing a BLR with a spherical symmetry and random velocity dispersion
class SourceBLRSph2 : public SourceBLR{
public:
	PosType SurfaceBrightness(PosType *y);
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
    PosType red;
    PosType mag_limit;
    PosType mstar;
    PosType log_phi;
    PosType alpha;
    PosType beta;
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





#endif /* SOURCE_H_ */
