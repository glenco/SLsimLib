/*
 * image_processing.h
 *
 *  Created on: Feb 28, 2010
 *      Author: R.B. Metcalf
 */

#ifndef IMAGE_PROCESSING_H_
#define IMAGE_PROCESSING_H_
 
#include "Tree.h"
#ifdef ENABLE_FFTW
#include "fftw3.h"
#endif

#include "utilities_slsim.h"

// forward declaration
struct Grid;

/** \ingroup Image
 * \brief Takes image structure and pixelizes the flux into regular pixel grid which then
 * can be exported as a fits file, smoothed, etc. like an image.
 *
 */
class PixelMap
{
public:
	PixelMap(const PixelMap& pmap, double res_ratio);
	PixelMap();
	PixelMap(const PixelMap& other);
	PixelMap(const PixelMap& pmap, const double* center, std::size_t Npixels);
	PixelMap(const double* center, std::size_t Npixels, double resolution);
	PixelMap(const double* center, std::size_t Nx, std::size_t Ny, double resolution);
	PixelMap(std::string fitsfilename,double resolution = -1);
	~PixelMap();
	
	PixelMap& operator=(PixelMap other);
	
	inline bool valid() const { return map.size(); };
	inline std::size_t size() const { return map.size(); };
	
	inline std::size_t getNx() const { return Nx; }
	inline std::size_t getNy() const { return Ny; }
	inline double getRangeX() const { return rangeX; }
	inline double getRangeY() const { return rangeY; }
	inline const double* getCenter() const{ return center; }
	inline double getResolution() const { return resolution; }
	
	void Clean();

  void AddImages(ImageInfo *imageinfo,int Nimages,float rescale = 1.);
  void AddImages(std::vector<ImageInfo> &imageinfo,int Nimages,float rescale = 1.);
  void AddUniformImages(ImageInfo *imageinfo,int Nimages,double value);
  
  void AddCurve(ImageInfo *curve,double value);
  void AddCurve(Kist<Point> *imagekist,PosType value);
  void AddCurve(std::vector<Point_2d> &curve,double value);
  
	void drawline(double x1[],double x2[],double value);
  void drawcircle(PosType r_center[],PosType radius,PosType value);
	void AddGrid(const Grid &grid,double value = 1.0);
  void AddGrid(const Grid &grid,LensingVariable val);

	void Renormalize(double factor);
	void AddValue(std::size_t i, double value);
	void AssignValue(std::size_t i, double value);
	void printASCII() const;
	void printASCIItoFile(std::string filename) const;
	void printFITS(std::string filename, bool verbose = false) const;
  void printFITS(std::string filename,std::vector<std::tuple<std::string,double,std::string>> &extra_header_info, bool verbose) const;

	void smooth(double sigma);

	inline double getValue(std::size_t i) const { return map[i]; }
	inline double & operator[](std::size_t i) { return map[i]; };
	inline double operator()(std::size_t i) const { return map[i]; };
	
	PixelMap& operator+=(const PixelMap& rhs);
	friend PixelMap operator+(const PixelMap&, const PixelMap&);

	PixelMap& operator-=(const PixelMap& rhs);
	friend PixelMap operator-(const PixelMap&, const PixelMap&);
    
	PixelMap& operator*=(const PixelMap& rhs);
	friend PixelMap operator*(const PixelMap&, const PixelMap&);

	PixelMap& operator*=(PosType b);
	friend PixelMap operator*(const PixelMap&, PosType b);
	
	std::valarray<double>& data() { return map; }
	
	bool agrees(const PixelMap& other) const;
	
	friend void swap(PixelMap&, PixelMap&);
  
  /// return average pixel value
  PosType ave() const;
  /// Total number of pixels
  size_t size(){return map.size();}
	
  void FindArc(PosType &radius,PosType *xc,PosType *arc_center,PosType &arclength,PosType &width
                         ,PosType threshold);
  
  /// get the index for a position, returns -1 if out of map, this version returns the 2D grid coordinates
  long find_index(PosType const x[],long &ix,long &iy);
  /// get the index for a position, returns -1 if out of map
  long find_index(PosType const x[]);
  
  /// get the index for a position, returns -1 if out of map, this version returns the 2D grid coordinates
  long find_index(PosType const x,PosType const y,long &ix,long &iy);
  /// get the index for a position, returns -1 if out of map
  long find_index(PosType const x,PosType const y);
  
  /// get the index for a position, returns -1 if out of map
  void find_position(PosType x[],std::size_t const index);
  /// get the index for a position, returns -1 if out of map
  void find_position(PosType x[],std::size_t const ix,std::size_t const iy);
  
  /// interpolate to point x[]
  PosType linear_interpolate(PosType x[]);
  
  /// draw a grid on the image that divides the each demension into N cells
  void drawgrid(int N,PosType value);
  void drawPoints(std::vector<Point *> points,PosType size,PosType value){
    if(size < resolution*3){
      size_t index;
      for(int i=0;i<points.size();++i){
        if(inMapBox(points[i]->x)){
          //index = Utilities::IndexFromPosition(x1,Nx,range,center);
          index = find_index(points[i]->x);
          map[index] = value;
        }
      }
    }else
      for(int i=0;i<points.size();++i) drawcircle(points[i]->x,0.01*rangeX,value);
    
  }
  void drawCurve(std::vector<Point *> points,PosType value){
    for(int i=0;i<points.size()-1;++i) drawline(points[i]->x,points[i+1]->x,value);
  }
  /// Draw a rectangle
  void drawBox(PosType p1[],PosType p2[],PosType value);
  
#ifdef ENABLE_FFTW

  /// Find the power spectrum of the map
  void PowerSpectrum(std::vector<PosType> &power_spectrum   /// output power spectrum
                     ,std::vector<PosType> &lvec            /// output l values of bands
                     ,bool overwrite = true                 /// if false add power to existing power_spectrum (used for averaging over many fields
                     ){
    
    if(power_spectrum.size() != lvec.size()) throw std::invalid_argument("these must be the same size");
    
    if(overwrite) Utilities::powerspectrum2d(map,map,Nx,Ny,rangeX,rangeY, lvec, power_spectrum);
    else{
      std::vector<PosType> tmp_power(power_spectrum.size());
      Utilities::powerspectrum2d(map,map,Nx,Ny,rangeX,rangeY, lvec, tmp_power);
      for(size_t ii=0;ii<power_spectrum.size();++ii) power_spectrum[ii] += tmp_power[ii];
    }
  }

  void AdaptiveSmooth(PosType value){
    
    std::valarray<double> tmp = Utilities::AdaptiveSmooth(data(),Nx,Ny,value);
    map = tmp;
  }

#endif

private:
	std::valarray<double> map;
  void AddGrid_(const PointList &list,LensingVariable val);

	std::size_t Nx;
	std::size_t Ny;
	double resolution,rangeX,rangeY,center[2];
	double map_boundary_p1[2],map_boundary_p2[2];

	double LeafPixelArea(IndexType i,Branch * branch1);
	void PointsWithinLeaf(Branch * branch1, std::list <unsigned long> &neighborlist);
	bool inMapBox(Branch * branch1) const;
	bool inMapBox(double * branch1) const;
  
};

typedef enum {Euclid_VIS,Euclid_Y,Euclid_J,Euclid_H,KiDS_u,KiDS_g,KiDS_r,KiDS_i,HST_ACS_I,CFHT_u,CFHT_g,CFHT_r,CFHT_i,CFHT_z} Telescope;

typedef enum {counts_x_sec, flux} unitType;

/** \ingroup Image
 * \brief It creates a realistic image from the output of a ray-tracing simulation.
 *
 * It translates pixel values in observed units (counts/sec), applies PSF and noise.
 * Input must be in photons/(cm^2*Hz).
 */
class Observation
{
public:
	Observation(Telescope tel_name);
	Observation(float diameter, float transmission, float exp_time, int exp_num, float back_mag, float ron, float seeing = 0.);
	Observation(float diameter, float transmission, float exp_time, int exp_num, float back_mag, float ron, std::string psf_file, float oversample = 1.);
	float getExpTime(){return exp_time;}
	int getExpNum(){return exp_num;}
	float getBackMag(){return back_mag;}
	float getDiameter(){return diameter;}
	float getTransmission(){return transmission;}
    /// read-out noise in electrons/pixel
	float getRon(){return ron;}
  /// seeing in arcsecs
	float getSeeing(){return seeing;}
	float getZeropoint(){return mag_zeropoint;}
    /// pixel size in radians
  float getPixelSize(){return pix_size;}
	std::valarray<double> getPSF(){return map_psf;}
  void setPSF(std::string psf_file, float os = 1.);
	PixelMap Convert (PixelMap &map, bool psf, bool noise,long *seed, unitType unit = counts_x_sec);
	PixelMap Convert_back (PixelMap &map);
    void setExpTime(float time){exp_time = time;}

private:
	float diameter;  // diameter of telescope (in cm)
	float transmission;  // total transmission of the instrument
	float mag_zeropoint;  // magnitude of a source that produces one count/sec in the image
	float exp_time;  // total exposure time (in sec)
	int exp_num;  // number of exposures
	float back_mag;  // sky (or background) magnitude in mag/arcsec^2
	float ron;  // read-out noise in electrons/pixel
	float seeing;  // full-width at half maximum of the gaussian smoothing
	std::valarray<double> map_psf;  // array of the point spread function
	float oversample; // psf oversampling factor
	double pix_size; // pixel size (in rad)
	bool telescope; // was the observation created from a default telescope?

	PixelMap AddNoise(PixelMap &pmap,long *seed);
	PixelMap PhotonToCounts(PixelMap &pmap);
	PixelMap ApplyPSF(PixelMap &pmap);

};

void pixelize(double *map,long Npixels,double range,double *center
		,ImageInfo *imageinfo,int Nimages,bool constant_sb,bool cleanmap
		,bool write_for_skymaker = false, std::string filename="");
void _SplitFluxIntoPixels(TreeHndl ptree,Branch *leaf,double *leaf_sb);
//void smoothmap(double *map_out,double *map_in,long Npixels,double range,double sigma);

namespace Utilities{
    void LoadFitsImages(std::string dir,const std::string& filespec,std::vector<PixelMap> & images,int maxN,double resolution = -1,bool verbose = false);
  void LoadFitsImages(std::string dir,std::vector<std::string> filespecs,std::vector<std::string> file_non_specs                                  ,std::vector<PixelMap> & images,std::vector<std::string> & names,int maxN,double resolution = -1,bool verbose = false);
    void ReadFileNames(std::string dir,const std::string filespec,std::vector<std::string> & filenames
                       ,bool verbose = false);
}

/** \brief Warning: Not tested yet. Class for doing adaptive smoothing using multiply resolution grids.
 */
class MultiGridSmoother{
public:
  MultiGridSmoother(double center[],std::size_t Nx,std::size_t Ny,double resolution);
  MultiGridSmoother(double center[],std::size_t Nx,double resolution);
  ~MultiGridSmoother(void){
    maps.clear();
  }
  
  /// resolution of finest grid from which interpolation is done
  PosType getHighestRes(){return maps[0].getResolution();}
  /// resolution of coarsest grid from which interpolation is done
  PosType getLowestRes(){return maps.back().getResolution();}
  
  /// Add particles to the map.  These do not need to be kept in memory after they are added.
  void add_particles(std::vector<PosType> x,std::vector<PosType> y);
  /// Output a map at the resolution of the map smoothed so that no superpixel as less than Nsmooth particles
  void output_map(PixelMap &map,int Nsmooth);
  void smooth(int Nsmooth,PixelMap &map);
  
private:
  void _smooth_(int k,size_t i,size_t j,int Nsmooth,PixelMap &map);
  std::vector<PixelMap> maps;
  std::vector<Utilities::Interpolator<PixelMap>> interpolators;
};
#endif
