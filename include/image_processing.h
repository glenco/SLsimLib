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
#include "source.h"

// forward declaration
struct Grid;
struct GridMap;
class Source;


/**
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
  PixelMap(PixelMap&& other);
	PixelMap(const PixelMap& pmap, const double* center, std::size_t Npixels);
	PixelMap(const double* center, std::size_t Npixels, double resolution);
	PixelMap(const double* center, std::size_t Nx, std::size_t Ny, double resolution);
	PixelMap(std::string fitsfilename,double resolution = -1);
  ~PixelMap(){
    map.resize(0);
  };
	
  PixelMap& operator=(const PixelMap &other);
  PixelMap& operator=(PixelMap &&other);
  
	inline bool valid() const { return map.size(); };
	inline std::size_t size() const { return map.size(); };
	
	inline std::size_t getNx() const { return Nx; }
	inline std::size_t getNy() const { return Ny; }
	inline double getRangeX() const { return rangeX; }
	inline double getRangeY() const { return rangeY; }
  //inline double* getCenter() const{ return &center[0]; }
  void const getCenter(Point_2d &c) const{ c[0]=center[0]; c[1]=center[1];}
  Point_2d getCenter() const{
    Point_2d c(center[0],center[1]);
    return c;
  }
	inline double getResolution() const { return resolution; }
	
  // Zero the whole map
  void Clean(){map *= 0;}

  void AddImages(ImageInfo *imageinfo,int Nimages,float rescale = 1.);
  void AddImages(std::vector<ImageInfo> &imageinfo,int Nimages,float rescale = 1.);
  /// Add an image from a the surface brightnesses of a Grid to the PixelMap
  void AddGridBrightness(Grid &grid);
  /// Add an image from a the surface brightnesses of a GridMap to the PixelMap
  void AddGridMapBrightness(const GridMap &grid);
  void AddUniformImages(ImageInfo *imageinfo,int Nimages,double value);
  PosType AddSource(Source &source);
  /// Add a source to the pixel map by oversamples the source so that oversample^2 points within each pixel are averaged
  PosType AddSource(Source &source,int oversample);
  
  /** \brief copy a PixelMap into this one.
   
   The size, resolution and center of the pixel maps are not changed and do
   not need to match.  The input pixel map is added while conserving the area integral
   of the map within the area of overlaping pixels.
   **/
  void copy_in(const PixelMap& pmap);
  
  /** \brief Replace overlaping pixel values with those of the input map.
   
   No attempt is made to interpolate or average the pixels of pmap.  No integration
   is done.  If the resolution of input map is higher than the reslution of the map then
   the pixels values will just be that of the last pixel visited while going through them.
   
   This can be used to sew tiles together into a larger map.
   **/
  void paste(const PixelMap& pmap);
  
  /** \brief copy a PixelMap that must be the same without creating a new one..
   
   This avoids calling a any constructor or destructor.
   **/
  void duplicate(const PixelMap& pmap);

  /// Adds source to map.  This version breaks pixels up into blocks and does them in seporate threads.
  template <typename T>
  PosType AddSource_parallel(T &source,int oversample){
    Point_2d s_center;
    source.getTheta(s_center);
    
    if( s_center[0] + source.getRadius() < map_boundary_p1[0] ) return 0.0;
    if( s_center[0] - source.getRadius() > map_boundary_p2[0] ) return 0.0;
    if( s_center[1] + source.getRadius() < map_boundary_p1[1] ) return 0.0;
    if( s_center[1] - source.getRadius() > map_boundary_p2[1] ) return 0.0;
    
    int nthreads = Utilities::GetNThreads();
    PosType totals[nthreads];
    std::vector<std::thread> thr;
    
    size_t block = map.size()/nthreads;
    for(int i = 0; i < nthreads ;++i){
      thr.push_back(std::thread(&PixelMap::addsource_<T>,this
                                   ,i*block,std::min((i+1)*block-1,map.size()-1),
                                oversample,std::ref(source),std::ref(totals[i])));
    }
 
    for(int ii=0;ii < nthreads;++ii) thr[ii].join();
    
    PosType total =0;
    for(int ii=0;ii < nthreads;++ii) total += totals[ii];
    
    return total;
  }

  void AddCurve(ImageInfo *curve,double value);
  void AddCurve(Kist<Point> *imagekist,PosType value);
  void AddCurve(std::vector<Point_2d> &curve,double value);
  
	void drawline(double x1[],double x2[],double value);
  void drawcircle(PosType r_center[],PosType radius,PosType value);
  void drawdisk(PosType r_center[],PosType radius,PosType value,int Nstrip);
	void AddGrid(const Grid &grid,double value = 1.0);
  void AddGrid(const Grid &grid,LensingVariable val);

	void Renormalize(double factor);
	void AddValue(std::size_t i, double value);
	void AssignValue(std::size_t i, double value);
	void printASCII() const;
	void printASCIItoFile(std::string filename) const;
	void printFITS(std::string filename, bool verbose = false);
  void printFITS(std::string filename,std::vector<std::tuple<std::string,double,std::string>> &extra_header_info, bool verbose);

	void smooth(double sigma);

	inline double getValue(std::size_t i) const { return map[i]; }
  inline double & operator[](std::size_t i) { return map[i]; };
  const double & operator[](std::size_t i) const { return map[i]; };
  inline double operator()(std::size_t i) const { return map[i]; };
  inline double operator()(std::size_t i,std::size_t j) const { return map[i + Nx*j]; };
	
  
	PixelMap& operator+=(const PixelMap& rhs);
	friend PixelMap operator+(const PixelMap&, const PixelMap&);

	PixelMap& operator-=(const PixelMap& rhs);
	friend PixelMap operator-(const PixelMap&, const PixelMap&);
    
	PixelMap& operator*=(const PixelMap& rhs);
	friend PixelMap operator*(const PixelMap&, const PixelMap&);

	PixelMap& operator*=(PosType b);
 	friend PixelMap operator*(const PixelMap&, PosType b);
  /// eliment wise multiplictioan
  PixelMap operator*=(const PixelMap &m) const;

	std::valarray<double>& data() { return map; }
	
  /// Check whether two PixelMaps agree in their physical dimensions.
	bool agrees(const PixelMap& other) const;
	
	friend void swap(PixelMap&, PixelMap&);
  static void swap(PixelMap&, PixelMap&);
  
  /// return average pixel value
  PosType ave() const;
  /// return sum of all pixel values
  PosType sum() const;
  /// Total number of pixels
  size_t size(){return map.size();}
	
  void FindArc(PosType &radius,PosType *xc,PosType *arc_center,PosType &arclength,PosType &width
                         ,PosType threshold);
  
  /// get the index for a position, returns -1 if out of map, this version returns the 2D grid coordinates
  long find_index(PosType const x[],long &ix,long &iy) const;
  /// get the index for a position, returns -1 if out of map
  long find_index(PosType const x[]) const;
  
  /// get the index for a position, returns -1 if out of map, this version returns the 2D grid coordinates
  long find_index(PosType const x,PosType const y,long &ix,long &iy) const;
  /// get the index for a position, returns -1 if out of map
  long find_index(PosType const x,PosType const y) const;
  
  /// get the index for a position, returns -1 if out of map
  void find_position(PosType x[],std::size_t const index) const;
  /// get the index for a position, returns -1 if out of map
  void find_position(PosType x[],std::size_t const ix,std::size_t const iy) const;
  
  /// interpolate to point x[]
  PosType linear_interpolate(PosType x[]);
  
  /// draw a grid on the image that divides the each demension into N cells
  void drawgrid(int N,PosType value);
  void drawPoints(std::vector<Point *> points,PosType size,PosType value);
  
  void drawPoints(std::vector<Point> points,PosType size,PosType value);
  
  void drawCurve(std::vector<Point *> points,PosType value){
    for(int i=0;i<points.size()-1;++i) drawline(points[i]->x,points[i+1]->x,value);
  }
  void drawCurve(std::vector<Point> points,PosType value){
    for(int i=0;i<points.size()-1;++i) drawline(points[i].x,points[i+1].x,value);
  }
  void drawPoints(std::vector<Point_2d> points,PosType size,PosType value);
  void drawCurve(std::vector<Point_2d> points,PosType value){
    for(int i=0;i<points.size()-1;++i) drawline(points[i].x,points[i+1].x,value);
  }
  /// Draw a rectangle
  void drawSquare(PosType p1[],PosType p2[],PosType value);
  void drawBox(PosType p1[],PosType p2[],PosType value,int Nstrip);
  
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

  /** \brief For a list of pixel indexes this will count and separated islands that are not connected.
   
   On return, 'pixel_index' is ordered into groups and the 'heads' list points to the first elemant 
   in each group plus the end of the list so that heads[i] to heads[i+1] is a group for 0 <= i <= ngroups.
   The number of groups is returned which is also heads.size() - 1
   */
  //int count_islands(std::list<size_t> &pixel_index,std::vector<std::list<size_t>::iterator> &heads) const;
  int count_islands(std::vector<size_t> &pixel_index) const;

  /// get a list of pixels above value
  size_t threshold(std::list<size_t> &pixel_index,PosType value){
    for(size_t i=0;i<map.size();++i) if(value < map[i]) pixel_index.push_back(i);
    return pixel_index.size();
  }
  
  /// reflects the image about the horizontal mid-line
  void flipY(){
    for(size_t i=0;i<Nx;++i){
      for(size_t j=0;j<(Ny-1)/2;++j){
        std::swap( map[i + Nx*j],map[i + Nx*(Ny-j-1)] );
      }
    }
  }
  /// reflects the image about the vertical mid-line
  void flipX(){
    for(size_t i=0;i<(Nx-1)/2;++i){
      for(size_t j=0;j<Ny;++j){
        std::swap( map[i + Nx*j],map[Nx-i-1 + Nx*j] );
      }
    }
  }
  /// rotate the image by 180deg or equivalently reflect it through the origin
  void doubleFlip(){
    flipX();
    flipY();
  }
  
  /// recenter the map without changing anything else
  void recenter(PosType c[2]);
  
  /** \brief convolve the image with a kernel.
   
   It is assumed that the size of the kernel is much smaller than the image and
   that the kernal has the same pixel size as the image.
   **/
  void convolve(PixelMap &kernel,long center_x = 0,long center_y = 0);
  
private:

  std::valarray<double> map;
	std::size_t Nx;
	std::size_t Ny;
	double resolution,rangeX,rangeY,center[2];
	double map_boundary_p1[2],map_boundary_p2[2];
  
  void AddGrid_(const PointList &list,LensingVariable val);

  double LeafPixelArea(IndexType i,Branch * branch1);
	void PointsWithinLeaf(Branch * branch1, std::list <unsigned long> &neighborlist);
	bool inMapBox(Branch * branch1) const;
	bool inMapBox(double * branch1) const;
  
  /// determines if pixels touch each other from i dimensional index
  bool pixels_are_neighbors(size_t i,size_t j) const;
  /** recursive function that finds all the pixels in reservoir beyond and including position 'group' that are attached to pixel current.
   On exit reserve is ordered so that pixels that are in the same group are in sequence and 'group' points to the element in 'reservoir' that is one past the group elements
   */
  void _count_islands_(size_t current,std::list<size_t> &reservoir
                       ,std::list<size_t>::iterator &group) const;
  
  //void addsource_(size_t i1,size_t i2,int oversample,Source source,PosType &total);
  //void addsource_(size_t i1,size_t i2,int oversample,Po,
  //                   PosType &total);

  template <typename T>
    void addsource_(size_t i1,size_t i2,int oversample,
                    T &source,
                    PosType &total){
    PosType tmp_res = resolution*1.0/oversample;
    PosType tmp = tmp_res*tmp_res;
    PosType bl = resolution /2 - 0.5*tmp_res;
    PosType y[2],x[2];
    
    total = 0;
    
    for(size_t index = i1 ;index <= i2; ++index){
      find_position(y,index);
      y[0] -= bl;
      y[1] -= bl;
      for(int i = 0 ; i < oversample ; ++i){
        x[0] = y[0] + i*tmp_res;
        for(int j=0; j < oversample;++j){
          x[1] = y[1] + j*tmp_res;
          map[index] += source.SurfaceBrightness(x)*tmp;
          total += source.SurfaceBrightness(x)*tmp;
        }
      }
    }
  }
  
};

typedef enum {Euclid_VIS,Euclid_Y,Euclid_J,Euclid_H,KiDS_u,KiDS_g,KiDS_r,KiDS_i,HST_ACS_I,CFHT_u,CFHT_g,CFHT_r,CFHT_i,CFHT_z} Telescope;

typedef enum {counts_x_sec, flux} unitType;

/** 
 * \brief It creates a realistic image from the output of a ray-tracing simulation.
 *
 * It translates pixel values in observed units (counts/sec), applies PSF and noise.
 * Input must be in ergs/(s*cm^2*Hz*hplanck).
 */
class Observation
{
public:
	Observation(Telescope tel_name,size_t Npix_x,size_t Npix_y);
	Observation(float diameter, float transmission, float exp_time, int exp_num, float back_mag, float ron
              ,size_t Npix_x,size_t Npix_y,float seeing = 0.);
	Observation(float diameter, float transmission, float exp_time, int exp_num, float back_mag, float ron ,std::string psf_file,size_t Npix_x,size_t Npix_y, float oversample = 1.);

  
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
  float getBackgroundNoise(float resolution, unitType unit = counts_x_sec);
	std::valarray<double> getPSF(){return map_psf;}
  void setPSF(std::string psf_file, float os = 1.);
  void setNoiseCorrelation(std::string nc_file);
	void Convert(PixelMap &map, bool psf, bool noise,long *seed, unitType unit = counts_x_sec);
  /// returns factor by which code image units need to be multiplied by to get flux units
  double flux_convertion_factor(){ return pow(10,-0.4*mag_zeropoint); }

	void Convert_back(PixelMap &map);
  void setExpTime(float time){exp_time = time;}

  size_t getNx(){ return Npix_x;}
  size_t getNy(){ return Npix_y;}

private:
  
  size_t Npix_x,Npix_y;
  std::vector<double> sqrt_noise_power;  // stores sqrt root of power noise spectrum
  size_t side_ncorr; // pixels on a side of input noise correlation function

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

	void AddNoise(PixelMap &pmap,long *seed);
  
	void PhotonToCounts(PixelMap &pmap);
	void ApplyPSF(PixelMap &pmap);
  void CorrelateNoise(PixelMap &pmap);

  //PixelMap noise_correlation;
};

void pixelize(double *map,long Npixels,double range,double *center
		,ImageInfo *imageinfo,int Nimages,bool constant_sb,bool cleanmap
		,bool write_for_skymaker = false, std::string filename="");
void _SplitFluxIntoPixels(TreeHndl ptree,Branch *leaf,double *leaf_sb);
//void smoothmap(double *map_out,double *map_in,long Npixels,double range,double sigma);

namespace Utilities{
    void LoadFitsImages(std::string dir,const std::string& filespec,std::vector<PixelMap> & images,int maxN,double resolution = -1,bool verbose = false);
    void LoadFitsImages(std::string dir,std::vector<std::string> filespecs,std::vector<std::string> file_non_specs                                  ,std::vector<PixelMap> & images,std::vector<std::string> & names,int maxN,double resolution = -1,bool verbose = false);
    void ReadFileNames(std::string dir,const std::string filespec
                       ,std::vector<std::string> & filenames
                       ,const std::string file_non_spec = " "
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
