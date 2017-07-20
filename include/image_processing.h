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

/** \ingroup Image
 * \brief Takes image structure and pixelizes the flux into regular pixel grid which then
 * can be exported as a fits file, smoothed, etc. like an image.
 *
 */
class Source;

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
  void const getCenter(Point_2d &c) const{ c[0]=center[0]; c[1]=center[1];}
  inline double getResolution() const { return resolution; }
  Point_2d getLowerLeft() const { return Point_2d(map_boundary_p1[0],map_boundary_p1[1]); }
  Point_2d getUpperRight() const { return Point_2d(map_boundary_p2[0],map_boundary_p2[1]); }
  
  inline double getValue(std::size_t i) const { return map[i]; }
  inline double & operator[](std::size_t i) { return map[i]; };
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
  
  std::valarray<double>& data() { return map; }
	
	void Clean();

  void AddImages(ImageInfo *imageinfo,int Nimages,float rescale = 1.);
  void AddImages(std::vector<ImageInfo> &imageinfo,int Nimages,float rescale = 1.);
  void AddGridBrightness(Grid &grid);
  void AddUniformImages(ImageInfo *imageinfo,int Nimages,double value);
  PosType AddSource(Source &source);
  /// Add a source to the pixel map by oversamples the source so that oversample^2 points within each pixel are averaged
  PosType AddSource(Source &source,int oversample);

  /// Adds source to map.  This version breaks pixels up into blocks and does them in seporate threads.
  template <typename T>
  PosType AddSource_parallel(T &source,int oversample){
    Point_2d s_center;
    source.getX(s_center);
    
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
	void printFITS(std::string filename, bool verbose = false) const;
  void printFITS(std::string filename,std::vector<std::tuple<std::string,double,std::string>> &extra_header_info, bool verbose) const;

	void smooth(double sigma);
	
	bool agrees(const PixelMap& other) const;
	
	friend void swap(PixelMap&, PixelMap&);
  void swap(PixelMap&, PixelMap&);
  
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
  void drawPoints(std::vector<Point *> points,PosType size,PosType value);
  void drawPoints(std::vector<Point> points,PosType size,PosType value);
  void drawCurve(std::vector<Point *> points,PosType value);
  void drawCurve(std::vector<Point> points,PosType value);
  void drawPoints(std::vector<Point_2d> points,PosType size,PosType value);
  void drawCurve(std::vector<Point_2d> points,PosType value);
  /// Draw a rectangle
  void drawSquare(PosType p1[],PosType p2[],PosType value);
  void drawBox(PosType p1[],PosType p2[],PosType value,int Nstrip);
  
#ifdef ENABLE_FFTW

  /// Find the power spectrum of the map
  void PowerSpectrum(std::vector<PosType> &power_spectrum   /// output power spectrum
                     ,std::vector<PosType> &lvec            /// output l values of bands
                     ,int zeropaddingfactor = 1
                     ,bool overwrite = true                 /// if false add power to existing power_spectrum (used for averaging over many fields
                     ,bool antialias =false
                     );
  
  void AdaptiveSmooth(PosType value){
    std::valarray<double> tmp = Utilities::AdaptiveSmooth(data(),Nx,Ny,value);
    map = tmp;
  }
  
  /// Calculate the deplection and shear from a convergence map using the single plane approximation
  void Shear(PixelMap &alpha1,PixelMap &alpha2,PixelMap &gamma1,PixelMap &gamma2,double zerosize);

#endif

  /** \brief For a list of pixel indexes this will count and separated islands that are not connected.
   
   On return, 'pixel_index' is ordered into groups and the 'heads' list points to the first elemant 
   in each group plus the end of the list so that heads[i] to heads[i+1] is a group for 0 <= i <= ngroups.
   The number of groups is returned which is also heads.size() - 1
   */
  int count_islands(std::list<size_t> &pixel_index,std::vector<std::list<size_t>::iterator> &heads) const;
  /// get a list of pixels above value
  size_t threshold(std::list<size_t> &pixel_index,PosType value){
    for(size_t i=0;i<map.size();++i) if(value < map[i]) pixel_index.push_back(i);
    return pixel_index.size();
  }
  /// check if two maps are the same in dimension and center
  bool compare(const PixelMap &m){
    return (Nx==m.Nx)*(Ny==m.Ny)*(rangeX==m.rangeX)*(rangeY==m.rangeY)
    *(center[0] == m.center[0])*(center[1] == m.center[1]);
  };

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

/** \ingroup Image
 * \brief It creates a realistic image from the output of a ray-tracing simulation.
 *
 * It translates pixel values in observed units (counts/sec), applies PSF and noise.
 * Input must be in ergs/(s*cm^2*Hz*hplanck).
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
  float getBackgroundNoise(float resolution, unitType unit = counts_x_sec);
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
  /// check is file exists
  bool file_exists (const std::string& name);
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

/**
 * \brief Computing deflection angles and shear by FFT from a rectangular map of the convergence / surface density
 */

template <typename V>  // V type must be like a vector - have resize() and [] operator
void ShearByFFT(
                       const V &convergence  /// vector or valarray of convergence or surface density
                       ,V &alpha1            /// output deflection
                       ,V &alpha2            /// output deflection
                       ,V &gamma1            /// output shear
                       ,V &gamma2            /// output shear
                       ,size_t nx            /// number of pixels in x direction
                       ,size_t ny            /// number of pixels in y direction
                       ,double length        /// length of x in units that will be used for derivatives
                       ,double zeropadding = 1.0  /// factor my which the
                       ){
#ifdef ENABLE_FFTW
  // initialize the quantities
  //int npix_filter = 0;   // filter the map if you want on a given number of pixels: CHECK IT NOT IMPLEMENTED YET
  
  // size of the new map in x and y directions, factor by which each size is increased
  int Nnx=int(zeropadding*nx);
  int Nny=int(zeropadding*ny);
  double Nboxlx = length*zeropadding;
  double Nboxly = length*zeropadding*ny/nx;
  
  std:: valarray<float> Nmap;
  try{
    Nmap.resize( Nnx*Nny );
  }catch(std::exception &e){
    std::cerr << "exception thrown in LensHaloMassMap::PreProcessFFTWMap(): " << e.what() << std::endl;
  }
  // assume locate in a rectangular map and build up the new one
  for( int j=0; j<Nny; j++ ){
    for( int i=0; i<Nnx; i++ ){
      Nmap[i+Nnx*j]=0;
      if(i>=int(Nnx/2-nx/2) && i<int(Nnx/2+nx/2) && j>=int(Nny/2-ny/2) && j<int(Nny/2+ny/2)){
        int ii = i-int(Nnx/2-nx/2);
        int jj = j-int(Nny/2-ny/2);
        
        if(ii>=nx || jj>=ny){
          std::cout << " 1 error mapping " << ii << "  " << jj << std::endl;
          exit(1);
        }
        if(ii<0 || jj<0){
          std::cout << " 2 error mapping " << ii << "  " << jj << std::endl;
          exit(1);
        }
        Nmap[i+Nnx*j]=convergence[ii+nx*jj];
      }
    }
  }
  
  double *dNmap=new double[Nnx*Nny];
  double *input=new double[Nnx*Nny];
  fftw_complex *fNmap=new fftw_complex[Nny*(Nnx/2+1)];
  fftw_complex *output=new fftw_complex[Nny*(Nnx/2+1)];
  for(int k=0;k<Nnx*Nny;k++) dNmap[k] = double(Nmap[k]);
  fftw_plan p;
  p=fftw_plan_dft_r2c_2d(Nny,Nnx,input,output,FFTW_ESTIMATE);
  
  for(int i=0;i<Nnx*Nny;i++) input[i] = dNmap[i];
  fftw_execute( p );
  for(int i=0; i<Nny*(Nnx/2+1);i++){
    fNmap[i][0] = output[i][0];
    fNmap[i][1] = output[i][1];
  }
  delete[] input;
  delete[] output;
  delete[] dNmap;
  fftw_destroy_plan(p);
  
  // fourier space
  // std:: cout << " allocating fourier space maps " << std:: endl;
  
  
  fftw_complex *fphi   = new fftw_complex[Nny*(Nnx/2+1)];
  // build modes for each pixel in the fourier space
  for( int i=0; i<Nnx/2+1; i++ ){
    double kx=double(i);
    kx=kx*2.*M_PI/Nboxlx;
    for( int j=0; j<Nny; j++ ){
      double ky=(j<Nny/2)?double(j):double(j-Nny);
      ky=ky*2.*M_PI/Nboxly;
      double k2 = kx*kx+ky*ky;
      //smooth if you want to IMPLEMENT
      // if(npix_filter>0){
      // fNmap[j+(Nnpixels/2+1)*i][0] = fNmap[j+(Nnpixels/2+1)*i][0]*exp(-k2/sigmag/sigmag/2.);
      // fNmap[j+(Nnpixels/2+1)*i][1] = fNmap[j+(Nnpixels/2+1)*i][1]*exp(-k2/sigmag/sigmag/2.);
      // }
      
      // fphi
      fphi[i+(Nnx/2+1)*j][0]= -2.*fNmap[i+(Nnx/2+1)*j][0]/k2;
      fphi[i+(Nnx/2+1)*j][1]= -2.*fNmap[i+(Nnx/2+1)*j][1]/k2;
      // null for k2 = 0 no divergence
      if(k2 == 0){
        fphi[i+(Nnx/2+1)*j][0] = 0.;
        fphi[i+(Nnx/2+1)*j][1] = 0.;
      }
    }
  }
  
  delete[] fNmap;
  
  
  fftw_complex *fft= new fftw_complex[Nny*(Nnx/2+1)];
  double *realsp = new double[Nnx*Nny];
  //fftw_plan pp = fftw_plan_dft_c2r_2d(Nny,Nnx,fft,realsp,FFTW_ESTIMATE);
  fftw_plan pp = fftw_plan_dft_c2r_2d(Nny,Nnx,fft,realsp,FFTW_MEASURE);
  
  // alpha1
  {
    
    // build modes for each pixel in the fourier space
    for( int i=0; i<Nnx/2+1; i++ ){
      double kx=double(i);
      kx=kx*2.*M_PI/Nboxlx;
      for( int j=0; j<Nny; j++ ){
        double ky=(j<Nny/2)?double(j):double(j-Nny);
        ky=ky*2.*M_PI/Nboxly;
        
        fft[i+(Nnx/2+1)*j][0] = -kx*fphi[i+(Nnx/2+1)*j][1];
        fft[i+(Nnx/2+1)*j][1] =  kx*fphi[i+(Nnx/2+1)*j][0];
      }
    }
    
    fftw_execute( pp );
    
    alpha1.resize(nx*ny);
    
    for( int j=Nny/2-ny/2; j<Nny/2+ny/2; j++ ){
      for( int i=Nnx/2-nx/2; i<Nnx/2+nx/2; i++ ){
        int ii = i-int(Nnx/2-nx/2);
        int jj = j-int(Nny/2-ny/2);
        
        alpha1[ii+nx*jj] = -1*float(realsp[i+Nnx*j]/Nnx/Nny);
      }
    }
  }
  
  // alpha2
  {
    
    // build modes for each pixel in the fourier space
    for( int i=0; i<Nnx/2+1; i++ ){
      double kx=double(i);
      kx=kx*2.*M_PI/Nboxlx;
      for( int j=0; j<Nny; j++ ){
        double ky=(j<Nny/2)?double(j):double(j-Nny);
        ky=ky*2.*M_PI/Nboxly;
        
        // alpha
        fft[i+(Nnx/2+1)*j][0] = -ky*fphi[i+(Nnx/2+1)*j][1];
        fft[i+(Nnx/2+1)*j][1] =  ky*fphi[i+(Nnx/2+1)*j][0];
        
      }
    }
    
    fftw_execute( pp );
    
    alpha2.resize(nx*ny);
    
    for( int j=Nny/2-ny/2; j<Nny/2+ny/2; j++ ){
      for( int i=Nnx/2-nx/2; i<Nnx/2+nx/2; i++ ){
        int ii = i-int(Nnx/2-nx/2);
        int jj = j-int(Nny/2-ny/2);
        
        alpha2[ii+nx*jj] = -1*float(realsp[i+Nnx*j]/Nnx/Nny);
      }
    }
  }
  // gamma1
  {
    
    // build modes for each pixel in the fourier space
    for( int i=0; i<Nnx/2+1; i++ ){
      double kx=double(i);
      kx=kx*2.*M_PI/Nboxlx;
      for( int j=0; j<Nny; j++ ){
        double ky=(j<Nny/2)?double(j):double(j-Nny);
        ky=ky*2.*M_PI/Nboxly;
        
        // gamma
        fft[i+(Nnx/2+1)*j][0] = 0.5*(kx*kx-ky*ky)*fphi[i+(Nnx/2+1)*j][0];
        fft[i+(Nnx/2+1)*j][1] = 0.5*(kx*kx-ky*ky)*fphi[i+(Nnx/2+1)*j][1];
        
      }
    }
    
    fftw_execute( pp );
    
    gamma1.resize(nx*ny);
    
    for( int j=Nny/2-ny/2; j<Nny/2+ny/2; j++ ){
      for( int i=Nnx/2-nx/2; i<Nnx/2+nx/2; i++ ){
        int ii = i-int(Nnx/2-nx/2);
        int jj = j-int(Nny/2-ny/2);
        
        gamma1[ii+nx*jj] = float( realsp[i+Nnx*j]/Nnx/Nny);
        
      }
    }
  }
  // gamma2
  {
    
    // build modes for each pixel in the fourier space
    for( int i=0; i<Nnx/2+1; i++ ){
      double kx=double(i);
      kx=kx*2.*M_PI/Nboxlx;
      for( int j=0; j<Nny; j++ ){
        double ky=(j<Nny/2)?double(j):double(j-Nny);
        ky=ky*2.*M_PI/Nboxly;
        
        // gamma
        fft[i+(Nnx/2+1)*j][0] = kx*ky*fphi[i+(Nnx/2+1)*j][0];
        fft[i+(Nnx/2+1)*j][1] = kx*ky*fphi[i+(Nnx/2+1)*j][1];
        
      }
    }
    
    fftw_execute( pp );
    
    gamma2.resize(nx*ny);
    
    for( int j=Nny/2-ny/2; j<Nny/2+ny/2; j++ ){
      for( int i=Nnx/2-nx/2; i<Nnx/2+nx/2; i++ ){
        int ii = i-int(Nnx/2-nx/2);
        int jj = j-int(Nny/2-ny/2);
        
        gamma2[ii+nx*jj] = float(-realsp[i+Nnx*j]/Nnx/Nny);
        
      }
    }
  }
  
  
  // std:: cout << " remapping the map in the original size " << std:: endl;
  delete[] fft;
  delete[] realsp;
  delete[] fphi;
#endif
}

#endif
