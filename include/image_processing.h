/*
 * image_processing.h
 *
 *  Created on: Feb 28, 2010
 *      Author: R.B. Metcalf
 */

#ifndef IMAGE_PROCESSING_H_
#define IMAGE_PROCESSING_H_

#include <complex>
#include <vector>
#include <tuple>

#include "point.h"
#include "image_info.h"
//#include "Tree.h"
//#include "utilities_slsim.h"
//#include "image_processing.h"
#include "source.h"

#ifdef ENABLE_FFTW
#include "fftw3.h"
#endif

// forward declaration
struct Grid;
struct GridMap;
//class Source;

/// These are partial units for the pixel map that can be used to ensure consistency.  For example, maps with different units cannot be added together.  default: ndef
enum class PixelMapUnits {
  ndef     // not defined
  ,surfb   // ergs / s / cm**2
  ,count_per_sec
  ,ADU   // Analogue-to-Digital Units
  ,mass
  ,mass_density
};

std::string to_string(PixelMapUnits unit);

/**
 * \brief Image structure that can be manipulated and exported to/from fits files.
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
  PixelMap(const PixelMap& pmap,long nx,long ny, std::size_t Npixels);
	PixelMap(const double* center, std::size_t Npixels, double resolution,PixelMapUnits u = PixelMapUnits::ndef);
	PixelMap(const double* center, std::size_t Nx, std::size_t Ny, double resolution,PixelMapUnits u = PixelMapUnits::ndef);
	PixelMap(std::string fitsfilename
           ,double resolution = -1,PixelMapUnits u = PixelMapUnits::ndef);
  ~PixelMap(){
    map.resize(0);
  };
	
  PixelMap& operator=(const PixelMap &other);
  PixelMap& operator=(PixelMap &&other);
  
  void ChangeUnits(PixelMapUnits u){units=u;}
  PixelMapUnits getUnits() const {return units;}
  
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
	
  /// returns right accention of center
  double getRA(){return RA;}
  /// returns declination of center
  double getDEC(){return DEC;}
  
  /// set the coordinates of center
  void setRAandDec(double RAin,double DECin){
    RA = RAin;
    DEC = DECin;
  }
  
  // Zero the whole map
  void Clean(){for(auto &a : map) a = 0;}

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
  void AddPointSource(const Point_2d &x,double flux);
  
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

  /// paste a PixelMap on with the lower left pixel match to [nx,ny] of this
  void paste(const PixelMap& pmap,long nx,long ny);

  /** \brief copy a PixelMap that must be the same without creating a new one..
   
   This avoids calling a any constructor or destructor.
   **/
  void duplicate(const PixelMap& pmap);

  /// Adds source to map.  This version breaks pixels up into blocks and does them in seporate threads.
  template <typename T>
  PosType AddSource_parallel(T &source,int oversample){
    Point_2d s_center;
    source.getTheta(s_center);
    
    if( (s_center[0] + source.getRadius() ) < map_boundary_p1[0] ) return 0.0;
    if( (s_center[0] - source.getRadius() ) > map_boundary_p2[0] ) return 0.0;
    if( (s_center[1] + source.getRadius() ) < map_boundary_p1[1] ) return 0.0;
    if( (s_center[1] - source.getRadius() ) > map_boundary_p2[1] ) return 0.0;
    
    int nthreads = Utilities::GetNThreads();
    PosType totals[nthreads];
    std::vector<std::thread> thr;
    
    size_t block = map.size()/nthreads;
    for(int i = 0; i < nthreads ;++i){
      thr.push_back(std::thread(&PixelMap::addsource_<T>,this
                      ,i*block,std::min((i+1)*block-1,map.size()-1)
                      ,oversample,std::ref(source)
                      ,std::ref(totals[i])));
    }
 
    for(int ii=0;ii < nthreads;++ii) thr[ii].join();
    
    PosType total =0;
    for(int ii=0;ii < nthreads;++ii) total += totals[ii];
    
    return total;
  }

  void AddCurve(ImageInfo *curve,double value);
  void AddCurve(Kist<Point> *imagekist,PosType value);
  void AddCurve(std::vector<Point_2d> &curve,double value);
  void AddCurve(std::vector<RAY> &curve,double value);
  
  /// simple line
	void drawline(double x1[],double x2[],double value,bool add);
  /// line by Bresenham's line algorithm
  void DrawLine(long x0,long x1,long y0,long y1,double value,bool add);
  void DrawLineGS(long x0,long x1,long y0,long y1,double value,bool add);
  void drawcircle(PosType r_center[],PosType radius,PosType value);
  void drawdisk(PosType r_center[],PosType radius,PosType value,int Nstrip);
	void AddGrid(const Grid &grid,double value = 1.0);
  void AddGrid(const Grid &grid,LensingVariable val);
  
	void Renormalize(double factor);
	void AddValue(std::size_t i, double value);
	void AssignValue(std::size_t i, double value);
	void printASCII() const;
	void printASCIItoFile(std::string filename) const;
	void printFITS(std::string filename,bool Xflip = false, bool verbose = false);
  void printFITS(std::string filename,std::vector< std::tuple<std::string,double,std::string> > &extra_header_info, bool verbose);
  
  /// This overides all header information and relaces it with the inputs. Meant for making a modified copy
  void printFITS(std::string filename  /// file to create
                ,std::vector<std::string> &headercards /// header information in cfitsio "card" format
                 );

	void smooth(double sigma);

	inline double getValue(std::size_t i) const { return map[i]; }
  inline double & operator[](std::size_t i) { return map[i]; };
  double operator[](std::size_t i) const { return map[i]; };
  inline double operator()(std::size_t i) const { return map[i]; };
  inline double operator()(std::size_t i,std::size_t j) const { return map[i + Nx*j]; };
  inline double & operator()(std::size_t i,std::size_t j) { return map[i + Nx*j]; };
  
	PixelMap& operator+=(const PixelMap& rhs);
  void operator+=(double f){map +=f;};
  //friend PixelMap operator+(const PixelMap&, const PixelMap&);
  PixelMap operator+(const PixelMap&) const;

	PixelMap& operator-=(const PixelMap& rhs);
  //friend PixelMap operator-(const PixelMap&, const PixelMap&);
  PixelMap operator-(const PixelMap&) const;

	PixelMap& operator*=(const PixelMap& rhs);
	//friend PixelMap operator*(const PixelMap&, const PixelMap&);
  PixelMap operator*(const PixelMap& a) const;
  PixelMap operator/(const PixelMap& a) const;

	PixelMap& operator*=(PosType b);
  PixelMap operator*(PosType b) const;

	std::valarray<double>& data() { return map; }
	
  /// Check whether two PixelMaps agree in their physical dimensions.
	bool agrees(const PixelMap& other) const;
	
	friend void swap(PixelMap&, PixelMap&);
  static void swap(PixelMap&, PixelMap&);
  
  /// return average pixel value
  PosType ave() const {return map.sum()/map.size();}
  /// return sum of all pixel values
  PosType sum() const {return map.sum();};
  /// Total number of pixels
  size_t size(){return map.size();}
  double max() const{ return map.max(); }
  double min() const{ return map.min(); }

  void FindArc(PosType &radius,PosType *xc,PosType *arc_center,PosType &arclength,PosType &width
                         ,PosType threshold);
  
  /// get the index for a position, returns -1 if out of map, this version returns the 2D grid coordinates
  long find_index(PosType const x[],long &ix,long &iy) const;
  /// get the index for a position, returns -1 if out of map
  long find_index(PosType const x[]) const;
  
  /// get the index for a position, returns -1 if out of map, this version returns the 2D grid coordinates
  long find_index(PosType x,PosType y,long &ix,long &iy) const;
  /// get the index for a position, returns -1 if out of map
  long find_index(PosType x,PosType y) const;
  
  /// get the index for a position, returns -1 if out of map
  void find_position(PosType x[],std::size_t const index) const;
  /// get the index for a position, returns -1 if out of map
  void find_position(PosType x[],std::size_t const ix,std::size_t const iy) const;
  
  /// rotate and scale the image while keeping pixels, resoluiton
  PixelMap rotate(
                  PosType theta  /// counter-clockwise rotation (radians)
                  ,PosType scale=1  /// scale <1 shrinks it
                  );
  
  /// interpolate to point x[]
  PosType linear_interpolate(PosType x[]);
  
  /// draw a grid on the image that divides the each demension into N cells
  void drawgrid(int N,PosType value);
  void drawPoints(std::vector<Point *> points,PosType size,PosType value);
  
  void drawPoints(std::vector<Point> points,PosType size,PosType value);
  
  void drawCurve(std::vector<Point *> points,PosType value){
    for(int i=0;i<points.size()-1;++i) drawline(points[i]->x,points[i+1]->x,value,false);
  }
  void drawCurve(std::vector<Point> points,PosType value){
    for(int i=0;i<points.size()-1;++i) drawline(points[i].x,points[i+1].x,value,false);
  }
  void drawPoints(std::vector<Point_2d> points,PosType size,PosType value);
  void drawCurve(std::vector<Point_2d> points,PosType value){
    for(int i=0;i<points.size()-1;++i) drawline(points[i].x,points[i+1].x,value,false);
  }
  /// Draw a rectangle
  void drawSquare(PosType p1[],PosType p2[],PosType value);
  void drawBox(PosType p1[],PosType p2[],PosType value,int Nstrip);

  /// Find the power spectrum of the map
  void PowerSpectrum(std::vector<PosType> &power_spectrum   /// output power spectrum
                     ,std::vector<PosType> &lvec            /// output l values of bands
                     ,bool overwrite = true                 /// if false add power to existing power_spectrum (used for averaging over many fields
                     ){
    
    if(power_spectrum.size() != lvec.size()) throw std::invalid_argument("these must be the same size");
    
    if(overwrite){
      Utilities::powerspectrum2d(map,Nx,Ny,rangeX,rangeY, lvec, power_spectrum);
    }else{
      std::vector<PosType> tmp_power(power_spectrum.size());
      Utilities::powerspectrum2d(map,Nx,Ny,rangeX,rangeY,lvec,tmp_power);
      for(size_t ii=0;ii<power_spectrum.size();++ii) power_spectrum[ii] += tmp_power[ii];
    }
  }

  /// Find the power spectrum of the map
  void PowerSpectrum(std::vector<PosType> &power_spectrum   /// output power spectrum
                     ,const std::vector<PosType> &lbins            /// output l values of bands
                     ,std::vector<PosType> &lave            /// output l values of bands
                     ,bool overwrite = true                 /// if false add power to existing power_spectrum (used for averaging over many fields
                     ){
    
    if(overwrite){
      Utilities::powerspectrum2dprebin(map,Nx,Ny,rangeX,rangeY,lbins,power_spectrum,lave);
    }else{
      if(power_spectrum.size() != lbins.size()-1) throw std::invalid_argument("these must be the same size");
      std::vector<PosType> tmp_power(power_spectrum.size());
      Utilities::powerspectrum2dprebin(map,Nx,Ny,rangeX,rangeY,lbins,tmp_power,lave);
      for(size_t ii=0;ii<power_spectrum.size();++ii) power_spectrum[ii] += tmp_power[ii];
    }
  }

  void AdaptiveSmooth(PosType value){
    std::valarray<double> tmp = Utilities::AdaptiveSmooth(data(),Nx,Ny,value);
    map = tmp;
  }
  
  /// returns a vector of  contour curves
  void find_contour(double level
                    ,std::vector<std::vector<Point_2d> > &points
                    ,std::vector<bool> &hits_edge
                    ) const;
  /** find all the points above level divided into seprated groups
   
   Groups with points are connected regions above level.
   Groups without points are regions surrounded by regions above level, i.e. holes.
   */
  
  void find_islands_holes(double level,
                    std::vector<std::vector<size_t> > &points
                    ) const;
  
  /** This applies a definition of lesning for a resolved source based on that of Sonnenfeld et al. 2023
   
   This is meant to be used on a signal-to-noise map of the lensed source only.
   
   The definition does the following tests:
   1) define a footprint at surface brightness level `pixel_threshold`
   2) romove disconected region that have signal-to-noise below `min_sn_per_image `
   3) if there are multiple images or a hole in the image `lens_TF =true`
   4) if there is one image with no hole increase the threshold and apply 3) until the it is classified as a
          lens or it reaches the maximum surface brightness level
   
   */
  void lens_definition(double min_sn_per_image                 /// signal-to-noise required for a seporate image (ex 10)
                       ,double pixel_threshold                 /// signal-to-noise threshold that defines the footprint (ex. 2)
                       ,int &Nimages                           /// the number of images
                       ,double &total_sig_noise_source         /// gives the total signal-to-noise of all images
                       ,std::vector<size_t> &maxima_indexes    /// index of maxima
                       ,std::vector<std::vector<size_t> > &image_points
                       ,bool &lens_TF                          /// whether it passes lens diffintion
                       ,double &level                          /// levels on which the multiple images are defined
                       ,size_t &n_pix_in_source                 /// number of pixels in footprint
                       ,bool verbose = false
                       );
 
  
  /// find maxima that are above minlevel
  std::vector<size_t> maxima(double minlevel) const;
  
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
  
  /// recenter the map without changing anything else.
  
  void recenter(PosType newcenter[2] /// in radians
                 );
  void recenter(Point_2d newcenter /// in radians
                 );
  
  /** \brief convolve the image with a kernel.
   
   It is assumed that the size of the kernel is much smaller than the image and
   that the kernal has the same pixel size as the image.
   **/
  void convolve(PixelMap &kernel,long center_x = 0,long center_y = 0);
 
  /** \brief Creates a PixelMap with a lower resolution.
   *  The value of the pixels are added for the new pixels.
   *   If n does not go into the orginial number of pixels evenly the right (top) redge is dropped.
   */
  PixelMap downsize(int n /// number of pixels each direction added into each new pixel
                    );


  /// add a heaader keyword that will appear in fits output
   void addheader(std::string label,long value,std::string comment){
     bool found = false;
     for(auto  &a : headers_long){
       if(std::get<0>(a) == label){
         std::get<1>(a)=value;
         std::get<2>(a)=comment;
         found = true;
         break;
       }
     }
     if(!found) headers_long.push_back(std::make_tuple(label,value,comment));
   }
   void addheader(std::string label,size_t value,std::string comment){
     bool found = false;
     for(auto  &a : headers_long){
       if(std::get<0>(a) == label){
         std::get<1>(a)=value;
         std::get<2>(a)=comment;
         found = true;
         break;
       }
     }
     if(!found) headers_long.push_back(std::make_tuple(label,value,comment));
   }
   void addheader(std::string label,float value,std::string comment){
     bool found = false;
     for(auto  &a : headers_float){
       if(std::get<0>(a) == label){
         std::get<1>(a)=value;
         std::get<2>(a)=comment;
         found=true;
         break;
       }
     }
     if(!found) headers_float.push_back(std::make_tuple(label,value,comment));
   }
  void addheader(std::string label,double value,std::string comment){
    bool found = false;
    for(auto  &a : headers_float){
      if(std::get<0>(a) == label){
        std::get<1>(a)=value;
        std::get<2>(a)=comment;
        found=true;
        break;
      }
    }
    if(!found) headers_float.push_back(std::make_tuple(label,value,comment));
  }
  void addheader(std::string label,std::string value,std::string comment){
    bool found = false;
    for(auto  &a : headers_string){
      if(std::get<0>(a) == label){
        std::get<1>(a)=value;
        std::get<2>(a)=comment;
        found=true;
        break;
      }
    }
    if(!found) headers_string.push_back(std::make_tuple(label,value,comment));
  }

private:
  std::vector<std::tuple<std::string,float,std::string> > headers_float;
  std::vector<std::tuple<std::string,long,std::string> > headers_long;
  std::vector<std::tuple<std::string,std::string,std::string> > headers_string;

  std::valarray<double> map;
	long Nx;
	long Ny;
  double resolution,rangeX,rangeY,center[2];
  double RA=0,DEC=0; // optional coordinates of center
	double map_boundary_p1[2],map_boundary_p2[2];
  PixelMapUnits units= PixelMapUnits::ndef;
  
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
  
  // find if pixel k is in the curve which must be in pixel units, meant to be used in conjunction with Utilities::find_boundaries<>
  bool incurve(long k,std::vector<Point_2d> &curve) const;
};

enum class Telescope {Euclid_VIS,Euclid_Y,Euclid_J,Euclid_H,KiDS_u,KiDS_g,KiDS_r,KiDS_i,HST_ACS_I,CFHT_u,CFHT_g,CFHT_r,CFHT_i,CFHT_z};

enum class UnitType {counts_x_sec, flux} ;


class Obs{
public:
  
  Obs(size_t Npix_xx,size_t Npix_yy  /// number of pixels in observation
      ,double pix_size               /// pixel size (in rad)
      ,int oversample          /// oversampling for input image
      ,float seeing = 0 // seeing in arcsec
  );
  
  virtual ~Obs(){};
  
  size_t getNxInput() const { return Npix_x_input;}
  size_t getNyInput() const { return Npix_y_input;}

  size_t getNxOutput() const { return Npix_x_output;}
  size_t getNyOutput() const { return Npix_y_output;}

  std::valarray<double> getPSF(){return map_psf;}
  //void setPSF(std::string psf_file);
  void setPSF(std::string psf_file,double resolution=0);
  void setPSF(PixelMap &psf_map);
  /// rotate and scale the psf from the original
  void rotatePSF(double theta   /// counter-clockwise rotation (radians)
                 ,double scale_x=1  /// scale <1 shrinks it
                 ,double scale_y=1  /// scale <1 shrinks it
  );
  
  /// add two PSFs to simulate stacking
  void coaddPSF(double f         /// relative weight of the PSFs, 1 being equal weight
                ,double theta1   /// ratation of first PSF
                ,double theta2   /// rotation of second PSF
                ,double scale_x  /// scale <1 shrinks it
                ,double scale_y  /// scale <1 shrinks it
                );
 
  void ApplyPSF(PixelMap &map_in,PixelMap &map_out);
  float getPixelSize() const {return pix_size;}
  void setNoiseCorrelation(std::string nc_file);
  
  // virtual methods

  virtual void AddNoise(PixelMap &pmap
                        ,PixelMap &error_map
                        ,Utilities::RandomNumbers_NR &ran,bool cosmics) = 0;
  
  virtual void Convert(PixelMap &map_in
                           ,PixelMap &map_out
                           ,PixelMap &error_map
                           ,bool psf
                           ,bool noise
                           ,Utilities::RandomNumbers_NR &ran,bool cosmics) = 0;
  
  virtual float getBackgroundNoise() const = 0;

  /// convert using stan
  virtual double mag_to_counts(double m) const = 0;
  virtual double counts_to_mag(double flux) const = 0;
  virtual double zeropoint() const = 0;
  virtual void setZeropoint(double zpoint) = 0;


protected:

  double pix_size; // pixel size (in rad)
  void CorrelateNoise(PixelMap &pmap);
  float seeing;  // full-width at half maximum of the gaussian smoothing
  
  // the number of pixels in the real image
  size_t Npix_x_output,Npix_y_output;
  // the number of pixels in the oversamples image
  size_t Npix_x_input,Npix_y_input;
 
  float psf_oversample; // psf oversampling factor
  void downsample(PixelMap &map_in,PixelMap &map_out) const;  // downsize from Npix_input to Npix_output
  
  PixelMap map_scratch;
 
private:
  double input_psf_pixel_size;
  size_t side_ncorr; // pixels on a side of input noise correlation function

  void fftpsf();  // FFT the psf for later use
  std::valarray<double> map_psf;  // array of the point spread function
  std::valarray<double> map_psfo;  // initial array of the point spread function
  
  std::vector<std::complex<double> > fft_psf;
  std::vector<std::complex<double> > fft_padded;
  std::vector<double> image_padded;
  std::vector<double> sqrt_noise_power;  // stores sqrt root of power noise spectrum

  // size of borders for psf convolution
  size_t nborder_x = 0;
  size_t nborder_y = 0;
  
  // size of padded images
  size_t n_x = 0;
  size_t n_y = 0;

  fftw_plan image_to_fft;
  fftw_plan fft_to_image;

  //PixelMap noise_correlation;
  std::vector<std::complex<double> > noise_fft_image;
  std::vector<double> noise_in_zeropad;
  fftw_plan p_noise_r2c;
  std::vector<double> noise_image_out;
  fftw_plan p_noise_c2r;
};

/**
 * \brief It creates a realistic image from the output of a ray-tracing simulation.
 *
 * It translates pixel values in observed units (counts/sec), applies PSF and noise.
 * Input must be in ergs/(s*cm^2*Hz).
 *
 *  see https://www.ucolick.org/~bolte/AY257/s_n.pdf
 */
class ObsVIS : public Obs{
private:
  
  // standard from magnitude to e- per sec
  double zero_point = 24.4;
  //double sigma_back_per_qsrttime = 0.00267 * sqrt(5.085000000000E+03);
  
  //double gain = 11160; // e-/ADU (Analog Digital Units)
  //double exp_num = 4;
  //double exp_time = 2260.;  // seconds
  //double l = 7103.43;
  //double dl = 3318.28;
  //double seeing = 0.18;
  
  // derived parameters;
  double sigma_background;  // background var in ergs / cm^2 / s / Hz
  //double sb_to_e;  // approximate convertion between ergs / cm^2 / s and e-

  // adds random cosmic rays to the noise map
  void cosmics(PixelMap &error_map
                ,double inv_sigma2 // for one dither
                ,int nc // number of cosmics to be added
                ,Utilities::RandomNumbers_NR &ran) const ;
  
public:
  
  // exposure times are set to wide survey expectations
  ObsVIS(size_t Npix_x,size_t Npix_y
         ,int oversample
         ,double resolution = 0.1*arcsecTOradians
         //,double t = 5.085000000000E+03  // observation time in seconds. default is for SC8
  );
  
  ObsVIS(size_t Npix_x
         ,size_t Npix_y
         ,const std::vector<double> &exposure_times  // in seconds
         ,int oversample
         );
  
  ObsVIS(size_t Npix_x
         ,size_t Npix_y
         ,const std::vector<double> &exposure_times  // in seconds
         ,int oversample
         ,double resolution
         ,double background_sigma
         );
  ~ObsVIS(){};
  
  /// add poisson noise to an image that is in units of electrons
  void AddPoisson(PixelMap &pmap
                         ,Utilities::RandomNumbers_NR &ran
                     );
  
  /// Applies  noise (read-out + Poisson) on an image, returns noise map
  void AddNoise(PixelMap &pmap
                 ,PixelMap &error_map
                 ,Utilities::RandomNumbers_NR &ran
                ,bool cosmic=true);

  void Convert(PixelMap &map_in
               ,PixelMap &map_out
               ,PixelMap &error_map  // this is sigma
               ,bool psf
               ,bool noise
               ,Utilities::RandomNumbers_NR &ran
               ,bool cosmic=true);
  
 
  double mag_to_counts(double m) const{
    if(m == 100) return 0;
    return pow(10,-0.4*(m + zero_point));
  }
  double counts_to_mag(double flux) const{
    if(flux <=0) return 100;
    return -2.5 * log10(flux) - zero_point;
  }

  double zeropoint() const {return zero_point;}
  void setZeropoint(double zpoint){zero_point=zpoint;}
 
  /// returns std of pixels in e-
  float getBackgroundNoise() const {
    double dt = Utilities::vec_sum(t_exp);// 3 * t1 + t2;
    
    return sigma_background / sqrt(dt);
  }
private:
  //double t1;
  //double t2;
  std::vector<double> t_exp;  // exposure times

};

/** 
 * \brief It creates a realistic image from the output of a ray-tracing simulation.
 *
 * It translates pixel values in observed units (counts/sec), applies PSF and noise.
 * Input must be in ergs/(s*cm^2*Hz).
 *
 *  see https://www.ucolick.org/~bolte/AY257/s_n.pdf
 */
class Observation : public Obs
{
public:
	Observation(Telescope tel_name,double exposure_time,int exposure_num
              ,size_t Npix_x,size_t Npix_y, float oversample);
	Observation(float diameter, float transmission, float exp_time, int exp_num, float back_mag, float read_out_noise
              ,size_t Npix_x,size_t Npix_y,double pix_size,float seeing = 0.);
	Observation(float diameter, float transmission, float exp_time, int exp_num, float back_mag, float read_out_noise ,std::string psf_file,size_t Npix_x,size_t Npix_y,double pix_size, float oversample = 1.);
  
  Observation(float zeropoint_mag, float exp_time, int exp_num, float back_mag, float read_out_noise, size_t Npix_x,size_t Npix_y,double pix_size,float seeing=0);
  Observation(float zeropoint_mag, float exp_time, int exp_num, float back_mag, float read_out_noise ,std::string psf_file,size_t Npix_x,size_t Npix_y,double pix_size, float oversample = 1.);
  
  ~Observation(){};
  
  float getExpTime() const {return exp_time;}
	int getExpNum() const {return exp_num;}
	float getBackMag() const {return back_mag;}
   /// read-out noise in electrons/pixel
	float getRon() const {return read_out_noise;}
  /// seeing in arcsecs
	float getSeeing() const {return seeing;}
	float getZeropoint() const {return mag_zeropoint;}
  void setZeropoint(double zpoint){mag_zeropoint=zpoint;}
    /// pixel size in radians
  float getBackgroundNoise(float resolution, UnitType unit = UnitType::counts_x_sec) const;
  float getBackgroundNoise() const {return 0;};

  void AddNoise(PixelMap &pmap,PixelMap &error_map,Utilities::RandomNumbers_NR &ran,bool dummy);

  void Convert(PixelMap &map_in
               ,PixelMap &map_out
               ,PixelMap &error_map
               ,bool psf
               ,bool noise
               ,Utilities::RandomNumbers_NR &ran
               ,bool cosmic=false
               );
 
  /// returns factor by which code image units need to be multiplied by to get flux units
  //double flux_convertion_factor(){ return pow(10,-0.4*mag_zeropoint); }

  void setExpTime(float time){exp_time = time;}
  void setPixelSize(float pixel_size){pix_size=pixel_size;}
 
  double mag_to_counts(double m) const {
    if(m == 100) return 0;
    return pow(10,-0.4*(m - mag_zeropoint));
  }
  double counts_to_mag(double flux) const{
    if(flux <=0) return 100;
    return -2.5 * log10(flux) + mag_zeropoint;
  }
  double zeropoint() const{
     return mag_zeropoint;
   }
private:
  
	//float diameter;  // diameter of telescope (in cm)
	//float transmission;  // total transmission of the instrument
	float mag_zeropoint;  // magnitude of a source that produces one count/sec in the image
	float exp_time;  // total exposure time (in sec)
	int exp_num;  // number of exposures
	float back_mag;  // sky (or background) magnitude in mag/arcsec^2
	float read_out_noise;  // read-out noise in electrons/pixel
  float gain;
  
	bool telescope; // was the observation created from a default telescope?
  float e_per_s_to_ergs_s_cm2;  // e- / s   for zero magnitudes
  float background_flux;  // e- / s / arcsec

  void set_up();
    
  void ToCounts(PixelMap &pmap);
  void ToSurfaceBrightness(PixelMap &pmap);
  void ToADU(PixelMap &pmap);
  
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
