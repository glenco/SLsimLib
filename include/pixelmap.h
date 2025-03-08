//
//  pixelmap.h
//  GLAMER
//
//  Created by Robert Benton Metcalf on 26/04/24.
//

#ifndef pixelmap_h
#define pixelmap_h

#include <vector>
#include <tuple>
#include <type_traits>

#include "utilities_slsim.h"
#include "point.h"
#include "cpfits.h"

struct Grid;
struct GridMap;
class Source;


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
template <typename T=double>
class PixelMap
{
  
  static_assert( std::is_same<T, float>::value || std::is_same<T, double>::value
                , "PixeMap can only be instantated with a double or float template argument");
  
public:
  PixelMap(const PixelMap<T>& pmap, double res_ratio);
  PixelMap();
  PixelMap(const PixelMap<T>& other);
  PixelMap(PixelMap<T>&& other);
  PixelMap(const PixelMap<T>& pmap, const double* center, std::size_t Npixels);
  PixelMap(const PixelMap<T>& pmap,long nx,long ny, std::size_t Npixels);
  PixelMap(const double* center, std::size_t Npixels, double resolution,PixelMapUnits u = PixelMapUnits::ndef);
  PixelMap(const double* center
           , std::size_t Nx, std::size_t Ny
           , double resolution
           ,PixelMapUnits u = PixelMapUnits::ndef);
  //PixelMap(std::string fitsfilename
  //         ,std::string extension = ""
  //         );
  PixelMap(std::string fitsfilename   /// file name of fits file to be read
          ,double my_res=-1         /// resolution (rad) of fits image if not given in fits file, use default or -1 otherwise
          ,PixelMapUnits u = PixelMapUnits::ndef
          ,std::string extension = ""  /// fits extension 1,2,,...
  );
  
  // for conversion from double to float
  template<typename OtherT>
  PixelMap(const PixelMap<OtherT>&);
  
  ~PixelMap(){
    map.resize(0);
  };
 
  
  PixelMap<T>& operator=(const PixelMap<T> &other);
  PixelMap<T>& operator=(PixelMap<T> &&other);
  
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
  
  // coordinates of lower left corner 
  inline Point_2d getLLBoundary() const{ return Point_2d(map_boundary_p1[0],map_boundary_p1[1]); }
  
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
  void AddUniformImages(ImageInfo *imageinfo,int Nimages,T value);
  PosType AddSource(Source &source);
  /// Add a source to the pixel map by oversamples the source so that oversample^2 points within each pixel are averaged
  PosType AddSource(Source &source,int oversample);
  void AddPointSource(const Point_2d &x,T flux);
  
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
  
  /** \brief convolve the image with a kernel.
   
   It is assumed that the kernal has the same pixel size as the image.
   
   The center of the kernel is at N/2+1.  If it has an odd number  of pixels
   in both directions this will be the center.
   **/
  PixelMap<T> convolve(const PixelMap<T> &kernel);
  PixelMap<T> convolve2(const PixelMap<T> &kernel);

  
  /** \brief cut out a part of the PixelMap
   
   Inputs are in pixels. xmax and xmin are one after the disired region
   so that the output size is (xmax-xmin) x (ymax-ymin)
   **/
  PixelMap<T> cutout(long xmin,long xmax,long ymin,long ymax);

  /** \brief copy a PixelMap that must be the same without creating a new one..
   
   This avoids calling a any constructor or destructor.
   **/
  void duplicate(const PixelMap& pmap);

  /// Adds source to map.  This version breaks pixels up into blocks and does them in seporate threads.
  template <typename S>
  PosType AddSource_parallel(S &source,int oversample){
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
      thr.push_back(std::thread(&PixelMap<T>::addsource_<S>,this
                      ,i*block,std::min((i+1)*block-1,map.size()-1)
                      ,oversample,std::ref(source)
                      ,std::ref(totals[i])));
    }
 
    for(int ii=0;ii < nthreads;++ii) thr[ii].join();
    
    PosType total =0;
    for(int ii=0;ii < nthreads;++ii) total += totals[ii];
    
    return total;
  }

  void AddCurve(ImageInfo *curve,T value);
  void AddCurve(Kist<Point> *imagekist,T value);
  void AddCurve(std::vector<Point_2d> &curve,T value);
  void AddCurve(std::vector<RAY> &curve,T value);
  
  /// simple line
  void drawline(double x1[],double x2[],T value,bool add);
  /// line by Bresenham's line algorithm
  void DrawLine(long x0,long x1,long y0,long y1,T value,bool add);
  void DrawLineGS(long x0,long x1,long y0,long y1,T value,bool add);
  void drawcircle(PosType r_center[],PosType radius,PosType value);
  void drawdisk(PosType r_center[],PosType radius,PosType value,int Nstrip);
  void AddGrid(const Grid &grid,T value = 1.0);
  void AddGrid(const Grid &grid,LensingVariable val);
  
  void Renormalize(T factor);
  void AddValue(std::size_t i, T value);
  void AssignValue(std::size_t i, T value);
  void printASCII() const;
  void printASCIItoFile(std::string filename) const;
  void printFITS(std::string filename,bool Xflip = false, bool ctype = false, bool verbose = false);
  void printFITS(std::string filename,std::vector< std::tuple<std::string,double,std::string> > &extra_header_info, bool verbose);
  
  /// This overides all header information and relaces it with the inputs. Meant for making a modified copy
  void printFITS(std::string filename  /// file to create
                ,std::vector<std::string> &headercards /// header information in cfitsio "card" format
                 );

  void smooth(double sigma);

  inline T getValue(std::size_t i) const { return map[i]; }
  inline T & operator[](std::size_t i) { return map[i]; };
  inline T operator[](std::size_t i) const { return map[i]; };
  inline T & operator()(std::size_t i) { return map[i]; };
  inline T operator()(std::size_t i) const { return map[i]; };
  inline T operator()(std::size_t i,std::size_t j) const { return map[i + Nx*j]; };
  inline T & operator()(std::size_t i,std::size_t j) { return map[i + Nx*j]; };
  
  PixelMap& operator+=(const PixelMap& rhs);
  void operator+=(float f){map +=f;};
  void operator+=(double f){map +=f;};
  //friend PixelMap operator+(const PixelMap&, const PixelMap&);
  PixelMap operator+(const PixelMap&) const;

  PixelMap& operator-=(const PixelMap& rhs);
  //friend PixelMap operator-(const PixelMap&, const PixelMap&);
  PixelMap operator-(const PixelMap&) const;
  void operator-=(float f){map -=f;};
  void operator-=(double f){map -=f;};

  PixelMap& operator*=(const PixelMap& rhs);
  //friend PixelMap operator*(const PixelMap&, const PixelMap&);
  PixelMap operator*(const PixelMap& a) const;
  PixelMap operator/(const PixelMap& a) const;

  PixelMap& operator*=(PosType b);
  PixelMap operator*(PosType b) const;

  std::valarray<T>& data() { return map; }
  
  /// Check whether two PixelMaps agree in their physical dimensions.
  bool agrees(const PixelMap& other) const;
  
  //friend void swap(PixelMap&, PixelMap&);
  static void swap(PixelMap<T>&, PixelMap<T>&);
  
  /// return average pixel value
  PosType ave() const {return map.sum()/map.size();}
  /// return sum of all pixel values
  PosType sum() const {return map.sum();};
  /// Total number of pixels
  size_t size(){return map.size();}
  T max() const{ return map.max(); }
  T min() const{ return map.min(); }

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
  PixelMap<T> rotate(
                  PosType theta  /// counter-clockwise rotation (radians)
                  ,T scale=1  /// scale <1 shrinks it
                  );
  
  /// interpolate to point x[]
  T linear_interpolate(PosType x[]);
  
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
                     );

  /// Find the power spectrum of the map
  void PowerSpectrum(std::vector<PosType> &power_spectrum   /// output power spectrum
                     ,const std::vector<PosType> &lbins            /// output l values of bands
                     ,std::vector<PosType> &lave            /// output l values of bands
                     ,bool overwrite = true                 /// if false add power to existing power_spectrum (used for averaging over many fields
                     );

  void AdaptiveSmooth(PosType value);
  
  /// returns a vector of  contour curves
  void find_contour(T level
                    ,std::vector<std::vector<Point_2d> > &points
                    ,std::vector<bool> &hits_edge
                    ) const;
  /** find all the points above level divided into seprated groups
   
   Groups with points are connected regions above level.
   Groups without points are regions surrounded by regions above level, i.e. holes.
   */
  
  void find_islands_holes(T level,
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
  void lens_definition(T min_sn_per_image                 /// signal-to-noise required for a seporate image (ex 10)
                       ,T pixel_threshold                 /// signal-to-noise threshold that defines the footprint (ex. 2)
                       ,int &Nimages                           /// the number of images
                       ,T &total_sig_noise_source         /// gives the total signal-to-noise of all images
                       ,std::vector<size_t> &maxima_indexes    /// index of maxima
                       ,std::vector<std::vector<size_t> > &image_points
                       ,bool &lens_TF                          /// whether it passes lens diffintion
                       ,T &level                          /// levels on which the multiple images are defined
                       ,size_t &n_pix_in_source                 /// number of pixels in footprint
                       ,bool verbose = false
                       );
 
  
  /// find maxima that are above minlevel
  std::vector<size_t> maxima(T minlevel) const;
  
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
  
  /* \brief convolve the image with a kernel.
   
   It is assumed that the size of the kernel is much smaller than the image and
   that the kernal has the same pixel size as the image.
   **/
  //void convolve(PixelMap &kernel,long center_x = 0,long center_y = 0);
 
  /** \brief Creates a PixelMap with a lower resolution.
   *  The value of the pixels are added for the new pixels.
   *   If n does not go into the orginial number of pixels evenly the right (top) redge is dropped.
   */
  PixelMap<T> downsize(int n /// number of pixels each direction added into each new pixel
                    );

  /** \brief Makes a PixelMap with resolution 1/n of the original with the values linearly interpolated.
   
   This normalizes so that the sum of the pixels is approximately constant.
   */
  PixelMap interpolate(int n);

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

  std::valarray<T> map;
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

  template <typename S>
    void addsource_(size_t i1,size_t i2,int oversample,
                    S &source,
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

template<typename T>
void swap(PixelMap<T>& x, PixelMap<T>& y)
{
  using std::swap;
  
  swap(x.map,y.map);
  
  swap(x.Nx, y.Nx);
  swap(x.Ny, y.Ny);
  swap(x.resolution, y.resolution);
  swap(x.rangeX, y.rangeX);
  swap(x.rangeY, y.rangeY);
  
  swap(x.center[0], y.center[0]);
  swap(x.center[1], y.center[1]);
  
  swap(x.map_boundary_p1[0], y.map_boundary_p1[0]);
  swap(x.map_boundary_p1[1], y.map_boundary_p1[1]);
  swap(x.map_boundary_p2[0], y.map_boundary_p2[0]);
  swap(x.map_boundary_p2[1], y.map_boundary_p2[1]);
}

#endif /* pixelmap_h */
