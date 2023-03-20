//
//  gridmap.h
//  GLAMER
//
//  Created by bmetcalf on 8/19/14.
//
//

#ifndef __GLAMER__gridmap__
#define __GLAMER__gridmap__

#include <iostream>
#include "lens.h"
#include "point.h"
#include "Tree.h"
#include "source.h"
#include <mutex>

/** 
 * \brief A simplified version of the Grid structure for making non-adaptive maps of the lensing quantities (kappa, gamma, etc...)
 *
 *  GripMap is faster and uses less memory than Grid.  It does not construct the tree structures for the points 
 *  and thus cannot be used for adaptive mapping or image finding.
 *
 *  The distance between the left (lower) most and right (upper) most ray is range so the resolution is range/(N-1).
 *  The lower left pixel is at center[]-0.5*range and the upper right is at center[]+0.5*range
 */

struct GridMap{
  
	GridMap(LensHndl lens,unsigned long N1d,const double center[2],double range);
  GridMap(LensHndl lens ,unsigned long Nx ,const PosType center[2] ,PosType rangeX ,PosType rangeY);
	~GridMap();
  
    /// reshoot the rays for example when the source plane has been changed
  GridMap ReInitialize(LensHndl lens);
  /**
   * \brief Recalculate surface brightness at every point without changing the positions of the gridmap or any lens properties.
   *
   *  Recalculate the surface brightness at all points on the gridmap.
   * This is useful when changing the source model while preserving changes in the grid.
   * Both i_tree and s_tree are both changed although only s_tree shows up here.
   *
   * returns the sum of the surface brightnesses
   */
  double RefreshSurfaceBrightnesses(SourceHndl source);
  /**
   Oversample some pixels where the usrface brightness is not smooth and update surface brighnesses to be the average inside the pixel.
   
   May be slow.
   */
  double AdaptiveRefreshSurfaceBrightnesses(Lens &lens,Source &source);
  
 /**
   * \brief Recalculate surface brightness just like GridMap::RefreshSurfaceBrightness but
   * the new source is added to any sources that were already there.
   *
   * returns total flux from the new source
   */
  double AddSurfaceBrightnesses(SourceHndl source);

  /// get the image point for a index number
  Point_2d image_point(size_t index){return i_points[index];}
  /// get the image point for a index number
  Point_2d source_point(size_t index){return s_points[index];}
 
  void ClearSurfaceBrightnesses();
  void assertNAN(); // check for nan in surface prightness
	size_t getNumberOfPoints() const {return Ngrid_init*Ngrid_init2;}
  
	/// return initial number of grid points in each direction
	int getInitNgrid(){return Ngrid_init;}
	/// return initial range of gridded region.  This is the distance from the first ray in a row to the last (unlike PixelMap)
	double getXRange(){return x_range;}
	double getYRange(){return x_range*axisratio;}
  /// resolution in radians, this is range / (N-1)
  double getResolution(){return x_range/(Ngrid_init-1);}
  
   /// make pixel map of lensing quantities at the resolution of the GridMap
  PixelMap writePixelMap(LensingVariable lensvar);
   /// fits output of lensing quantities at the resolution of the GridMap
  void writeFits(LensingVariable lensvar,std::string filensame);

  void writePixelMapUniform(PixelMap &map,LensingVariable lensvar);
  void writeFitsUniform(const PosType center[],size_t Nx,size_t Ny,LensingVariable lensvar,std::string filename);
  PixelMap writePixelMapUniform(const PosType center[],size_t Nx,size_t Ny,LensingVariable lensvar);

  /// this will make a fits map of the grid as is.
  void writeFitsUniform(
                        LensingVariable lensvar    ///< quantity to be output
                        ,std::string filename     ///< name of output fits file
                        ){
    PixelMap map = writePixelMap(lensvar);
    map.printFITS(filename);
  }
  
  /// returns a PixelMap with the flux in pixels at a resolution of res times the original resolution
  PixelMap getPixelMapFlux(int res) const;
  /// update a PixelMap with the flux in pixels at a resolution of res times the original resolution.
  /// The map must have precisely the right size and center to match or an exception will be thrown.
  /// Constructing the map with PixelMap getPixelMapFlux(int res) will insure that it does.
  void getPixelMapFlux(PixelMap &map) const;
  
  /// returns the area (radians^2) of the region with negative magnification at resolution of fixed grid
  PosType EinsteinArea() const;

  /// flux weighted magnification with current surface brightness averaged on the image plane
  PosType magnification() const;
  PosType magnification2() const;
  /// returns centroid of flux on the grid
  Point_2d centroid() const;
  
  Point_2d getCenter(){return center;}
  
  Point * operator[](size_t i){return i_points + i;};
  
  GridMap(GridMap &&grid){
    *this = std::move(grid);
  }
  
  GridMap & operator=(GridMap &&grid){
    Ngrid_init = grid.Ngrid_init;
    Ngrid_init2 = grid.Ngrid_init2;
    pointID = grid.pointID;
    axisratio = grid.axisratio;
    x_range = grid.x_range;
    
    i_points = grid.i_points;
    grid.i_points = nullptr;
    s_points = grid.s_points;
    grid.s_points = nullptr;
    
    center = grid.center;
    
    point_factory = std::move(grid.point_factory);

    return *this;
  }

  struct Triangle{
    Triangle(size_t i,size_t j,size_t k){
      index[0] = i;
      index[1] = j;
      index[2] = k;
    }
    size_t index[3];
    size_t & operator[](int i){return index[i];}
  };
  
    
  /*** \brief Returns a list of  RAYs from a set of source positions.
   
   The image positions are found in parallel by the triangle method.  The order of the
   output rays will be the same as the sources with multiple images consecutive.
   */
  std::list<RAY> find_images(std::vector<Point_2d> &ys) const{
    int nthreads = Utilities::GetNThreads();
 
    long N =  ys.size();
    long n = (int)(N/nthreads + 1);
    
    std::vector<std::list<RAY>> v_of_lists(nthreads);
    
    std::vector<std::thread> thr;
    long m = 0,i=0;
    //for(int i = 0; i < nthreads ;++i){
    for(int i=0 ; i<nthreads ; ++i ){
 
      if(m > N-n ) n = N-m;
      thr.push_back(std::thread(
                                &GridMap::_find_images_
                                ,this
                                ,ys.data() + m
                                ,n
                                ,std::ref(v_of_lists[i])
                                )
                    );
      
      m += n;
    }
 
    for(auto &t : thr) t.join();
  
    //join ray lists
    for(int i = 1 ; i<nthreads ; ++i) v_of_lists[0].splice(v_of_lists[0].end(),v_of_lists[i]);
  
    return v_of_lists[0];
  }

  /** find all images by triangle method
   */
  void find_images(Point_2d y
                   ,std::vector<Point_2d> &image_points  /// positions of the images limited by resolution of the gridmap
                   ,std::vector<Triangle> &triangles     /// index's of the points that form the triangles that the images are in
                   ) const {
    
    image_points.clear();
    triangles.clear();
    
    size_t k;
    size_t k1;
    size_t k2;
    size_t k3;
    
    int sig1,sig2,sig3,sig_sum;

    for(size_t i=0 ; i< Ngrid_init-1 ; i++){
      for(size_t j=0 ; j< Ngrid_init2-1 ; j++){
        k = i + j*Ngrid_init;
        k1 = k + Ngrid_init + 1;
        
        k2 = k + 1;
        k3 = k + Ngrid_init;

        sig1 = sign( (y-s_points[k])^(s_points[k1]-s_points[k]) );
        sig2 = sign( (y-s_points[k1])^(s_points[k2]-s_points[k1]) );
        sig3 = sign( (y-s_points[k2])^(s_points[k]-s_points[k2]) );
        
        sig_sum = sig1 + sig2 + sig3;
        if(abs(sig_sum) == 3){ // inside
          image_points.push_back( ( i_points[k] + i_points[k1] + i_points[k2] )/3  );
          triangles.push_back(Triangle(k,k1,k2));
        }else if(abs(sig_sum) == 2){ // on edge
          if(sig_sum > 0){
            if(sig1 == 0){
              image_points.push_back( ( i_points[k] + i_points[k1])/2  );
              triangles.push_back(Triangle(k,k1,k2));
            }else if(sig2 == 0){
              image_points.push_back( ( i_points[k1] + i_points[k2])/2  );
              triangles.push_back(Triangle(k,k1,k2));
            }else{
              image_points.push_back( ( i_points[k2] + i_points[k])/2  );
              triangles.push_back(Triangle(k,k1,k2));
            }
          }
        }else if (sig1 == 0 && sig3 == 0){ // a vertex
          image_points.push_back( i_points[k] );
        }
       
        sig1 = sign( (y-s_points[k])^(s_points[k1]-s_points[k]) );
        sig2 = sign( (y-s_points[k1])^(s_points[k3]-s_points[k1]) );
        sig3 = sign( (y-s_points[k3])^(s_points[k]-s_points[k3]) );
        
        sig_sum = sig1 + sig2 + sig3;
        if(abs(sig_sum) == 3){ // inside
          image_points.push_back( ( i_points[k] + i_points[k1] + i_points[k3] )/3  );
          triangles.push_back(Triangle(k,k1,k3));
        }else if(abs(sig_sum) == 2){ // on edge
          if(sig_sum > 0){
            if(sig1 == 0){
              image_points.push_back( ( i_points[k] + i_points[k1] )/2  );
              triangles.push_back(Triangle(k,k1,k3));
            }else if(sig2 == 0){
              image_points.push_back( ( i_points[k1] + i_points[k3] )/2  );
              triangles.push_back(Triangle(k,k1,k3));
            }else{
              image_points.push_back( ( i_points[k3] + i_points[k] )/2  );
              triangles.push_back(Triangle(k,k1,k3));
            }
          }
        }
      }
    }
    
    return;
  }
  
  /** \brief add flux to the rays that are nearest to the source on the source plane for each image
  *
   * This uses GridMap::find_images to find the images.  It then finds the point that is closest to the source position.
   *  The flux is added to one point per image.  The total flux added is returned.  No further refinement of the grid is done
   *  so it is limited by the resolution of the GridMap.  Some spurious low magnification images can be found.
   */
  double AddPointSource(const Point_2d &y,double flux){
    std::vector<Point_2d> images;
    std::vector<Triangle> tri;

    find_images(y,images,tri);
    
    int n = images.size();
    
    double total_flux = 0;
    for(int i = 0 ; i<n ; ++i){
      
      int closest = 0;
      double d = (s_points[ tri[i][0] ] - y).length_sqr();
      double tmp = (s_points[ tri[i][1] ] - y).length_sqr();
      if(tmp < d){ d=tmp; closest = 1;}
      tmp = (s_points[ tri[i][2] ] - y).length_sqr();
      if(tmp < d){ d=tmp; closest = 2;}
 
      size_t ii = tri[i][closest];
      total_flux += flux / fabs(i_points[ii].invmag());
      d = flux / fabs(i_points[ii].invmag()) / getResolution() / getResolution() ;
      i_points[ii].surface_brightness += d;
      s_points[ii].surface_brightness += d;
    }
    
    return total_flux;
  }
  
  void find_crit(std::vector<std::vector<Point_2d> > &points
                 ,std::vector<bool> &hits_boundary
                 ,std::vector<CritType> &crit_type
                 ){
    
    points.resize(0);
    
    size_t N = Ngrid_init * Ngrid_init2;
    std::vector<bool> bitmap(N);
    size_t count = 0;
    double eigenv[2];
 
    // find tangential critical curves
    for(size_t k = 0 ; k < N ; ++k){
      i_points[k].A.eigenv(eigenv);
      if(eigenv[1] < 0){
        bitmap[k] = true;
        ++count;
      } else {
        bitmap[k] = false;
      }
    }
    if(count>0) find_boundaries(bitmap,points,hits_boundary);
    crit_type.resize(points.size());
    for(CritType &b : crit_type) b = CritType::tangential;

    // find radial critical curves
    count=0;
    for(size_t k = 0 ; k < N ; ++k){
      i_points[k].A.eigenv(eigenv);
      if(eigenv[0] < 0 && eigenv[1] < 0){
        bitmap[k] = true;
        ++count;
      }else{
        bitmap[k] = false;
      }
    }
    if(count>0) find_boundaries(bitmap,points,hits_boundary,true);
    
    int k = crit_type.size();
    crit_type.resize(points.size());
    for(int i=k ; i<crit_type.size() ; ++i) crit_type[i] = CritType::radial;
  }
  
  /** finds ordered boundaries to regions where bitmap == true

   This can be used to find critical curves or contours.
   `bitmap` should be the same size as the `Gridmap`
   If the boundary curve  touches the edge of the `GridMap` it will be indicated in `hits_boundary` as
   `true`.
   
   Boundaries will never cross or lead off the grid.  On the edges they will leave the edge pixels out even if they should be in.  This is a technical comprimise.
  */
  void find_boundaries(std::vector<bool> &bitmap  // = true inside
                       ,std::vector<std::vector<Point_2d> > &points
                       ,std::vector<bool> &hits_edge
                       //,std::vector<CritType> &crit_type
                       ,bool add_to_vector=false
                       ){
    
    size_t nx = Ngrid_init;
    size_t ny = Ngrid_init2;
    size_t n = nx*ny;
    size_t ncells = nx*ny;
    
    std::vector<bool> not_used(bitmap.size(),true);
    
    assert(bitmap.size()==ncells);
    
    // pad edge of field with bitmap=false
    for(size_t i=0 ; i<nx ; ++i) bitmap[i]=false;
    size_t j = nx*(ny-1);
    for(size_t i=0 ; i<nx ; ++i) bitmap[i + j]=false;
    for(size_t i=0 ; i<ny ; ++i) bitmap[i*nx]=false;
    j = nx-1;
    for(size_t i=0 ; i<ny ; ++i) bitmap[j + i*nx]=false;

    std::list<std::list<Point_2d>> contours;
    if(!add_to_vector){
      hits_edge.resize(0);
      //crit_type.resize(0);
    }
    
    bool done = false;
    size_t kfirst_in_bound = 0;
    while(!done){
      // find first cell in edge
      size_t k=0;
      int count;
      for( k = kfirst_in_bound + 1; k < n - nx ; ++k){
        if(k % nx != nx-1){ // one less cells than points
          count = 0;
          if(bitmap[k]) ++count;
          if(bitmap[k+1]) ++count;;
          if(bitmap[k + Ngrid_init]) ++count;;
          if(bitmap[k + Ngrid_init + 1]) ++count;
          if(count > 0
             && count < 4
             && not_used[k]
             ) break;
        }
      }
      
      kfirst_in_bound = k;
      
      if(k == n-nx){
        done=true;
      }else{ // found an edge
        
        contours.resize(contours.size() + 1);
        std::list<Point_2d> &contour = contours.back();
        hits_edge.push_back(false);
        
        /* find type of critical curve
        if(i_points[k].inverted()
           || i_points[k+1].inverted()
           || i_points[k+nx].inverted()
           || i_points[k+nx+1].inverted()
           ){
          crit_type.push_back(CritType::radial);
        }else{
          crit_type.push_back(CritType::tangential);
        }
        */
        
        int face_in=0;
        int type;
        size_t n_edge = 0;
        // follow edge until we return to the first point
        while(k != kfirst_in_bound || n_edge==0){
          
          if(k%nx == 0 || k%nx == nx-2) hits_edge.back() = true;
          if(k/nx == 0 || k/nx == ny-2) hits_edge.back() = true;
          
          not_used[k] = false;
          
          ++n_edge;
          type = 0;
          // find type of cell
          if(bitmap[k]) type +=1;
          if(bitmap[k+1]) type += 10;
          if(bitmap[k + Ngrid_init]) type += 100;
          if(bitmap[k + Ngrid_init + 1]) type += 1000;
          
          if(type == 0 || type == 1111){  // all in or all out
            throw std::runtime_error("off edge!!");
          }else if(type == 1 || type == 1110){ // lower left only
            
            if(face_in==0){
              contour.push_back( (i_points[k] + i_points[k+1]) / 2 );
              face_in=1;
              k -= nx;
            }else{
              contour.push_back( (i_points[k] + i_points[k+nx]) / 2 );
              face_in=2;
              k -= 1;
            }
            
          }else if(type == 10 || type == 1101){ // lower right only
            
            if(face_in==2){
              contour.push_back( (i_points[k] + i_points[k+1]) / 2 );
              face_in=1;
              k -= nx;
            }else{
              contour.push_back( (i_points[k+nx] + i_points[k+1]) / 2 );
              face_in=0;
              k += 1;
            }
            
          }else if(type == 100 || type == 1011){ // upper left only
            
            if(face_in==0){
              contour.push_back( (i_points[k+nx] + i_points[k+nx+1]) / 2 );
              face_in=3;
              k += nx;
            }else{
              contour.push_back( (i_points[k] + i_points[k+nx]) / 2 );
              face_in=2;
              k -= 1;
            }
            
          }else if(type == 1000 || type == 111){ // upper right only
            
            if(face_in==1){
              contour.push_back( (i_points[k+1] + i_points[k+nx+1]) / 2 );
              face_in=0;
              k += 1;
            }else{
              contour.push_back( (i_points[k+nx] + i_points[k+nx+1]) / 2 );
              face_in=3;
              k += nx;
            }
            
          }else if(type == 11 || type == 1100){ // lower two
            
            if(face_in==0){
              contour.push_back( (i_points[k+1] + i_points[k+nx+1]) / 2 );
              k += 1;
            }else{
              contour.push_back( (i_points[k] + i_points[k+nx]) / 2 );
              k -= 1;
            }
            
          }else if(type == 1010 || type == 101){ // right two
            
            if(face_in==1){
              contour.push_back( (i_points[k] + i_points[k+1]) / 2 );
              k -= nx;
            }else{
              contour.push_back( (i_points[k+nx] + i_points[k+nx+1]) / 2 );
              k += nx;
            }
            
          }else if(type == 1001){ // lower left upper right
            
            if(face_in==0){
              contour.push_back( (i_points[k+nx] + i_points[k+nx+1]) / 2 );
              face_in=3;
              k += nx;
            }else if(face_in==1){
              contour.push_back( (i_points[k] + i_points[k+nx]) / 2 );
              face_in=2;
              k -= 1;
            }else if(face_in==2){
              contour.push_back( (i_points[k] + i_points[k+1]) / 2 );
              face_in=1;
              k -= nx;
            }else{
              contour.push_back( (i_points[k+nx] + i_points[k+1]) / 2 );
              face_in=0;
              k += 1;
            }
            
          }else if(type == 110){ // upper left lower right
            
            if(face_in==0){
              contour.push_back( (i_points[k] + i_points[k+1]) / 2 );
              face_in=1;
              k -= nx;
            }else if(face_in==1){
              contour.push_back( (i_points[k+1] + i_points[k+nx+1]) / 2 );
              face_in=0;
              k += 1;
            }else if(face_in==2){
              contour.push_back( (i_points[k + nx] + i_points[k+nx+1]) / 2 );
              face_in=3;
              k += nx;
            }else{
              contour.push_back( (i_points[k] + i_points[k+nx]) / 2 );
              face_in=2;
              k -= 1;
            }
          }
        }
      }
    }
    
    int offset = 0;
    if(!add_to_vector){
      points.resize(contours.size());
    }else{
      offset = points.size();
      points.resize(points.size() + contours.size());
    }
    
    // copy lists of points into vectors
    int i=0;
    for(auto &c: contours){
      points[offset+i].resize(c.size());
      size_t j=0;
      for(auto &p: c){
        points[offset+i][j] = p;
        ++j;
      }
      ++i;
    }
  }
  
private:
  
  // cluge to make compatible with old method of producing points
  Point * NewPointArray(size_t N){
    Point * p = point_factory(N);
    p[0].head = N;
    for(size_t i=1; i < N ; ++i) p[i].head = 0;
    return p;
  }
  
  void xygridpoints(Point *points,double range,const double *center,long Ngrid
                    ,short remove_center);
  
	/// one dimensional size of initial grid
	int Ngrid_init;
  int Ngrid_init2;
  
  unsigned long pointID;
  PosType axisratio;
  PosType x_range;
  void writePixelMapUniform_(Point* points,size_t size,PixelMap *map,LensingVariable val);
  
  Point *i_points;
  Point *s_points;
  Point_2d center;
  
  bool to_refine(long i,long j,double total,double f) const ;
  static std::mutex grid_mutex;
  
  MemmoryBank<Point> point_factory;
  
  void _find_images_(Point_2d *ys,long Nys,std::list<RAY> &rays) const{
    
    std::vector<Point_2d> x;
    std::vector<Triangle> triangles;
    rays.resize(0);
    auto itr = rays.begin();
    for(long k=0 ; k<Nys ; ++k){
      find_images(ys[k],x,triangles);
      rays.emplace_back();
      RAY &ray = rays.back();
      for(int i=0 ; i < x.size() ; ++i){
        ray.x=x[i];
        ray.y = ys[k];
        ray.A *= 0;
        for(int j=0 ; j<3 ; ++j){
          ray.A = ray.A  + i_points[triangles[i][j]].A;
        }
        ray.A /= 3;
      }
    }
    
    return;
  }

};

#endif // defined(__GLAMER__gridmap__)
