//
//  gridmap.cpp
//  GLAMER
//
//  Created by bmetcalf on 8/19/14.
//
//

#include <utility>
#include <mutex>
#include <thread>
#include "gridmap.h"
#include "Tree.h"

std::mutex GridMap::grid_mutex;

/** 
 * \brief Constructor for initializing rectangular grid.
 *
 * Cells of grid will always be square with initial resolution rangeX/(Nx-1).
 * The Y range may not be exactly rangeY, but will be the nearest value that
 * is a whole number of cells.
 *
 * Note: Deflection solver must be specified before creating a GridMap.
 */
GridMap::GridMap(
                 LensHndl lens       /// lens model for initializing grid
                 ,unsigned long Nx   /// Initial number of grid points in X dimension.
                 ,const PosType my_center[2]  /// Center of grid.
                 ,PosType rangeX   /// Full width of grid in x direction in whatever units will be used.
                 ,PosType rangeY  /// Full width of grid in y direction in whatever units will be used.
): Ngrid_init(Nx),axisratio(rangeY/rangeX),x_range(rangeX){
  
  pointID = 0;
  
  center[0] = my_center[0];
  center[1] = my_center[1];

  assert(Nx > 0);
  assert(rangeX > 0 && rangeY >0);
  
  if(Nx <= 0){ERROR_MESSAGE();
    std::cout << "cannot make GridMap with no points" << std::endl;
    throw std::runtime_error("");
  }
  if(rangeX <= 0 || rangeY <= 0 ){ERROR_MESSAGE();
    throw std::runtime_error("");
    std::cout << "cannot make GridMap with no range" << std::endl;
  }
  
  //if(Ngrid_init % 2 == 1 ) ++Ngrid_init;
  
  Ngrid_init2 = (int)(Ngrid_init*axisratio);
  //if(Ngrid_init2 % 2 == 1) ++Ngrid_init2;
  
  i_points = NewPointArray(Ngrid_init*Ngrid_init2);
 
  int i;
  // set grid of positions
  for(int ii=0;ii<Ngrid_init;++ii){
    for(int jj=0;jj<Ngrid_init2;++jj){
      
      i = ii + jj*Ngrid_init;
      
      i_points[i].id=pointID;
      ++pointID;
      
      i_points[i].x[0] = center[0] + rangeX*ii/(Ngrid_init-1) - 0.5*rangeX;
      i_points[i].x[1] = center[1] + rangeX*jj/(Ngrid_init-1) - 0.5*rangeY;
      i_points[i].gridsize=rangeX/(Ngrid_init-1);
    }
  }
  
  assert(i == Ngrid_init*Ngrid_init2-1);
  
  s_points = NewPointArray(Ngrid_init*Ngrid_init2);
  LinkToSourcePoints(i_points.data(),s_points.data(),Ngrid_init*Ngrid_init2);
  {
    std::lock_guard<std::mutex> hold(grid_mutex);
    lens->rayshooterInternal(Ngrid_init*Ngrid_init2,i_points.data());
  }
}

/** 
 * \brief Constructor for initializing square grid.
 *
 * Note: Deflection solver must be specified before creating a GridMap.
 */
GridMap::GridMap(
                 LensHndl lens               /// lens model for initializing grid
                 ,unsigned long N1d          /// Initial number of grid points in each dimension.
                 ,const PosType my_center[2]    /// Center of grid.
                 ,PosType range              /// Full width of grid in whatever units will be used.
): Ngrid_init(N1d),Ngrid_init2(N1d),axisratio(1.0),x_range(range){
  
  pointID = 0;
  
  if(N1d < 1) throw std::runtime_error("GridMap size is < 1!");
  if(range <= 0) throw std::runtime_error("GridMap range is <= 0!");
  
  center[0] = my_center[0];
  center[1] = my_center[1];
  
  if(N1d <= 0){ERROR_MESSAGE(); std::cout << "cannot make GridMap with no points" << std::endl; exit(1);}
  if(range <= 0){ERROR_MESSAGE(); std::cout << "cannot make GridMap with no range" << std::endl; exit(1);}
  
  i_points = NewPointArray(Ngrid_init*Ngrid_init);
  xygridpoints(i_points.data(),range,center.x,Ngrid_init,0);
  s_points = NewPointArray(Ngrid_init*Ngrid_init);
  LinkToSourcePoints(i_points.data(),s_points.data(),Ngrid_init*Ngrid_init);
    
  {
    std::lock_guard<std::mutex> hold(grid_mutex);
    lens->rayshooterInternal(Ngrid_init*Ngrid_init,i_points.data());
  }
}

GridMap::GridMap(
                 unsigned long N1d          /// Initial number of grid points in each dimension.
                 ,const PosType my_center[2]    /// Center of grid.
                 ,PosType range              /// Full width of grid in whatever units will be used.
): Ngrid_init(N1d),Ngrid_init2(N1d),axisratio(1.0),x_range(range){
  
  pointID = 0;
  
  if(N1d < 1) throw std::runtime_error("GridMap size is < 1!");
  if(range <= 0) throw std::runtime_error("GridMap range is <= 0!");
  
  center[0] = my_center[0];
  center[1] = my_center[1];
  
  if(N1d <= 0){ERROR_MESSAGE(); std::cout << "cannot make GridMap with no points" << std::endl; exit(1);}
  if(range <= 0){ERROR_MESSAGE(); std::cout << "cannot make GridMap with no range" << std::endl; exit(1);}
  
  i_points = NewPointArray(Ngrid_init*Ngrid_init);
  xygridpoints(i_points.data(),range,center.x,Ngrid_init,0);
  s_points = NewPointArray(Ngrid_init*Ngrid_init);
  xygridpoints(s_points.data(),range,center.x,Ngrid_init,0);
  LinkToSourcePoints(i_points.data(),s_points.data(),Ngrid_init*Ngrid_init);
}

GridMap::~GridMap(){}

GridMap GridMap::ReInitialize(LensHndl lens){
  
  GridMap newgrid(lens,Ngrid_init,center.x,x_range,getYRange());
  
  return newgrid;
  /*
  {
    std::lock_guard<std::mutex> hold(grid_mutex);
    lens->rayshooterInternal(Ngrid_init*Ngrid_init,i_points);
  }
  ClearSurfaceBrightnesses();
   */
}

void GridMap::deLens(){
  size_t n = Ngrid_init*Ngrid_init2;
  for(size_t i = 0 ; i < n ; ++i){
    i_points[i].A.setToI();
    s_points[i].A.setToI();
    i_points[i].dt = 0;
    s_points[i].dt = 0;
    s_points[i].x[0] = i_points[i].x[0];
    s_points[i].x[1] = i_points[i].x[1];
  }
}

double GridMap::RefreshSurfaceBrightnesses(Source* source){
  PosType total=0,tmp;
  
  double res2 = pow(getResolution(),2);
  for(size_t i=0;i <s_points[0].head;++i){
    tmp = source->SurfaceBrightness(s_points[i].x);
    s_points[i].surface_brightness = s_points[i].image->surface_brightness
    = tmp;
    total += tmp * res2;
    s_points[i].in_image = s_points[i].image->in_image = NO;
  }
  
  return total;
}

double GridMap::AdaptiveRefreshSurfaceBrightnesses(Lens &lens,Source &source){
  PosType f=1.0e-4;
  
  double resolution = getResolution();
  LinkedPoint point;

  double total_flux = RefreshSurfaceBrightnesses(&source)/resolution/resolution;
  double original_total_flux = total_flux; // ????
  if(total_flux == 0) return 0;
  
  for(long j=0 ; j < Ngrid_init2 ; ++j){
    size_t k = j*Ngrid_init;
    for(long i=0 ; i < Ngrid_init ; ++i,++k){

      if(to_refine(i,j,total_flux,f)){
 
        Point_2d ll = *(i_points[k].image);
        ll[0] -= resolution/2;
        ll[1] -= resolution/2;

        double new_flux = i_points[k].surface_brightness;
        double old_flux=0;
        int n = 1;

        total_flux -= new_flux;
        while ( fabs(old_flux-new_flux) > 1.0e-3*new_flux ){
          old_flux = new_flux;
          n *= 3;
          
          double res = resolution/n;
          new_flux = 0;

          for(int u=0 ; u<n ; ++u){
            point.x[0] = ll[0] + (0.5 + u)*res;
            for(int v=0 ; v<n ; ++v){
               point.x[1] = ll[1] + (0.5 + v)*res;

               lens.rayshooterInternal(1,&point);
               new_flux += source.SurfaceBrightness(point.image->x);
            }
          }
          
          new_flux /= n*n;
        }
        
        //std::cout << " change in pixel flux = " << (new_flux-original) / original << std::endl;
        i_points[k].surface_brightness = new_flux;
        total_flux += new_flux;
      }
      
    }
  }
  
  //std::cout << " change in total flux = " << (total_flux-original_total_flux) / original_total_flux << std::endl;

  assert( (total_flux-original_total_flux) <  0.1*original_total_flux);
  return total_flux * resolution * resolution;
}

bool GridMap::to_refine(long i,long j,double total,double f) const {
  double flux = i_points[i + j*Ngrid_init].surface_brightness;
  
  for(size_t ii=MAX<size_t>(i-1,0) ; ii < MIN<size_t>(i+2,Ngrid_init) ;++ii){
    for(size_t jj=MAX<size_t>(j-1,0) ; jj < MIN<size_t>(j+2,Ngrid_init2) ;++jj){
      size_t kk = ii + jj*Ngrid_init;
      if( fabs(flux - i_points[kk].surface_brightness ) / total > f ) return true;
    }
  }
  return false;
}

double GridMap::AddSurfaceBrightnesses(Source* source){
  PosType total=0,tmp;
  
  for(size_t i=0;i <s_points[0].head;++i){
    tmp = source->SurfaceBrightness(s_points[i].x);
    s_points[i].surface_brightness += tmp;
    s_points[i].image->surface_brightness += tmp;
    total += tmp;
    s_points[i].in_image = s_points[i].image->in_image = NO;
  }
  
  return total * pow(getResolution(),2);
}

void GridMap::ClearSurfaceBrightnesses(){
  
  for(size_t i=0;i <s_points[0].head;++i){
    s_points[i].surface_brightness = s_points[i].image->surface_brightness
    = 0.0;
    s_points[i].in_image = s_points[i].image->in_image = NO;
  }
}

void GridMap::assertNAN(){
  for(size_t i=0;i <s_points[0].head;++i){
    assert(!isnan(s_points[i].surface_brightness));
  }
  for(size_t i=0;i <i_points[0].head;++i){
     assert(!isnan(i_points[i].surface_brightness));
   }
}


void GridMap::xygridpoints(Point *i_points,PosType range,const PosType *center,long Ngrid_1d,short remove_center){
  /* make a new rectolinear grid of points on the image plane **/
  /* and link them to points on the source plane **/
  /* remove_center = 0 include center point of grid */
  /*              != 1 leave out center point of grid if Ngrid_1d is odd, Ngrid_1d*Ngrid_1d-1 points outputted */
  /* warning: memory for i_points must be allocated before entering */
  long i,j;
  
  if(remove_center && (Ngrid_1d%2 == 1)){
    for(i=0,j=0;i<Ngrid_1d*Ngrid_1d;++i){
      
      if( (2*(i/Ngrid_1d)/(Ngrid_1d-1) == 1) && (i%Ngrid_1d == Ngrid_1d/2+1) ) j=1;
      i_points[i-j].id=pointID;
      ++pointID;
      Utilities::PositionFromIndex(i,i_points[i-j].x,Ngrid_1d,range,center);
      i_points[i-j].gridsize=range/(Ngrid_1d-1);
    }
    
  }else{
    for(i=0;i<Ngrid_1d*Ngrid_1d;++i){
      i_points[i].id=pointID;
      ++pointID;
      Utilities::PositionFromIndex(i,i_points[i].x,Ngrid_1d,range,center);
      i_points[i].gridsize=range/(Ngrid_1d-1);
    }
  }
  
  return;
}

PosType GridMap::EinsteinArea() const{
  size_t count = 0;
  size_t N = Ngrid_init*Ngrid_init2;
  for(size_t i=0;i<N;++i){
    if(i_points[i].invmag() < 0) ++count;
  }
  
  return count*x_range*x_range/Ngrid_init/Ngrid_init;
}
// discontinued because it is unstable when there are very demagnified regiona
////PosType GridMap::magnification() const{
//
////  double mag = 0,flux = 0;
//  size_t N = Ngrid_init*Ngrid_init2;
//  for(size_t i=0;i<N;++i){
//    mag += i_points[i].surface_brightness*fabs(i_points[i].invmag());
//    flux += i_points[i].surface_brightness;
//  }
//  return flux/mag;
//}

PosType GridMap::magnificationFlux(Source &source) const{

  double magnified_flux = 0,unmagnified_flux = 0;
  size_t N = Ngrid_init*Ngrid_init2;
  for(size_t i=0;i<N;++i){
    magnified_flux += source.SurfaceBrightness(s_points[i].x);
    unmagnified_flux += source.SurfaceBrightness(i_points[i].x);
  }
  return magnified_flux / unmagnified_flux ;
}

double GridMap::magnificationTr() const {
  
  size_t k;
  size_t k1;
  size_t k2;
  size_t k3;
  
  double sb;
  double flux_source=0.0,flux_image=0.0;
  for(size_t i=0 ; i< Ngrid_init-1 ; i++){
    for(size_t j=0 ; j< Ngrid_init2-1 ; j++){
      
      k = i + j*Ngrid_init;
      assert( s_points[k].surface_brightness == i_points[k].surface_brightness);
      k1 = k + Ngrid_init + 1;

      Point_2d cross = s_points[k1] - s_points[k];

      k2 = k + 1;
      k3 = k + Ngrid_init;

      // 9 X surface brightness at centroid of lower triangle
      sb = ( 2*(s_points[k].surface_brightness + s_points[k1].surface_brightness)
            + 4*s_points[k2].surface_brightness  + s_points[k3].surface_brightness );

      
      if(sb !=0.0){
        flux_source += sb * fabs( cross^(s_points[k2] - s_points[k]) );
      }
      
      flux_image += sb;
      
      // 9 X surface brightness at centroid of upper triangle
       sb = ( 2*(s_points[k].surface_brightness + s_points[k1].surface_brightness)
            + 4*s_points[k3].surface_brightness  + s_points[k2].surface_brightness );
      
      if(sb !=0.0){
        flux_source += sb * fabs( cross^(s_points[k3] - s_points[k]) );
      }
      
      flux_image += sb;
      
    }
  }
  
  return flux_image / flux_source * getResolution() * getResolution();
}

double GridMap::magnificationTr(std::vector<size_t> &ps) const {
  
  size_t k1;
  size_t k2;
  size_t k3;
  
  long nx = Ngrid_init-1;
  long ny = Ngrid_init2-1;
 
  double sb;
  double flux_source=0.0,flux_image=0.0;
  for(size_t k : ps){
      
    if( k % Ngrid_init < nx && k / Ngrid_init < ny ){
      assert( s_points[k].surface_brightness == i_points[k].surface_brightness);
      k1 = k + Ngrid_init + 1;

      Point_2d cross = s_points[k1] - s_points[k];

      k2 = k + 1;
      k3 = k + Ngrid_init;

      // 9 X surface brightness at centroid of lower triangle
      sb = ( 2*(s_points[k].surface_brightness + s_points[k1].surface_brightness)
            + 4*s_points[k2].surface_brightness  + s_points[k3].surface_brightness );

      
      if(sb !=0.0){
        flux_source += sb * fabs( cross^(s_points[k2] - s_points[k]) );
      }
      
      flux_image += sb;
      
      // 9 X surface brightness at centroid of upper triangle
       sb = ( 2*(s_points[k].surface_brightness + s_points[k1].surface_brightness)
            + 4*s_points[k3].surface_brightness  + s_points[k2].surface_brightness );
      
      if(sb !=0.0){
        flux_source += sb * fabs( cross^(s_points[k3] - s_points[k]) );
      }
      
      flux_image += sb;
      
    }
  }
  
  return flux_image / flux_source * getResolution() * getResolution();
}

double GridMap::AreaCellOnSourcePlane(size_t k) const{
  
  assert(k % Ngrid_init != Ngrid_init-1);
  assert(k / Ngrid_init != Ngrid_init2-1);
  long k1 = k + Ngrid_init + 1;

  Point_2d cross = s_points[k1] - s_points[k];

  long k2 = k + 1;
  long k3 = k + Ngrid_init;
  
  return ( fabs( cross^(s_points[k2] - s_points[k]) )/2 + fabs( cross^(s_points[k3] - s_points[k]) )/2 );
}

double GridMap::AddPointSource(const Point_2d &y,double flux){
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

void GridMap::find_magnification_contour(
  std::vector<std::vector<Point_2d> > &curves
  ,std::vector<bool> &hits_boundary
  ,double invmag
  ){
  curves.resize(0);
  
  size_t N = Ngrid_init * Ngrid_init2;
  std::vector<bool> bitmap(N);
  size_t count = 0;
  // find tangential critical curves
  for(size_t k = 0 ; k < N ; ++k){
    if(i_points[k].invmag() < invmag){
      bitmap[k] = true;
      ++count;
    } else {
      bitmap[k] = false;
    }
  }
  std::vector<std::vector<long> > indexes;
  if(count>0){
    Utilities::find_boundaries<Point_2d>(bitmap,Ngrid_init,curves,hits_boundary,false);
  }
  // rescale from pixel units to those of grid
  double resolution = getResolution();
  for(int i=0; i<curves.size() ; ++i){
    for(Point_2d &p : curves[i]) p = p * resolution + i_points[0];
  }
}




void GridMap::find_crit(std::vector<std::vector<Point_2d> > &curves
               ,std::vector<bool> &hits_boundary
               ,std::vector<CritType> &crit_type                      
                        ){
  
  curves.resize(0);
  
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
  //if(count>0) find_boundaries(bitmap,points,hits_boundary);
  std::vector<std::vector<long> > indexes;
  if(count>0){
    Utilities::find_boundaries<Point_2d>(bitmap,Ngrid_init,curves,hits_boundary,false);
    
    Utilities::find_islands(bitmap,Ngrid_init,indexes,hits_boundary);
  }
  
  crit_type.resize(curves.size());
  for(CritType &b : crit_type) b = CritType::tangential;
  long ntange = curves.size();  // number of tangential curves
  
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
  
  if(count>0){
    Utilities::find_boundaries<Point_2d>(bitmap,Ngrid_init,curves,hits_boundary,true);
  }
  
  long m=curves.size();
  crit_type.resize(m);
  for(long i=ntange ; i<m ; ++i) crit_type[i] = CritType::radial;
  
  // reorder them so that radial curves follow the tangential curves they are within
  for(long j=ntange ; j<m ; ++j){
    // pixel in radial critical curve
    for(int i=0 ; i<j ;++i){
      if(crit_type[i] == CritType::tangential 
        && Utilities::inCurve(curves[j][0] ,curves[i] )
         ){
        for(long k=j ; k>i+1 ; --k){
          std::swap(curves[k],curves[k-1]);
          std::swap(crit_type[k],crit_type[k-1]);
          {
            bool tmp = hits_boundary[k];
            hits_boundary[k] = hits_boundary[k-1];
            hits_boundary[k-1] = tmp;
          }
          //std::swap(hits_boundary[k],hits_boundary[k-1]);
        }
        
        break;
      }
    }
  }
  
  // rescale from pixel units to those of grid
  double resolution = getResolution();
  for(int i=0; i<curves.size() ; ++i){
    for(Point_2d &p : curves[i]) p = p * resolution + i_points[0];
  }
  
  // if radial caustic has not been found, estimate a pseudo caustic as the convex hull of negative magnification images
  int ii_tan=0;
  for(int j=0 ; j<curves.size() ; ++j){
    if(crit_type[j] == CritType::tangential){
      if(j==curves.size()-1 || crit_type[j+1] == CritType::tangential ){ // has no radial critical curve
        
        std::vector<Point_2d> psudo(indexes[ii_tan].size());
        
        for(size_t i=0 ; i<indexes[ii_tan].size() ; ++i) psudo[i] = s_points[ indexes[ii_tan][i] ];
        
        std::vector<size_t> hull_index;
        Utilities::convex_hull(psudo,hull_index);
        
        std::vector<Point_2d> v(hull_index.size());
        curves.push_back(v);
        
        for(long i=0 ; i< hull_index.size() ; ++i) curves.back()[i] = i_points[ indexes[ii_tan][ hull_index[i] ]  ];
        
        //write_csv("test_pseud.csv", curves.back());
        
        crit_type.push_back(CritType::pseudo);
        hits_boundary.push_back(false);
        
        for(int k=curves.size()-1 ; k>j+1 ; --k){
          std::swap(curves[k],curves[k-1]);
          std::swap(crit_type[k],crit_type[k-1]);
          {
            bool tmp = hits_boundary[k];
            hits_boundary[k] = hits_boundary[k-1];
            hits_boundary[k-1] = tmp;
          }
        }
      } // radial missing
      ++ii_tan;
    } // is a tangent
  } // loop through all
  assert(2*ntange <= curves.size());
}

//PosType GridMap::magnification2() const{
//  double mag = 0,flux = 0;
//  size_t N = Ngrid_init*Ngrid_init2;
//  for(size_t i=0;i<N;++i){
//    mag += i_points[i].surface_brightness;
//    flux += i_points[i].surface_brightness/fabs(i_points[i].invmag());
//  }
//  return flux/mag;
//}

//void GridMap::find_boundaries(std::vector<bool> &bitmap  // = true inside
//                     ,std::vector<std::vector<Point_2d> > &points
//                     ,std::vector<bool> &hits_edge
//                     ,bool add_to_vector
//                     ){
//
//  size_t nx = Ngrid_init;
//  size_t ny = Ngrid_init2;
//  size_t n = nx*ny;
//
//  std::vector<bool> not_used(bitmap.size(),true);
//
//  assert(bitmap.size()==n);
//
//  // pad edge of field with bitmap=false
//  for(size_t i=0 ; i<nx ; ++i) bitmap[i]=false;
//  size_t j = nx*(ny-1);
//  for(size_t i=0 ; i<nx ; ++i) bitmap[i + j]=false;
//  for(size_t i=0 ; i<ny ; ++i) bitmap[i*nx]=false;
//  j = nx-1;
//  for(size_t i=0 ; i<ny ; ++i) bitmap[j + i*nx]=false;
//
//  std::list<std::list<Point_2d>> contours;
//  if(!add_to_vector){
//    hits_edge.resize(0);
//  }
//
//  bool done = false;
//  long kfirst_in_bound = -1;
//  while(!done){
//    // find first cell in edge
//    size_t k=0;
//    int type;
//    for( k = kfirst_in_bound + 1 ; k < n - nx ; ++k){
//      if(k % nx != nx-1){ // one less cells than points
//        type = 0;
//        if(bitmap[k] ) type +=1;
//        if(bitmap[k+1]) type += 10;
//        if(bitmap[k + nx]) type += 100;
//        if(bitmap[k + nx + 1]) type += 1000;
//
//        if(type > 0
//           && type != 1111
//           && not_used[k]
//           ) break;
//      }
//    }
//
//    kfirst_in_bound = k;
//
//    if(k == n-nx){
//      done=true;
//    }else{ // found an edge
//
//      contours.resize(contours.size() + 1);
//      std::list<Point_2d> &contour = contours.back();
//      hits_edge.push_back(false);
//
//      int type;
//      int face_in=0;
//      size_t n_edge = 0;
//
//      // follow edge until we return to the first point
//      while(k != kfirst_in_bound || n_edge==0){
//
//        if(n_edge >= n){  // infinite loop, output debugging data
//          std::cerr << "Too many points in GridMap::find_boundaries()." << std::endl;
//          std::cerr << "kfirst_in_bound " << kfirst_in_bound << std::endl;
//          std::cerr << "  countour is output to boundary_error_file.csv and bitmap_error_file.csv" << std::endl;
//          {
//            std::ofstream file("bitmap_error_file.csv");
//            file << "in,x,y" << std::endl;
//            for(size_t i=0 ; i<n ; ++i){
//              file << bitmap[i] << "," << i_points[i][0] << "," << i_points[i][1] << std::endl;
//            }
//          }
//
//          {
//            std::ofstream file("boundary_error_file.csv");
//            file << "contour,x,y" << std::endl;
//            int i = 0;
//            for(auto &v : contours){
//              for(Point_2d &p : v){
//                file << i << "," << p[0] << "," << p[1] << std::endl;
//              }
//              ++i;
//            }
//          }
//          throw std::runtime_error("caught in loop.");
//        }
//
//        if(k%nx == 0 || k%nx == nx-2) hits_edge.back() = true;
//        if(k/nx == 0 || k/nx == ny-2) hits_edge.back() = true;
//
//        not_used[k] = false;
//
//        ++n_edge;
//        type = 0;
//        // find type of cell
//        if(bitmap[k] ) type +=1;
//        if(bitmap[k+1]) type += 10;
//        if(bitmap[k + nx]) type += 100;
//        if(bitmap[k + nx + 1]) type += 1000;
//
//        if(type == 0 || type == 1111){  // all in or all out
//          throw std::runtime_error("off edge!!");
//        }else if(type == 1 || type == 1110){ // lower left only
//
//          if(face_in==0){
//            contour.push_back( (i_points[k] + i_points[k+1]) / 2 );
//            face_in=1;
//            k -= nx;
//          }else{
//            contour.push_back( (i_points[k] + i_points[k+nx]) / 2 );
//            face_in=2;
//            k -= 1;
//          }
//
//        }else if(type == 10 || type == 1101){ // lower right only
//
//          if(face_in==2){
//            contour.push_back( (i_points[k] + i_points[k+1]) / 2 );
//            face_in=1;
//            k -= nx;
//          }else{
//            contour.push_back( (i_points[k+nx+1] + i_points[k+1]) / 2 );
//            face_in=0;
//            k += 1;
//          }
//
//        }else if(type == 100 || type == 1011){ // upper left only
//
//          if(face_in==0){
//            contour.push_back( (i_points[k+nx] + i_points[k+nx+1]) / 2 );
//            face_in=3;
//            k += nx;
//          }else{
//            contour.push_back( (i_points[k] + i_points[k+nx]) / 2 );
//            face_in=2;
//            k -= 1;
//          }
//
//        }else if(type == 1000 || type == 111){ // upper right only
//
//          if(face_in==1){
//            contour.push_back( (i_points[k+1] + i_points[k+nx+1]) / 2 );
//            face_in=0;
//            k += 1;
//          }else{
//            contour.push_back( (i_points[k+nx] + i_points[k+nx+1]) / 2 );
//            face_in=3;
//            k += nx;
//          }
//
//        }else if(type == 11 || type == 1100){ // lower two
//
//          if(face_in==0){
//            contour.push_back( (i_points[k+1] + i_points[k+nx+1]) / 2 );
//            k += 1;
//          }else{
//            contour.push_back( (i_points[k] + i_points[k+nx]) / 2 );
//            face_in = 2;
//            k -= 1;
//          }
//
//        }else if(type == 1010 || type == 101){ // right two
//
//          if(face_in==1){
//            contour.push_back( (i_points[k] + i_points[k+1]) / 2 );
//            k -= nx;
//          }else{
//            contour.push_back( (i_points[k+nx] + i_points[k+nx+1]) / 2 );
//            face_in = 3;
//            k += nx;
//          }
//
//        }else if(type == 1001){ // lower left upper right
//
//          if(face_in==0){
//            contour.push_back( (i_points[k+nx] + i_points[k+nx+1]) / 2 );
//            face_in=3;
//            k += nx;
//          }else if(face_in==1){
//            contour.push_back( (i_points[k] + i_points[k+nx]) / 2 );
//            face_in=2;
//            k -= 1;
//          }else if(face_in==2){
//            contour.push_back( (i_points[k] + i_points[k+1]) / 2 );
//            face_in=1;
//            k -= nx;
//          }else{
//            contour.push_back( (i_points[k+nx+1] + i_points[k+1]) / 2 );
//            face_in=0;
//            k += 1;
//          }
//
//        }else if(type == 110){ // upper left lower right
//
//          if(face_in==0){
//            contour.push_back( (i_points[k] + i_points[k+1]) / 2 );
//            face_in=1;
//            k -= nx;
//          }else if(face_in==1){
//            contour.push_back( (i_points[k+1] + i_points[k+nx+1]) / 2 );
//            face_in=0;
//            k += 1;
//          }else if(face_in==2){
//            contour.push_back( (i_points[k + nx] + i_points[k+nx+1]) / 2 );
//            face_in=3;
//            k += nx;
//          }else{
//            contour.push_back( (i_points[k] + i_points[k+nx]) / 2 );
//            face_in=2;
//            k -= 1;
//          }
//        }
//      }
//
//      // the diamand with a hole
//      //if(n_edge == 12 && bitmap[k + nx + 1] ){
//      //  n_edge=0;
//      //  face_in = 2;
//      //}
//
//    }
//  }
//
//  int offset = 0;
//  if(!add_to_vector){
//    points.resize(contours.size());
//  }else{
//    offset = points.size();
//    points.resize(points.size() + contours.size());
//  }
//
//  // copy lists of points into vectors
//  int i=0;
//  for(auto &c: contours){
//    points[offset+i].resize(c.size());
//    size_t j=0;
//    for(auto &p: c){
//      points[offset+i][j] = p;
//      ++j;
//    }
//    ++i;
//  }
//}


//void GridMap::find_crit_boundary(std::vector<std::vector<Point_2d> > &points
//               ,std::vector<bool> &hits_boundary
//               ) const{
//  
//}


Point_2d GridMap::centroid() const{
  double flux = 0;
  Point_2d centroid(0,0);
  
  size_t N = Ngrid_init*Ngrid_init2;
  for(size_t i=0;i<N;++i){
    centroid += i_points[i]*i_points[i].surface_brightness;
    flux += i_points[i].surface_brightness;
  }
  return centroid/flux;
}

std::list<RAY> GridMap::find_images(std::vector<Point_2d> &ys
                                    ,std::vector<int> &multiplicity
                                    ) const{
  int nthreads = Utilities::GetNThreads();

  long N =  ys.size();
  long n = (int)(N/nthreads + 1);
  
  multiplicity.resize(N);
  std::vector<std::list<RAY>> v_of_lists(nthreads);
  
  std::vector<std::thread> thr;
  long m = 0;
  for(int i=0 ; i<nthreads ; ++i ){

    if(m > N-n ) n = N-m;
    thr.push_back(std::thread(
                              &GridMap::_find_images_
                              ,this
                              ,ys.data() + m
                              ,multiplicity.data() + m
                              ,n
                              ,std::ref(v_of_lists[i])
                              )
                  );
    
    m += n;
  }

  for(auto &t : thr) t.join();

  //join ray lists
  for(int i = 1 ; i<nthreads ; ++i)
    v_of_lists[0].splice(v_of_lists[0].end(),v_of_lists[i]);
  
  return v_of_lists[0];
}

void GridMap::find_images(Point_2d y
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
          }else if(sig2 == 0){
            image_points.push_back( ( i_points[k1] + i_points[k2])/2  );
          }else{
            image_points.push_back( ( i_points[k2] + i_points[k])/2  );
          }
          triangles.push_back(Triangle(k,k1,k2));
        }
      }else if (sig1 == 0 && sig3 == 0){ // a vertex
        image_points.push_back( i_points[k] );
        triangles.push_back(Triangle(k,k1,k2));
      }
     
      //sig1 = sign( (y-s_points[k])^(s_points[k1]-s_points[k]) );
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
          }else if(sig2 == 0){
            image_points.push_back( ( i_points[k1] + i_points[k3] )/2  );
          }else{
            image_points.push_back( ( i_points[k3] + i_points[k] )/2  );
          }
          triangles.push_back(Triangle(k,k1,k3));
        }
      }
    }
  }
  
  return;
}

void GridMap::limited_image_search(Point_2d &y
                 ,std::vector<size_t> &cell_numbers  /// positions of the images limited by resolution of the gridmap
                 ,std::vector<Triangle> &triangles     /// index's of the points that form the triangles that the images are in
) const {
  
  triangles.clear();
  
  size_t k1;
  size_t k2;
  size_t k3;
  
  int sig1,sig2,sig3,sig_sum;
  
  for(size_t k : cell_numbers){
    k1 = k + Ngrid_init + 1;
    
    k2 = k + 1;
    k3 = k + Ngrid_init;
    
    sig1 = sign( (y-s_points[k])^(s_points[k1]-s_points[k]) );
    sig2 = sign( (y-s_points[k1])^(s_points[k2]-s_points[k1]) );
    sig3 = sign( (y-s_points[k2])^(s_points[k]-s_points[k2]) );
    
    sig_sum = sig1 + sig2 + sig3;
    if(abs(sig_sum) == 3){ // inside
      triangles.push_back(Triangle(k,k1,k2));
    }else if(abs(sig_sum) == 2){ // on edge
      if(sig_sum > 0){
        triangles.push_back(Triangle(k,k1,k2));
      }
    }else if (sig1 == 0 && sig3 == 0){ // a vertex
      triangles.push_back(Triangle(k,k1,k2));
    }
    
    //sig1 = sign( (y-s_points[k])^(s_points[k1]-s_points[k]) );
    sig2 = sign( (y-s_points[k1])^(s_points[k3]-s_points[k1]) );
    sig3 = sign( (y-s_points[k3])^(s_points[k]-s_points[k3]) );
    
    sig_sum = sig1 + sig2 + sig3;
    if(abs(sig_sum) == 3){ // inside
      triangles.push_back(Triangle(k,k1,k3));
    }else if(abs(sig_sum) == 2){ // on edge
      if(sig_sum > 0){
        triangles.push_back(Triangle(k,k1,k3));
      }
    }
  }
  
  return;
}

void GridMap::_find_images_(Point_2d *ys,int *multiplicity,long Nys,std::list<RAY> &rays) const{
  
  std::vector<Point_2d> x;
  std::vector<Triangle> triangles;
  rays.resize(0);
  //auto itr = rays.begin();
  for(long k=0 ; k<Nys ; ++k){
    find_images(ys[k],x,triangles);
    multiplicity[k] = x.size();
    for(int i=0 ; i < x.size() ; ++i){
      rays.emplace_back();
      rays.back().x = x[i];
      rays.back().y = ys[k];
      rays.back().A *= 0;
      for(int j=0 ; j<3 ; ++j){
        rays.back().A = rays.back().A  + i_points[triangles[i][j]].A;
      }
      rays.back().A /= 3;
    }
  }
  
  return;
}

bool  GridMap::incurve(long k,std::vector<Point_2d> &curve) const{
  int n=0;
  long i = k % Ngrid_init , j = k / Ngrid_init;
  for(Point_2d &p : curve){
     if( long(p[0]+0.5) > i && (j - long(p[1]+0.5) ) == 0 ) ++n;
  }
  
  return n%2 == 1;
}

