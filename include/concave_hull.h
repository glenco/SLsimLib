//
//  concave_hull.h
//  GLAMER
//
//  Created by bmetcalf on 23/02/16.
//
//

#ifndef concave_hull_h
#define concave_hull_h

#include <set>
#include "point.h"
#include "geometry.h"
#include "image_processing.h"
#include "simpleTreeVec.h"

namespace Utilities{

template<typename T,typename P>
Point_2d subtract(T& p1,P& p2){
  return Point_2d(p1[0] - p2[0], p1[1] - p2[1]);
}
template<typename P>
double crossD(P &O,P &A,P &B){
  return (A[0] - O[0]) * (B[1] - O[1]) - (A[1] - O[1]) * (B[0] - O[0]);
}


/// removes the intersections of the curve
template <typename T>
size_t RemoveIntersections(std::vector<T> &curve){
  
  if(curve.size() <=3) return 0;
  
  size_t N = curve.size(),count = 0;
  T tmp;
  
  curve.push_back(curve[0]);
  
  for(size_t i=0;i<N-2;++i){
    for(size_t j=i+2;j<N;++j){
      if(Utilities::Geometry::intersect(curve[i].x,curve[i+1].x,curve[j].x,curve[j+1].x)){
        
        size_t k=i+1,l=j;
        while(k < l){
          std::swap(curve[k] , curve[l]);
          ++k;
          --l;
        }
        ++count;
      }
    }
  }
  
  //assert(curve[0]==curve.back());
  curve.pop_back();
  
  return count;
}

/// removes the intersections while removing interior loops
///  The input curve needs to be ordered already.  No points in the
///  input curve will be outside the output hull.
///  Will fail if there are overlapping segments on the hull.
std::vector<Point_2d> TighterHull(const std::vector<Point_2d> &v);

//template <typename T>
//std::vector<T> TightHull(const std::vector<T> &curve){
//
//  if(curve.size() <=3) return curve;
//
//  size_t N = curve.size(),count = 0;
//  T tmp;
//
//  // find left most point
//  long init = 0;
//  double pmin = curve[0][0];
//  for(int i=1 ; i<N ; ++i){
//    if(curve[i][0] < pmin){
//      init = i;
//      pmin = curve[i][0];
//    }
//  }
//
//  // orientation == 1 clockwise
//  Geometry::CYCLIC cyc(N);
//  int orientation = sign( (curve[cyc[init-1]] - curve[init])^(curve[cyc[init+1]] - curve[init])  );
//  if(orientation == 0){
//    orientation = sign( curve[cyc[init+1]][1] - curve[init][1]);
//  }
//
//  orientation *= -1; // orientation == 1 is counter clockwise
//
//  // make a copy where first one is the leftmost
//  std::vector<T> hull(N);
//  for(size_t i=0;i<N;++i){
//    hull[i] = curve[ cyc[ orientation * i + init] ];
//  }
//
//  for(size_t i=0;i<N-2;++i){
//
//    for(size_t j=i+2;j<N;++j){
//      if(Utilities::Geometry::intersect(hull[i].x,hull[ (i+1) % N ].x,hull[j].x,hull[ (j+1) % N ].x)){
//
//
//        long ii = cyc[i-1];
//        //if(i==0) ii = N-1;
//        //else ii = i-1;
//
//        Point_2d x =  hull[i] - hull[ii];
//        x.unitize();
//
//        if(
//           //( x^(hull[ (j+1) % N ] - hull[i]) ) / (hull[ (j+1) % N ] - hull[i]).length()
//           atan2( x^(hull[ (j+1) % N ] - hull[i]) , x*(hull[ (j+1) % N ] - hull[i]) )
//           <
//           atan2( x^(hull[j] - hull[i]) , x*(hull[j] - hull[i]) )
//           //( x^(hull[j] - hull[i]) ) / (hull[j] - hull[i]).length()
//           ){
//
//          // inner loop - remove all points between i and j+1
//
//          int n = j - i;
//          int k = j+1;
//          while(k < N){
//            //std::swap(hull[k],hull[k-n]);
//            hull[k-n] = hull[k];
//            ++k;
//          }
//          N -= n;
//        }else{
//
//          // outer loop - put them in reverse order
//          size_t k=i+1,l=j;
//          while(k < l){
//            std::swap(hull[k] , hull[l]);
//            ++k;
//            --l;
//          }
//        }
//        ++count;
//      }
//    }
//  }
//
//  //assert(hull[0]==hull.back());
//  hull.resize(N);
//
//  assert(N>2);
//  return hull;
//}

/// returns random point within a trinagle
Point_2d RandomInTriangle(const Point_2d &x1,
                          const Point_2d &x2,
                          const Point_2d &x3,
                          Utilities::RandomNumbers_NR &ran);

/// return a point within a convex polygon
Point_2d RandomInConvexPoly(const std::vector<Point_2d> &pp,
                            Utilities::RandomNumbers_NR &ran);

/// return a point within a polygon that doesn't need to be convex
Point_2d RandomInPoly(std::vector<Point_2d> &pp,
                      Utilities::RandomNumbers_NR &ran);

/** \brief finds ordered boundaries to regions where bitmap == true

 This can be used to find critical curves or contours.
 `bitmap` should be the same size as the `Gridmap`
 If the boundary curve  touches the edge of the `GridMap` it will be indicated in `hits_boundary` as
 `true`.
 
 Boundaries will never cross or lead off the grid.  On the edges they will leave the edge pixels out even if they should be in.  This is a technical comprimise.
 
 Output points are in pixel units with (0,0) being pioint (0,0)
*/
template <typename P>
void find_boundaries(std::vector<bool> &bitmap  // = true inside
                     ,long nx  // number of pixels in x direction
                     ,std::vector<std::vector<P> > &points
                     ,std::vector<bool> &hits_edge
                     ,bool add_to_vector=false
                     ){
  
  size_t n = bitmap.size();
  long ny = n/nx;
  
  if(n != nx*ny){
    std::cerr << "Wrong sizes in Utilities::find_boundaries." << std::endl;
    throw std::invalid_argument("invalid size");
  }
  
  std::vector<bool> not_used(n,true);
  
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
  }
  
  bool done = false;
  long kfirst_in_bound = -1;
  while(!done){
    // find first cell in edge
    size_t k=0;
    int type;
    for( k = kfirst_in_bound + 1 ; k < n - nx ; ++k){
      if(k % nx != nx-1){ // one less cells than points
        type = 0;
        if(bitmap[k] ) type +=1;
        if(bitmap[k+1]) type += 10;
        if(bitmap[k + nx]) type += 100;
        if(bitmap[k + nx + 1]) type += 1000;

        if(type > 0
           && type != 1111
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
      
      
      int type;
      int face_in=0;
      size_t n_edge = 0;
      
      // follow edge until we return to the first point
      while(k != kfirst_in_bound || n_edge==0){
        
        if(n_edge >= n){  // infinite loop, output debugging data
          std::cerr << "Too many points in Utilities::find_boundaries()." << std::endl;
          std::cerr << "kfirst_in_bound " << kfirst_in_bound << std::endl;
          std::cerr << "  countour is output to boundary_error_file.csv and bitmap_error_file.csv" << std::endl;
          {
            std::ofstream file("bitmap_error_file.csv");
            file << "in,x,y" << std::endl;
            for(size_t i=0 ; i<n ; ++i){
              file << bitmap[i] << "," << i%nx << "," << i/nx << std::endl;
            }
          }
          
          {
            std::ofstream file("boundary_error_file.csv");
            file << "contour,x,y" << std::endl;
            int i = 0;
            for(auto &v : contours){
              for(Point_2d &p : v){
                file << i << "," << p[0] << "," << p[1] << std::endl;
              }
              ++i;
            }
          }
          throw std::runtime_error("caught in loop.");
        }
        
        if(k%nx == 0 || k%nx == nx-2) hits_edge.back() = true;
        if(k/nx == 0 || k/nx == ny-2) hits_edge.back() = true;
        
        not_used[k] = false;
        
        ++n_edge;
        type = 0;
        // find type of cell
        if(bitmap[k]) type +=1;
        if(bitmap[k+1]) type += 10;
        if(bitmap[k + nx]) type += 100;
        if(bitmap[k + nx + 1]) type += 1000;
        
        if(type == 0 || type == 1111){  // all in or all out
          throw std::runtime_error("off edge!!");
        }else if(type == 1 || type == 1110){ // lower left only
          
          if(face_in==0){
            contour.push_back( Point_2d( k%nx + 0.5
                                       , k/nx ) );
            //contour.push_back( (i_points[k] + i_points[k+1]) / 2 );
            face_in=1;
            k -= nx;
          }else{
            contour.push_back( Point_2d( k%nx
                                       , k/nx + 0.5  ) );
            //contour.push_back( (i_points[k] + i_points[k+nx]) / 2 );
            face_in=2;
            k -= 1;
          }
          
        }else if(type == 10 || type == 1101){ // lower right only
          
          if(face_in==2){
            contour.push_back( Point_2d( k%nx + 0.5
                                       , k/nx  ) );
            //contour.push_back( (i_points[k] + i_points[k+1]) / 2 );
            face_in=1;
            k -= nx;
          }else{
            contour.push_back( Point_2d( k%nx + 1
                                       , k/nx + 0.5) );
            //contour.push_back( (i_points[k+nx+1] + i_points[k+1]) / 2 );
            face_in=0;
            k += 1;
          }
          
        }else if(type == 100 || type == 1011){ // upper left only
          
          if(face_in==0){
            contour.push_back( Point_2d( k%nx + 0.5
                                       , k/nx + 1 ) );
            //contour.push_back( (i_points[k+nx] + i_points[k+nx+1]) / 2 );
            face_in=3;
            k += nx;
          }else{
            contour.push_back( Point_2d( k%nx
                                       , k/nx + 0.5 ) );
            //contour.push_back( (i_points[k] + i_points[k+nx]) / 2 );
            face_in=2;
            k -= 1;
          }
          
        }else if(type == 1000 || type == 111){ // upper right only
          
          if(face_in==1){
            contour.push_back( Point_2d( k%nx + 1
                                       , k/nx + 0.5 ) );
             //contour.push_back( (i_points[k+1] + i_points[k+nx+1]) / 2 );
            face_in=0;
            k += 1;
          }else{
            contour.push_back( Point_2d( k%nx + 0.5
                                       , k/nx + 1 ) );
            //contour.push_back( (i_points[k+nx] + i_points[k+nx+1]) / 2 );
            face_in=3;
            k += nx;
          }
          
        }else if(type == 11 || type == 1100){ // lower two
          
          if(face_in==0){
            contour.push_back( Point_2d( k%nx + 1
                                       , k/nx + 0.5 ) );
            //contour.push_back( (i_points[k+1] + i_points[k+nx+1]) / 2 );
            k += 1;
          }else{
            contour.push_back( Point_2d( k%nx
                                       , k/nx + 0.5 ) );
            //contour.push_back( (i_points[k] + i_points[k+nx]) / 2 );
            face_in = 2;
            k -= 1;
          }
          
        }else if(type == 1010 || type == 101){ // right two
          
          if(face_in==1){
            contour.push_back( Point_2d( k%nx + 0.5
                                       , k/nx ) );
            //contour.push_back( (i_points[k] + i_points[k+1]) / 2 );
            k -= nx;
          }else{
            contour.push_back( Point_2d( k%nx + 0.5
                                       , k/nx + 1 ) );
            //contour.push_back( (i_points[k+nx] + i_points[k+nx+1]) / 2 );
            face_in = 3;
            k += nx;
          }
          
        }else if(type == 1001){ // lower left upper right
          
          if(face_in==0){
            contour.push_back( Point_2d( k%nx + 0.5
                                       , k/nx + 1 ) );
            //contour.push_back( (i_points[k+nx] + i_points[k+nx+1]) / 2 );
            face_in=3;
            k += nx;
          }else if(face_in==1){
            contour.push_back( Point_2d( k%nx
                                       , k/nx + 0.5 ) );
            //contour.push_back( (i_points[k] + i_points[k+nx]) / 2 );
            face_in=2;
            k -= 1;
          }else if(face_in==2){
            contour.push_back( Point_2d( k%nx + 0.5
                                       , k/nx ) );
            //contour.push_back( (i_points[k] + i_points[k+1]) / 2 );
            face_in=1;
            k -= nx;
          }else{
            contour.push_back( Point_2d( k%nx + 1
                                       , k/nx + 0.5 ) );
            //contour.push_back( (i_points[k+nx + 1] + i_points[k+1]) / 2 );
            face_in=0;
            k += 1;
          }
          
        }else if(type == 110){ // upper left lower right
          
          if(face_in==0){
            contour.push_back( Point_2d( k%nx + 0.5
                                       , k/nx  ) );
            //contour.push_back( (i_points[k] + i_points[k+1]) / 2 );
            face_in=1;
            k -= nx;
          }else if(face_in==1){
            contour.push_back( Point_2d( k%nx + 1
                                       , k/nx + 0.5 ) );
            //contour.push_back( (i_points[k+1] + i_points[k+nx+1]) / 2 );
            face_in=0;
            k += 1;
          }else if(face_in==2){
            contour.push_back( Point_2d( k%nx + 0.5
                                       , k/nx + 1 ) );
            //contour.push_back( (i_points[k + nx] + i_points[k+nx+1]) / 2 );
            face_in=3;
            k += nx;
          }else{
            contour.push_back( Point_2d( k%nx
                                       , k/nx + 0.5 ) );
            //contour.push_back( (i_points[k] + i_points[k+nx]) / 2 );
            face_in=2;
            k -= 1;
          }
        }
      }
      
      // special case of the diamand with a hole
      //if(n_edge == 12 && bitmap[k + nx + 1] ){
      //  n_edge=0;
      //  face_in = 2;
      //}
      
    }
  }
  
  long offset = 0;
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

/// find the indexes of areas with bitmap[]=true broken up into diconnected islands
void find_islands(std::vector<bool> &bitmap  // = true inside
                  ,long nx  // number of pixels in x direction
                  ,std::vector<std::vector<long> > &indexes
                  ,std::vector<bool> &hits_edge
                  ,bool add_to_vector=false
                  );

/// this returns area within the curve x average kappa iwithin the curve
double interior_mass(const std::vector<Point_2d> &alpha
                     ,const std::vector<Point_2d> &x);

/// this returns area within the curve x average kappa iwithin the curve
double interior_mass(const std::vector<RAY> &rays);

/// Returns a vector of points on the convex hull in counter-clockwise order.
template<typename T>
std::vector<T> convex_hull(const std::vector<T> &PP)
{
 
  if(PP.size() <= 3){
    return PP;
  }

  std::vector<T> P=PP;

  size_t n = P.size();
  size_t k = 0;
  std::vector<T> hull(2*n);
  
  // Sort points lexicographically
  std::sort(P.begin(), P.end(),
            [](T p1,T p2){
    if(p1[0]==p2[0]) return p1[1] < p2[1];
    return p1[0] < p2[0];});
  
  
  // Build lower hull
  for (size_t i = 0; i < n; i++) {
    while (k >= 2 && crossD(hull[k-2], hull[k-1], P[i]) <= 0){
      k--;
    }
    hull[k++] = P[i];
  }
  
  // Build upper hull
  for (long i = n-2, t = k+1; i >= 0; i--) {
    while (k >= t && crossD(hull[k-2], hull[k-1], P[i]) <= 0){
      k--;
      assert(k > 1);
    }
    hull[k++] = P[i];
  }
  
  
  hull.resize(k);
  hull.pop_back();
  
  return;
}

/// Returns a vector of points on the convex hull in counter-clockwise order.
///
template<typename T>
void convex_hull(const std::vector<T> &P,std::vector<size_t> &hull)
{
  
  size_t n = P.size();
  
  if(n <= 3){
    hull.resize(n);
    size_t i = 0;
    for(size_t &d : hull) d = i++;
    return;
  }
  
  //std::vector<T> hull(n + 1);
  hull.resize(n + 1);
  size_t k = 0;
  
  // Sort points lexicographically
  std::vector<size_t> sort_index(P.size());
  {
    int i=0;
    for(size_t &a :  sort_index) a = i++;
  }
  std::sort(sort_index.begin(), sort_index.end(),
            [&P](size_t i1,size_t i2){
    if(P[i1][0]==P[i2][0]) return P[i1][1] < P[i2][1];
    return P[i1][0] < P[i2][0];});
  
  
  // Build lower hull
  for (size_t i = 0; i < n; i++) {
    while (k >= 2 && crossD(P[hull[k-2]], P[hull[k-1]], P[sort_index[i]]) <= 0){
      k--;
    }
    hull[k++] = sort_index[i];
  }
  
  // Build upper hull
  for (long i = n-2, t = k+1; i >= 0; i--) {
    while (k >= t && crossD(P[hull[k-2]], P[hull[k-1]],P[sort_index[i]]) <= 0){
      k--;
      assert(k > 1);
    }
    hull[k++] = sort_index[i];
  }
    
  hull.resize(k);
  hull.pop_back();
  
  return;
}


struct Edge{
  Edge(size_t i,double l):length(l),index(i){};
  double length = 0;
  size_t index = 0;
};
/** \brief Creates the concave hull of a group of 2 dimensional points
 by the shrink-wrap algorithm.
 
 The type of the input vector points must have an operator [].
 If the input vector is the same as the output vector it will be replaced,
 and the function will still work.
 
 It is guaranteed that the resulting hull will surround the all the points.  Any edge that is greater than scale will be refined until it is either smaller than scale or it cannot be refined further.  As a result some edges might be larger than scale and some smaller.
 
 This should be a NlogN algorithm.
 
 The algorithm:  1) The convex hull is found.  2) The longest edge is found
 3) all the points that are not in the hull are tested to see if they are within the rays extending from the end point perpendicular to the edge.  4) Of the points that are the one that makes the smallest area triangle with the end points is chosen and added 5) go back to 3 if there are edges that are larger than scale and new points exist to be added 6) remove all intersections in the hull
 
 */
template<typename T>
void concave(std::vector<T> &init_points
             ,std::vector<T> &hull_out,double scale)
{
  
  //typedef typename InputIt::value_type point;
  
  bool TEST = false;
  
  // find the convex hull
  std::vector<T> hull = convex_hull(init_points);
  
  if(init_points.size() == hull.size()) return;
  
  std::list<Edge> edges;
  std::list<T> leftovers;
  hull.push_back(hull[0]);
  
  double tmp;
  
  // find pairs of hull points that are further appart than scale
  for(int i=0;i<hull.size()-1;++i){
    tmp = (subtract(hull[i],hull[i+1])).length();
    if(tmp > scale) edges.emplace_back(i,tmp);// .push_back(std::pair<size_t,double>(i,tmp));
  }
  
  if(edges.size() == 0) return;
  
  if(TEST){
    PosType cent[] ={0.5,0.5};
    PixelMap map(cent,256, 1.0/256);
    
    map.drawgrid(10,0.5);
    
    map.drawPoints(init_points,0.01,1.0);
    map.drawCurve(hull,2);
    
    map.printFITS("!test_concave.fits");
  }
  
  // sort edges by length
  edges.sort([](const Edge &p1,const Edge &p2){return p1.length > p2.length;});
  assert(edges.front().length >= edges.back().length);
  
  
  {   // make a list of points that are not already in the convex hull
    
    hull.pop_back();
    std::vector<size_t> index(hull.size());
    size_t i = 0;
    for(size_t &ind: index) ind = i++;
    std::sort(index.begin(),index.end(),
              [&hull](size_t p1,size_t p2){
      if(hull[p1][0] == hull[p2][0]) return hull[p1][1] < hull[p2][1];
      return hull[p1][0] < hull[p2][0];});
    
    size_t j = 0;
    
    for(auto pit = init_points.begin(); pit != init_points.end() ; ++pit){
      
      if((hull[index[j]][0] == (*pit)[0])*(hull[index[j]][1] == (*pit)[1])){
        ++j;
      }else{
        leftovers.push_back(*pit);
      }
      
    }
    assert(hull.size() + leftovers.size() == init_points.size());
    hull.push_back(hull[0]);
  }
  
  typename std::list<T>::iterator p,nextpoint;
  //auto p = leftovers.begin();
  
  double minarea,area,co1;//,co2;
  size_t i;
  Point_2d vo,v1,v2;
  
  //int count = 0;
  while(edges.front().length > scale ){
    
    if(leftovers.size() == 0) break;
    
    // find the point which if added to this edge would change the area least
    minarea = HUGE_VAL;
    i = edges.front().index;
    
    for(p = leftovers.begin() ; p != leftovers.end() ; ++p){
      
      v1 = subtract(*p,hull[i]);
      
      if( edges.front().length > v1.length()){
        
        vo = subtract(hull[i+1], hull[i]);
        //v2 = subtract(*p,hull[i+1]);
        area = vo^v1;
        co1 = v1*vo;
        //co2 = v2*vo;
        
        if( co1 > 0 && (area > 0)*(area < minarea) ){
          //      if(co1 > 0 && area > 0 && index_length.front().second > v1.length()){
          minarea = area;
          nextpoint = p;
        }
      }
    }
    
    if(minarea == HUGE_VAL){
      // if there is no acceptable point for this edge continue with the second longest edge
      edges.pop_front();
      continue;
    }
    
    // insert new point into hull
    T tmp = *nextpoint;
    leftovers.erase(nextpoint);
    hull.insert(hull.begin() + i + 1,tmp);
    
    // update index_length list
    Edge new1 = edges.front();
    edges.pop_front();
    new1.length = (subtract(hull[i],hull[i+1])).length();
    Edge new2(i+1,(subtract(hull[i+1],hull[i+2])).length());
    
    for(auto &p : edges) if(p.index > i) ++(p.index);
    
    auto p = edges.begin();
    if(new1.length > scale){
      for(p = edges.begin() ; p != edges.end() ; ++p){
        if((*p).length < new1.length){
          edges.insert(p,new1);
          break;
        }
      }
      if(p==edges.end()) edges.push_back(new1);
    }
    if(new2.length > scale){
      for(p = edges.begin(); p != edges.end() ; ++p){
        if((*p).length < new2.length){
          edges.insert(p,new2);
          break;
        }
      }
      if(p==edges.end()) edges.push_back(new2);
    }
    
    
    if(TEST){
      PosType cent[] ={0.5,0.5};
      PixelMap map(cent,256, 1.0/256);
      
      map.drawgrid(10,0.5);
      
      map.drawPoints(init_points,0.03,1.0);
      map.drawCurve(hull,2);
      
      map.printFITS("!test_concave.fits");
    }
  }
  
  hull.pop_back();
  
  //Utilities::RemoveIntersections(hull);
  
  std::swap(hull,hull_out);
  
  return;
}

template<typename T>
std::vector<T> concave2(std::vector<T> &init_points,double scale)
{
  
  //typedef typename InputIt::value_type point;
  
  bool TEST = false;
  
  // find the convex hull
  std::vector<T> hull = convex_hull(init_points);
  
  if(init_points.size() == hull.size()) return hull;
  
  std::list<Edge> edges;
  std::list<T> leftovers;
  hull.push_back(hull[0]);
  
  double tmp;
  
  // find pairs of hull points that are further appart than scale
  for(int i=0;i<hull.size()-1;++i){
    tmp = (subtract(hull[i],hull[i+1])).length();
    if(tmp > scale) edges.emplace_back(i,tmp);// .push_back(std::pair<size_t,double>(i,tmp));
  }
  
  if(edges.size() == 0) return hull;
  
  if(TEST){
    PosType cent[] ={0.5,0.5};
    PixelMap map(cent,256, 1.0/256);
    
    map.drawgrid(10,0.5);
    
    map.drawPoints(init_points,0.01,1.0);
    map.drawCurve(hull,2);
    
    map.printFITS("!test_concave.fits");
  }
  
  // sort edges by length
  edges.sort([](const Edge &p1,const Edge &p2){return p1.length > p2.length;});
  assert(edges.front().length >= edges.back().length);
  
  
  {   // make a list of points that are not already in the convex hull
    
    hull.pop_back();
    std::vector<size_t> index(hull.size());
    size_t i = 0;
    for(size_t &ind: index) ind = i++;
    std::sort(index.begin(),index.end(),
              [&hull](size_t p1,size_t p2){
      if(hull[p1][0] == hull[p2][0]) return hull[p1][1] < hull[p2][1];
      return hull[p1][0] < hull[p2][0];});
    
    size_t j = 0;
    
    for(auto pit = init_points.begin(); pit != init_points.end() ; ++pit){
      
      if((hull[index[j]][0] == (*pit)[0])*(hull[index[j]][1] == (*pit)[1])){
        ++j;
      }else{
        leftovers.push_back(*pit);
      }
      
    }
    assert(hull.size() + leftovers.size() == init_points.size());
    hull.push_back(hull[0]);
  }
  
  typename std::list<T>::iterator p,nextpoint;
  //auto p = leftovers.begin();
  
  double minarea,area,co1;
  size_t i;
  Point_2d vo,v1,v2;
  
  while(edges.front().length > scale ){
    
    if(leftovers.size() == 0) break;
    
    // find the point which if added to this edge would change the area least
    minarea = HUGE_VAL;
    i = edges.front().index;
    
    for(p = leftovers.begin() ; p != leftovers.end() ; ++p){
      
      v1 = subtract(*p,hull[i]);
      
      if( edges.front().length > v1.length()){
        
        vo = subtract(hull[i+1], hull[i]);
        //v2 = subtract(*p,hull[i+1]);
        area = vo^v1;
        co1 = v1*vo;
        //co2 = v2*vo;
        
        if( co1 > 0 && (area > 0)*(area < minarea) ){
          //      if(co1 > 0 && area > 0 && index_length.front().second > v1.length()){
          minarea = area;
          nextpoint = p;
        }
      }
    }
    
    if(minarea == HUGE_VAL){
      // if there is no acceptable point for this edge continue with the second longest edge
      edges.pop_front();
      continue;
    }
    
    // insert new point into hull
    T tmp = *nextpoint;
    leftovers.erase(nextpoint);
    hull.insert(hull.begin() + i + 1,tmp);
    
    // update index_length list
    Edge new1 = edges.front();
    edges.pop_front();
    new1.length = (subtract(hull[i],hull[i+1])).length();
    Edge new2(i+1,(subtract(hull[i+1],hull[i+2])).length());
    
    for(auto &p : edges) if(p.index > i) ++(p.index);
    
    auto p = edges.begin();
    if(new1.length > scale){
      for(p = edges.begin() ; p != edges.end() ; ++p){
        if((*p).length < new1.length){
          edges.insert(p,new1);
          break;
        }
      }
      if(p==edges.end()) edges.push_back(new1);
    }
    if(new2.length > scale){
      for(p = edges.begin(); p != edges.end() ; ++p){
        if((*p).length < new2.length){
          edges.insert(p,new2);
          break;
        }
      }
      if(p==edges.end()) edges.push_back(new2);
    }
    
    
    if(TEST){
      PosType cent[] ={0.5,0.5};
      PixelMap map(cent,256, 1.0/256);
      
      map.drawgrid(10,0.5);
      
      map.drawPoints(init_points,0.03,1.0);
      map.drawCurve(hull,2);
      
      map.printFITS("!test_concave.fits");
    }
  }
  
  hull.pop_back();
  
  //Utilities::RemoveIntersections(hull);
  
  return hull;
}

template <typename Ptype>
bool segments_cross(const Ptype &a1,const Ptype &a2
                    ,const Ptype &b1,const Ptype &b2){
  Ptype db= b2 - b1;
  Ptype da= a2 - a1;
  Ptype d1= b1 - a1;
  
  double tmp = (db^da);
  if(tmp==0) return false; // parallel case
  
  double B = (da^d1) / tmp;
  if( (B>1) || (B<0)) return false;
  B = (db^d1) / tmp;
  if( (B>1) || (B<0)) return false;
  return true;
}

template <typename Ptype>
bool inCurve(const Ptype &x,const std::vector<Ptype> &H){
  
  size_t n = H.size();
  if(n <=2) return false;
  long w=0;
  Ptype dH,dD;
  double B,A;
  for(size_t i=0 ; i<n-1 ; ++i){
    dH = H[i+1]-H[i];
    if(dH[1] == 0){  // // horizontal segment
      if(x[1] == H[i][1]){
        B = (x[0]-H[i][0])/dH[0];
        if( B>=0 && B<=1 ) return true;  // point on boundary
      }
    }else{
      dD = x-H[i];
      A = dD[1]/dH[1];
      if( A > 1 || A <= 0) continue;
      B = (dH^dD)/dH[1];
      if(B==0){               // on the boundary
        return true;          // boundaries are always in
      }else if(B>0){
        if(dH[1] > 0) ++w;
        else --w;
      }
    }
  }
  
  return w;
}

template <typename Ptype>
bool inhull(PosType x[],const std::vector<Ptype> &H){
  
  size_t n = H.size();
  if(n <=2) return false;
  long w=0;
  Ptype dH,dD;
  double B,A;
  for(size_t i=0 ; i<n-1 ; ++i){
    dH = H[i+1]-H[i];
    if(dH[1] == 0){  // // horizontal segment
      if(x[1] == H[i][1]){
        B = (x[0]-H[i][0])/dH[0];
        if( B>=0 && B<=1 ) return true;  // point on boundary
      }
    }else{
      dD[0] = x[0]-H[i][0];
      dD[1] = x[1]-H[i][1];
      A = dD[1]/dH[1];
      if( A > 1 || A <= 0) continue;
      B = (dH^dD)/dH[1];
      if(B==0){               // on the boundary
        return true;          // boundaries are always in
      }else if(B>0){
        if(dH[1] > 0) ++w;
        else --w;
      }
    }
  }
  
  return abs(w);
}

template <>
inline bool inhull<Point *>(PosType x[],const std::vector<Point *> &H){
  
  size_t n = H.size();
  if(n <=2) return false;
  long w=0;
  Point_2d dH,dD;
  double B,A;
  for(size_t i=0 ; i<n-1 ; ++i){
    dH = *H[i+1]-*H[i];
    if(dH[1] == 0){  // // horizontal segment
      if(x[1] == H[i]->x[1]){
        B = (x[0]-H[i]->x[0])/dH[0];
        if( B>=0 && B<=1 ) return true;  // point on boundary
      }
    }else{
      dD[0] = x[0]-H[i]->x[0];
      dD[1] = x[1]-H[i]->x[1];
      A = dD[1]/dH[1];
      if( A > 1 || A <= 0) continue;
      B = (dH^dD)/dH[1];
      if(B==0){               // on the boundary
        return true;          // boundaries are always in
      }else if(B>0){
        if(dH[1] > 0) ++w;
        else --w;
      }
    }
  }
  
  return w;
}

/// finds in x is within the curve discribed by the H[].x points ie image points
template <>
inline bool inhull<RAY>(PosType x[],const std::vector<RAY> &H){
  
  size_t n = H.size();
  if(n <=2) return false;
  long w=0;
  Point_2d dH,dD;
  double B,A;
  for(size_t i=0 ; i<n-1 ; ++i){
    dH = H[i+1].x-H[i].x;
    if(dH[1] == 0){  // // horizontal segment
      if(x[1] == H[i].x[1]){
        B = (x[0]-H[i].x[0])/dH[0];
        if( B>=0 && B<=1 ) return true;  // point on boundary
      }
    }else{
      dD[0] = x[0]-H[i].x[0];
      dD[1] = x[1]-H[i].x[1];
      A = dD[1]/dH[1];
      if( A > 1 || A <= 0) continue;
      B = (dH^dD)/dH[1];
      if(B==0){               // on the boundary
        return true;          // boundaries are always in
      }else if(B>0){
        if(dH[1] > 0) ++w;
        else --w;
      }
    }
  }
  
  return w;
}

template <>
inline bool inhull<PosType *>(PosType x[],const std::vector<PosType *> &H){
  
  size_t n = H.size();
  if(n <=2) return false;
  long w=0;
  Point_2d dH,dD;
  double B,A;
  for(size_t i=0 ; i<n-1 ; ++i){
    dH[0] = H[i+1][0]-H[i][0];
    dH[1] = H[i+1][1]-H[i][1];
    if(dH[1] == 0){  // // horizontal segment
      if(x[1] == H[i][1]){
        B = (x[0]-H[i][0])/dH[0];
        if( B>=0 && B<=1 ) return true;  // point on boundary
      }
    }else{
      dD[0] = x[0]-H[i][0];
      dD[1] = x[1]-H[i][1];
      A = dD[1]/dH[1];
      if( A > 1 || A <= 0) continue;
      B = (dH^dD)/dH[1];
      if(B==0){               // on the boundary
        return true;          // boundaries are always in
      }else if(B>0){
        if(dH[1] > 0) ++w;
        else --w;
      }
    }
  }
  
  return w;
}


/*** \brief Calculate the k nearest neighbors concave hull.
 
 This algorithem is guarenteed to find a curve that serounds an island of points, but not all islands unless check==true.
 
 If it at first fails with the k input, k will increase.  The final k relaces the input k.
 
 check determines whether a final check is done to determine whether all the raminaing points are inside the hull.  Otherwise it is possible to get multiple disconnected clusters and the hull will surround only one.
 */

template <typename Ptype>
std::vector<Ptype> concaveK(std::vector<Ptype> &points,int &k,bool check=true)
{
  //std::cout << "finding hull .... ";
  if(points.size() <= 3){
    return points;
  }
  
  if(k  < 3) k =3;
  
  size_t npoints = points.size();
  
  //  {
  //    tree.pop(7);
  //    RandomNumbers_NR ran(123);
  //    std::vector<double> radii;
  //    std::vector<size_t> neighbors;
  //    Point_2d point;
  //    for(int i = 0 ; i<100; ++i){
  //      point[0] = (1-2*ran());
  //      point[1] = (1-2*ran());
  //      tree.NearestNeighbors(point.x,k,radii,neighbors);
  //    }
  //    for(auto p : points){
  //      tree.NearestNeighbors(p.x,k,radii,neighbors);
  //    }
  //  }
  std::vector<Ptype> hull;
  
  double xmin = points[0][0];
  size_t first_point_index = 0;
  for(size_t i=0 ; i<npoints ; ++i){
    if(points[i][0] < xmin){
      xmin = points[i][0];
      first_point_index = i;
    }
  }
  
  Ptype firstpoint = points[first_point_index];
  
  std::vector<size_t> remaining_index;
  bool segmented = true;
  while(segmented){
    bool found=false;
    while(!found){  // if hull leaves points out repeat with larger k
      hull.resize(0);
      TreeSimpleVec<Ptype> tree(points.data(),npoints,2);
      
      if(k>=points.size()){
        hull = convex_hull(points);
        return hull;
      }
      
      hull.push_back(firstpoint);
      // point is not popper frum the free so that it can be found at the end
      
      std::vector<double> radii;
      std::vector<size_t> neighbors;
      Ptype hull_end = hull.back();
      
      Ptype v1;
      Ptype last_point = Ptype(hull[0][0],hull[0][1]-1);
      Ptype v2;
      
      size_t new_index;
      double theta_max=0;
      Ptype new_point;
      std::vector<double> thetas;
      std::vector<int> sorted_index(k);
      
      while((hull[0] != hull.back() && found) || hull.size()==1 ){
        
        tree.NearestNeighbors(hull_end.x,k,radii,neighbors);
        
        long Nneighbors = neighbors.size(); // incase there are not k left
        thetas.resize(Nneighbors);
        sorted_index.resize(Nneighbors);
        
        v1 = hull_end - last_point;
        theta_max = 0;
        for(int i=0; i<Nneighbors ; ++i){
          v2 = points[neighbors[i]] - hull_end;
          double cross = v1^v2,dot = (v1*v2);
          if(cross == 0 && dot <= 0) thetas[i] = -PI;  // prevents backtracking at beginning
          else thetas[i] = atan2( cross , dot );
          sorted_index[i] = i;
        }
        
        std::sort(sorted_index.begin(),sorted_index.end()
                  ,[&thetas](int i,int j){return thetas[i] > thetas[j];} );
        
        int trial=0;
        bool intersect = true;
        while(intersect && trial < Nneighbors){
          
          new_index = neighbors[ sorted_index[trial] ];
          new_point = points[ new_index ];
          
          if(hull.size() <= 3) intersect = false;
          
          // check that new edge doesn't cross earlier edge
          for(long i = hull.size() - 2 ; i > 1 ; --i){
            intersect = segments_cross(hull[i],hull[i-1]
                                       ,hull.back(),new_point);
            
            if(intersect){
              break;
            }
          }
          
          if(new_index == first_point_index && !intersect){
            // check that neighbors are inside hull that would be closed
            // this is to prevent pre-mature closing of the loop
            hull.push_back(new_point);
            for(size_t i : neighbors){
              if(i != new_index && !inCurve(points[i],hull)){
                intersect = true;
                break;
              }
            }
            hull.pop_back();  /// will be added back later
          }
          
          ++trial;
        }
        
        if(intersect){
          //          std::cout << "Intersection ..." << k << std::endl;
          //
          //          std::ofstream logfile("testpoint.csv");
          //          for(auto &p : points){
          //            logfile << p[0] <<","<< p[1] << std::endl;
          //          }
          //          logfile.close();
          //
          //          logfile.open("testhull.csv");
          //          for(auto &p : hull){
          //            logfile << p[0] <<","<< p[1] << std::endl;
          //          }
          //          logfile << new_point[0] <<","<< new_point[1] << std::endl;
          //          logfile.close();
          
          k *= 2;
          found = false;
        }else{
          hull.push_back(new_point);
          last_point = hull_end;
          hull_end = hull.back();
          //tree.print();
          tree.pop(new_index);
          //tree.print();
          found = true;
          
          //          {
          //            Ptype center;
          //            PixelMap map(center.x, 512, 3. / 512.);
          //
          //            map.drawCurve(hull,1);
          //            map.drawPoints(points,0,2);
          //            map.printFITS("!concavek_test"+ std::to_string(f) +".fits");
          //            std::cout << "Test plot : concavek_test" << f++ << ".fits" << std::endl;
          //          }
        }
      }
      remaining_index = tree.get_index();
    }
    // test if all remaining points are in the hull
    
    
    segmented = false;
    if(check){
      for(auto i : remaining_index){
        if(!inCurve(points[i],hull)){
          segmented = true;
          k *= 2;
          //          std::cout << "point outside ... ";
          //          std::ofstream logfile("testpoint.csv");
          //          for(auto &p : points){
          //            logfile << p[0] <<","<< p[1] << std::endl;
          //          }
          //          logfile.close();
          //
          //          logfile.open("testhull.csv");
          //          for(auto &p : hull){
          //            logfile << p[0] <<","<< p[1] << std::endl;
          //          }
          //          logfile.close();
          //
          
          break;
        }
      }
    }
  }
  
  //std::cout << "found hull" << std::endl;
  //hull.pop_back();
  
  return hull;
}


template <typename Ptype>
void testconcaveK(){
  std::vector<Ptype> points(200);
  RandomNumbers_NR ran(2312);
  
  for(Ptype &p : points){
    p[0] = (1-ran());
    p[1] = (1-2*ran());
  }
  
  double x=-1;
  for(int i=0 ; i< points.size()/4 ; ++i){
    points[i][0] = 1;
    points[i][1] = x;
    x += 1/25.;
  }
  x=-1;
  for(int i=points.size()/4 ; i< 2*points.size()/4 ; ++i){
    points[i][0] = -1;
    points[i][1] = x;
    x += 1/25.;
  }
  x=-1;
  for(int i=2*points.size()/4 ; i< 3*points.size()/4 ; ++i){
    points[i][0] = x;
    points[i][1] = 1;
    x += 1/25.;
  }
  x=-1;
  for(int i=3*points.size()/4 ; i< points.size() ; ++i){
    points[i][0] = x;
    points[i][1] = -1;
    x += 1/25.;
  }
  
  //  for(int i=0 ; i< points.size()/2 ; ++i){
  //    points[i][0] = (1-ran()) - 2;
  //    points[i][1] = (1-2*ran());
  //  }
  //
  //
  //  points[0][0] = 0; points[0][1] = -1;
  //  points[1][0] = -1; points[1][1] = 0;
  //  points[2][0] = 0; points[2][1] = 1;
  //  points[3][0] = 1; points[3][1] = 0;
  //
  //  points[4][0] = -1; points[4][1] = 1;
  //  points[5][0] = 1; points[5][1] = 1;
  
  int k = 5;
  std::vector<Ptype> hull = concaveK(points,k);
  Ptype center;
  PixelMap map(center.x, 512, 3. / 512.);
  
  map.drawCurve(hull,1);
  map.drawPoints(points,0,2);
  map.printFITS("!concavek_test.fits");
  std::cout << "Test plot : concavek_test.fits" << std::endl;
}

/// return the interesetion point  of line defined by v1,v2 and w1,w2
Point_2d line_intersection(const Point_2d &v1,const Point_2d &v2,
                           const Point_2d &w1,const Point_2d &w2);


/// returns true if a circle of radius r around the point x intersects with the curve v.  Does not include case where one compleatly encloses the other.
bool circleIntersetsCurve(const Point_2d &x,double r,const std::vector<Point_2d> &v);

/** \brief Returns true if there is any overlap between a curve and a circle.  This includea cases where one compleatly encloses the other.
 */
bool circleOverlapsCurve(const Point_2d &x,double r,const std::vector<Point_2d> &v);

/** \brief Find a curve that is made up of segments from v and w that surrounds them and does not self intersect
 
 v and w must be non-self intersecting
 If they do not intersect and one is not inside the other an empty vector is returned
 */
std::vector<Point_2d> envelope(const std::vector<Point_2d> &v
                               ,const std::vector<Point_2d> &w);


template <typename R>
Point_2d RandomPointWithinCurve(const std::vector<Point_2d> &curve,R &ran){
  
  if(curve.size()==0) throw std::runtime_error("bad curve");
    
  Point_2d p1,p2,center;
  p1=p2=curve[0];
  for(const Point_2d &p : curve){
    if(p[0] < p1[0]){
      p1[0]=p[0];
    }else if(p[0] > p2[0]){
      p2[0]=p[0];
    }
    if(p[1] < p1[1]){
      p1[1]=p[1];
    }else if(p[1] > p2[1]){
      p2[1]=p[1];
    }
  }
  
  if( p1==p2) return p1;
  Point_2d p;
  double area;
  do{
    p[0] = p1[0] + ran()*( p2[0]- p1[0] );
    p[1] = p1[1] + ran()*( p2[1]- p1[1] );

  }while( Utilities::windings(p,curve,&area) == 0 );
  
  return p;
}
/** \brief Return a point that is either within the curve or within a distance r of the curve
 */
template <typename R>
Point_2d RandomPointTouchingCurve(const std::vector<Point_2d> &curve
                                  ,double r
                                  ,R &ran){
  
  if(curve.size()==0) throw std::runtime_error("bad curve");
    
  Point_2d p1,p2,center;
  p1=p2=curve[0];
  for(const Point_2d &p : curve){
    if(p[0] < p1[0]){
      p1[0]=p[0];
    }else if(p[0] > p2[0]){
      p2[0]=p[0];
    }
    if(p[1] < p1[1]){
      p1[1]=p[1];
    }else if(p[1] > p2[1]){
      p2[1]=p[1];
    }
  }

  p1[0] -= r;
  p1[1] -= r;
  p2[0] += r;
  p2[1] += r;

  Point_2d p;
  double area;
  do{
    p[0] = p1[0] + ran()*( p2[0]- p1[0] );
    p[1] = p1[1] + ran()*( p2[1]- p1[1] );
  }while( Utilities::windings(p.x,curve,&area) == 0
         && Utilities::circleOverlapsCurve(p,r,curve) == false
         );
  
  return p;
}

}
#endif /* concave_hull_h */
