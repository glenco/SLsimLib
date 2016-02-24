//
//  concave_hull.h
//  GLAMER
//
//  Created by bmetcalf on 23/02/16.
//
//

#ifndef concave_hull_h
#define concave_hull_h

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
            tmp = curve[k];
            curve[k] = curve[l];
            curve[l] = tmp;
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

  
  /// Returns a vector of points on the convex hull in counter-clockwise order.
  template<typename T>
  void convex_hull(std::vector<T> &P,std::vector<T> &hull_out)
  {
    
    if(P.size() <= 3){
      hull_out = P;
      return;
    }
    
    std::vector<T> hull;
    
    size_t n = P.size();
    size_t k = 0;
    hull.resize(2*n);
    
    // Sort points lexicographically
    std::sort(P.begin(), P.end(),
              [](T &p1,T &p2){
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
    
    std::swap(hull,hull_out);
    
    return;
  }
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
  
  std::vector<T> hull;
  
  // find the convex hull
  convex_hull(init_points,hull);
  
  if(init_points.size() == hull.size()) return;
  
  std::list<std::pair<size_t,double>> index_length;
  
  std::list<T> leftovers;
  hull.push_back(hull[0]);
  
  double tmp;
  
  for(int i=0;i<hull.size()-1;++i){
    tmp = (subtract(hull[i],hull[i+1])).length();
    if(tmp > scale) index_length.push_back(std::pair<size_t,double>(i,tmp));
  }
  
  if(index_length.size() == 0) return;
  
  if(TEST){
    PosType cent[] ={0.5,0.5};
    PixelMap map(cent,256, 1.0/256);
    
    map.drawgrid(10,0.5);
    
    map.drawPoints(init_points,0.01,1.0);
    map.drawCurve(hull,2);
    
    
    map.printFITS("!test_concave.fits");
  }
  
  
  
  
  // sort edges by length
  index_length.sort([](const std::pair<size_t,double> &p1,const std::pair<size_t, double> &p2){return p1.second > p2.second;});
  assert(index_length.front().second >= index_length.back().second);
  
  
  {   // make a list of points that are not already in the convex hull
    
    hull.pop_back();
    std::vector<size_t> index(hull.size());
    size_t i = 0;
    for(size_t &ind: index) ind = i++;
    std::sort(index.begin(),index.end(),
              [&hull](size_t &p1,size_t &p2){
                if(hull[p1][0] == hull[p2][0]) return hull[p1][1] < hull[p2][1];
                return hull[p1][0] < hull[p2][0];});
    
    size_t j = 0;
    
    for(auto pit = init_points.begin(); pit != init_points.end() ; ++pit){
      
      if((hull[index[j]][0] == (*pit)[0])*(hull[index[j]][1] == (*pit)[1])){
        ++j;
      }else{
        //leftovers.push_back(Point_2d((*pit)[0],(*pit)[1]));
        leftovers.push_back(*pit);
      }
      
    }
    assert(hull.size() + leftovers.size() == init_points.size());
    hull.push_back(hull[0]);
  }
  
  typename std::list<T>::iterator p,nextpoint;
  //auto p = leftovers.begin();
  
  double minarea,area,co;
  size_t i;
  Point_2d vo,v1;
  
  while(index_length.front().second > scale){
    
    if(leftovers.size() == 0) break;
    
    // find the point which if added to this edge would change the area least
    minarea = HUGE_VAL;
    i = index_length.front().first;
    
    for(p = leftovers.begin() ; p != leftovers.end() ; ++p){
      
      vo = subtract(hull[i+1], hull[i]);
      v1 = subtract(*p,hull[i]);
      
      area = vo^v1;
      co = v1*vo;
      
      if(co > 0 && area > 0 && index_length.front().second > v1.length()){
        if(area < minarea){
          minarea = area;
          nextpoint = p;
        }
      }
    }
    
    if(minarea == HUGE_VAL){
      // if there is no acceptable point for this edge continue with the second longest edge
      index_length.pop_front();
      continue;
    }
    
    // insert new point into hull
    T tmp = *nextpoint;
    leftovers.erase(nextpoint);
    hull.insert(hull.begin() + i + 1,tmp);
    
    // update index_length list
    std::pair<size_t,double> new1 = index_length.front();
    index_length.pop_front();
    new1.second = (subtract(hull[i],hull[i+1])).length();
    std::pair<size_t,double> new2(i+1,(subtract(hull[i+1],hull[i+2])).length());
    
    for(auto &p : index_length) if(p.first > i) ++(p.first);
    
    auto p = index_length.begin();
    for(; p != index_length.end() ; ++p){
      if((*p).second < new1.second){
        index_length.insert(p,new1);
        break;
      }
    }
    if(p==index_length.end()) index_length.push_back(new1);
    
    for(p = index_length.begin(); p != index_length.end() ; ++p){
      if((*p).second < new2.second){
        index_length.insert(p,new2);
        break;
      }
    }
    if(p==index_length.end()) index_length.push_back(new2);
    
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
  
  Utilities::RemoveIntersections(hull);
  
  std::swap(hull,hull_out);
  
  return;
}

}
#endif /* concave_hull_h */
