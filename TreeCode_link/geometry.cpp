//
//  geometry.cpp
//  GLAMER
//
//  Created by bmetcalf on 7/24/14.
//
//

#include "geometry.h"
#include "Tree.h"
#include "point.h"

bool Utilities::Geometry::intersect(const PosType a1[],const PosType a2[],const PosType b1[],const PosType b2[]){
  
  if(a1 == b1 || a1 == b2) return false;
  if(a2 == b1 || a2 == b2) return false;
  
  int o1 = Utilities::Geometry::orientation(a1, a2, b1);
  int o2 = Utilities::Geometry::orientation(a1, a2, b2);
  int o3 = Utilities::Geometry::orientation(b1, b2, a1);
  int o4 = Utilities::Geometry::orientation(b1, b2, a2);
  
  // General case
  if (o1 != o2 && o3 != o4)
    return true;
  
  // Special Cases
  // p1, q1 and p2 are colinear and p2 lies on segment p1q1
  if (o1 == 0 && Utilities::Geometry::onSegment(a1, b1, a2)) return true;
  
  // p1, q1 and p2 are colinear and q2 lies on segment p1q1
  if (o2 == 0 && Utilities::Geometry::onSegment(a1, b2, a2)) return true;
  
  // p2, q2 and p1 are colinear and p1 lies on segment p2q2
  if (o3 == 0 && Utilities::Geometry::onSegment(b1, a1, b2)) return true;
  
  // p2, q2 and q1 are colinear and q1 lies on segment p2q2
  if (o4 == 0 && Utilities::Geometry::onSegment(b1, a2, b2)) return true;
  
  return false; // Doesn't fall in any of the above cases
}

int Utilities::Geometry::intersect(const std::vector<Point_2d> &curve){
  if(curve.size() < 4) return 0;
  
  int intersections = 0;
  for(size_t i=0;i<curve.size()-2;++i){
    for(size_t j=i+2;j<curve.size()-1;++j){
      intersections += Utilities::Geometry::intersect(curve[i].x,curve[i+1].x,curve[j].x,curve[j+1].x);
    }
    intersections += Utilities::Geometry::intersect(curve[i].x,curve[i+1].x
                                                    ,curve.back().x,curve[0].x);
  }
  
  return intersections;
}

int Utilities::Geometry::intersections(const std::vector<Point_2d> &curve,std::vector<Point_2d> &intersections,std::vector<std::pair<int,int> > &segments){
  
  intersections.clear();
  segments.clear();
  
  if(curve.size() < 4) return 0;
  
  int n = 0;
  for(size_t i=0;i<curve.size()-2;++i){
    for(size_t j=i+2;j<curve.size()-1;++j){
      if(Utilities::Geometry::intersect(curve[i].x,curve[i+1].x,curve[j].x,curve[j+1].x)){
        segments.emplace_back(i,j);
        Point_2d dp = curve[i+1] - curve[i];
        Point_2d dq = curve[j+1] - curve[j];
        double u = ( (curve[i+1]^curve[i]) - (dp^curve[j]) ) / (dp^dq) ;
        intersections.push_back( dq*u + curve[j] );
        ++n;
      }
    }
  }
  
  return n;
}


int Utilities::Geometry::orientation(const PosType p[],const PosType q[],const PosType r[])
{
  // See 10th slides from following link for derivation of the formula
  double val = (q[1] - p[1]) * (r[0] - q[0]) -
  (q[0] - p[0]) * (r[1] - q[1]);
  
  if (val == 0.0) return 0;  // colinear
  
  return (val > 0)? 1: 2; // clock or counterclock wise
}
/** \brief Given three colinear points p, q, r, the function checks if
 point q lies on line segment 'pr', but not at p or r
 */
bool Utilities::Geometry::onSegment(const PosType p[], const PosType q[], const PosType r[])
{
  if (q[0] < MAX(p[0], r[0]) && q[0] > MIN(p[0], r[0]) &&
      q[1] < MAX(p[1], r[1]) && q[1] > MIN(p[1], r[1]))
    return true;
  
  return false;
}

double Utilities::Geometry::AngleBetween2d(double v1[],double v2[]){
  double y = (v1[0] * v2[1]) - (v2[0] * v1[1]);
  double x = (v1[0] * v2[0]) + (v1[1] * v2[1]);

  if(y == 0 && x < 0 ) return PI;
  return atan2(y, x);
}

//std::vector<Point_2d> Utilities::Geometry::MagicHull(const std::vector<Point_2d> &points){
//
//  std::vector<Point_2d> hull;
//  if(points.size() == 0) return hull;
//
//  // find left most point
//  long init = 0,Npoints=points.size();
//  double pmin = points[0][0];
//  for(int i=1 ; i<Npoints ; ++i){
//    if(points[i][0] < pmin){
//      init = i;
//      pmin = points[i][0];
//    }
//  }
//
//  CYCLIC cyc(Npoints);
//  int orientation = sign( (points[cyc[init-1]] - points[init])^(points[cyc[init+1]] - points[init])  );
//
//  hull.reserve(3*Npoints);
//
//  long k=init,count=0;
//  int direction = 1;
//  Point_2d intersection_point;
//  long next_segment;
//  Point_2d current_point = points[init];
//  hull.push_back(points[init]);
//
//  k = cyc[init + 1];
//  long current_segment = init;
//  while(count < 3*Npoints ) {
//
//    int Nintersect = 0;
//
//    {
//      long j=0,jp;
//      double minimum=HUGE_VALF;
//      // find closest itersection if any
//      while(j<Npoints){
//        jp = cyc[j+1];
//
//        if(Utilities::Geometry::intersect(current_point.x,points[k].x,points[j].x,points[jp].x) ){
//
//          Point_2d dp = points[k] - current_point;
//          Point_2d dq = points[jp] - points[j];
//          //double v = ( (points[k]^current_point) - (dp^points[j]) ) / (dp^dq) ;
//          double u = ( dq^(points[j]-current_point) ) / (dq^dp) ;
//
//          Point_2d p = dp*u + current_point;
//
//          if (abs(u) > 1.0e-10){ // requires the current_point not be on the next edge
//            if ( u < minimum ){
//              intersection_point = p;
//              minimum = u;
//              next_segment = j;
//            }
//            ++Nintersect;
//          }
//          //else{
//          //  std::cout << u/dp.length() << "," << k <<","<< j << "," << jp << "," << current_segment << std::endl;
//          //}
//        }
//        ++j;
//      }
//    }
//
//    if(Nintersect==0){      // no intersection
//      if(k==init) break;
//
//      hull.push_back(points[k]);
//      current_segment = k;
//      k = cyc[k + direction];
//    }else{
//
//      hull.push_back(intersection_point);
//      if(orientation == sign( (points[next_segment] - intersection_point)^(points[ k ] - intersection_point) ) ){
//        k = cyc[next_segment + 1];
//        direction = 1;
//      }else{
//        k = next_segment;
//        direction = -1;
//      }
//    }
//    current_point = hull.back();
//    ++count;
//  }
//
//  hull.shrink_to_fit();
//
//  return hull;
//}

/** test code for MagicHull
 
 {
   //Utilities::RandomNumbers_NR ran(128234);
   //Utilities::RandomNumbers_NR ran(2328234);
   //Utilities::RandomNumbers_NR ran(231277);
   //Utilities::RandomNumbers_NR ran(123477);
   //Utilities::RandomNumbers_NR ran(1233223477);
   Utilities::RandomNumbers_NR ran(1233127);

   std::vector<Point_2d> points(10);
   for(Point_2d &p : points){
     p[0] = ran();
     p[1] = ran();
   }
   //points.push_back(points[0]);
   
   /////////////////////////////////////////////////////////////////////////////////////
   
   //std::vector<Point_2d> hull = points;
   //std::vector<Point_2d> intersections;
   //std::vector<std::pair<int,int> > segments;
   
   //Utilities::Geometry::intersections(points,intersections,segments);
   
   assert(points[0] != points.back());
   // orientation == 1 is clockwise
 
   std::vector<Point_2d> hull = Utilities::Geometry::MagicHull(points);
   
   std::ofstream file("original.csv");
   file << "x,y" << std::endl;
   for(Point_2d &p : points) file << p[0] <<"," << p[1] << std::endl;
   
   std::ofstream file2("hull.csv");
   file2 << "x,y" << std::endl;
   for(Point_2d &p : hull) file2 << p[0] <<"," << p[1] << std::endl;

   exit(0);
 }
 */
