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

int Utilities::Geometry::incurve(PosType x[],std::vector<double *> &curve){
  int number = 0;
  size_t i;
  
  Point point;
  for(i=0;i<curve.size()-1;++i){
    
    if( (x[1] >= curve[i][1])*(x[1] <= curve[i+1][1]) ){
      if(Utilities::Geometry::orientation(curve[i], x, curve[i+1]) <= 1) ++number;
    }else if( (x[1] <= curve[i][1])*(x[1] > curve[i+1][1]) ){
      if(Utilities::Geometry::orientation(curve[i], x, curve[i+1]) == 2) --number;
    }
    
  }
  
  if( (x[1] >= curve[i][1])*(x[1] <= curve[0][1]) ){
    if(Utilities::Geometry::orientation(curve[i], x, curve[0]) <= 1) ++number;
  }else if( (x[1] <= curve[i][1])*(x[1] > curve[0][1]) ){
    if(Utilities::Geometry::orientation(curve[i], x, curve[0]) == 2) --number;
  }
  
  return number == 0 ? 0 : 1;
}


