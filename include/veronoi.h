//
//  veronoi.h
//  GLAMER
//
//  Created by Ben Metcalf on 27/02/15.
//
//

#ifndef GLAMER_veronoi_h
#define GLAMER_veronoi_h

#include "point.h"

class Point_Nd{
  
  Point_Nd(int Nin){
    dim = Nin;
    x.resize(Nin);
  }
  ~Point_Nd(){x.clear();};
  
  Point_Nd(const Point_Nd &p){
    for(int i=0;i<dim;++i) x[i]=p.x[i];
  }
  Point_Nd & operator=(const Point_Nd &p){
    if(this == &p) return *this;
    for(int i=0;i<dim;++i) x[i]=p.x[i];
    return *this;
  }
  Point_Nd & operator=(const Point &p){
    for(int i=0;i<dim;++i) x[i]=p.x[i];
    return *this;
  }
  
  Point_Nd & operator+=(const Point &p){
    for(int i=0;i<dim;++i) x[i]+=p.x[i];
    return *this;
  }
  Point_Nd  operator+(const Point_Nd &p) const{
    Point_Nd tmp(p.dim);
    for(int i=0;i<dim;++i) tmp.x[i] = x[i] + p.x[i];
    return tmp;
  }
  Point_Nd  operator-(const Point_Nd &p) const{
    Point_Nd tmp(p.dim);
    for(int i=0;i<dim;++i) tmp.x[i] = x[i] - p.x[i];
    return tmp;
  }
  Point_Nd & operator+=(const Point_Nd &p){
    for(int i=0;i<dim;++i) x[i]+=p.x[i];
    return *this;
  }
  Point_Nd & operator/=(PosType value){
    for(int i=0;i<dim;++i) x[i]/=value;
    return *this;
  }
  Point_Nd & operator/(PosType value){
    for(int i=0;i<dim;++i) x[i]/=value;
    return *this;
  }
  Point_Nd & operator*=(PosType value){
    for(int i=0;i<dim;++i) x[i]*=value;
    return *this;
  }
  /// scalar product
  PosType operator*(const Point_Nd &p){
    PosType sum =0;
    for(int i=0;i<dim;++i) sum += x[i]*p.x[i];
    return 0.0;
  }
  /// outer product
  //PosType operator^(const Point_Nd &p){
  //  return x[0]*p.x[1] - x[1]*p.x[0];
  //}
  
  /// length
  PosType length(){
    return sqrt((*this)*(*this));
  }
  
  PosType & operator[](size_t i){return x[i];}
private:
  std::vector<PosType> x;
  int dim;
};


class Simplex_Nd{
  Simplex_Nd(int D,Point_Nd *points){
    vecs.resize(D+1);
    for(int i=0;i<D+1;++i) vecs[i] = &points[i];
  }
  
  Point_Nd & RandomPointWithin();
  double Volume(){return volume;}
private:
  std::vector<Point_Nd *> vecs;
  PosType volume;
  void CalcVolume();
};

#endif
