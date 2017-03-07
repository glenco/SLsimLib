/*
 * profiles.h
 *
 *  Created on: 2017
 *      Author: bmetcalf
 */

#ifndef PROFILES_H_
#define PROFILES_H_

#include "point.h"
#include "geometry.h"
#include "utilities_slsim.h"


namespace Profiles {
  
  /// cubic b-spline kernel in different dimensions normalized to 1
  template <int d>
  double Bspline(double q){
  }
  template <>
  inline double Bspline<1>(double q){
    if(q >= 2) return 0;
    double q2 = 2-q,q1 = 1-q;
    if(q <= 1) return 1.5*(0.25*q2*q2*q2 - q1*q1*q1);
    return q2*q2*q2/6;
  }
  template <>
  inline double Bspline<2>(double q){
    if(q >= 2) return 0;
    double q2 = 2-q,q1 = 1-q;
    if(q <= 1) return (0.25*q2*q2*q2 - q1*q1*q1)*0.3183098861837907;
    return q2*q2*q2*0.07957747154594767;
  }
  template <>
  inline double Bspline<3>(double q){
    if(q >= 2) return 0;
    double q2 = 2-q,q1 = 1-q;
    if(q <= 1) return (0.25*q2*q2*q2 - q1*q1*q1)/pi;
    return q2*q2*q2*0.07957747154594767;
  }
  
  /** \breaf Two dimensional profile of a three dimensional NFW profile trincated in 3 D.
   
   This is the surface density in units of M rs^2.
   **/
  class TNFW2D{
    
  public:
    TNFW2D(double cons):c(cons){
      c2=c*c;
      cp1 = c+1;
      To = sqrt((c-1)/(cp1));
      
      M = (log(cp1) - c/(cp1) )/2/pi;
    };
    
    double operator()(double x  /// x=r/rs
                      ){
      
      if(x>c) return 0.0;
      double x1 = sqrt(fabs(x*x - 1 ));
      double a = sqrt(c2 - x*x )/ x1 ;
      double T;
      if(x>1) T = (atan(a) - atan(a/c))/x1 ;
        else if(x<1) T = 0.5*log( (1-a)*(1+a/c)/(1-a/c)/(1+a) )/x1 ;
          else T = To;
        
      return ( (sqrt(c2 - x*x )/cp1 - T)/(x*x -1) )/M;
    };
    
  private:
    double c;
    double c2;
    double cp1;
    double M;
    double To;
  };
  
  class BsplineGEN{
  public:
    BsplineGEN(long seed);
    
    /// returns a vector of positions with length between 0 and 2
    void draw(std::vector<Point_3d> &v);
    void draw(Point_3d &v);
    
  private:
    Utilities::RandomNumbers_NR ran;
    std::vector<double> X;
    std::vector<double> F;
    const int N = 1000;
    double dx;
    
    double mass_frac(double q);
  };
  
  /** \brief class for generating positions in proportion to mass in an NFW profiles
   */
  
  using Utilities::Geometry::SphericalPoint;
  
  class NFWgenerator{
  public:
    NFWgenerator(Utilities::RandomNumbers_NR &ran_in,double max_cons);
    /// returns a vector of points drawn from a spherical halo
    void drawSpherical(std::vector<Point_3d> &points  /// output points
                       ,double cons                   /// concentration
                       ,double Rvir                   /// maximum elliptical radius
    );
    ///  STILL UNDER CONSTRUCTION returns a vector of points drawn from a triaxial halo,
    void drawTriAxial(std::vector<Point_3d> &points  /// output points
                      ,double cons                   /// concentration
                      ,double Rvir                   /// maximum elliptical radius
                      ,double f1                     /// axis ratio 1 to 3
                      ,double f2                     /// axis ratio 2 to 3
                      ,SphericalPoint v              /// direction of axis 3
    );
  private:
    Utilities::RandomNumbers_NR &ran;
    double dx;
    std::vector<double> X;
    std::vector<double> F;
    const int N = 1000;
  };

  /// projected NFW profile, not truncated
  class NFW2d{
  public:
    
    NFW2d():func(gf,1.0e-3,20,100){}
    
    double operator()(double x){
      x = (x > 1.0e-3) ? x : 1.0e-3;
      double y;
      func(x,y);
      return y;
    }
  private:
    
    Utilities::LogLookUpTable<double> func;
    
    static double gf(double x){
      double ans;
      
      if(x<1e-5) x=1e-5;
      ans=log(x/2);
      if(x==1.0){ ans += 1.0; return ans;}
      if(x>1.0){  ans +=  2*atan(sqrt((x-1)/(x+1)))/sqrt(x*x-1); return ans;}
      if(x<1.0){  ans += 2*atanh(sqrt((1-x)/(x+1)))/sqrt(1-x*x); return ans;}
      return 0.0;
    }
  };
  
}
#endif /* PROFILES_H_ */
