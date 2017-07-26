//
//  profiles.cpp
//  GLAMER
//
//  Created by bmetcalf on 28/02/17.
//
//

#include "profiles.h"
#include "geometry.h"
#include "utilities_slsim.h"
//#include "Tree.h"

namespace Profiles{
  
  Utilities::LogLookUpTable<double> NFW2D::flookup(NFW2D::ffunction,1.0e-3,20,500);
  Utilities::LogLookUpTable<double> NFW2D::glookup(NFW2D::ggfunction,1.0e-3,20,500);
  
  double NFW2D::ffunction(PosType x){
    
    if(x==0) x=1e-5;
    if(x==1.0){ return 1.0/3.0/pi;}
    if(x>1.0){  return 0.5*(1-2*atan(sqrt((x-1)/(x+1)))/sqrt(x*x-1))/(x*x-1)/pi; }
    if(x<1.0){  return 0.5*(1-2*atanh(sqrt((1-x)/(x+1)))/sqrt(1-x*x))/(x*x-1)/pi;}
    return 0.0;
  }
  double NFW2D::ggfunction(PosType x){
    
    if(x<1e-5) x=1e-5;
    PosType ans=log(x/2);
    if(x==1.0){ ans += 1.0; return ans;}
    if(x>1.0){  ans += 2*atan(sqrt((x-1)/(x+1)))/sqrt(x*x-1); return ans;}
    if(x<1.0){  ans += 2*atanh(sqrt((1-x)/(x+1)))/sqrt(1-x*x); return ans;}
    return 0.0;
  }


  BsplineGEN::BsplineGEN(long seed):ran(seed){
    X.resize(N);
    F.resize(N);
    X[0] = F[0] = 0.0;
    dx = 2.0/(N-1);
    for(int i=1;i<N;++i){
      X[i] = i*dx;
      F[i] = mass_frac(X[i]);
    }
  }
  
  void BsplineGEN::draw(std::vector<Point_3d> &v){
    
    for(auto &p : v) draw(p);
  };
  
  void BsplineGEN::draw(Point_3d &p){
    
    double f = ran();
    size_t i = Utilities::locate(F,f);
    double x = X[i] + dx*(f - F[i])/(F[i+1] - F[i]);
    
    double theta = 2*pi*ran();
    double z = 2*ran() - 1;
    double t = sqrt(1-z*z);
    p[0] = t*cos(theta);
    p[1] = t*sin(theta);
    p[2] = z;
    
    p *= x;
  };
  
  
  double BsplineGEN::mass_frac(double q){
    
    if(q > 2) return 1.0;
    PosType q2=q*q,q3=q2*q;
    
    if(q > 1) return -0.0666667 + q3*( 8./3. - 3.*q + 6.*q2/5. - q3/6.);
    
    return 4*q3*(0.333333 - 0.3*q2 + 0.125*q3);
  };

  /** \brief class for generating positions in proportion to mass in an NFW profiles
   */
  NFWgenerator::NFWgenerator(Utilities::RandomNumbers_NR &ran_in,double max_cons)
  :ran(ran_in)
  {
    X.resize(N);
    F.resize(N);
    X[0] = F[0] = 0.0;
    dx = max_cons/(N-1);
    for(int i=1;i<N;++i){
      X[i] = i*dx;
      F[i] = log(1+X[i]) - X[i]/(1+X[i]) ;
    }
  }
  
  /// returns a vector of points drawn from a spherical halo
  void  NFWgenerator::drawSpherical(std::vector<Point_3d> &points  /// output points
                                    ,double cons                   /// concentration
                                    ,double Rvir                   /// maximum elliptical radius
  ){
    
    double Fmax = log(1+cons) - cons/(1+cons);
    double rs = Rvir/cons;
    for(auto &p : points){
      double f = Fmax*ran();
      size_t i = MIN<size_t>(Utilities::locate(F,f),N-2);
      double x = X[i] + dx*(f - F[i])/(F[i+1] - F[i]);
      double theta = 2*pi*ran();
      
      p[0] = 2*ran() - 1;
      double tmp = sqrt(1-p[0]*p[0]);
      p[1] = tmp*cos(theta);
      p[2] = tmp*sin(theta);
      
      p *= x*rs;
    }
  }
  
  ///  STILL UNDER CONSTRUCTION returns a vector of points drawn from a triaxial halo,
  void  NFWgenerator::drawTriAxial(std::vector<Point_3d> &points  /// output points
                                   ,double cons                   /// concentration
                                   ,double Rvir                   /// maximum elliptical radius
                                   ,double f1                     /// axis ratio 1 to 3
                                   ,double f2                     /// axis ratio 2 to 3
                                   ,Utilities::Geometry::SphericalPoint v              /// direction of axis 3
  ){
    
    double a3 = 1.0/pow(f1*f2,1.0/3.0);
    double a1 = f1*a3,a2 = f2*a3;
    Utilities::Geometry::Quaternion rot = Utilities::Geometry::Quaternion::q_z_rotation( v.phi )*Utilities::Geometry::Quaternion::q_y_rotation( v.theta );
    
    double Fmax = log(1+cons) - cons/(1+cons);
    double rs = Rvir/cons;
    for(auto p : points){
      
      p[0] = ran.gauss();
      p[1] = ran.gauss();
      p[2] = ran.gauss();
      
      double f = Fmax*ran();
      size_t i = Utilities::locate(F,f);
      double x = X[i] + dx*(f - F[i])/(F[i+1] - F[i]);
      
      p *= x*rs/p.length();
      
      p[0] *= a1;
      p[1] *= a2;
      p[2] *= a3;
      
      p = rot.Rotate(p);
    }
  }

}
