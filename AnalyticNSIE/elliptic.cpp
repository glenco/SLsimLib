//
//  elliptic.cpp
//  GLAMER
//
//  Created by bmetcalf on 8/5/14.
//
//

#include "elliptic.h"
#include "utilities_slsim.h"

PosType Elliptic::DALPHAXDM::operator()(PosType logm){
  
  double m=exp(logm);
  double ap = m*m*a2 + lambda,bp = m*m*b2 + lambda;
  double p2 = x[0]*x[0]/ap/ap/ap/ap + x[1]*x[1]/bp/bp/bp/bp;  // actually the inverse of equation (5) in Schramm 1990
  
  //return m*isohalo->kappa(x)/(ap*ap*ap*bp*p2);
  PosType alpha[2]={0,0},tmp[2] = {m*(isohalo->get_Rmax()),0};
  KappaType kappa=0,gamma[2]={0,0},phi;
  
  isohalo->force_halo(alpha,&kappa,gamma,&phi,tmp);
  assert(kappa >= 0.0);
  //return m*kappa/(x[0]*x[0] + x[1]*x[1]);
  return m*m*kappa/(ap*ap*ap*bp*p2);
}

PosType Elliptic::DALPHAYDM::operator()(PosType m){
  
  double ap = m*m*a2 + lambda,bp = m*m*b2 + lambda;
  double p2 = x[0]*x[0]/ap/ap/ap/ap + x[1]*x[1]/bp/bp/bp/bp;  // actually the inverse of equation (5) in Schramm 1990
  
  PosType alpha[2]={0,0},tmp[2] = {m*(isohalo->get_Rmax()),0};
  KappaType kappa=0,gamma[2]={0,0},phi;
  
  isohalo->force_halo(alpha,&kappa,gamma,&phi,tmp);
  //return m*kappa/(x[0]*x[0] + x[1]*x[1]);
  return m*kappa/(ap*bp*bp*bp*p2);
}

void Elliptic::alpha(PosType x[],PosType alpha[]){
  
  if(x[0] == 0.0 && x[1] == 0.0){
    alpha[0] = alpha[1] =0.0;
    return;
  }
  std::cerr << "Warning: Elliptic::alpha() seems to give a value that is a factor ~ 1.27 too large." << std::endl;
  double a2=a*a,b2 = b*b;
  double tmp = (a2 + b2 - x[0]*x[0] - x[1]*x[1]);
  double c = cos(inclination),s = sin(inclination);
  double xtmp[2] = {x[0]*c - x[1]*s
                   ,x[0]*s + x[1]*c};
  
  //std::cout << "xtmp = " << xtmp[0] << " " << xtmp[1] << std::endl;
  
  double lambda = (tmp + sqrt(tmp*tmp + 4*(x[0]*x[0]*b2 + x[1]*x[1]*a2 - a2*b2 )))/2;
  //double lambda = (tmp - sqrt(tmp*tmp + 4*(x[0]*x[0]*b2 + x[1]*x[1]*a2 - a2*b2 )))/2;
  assert(lambda == lambda);
  //std::cout << "lambda = " << lambda << std::endl;
  
  PosType mo = sqrt(xtmp[0]*xtmp[0]/a2 + xtmp[1]*xtmp[1]/b2);

  //std::cout << "mo = " << mo << std::endl;

  Elliptic::DALPHAXDM funcX(lambda,a2,b2,xtmp,isohalo);
  //alpha[0] = -8*a*b*xtmp[0]*Utilities::nintegrate<Elliptic::DALPHAXDM,PosType>(funcX,0,MIN(mo,1.0),1.0e-6)/pi;
  alpha[0] = -8*a*b*xtmp[0]*Utilities::nintegrate<Elliptic::DALPHAXDM,PosType>(funcX,log(1.0e-7),log(MIN(mo,1.0)),1.0e-6)/pi;
  
  
  Elliptic::DALPHAYDM funcY(lambda,a2,b2,xtmp,isohalo);
  alpha[1] = -8*a*b*xtmp[1]*Utilities::nintegrate<Elliptic::DALPHAYDM,PosType>(funcY,0,MIN(mo,1.0),1.0e-6)/pi;
  
  //std::cout << "alpha = " << alpha[0] << " " << alpha[1] << std::endl;

  /// rotate the
  tmp = alpha[0];
  alpha[0] = alpha[0]*c + alpha[1]*s;
  alpha[1] = -tmp*s + alpha[1]*c;
};

