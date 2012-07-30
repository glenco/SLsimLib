/*
 * fits.h
 *
 *  Created on: Jun 19, 2012
 *      Author: mpetkova
 */

#ifdef WITH_MOKA

#ifndef FITS_H_
#define FITS_H_

#include <iostream>
#include <valarray>
#include <vector>

struct LensHalo{
  double zl,m,zs,DL,DLS,DS,c,cS,fsub,mstar,minsubmass;   // lens and source properties
  int nsub;                                              // subhalo number
  double ea,eb,ec;                                       // axes 
  double Ermed,Ereff;                                    // median and effective Einstein radius
  double beta;                                           // inner slope of the main halo density profile
  double boxlMpc,boxlarcsec;                             // boxsize
  double omegam,omegal,h,wq;                             // cosmology
  int npix;                                              // number of pixels
};

void getDims(std::string fn
		,int *nx
		,int *ny);

void readImage(std::string fn
		,std::valarray<float> *convergence
		,std::valarray<float> *alpha1
		,std::valarray<float> *alpha2
		,std::valarray<float> *gamma1
		,std::valarray<float> *gamma2
	       ,struct LensHalo *LH);

void writeImage(std::string filename
		,std::valarray<float> convergence
		,std::valarray<float> gamma1
		,std::valarray<float> gamma2
		,std::valarray<float> gamma3
		,int nx
		,int ny
		,struct LensHalo LH);


void make_friendship(int ii,int ji,int np,std:: vector<int> &friends, std:: vector<double> &pointdist);

int fof(double l,std:: vector<double> xci, std:: vector<double> yci, std:: vector<int> &groupid);

#endif /* FITS_H_ */

#endif
