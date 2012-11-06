/*
 * fits.h
 *
 *  Created on: Jun 19, 2012
 *      Author: mpetkova
 */


#ifndef FITS_H_
#define FITS_H_

#include "standard.h"

//TODO Carlo/Margarita: What is this structure for?  It would be better to rename it because we have other structures that start with Lens and this is not related.
struct LensHalo{
	/// lens and source properties
  double zl,m,zs,DL,DLS,DS,c,cS,fsub,mstar,minsubmass;
  /// subhalo number
  int nsub;
  /// axes
  double ea,eb,ec;
  /// median Einstein radius
  double Ermed;
  /// effective Einstein radius
  double Ereff;
  /// inner slope of the main halo density profile
  double beta;
  /// boxsize
  double boxlMpc,boxlarcsec;
  /// cosmology
  double omegam,omegal,h,wq;
  /// number of pixels
  int npix;
};

//TODO Could these be made into methods for a class?  It would help to organize things more clearly.
void getDims(std::string fn
		,int *nx
		,int *ny);

void readImage(std::string fn
		,std::valarray<float> *convergence
		,std::valarray<float> *alpha1
		,std::valarray<float> *alpha2
		,std::valarray<float> *gamma1
		,std::valarray<float> *gamma3
		,LensHalo *LH);

void writeImage(std::string filename
		,std::valarray<float> convergence
		,std::valarray<float> gamma1
		,std::valarray<float> gamma2
		,std::valarray<float> gamma3
		,int nx
		,int ny
		,LensHalo *LH);

void writeImage(std::string filename
		,std::valarray<float> convergence
		,std::valarray<float> gamma1
		,std::valarray<float> gamma2
		,std::valarray<float> gamma3
		,int nx
		,int ny);


void make_friendship(int ii,int ji,int np,std:: vector<int> &friends, std:: vector<double> &pointdist);

int fof(double l,std:: vector<double> xci, std:: vector<double> yci, std:: vector<int> &groupid);

#endif /* FITS_H_ */

