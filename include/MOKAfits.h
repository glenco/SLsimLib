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

void getDims(std::string fn
		,int *nx
		,int *ny);

void readImage(std::string fn
		,std::valarray<float> *convergence
		,std::valarray<float> *alpha1
		,std::valarray<float> *alpha2
		,std::valarray<float> *gamma1
		,std::valarray<float> *gamma2
		,double *boxl
		,double *boxlMpc
		,double *zlens
		,double *zsource
		,double *omegam
		,double *omegal
		,double *h
		,double *DL);

void writeImage(std::string filename
		,std::valarray<float> convergence
		,std::valarray<float> gamma1
		,std::valarray<float> gamma2
		,std::valarray<float> gamma3
		,int nx
		,int ny
		,double boxl
		,double boxlMpc
		,double zlens
		,double zsource
		,double omegam
		,double omegal
		,double h
		,double DL);

#endif /* FITS_H_ */

#endif
