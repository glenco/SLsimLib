/*
 * profile.h
 *
 *  Created on: May 23, 2012
 *      Author: cgiocoli
 */

#ifndef PROFILE_H_
#define PROFILE_H_

#include <vector>
#include <assert.h>
#include <valarray>
#include <iostream>

double * estprof(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, double dr0, double xmax);
double * estsigmaprof(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, double dr0, double xmax, double* qm);
double * estcprof(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, double dr0, double xmax);
double * estsigmacprof(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, double dr0, double xmax, double* qm);

#endif /* PROFILE_H_ */