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

double *estprof(std:: valarray<double> q,int nx,int ny, std:: valarray<double> r, double dr0, 
		double xmax,std:: vector<int> &vi, std:: vector<int> &vj,int ngal);
double *estsigmaprof(std:: valarray<double> q,int nx,int ny, std:: valarray<double> r, double dr0, 
		     double xmax, std:: vector<int> &vi, std:: vector<int> &vj, int ngal, double* qm);
double *estcprof(std:: valarray<double> q,int nx,int ny, std:: valarray<double> r, double dr0, 
		 double xmax, std:: vector<int> &vi, std:: vector<int> &vj, int ngal);
double *estsigmacprof(std:: valarray<double> q,int nx,int ny, std:: valarray<double> r, double dr0, 
		      double xmax, std:: vector<int> &vi, std:: vector<int> &vj, int ngal, double* qm);

#endif /* PROFILE_H_ */
