/*
 * fitlens.h
 *
 *  Created on: Jan 13, 2011
 *      Author: RB Metcalf
 */

#ifndef FITLENS_H_
#define FITLENS_H_

#include <mutex>
#include <analytic_lens.h>


void find_lens(int Nimages,int Nsources,int *pairing,double **xob,double *xg,double beta
               ,int N,int *degen,double *mod,double **v,double **dx);
double modfind(double theta);

void deflection_total(double *ximage,double *angle,int fulllens,double sdeg,double *ang_lens
			,int mag,double *Amag);
double deflection_model(double beta,double *mod,double *x,double *y,double *mag,int N
		    ,int Nlenses,double Re2,double *x2);
int find_image_number(double *yo,double *x_center,double *mod,int Nmod,int Nlenses,double Re2,double *x2);
double finiteMag(double radsource,double *xo,double *mod,int Nmod,int Nlenses,double Re2,double *x2t);

double minEllip(double *q);
double minaxis(double thetaX);
void RotateModel(double thetaX,double *mod,int N,int Nsources);
double regularize(int Nmax,int Nmin,int N,int Nsources,int degen,double *mod
		    ,double **v,double *modo);

#endif
