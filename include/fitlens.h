/*
 * fitlens.h
 *
 *  Created on: Jan 13, 2011
 *      Author: RB Metcalf
 */

#ifndef FITLENS_H_
#define FITLENS_H_

#include <analytic_lens.h>

double ElliptisizeLens(int Nimages,int Nsources,int Nlenses,int *pairing,double **xob
		       ,double *xc,double **xg,double sigG,double beta,int Nmod
		       ,double *mod,double **dx,double *re2,double *q);
void find_lens(int Nimages,int Nsources,int *pairing,double **xob,double *xg,double beta
		 ,int N,int *degen,double *mod,double **v,double **dx);
double minEllip(double *q);
double find_axis(double *mod,int Nmod);
void RotateModel(double thetaX,double *mod,int N,int Nsources);
void find_lens(int Nimages,int Nsources,int *pairing,double **xob,double *xg,double beta
		 ,int N,int *degen,double *mod,double **v,double **dx);
double regularize(int Nmax,int Nmin,int N,int Nsources,int degen,double *mod
		    ,double **v,double *modo);
double modfind(double theta);
//void find_crit(double *mod,int Nmod,int Nlenses,double Re2,double *x2);
double deflect_translated(double beta,double *mod,double *x,double *y,double *mag,int N
   ,int Nlenses,double Re2,double *x2);
void deflection_total(double *ximage,double *angle,int fulllens,double sdeg,double *ang_lens
			,int mag,double *Amag);
double minaxis(double thetaX);
double deflection_model(double beta,double *mod,double *x,double *y,double *mag,int N
		    ,int Nlenses,double Re2,double *x2);
int find_image_number(double *yo,double *x_center,double *mod,int Nmod,int Nlenses,double Re2,double *x2);
double finiteMag(double radsource,double *xo,double *mod,int Nmod,int Nlenses,double Re2,double *x2t);

void FindLensSimple(AnaLens *lens,int Nimages,Point *image_positions,double *y,double **dx_sub);
void FindLensSimple(AnaLens *lens,ImageInfo *imageinfo ,int Nimages ,double *y,double **dx_sub);

#endif
