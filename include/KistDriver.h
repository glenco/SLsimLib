/*
 * KistDriver.h
 *
 *  Created on: Nov 16, 2010
 *      Author: bmetcalf
 */

#ifndef KIST_DRIVER_H_
#define KIST_DRIVER_H_

#include <Tree.h>
#include <point.h>
#include <Kist.h>
#include <image_info.h>

void FindAllBoxNeighborsKist(TreeHndl tree,Point *point,KistHndl neighbors);
void _FindAllBoxNeighborsKist(TreeHndl tree,Branch *leaf,KistHndl neighbors);
void _FindAllBoxNeighborsKist_iter(TreeHndl tree,Branch *leaf,KistHndl neighbors);
void PointsWithinEllipKist(TreeHndl tree,double *ray,float rmax,float rmin,float posangle,KistHndl neighborkist);
double PointsWithinKist(TreeHndl tree,double *ray,float rmax,KistHndl neighborkist,short markpoints);
void _PointsWithinKist(TreeHndl tree,double *ray,float *rmax,KistHndl neighborkist
		,short markpoints,double *maxgridsize);
void PointsWithinKist_iter(TreeHndl tree,double *ray,float rmin,float rmax,KistHndl neighborkist);
Point *NearestNeighborKist(TreeHndl tree,double *ray,int Nneighbors,KistHndl neighborkist);
void SetInImage(KistHndl kist,Boo value);
void DirtyFoF(ImageInfo *imageinfo ,int *Nimages ,double linkinglength ,int MaxNimages );
void _DirtyFoF(KistHndl neighbors,KistHndl wholekist,double linkinglength);
void DirtyDivider(ImageInfo *imageinfo,int *Nimages ,int MaxNimages ,int Ngroup);
void _DirtyDivider(KistHndl neighbors,KistHndl wholekist,int Ngroup);

#endif
