/*
 * KistDriver.h
 *
 *  Created on: Nov 16, 2010
 *      Author: bmetcalf
 */

void FindAllBoxNeighborsKist(TreeHndl tree,Point *point,KistHndl neighbors);
void _FindAllBoxNeighborsKist(TreeHndl tree,Branch *leaf,KistHndl neighbors);
void _FindAllBoxNeighborsKist_iter(TreeHndl tree,Branch *leaf,KistHndl neighbors);
void PointsWithinEllipKist(TreeHndl tree,double *ray,float rmax,float rmin,float posangle,KistHndl neighborkist);
void PointsWithinKist(TreeHndl tree,double *ray,float rmax,KistHndl neighborkist,short markpoints);
void _PointsWithinKist(TreeHndl tree,double *ray,float *rmax,KistHndl neighborkist
		,short markpoints);
Point *NearestNeighborKist(TreeHndl tree,double *ray,int Nneighbors,KistHndl neighborkist);
