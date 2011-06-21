/*
 * peak_refinement.c
 *
 *  Created on: Apr 8, 2011
 *      Author: bmetcalf
 *
 *      Written for project with S.Hilbert for adaptive ray tracing through the Millenium II simulation.
 *
 *      Uses a defection solver provided by the user.  Adaptively finds regions with high kappa and refines
 *      them to higher resolution
 */

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <Tree.h>
#include <Kist.h>
#include <point.h>
#include <KistDriver.h>
#include <peak_refinement.h>


BeamHndl find_peaks(double center[],double range,unsigned long Ngrid,double rEinsteinMin,double kappa_max
		,unsigned long *Nbeams,void (*rayshooter)(unsigned long N,Point *i_point)){

	/*       center[2] - center of the grid region
	 *       range  - range of the grid region in whatever units you like
	 *       Ngrid - the initial 1-D number of grid points ie grid resolution
	 *               the actual number must be a power of 2, if the input is not it will be
	 *               changed to the next higher power of 2
	 *       rEinsteinMin - the Einstein radius of the smallest lens that should be resolved, sets resolution target
	 *       kappa_max - highest kappa to be refined to, 1 or 2 is probably good enough
	 *       Nbeams - number of final grid points
	 *       rayshooter - user supplied function that does the
	 *       	         N - number of points to be shot
	 *       			 i_points - an array of points i_point[i].x[] is the position on the image plane
	 *                              i_points[i].kappa should be calculated and the source positions put into
	 *                              i_points[i].image->x,  same units as range
	 *
	 *
	 *     returns a pointer to an array of beams, beams use less memory than points
	 */

	Point **i_points,*s_points,*dummy;
	double res_target = 0,threshold;
	long Nnewpoints = 0,Ntemp;
	unsigned long i;
	ImageInfo *imageinfo;
	TreeHndl i_tree,s_tree;

	if( (Ngrid & (Ngrid-1)) != 0)
		Ngrid = pow(2, (int)( log(Ngrid)/log(2.) ) + 1 );

	//printf("Ngrid = %li\n",Ngrid);

	// set up grid
	i_points = &dummy;   // keep this on the stack

	(*i_points) = NewPointArray(Ngrid*Ngrid,True);
	xygridpoints((*i_points),range,center,Ngrid,0);
	s_points=LinkToSourcePoints((*i_points),Ngrid*Ngrid);

	rayshooter(Ngrid*Ngrid,(*i_points));

	// build trees
	i_tree=BuildTree((*i_points),Ngrid*Ngrid);
	s_tree=BuildTree(s_points,Ngrid*Ngrid);

	imageinfo=NewImageInfo(1);

	for(i = 0; i< Ngrid*Ngrid; ++i) InsertAfterCurrentKist(imageinfo->imagekist,&(*i_points)[i]);

	// increase threshold while increasing angular resolution
	for(threshold =  rEinsteinMin*Ngrid/range/4 ; threshold <= kappa_max ; threshold *= 3){

		res_target = rEinsteinMin/threshold/2;  // keeps resolution below size of smallest lens

		imageinfo->gridrange[2] = 1.0e99;
		imageinfo->gridrange[0] = imageinfo->gridrange[1]  = 0.0;

		// take out points that are no longer above threshold
		MoveToTopKist(imageinfo->imagekist);
		Ntemp = imageinfo->imagekist->Nunits;
		for(i=0;i<Ntemp;++i){
			if(getCurrentKist(imageinfo->imagekist)->kappa < threshold){

				getCurrentKist(imageinfo->imagekist)->in_image = False;
				if(AtTopKist(imageinfo->imagekist)){
					TakeOutCurrentKist(imageinfo->imagekist);
				}else{
					TakeOutCurrentKist(imageinfo->imagekist);
					MoveDownKist(imageinfo->imagekist);
				}
			}else{
				getCurrentKist(imageinfo->imagekist)->in_image = True;

				if(getCurrentKist(imageinfo->imagekist)->gridsize > imageinfo->gridrange[1])
					imageinfo->gridrange[1] = getCurrentKist(imageinfo->imagekist)->gridsize;

				if(getCurrentKist(imageinfo->imagekist)->gridsize < imageinfo->gridrange[2])
					imageinfo->gridrange[2] = getCurrentKist(imageinfo->imagekist)->gridsize;

				//printf("    gridsize = %e\n",getCurrentKist(imageinfo->imagekist)->gridsize);
				MoveDownKist(imageinfo->imagekist);
			}
		}
		assert(imageinfo->imagekist->Nunits <= Ntemp);

		//printf("restarget = %e gridrange[2] = %e  gridrange[1] = %e\n",res_target,imageinfo->gridrange[2],imageinfo->gridrange[1]);

		//printf("first pass in image %li\n",imageinfo->imagekist->Nunits);
		findborders4(i_tree,imageinfo);

		//printf("restarget = %e gridrange[2] = %e  gridrange[1] = %e\n",res_target,imageinfo->gridrange[2],imageinfo->gridrange[1]);
		Nnewpoints = -i_tree->pointlist->Npoints;
		refine_grid_kist(i_tree,s_tree,imageinfo,1,res_target,2,True,False,i_points);
		Nnewpoints += i_tree->pointlist->Npoints;

		//printf("Nnewpoints = %li\n",Nnewpoints);
		// **** do rayshooting ****
		rayshooter(Nnewpoints,(*i_points));

		while(Nnewpoints > 0){

			// add new points that are above the threshold to image
			for(i=0;i<Nnewpoints;++i){
				if((*i_points)[i].kappa > threshold){

					InsertAfterCurrentKist(imageinfo->imagekist,&(*i_points)[i]);
					MoveDownKist(imageinfo->imagekist);
					getCurrentKist(imageinfo->imagekist)->in_image = True;

					if(getCurrentKist(imageinfo->imagekist)->gridsize > imageinfo->gridrange[1])
						imageinfo->gridrange[1] = getCurrentKist(imageinfo->imagekist)->gridsize;
					if(getCurrentKist(imageinfo->imagekist)->gridsize < imageinfo->gridrange[2])
						imageinfo->gridrange[2] = getCurrentKist(imageinfo->imagekist)->gridsize;

				}else{
					(*i_points)[i].in_image = False;
				}
			}

			// find borders
			findborders4(i_tree,imageinfo);

			// refine all image points and outer border

			Nnewpoints = -i_tree->pointlist->Npoints;
			refine_grid_kist(i_tree,s_tree,imageinfo,1,res_target,2,True,False,i_points);
			Nnewpoints += i_tree->pointlist->Npoints;
			//printf("Nnewpoints = %li\n",Nnewpoints);

			// **** do rayshooting  ****
			rayshooter(Nnewpoints,(*i_points));
		}
	}

	// need some way of returning all the points

	*Nbeams = i_tree->pointlist->Npoints;
	BeamHndl beams = (Beam *)malloc(*Nbeams*sizeof(Beam));
	copyPointToBeam(i_tree->pointlist,beams);

	// free memory in grids, trees and images
	//free(imageinfo->points);
	freeImageInfo(imageinfo,1);
	freeTree(i_tree);
	freeTree(s_tree);

	return beams;
}

void copyPointToBeam(ListHndl pointlist,Beam *beam){
	unsigned long i;

	MoveToTopList(pointlist);
	i=0;
	do{
		beam[i].image[0] = pointlist->current->x[0];
		beam[i].image[1] = pointlist->current->x[1];

		beam[i].source[0] = pointlist->current->image->x[0];
		beam[i].source[1] = pointlist->current->image->x[1];

		beam[i].kappa = pointlist->current->kappa;
		beam[i].gamma[0] = pointlist->current->gamma[0];
		beam[i].gamma[1] = pointlist->current->gamma[1];

		++i;
	}while(MoveDownList(pointlist));

	return;
}

void copyBeamToPoint(Beam *beam,Point *pointarr,unsigned long N){
	unsigned long i;

	for(i=0;i<N;++i){
		pointarr[i].x[0] = beam[i].source[0];
		pointarr[i].x[1] = beam[i].source[1];

		pointarr[i].image->x[0] = beam[i].image[0];
		pointarr[i].image->x[1] = beam[i].image[1];

		pointarr[i].kappa = beam[i].kappa;
		pointarr[i].gamma[0] = beam[i].gamma[0];
		pointarr[i].gamma[1] = beam[i].gamma[1];
	}

	return;
}
