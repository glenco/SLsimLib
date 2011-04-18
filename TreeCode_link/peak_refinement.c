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
#include <Tree.h>
#include <Kist.h>
#include <point.h>
#include <KistDriver.h>


void find_peaks(double center[],double range,unsigned long Ngrid,void (*rayshooter)(unsigned long N,Point *i_points)){
	Point *i_points = 0,*s_points;
	double res_target = 0,threshold;
	long Nnewpoints = 0;
	unsigned long i;
	const float rEinsteinMin = 1.0;
	ImageInfo *imageinfo;
	TreeHndl i_tree,s_tree;

	// set up grid
	i_points = NewPointArray(Ngrid*Ngrid,True);
	xygridpoints(i_points,range,center,Ngrid,0);
	s_points=LinkToSourcePoints(i_points,Ngrid*Ngrid);

	rayshooter(Ngrid*Ngrid,i_points);

	// build trees
	i_tree=BuildTree(i_points,Ngrid*Ngrid);
	s_tree=BuildTree(s_points,Ngrid*Ngrid);

	imageinfo=NewImageInfo(1);

	// increase threshold while increasing angular resolution
	for(threshold = 0.01 ; threshold < 1 ; threshold += 0.1){

		res_target = rEinsteinMin/threshold/2;  // keeps resolution below size of smallest lens

		imageinfo->gridrange[2] = 1.0e99;
		imageinfo->gridrange[0] = imageinfo->gridrange[1]  = 0.0;

		// take out points that are no longer above threshold
		MoveToTopKist(imageinfo->imagekist);
		do{
			if(getCurrentKist(imageinfo->imagekist)->kappa < threshold){
				getCurrentKist(imageinfo->imagekist)->in_image = False;
				TakeOutCurrentKist(imageinfo->imagekist);
			}else{
				getCurrentKist(imageinfo->imagekist)->in_image = True;

				if(getCurrentKist(imageinfo->imagekist)->gridsize < imageinfo->gridrange[1])
					imageinfo->gridrange[1] = getCurrentKist(imageinfo->imagekist)->gridsize;
				if(getCurrentKist(imageinfo->imagekist)->gridsize > imageinfo->gridrange[2])
					imageinfo->gridrange[2] = getCurrentKist(imageinfo->imagekist)->gridsize;
			}
		}while(MoveDownKist(imageinfo->imagekist));
		findborders4(i_tree,imageinfo);

		Nnewpoints = -i_tree->pointlist->Npoints;
		refine_grid2(i_tree,s_tree,imageinfo,0,res_target,2,True,False,i_points);
		Nnewpionts += i_tree->pointlist->Npoints;
		// **** do rayshooting ****
		rayshooter(Nnewpoints,i_points);

		while(Nnewpoints > 0){

			// add new points that are above the threshold to image
			for(i=0;i<Nnewpoints;++i){
				if(i_points[i].kappa > threshold){
					InsertAfterCurrentKist(imageinfo->imagekist,&i_points[i]);
					MoveDownKist(imageinfo->imagekist);
					getCurrentKist(imageinfo->imagekist)->in_image = True;

					if(getCurrentKist(imageinfo->imagekist)->gridsize < imageinfo->gridrange[1])
						imageinfo->gridrange[1] = getCurrentKist(imageinfo->imagekist)->gridsize;
					if(getCurrentKist(imageinfo->imagekist)->gridsize > imageinfo->gridrange[2])
						imageinfo->gridrange[2] = getCurrentKist(imageinfo->imagekist)->gridsize;

				}else{
					i_points[i].in_image = False;
				}
			}

			// find borders
			findborders4(i_tree,imageinfo);

			// refine all image points and outer border

			Nnewpoints = -i_tree->pointlist->Npoints;
			refine_grid2(i_tree,s_tree,imageinfo,0,res_target,2,True,False,i_points);
			Nnewpionts += i_tree->pointlist->Npoints;
			// **** do rayshooting  ****
			rayshooter(Nnewpoints,i_points);
		}
	}

	// need some way of returning all the points

	i_points = NewPointArray(i_tree->pointlist->Npoints,True);

	MoveToTopList(i_tree->pointlist);
	for(i=0 ; i<i_tree->pointlist->Npoints ; ++i){

		output[i].image_x[0] = i_tree->pointlist->current->x[0];
		output[i].image_x[1] = i_tree->pointlist->current->x[1];

		output[i].source_x[0] = i_tree->pointlist->current->image->x[0];
		output[i].source_x[1] = i_tree->pointlist->current->image->x[1];

		output[i].kappa = i_tree->pointlist->current->kappa;
		output[i].gamma[0] = i_tree->pointlist->current->gamma[0];
		output[i].gamma[1] = i_tree->pointlist->current->gamma[1];

		MoveDownList(i_tree->pointlist);
	}


	// free memory in grids, trees and images
	free(imageinfo->points);
	freeImageInfo(imageinfo,1);
	freeTree(i_tree);
	freeTree(s_tree);

	return;
}
