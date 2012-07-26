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


#include <slsimlib.h>

using namespace std;

/** \ingroup ImageFinding
 *
 *  \brief Refines the grid based on the convergence so that high density regions have high resolution.
 *
 *  No source is used in this process.  The code acts as if all the mass is in SIS halos iteratively increasing the
 *  resolution and kappa threshhold until the desired resolution is found.
 *
 */
short find_peaks(
		LensHndl lens         /// Lens model
		,GridHndl grid        /// Grid to be refined.  It must be initialized.
		,double rEinsteinMin  /// the Einstein radius of the smallest lens that should be resolved, sets resolution target
		,double kappa_max     /// highest kappa to be refined to, 1 or 2 is probably good enough
		,ImageInfo *imageinfo /// the image
		,int *Nimages		/// number of peaks
		){


	//Point **i_points,*s_points,*dummy;
	double res_target = 0,threshold;
	long Nnewpoints = 0,Ntemp;
	unsigned long i;
	//ImageInfo *imageinfo = new ImageInfo;
	KistHndl newpointskist = new Kist;

	if(grid->getInitRange() != grid->NumberOfPoints() ) grid->ReInitializeGrid(lens);

	// Add all points to imageinfo
	MoveToTopList(grid->i_tree->pointlist);
	do{
		imageinfo->imagekist->InsertAfterCurrent(grid->i_tree->pointlist->current);
	}while(MoveDownList(grid->i_tree->pointlist));


	// increase threshold while increasing angular resolution
	for(threshold =  rEinsteinMin*(grid->getNgrid())/(grid->getInitRange())/4 ; threshold <= kappa_max ; threshold *= 3){
		cout << "threshold " << threshold << endl;
		res_target = rEinsteinMin/threshold/2;  // keeps resolution below size of smallest lens

		imageinfo->gridrange[2] = 1.0e99;
		imageinfo->gridrange[0] = imageinfo->gridrange[1]  = 0.0;

		// take out points that are no longer above threshold
		imageinfo->imagekist->MoveToTop();
		Ntemp = imageinfo->imagekist->Nunits();
		for(i=0;i<Ntemp;++i){
			if(imageinfo->imagekist->getCurrent()->kappa < threshold){

				imageinfo->imagekist->getCurrent()->in_image = FALSE;
				if(imageinfo->imagekist->AtTop()){
					imageinfo->imagekist->TakeOutCurrent();
				}else{
					imageinfo->imagekist->TakeOutCurrent();
					imageinfo->imagekist->Down();
				}
			}else{
				imageinfo->imagekist->getCurrent()->in_image = TRUE;

				if(imageinfo->imagekist->getCurrent()->gridsize > imageinfo->gridrange[1])
					imageinfo->gridrange[1] = imageinfo->imagekist->getCurrent()->gridsize;

				if(imageinfo->imagekist->getCurrent()->gridsize < imageinfo->gridrange[2])
					imageinfo->gridrange[2] = imageinfo->imagekist->getCurrent()->gridsize;

				//printf("    gridsize = %e\n",getCurrentKist(imageinfo->imagekist)->gridsize);
				imageinfo->imagekist->Down();
			}
		}
		assert(imageinfo->imagekist->Nunits() <= Ntemp);

		//printf("restarget = %e gridrange[2] = %e  gridrange[1] = %e\n",res_target,imageinfo->gridrange[2],imageinfo->gridrange[1]);

		//printf("first pass in image %li\n",imageinfo->imagekist->Nunits);
		findborders4(grid->i_tree,imageinfo);

		//printf("restarget = %e gridrange[2] = %e  gridrange[1] = %e\n",res_target,imageinfo->gridrange[2],imageinfo->gridrange[1]);

		Nnewpoints = refine_grid_kist(lens,grid,imageinfo,1,res_target,2,true,newpointskist);


		while(newpointskist->Nunits() > 0){

			// add new points that are above the threshold to image

			while(newpointskist->Nunits() > 0){
				if(newpointskist->getCurrent()->kappa > threshold){

					InsertAfterCurrentKist(imageinfo->imagekist,newpointskist->TakeOutCurrent());
					imageinfo->imagekist->Down();
					imageinfo->imagekist->getCurrent()->in_image = TRUE;

					if(imageinfo->imagekist->getCurrent()->gridsize > imageinfo->gridrange[1])
						imageinfo->gridrange[1] = imageinfo->imagekist->getCurrent()->gridsize;
					if(imageinfo->imagekist->getCurrent()->gridsize < imageinfo->gridrange[2])
						imageinfo->gridrange[2] = imageinfo->imagekist->getCurrent()->gridsize;

				}else{
					newpointskist->TakeOutCurrent()->in_image = FALSE;
				}
			}

			// find borders
			findborders4(grid->i_tree,imageinfo);

			// refine all image points and outer border

			Nnewpoints = refine_grid_kist(lens,grid,imageinfo,1,res_target,2,true,newpointskist);
			//printf("Nnewpoints = %li\n",Nnewpoints);

		}
	}

	int NimageMax = 1000;

	divide_images_kist(grid->i_tree,imageinfo,Nimages,NimageMax);

	// need some way of returning all the points

	// free memory in grids, trees and images
	//free(imageinfo->points);
	//delete imageinfo;
	delete newpointskist;

	return 1;
}

/** \ingroup ImageFinding
 *
 *  \brief Refines the grid based on the flux from implanted sources.
 *
 *
 */
// TODO Add function declaration to header file.  Maybe move this function to another file.
short refine_on_implanted_source(
		MultiLens &lens        /// MultiLens model
		,GridHndl grid         /// Grid to be refined.  It must be initialized.
		,double *theta         /// position on the sky
		,double radius         /// size of region to look for image
		,double res_target     /// the final grid resolution that is required
		,ImageInfo *imageinfo  /// the output image
		,int *Nimages		   /// number of images
		,bool kappa_off        /// turn off convergence and shear calculation to save time
		){

	long Nnewpoints = 0;
	KistHndl newpointskist = new Kist;

	imageinfo->imagekist->Empty();
	imageinfo->gridrange[2] = 1.0e99;
	imageinfo->gridrange[0] = imageinfo->gridrange[1]  = 0.0;

	PointsWithinKist(grid->i_tree,theta,radius,newpointskist,0);

	newpointskist->MoveToTop();
	do{
		// re-shoot rays to add in surface brightness from implanted sources
		lens.rayshooterInternal(1,imageinfo->imagekist->getCurrent(),kappa_off);

		if(newpointskist->getCurrent()->surface_brightness != 0){
			imageinfo->imagekist->InsertAfterCurrent(newpointskist->getCurrent());
			newpointskist->getCurrent()->in_image = TRUE;

			imageinfo->imagekist->Down();
			if(imageinfo->imagekist->getCurrent()->gridsize > imageinfo->gridrange[1])
				imageinfo->gridrange[1] = imageinfo->imagekist->getCurrent()->gridsize;

			if(imageinfo->imagekist->getCurrent()->gridsize < imageinfo->gridrange[2])
				imageinfo->gridrange[2] = imageinfo->imagekist->getCurrent()->gridsize;

		}
	}while(newpointskist->Down());

	// if there are no points
	if(newpointskist->Nunits() == 0){
		NearestNeighborKist(grid->i_tree,theta,8,imageinfo->imagekist);
	}

	// TODO Deal with case where the image is not found initially so that multiple refinements can be made before reaching the scale.

	newpointskist->Empty();

	findborders4(grid->i_tree,imageinfo);

	// refine grid to wanted resolution
	Nnewpoints = 1;
	while(Nnewpoints){
		Nnewpoints = refine_grid_kist(&lens,grid,imageinfo,1,res_target,2,true,newpointskist);
		if(Nnewpoints > 0){
			newpointskist->MoveToTop();
			do{
				if(newpointskist->getCurrent()->surface_brightness != 0){
					imageinfo->imagekist->InsertAfterCurrent(newpointskist->getCurrent());
					newpointskist->getCurrent()->in_image = TRUE;

					imageinfo->imagekist->Down();
					if(imageinfo->imagekist->getCurrent()->gridsize > imageinfo->gridrange[1])
						imageinfo->gridrange[1] = imageinfo->imagekist->getCurrent()->gridsize;

					if(imageinfo->imagekist->getCurrent()->gridsize < imageinfo->gridrange[2])
						imageinfo->gridrange[2] = imageinfo->imagekist->getCurrent()->gridsize;
				}
			}while(newpointskist->Down());

			findborders4(grid->i_tree,imageinfo);
		}
	}

	imageinfo->imagekist->MoveToTop();
	do{
		imageinfo->imagekist->getCurrent()->in_image = FALSE;
	}while(imageinfo->imagekist->Down());

	delete newpointskist;

	return 1;
}
