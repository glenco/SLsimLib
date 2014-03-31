/*
 * peak_refinement.c
 *
 *  Created on: Apr 8, 2011
 *      Author: bmetcalf
 *
 *      Written for project with S.Hilbert for adaptive ray tracing through the Millennium II simulation.
 *
 *      Uses a defection solver provided by the user.  Adaptively finds regions with high kappa and refines
 *      them to higher resolution
 */


#include "slsimlib.h"

using namespace std;

namespace ImageFinding{

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
		,PosType rEinsteinMin  /// the Einstein radius of the smallest lens that should be resolved, sets resolution target
		,PosType kappa_max     /// highest kappa to be refined to, 1 or 2 is probably good enough
		,ImageInfo *imageinfo /// the image
		,int *Nimages		/// number of peaks
		){


	//Point **i_points,*s_points,*dummy;
	PosType res_target = 0,threshold;
	long Nnewpoints = 0,Ntemp;
	unsigned long i;
	//ImageInfo *imageinfo = new ImageInfo;
	Kist<Point> * newpointskist = new Kist<Point>;

	if(grid->getInitRange() != grid->getNumberOfPoints() ) grid->ReInitializeGrid(lens);

	// Add all points to imageinfo
	MoveToTopList(grid->i_tree->pointlist);
	do{
		imageinfo->imagekist->InsertAfterCurrent(grid->i_tree->pointlist->current);
	}while(MoveDownList(grid->i_tree->pointlist));


	// increase threshold while increasing angular resolution
	for(threshold =  rEinsteinMin*(grid->getInitNgrid())/(grid->getInitRange())/4 ; threshold <= kappa_max ; threshold *= 3){
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

		Nnewpoints = ImageFinding::refine_grid_kist(lens,grid,imageinfo,1,res_target,2,newpointskist);


		while(newpointskist->Nunits() > 0){

			// add new points that are above the threshold to image

			while(newpointskist->Nunits() > 0){
				if(newpointskist->getCurrent()->kappa > threshold){

					imageinfo->imagekist->InsertAfterCurrent(newpointskist->TakeOutCurrent());
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

			Nnewpoints = ImageFinding::refine_grid_kist(lens,grid,imageinfo,1,res_target,2,newpointskist);
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
} // end of namespace ImageFinding
