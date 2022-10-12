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

/** Finding
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
    ,std::vector<ImageInfo> &imageinfo /// the image
		,int *Nimages		/// number of peaks
		){
	//Point **i_points,*s_points,*dummy;
	PosType res_target = 0,threshold;
	long Nnewpoints = 0,Ntemp;
	unsigned long i;
	//ImageInfo *imageinfo = new ImageInfo;
	Kist<Point> * newpointskist = new Kist<Point>;
      
      if(grid->getInitRange() != grid->getNumberOfPoints() ){
        *grid = grid->ReInitialize(lens);
      }

	// Add all points to imageinfo
	//MoveToTopList(grid->i_tree->pointlist);
      PointList::iterator i_tree_pl_current;
      i_tree_pl_current.current = (grid->i_tree->pointlist.Top());
	do{
		imageinfo[0].imagekist->InsertAfterCurrent(*i_tree_pl_current);
	}while(--i_tree_pl_current);


	// increase threshold while increasing angular resolution
	for(threshold =  rEinsteinMin*(grid->getInitNgrid())/(grid->getInitRange())/4 ; threshold <= kappa_max ; threshold *= 3){
		cout << "threshold " << threshold << endl;
		res_target = rEinsteinMin/threshold/2;  // keeps resolution below size of smallest lens

		imageinfo[0].gridrange[2] = 1.0e99;
		imageinfo[0].gridrange[0] = imageinfo[0].gridrange[1]  = 0.0;

		// take out points that are no longer above threshold
		imageinfo[0].imagekist->MoveToTop();
		Ntemp = imageinfo[0].imagekist->Nunits();
		for(i=0;i<Ntemp;++i){
			if(imageinfo[0].imagekist->getCurrent()->kappa() < threshold){

				imageinfo[0].imagekist->getCurrent()->in_image = NO;
				if(imageinfo[0].imagekist->AtTop()){
					imageinfo[0].imagekist->TakeOutCurrent();
				}else{
					imageinfo[0].imagekist->TakeOutCurrent();
					imageinfo[0].imagekist->Down();
				}
			}else{
				imageinfo[0].imagekist->getCurrent()->in_image = YES;

				if(imageinfo[0].imagekist->getCurrent()->gridsize > imageinfo[0].gridrange[1])
					imageinfo[0].gridrange[1] = imageinfo[0].imagekist->getCurrent()->gridsize;

				if(imageinfo[0].imagekist->getCurrent()->gridsize < imageinfo[0].gridrange[2])
					imageinfo[0].gridrange[2] = imageinfo[0].imagekist->getCurrent()->gridsize;

				//printf("    gridsize = %e\n",getCurrentKist(imageinfo[0].imagekist)->gridsize);
				imageinfo[0].imagekist->Down();
			}
		}
		assert(imageinfo[0].imagekist->Nunits() <= Ntemp);

		//printf("restarget = %e gridrange[2] = %e  gridrange[1] = %e\n",res_target,imageinfo[0].gridrange[2],imageinfo[0].gridrange[1]);

		//printf("first pass in image %li\n",imageinfo[0].imagekist->Nunits);
    bool touches_edge;
		findborders4(grid->i_tree,imageinfo.data(),touches_edge);

		//printf("restarget = %e gridrange[2] = %e  gridrange[1] = %e\n",res_target,imageinfo[0].gridrange[2],imageinfo[0].gridrange[1]);

		Nnewpoints = IF_routines::refine_grid_kist(lens,grid,imageinfo.data(),1,res_target,2,newpointskist);


		while(newpointskist->Nunits() > 0){

			// add new points that are above the threshold to image

			while(newpointskist->Nunits() > 0){
				if(newpointskist->getCurrent()->kappa() > threshold){

					imageinfo[0].imagekist->InsertAfterCurrent(newpointskist->TakeOutCurrent());
					imageinfo[0].imagekist->Down();
					imageinfo[0].imagekist->getCurrent()->in_image = YES;

					if(imageinfo[0].imagekist->getCurrent()->gridsize > imageinfo[0].gridrange[1])
						imageinfo[0].gridrange[1] = imageinfo[0].imagekist->getCurrent()->gridsize;
					if(imageinfo[0].imagekist->getCurrent()->gridsize < imageinfo[0].gridrange[2])
						imageinfo[0].gridrange[2] = imageinfo[0].imagekist->getCurrent()->gridsize;

				}else{
					newpointskist->TakeOutCurrent()->in_image = NO;
				}
			}

			// find borders
			findborders4(grid->i_tree,imageinfo.data(),touches_edge);

			// refine all image points and outer border

			Nnewpoints = IF_routines::refine_grid_kist(lens,grid,imageinfo.data(),1,res_target,2,newpointskist);
			//printf("Nnewpoints = %li\n",Nnewpoints);

		}
	}

	divide_images_kist(grid->i_tree,imageinfo,Nimages);

	// need some way of returning all the points

	// free memory in grids, trees and images
	//free(imageinfo->points);
	//delete imageinfo;
	delete newpointskist;

	return 1;
}
} // end of namespace ImageFinding
