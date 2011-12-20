/*
 * pixelize.c
 *
 *  Created on: Feb 27, 2010
 *      Author: R.B. Metcalf
 */
/*#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <omp.h>
#include <nr.h>
#include <nrutil.h>
#include <Kist.h>
#include <Tree.h>
#include <tree_maintenance.h>
#include "image_processing.h"*/

#include <slsimlib.h>

double area=0;

/** \ingroup Image
 * \brief Takes images and pixelizes the flux into regular pixels.
 *
 * The routine constructs a tree structure of the pixel points for fast calculation of location
 * of each point in the images.  This tree is stored internally and cannot be freed.
 */
void pixelize(
		double *map    /// Output map in one dimensional array. It is always square. map[0...Npixels*Npixels-1]
		,long Npixels   /// Number of pixels in one dimension of map.
		,double range   /// One dimensional range of map in whatever units the point positions are in (generally Mpc on the lens plane.)
		,double *center /// The location of the center of the map
		,ImageInfo *imageinfo  /// An array of ImageInfo-s.  There is no reason to separate images for this routine
		,int Nimages           /// Number of images on input.
		,bool constant_sb  /// true - all images will have surface brightness = 1,
		                      /// false - surface brightness is taken from surface_brighness in  the image points
		,bool cleanmap     ///  true - erases previous pixel map, false - adds new flux to map
		){

	long ix;
	double sb=1.0,resolution=0;
	unsigned long i,ii;
	static unsigned long count=0;
	Point *points;
	static TreeHndl ptree;

	std::printf("%d %g %g %g\n", Npixels, range, center[0], center[1]);

	if( (Npixels & (Npixels-1)) != 0){
		ERROR_MESSAGE();
		std::printf("ERROR: pixelsize, Npixels is not a power of 2\n");
		exit(1);
	}

	++count;
	resolution=range/(Npixels-1);

	// initialize pixel tree
	if(count==1){
		points=NewPointArray(Npixels*Npixels,true);
		xygridpoints(points,range,center,Npixels,false);
		ptree=BuildTree(points,Npixels*Npixels);
	}

	if(cleanmap){

		if(count > 1){
			// rebuild pixel tree
			emptyTree(ptree);
			points=NewPointArray(Npixels*Npixels,true);
			xygridpoints(points,range,center,Npixels,0);
			FillTree(ptree,points,Npixels*Npixels);
			//ptree=BuildTree(points,Npixels*Npixels);
		}

		MoveToTopList(ptree->pointlist);
		for(i=0 ; i < ptree->pointlist->Npoints ; ++i){
			ptree->pointlist->current->surface_brightness = 0.0;
			MoveDownList(ptree->pointlist);
		}
		for(i=0 ; i < Npixels*Npixels ; ++i) map[i]=0.0;
	}


	sb = 1;
	for(ii=0;ii<Nimages;++ii){
		MoveToTopKist(imageinfo[ii].imagekist);

		do{

			if(!constant_sb) sb = getCurrentKist(imageinfo[ii].imagekist)->surface_brightness;

			assert(getCurrentKist(imageinfo[ii].imagekist)->leaf);
			moveTop(ptree);
			_SplitFluxIntoPixels(ptree,getCurrentKist(imageinfo[ii].imagekist)->leaf,&sb);

/*
			getCurrentKist(imageinfo[ii].imagekist)->gridsize = getCurrentKist(imageinfo[ii].imagekist)->leaf->boundary_p2[0]
		                             - getCurrentKist(imageinfo[ii].imagekist)->leaf->boundary_p1[0];
			if( boxinbox(getCurrentKist(imageinfo[ii].imagekist)->leaf,ptree->top) ){


				moveTop(ptree);
				_FindBox(ptree,getCurrentKist(imageinfo[ii].imagekist)->x);
				//ptree->current->points->surface_brightness += sb;// *pow(getCurrentKist(imageinfo[ii].imagekist)->gridsize,2);

				if(boxinbox(getCurrentKist(imageinfo[ii].imagekist)->leaf,ptree->current) == true){
					// entire cell is inside pixel
					//ptree->current->points->surface_brightness += sb*pow(getCurrentKist(imageinfo[ii].imagekist)->gridsize/resolution,2);
					ptree->current->points->surface_brightness = sb;

				}else{  // image cell is in multiple pixels

					moveUp(ptree);

					// find a super-pixel that contains all of the cell
					while(boxinbox(getCurrentKist(imageinfo[ii].imagekist)->leaf,ptree->current) == false) moveUp(ptree);

					assert((ptree->current->boundary_p1[0]-getCurrentKist(imageinfo[ii].imagekist)->leaf->boundary_p1[0]) < 0);
					assert((ptree->current->boundary_p1[1]-getCurrentKist(imageinfo[ii].imagekist)->leaf->boundary_p1[1]) < 0);
					assert((ptree->current->boundary_p2[0]-getCurrentKist(imageinfo[ii].imagekist)->leaf->boundary_p2[0]) > 0);
					assert((ptree->current->boundary_p2[1]-getCurrentKist(imageinfo[ii].imagekist)->leaf->boundary_p2[1]) > 0);

					_SplitFluxIntoPixels(ptree,getCurrentKist(imageinfo[ii].imagekist)->leaf,&sb);
				}
				// TODO fix above so that flux is distributed in pixels correctly
				//ptree->current->points->surface_brightness = getCurrentKist(imageinfo[ii].imagekist)->surface_brightness;
			}
*/

		}while(MoveDownKist(imageinfo[ii].imagekist));
	}

	MoveToTopList(ptree->pointlist);
	for(i=0;i<ptree->pointlist->Npoints;++i){
		if(ptree->pointlist->current->surface_brightness > 0.0){
			ix = IndexFromPosition(ptree->pointlist->current->x,Npixels,range,center);
			if(ix > -1) map[ix] =  ptree->pointlist->current->surface_brightness/resolution/resolution;
		}
		MoveDownList(ptree->pointlist);
	}

	return;
}


/** \ingroup LowLevel
 *  Recursively determine which pixels a leaf intersects with
 *  and add flux to that pixel.  Used in pixelizer().
 */
void _SplitFluxIntoPixels(TreeHndl ptree,Branch *leaf,double *leaf_sb){

	area = BoxIntersection(leaf,ptree->current);

	if(area > 0.0){
		if(ptree->current->child1 == NULL && ptree->current->child2 == NULL){
			assert(ptree->current->npoints == 1);
			ptree->current->points->surface_brightness += (*leaf_sb)*area;
			//ptree->current->points->surface_brightness = *leaf_sb;

			return;
		}

		if(ptree->current->child1 != NULL){
			moveToChild(ptree,1);
			_SplitFluxIntoPixels(ptree,leaf,leaf_sb);
			moveUp(ptree);
		}

		if(ptree->current->child2 != NULL){
			moveToChild(ptree,2);
			_SplitFluxIntoPixels(ptree,leaf,leaf_sb);
			moveUp(ptree);
		}
	}

	return;
}
/** \ingroup Image
 *
 * \brief Smoothes a map with a Gaussian kernel. Needs to be tested.
 */
void smoothmap(double *map_out,double *map_in,long Npixels,double range,double sigma){
	double sum=0,**mask;
	unsigned long i=0;
	long j,k,ix,iy;
	int Nmask;

	Nmask=2*(int)(3*sigma*Npixels/range + 1);
	if(Nmask < 4 ) std::printf("WARNING: pixels are large compare to psf Nmask=%i\n",Nmask);

	// set up mask
	mask=dmatrix(0,Nmask,0,Nmask);
	for(j=-Nmask/2,sum=0;j<=Nmask/2;++j){
		for(k=-Nmask/2;k<=Nmask/2;++k){
			mask[j+Nmask/2][k+Nmask/2]= exp(-(pow(j*range/(Npixels-1),2)
					                        + pow(k*range/(Npixels-1),2))/2/pow(sigma,2) );
			sum+=mask[j+Nmask/2][k+Nmask/2];
		}
	}
	for(j=-Nmask/2;j<=Nmask/2;++j) for(k=-Nmask/2;k<=Nmask/2;++k) mask[j+Nmask/2][k+Nmask/2]/=sum;

	for(i=0;i<Npixels*Npixels;++i){
		for(j=0;j<=Nmask;++j){
			ix=i%Npixels + j-Nmask/2;
			if( (ix>-1)*(ix<Npixels) ){
				for(k=0;k<=Nmask;++k){
					iy=i/Npixels + k-Nmask/2;
					if( (iy>-1)*(iy<Npixels) ){

						map_out[ix+Npixels*iy] += mask[ix][iy]*map_in[i];
					}
				}
			}
		}
	}

	free_dmatrix(mask,0,Nmask,0,Nmask);

	return ;
}
