/*
 * divide_images.c
 *
 *  Created on: Nov 11, 2010
 *      Author: R.B. Metcalf
 */
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <Tree.h>
#include <Kist.h>
#include <KistDriver.h>
#include <divide_images.h>

void find_divide_images(TreeHndl i_tree,TreeHndl s_tree
		,double *source_x,double source_r
		,ImageInfo *imageinfo,int *Nimages,int Nimagesmax){
	/* finds the points that are within the circular source and divides
	 * the images.
	 *
	 */

	PointsWithinKist(s_tree,source_x,source_r,imageinfo->imagekist,0);

	// no image points found
	if(imageinfo->imagekist->Nunits == 0){
		*Nimages = 0;
		return;
	}

	// move from source plane to image plane
	TranformPlanesKist(imageinfo->imagekist);

	if(imageinfo->imagekist->Nunits == 1){
		*Nimages = 1;
		imageinfo->Npoints = 1;
		return;
	}

	divide_images_kist(i_tree,imageinfo,Nimages,Nimagesmax);
	return;
}

void divide_images(TreeHndl i_tree,ImageInfo *imageinfo
		,int *Nimages,int Nimagesmax){
	/* divide_images
	 *
	 * Divides the image points up into separate images that are linked by cell
	 * neighbors.
	 *
	 * Should scale like NlogN  instead of N^2.  This is achieved with a recursive
	 *   algorithm.
	 *
	 * on entering:
	 *     imageinfo->imagekist must contain all the point in all the images in any order.
	 *     The flags in_image == False for all points in i_tree and their images that
	 *        are not in the image.  Not required that points in the image be
	 *        flagged.
	 * on exit:
	 * 	   imagelist is reordered so that all points in an image are contiguous
	 *     Nimages is updataed
	 *     imageinfo[i].points are not changed
	 *     imageinfo[i].Npoints is set to number of points in ith image
	 *	   image point flags in_image == True
	 *	   the area and area_error of each image is calculated
	 */
	unsigned long i,j,Ntemp;
	KistHndl new_imagekist = NewKist();
	double tmp = 0;

	if(imageinfo->imagekist->Nunits < 2){
		*Nimages = 1;
		imageinfo->Npoints = imageinfo->imagekist->Nunits;

		return ;
	}
	// mark points in tree as in image

	assert(imageinfo->imagekist->top->data);
	MoveToTopKist(imageinfo->imagekist);
	do{
		getCurrentKist(imageinfo->imagekist)->in_image = True;
		getCurrentKist(imageinfo->imagekist)->image->in_image = True;
	}while(MoveDownKist(imageinfo->imagekist));

	assert(imageinfo->imagekist->top->data);

/*	MoveToTopList(imagelist);
	do{
		imagelist->current->image->in_image = True;
		imagelist->current->image->image->in_image = True;
	}while(MoveDownList(imagelist));
*/

	i=0;
	do{
		//imageinfo[i].points = getCurrentKist(imageinfo->imagekist);

		getCurrentKist(imageinfo->imagekist)->in_image = False;
		getCurrentKist(imageinfo->imagekist)->image->in_image = False;
		imageinfo[i].Npoints = 1;

		assert(imageinfo->imagekist->top->data);
		partition_images(getCurrentKist(imageinfo->imagekist),&(imageinfo[i].Npoints),i_tree);
		assert(imageinfo->imagekist->top->data);

		// take out points that got un-marked in partition_images

		Ntemp = imageinfo->imagekist->Nunits;
		imageinfo[i].area = imageinfo[i].area_error = 0.0;
		for( j=0,MoveToTopKist(imageinfo->imagekist) ; j < Ntemp ; ++j,MoveDownKist(imageinfo->imagekist) ){
			if(getCurrentKist(imageinfo->imagekist)->in_image == 0){

				// calculates area of image
				tmp = pow(getCurrentKist(imageinfo->imagekist)->gridsize,2 );
				imageinfo[i].area += tmp;
			    if(imageinfo[i].area_error < tmp) imageinfo[i].area_error = tmp;

				InsertAfterCurrentKist(new_imagekist,TakeOutCurrentKist(imageinfo->imagekist));
				MoveDownKist(new_imagekist);
			}
		}
		imageinfo[i].area_error /= imageinfo[i].area;

		++i;
	}while(imageinfo->imagekist->Nunits > 0 && i < Nimagesmax);

	*Nimages = i;

	// if more than Nimagemax images add points to last image
	if(i == Nimagesmax){
		MoveToTopKist(imageinfo->imagekist);
		do{
			InsertAfterCurrentKist(new_imagekist,TakeOutCurrentKist(imageinfo->imagekist));
			MoveDownKist(new_imagekist);
			++(imageinfo[i-1].Npoints);
		}while(MoveDownKist(imageinfo->imagekist));
	}

	assert(imageinfo->imagekist->Nunits == 0);

	imageinfo->imagekist->Nunits = new_imagekist->Nunits;
	imageinfo->imagekist->top = new_imagekist->top;
	imageinfo->imagekist->bottom = new_imagekist->bottom;
	imageinfo->imagekist->current = new_imagekist->current;

	free(new_imagekist);


	MoveToTopKist(imageinfo->imagekist);
	do{
		getCurrentKist(imageinfo->imagekist)->in_image = True;
		getCurrentKist(imageinfo->imagekist)->image->in_image = True;
	}while(MoveDownKist(imageinfo->imagekist));


	return;
}
void divide_images_kist(TreeHndl i_tree,ImageInfo *imageinfo,int *Nimages,int Nimagesmax){
	/* divide_images
	 *
	 * Divides the image points up into separate images that are linked by cell
	 * neighbors.
	 *
	 * Should scale like NlogN  instead of N^2.  This is achieved with a recursive
	 *   algorithm.
	 *
	 * on entering:
	 *     imageinfo->imagekist must contain all the point in all the images in any order.
	 *     The flags in_image == False for all points in i_tree and their images that
	 *        are not in the image.  Not required that points in the image be
	 *        flagged.
	 * on exit:
	 * 	   imagelist[i] contains all points in that image
	 *     Nimages is updataed
	 *     imageinfo[i].points are not changed
	 *     imageinfo[i].Npoints is set to number of points in ith image
	 *	   image point flags in_image == True
	 *	   the area and area_error of each image are calculated
	 */
	unsigned long i,j,Ntemp,Ntest;
	KistHndl new_imagekist = NewKist();
	double tmp = 0;

	if(imageinfo->imagekist->Nunits < 2){
		*Nimages = 1;
		imageinfo->Npoints = imageinfo->imagekist->Nunits;

		return ;
	}

	// mark points in tree as in image and transfer points to temporary kist
	assert(imageinfo->imagekist->top->data);
	MoveToTopKist(imageinfo->imagekist);
	Ntest = imageinfo->imagekist->Nunits;
	while(imageinfo->imagekist->Nunits > 0){
		getCurrentKist(imageinfo->imagekist)->in_image = True;
		getCurrentKist(imageinfo->imagekist)->image->in_image = True;
		InsertAfterCurrentKist(new_imagekist,TakeOutCurrentKist(imageinfo->imagekist));
		MoveDownKist(new_imagekist);
	}

	i=0;
	do{

		//printf("   new_imagekist %li\n",new_imagekist->Nunits);
		imageinfo[i].area = partition_images2(getCurrentKist(new_imagekist),imageinfo[i].imagekist,i_tree);
		assert(imageinfo[i].imagekist->Nunits <= new_imagekist->Nunits);  // check that no more than

		imageinfo[i].Npoints = imageinfo[i].imagekist->Nunits;
		// take out points that got un-marked in partition_images

		Ntemp = new_imagekist->Nunits;
		imageinfo[i].area = imageinfo[i].area_error = 0.0;

		for( j=0,MoveToTopKist(new_imagekist) ; j < Ntemp ; ++j){
			if(getCurrentKist(new_imagekist)->in_image == False){

				// calculates area of image
				tmp = pow(getCurrentKist(new_imagekist)->gridsize,2 );
				imageinfo[i].area += tmp;
			    if(imageinfo[i].area_error < tmp) imageinfo[i].area_error = tmp;

			    if(AtTopKist(new_imagekist)){
			    	TakeOutCurrentKist(new_imagekist);
			    }else{
			    	TakeOutCurrentKist(new_imagekist);
			    	MoveDownKist(new_imagekist);
			    }
			}else{
				MoveDownKist(new_imagekist);
			}

			//printf("%li %li\n",new_imagekist->Nunits,imageinfo[i].imagekist->Nunits);
		}
		imageinfo[i].area_error /= imageinfo[i].area;

		assert((Ntemp - new_imagekist->Nunits - imageinfo[i].imagekist->Nunits) == 0);
		assert(AtBottomKist(new_imagekist));
		++i;
	}while(new_imagekist->Nunits > 0 && i < Nimagesmax);

	*Nimages = i;

	// if more than Nimagemax images add points to last image
	if(i == Nimagesmax){
		MoveToBottomKist(imageinfo[i-1].imagekist);
		do{
			InsertAfterCurrentKist(imageinfo[i-1].imagekist,TakeOutCurrentKist(new_imagekist));
			MoveDownKist(imageinfo[i-1].imagekist);
			++(imageinfo[i-1].Npoints);
		}while(new_imagekist->Nunits > 0);
	}

	assert(new_imagekist->Nunits == 0);
	freeKist(new_imagekist);

	// mark all image points
	for(i=0;i<*Nimages;++i){
		MoveToTopKist(imageinfo[i].imagekist);
		do{
			getCurrentKist(imageinfo[i].imagekist)->in_image = True;
			getCurrentKist(imageinfo[i].imagekist)->image->in_image = True;
			--Ntest;
		}while(MoveDownKist(imageinfo[i].imagekist));
	}

	assert(Ntest == 0);
	return;
}

void partition_images(Point *point,unsigned long *N_in_image,TreeHndl i_tree){
	/* recursive function that un-marks all the points that are attached
	 * to point by cell neighbors that were previously marked in_image==True
	 */

	assert(point);
	assert(i_tree);

	KistHndl neighbors = NewKist();

	// find all the neighbors of point
	FindAllBoxNeighborsKist(i_tree,point,neighbors);

	// remove neighbors that are not marked as in image
	//  It reduces the memory usage to do this before recursion

	for(MoveToTopKist(neighbors) ; neighbors->Nunits > 0 ;){
		if(getCurrentKist(neighbors)->in_image == False){
			TakeOutCurrentKist(neighbors);
		}else{
			if(!MoveDownKist(neighbors)) break;
		}
	}

	while(neighbors->Nunits > 0){
		if(getCurrentKist(neighbors)->in_image == True){
			getCurrentKist(neighbors)->in_image = False;
			getCurrentKist(neighbors)->image->in_image = False;
			++(*N_in_image);
			partition_images(getCurrentKist(neighbors),N_in_image,i_tree);
		}
		TakeOutCurrentKist(neighbors);
	}

	//printf("exiting partition_images\n");

	assert(neighbors->Nunits == 0);
	free(neighbors);

	return;
}

double partition_images2(Point *point,KistHndl imagekist,TreeHndl i_tree){
/* finds all the points with in_image = True that are connected to point
 *    by cell neighbors of cell neighbors.  The resulting kist of points
 *    is left in imagekist.  The in_image marks are NOT returned to their original
 *    values.  The ones that are put into imagekist are changed to in_inage = False
 */
	assert(point);
	assert(i_tree);
	assert(point->in_image == True);

	double area = 0.0;
	KistHndl neighbors = NewKist();

	EmptyKist(imagekist);

	InsertAfterCurrentKist(imagekist,point);
	do{

		if(getCurrentKist(imagekist)->in_image == True){
			getCurrentKist(imagekist)->in_image = False;

			area += pow(getCurrentKist(imagekist)->gridsize,2);

			// find all the neighbors of point
			FindAllBoxNeighborsKist(i_tree,getCurrentKist(imagekist),neighbors);

			// remove neighbors that are not marked as in image
			for(MoveToTopKist(neighbors) ; neighbors->Nunits > 0 ;){
				if(getCurrentKist(neighbors)->in_image == False){
					TakeOutCurrentKist(neighbors);
				}else{
					InsertAfterCurrentKist(imagekist,TakeOutCurrentKist(neighbors));
				}
			}

		}else{
			TakeOutCurrentKist(imagekist);
		}
	}while(MoveDownKist(imagekist));

	//assert(AreDataUniqueKist(imagekist));

	assert(neighbors->Nunits == 0);
	freeKist(neighbors);

	return area;
}
