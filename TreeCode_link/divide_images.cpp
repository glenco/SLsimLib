/*
 * divide_images.c
 *
 *  Created on: Nov 11, 2010
 *      Author: R.B. Metcalf
 */
#include "slsimlib.h"

/** \ingroup ImageFindingL2
 * \brief finds the points that are within the circular source and divides
 * the images.
 *
 */
void find_divide_images(TreeHndl i_tree,TreeHndl s_tree
	,double *source_x,double source_r
	,ImageInfo *imageinfo,int *Nimages,int Nimagesmax){

	PointsWithinKist(s_tree,source_x,source_r,imageinfo->imagekist,0);

	// no image points found
	if(imageinfo->imagekist->Nunits() == 0){
		*Nimages = 0;
		return;
	}

	// move from source plane to image plane
	TranformPlanesKist(imageinfo->imagekist);

	if(imageinfo->imagekist->Nunits() == 1){
		*Nimages = 1;
		return;
	}

	divide_images_kist(i_tree,imageinfo,Nimages,Nimagesmax);
	return;
}

/** \ingroup functions
 *
 *	Returns the genus of an image by counting the number of disconnected outer borders.
 *	The in_image flag must be set to FALSE for all points on the grid and are returned to
 *	that value.  The outer border of the image must be found first;
 */
int ImageGenus(TreeHndl i_tree,ImageInfo *imageinfo){
	assert(imageinfo->outerborder);

	KistHndl kist = imageinfo->outerborder;
	KistHndl tmp_kist = new Kist,new_kist = new Kist;
	int number,Nimagesmax=100;
	unsigned long Ntemp,j;


	if(kist->Nunits() < 2){
		return 0;
	}

	// mark points in tree as in image and transfer points to temporary temporary new_imagekist

	kist->MoveToTop();
	do{
		getCurrentKist(kist)->in_image = TRUE;
		getCurrentKist(kist)->image->in_image = TRUE;
		InsertAfterCurrentKist(new_kist,kist->getCurrent());
		MoveDownKist(new_kist);
	}while(kist->Down());

	number=0;
	do{

		//printf("   new_imagekist %li\n",new_imagekist->Nunits());
		partition_images_kist(getCurrentKist(new_kist),tmp_kist,i_tree);

		// take out points that got un-marked in partition_images

		Ntemp = new_kist->Nunits();

		for( j=0,MoveToTopKist(new_kist) ; j < Ntemp ; ++j){
			if(getCurrentKist(new_kist)->in_image == FALSE){

			    if(AtTopKist(new_kist)){
			    	TakeOutCurrentKist(new_kist);
			    }else{
			    	TakeOutCurrentKist(new_kist);
			    	MoveDownKist(new_kist);
			    }
			}else{
				MoveDownKist(new_kist);
			}
		}

		assert(AtBottomKist(new_kist));
		++number;
	}while(new_kist->Nunits() > 0 && number < Nimagesmax);

	assert(new_kist->Nunits() == 0);
	delete new_kist;
	delete tmp_kist;

	return number-1;
}


/* \ingroup ImageFindingL2
 *  divide_images
 *
 * \brief `Reorders the image points up into separate images that are linked by cell
 * neighbors.
 *
 * Should scale like NlogN  instead of N^2.  This is achieved with a recursive
 *   algorithm.
 *
 * on entering:
 *     imageinfo->imagekist must contain all the point in all the images in any order.
 *     The flags in_image == FALSE for all points in i_tree and their images that
 *        are not in the image.  Not required that points in the image be
 *        flagged.
 * on exit:
 * 	   imagelist is reordered so that all points in an image are contiguous
 *     Nimages is updataed
 *     imageinfo[i].points are not changed
 *     imageinfo[i].Npoints is set to number of points in ith image
 *	   image point flags in_image == TRUE
 *	   the area and area_error of each image is calculated
 *
void divide_images(TreeHndl i_tree,ImageInfo *imageinfo
		,int *Nimages,int Nimagesmax){
	unsigned long i,j,Ntemp;
	KistHndl new_imagekist = new Kist;
	double tmp = 0;

	if(imageinfo->imagekist->Nunits() < 2){
		*Nimages = 1;
		//imageinfo->Npoints = imageinfo->imagekist->Nunits();

		return ;
	}
	// mark points in tree as in image

	assert(imageinfo->imagekist->top->data);
	MoveToTopKist(imageinfo->imagekist);
	do{
		getCurrentKist(imageinfo->imagekist)->in_image = TRUE;
		getCurrentKist(imageinfo->imagekist)->image->in_image = TRUE;
	}while(MoveDownKist(imageinfo->imagekist));

	assert(imageinfo->imagekist->top->data);

	i=0;
	do{
		//imageinfo[i].points = getCurrentKist(imageinfo->imagekist);

		getCurrentKist(imageinfo->imagekist)->in_image = FALSE;
		getCurrentKist(imageinfo->imagekist)->image->in_image = FALSE;
		imageinfo[i].Npoints = 1;

		assert(imageinfo->imagekist->top->data);

		partition_images(getCurrentKist(imageinfo->imagekist),&(imageinfo[i].Npoints),i_tree);
		assert(imageinfo->imagekist->top->data);

		// take out points that got un-marked in partition_images

		Ntemp = imageinfo->imagekist->Nunits();
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
	}while(imageinfo->imagekist->Nunits() > 0 && i < Nimagesmax);

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

	assert(imageinfo->imagekist->Nunits() == 0);

	new_imagekist->MoveToTop();
	while(new_imagekist->Nunits() > 0){
		imageinfo->imagekist->InsertAfterCurrent(new_imagekist->TakeOutCurrent());
		imageinfo->imagekist->Down();
	}

	delete new_imagekist;

	MoveToTopKist(imageinfo->imagekist);
	do{
		getCurrentKist(imageinfo->imagekist)->in_image = TRUE;
		getCurrentKist(imageinfo->imagekist)->image->in_image = TRUE;
	}while(MoveDownKist(imageinfo->imagekist));


	return;
}*/

/** \ingroup ImageFindingL2
 *
 *  divide_images_kist
 *
 * \brief Divides the image points up into separate images that are linked by cell
 * neighbors.
 *
 * Should scale like NlogN  instead of N^2.  This is achieved with a recursive
 *   algorithm.
 *
 * on entering:
 *     imageinfo->imagekist must contain all the point in all the images in any order.
 *     The flags in_image == FALSE for all points in i_tree and their images that
 *        are not in the image.  Not required that points in the image be
 *        flagged.
 * on exit:
 *     Nimages is updataed
 *     imageinfo[i].imagekist is set to the points in ith image
 *	   image point flags in_image == TRUE
 *	   the area and area_error of each image are calculated
 *
 *	   If more than Nimagesmax images are found imageinfo the remaining points are put into the nearest images.
 *
 *	   imageinfo[i].ShouldNotRefile = 0 is set for every image.
 *
 *	   Calculates images' geometric centers and area's.
 *
 */
void divide_images_kist(
	    TreeHndl i_tree
	    ,ImageInfo *imageinfo
	    ,int *Nimages
	    ,int Nimagesmax
	    ){
	unsigned long i,j,Ntemp,Ntest;
	KistHndl new_imagekist = new Kist;
	double tmp = 0;

	if(imageinfo->imagekist->Nunits() < 2){
		*Nimages = 1;
		return ;
	}

	// mark points in tree as in image and transfer points to temporary temporary new_imagekist

	MoveToTopKist(imageinfo->imagekist);
	Ntest = imageinfo->imagekist->Nunits();
	while(imageinfo->imagekist->Nunits() > 0){
		getCurrentKist(imageinfo->imagekist)->in_image = TRUE;
		getCurrentKist(imageinfo->imagekist)->image->in_image = TRUE;
		InsertAfterCurrentKist(new_imagekist,TakeOutCurrentKist(imageinfo->imagekist));
		MoveDownKist(new_imagekist);
	}

	i=0;
	do{

		//printf("   new_imagekist %li\n",new_imagekist->Nunits());
		imageinfo[i].area = partition_images_kist(getCurrentKist(new_imagekist),imageinfo[i].imagekist,i_tree);

/*		if(imageinfo[i].imagekist->Nunits() > new_imagekist->Nunits()){
			bool test1,test2;
			test1 = new_imagekist->AreDataUnique();
			test2 = imageinfo[i].imagekist->AreDataUnique();
			assert(imageinfo[i].imagekist->Nunits() <= new_imagekist->Nunits());  // check that no more than
		}
*/

		// take out points that got un-marked in partition_images

		Ntemp = new_imagekist->Nunits();
		imageinfo[i].area = imageinfo[i].area_error = 0.0;

		for( j=0,MoveToTopKist(new_imagekist) ; j < Ntemp ; ++j){
			if(getCurrentKist(new_imagekist)->in_image == FALSE){

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

			//printf("%li %li\n",new_imagekist->Nunits(),imageinfo[i].imagekist->Nunits());
		}
		imageinfo[i].area_error /= imageinfo[i].area;

		// calculate image centroid
		imageinfo[i].centroid[0] = 0;
		imageinfo[i].centroid[1] = 0;
		MoveToTopKist(imageinfo[i].imagekist);
		do{
			tmp = pow(getCurrentKist(imageinfo[i].imagekist)->gridsize,2 );
			imageinfo[i].centroid[0] += tmp*getCurrentKist(imageinfo[i].imagekist)->x[0];
			imageinfo[i].centroid[1] += tmp*getCurrentKist(imageinfo[i].imagekist)->x[1];
		}while(MoveDownKist(imageinfo[i].imagekist));
		imageinfo[i].centroid[0] /= imageinfo[i].area;
		imageinfo[i].centroid[1] /= imageinfo[i].area;

		assert((Ntemp - new_imagekist->Nunits() - imageinfo[i].imagekist->Nunits()) == 0);
		assert(AtBottomKist(new_imagekist));
		++i;
	}while(new_imagekist->Nunits() > 0 && i < Nimagesmax);

	*Nimages = i;


	// If there are more than NimageMax images, put the extra points into the closest image.
	if(i == Nimagesmax){
		double r2,rmin;

		while(new_imagekist->Nunits() > 0){
			for(i=0,rmin=1.0e200,j=0;i<*Nimages;++i){
				r2 = pow(getCurrentKist(new_imagekist)->x[0] - imageinfo[i].centroid[0],2)
					+ pow(getCurrentKist(new_imagekist)->x[1] - imageinfo[i].centroid[1],2);
				if(rmin > r2 ){
					rmin = r2;
					j = i;
				}
			}
			tmp = pow(getCurrentKist(new_imagekist)->gridsize,2);
			imageinfo[j].centroid[0] = ( imageinfo[j].centroid[0]*imageinfo[j].area
						+ tmp*getCurrentKist(new_imagekist)->x[0] )/(imageinfo[j].area + tmp);
			imageinfo[j].centroid[1] = ( imageinfo[j].centroid[1]*imageinfo[j].area
						+ tmp*getCurrentKist(new_imagekist)->x[1] )/(imageinfo[j].area + tmp);
			imageinfo[j].area += tmp;
			InsertAfterCurrentKist(imageinfo[j].imagekist,TakeOutCurrentKist(new_imagekist));
		}
	}

	assert(new_imagekist->Nunits() == 0);
	delete new_imagekist;

	// mark all image points
	for(i=0;i<*Nimages;++i){
		imageinfo[i].ShouldNotRefine = 0;
		MoveToTopKist(imageinfo[i].imagekist);
		do{
			getCurrentKist(imageinfo[i].imagekist)->in_image = TRUE;
			getCurrentKist(imageinfo[i].imagekist)->image->in_image = TRUE;
			--Ntest;
		}while(MoveDownKist(imageinfo[i].imagekist));
	}

	assert(Ntest == 0);
	return;
}

/** \ingroup ImageFindingL2
 *
 * \brief recursive function that un-marks all the points that are attached
 * to point by cell neighbors that were previously marked in_image==TRUE
 */

void partition_images(Point *point,unsigned long *N_in_image,TreeHndl i_tree){

	assert(point);
	assert(i_tree);

	KistHndl neighbors = new Kist;

	// find all the neighbors of point
	FindAllBoxNeighborsKist(i_tree,point,neighbors);

	// remove neighbors that are not marked as in image
	//  It reduces the memory usage to do this before recursion

	for(MoveToTopKist(neighbors) ; neighbors->Nunits() > 0 ;){
		if(getCurrentKist(neighbors)->in_image == FALSE){
			TakeOutCurrentKist(neighbors);
		}else{
			if(!MoveDownKist(neighbors)) break;
		}
	}

	while(neighbors->Nunits() > 0){
		if(getCurrentKist(neighbors)->in_image == TRUE){
			getCurrentKist(neighbors)->in_image = FALSE;
			getCurrentKist(neighbors)->image->in_image = FALSE;
			++(*N_in_image);
			partition_images(getCurrentKist(neighbors),N_in_image,i_tree);
		}
		TakeOutCurrentKist(neighbors);
	}

	//printf("exiting partition_images\n");

	assert(neighbors->Nunits() == 0);
	delete neighbors;

	return;
}

/** \ingroup ImageFindingL2
 *
 *  finds all the points in i_tree with in_image = TRUE that are connected to point
 *    by cell neighbors of cell neighbors.  The resulting kist of points
 *    is left in imagekist.  The in_image marks are NOT returned to their original
 *    values.  The ones that are put into imagekist are changed to in_image = FALSE
 */
double partition_images_kist(Point *point,KistHndl imagekist,TreeHndl i_tree){
	assert(point);
	assert(i_tree);
	assert(point->in_image == TRUE);

	double area = 0.0;
	KistHndl neighbors = new Kist;

	EmptyKist(imagekist);

	InsertAfterCurrentKist(imagekist,point);
	do{

		if(getCurrentKist(imagekist)->in_image == TRUE){
			getCurrentKist(imagekist)->in_image = FALSE;

			area += pow(getCurrentKist(imagekist)->gridsize,2);

			// find all the neighbors of point
			FindAllBoxNeighborsKist(i_tree,getCurrentKist(imagekist),neighbors);

			// remove neighbors that are not marked as in image
			for(MoveToTopKist(neighbors) ; neighbors->Nunits() > 0 ;){
				if(getCurrentKist(neighbors)->in_image == FALSE){
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

	assert(neighbors->Nunits() == 0);
	delete neighbors;

	return area;
}
