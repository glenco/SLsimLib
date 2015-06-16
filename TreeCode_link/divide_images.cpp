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
	,PosType *source_x,PosType source_r
  ,std::vector<ImageInfo> &imageinfo,int *Nimages){

	s_tree->PointsWithinKist(source_x,source_r,imageinfo[0].imagekist,0);

	// no image points found
	if(imageinfo[0].imagekist->Nunits() == 0){
		*Nimages = 0;
		return;
	}

	// move from source plane to image plane
	imageinfo[0].imagekist->TranformPlanes();

	if(imageinfo[0].imagekist->Nunits() == 1){
		*Nimages = 1;
		return;
	}

	divide_images_kist(i_tree,imageinfo,Nimages);
	return;
}

/** \ingroup functions
 *
 *	Returns the genus of an image by counting the number of disconnected outer borders.
 *	The in_image flag must be set to NO for all points on the grid and are returned to
 *	that value.  The outer border of the image must be found first;
 */
int ImageGenus(TreeHndl i_tree,ImageInfo *imageinfo){
	assert(imageinfo->outerborder);

	Kist<Point> * kist = imageinfo->outerborder;
	Kist<Point> * tmp_kist = new Kist<Point>;
	Kist<Point> * new_kist = new Kist<Point>;
	int number,Nimagesmax=100;
	unsigned long Ntemp,j;


	if(kist->Nunits() < 2){
		return 0;
	}

	// mark points in tree as in image and transfer points to temporary temporary new_imagekist

	kist->MoveToTop();
	do{
		kist->getCurrent()->in_image = YES;
		kist->getCurrent()->image->in_image = YES;
		new_kist->InsertAfterCurrent(kist->getCurrent());
		new_kist->Down();
	}while(kist->Down());

	number=0;
	do{

		//printf("   new_imagekist %li\n",new_imagekist.Nunits());
		partition_images_kist(new_kist->getCurrent(),tmp_kist,i_tree);

		// take out points that got un-marked in partition_images

		Ntemp = new_kist->Nunits();

		for( j=0,new_kist->MoveToTop() ; j < Ntemp ; ++j){
			if(new_kist->getCurrent()->in_image == NO){

			    if(new_kist->AtTop()){
			    	new_kist->TakeOutCurrent();
			    }else{
			    	new_kist->TakeOutCurrent();
			    	new_kist->Down();
			    }
			}else{
				new_kist->Down();
			}
		}

		assert(new_kist->AtBottom() || new_kist->OffBottom());
		++number;
	}while(new_kist->Nunits() > 0 && number < Nimagesmax);
	if(number == Nimagesmax)
		new_kist->Empty();


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
 *     The flags in_image == NO for all points in i_tree and their images that
 *        are not in the image.  Not required that points in the image be
 *        flagged.
 * on exit:
 * 	   imagelist is reordered so that all points in an image are contiguous
 *     Nimages is updataed
 *     imageinfo[i].points are not changed
 *     imageinfo[i].Npoints is set to number of points in ith image
 *	   image point flags in_image == YES
 *	   the area and area_error of each image is calculated
 *
void divide_images(TreeHndl i_tree,ImageInfo *imageinfo
		,int *Nimages,int Nimagesmax){
	unsigned long i,j,Ntemp;
	Kist<Point> * new_imagekist = new Kist<Point>;
	PosType tmp = 0;

	if(imageinfo->imagekist->Nunits() < 2){
		*Nimages = 1;
		//imageinfo->Npoints = imageinfo->imagekist->Nunits();

		return ;
	}
	// mark points in tree as in image

	assert(imageinfo->imagekist->top->data);
	MoveToTopKist(imageinfo->imagekist);
	do{
		getCurrentKist(imageinfo->imagekist)->in_image = YES;
		getCurrentKist(imageinfo->imagekist)->image->in_image = YES;
	}while(MoveDownKist(imageinfo->imagekist));

	assert(imageinfo->imagekist->top->data);

	i=0;
	do{
		//imageinfo[i].points = getCurrentKist(imageinfo->imagekist);

		getCurrentKist(imageinfo->imagekist)->in_image = NO;
		getCurrentKist(imageinfo->imagekist)->image->in_image = NO;
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

	new_imagekist.MoveToTop();
	while(new_imagekist.Nunits() > 0){
		imageinfo->imagekist->InsertAfterCurrent(new_imagekist->TakeOutCurrent());
		imageinfo->imagekist->Down();
	}

	MoveToTopKist(imageinfo->imagekist);
	do{
		getCurrentKist(imageinfo->imagekist)->in_image = YES;
		getCurrentKist(imageinfo->imagekist)->image->in_image = YES;
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
 *     The flags in_image == NO for all points in i_tree and their images that
 *        are not in the image.  Not required that points in the image be
 *        flagged.
 * on exit:
 *     Nimages is updataed
 *     imageinfo[i].imagekist is set to the points in ith image
 *	   image point flags in_image == YES
 *	   the area and area_error of each image are calculated
 *
 *
 *	   imageinfo[i].ShouldNotRefile = 0 is set for every image.
 *
 *	   Calculates images' geometric centers and area's.
 *
 */
void divide_images_kist(
	    TreeHndl i_tree
      ,std::vector<ImageInfo> &imageinfo
	    ,int *Nimages
	    ){
	unsigned long i,j,Ntemp,Ntest;
	Kist<Point> new_imagekist;
	PosType tmp = 0;

  assert(imageinfo.size() > 0);
	if(imageinfo[0].imagekist->Nunits() < 2){
		*Nimages = 1;
		return ;
	}

	  // mark points in tree as in image and transfer points to temporary temporary new_imagekist

	imageinfo[0].imagekist->MoveToTop();
	Ntest = imageinfo[0].imagekist->Nunits();
	while(imageinfo[0].imagekist->Nunits() > 0){
		imageinfo[0].imagekist->getCurrent()->in_image = YES;
		imageinfo[0].imagekist->getCurrent()->image->in_image = YES;
		new_imagekist.InsertAfterCurrent(imageinfo[0].imagekist->TakeOutCurrent());
		new_imagekist.Down();
	}

	i=0;
	do{

    if(i > imageinfo.size()-1) imageinfo.resize(i+5);
		//printf("   new_imagekist %li\n",new_imagekist.Nunits());
		imageinfo[i].area = partition_images_kist(new_imagekist.getCurrent(),imageinfo[i].imagekist,i_tree);

/*		if(imageinfo[i].imagekist->Nunits() > new_imagekist.Nunits()){
			bool test1,test2;
			test1 = new_imagekist.AreDataUnique();
			test2 = imageinfo[i].imagekist->AreDataUnique();
			assert(imageinfo[i].imagekist->Nunits() <= new_imagekist.Nunits());  // check that no more than
		}
*/

		// take out points that got un-marked in partition_images

		Ntemp = new_imagekist.Nunits();
		imageinfo[i].area = imageinfo[i].area_error = 0.0;

		for( j=0,new_imagekist.MoveToTop() ; j < Ntemp ; ++j){
			if(new_imagekist.getCurrent()->in_image == NO){

				// calculates area of image
				tmp = pow(new_imagekist.getCurrent()->gridsize,2 );
				imageinfo[i].area += tmp;
			    if(imageinfo[i].area_error < tmp) imageinfo[i].area_error = tmp;

			    if(new_imagekist.AtTop()){
			    	new_imagekist.TakeOutCurrent();
			    }else{
			    	new_imagekist.TakeOutCurrent();
			    	new_imagekist.Down();
			    }
			}else{
				new_imagekist.Down();
			}

			//printf("%li %li\n",new_imagekist.Nunits(),imageinfo[i].imagekist->Nunits());
		}
		imageinfo[i].area_error /= imageinfo[i].area;

		// calculate image centroid
		imageinfo[i].centroid[0] = 0;
		imageinfo[i].centroid[1] = 0;
    imageinfo[i].gridrange[0] = imageinfo[i].gridrange[2] = imageinfo[i].imagekist->getCurrent()->gridsize;
    
		imageinfo[i].imagekist->MoveToTop();
		do{
			tmp = imageinfo[i].imagekist->getCurrent()->gridsize;
			imageinfo[i].centroid[0] += tmp*tmp*imageinfo[i].imagekist->getCurrent()->x[0];
			imageinfo[i].centroid[1] += tmp*tmp*imageinfo[i].imagekist->getCurrent()->x[1];
      imageinfo[i].gridrange[0] = MAX(imageinfo[i].gridrange[0],tmp);
      imageinfo[i].gridrange[2] = MIN(imageinfo[i].gridrange[2],tmp);
		}while(imageinfo[i].imagekist->Down());
		imageinfo[i].centroid[0] /= imageinfo[i].area;
		imageinfo[i].centroid[1] /= imageinfo[i].area;
    imageinfo[i].gridrange[1] = imageinfo[i].gridrange[2];

		assert((Ntemp - new_imagekist.Nunits() - imageinfo[i].imagekist->Nunits()) == 0);
		assert(new_imagekist.OffBottom());
		++i;
	}while(new_imagekist.Nunits() > 0 );

	*Nimages = i;


	// If there are more than NimageMax images, put the extra points into the closest image.
	/*if(i == Nimagesmax){
		PosType r2,rmin;

		while(new_imagekist.Nunits() > 0){
			for(i=0,rmin=1.0e200,j=0;i<*Nimages;++i){
				r2 = pow(new_imagekist.getCurrent()->x[0] - imageinfo[i].centroid[0],2)
					+ pow(new_imagekist.getCurrent()->x[1] - imageinfo[i].centroid[1],2);
				if(rmin > r2 ){
					rmin = r2;
					j = i;
				}
			}
			tmp = pow(new_imagekist.getCurrent()->gridsize,2);
			imageinfo[j].centroid[0] = ( imageinfo[j].centroid[0]*imageinfo[j].area
						+ tmp*new_imagekist.getCurrent()->x[0] )/(imageinfo[j].area + tmp);
			imageinfo[j].centroid[1] = ( imageinfo[j].centroid[1]*imageinfo[j].area
						+ tmp*new_imagekist.getCurrent()->x[1] )/(imageinfo[j].area + tmp);
			imageinfo[j].area += tmp;
			imageinfo[j].imagekist->InsertAfterCurrent(new_imagekist.TakeOutCurrent());
		}
	}*/

	assert(new_imagekist.Nunits() == 0);
  
	// mark all image points
	for(i=0;i<*Nimages;++i){
		imageinfo[i].ShouldNotRefine = 0;
		imageinfo[i].imagekist->MoveToTop();
		do{
			imageinfo[i].imagekist->getCurrent()->in_image = YES;
			imageinfo[i].imagekist->getCurrent()->image->in_image = YES;
			--Ntest;
		}while(imageinfo[i].imagekist->Down());
	}

  imageinfo.resize(*Nimages);
  
	assert(Ntest == 0);
	return;
}

/** \ingroup ImageFindingL2
 *
 * \brief recursive function that un-marks all the points that are attached
 * to point by cell neighbors that were previously marked in_image==YES
 */

void partition_images(Point *point,unsigned long *N_in_image,TreeHndl i_tree){

	assert(point);
	assert(i_tree);

	Kist<Point> * neighbors = new Kist<Point>;

	// find all the neighbors of point
	i_tree->FindAllBoxNeighborsKist(point,neighbors);

	// remove neighbors that are not marked as in image
	//  It reduces the memory usage to do this before recursion

	for(neighbors->MoveToTop() ; neighbors->Nunits() > 0 ;){
		if(neighbors->getCurrent()->in_image == NO){
			neighbors->TakeOutCurrent();
		}else{
			if(!(neighbors->Down())) break;
		}
	}

	while(neighbors->Nunits() > 0){
		if(neighbors->getCurrent()->in_image == YES){
			neighbors->getCurrent()->in_image = NO;
			neighbors->getCurrent()->image->in_image = NO;
			++(*N_in_image);
			partition_images(neighbors->getCurrent(),N_in_image,i_tree);
		}
		neighbors->TakeOutCurrent();
	}

	//printf("exiting partition_images\n");

	assert(neighbors->Nunits() == 0);
	delete neighbors;

	return;
}

/** \ingroup ImageFindingL2
 *
 *  finds all the points in i_tree with in_image = YES that are connected to point
 *    by cell neighbors of cell neighbors.  The resulting kist of points
 *    is left in imagekist.  The in_image marks are NOT returned to their original
 *    values.  The ones that are put into imagekist are changed to in_image = NO
 */
PosType partition_images_kist(Point *point,Kist<Point> * imagekist,TreeHndl i_tree){
	assert(point);
	assert(i_tree);
	assert(point->in_image == YES);

	PosType area = 0.0;
	Kist<Point> * neighbors = new Kist<Point>;

	imagekist->Empty();

	imagekist->InsertAfterCurrent(point);
	do{

		if(imagekist->getCurrent()->in_image == YES){
			imagekist->getCurrent()->in_image = NO;

			area += pow(imagekist->getCurrent()->gridsize,2);

			// find all the neighbors of point
			i_tree->FindAllBoxNeighborsKist(imagekist->getCurrent(),neighbors);

			// remove neighbors that are not marked as in image
			for(neighbors->MoveToTop() ; neighbors->Nunits() > 0 ;){
				if(neighbors->getCurrent()->in_image == NO){
					neighbors->TakeOutCurrent();
				}else{
					imagekist->InsertAfterCurrent(neighbors->TakeOutCurrent());
				}
			}

		}else{
			imagekist->TakeOutCurrent();
		}
	}while(imagekist->Down());

	//assert(AreDataUniqueKist(imagekist));

	assert(neighbors->Nunits() == 0);
	delete neighbors;

	return area;
}
