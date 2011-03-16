/*
 * map_images.c
 *
 *  Created on: Oct 6, 2010
 *      Author: bmetcalf
 */
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <nrutil.h>
#include "Tree.h"
#include "../AnalyticNSIE/analytic_lens.h"
#include "../Kist/Kist.h"
#include "KistDriver.h"
#include "divide_images.h"
#include "map_images.h"

const float mumin = 0.3;  // actually the sqrt of the minimum magnification
const int Ngrid_block = 3;
const float FracResTarget = 3.0e-4;
//const float FracResTarget = 1.0e-3;
const float target_all = 1.0e-2;
const float smallest = 0.1;

void map_images(AnaLens *lens,TreeHndl s_tree,TreeHndl i_tree
		,int *Nimages,ImageInfo *imageinfo,int NimageMax,double initial_size
		  ,Boolean splitimages,ExitCriterion criterion,Boolean kappa_off){

	/*
	 * map_images is intended for mapping images of sources more complicated than simple circles
	 *
	 * returns finite refined images
	 *  it starts with a large source and reduces down to the right size refining at each step
	 *  should not miss any image larger than ~ munin*r_source linear size
	 *  if r_souce < 0 the source is a point source and fabs(r_source) is the resolution
	 *  if splitimage==TRUE each image is refined to target accuracy, otherwise all images are treated as one
	 *
	 *
	 *  kappa_off = True - turns off calculation of kappa and gamma
	 *              False - turns on calculation of kappa and gamma
	 */

	assert(lens);
	assert(s_tree);
	assert(i_tree);
	assert(imageinfo->imagekist);

	long Nsizes;
	unsigned long Nimagepoints;
	double rtemp,tmp,y[2];
	static double oldy[2],oldr=0;
	short moved;
	int i,j;
	//Point *i_points,*s_points;
	//Point *point;
	time_t to;
	//ListHndl tmp_border_pointer;
	static int oldNimages=0;

	// flush in-source-markers and surface brightnesses
	//sublist=NewList();

	if(oldr==0) oldr=lens->source_r;
	if((oldy[0]==lens->source_x[0])*(oldy[1]==lens->source_x[1])*(oldr > lens->source_r)) initial_size=oldr;

	if(lens->source_r <= 0.0){ERROR_MESSAGE(); printf("ERROR: find_images, point source must have a resolution target\n"); exit(1);}

	EmptyKist(imageinfo->imagekist);
	// reset in_image for points in tree incase they haven't
	ClearAllMarkes(i_tree);
	ClearAllMarkes(s_tree);

	// starting with a larger source size make sure all the grid sizes are small enough to find it
	time(&to);

	//if(initial_size == 0) initial_size=initialgridsize;
    Nsizes=(int)(log(initial_size/sqrt(pi)/fabs(smallest*lens->source_r))/log(Ngrid_block) ) + 1.0; // round up


    //////////////////////////////////////////
    // telescope source size down to target
    //////////////////////////////////////////
	//if(!( (oldy[0]==lens->source_x[0])*(oldy[1]==lens->source_x[1])*(oldr < lens->source_r) )){

		//rtelescope = lens->source_r;
    for(rtemp = initial_size*pow(Ngrid_block,Nsizes)
    		;rtemp >= 0.999*Ngrid_block*fabs(smallest*lens->source_r)
    		;rtemp /= Ngrid_block ){

    	do{
    		moved=image_finder(lens->source_x,rtemp,s_tree,i_tree
    				,Nimages,imageinfo,NimageMax,&Nimagepoints,-1,0);

    		assert(imageinfo->Npoints > 1);
    		assert(*Nimages);
     	}while(refine_grid(i_tree,s_tree,imageinfo,*Nimages,rtemp*mumin/Ngrid_block,2,kappa_off));

    	assert(imageinfo->Npoints > 1);
    	assert(*Nimages);

    	//printf("rtemp = %e = %e x r_source\n",rtemp,rtemp/lens->source_r);
    }

    // To conserve space, erase image point information.  It will be re-found later
    free(imageinfo->points);
    imageinfo->points = NULL;

	//////////////////////////////////////////////////////////////////////////////////
	// target source size has been reached, do full image decomposition and refinement
	/////////////////////////////////////////////////////////////////////////////////

	printf("telescoped\n");

	// do an initial uniform refinement to make sure there are enough point in
	//  the images
/*	do{
		time(&t3);
*/

	//splitter(imageinfo,300,Nimages);

			// mark image points in tree
//		PointsWithin(s_tree,lens->source_x,lens->source_r,imagelist,1);
//		moved=image_finder(lens->source_x,fabs(lens->source_r),s_tree,i_tree
//				,Nimages,imageinfo,&Nimagepoints,0,1);
/*
		if(*Nimages < 1) printf("  Nimages=%i i=%i\n",*Nimages,i);

		time(&now);
		++i;
	}while( refine_grid(i_tree,s_tree,imageinfo,*Nimages,1.0e-1,1,kappa_off)
				|| moved );
*/
//	PointsWithin(s_tree,lens->source_x,lens->source_r,imagelist,-1);

	assert(*Nimages);
	assert(imageinfo->Npoints);

	// find kist of image points and divide into images
	//printf("finding images\n");
	find_divide_images(i_tree,s_tree,lens->source_x,lens->source_r
			,imageinfo,Nimages,NimageMax);

	//printf("images identified\n");
	//printf("Nimages = %i\n",*Nimages);


	/////////////////////////////////////////////
	// link image points lists for each image
	// and calculate surface brightness at each point
	/////////////////////////////////////////////

/*	for(i=0,Nold=0;i<(*Nimages);++i) Nold += imageinfo[i].Npoints;
	for(i=0;i<Nold;++i){
		// make a copy of image point
		//InsertPointAfterCurrent(imagelist,imageinfo[0].points[i].x,imageinfo[0].points[i].id
		//		,imageinfo[0].points[i].image);
		InsertPointAfterCurrent(imagelist,&(imageinfo[0].points[i]));
		MoveDownList(imagelist);
		//PointCopyData(imagelist->current,&(imageinfo[0].points[i]));
		imagelist->current->in_image = 1;  // Initialized value for refine_grid_on_image
	}
*/

	// calculate surface brightnesses and flux of each image
	MoveToTopKist(imageinfo->imagekist);
	for(i=0 ; i < *Nimages ; ++i){
		imageinfo[i].area = 0.0;
		//printf("%li points in image %i\n",imageinfo[i].Npoints,i);
		for(j = 0 ; j < imageinfo[i].Npoints ; ++j,MoveDownKist(imageinfo->imagekist)){

			y[0] = getCurrentKist(imageinfo->imagekist)->image->x[0] - lens->source_x[0];
			y[1] = getCurrentKist(imageinfo->imagekist)->image->x[1] - lens->source_x[1];
			getCurrentKist(imageinfo->imagekist)->surface_brightness = lens->source_sb_func(y);
			getCurrentKist(imageinfo->imagekist)->in_image = 1;  // Initialized value for refine_grid_on_image

			imageinfo[i].area += pow(getCurrentKist(imageinfo->imagekist)->gridsize,2)
			                     *getCurrentKist(imageinfo->imagekist)->surface_brightness;
		}
	}

	//PrintList(imagelist);

	//printf("Nimages = %i  Npoints_total = %li\n",*Nimages,imagelist->Npoints);
	//for(i=0;i<(*Nimages);++i) printf("    Npoints = %li\n",imageinfo[i].Npoints);
	// refine images based on flux in each pixel
	i=0;
	while( refine_grid_on_image(lens,i_tree,s_tree,imageinfo,*Nimages
			,FracResTarget,criterion,kappa_off) ) ++i;

//	while( refine_grid_on_image(lens,i_tree,s_tree,imageinfo,imagelist,*Nimages,1.0e-9,criterion,kappa_off) ) ++i;

	//printf("i=%i Nold=%li\n",i,Nold);
	//printf("%li\n",imagelist->Npoints);

	oldy[0]=lens->source_x[0];
	oldy[1]=lens->source_x[1];
	oldr=lens->source_r;

	oldNimages=*Nimages;

	// find image centroid
	MoveToTopKist(imageinfo->imagekist);
	for(i=0;i<*Nimages;++i){
		tmp=0.0;
		imageinfo[i].centroid[0] = 0.0;
		imageinfo[i].centroid[1] = 0.0;
		for(j = 0 ; j < imageinfo[i].Npoints ; ++j,MoveDownKist(imageinfo->imagekist) ){
			tmp += pow(getCurrentKist(imageinfo->imagekist)->gridsize,2)*getCurrentKist(imageinfo->imagekist)->surface_brightness;
			imageinfo[i].centroid[0] += getCurrentKist(imageinfo->imagekist)->x[0]*pow(getCurrentKist(imageinfo->imagekist)->gridsize,2)
					*getCurrentKist(imageinfo->imagekist)->surface_brightness;
			imageinfo[i].centroid[1] += getCurrentKist(imageinfo->imagekist)->x[1]*pow(getCurrentKist(imageinfo->imagekist)->gridsize,2)
					*getCurrentKist(imageinfo->imagekist)->surface_brightness;
		}
		if(imageinfo[i].Npoints > 0 ){
			imageinfo[i].centroid[0] /= tmp;
			imageinfo[i].centroid[1] /= tmp;
		}

		//printf("  %i  centroid = %e %e N = %li\n",i,imageinfo[i].centroid[0],imageinfo[i].centroid[1]
		//                                    ,imageinfo[i].Npoints);
	}
	assert(AtBottomKist(imageinfo->imagekist));

	return ;
}

int refine_grid_on_image(AnaLens *lens,TreeHndl i_tree,TreeHndl s_tree,ImageInfo *imageinfo
		,unsigned long Nimages,double res_target,ExitCriterion criterion
		,Boolean kappa_off){

	/* refines grid according to criterion based on flux in each grid cell
	 *
	 *  points with in_image == 1 are tested for refinement
	 *
	 * criterion = 0 stops refining when error in total area reaches res_target
	 * 	         = 1 stops refining when each image reaches error limit or is smaller than res_target
	 *           = 2 stops refining when grid resolution is smaller than res_target in all images
	 */

	//printf("entering refine_grid\n");

  if(Nimages < 1) return 0;

  int i,j,k,number_of_refined,count; /* Ngrid_block must be odd */
  double total_area,y[2],flux=0,tmp=0,sbmax,sbmin;
  Point *i_points,*s_points;
  unsigned long Ncells,Nold;
  KistHndl nearest=NewKist();

  for(i=0,total_area=0;i<Nimages;++i) total_area +=
		  imageinfo[i].area;
  for(i=0,Nold=0;i<Nimages;++i) Nold += imageinfo[i].Npoints;
  assert(total_area > 0.0);

  number_of_refined=0;
  MoveToTopKist(imageinfo->imagekist);
  for(i=0,Ncells=0;i<Nimages;++i){
	assert(imageinfo[i].area > 0.0);
  	for(j = 0 ; j < imageinfo[i].Npoints ; ++j,MoveDownKist(imageinfo->imagekist) ){

    // count number of gridcells to be refined
    // cells in image

  		//if(getCurrentKist(imageinfo->imagekist)->in_image == 1){

  			flux = pow(getCurrentKist(imageinfo->imagekist)->gridsize,2)*getCurrentKist(imageinfo->imagekist)->surface_brightness;
  			getCurrentKist(imageinfo->imagekist)->in_image = 0;

  			if(criterion == TotalArea && flux > res_target*total_area){
  				getCurrentKist(imageinfo->imagekist)->in_image = 1;
  				++Ncells;
  			}

  			if(criterion == EachImage && flux > res_target*imageinfo[i].area
  				&& flux > target_all*res_target*total_area ){
  				getCurrentKist(imageinfo->imagekist)->in_image = 1;
  				++Ncells;
  			}

  			if(criterion == Resolution &&
  					getCurrentKist(imageinfo->imagekist)->gridsize > res_target ){
  				getCurrentKist(imageinfo->imagekist)->in_image = 1;
  				++Ncells;
  			}

  			if(criterion == FillHoles ){

  				if(flux > res_target*imageinfo[i].area
  						&& flux > target_all*res_target*total_area
  				){
  					getCurrentKist(imageinfo->imagekist)->in_image = 1;
  					++Ncells;

  				}else{  // fill the holes

  					FindAllBoxNeighborsKist(i_tree,getCurrentKist(imageinfo->imagekist),nearest);
  					//PrintList(nearest);
  					sbmin = sbmax = 0;
  					MoveToTopKist(nearest);
  					do{
  						y[0] = getCurrentKist(nearest)->image->x[0] - lens->source_x[0];
  						y[1] = getCurrentKist(nearest)->image->x[1] - lens->source_x[1];
  						tmp = lens->source_sb_func(y)*pow(getCurrentKist(nearest)->gridsize,2);

  						sbmax = MAX(sbmax,tmp);
  						//sbmin = MIN(sbmin,tmp);

  						if( tmp > res_target*imageinfo[i].area
  						 && tmp > target_all*res_target*total_area){
  							getCurrentKist(imageinfo->imagekist)->in_image = 1;
   							++Ncells;
  							break;
  						}
  					}while( MoveDownKist(nearest) );

  					/*
  					//if( sbmax > 0.0 && sbmax > res_target*imageinfo[i].area ){
  					if( sbmax > res_target*imageinfo[i].area
  							&& sbmax > target_all*res_target*total_area){

  						getCurrentKist(imageinfo->imagekist)->in_image = 1;
  						++Ncells;
  					}
	               */
  				}
  			//}
  		}
  	}
  }

  freeKist(nearest);
  if(Ncells == 0) return 0;

  if(Ncells > 0){

	  // make space for new points to go into tree
 	  i_points = NewPointArray((Ngrid_block*Ngrid_block-1)*Ncells,True);

 	  MoveToTopKist(imageinfo->imagekist);
	  for(i=0,Ncells=0;i<Nimages;++i){
		  for(j = 0,Nold = 0; j < imageinfo[i].Npoints ; ++j){

			  if(getCurrentKist(imageinfo->imagekist)->in_image){

				  ++count;
				  xygridpoints(&i_points[(Ngrid_block*Ngrid_block-1)*Ncells]
				      ,getCurrentKist(imageinfo->imagekist)->gridsize*(Ngrid_block-1)/Ngrid_block
				      ,getCurrentKist(imageinfo->imagekist)->x,Ngrid_block,1);

				  // reduces size of actual gridpoint in the tree here
				  getCurrentKist(imageinfo->imagekist)->gridsize /= Ngrid_block;

				  // link new points into list
				  for(k=0;k < Ngrid_block*Ngrid_block-1 ;++k){

					  // put point into image imageinfo->imagekist
					  InsertAfterCurrentKist(imageinfo->imagekist,&(i_points[(Ngrid_block*Ngrid_block-1)*Ncells+k]));
					  MoveDownKist(imageinfo->imagekist);

					  // mark for future calculation of surface brightness
					  getCurrentKist(imageinfo->imagekist)->surface_brightness = -1;
					  // mark for future possible refinement
					  getCurrentKist(imageinfo->imagekist)->in_image = 1;

					  ++number_of_refined;
				  }

				  ++Ncells;

				  Nold += Ngrid_block*Ngrid_block-1;
				  // add points to inner edge list

			  }
			  MoveDownKist(imageinfo->imagekist);
		  }
		  imageinfo[i].Npoints += Nold;
	  }

      s_points=LinkToSourcePoints(i_points,(Ngrid_block*Ngrid_block-1)*Ncells);

      // lens the added points
      rayshooterInternal((Ngrid_block*Ngrid_block-1)*Ncells,i_points,i_tree,kappa_off);

      // add points to trees
      AddPointsToTree(i_tree,i_points,Ncells*(Ngrid_block*Ngrid_block-1));
      AddPointsToTree(s_tree,s_points,Ncells*(Ngrid_block*Ngrid_block-1));

      // Go through list of image points and copy the information
      //  just calculated into the new points and re-calculate the area/flux
      //  for each image
     MoveToTopKist(imageinfo->imagekist);
     for(i=0,Ncells=0;i<Nimages;++i){
    	 imageinfo[i].area = 0.0;
    	 for(j = 0 ; j < imageinfo[i].Npoints ; ++j,MoveDownKist(imageinfo->imagekist) ){
    		 if(getCurrentKist(imageinfo->imagekist)->surface_brightness == -1){
    			 //PointCopyData(getCurrentKist(imageinfo->imagekist),getCurrentKist(imageinfo->imagekist)->image);
    			 y[0] = getCurrentKist(imageinfo->imagekist)->image->x[0] - lens->source_x[0];
    			 y[1] = getCurrentKist(imageinfo->imagekist)->image->x[1] - lens->source_x[1];
    			 getCurrentKist(imageinfo->imagekist)->surface_brightness = lens->source_sb_func(y);
    		 }
    		 imageinfo[i].area += pow(getCurrentKist(imageinfo->imagekist)->gridsize,2)
    				 *getCurrentKist(imageinfo->imagekist)->surface_brightness;
    	 }
     }
     assert(AtBottomKist(imageinfo->imagekist));
  }

  return number_of_refined;
}
