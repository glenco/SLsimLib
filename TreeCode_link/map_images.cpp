/*
 * map_images.c
 *
 *  Created on: Oct 6, 2010
 *      Author: bmetcalf
 */
/*#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <nrutil.h>
#include <Tree.h>
#include <analytic_lens.h>
#include <Kist.h>
#include <KistDriver.h>
#include <divide_images.h>
#include <map_images.h>
#include <tree_maintenance.h>*/

#include <slsimlib.h>

const float mumin = 0.3;  // actually the sqrt of the minimum magnification
const int Ngrid_block = 3;
const float FracResTarget = 3.0e-4;
//const float FracResTarget = 1.0e-3;
const float target_all = 1.0e-2;
const float smallest = 0.1;

/** \ingroup ImageFinding
 *  \brief Find images and refine them based on their surface brightness distribution.
 *
 *  Uses find_images_kist() to initially find and refine images and then uses a surface brightness
 *  based criterion the refine the most important parts of the lens.
 *
 *  map_images is intended for mapping images of sources more complicated than simple circles.
 *
 */
void map_images(
	    /// model
		LensHndl lens,SourceBLR *source
		 /// Tree of grid points
		,GridHndl grid
		,int *Nimages /// number of images found
		,ImageInfo *imageinfo /// information on each image
		,int NimageMax   /// maximum number of images allowed
		,double initial_size   /// Initial size of source for telescoping, 0 to start from the initial grid size.
		,bool splitimages  /// TRUE each image is refined to target accuracy, otherwise all images are treated as one
		,ExitCriterion criterion  /// see data type
		,bool kappa_off  /// turns off calculation of surface density, shear, magnification and time delay
		){

	assert(lens);
	assert(grid->s_tree);
	assert(grid->i_tree);
	assert(imageinfo->imagekist);

	unsigned long Nimagepoints;
	double tmp,y[2],area_tot;
	static double oldy[2],oldr=0;
	//short moved;
	long i,j;
	//Point *i_points,*s_points;
	//Point *point;
	//time_t to;
	//ListHndl tmp_border_pointer;
	static int oldNimages=0;
	//Point **dummy_pnt;

	// flush in-source-markers and surface brightnesses
	//sublist=NewList();

	if(oldr==0) oldr=source->source_r;
	if((oldy[0]==source->source_x[0])*(oldy[1]==source->source_x[1])*(oldr > source->source_r)) initial_size=oldr;

	if(source->source_r <= 0.0){ERROR_MESSAGE(); std::printf("ERROR: find_images, point source must have a resolution target\n"); exit(1);}


	// do an initial refinement to find all images and refine grid
	// the source->source_r_in is used as a characteristic small size for the source
	assert(source->source_r_in > 0);
	find_images_kist(lens,source->source_x,source->source_r_in,grid,Nimages
			  ,imageinfo,NimageMax,&Nimagepoints,0,false,0,false,true);


	// find kist of image points and divide into images
	find_divide_images(grid->i_tree,grid->s_tree,source->source_x,source->source_r,imageinfo,Nimages,NimageMax);

	//printf("images identified\n");
	//printf("Nimages = %i\n",*Nimages);


	/////////////////////////////////////////////
	// link image points lists for each image
	// and calculate surface brightness at each point
	/////////////////////////////////////////////

	// ****** calculate surface brightnesses and flux of each image  ******
	area_tot=0;
//#pragma omp parallel for schedule(dynamic) private(i, j, y) reduction(+:area_tot)
	for(i=0 ; i < *Nimages ; ++i){
		MoveToTopKist(imageinfo[i].imagekist);
		imageinfo[i].area = 0.0;
		//printf("%li points in image %i\n",imageinfo[i].Npoints,i);
		for(j = 0 ; j < imageinfo[i].imagekist->Nunits() ; ++j,MoveDownKist(imageinfo[i].imagekist)){

			y[0] = getCurrentKist(imageinfo[i].imagekist)->image->x[0] - source->source_x[0];
			y[1] = getCurrentKist(imageinfo[i].imagekist)->image->x[1] - source->source_x[1];
			getCurrentKist(imageinfo[i].imagekist)->surface_brightness = (source->source_sb_func)(y);
			//printf("%e  %e\n",getCurrentKist(imageinfo[i].imagekist)->surface_brightness,source->source_sb_func(y));
			//assert(getCurrentKist(imageinfo[i].imagekist)->surface_brightness >= 0.0);
			if(getCurrentKist(imageinfo[i].imagekist)->surface_brightness > 0.0){
			//	printf(" %li %li sb = %e",i,j,getCurrentKist(imageinfo[i].imagekist)->surface_brightness);
				getCurrentKist(imageinfo[i].imagekist)->in_image = true;  // Initialized value for refine_grid_on_image
				imageinfo[i].area += pow(getCurrentKist(imageinfo[i].imagekist)->gridsize,2)
			                     *getCurrentKist(imageinfo[i].imagekist)->surface_brightness;
			}else{
				getCurrentKist(imageinfo[i].imagekist)->in_image = false;
			}
		}
		area_tot += imageinfo[i].area;
		//assert(imageinfo[i].area > 0);
	}
	if(area_tot == 0.0) return;

	/*
	 ******* refine images based on flux in each pixel ******
	 */
	 MoveToTopList(grid->i_tree->pointlist);
	 i=0;
	 do{ assert(grid->i_tree->pointlist->current->leaf); ++i;}while(MoveDownList(grid->i_tree->pointlist));

	i=0;
	if(area_tot != 0.0) while( refine_grid_on_image(lens,source,grid,imageinfo,*Nimages
			,FracResTarget,criterion,kappa_off) ){ ++i;

			/***************** test lines  **************************
			MoveToTopList(grid->i_tree->pointlist);
			do{ assert(grid->i_tree->pointlist->current->leaf); }while(MoveDownList(grid->i_tree->pointlist));
			/***************** **************************/
			}

	//printf("i=%i Nold=%li\n",i,Nold);
	//printf("%li\n",imagelist->Npoints);

	oldy[0]=source->source_x[0];
	oldy[1]=source->source_x[1];
	oldr=source->source_r;

	oldNimages=*Nimages;

	// find image centroid
#pragma omp parallel for schedule(dynamic) private(i, j, tmp)
	for(i=0;i<*Nimages;++i){
		MoveToTopKist(imageinfo[i].imagekist);
		tmp=0.0;
		imageinfo[i].centroid[0] = 0.0;
		imageinfo[i].centroid[1] = 0.0;
		for(j = 0 ; j < imageinfo[i].imagekist->Nunits() ; ++j,MoveDownKist(imageinfo[i].imagekist) ){
			tmp += pow(getCurrentKist(imageinfo[i].imagekist)->gridsize,2)*getCurrentKist(imageinfo[i].imagekist)->surface_brightness;
			imageinfo[i].centroid[0] += getCurrentKist(imageinfo[i].imagekist)->x[0]*pow(getCurrentKist(imageinfo[i].imagekist)->gridsize,2)
					*getCurrentKist(imageinfo[i].imagekist)->surface_brightness;
			imageinfo[i].centroid[1] += getCurrentKist(imageinfo[i].imagekist)->x[1]*pow(getCurrentKist(imageinfo[i].imagekist)->gridsize,2)
					*getCurrentKist(imageinfo[i].imagekist)->surface_brightness;
			getCurrentKist(imageinfo[i].imagekist)->in_image = false;  // re-set marks
		}
		if(imageinfo[i].Npoints > 0 ){
			imageinfo[i].centroid[0] /= tmp;
			imageinfo[i].centroid[1] /= tmp;
		}

		//printf("  %i  centroid = %e %e N = %li\n",i,imageinfo[i].centroid[0],imageinfo[i].centroid[1]
		//                                    ,imageinfo[i].Npoints);
	}
	assert(AtBottomKist(imageinfo[i].imagekist));

	return ;
}

int refine_grid_on_image(LensHndl lens,SourceHndl source,GridHndl grid,ImageInfo *imageinfo
		,unsigned long Nimages,double res_target,ExitCriterion criterion
		,bool kappa_off){

	/* refines grid according to criterion based on flux in each grid cell
	 *
	 *  points with in_image == 1 are tested for refinement
	 *
	 * criterion = TotalArea    stops refining when error in total area reaches res_target
	 * 	         = EachImage    stops refining when each image reaches error limit or is smaller than res_target
	 *           = Resolution   stops refining when grid resolution is smaller than res_target in all images
	 *           = FillHoles    does the same as EachImage but also refines the neighbors to the cells that fullfill
	 *                            the EachImage criterion
	 */

	//printf("entering refine_grid\n");

  if(Nimages < 1) return 0;

  int i,j,k,number_of_refined,count; /* Ngrid_block must be odd */
  double total_area,y[2],flux=0,tmp=0,sbmax,sbmin;
  Point *i_points,*s_points;
  unsigned long Ncells,Nold,Ntmp;

  for(i=0,total_area=0;i<Nimages;++i) total_area += imageinfo[i].area;
  for(i=0,Nold=0;i<Nimages;++i) Nold += imageinfo[i].imagekist->Nunits();
  if(total_area == 0.0) return 0;

  KistHndl nearest=new Kist;

  number_of_refined=0;
  for(i=0,Ncells=0;i<Nimages;++i){
	MoveToTopKist(imageinfo[i].imagekist);
	if(imageinfo[i].area > 0.0){
		for(j = 0 ; j < imageinfo[i].imagekist->Nunits() ; ++j,MoveDownKist(imageinfo[i].imagekist) ){

			// count number of gridcells to be refined
			// cells in image

			//if(getCurrentKist(imageinfo[i].imagekist)->in_image == 1){

  			flux = pow(getCurrentKist(imageinfo[i].imagekist)->gridsize,2)*getCurrentKist(imageinfo[i].imagekist)->surface_brightness;
  			getCurrentKist(imageinfo[i].imagekist)->in_image = 0;

  			if(criterion == TotalArea && flux > res_target*total_area){
  				getCurrentKist(imageinfo[i].imagekist)->in_image = 1;
  				++Ncells;
  			}

  			if(criterion == EachImage && flux > res_target*imageinfo[i].area
  				&& flux > target_all*res_target*total_area ){
  				getCurrentKist(imageinfo[i].imagekist)->in_image = 1;
  				++Ncells;
  			}

  			if(criterion == Resolution &&
  					getCurrentKist(imageinfo[i].imagekist)->gridsize > res_target ){
  				getCurrentKist(imageinfo[i].imagekist)->in_image = 1;
  				++Ncells;
  			}

  			if(criterion == FillHoles ){

  				if(flux > res_target*imageinfo[i].area
  						&& flux > target_all*res_target*total_area
  				){
  					getCurrentKist(imageinfo[i].imagekist)->in_image = 1;
  					++Ncells;

  				}else{  // fill the holes

  					FindAllBoxNeighborsKist(grid->i_tree,getCurrentKist(imageinfo[i].imagekist),nearest);
  					//PrintList(nearest);
  					sbmin = sbmax = 0;
  					MoveToTopKist(nearest);
  					do{
  						y[0] = getCurrentKist(nearest)->image->x[0] - source->source_x[0];
  						y[1] = getCurrentKist(nearest)->image->x[1] - source->source_x[1];
  						tmp = (source->source_sb_func)(y)*pow(getCurrentKist(nearest)->gridsize,2);

  						sbmax = MAX(sbmax,tmp);
  						//sbmin = MIN(sbmin,tmp);

  						if( tmp > res_target*imageinfo[i].area
  						 && tmp > target_all*res_target*total_area){
  							getCurrentKist(imageinfo[i].imagekist)->in_image = 1;
   							++Ncells;
  							break;
  						}
  					}while( MoveDownKist(nearest) );

  					/*
  					//if( sbmax > 0.0 && sbmax > res_target*imageinfo[i].area ){
  					if( sbmax > res_target*imageinfo[i].area
  							&& sbmax > target_all*res_target*total_area){

  						getCurrentKist(imageinfo[i].imagekist)->in_image = 1;
  						++Ncells;
  					}
	               */
  				}
  			//}
  			}
		}
	}else{
		// case of image with no flux in it
		MoveToTopKist(imageinfo[i].imagekist);
		do{ getCurrentKist(imageinfo[i].imagekist)->in_image = 0;}while(MoveDownKist(imageinfo[i].imagekist));
	}
  }

  delete nearest;
  if(Ncells == 0) return 0;

  if(Ncells > 0){

	  // make space for new points to go into tree
 	  i_points = NewPointArray((Ngrid_block*Ngrid_block-1)*Ncells,true);

	  for(i=0,Ncells=0 ; i < Nimages ; ++i){
		  Ntmp = imageinfo[i].imagekist->Nunits();
	 	  MoveToTopKist(imageinfo[i].imagekist);
		  for(j = 0,Nold = 0; j < Ntmp ; ++j){

			  if(getCurrentKist(imageinfo[i].imagekist)->in_image){

				  ++count;
				  xygridpoints(&i_points[(Ngrid_block*Ngrid_block-1)*Ncells]
				      ,getCurrentKist(imageinfo[i].imagekist)->gridsize*(Ngrid_block-1)/Ngrid_block
				      ,getCurrentKist(imageinfo[i].imagekist)->x,Ngrid_block,1);

				  // reduces size of actual gridpoint in the tree here
				  getCurrentKist(imageinfo[i].imagekist)->gridsize /= Ngrid_block;

				  // link new points into list
				  for(k=0;k < Ngrid_block*Ngrid_block-1 ;++k){

					  // put point into image imageinfo[i].imagekist
					  InsertAfterCurrentKist(imageinfo[i].imagekist,&(i_points[(Ngrid_block*Ngrid_block-1)*Ncells+k]));
					  MoveDownKist(imageinfo[i].imagekist);

					  // mark for future calculation of surface brightness
					  getCurrentKist(imageinfo[i].imagekist)->surface_brightness = -1;
					  // mark for future possible refinement
					  getCurrentKist(imageinfo[i].imagekist)->in_image = 1;

					  ++number_of_refined;
				  }

				  ++Ncells;

				  Nold += Ngrid_block*Ngrid_block-1;
				  // add points to inner edge list

			  }
			  MoveDownKist(imageinfo[i].imagekist);
		  }
		  imageinfo[i].Npoints += Nold;
	  }

      s_points = LinkToSourcePoints(i_points,(Ngrid_block*Ngrid_block-1)*Ncells);

      // lens the added points
      lens->rayshooterInternal((Ngrid_block*Ngrid_block-1)*Ncells,i_points,kappa_off);

		/***************** test lines  **************************
 	 do{ assert(grid->i_tree->pointlist->current->leaf);}while(MoveDownList(grid->i_tree->pointlist));
  	 moveTop(grid->i_tree);
  	 do{
  		 if(grid->i_tree->current->number == 4829221){
  			 printf("STOP!\n");
  		 }
  		 assert(grid->i_tree->current->points->next || grid->i_tree->current->points->prev);
  	 }while(TreeWalkStep(grid->i_tree,true));
		/*****************  **************************/

      // add points to trees
      AddPointsToTree(grid->i_tree,i_points,Ncells*(Ngrid_block*Ngrid_block-1));

		/***************** test lines  **************************
  	 do{ assert(grid->i_tree->pointlist->current->leaf);}while(MoveDownList(grid->i_tree->pointlist));
  	 moveTop(grid->i_tree);
  	 do{ assert(grid->i_tree->current->points->next || grid->i_tree->current->points->prev); }while(TreeWalkStep(grid->i_tree,true));
		/***************** **************************/

      AddPointsToTree(grid->s_tree,s_points,Ncells*(Ngrid_block*Ngrid_block-1));

      // Go through list of image points and copy the information
      //  just calculated into the new points and re-calculate the area/flux
      //  for each image
     for(i=0,Ncells=0;i<Nimages;++i){
    	 imageinfo[i].area = 0.0;
    	 MoveToTopKist(imageinfo[i].imagekist);
    	 for(j = 0 ; j < imageinfo[i].imagekist->Nunits() ; ++j,MoveDownKist(imageinfo[i].imagekist) ){
    		 if(getCurrentKist(imageinfo[i].imagekist)->surface_brightness == -1){
    			 //PointCopyData(getCurrentKist(imageinfo[i].imagekist),getCurrentKist(imageinfo[i].imagekist)->image);
    			 y[0] = getCurrentKist(imageinfo[i].imagekist)->image->x[0] - source->source_x[0];
    			 y[1] = getCurrentKist(imageinfo[i].imagekist)->image->x[1] - source->source_x[1];
    			 getCurrentKist(imageinfo[i].imagekist)->surface_brightness = (source->source_sb_func)(y);
    		 }
    		 imageinfo[i].area += pow(getCurrentKist(imageinfo[i].imagekist)->gridsize,2)
    				 * getCurrentKist(imageinfo[i].imagekist)->surface_brightness;
    	 }
    	 assert(AtBottomKist(imageinfo[i].imagekist));
     }
  }

  return number_of_refined;
}
