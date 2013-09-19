

#include "slsimlib.h"

static const int NpointsRequired = 50;  // number of points required to be within an image
static const int Ngrid_block = 3;       // each cell is divided into Ngrid_block^2 subcells

//static const float mumin = 0.1;  // actually the sqrt of the minimum magnification
static const float mumin = 0.3;  // actually the sqrt of the minimum magnification
//static const float mumin = 0.5;
//static const float FracResTarget = 4.0e-4;
static const float FracResTarget = 1.0e-4;

double initialgridsize=0;


/**
 * This routine has been replaced by refine_grid_kist().  It is kept only for compatibility
 * with find_crit().
 *
 *  criterion = 0 stops refining when error in total area reaches res_target
 * 	         = 1 stops refining when each image reaches error limit or is smaller than res_target
 *           = 2 stops refining when grid resolution is smaller than res_target in all images
 */
//int refine_grid(LensHndl lens,TreeHndl i_tree,TreeHndl s_tree,OldImageInfo *imageinfo
//		,unsigned long Nimages,double res_target,short criterion,bool kappa_off){
int refine_grid(LensHndl lens,GridHndl grid,OldImageInfo *imageinfo
			,unsigned long Nimages,double res_target,short criterion,bool kappa_off,bool batch){


	//printf("entering refine_grid\n");

  if(Nimages < 1) return 0;

  int i,j,number_of_refined,count; /* Ngrid_block must be odd */
  double rmax,total_area;
  Point* point;
  short pass=0;
  long Ncells,Ncells_o;
  std::vector<Point*> points_to_refine;

  for(i=0,total_area=0;i<Nimages;++i) total_area += imageinfo[i].area;

  number_of_refined=0; 
  for(i=0;i<Nimages;++i){
    count=0;

    if(criterion == 0) pass = imageinfo[i].area*imageinfo[i].area_error/total_area > res_target;
    if(criterion == 1) pass = (imageinfo[i].area_error > res_target)*(imageinfo[i].area > 1.0e-5*total_area);
    if(criterion == 2) pass = imageinfo[i].gridrange[1] > res_target;

    // make sure no border point has a lower res than any image point

    //printf("imageinfo[%i].area_error=%e N=%i\n",i,imageinfo[i].area_error,imageinfo[i].Npoints);
    if( pass || imageinfo[i].gridrange[0]>1.01*imageinfo[i].gridrange[1]){

    	/* largest gridsize in image */
    	rmax=MAX(imageinfo[i].gridrange[1],imageinfo[i].gridrange[0]);

      // count number of gridcells to be refined
      // cells in image
    	for(j=0,Ncells=0;j<imageinfo[i].Npoints;++j) if( imageinfo[i].points[j].gridsize
    			> 1.01*rmax/Ngrid_block) ++Ncells;

    	// cells on border
		MoveToTopKist(imageinfo[i].outerborder);
    	for(j=0;j<imageinfo[i].outerborder->Nunits();++j){
    		if( getCurrentKist(imageinfo[i].outerborder)->gridsize > 1.01*rmax/Ngrid_block ){
    			// border point is marked to prevent refining more than once
    			//   it will be unmarked by the end of refine grid
    			getCurrentKist(imageinfo[i].outerborder)->in_image = TRUE;
    			++Ncells;
    		}
    		MoveDownKist(imageinfo[i].outerborder);
    	}

		Ncells_o=Ncells;

       // loop through points in ith image
    	for(j=0,Ncells=0;j<imageinfo[i].Npoints;++j){

    		if( imageinfo[i].points[j].gridsize > 1.01*rmax/Ngrid_block){  /* only refine largest grid size in image*/

    			// get real point in tree
    			point = imageinfo[i].points[j].image->image;
    			assert(point->gridsize > 0);

    			if(batch){
    				points_to_refine.push_back(point);
    			}else{
    				grid->RefineLeaf(lens,point,kappa_off);
        			// TODO: This seems like it shouldn't have been there. imageinfo[i].points[j].gridsize /= Ngrid_block;
    			}
    			++count;

    			++Ncells;
    		}
    	}

      /* loop through outer border of ith image */

      MoveToTopKist(imageinfo[i].outerborder);
      for(j=0;j<imageinfo[i].outerborder->Nunits();++j){
    	  if( getCurrentKist(imageinfo[i].outerborder)->gridsize > 1.01*rmax/Ngrid_block){/* only refine largest grid size in image*/

    		  point = getCurrentKist(imageinfo[i].outerborder);
    		  assert(point->gridsize > 0);

    		  if(point->in_image){ /* point has not been refined yet */
    			  ++count;

      			if(batch){
      				points_to_refine.push_back(point);
      			}else{
      				grid->RefineLeaf(lens,point,kappa_off);
      			}
      			++Ncells;
    			  point->in_image = FALSE;  // unmark so that it wouldn't double refine
    		  }
    		  //imageinfo[i].outerborderlist->current->gridsize /= Ngrid_block;/**/
    	  }
    	  MoveDownKist(imageinfo[i].outerborder);
      }

    }

    if(count > 0) ++number_of_refined;
  } /* end of image loop */

  if(batch){
  	grid->RefineLeaves(lens,points_to_refine,kappa_off);
  	points_to_refine.clear();
  }

   return number_of_refined;
}



/**
 * \brief	refines only inner and outer edges of all images
 * <pr>
 *    Does NOT update borders or image points.  These need to be re-found elsewhere.
 *
 * criterion = 0 stops refining when each image reaches error limit
 *           = 1 stops refining when grid resolution is smaller than res_target in all images
 *           = 2 stop when area of a cell reaches res_target * area of all images
 *
 * </pr>
 */
long refine_edges(
		LensHndl lens
		,GridHndl grid
		,ImageInfo *imageinfo
		,unsigned long Nimages
		,double res_target
		,short criterion
		,bool kappa_off
		,Kist<Point> * newpointskist  /// returns a Kist of the points that were added to the grid on this pass, if == NULL will not be added
		,bool batch){
	 //printf("entering refine_edges\n");

	if(newpointskist) newpointskist->Empty();

	if(Nimages < 1) return 0;

	long i,j,Ncells=0,Ncells_o=0,count=0;
	Point *point;
	double area_total=0;
	std::vector<Point *> points_to_refine;

	// count border points
	if( criterion == 2 ) for(i=0,area_total=0.0;i<Nimages;++i) area_total += imageinfo[i].area;

	for(i=0,Ncells=0;i<Nimages;++i){

		MoveToTopKist(imageinfo[i].outerborder);
		for(j=0;j<imageinfo[i].outerborder->Nunits();++j){
			if( ( criterion==0 && pow(getCurrentKist(imageinfo[i].outerborder)->gridsize,2)/imageinfo[i].area > res_target )
							|| ( criterion==1 && getCurrentKist(imageinfo[i].outerborder)->gridsize > res_target)
							|| ( criterion==2 && pow(getCurrentKist(imageinfo[i].outerborder)->gridsize,2)/area_total > res_target) ){
				++Ncells;
				getCurrentKist(imageinfo[i].outerborder)->in_image = TRUE;  // Temporarily mark point so they are not double refined
			}
			MoveDownKist(imageinfo[i].outerborder);
		}
		MoveToTopKist(imageinfo[i].innerborder);
		for(j=0;j<imageinfo[i].innerborder->Nunits();++j){
			if( criterion==0 && pow(getCurrentKist(imageinfo[i].innerborder)->gridsize,2)/imageinfo[i].area > res_target ) ++Ncells;
			if( criterion==1 && getCurrentKist(imageinfo[i].innerborder)->gridsize > res_target) ++Ncells;
			if( criterion==2 && pow(getCurrentKist(imageinfo[i].innerborder)->gridsize,2)/area_total > res_target ) ++Ncells;
			MoveDownKist(imageinfo[i].innerborder);
		}

		//printf("       %i Nouter=%i Ninner=%i\n",i,imageinfo[i].outerborder->Npoints,imageinfo[i].innerborderkist->Npoints);
	}

	if(Ncells==0) return 0;

	Ncells_o=Ncells;

	for(i=0,Ncells=0;i<Nimages;++i){
			/* loop through outer border of ith image */

		MoveToTopKist(imageinfo[i].outerborder);
		for(j=0;j<imageinfo[i].outerborder->Nunits();++j){

			if( ( criterion==0 && pow(getCurrentKist(imageinfo[i].outerborder)->gridsize,2)/imageinfo[i].area > res_target )
				|| ( criterion==1 && getCurrentKist(imageinfo[i].outerborder)->gridsize > res_target)
				|| ( criterion==2 && pow(getCurrentKist(imageinfo[i].outerborder)->gridsize,2)/area_total > res_target) ){

				point = getCurrentKist(imageinfo[i].outerborder);
    			assert(point->gridsize > 0);

				if(point->in_image){ // point has not been refined yet
					++count;

					if(batch) points_to_refine.push_back(point);
					else{
						Point *i_points = grid->RefineLeaf(lens,point,kappa_off);
						if(newpointskist && i_points != NULL){
							for(unsigned int k=0;k < i_points->head; ++k) newpointskist->InsertAfterCurrent(&i_points[k]);
						}
					}

					point->in_image = FALSE;

					++Ncells;
				}
			}
			MoveDownKist(imageinfo[i].outerborder);
		}

		MoveToTopKist(imageinfo[i].innerborder);
		for(j=0;j<imageinfo[i].innerborder->Nunits();++j){

			if( ( criterion==0 && pow(getCurrentKist(imageinfo[i].innerborder)->gridsize,2)/imageinfo[i].area > res_target )
					|| (criterion==1 && getCurrentKist(imageinfo[i].innerborder)->gridsize > res_target)
					|| (criterion==2 && pow(getCurrentKist(imageinfo[i].innerborder)->gridsize,2)/area_total > res_target) ){

				//point = getCurrentKist(imageinfo[i].innerborder);
    			//assert(point->gridsize > 0);

				//if( getCurrentKist(imageinfo[i].innerborderkist)->gridsize == point->gridsize){ /* point has not been refined yet */
    			++count;

    			if(batch) points_to_refine.push_back(getCurrentKist(imageinfo[i].innerborder));
    			else{
    				Point *i_points = grid->RefineLeaf(lens,getCurrentKist(imageinfo[i].innerborder),kappa_off);
    				if(newpointskist && i_points != NULL)
    					for(unsigned int k=0;k < i_points->head; ++k) newpointskist->InsertAfterCurrent(&i_points[k]);
    			}
       			++Ncells;
			}
			MoveDownKist(imageinfo[i].innerborder);
		}

	}

	if(batch){
		Point *i_points = grid->RefineLeaves(lens,points_to_refine,kappa_off);
		if(newpointskist && i_points != NULL){
			for(unsigned int k=0;k < i_points->head; ++k) newpointskist->InsertAfterCurrent(&i_points[k]);
		}
	  	points_to_refine.clear();
	}

	return count;
}

/**
 * <pr>
 * 	refines only inner and outer edges of all images
 *    then update borders, area and area_error so that the
 *    images do not need to be found every time, but they should
 *    be well defined before using this routine.  It will not break
 *    already existing images up into sub-images and will not merge
 *    images that are found to be connected at higher resolution.
 *
 *    Note:   in_image marks must be set on entry
 *
 *    on exit:
 *       new image points are not added to imageinfo->points so they will be out of date
 *       but they are added to imageinfo[i].imagekist
 *       *image_overlap = true if the images merge, but not guaranteed when they separate.
 *
 * criterion = 0 stops refining when each image reaches error limit or is smaller than res_target
 *           = 1 stops refining when grid resolution is smaller than res_target in all images
 *           = 2 stop when area of a cell reaches res_target * area of all images
 *
 * </pr>
 */
long refine_edges2(LensHndl lens,double *y_source,double r_source,GridHndl grid
		,ImageInfo *imageinfo,bool *image_overlap,unsigned long Nimages,double res_target
		,short criterion,bool kappa_off,bool batch){

	 //printf("entering refine_edges2\n");

	if(Nimages < 1) return 0;

	long i,j,k,n,Ncells=0,Ncells_o=0,count=0,Npoints=0;
	Point *i_points,*point;
	Kist<Point> * neighborkist = new Kist<Point>;
	double tmp_area=0,area_total=0;
	bool addinner;
	std::vector<Point *> points_to_refine;

	// count border points
	if( criterion==2) for(i=0,area_total = 0.0;i<Nimages;++i) area_total += imageinfo[i].area;

	*image_overlap=false;

	// loop through outer border of ith image

	for(i=0;i<Nimages;++i){

        if(!(imageinfo[i].ShouldNotRefine)){
		// count border points that needs to be refined
		MoveToTopKist(imageinfo[i].outerborder);
		for(j=0,Ncells=0;j<imageinfo[i].outerborder->Nunits();++j){
			if( ( criterion==0 && pow( getCurrentKist(imageinfo[i].outerborder)->gridsize,2)/imageinfo[i].area > res_target )
					|| ( criterion==1 && getCurrentKist(imageinfo[i].outerborder)->gridsize > res_target)
					|| ( criterion==2 && pow( getCurrentKist(imageinfo[i].outerborder)->gridsize,2)/area_total > res_target)){
				++Ncells;
				getCurrentKist(imageinfo[i].outerborder)->in_image = TRUE; // temporary mark
			}
			MoveDownKist(imageinfo[i].outerborder);
		}
		MoveToTopKist(imageinfo[i].innerborder);
		for(j=0;j<imageinfo[i].innerborder->Nunits();++j){
			if( criterion==0 && pow(getCurrentKist(imageinfo[i].innerborder)->gridsize,2)/imageinfo[i].area > res_target ) ++Ncells;
			if( criterion==1 && getCurrentKist(imageinfo[i].innerborder)->gridsize > res_target) ++Ncells;
			if( criterion==2 && pow(getCurrentKist(imageinfo[i].innerborder)->gridsize,2)/area_total > res_target ) ++Ncells;
			MoveDownKist(imageinfo[i].innerborder);
		}
		//printf("Ncells=%i total \n",Ncells);
		//printf("       %i Nouter=%i Ninner=%i\n",i,imageinfo[i].outerborder->Npoints,imageinfo[i].innerborderkist->Nunits());

		EmptyKist(neighborkist);

		if(Ncells > 0){

			//i_points=NewPointArray((Ngrid_block*Ngrid_block-1)*Ncells,true);
			Ncells_o=Ncells;

			MoveToTopKist(imageinfo[i].outerborder);
			for(j=0,Ncells=0;j<imageinfo[i].outerborder->Nunits();++j){

			//if( imageinfo[i].outerborder->current->gridsize > res_target){
				if( ( criterion==0 && pow( getCurrentKist(imageinfo[i].outerborder)->gridsize,2)/imageinfo[i].area > res_target )
						|| ( criterion==1 && getCurrentKist(imageinfo[i].outerborder)->gridsize > res_target)
						|| ( criterion==2 && pow( getCurrentKist(imageinfo[i].outerborder)->gridsize,2)/area_total > res_target)){

					point = getCurrentKist(imageinfo[i].outerborder);
	    			assert(point->gridsize > 0);

	                // point has not been refined yet
					if(point->in_image){
						++count;

						if(batch){
							points_to_refine.push_back(point);
						}else{
							i_points = grid->RefineLeaf(lens,point,kappa_off);
							sort_out_points(i_points,&imageinfo[i],r_source,y_source);
							/*
							//  sort new points into in and out of image
							//    and add them to inner and outer borders
							for(j=0;j<i_points->head;++j){
								if( sqrt(pow(i_points[j].image->x[0]-y_source[0],2)
									+ pow(i_points[j].image->x[1]-y_source[1],2)) < r_source){

									// mark points
									i_points[j].in_image = TRUE;
									i_points[j].image->in_image = TRUE;

									InsertAfterCurrentKist(imageinfo[i].innerborder,&(i_points[j]));
									InsertAfterCurrentKist(imageinfo[i].imagekist,&(i_points[j]));
									MoveDownKist(imageinfo[i].innerborder);

									imageinfo[i].area += pow(i_points[j].gridsize,2);
								}else{

									// un-mark points
									i_points[j].in_image = FALSE;
									i_points[j].image->in_image = FALSE;
								}
							}*/
						}

						point->in_image = FALSE;
					}else *image_overlap=true;

				}
				MoveDownKist(imageinfo[i].outerborder);
			}
			//printf("Ncells=%i in outer border the second time \n",Ncells);

			MoveToTopKist(imageinfo[i].innerborder);
			tmp_area=imageinfo[i].area;
			for(j=0;j<imageinfo[i].innerborder->Nunits();++j){

			//if( getCurrentKist(imageinfo[i].innerborderkist)->gridsize > res_target){
				if( ( criterion == 0 && pow(getCurrentKist(imageinfo[i].innerborder)->gridsize,2)/tmp_area > res_target )
						|| (criterion == 1 && getCurrentKist(imageinfo[i].innerborder)->gridsize > res_target)
						|| (criterion == 2 && pow(getCurrentKist(imageinfo[i].innerborder)->gridsize,2)/area_total > res_target)){

					//point=getCurrentKist(imageinfo[i].innerborderkist)->leaf->points;
					point = getCurrentKist(imageinfo[i].innerborder);
	    			assert(point->gridsize > 0);

					// point has not been refined yet
					//if( getCurrentKist(imageinfo[i].innerborderkist)->gridsize == point->gridsize){
						++count;

						if(batch){
							points_to_refine.push_back(point);
							// subtract area of new points from image area
							imageinfo[i].area -= (pow(grid->getNgrid_block(),2)-1)*pow(point->gridsize/grid->getNgrid_block(),2);
							if(imageinfo[i].gridrange[2] == point->gridsize/grid->getNgrid_block())
								imageinfo[i].gridrange[2] = point->gridsize/grid->getNgrid_block();
						}else{
							i_points = grid->RefineLeaf(lens,point,kappa_off);
							sort_out_points(i_points,&imageinfo[i],r_source,y_source);
							//  sort new points into in and out of image
							//    and add them to inner and outer borders
							/*
							for(j=0;j<i_points->head;++j){
								if( sqrt(pow(i_points[j].image->x[0]-y_source[0],2)
										+ pow(i_points[j].image->x[1]-y_source[1],2)) < r_source){

									// mark points
									i_points[j].in_image = TRUE;
									i_points[j].image->in_image = TRUE;

									InsertAfterCurrentKist(imageinfo[i].innerborder,&(i_points[j]));
									InsertAfterCurrentKist(imageinfo[i].imagekist,&(i_points[j]));
									MoveDownKist(imageinfo[i].innerborder);
									//PointCopyData(imageinfo[i].innerborder->current,&(i_points[j]));

									imageinfo[i].area += pow(i_points[j].gridsize,2);
								}else{

									// un-mark points
									i_points[j].in_image = FALSE;
									i_points[j].image->in_image = FALSE;

								}
							}
							 */

							// subtract area of new points from image area
							imageinfo[i].area -= (pow(grid->getNgrid_block(),2)-1)*pow(point->gridsize,2);
							if(imageinfo[i].gridrange[2] > point->gridsize)
								imageinfo[i].gridrange[2] = point->gridsize;
						}

						++Ncells;
				}
				MoveDownKist(imageinfo[i].innerborder);
			}


			if(batch){
				i_points = grid->RefineLeaves(lens,points_to_refine,kappa_off);
				if(i_points) sort_out_points(i_points,&imageinfo[i],r_source,y_source);
				points_to_refine.clear();
			}

			if(Ncells) imageinfo[i].gridrange[0]/=grid->getNgrid_block(); // maximum grid size in outerborder

			/*
			 *     Weed out of the borders the points that are no longer on the border
			 */

			EmptyKist(imageinfo[i].outerborder);

			Npoints=imageinfo[i].innerborder->Nunits();
			MoveToTopKist(imageinfo[i].innerborder);
			for(j=0;j<Npoints;++j){
				addinner=false;

				// update leaf pointer of inner border point if necessary
				/*** don't think this is necessary anymore
				if(getCurrentKist(imageinfo[i].innerborder)->leaf->npoints > 1){
					grid->i_tree->current = getCurrentKist(imageinfo[i].innerborder)->leaf;
					_FindBox(grid->i_tree,getCurrentKist(imageinfo[i].innerborder)->x);
					getCurrentKist(imageinfo[i].innerborder)->leaf = grid->i_tree->current;
				}*/

				grid->i_tree->FindAllBoxNeighborsKist(getCurrentKist(imageinfo[i].innerborder),neighborkist);

				MoveToTopKist(neighborkist);
				for(k=0;k < neighborkist->Nunits() ;++k){
					if(getCurrentKist(neighborkist)->in_image == FALSE){
						addinner=true;

						MoveToTopKist(imageinfo[i].outerborder);
						for(n=0;n<imageinfo[i].outerborder->Nunits();++n){
							if(getCurrentKist(imageinfo[i].outerborder) == getCurrentKist(neighborkist)) break;
								MoveDownKist(imageinfo[i].outerborder);
							}
							// add to outer border
						if(n==imageinfo[i].outerborder->Nunits()){
							InsertBeforeCurrentKist(imageinfo[i].outerborder,getCurrentKist(neighborkist));
							MoveUpKist(imageinfo[i].outerborder);
						}
					}
					MoveDownKist(neighborkist);
				}

				if(addinner) MoveDownKist(imageinfo[i].innerborder);
				else{
					TakeOutCurrentKist(imageinfo[i].innerborder);
					if(!AtTopKist(imageinfo[i].innerborder))
						MoveDownKist(imageinfo[i].innerborder);
				}
			}
		}
        }
	}  // end loop through images

	delete neighborkist;

	return count;
}

//  sort new points into in and out of image and add them to inner and outer borders
// Also adds the area of new image points to the image area.
void sort_out_points(Point *i_points,ImageInfo *imageinfo,double r_source,double y_source[]){
	for(unsigned long j=0;j<i_points->head;++j){
		if( sqrt(pow(i_points[j].image->x[0]-y_source[0],2)
			+ pow(i_points[j].image->x[1]-y_source[1],2)) < r_source){

			// mark points
			i_points[j].in_image = TRUE;
			i_points[j].image->in_image = TRUE;

			InsertAfterCurrentKist(imageinfo->innerborder,&(i_points[j]));
			InsertAfterCurrentKist(imageinfo->imagekist,&(i_points[j]));
			MoveDownKist(imageinfo->innerborder);

			imageinfo->area += pow(i_points[j].gridsize,2);
		}else{

			// un-mark points
			i_points[j].in_image = FALSE;
			i_points[j].image->in_image = FALSE;
		}
	}
}

/** Finds inner and outer border of image by bordering box neighbor method.
 * Maintained only for find_crit().  Not very efficient.
 * Should be replaced with findborder4() in any new implementation.
 */

void findborders2(TreeHndl i_tree,OldImageInfo *imageinfo){
	int i;
	unsigned long l,m,j;
	bool addinner;

	//printf("beginning findborders2\n");
	//checkTree(i_tree);

	//point=(Point *)malloc(sizeof(Point));
	EmptyKist(imageinfo->innerborder);
	EmptyKist(imageinfo->outerborder);

	imageinfo->gridrange[2] = 1.0e99; /* minimum grid size in image */
	imageinfo->gridrange[0] = 0.0; /* maximum grid size in outerborder */
	imageinfo->gridrange[1] = 0.0;      /* maximum grid size in image */

	if(imageinfo->Npoints < 1) return;

	Kist<Point> * neighborkist = new Kist<Point>;

	for(j=0;j<imageinfo->Npoints;++j){

		if(imageinfo->gridrange[1] < imageinfo->points[j].gridsize)
			imageinfo->gridrange[1] = imageinfo->points[j].gridsize;
		if(imageinfo->gridrange[2] > imageinfo->points[j].gridsize)
			imageinfo->gridrange[2] = imageinfo->points[j].gridsize;

		addinner=false;

		i_tree->FindAllBoxNeighborsKist(&(imageinfo->points[j]),neighborkist);

		MoveToTopKist(neighborkist);
		for(i=0;i<neighborkist->Nunits();++i){

			for(l=0;l<imageinfo->Npoints;++l) if( getCurrentKist(neighborkist)->id
					== imageinfo->points[l].id) break;

			if(l==imageinfo->Npoints){  // point is a neighbor
				addinner=true;
				// check if point is already in list
				MoveToTopKist(imageinfo->outerborder);
				for(m=0;m<imageinfo->outerborder->Nunits();++m){
					if( getCurrentKist(imageinfo->outerborder) == getCurrentKist(neighborkist) ) break;
					MoveDownKist(imageinfo->outerborder);
				}
				if(m==imageinfo->outerborder->Nunits()){
				// add point to outerborder
					InsertAfterCurrentKist(imageinfo->outerborder,getCurrentKist(neighborkist));
					MoveDownKist(imageinfo->outerborder);
				}
			}
			MoveDownKist(neighborkist);
		}

    // go to eight points of the compass
		/*for(i=-1;i<=1;++i){
			ray[0]=imageinfo->points[j].x[0]+1.001*i*imageinfo->points[j].gridsize/2;
			for(k=-1;k<=1;++k){
				if(i != 0 || k != 0){

					ray[1]=imageinfo->points[j].x[1]+1.001*k*imageinfo->points[j].gridsize/2;
					//printf("ray =  %e %e\n",ray[0],ray[1]);
					FindBoxPoint(i_tree,ray,point);
					//printf("ray =  %e %e point = %e %e\n",ray[0]
					//            ,ray[1],point->x[0],point->x[1]);

					//PrintPoint(point);
					if(point->id == imageinfo->points[j].id){  // error check
					    ERROR_MESSAGE();
						printf("ERROR: findborders2 FindBoxPoint returned same point j=%i i=%i k=%i\n"
								,j,i,k);
						exit(1);
					}

					for(l=0;l<imageinfo->Npoints;++l) if(point->id == imageinfo->points[l].id) break;

					if(l==imageinfo->Npoints){  // point is a neighbor

						addinner=true;
						// check if point is already in list
						MoveToTopKist(imageinfo->outerborder);
						for(m=0;m<imageinfo->outerborder->Npoints;++m){
							if(imageinfo->outerborder->current->id == point->id) break;
							MoveDownKist(imageinfo->outerborder);
						}
						if(m==imageinfo->outerborder->Npoints){
        				// add point to outerborder
							InsertAfterCurrent(imageinfo->outerborder,point->x,point->id,point->image);
							MoveDownKist(imageinfo->outerborder);
							PointCopyData(imageinfo->outerborder->current,point);
						}

						// if neighbor is on a smaller grid put two neighbors in outer border
						if(point->gridsize < imageinfo->points[j].gridsize){
							if(i==0 || k==0){
								ray[0]=point->x[0]+1.001*k*point->gridsize/2;
								ray[1]=point->x[1]+1.001*i*point->gridsize/2;
								FindBoxPoint(i_tree,ray,point);
								for(l=0;l<imageinfo->Npoints;++l) if(point->id == imageinfo->points[l].id) break;
								if(l==imageinfo->Npoints){  // point is a neighbor
				        			// check if point is already in list
									MoveToTopKist(imageinfo->outerborder);
									for(m=0;m<imageinfo->outerborder->Npoints;++m){
										if(imageinfo->outerborder->current->id == point->id) break;
										MoveDownKist(imageinfo->outerborder);
									}
									if(m==imageinfo->outerborder->Npoints){
										// add point to outerborder
										InsertAfterCurrent(imageinfo->outerborder,point->x,point->id,point->image);
										MoveDownKist(imageinfo->outerborder);
										PointCopyData(imageinfo->outerborder->current,point);
									}
								}

								ray[0]=point->x[0]-2.001*k*point->gridsize;
								ray[1]=point->x[1]-2.001*i*point->gridsize;
								FindBoxPoint(i_tree,ray,point);
								for(l=0;l<imageinfo->Npoints;++l) if(point->id == imageinfo->points[l].id) break;
								if(l==imageinfo->Npoints){  // point is a neighbor
									// check if point is already in list
									MoveToTopKist(imageinfo->outerborder);
									for(m=0;m<imageinfo->outerborder->Npoints;++m){
										if(imageinfo->outerborder->current->id == point->id) break;
										MoveDownKist(imageinfo->outerborder);
									}
									if(m==imageinfo->outerborder->Npoints){
										// add point to outerborder
										InsertAfterCurrent(imageinfo->outerborder,point->x,point->id,point->image);
										MoveDownKist(imageinfo->outerborder);
										PointCopyData(imageinfo->outerborder->current,point);
									}
								}


							}
						}

					}

				}
			}
		}*/
		if(addinner){
			// add point to innerborderkist
			InsertAfterCurrentKist(imageinfo->innerborder,imageinfo->points[j].image->image);
			MoveDownKist(imageinfo->innerborder);
		}
	}


	//checkTree(i_tree);
	//printf("exiting findborders2\n");
	delete neighborkist;

	return;
}

/**
 * Finds inner and outer borders of an image using
 * bordering box method
 *   same as findborder2() but it uses the in_image markers
 *   which makes it faster
 *   Note:  markers in_image must be set
 *
 *   Maintained only for find_crit().  Uses imageinfo->points[] array instead of imageinfo->imagekist
 * Should be replaced with findborder4() in any new implementation.
 */
void findborders3(TreeHndl i_tree,OldImageInfo *imageinfo){
	int i;
	unsigned long m,j;
	bool addinner;

	//printf("beginning findborders2\n");
	//checkTree(i_tree);

	//point=(Point *)malloc(sizeof(Point));
	EmptyKist(imageinfo->innerborder);
	EmptyKist(imageinfo->outerborder);

	imageinfo->gridrange[2] = 1.0e99; /* minimum grid size in image */
	imageinfo->gridrange[0] = 0.0; /* maximum grid size in outerborder */
	imageinfo->gridrange[1] = 0.0;      /* maximum grid size in image */

	if(imageinfo->Npoints < 1) return;

	Kist<Point> * neighborkist = new Kist<Point>;

	for(j=0;j<imageinfo->Npoints;++j){

		if(imageinfo->gridrange[1] < imageinfo->points[j].gridsize)
			imageinfo->gridrange[1] = imageinfo->points[j].gridsize;
		if(imageinfo->gridrange[2] > imageinfo->points[j].gridsize)
			imageinfo->gridrange[2] = imageinfo->points[j].gridsize;

		addinner=false;

		i_tree->FindAllBoxNeighborsKist(&(imageinfo->points[j]),neighborkist);

		MoveToTopKist(neighborkist);
		for(i=0;i<neighborkist->Nunits();++i){

			if( getCurrentKist(neighborkist)->in_image == FALSE){  // point is a neighbor
				addinner=true;
				// check if point is already in list
				MoveToTopKist(imageinfo->outerborder);
				for(m=0;m<imageinfo->outerborder->Nunits();++m){
					if( getCurrentKist(imageinfo->outerborder) == getCurrentKist(neighborkist) ) break;
					MoveDownKist(imageinfo->outerborder);
				}
				if(m==imageinfo->outerborder->Nunits()){
				// add point to outerborder
					InsertAfterCurrentKist(imageinfo->outerborder,getCurrentKist(neighborkist));
					MoveDownKist(imageinfo->outerborder);
				}
			}
			MoveDownKist(neighborkist);
		}

		if(addinner){
			// add point to innerborderkist
			InsertAfterCurrentKist(imageinfo->innerborder,imageinfo->points[j].image->image);  // need to put in real point
			MoveDownKist(imageinfo->innerborder);
		}
	}

	delete neighborkist;

	return;
}

void xygridpoints(Point *i_points,double range,double *center,long Ngrid_1d,short remove_center){
  /* make a new rectolinear grid of points on the image plane **/
  /* and link them to points on the source plane **/
  /* remove_center = 0 include center point of grid */
  /*              != 1 leave out center point of grid if Ngrid_1d is odd, Ngrid_1d*Ngrid_1d-1 points outputted */
  /* warning: memory for i_points must be allocated before entering */
  long i,j;
  static long id=0;

  if(id==0) initialgridsize=range/(Ngrid_1d-1);

  if(remove_center && (Ngrid_1d%2 == 1)){
    /*i_points=NewPointArray(Ngrid_1d*Ngrid_1d-1);*/
    for(i=0,j=0;i<Ngrid_1d*Ngrid_1d;++i){

      if( (2*(i/Ngrid_1d)/(Ngrid_1d-1) == 1) && (i%Ngrid_1d == Ngrid_1d/2+1) ) j=1;
      i_points[i-j].id=id;
      ++id;
      Utilities::PositionFromIndex(i,i_points[i-j].x,Ngrid_1d,range,center);
      //i_points[i-j].x[0] = center[0] + range*( 1.0*(i%Ngrid_1d)/(Ngrid_1d-1) - 0.5 );
      //i_points[i-j].x[1] = center[1] + range*( 1.0*(i/Ngrid_1d)/(Ngrid_1d-1) - 0.5 );
      i_points[i-j].gridsize=range/(Ngrid_1d-1);
    }

  }else{
    /*i_points=NewPointArray(Ngrid_1d*Ngrid_1d);*/
    for(i=0;i<Ngrid_1d*Ngrid_1d;++i){
      i_points[i].id=id;
      ++id;
      Utilities::PositionFromIndex(i,i_points[i].x,Ngrid_1d,range,center);
      //i_points[i].x[0] = center[0] + range*( 1.0*(i%Ngrid_1d)/(Ngrid_1d-1) - 0.5 );
      //i_points[i].x[1] = center[1] + range*( 1.0*(i/Ngrid_1d)/(Ngrid_1d-1) - 0.5 );
      i_points[i].gridsize=range/(Ngrid_1d-1);
    }
  }

  return;
}

void combineCloseImages(double linkinglength,ImageInfo *imageinfo,int *Nimages
		,int *NewNimages,int NimagesMax){
	unsigned long i,j,k;

	assert(*Nimages < NimagesMax);

	*NewNimages=*Nimages;
	for(i=0;i<(*NewNimages);++i){

		for(j=i+1;j<(*NewNimages);++j){

			if(sqrt( pow(imageinfo[i].centroid[0] - imageinfo[j].centroid[0],2)
				   + pow(imageinfo[i].centroid[1] - imageinfo[j].centroid[1],2) )
					< linkinglength){

				// update image center

				imageinfo[i].centroid[0] = (imageinfo[i].area*imageinfo[i].centroid[0]
				                              + imageinfo[j].area*imageinfo[j].centroid[0])/
						(imageinfo[i].area + imageinfo[j].area);
				imageinfo[i].centroid[1] = (imageinfo[i].area*imageinfo[i].centroid[1]
				                              + imageinfo[j].area*imageinfo[j].centroid[1])/
						(imageinfo[i].area + imageinfo[j].area);

				// combine points in image and on borders
				MoveToBottomKist(imageinfo[j].imagekist);
				while(imageinfo[j].imagekist->Nunits() > 0) InsertAfterCurrentKist(imageinfo[i].imagekist,TakeOutCurrentKist(imageinfo[j].imagekist) );
				MoveToBottomKist(imageinfo[j].innerborder );
				while(imageinfo[j].innerborder ->Nunits() > 0) InsertAfterCurrentKist(imageinfo[i].innerborder ,TakeOutCurrentKist(imageinfo[j].innerborder) );
				MoveToBottomKist(imageinfo[j].outerborder );
				while(imageinfo[j].outerborder ->Nunits() > 0) InsertAfterCurrentKist(imageinfo[i].outerborder ,TakeOutCurrentKist(imageinfo[j].outerborder) );

				// images that are smaller than the target error are not included
				if(DMIN(imageinfo[i].area,imageinfo[j].area)
					/DMAX(imageinfo[i].area,imageinfo[j].area) < 2*FracResTarget) --(*Nimages);

				imageinfo[i].area += imageinfo[j].area;

				// move image j to end of list
				for(k=j;k<*NewNimages-1;++k) SwapImages(&(imageinfo[k]),&(imageinfo[k+1]));
				j=i;
				--(*NewNimages);
			}
		}
	}

	return ;
}

void SwapImages(OldImageInfo *image1,OldImageInfo *image2){
	Point *point;
	unsigned long Npoints,i;
	double tmp;
	Kist<Point> * list;

	point = image1->points;
	image1->points = image2->points;
	image2->points = point;

	Npoints = image1->Npoints;
	image1->Npoints = image2->Npoints;
	image2->Npoints = Npoints;
	Npoints = image1->ShouldNotRefine;
	image1->ShouldNotRefine = image2->ShouldNotRefine;
	image2->ShouldNotRefine = Npoints;

	tmp = image1->area;
	image1->area = image2->area;
	image2->area = tmp;

	tmp = image1->area_error;
	image1->area_error = image2->area_error;
	image2->area_error = tmp;

	for(i=0;i<3;++i){
		tmp = image1->gridrange[i];
		image1->gridrange[i] = image2->gridrange[i];
		image2->gridrange[i] = tmp;
	}

	for(i=0;i<2;++i){
		tmp = image1->centroid[i];
		image1->centroid[i] = image2->centroid[i];
		image2->centroid[i] = tmp;
	}

	//list = image1->imagekist;
	//image1->imagekist = image2->imagekist;
	//image2->imagekist = list;

	list = image1->innerborder;
	image1->innerborder = image2->innerborder;
	image2->innerborder = list;

	list = image1->outerborder;
	image1->outerborder = image2->outerborder;
	image2->outerborder = list;

	return ;
}
