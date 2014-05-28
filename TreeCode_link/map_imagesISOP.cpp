/*
 * map_images.c
 *
 *  Created on: Oct 6, 2010
 *      Author: bmetcalf
 */
#include "slsimlib.h"
#include "isop.h"
const bool verbose = false;
const PosType FracResTarget = 3.0e-5;

/**
 *  \brief Find images and refine them based on their surface brightness distribution.
 *
 *  Uses ImageFinding::find_images_kist() to initially find and refine images and then using 
 *  IntegrateFluxInCell().  If cells cannot be integrated refiments are madebased on a surface brightness
 *  based criterion to refine the most important parts of the lens.
 *
 *  map_images is intended for mapping images of sources more complicated than simple circles.
 *
 */
void ImageFinding::map_imagesISOP(
	    /// model
		LensHndl lens
		,Source *source
		,GridHndl grid          /// Tree of grid points
		,int *Nimages           /// number of images found
		,ImageInfo *imageinfo   /// information on each image
		,int NimageMax          /// Size of imageinfo array on entry.  This could increase if more images are found
		,PosType rmax            /// Maximum size of source on souce plane.  The entire source must be within this distance from
		                        ///  source->getX()[]
		,PosType res_min        /// requred resolution of image, typically the pixel size of the final image
		,PosType initial_size    /// Initial size of source for telescoping, 0 to start from the initial grid size.
		                        /// If < 0 no telescoping is used and only the already existing points are used to
		                        /// to initiate the image finding.
		,ExitCriterion criterion  /// see data type
		,bool divide_images    /// if true will divide images and apply the exit criterion to them separately.
    ,bool verbose
                                  ){

	assert(lens);
	assert(grid->s_tree);
	assert(grid->i_tree);
	assert(imageinfo->imagekist);

	unsigned long Nimagepoints,Ntmp;
	PosType tmp,area_tot,flux;
	long i,j;
	static int oldNimages=0;
	bool go;
	PosType center[2],y[2],sb,xx[2],source_flux = 0;
	int Nsources;

	assert(res_min > 0);
	assert(rmax > 0);

  if(source->getRadius() <= 0.0){ERROR_MESSAGE(); printf("ERROR: find_images, point source must have a resolution target\n"); exit(1);}

	if(initial_size == 0) initial_size = 3*(grid->i_tree->top->boundary_p2[0]
                                          - grid->i_tree->top->boundary_p1[0])/grid->getInitNgrid();

	unsigned long tmp_count = 0;

  if(verbose) std::cout << "number of grid points before ImageFinding::find_images_kist: "
    << grid->getNumberOfPoints() << std::endl;
    
  ImageFinding::find_images_kist(lens,source->getX(),source->getRadius(),grid,Nimages
                                 ,imageinfo,NimageMax,&Nimagepoints,0.0,false,0,verbose);
  
  if(verbose) std::cout << "number of grid points after ImageFinding::find_images_kist: "<< grid->getNumberOfPoints() << std::endl;
  Nsources = 1;

	if(verbose) printf("total number of points after telescope: %li\n",grid->getNumberOfPoints());

	// Set all in_image flags to false.  This should not be necessary.  !!!!!!
	grid->ClearAllMarks();
	tmp = grid->RefreshSurfaceBrightnesses(source);

	if(tmp == 0.0){  // no flux was found
	  imageinfo->imagekist->Empty();
	  *Nimages = 0;
	  for(i=0; i < NimageMax ; ++i) imageinfo[i].area = 0.0;
	  return;
	}
	assert(tmp > 0.0 || imageinfo->imagekist->Nunits() == 0);

	// take out points with no flux
	if(verbose) printf("before taking out zero surface brightness points: %li\n",imageinfo->imagekist->Nunits());
	Ntmp = imageinfo->imagekist->Nunits();
	for(i=0,imageinfo->imagekist->MoveToTop(),area_tot=0; i < Ntmp ; ++i ){

		flux = imageinfo->imagekist->getCurrent()->surface_brightness * pow(imageinfo->imagekist->getCurrent()->gridsize,2);

		area_tot += flux;

		if(imageinfo->imagekist->getCurrent()->surface_brightness > 0){

			imageinfo->imagekist->getCurrent()->in_image = TRUE;
			imageinfo->imagekist->getCurrent()->image->in_image = TRUE;
			imageinfo->imagekist->Down();
		}else{
			imageinfo->imagekist->getCurrent()->in_image = FALSE;
			imageinfo->imagekist->getCurrent()->image->in_image = FALSE;
			if(imageinfo->imagekist->AtTop()) go = false; else go = true;
			imageinfo->imagekist->TakeOutCurrent();
			if(go) imageinfo->imagekist->Down();
		}

	}
  
	if(verbose) printf("after taking out zero surface brightness points: %li\n"
                     ,imageinfo->imagekist->Nunits());

	if(imageinfo->imagekist->Nunits() == 0){
		*Nimages = 0;
		for(i=0; i < NimageMax ; ++i) imageinfo[i].area = 0.0;
		return;
	}

	// divide up images
	if(divide_images) divide_images_kist(grid->i_tree,imageinfo,Nimages,NimageMax);
	else *Nimages = 1;
  
	if(verbose) printf("number of images after first division is %i\n",*Nimages);

	/////////////////////////////////////////////
	// link image points lists for each image
	// and calculate surface brightness at each point
	/////////////////////////////////////////////

  Boo outcome;
  long count=0;
  Point *current;

	// ****** calculate surface brightnesses and flux of each image   ******
	for(i=0,area_tot=0.0; i < *Nimages ; ++i){

		imageinfo[i].ShouldNotRefine = 0;
		imageinfo[i].uniform_mag = unchecked;

		findborders4(grid->i_tree,&(imageinfo[i]));

		imageinfo[i].imagekist->MoveToTop();
		imageinfo[i].area = 0.0;
		do{
      
      current = imageinfo[i].imagekist->getCurrent();
      ImageFinding::IntegrateFluxInCell(current,*source,1.0e-2,outcome);
      if(outcome != TRUE || current->gridsize > res_min ){
        ++count;
        imageinfo[i].imagekist->getCurrent()->flag = false;
      }else{
        imageinfo[i].imagekist->getCurrent()->flag = true;
      }
      
			imageinfo[i].area += pow(imageinfo[i].imagekist->getCurrent()->gridsize,2)
         *(imageinfo[i].imagekist->getCurrent()->surface_brightness);

		}while(imageinfo[i].imagekist->Down());

		area_tot += imageinfo[i].area;
		//printf("   %i area = %e\n",i,imageinfo[i].area);
		assert(imageinfo[i].area >= 0);
	}
  
  if(count){

	/********************************************************
	 ******* refine images based on flux in each pixel ******
	 *******************************************************/
    i=0;
    while( 0 < ImageFinding::refine_grid_on_imageISOP(lens,source,grid,imageinfo,Nimages
                                                      ,Nsources,NimageMax,FracResTarget,res_min
                                                ,criterion,divide_images) ) ++i;
      
  }

	// find image centroid
	for(i=0;i<*Nimages;++i){
		// reset these flags
		imageinfo[i].ShouldNotRefine = 0;
		imageinfo[i].uniform_mag = unchecked;

		imageinfo[i].imagekist->MoveToTop();
		imageinfo[i].centroid[0] = 0.0;
		imageinfo[i].centroid[1] = 0.0;
    imageinfo[i].area = 0.0;
		for(j = 0 ; j < imageinfo[i].imagekist->Nunits() ; ++j,imageinfo[i].imagekist->Down() ){
			imageinfo[i].area += pow(imageinfo[i].imagekist->getCurrent()->gridsize,2)
               *imageinfo[i].imagekist->getCurrent()->surface_brightness;
			imageinfo[i].centroid[0] += imageinfo[i].imagekist->getCurrent()->x[0]
          *pow(imageinfo[i].imagekist->getCurrent()->gridsize,2)
					*imageinfo[i].imagekist->getCurrent()->surface_brightness;
			imageinfo[i].centroid[1] += imageinfo[i].imagekist->getCurrent()->x[1]
          *pow(imageinfo[i].imagekist->getCurrent()->gridsize,2)
					*imageinfo[i].imagekist->getCurrent()->surface_brightness;

			imageinfo[i].imagekist->getCurrent()->in_image = FALSE;  // re-set marks
			imageinfo[i].imagekist->getCurrent()->image->in_image = FALSE;  // re-set marks
      imageinfo[i].imagekist->getCurrent()->flag = false;
      
      
		}
		if(imageinfo[i].getNimagePoints() > 0 ){
			imageinfo[i].centroid[0] /= imageinfo[i].area;
			imageinfo[i].centroid[1] /= imageinfo[i].area;
		}

		//printf("  %i  centroid = %e %e N = %li\n",i,imageinfo[i].centroid[0],imageinfo[i].centroid[1]
		//                                    ,imageinfo[i].Npoints);
	}

	return ;
}

/**
 * \brief refines grid according to criterion based on flux in each grid cell
 *
 *  points with in_image == 1 are tested for refinement
 *
 * criterion = TotalArea    stops refining when error in total area reaches res_target
 * 	         = EachImage    stops refining when each image reaches error limit or is smaller than res_target
 *           = Resolution   stops refining when grid resolution is smaller than res_target in all images
 *           = FillHoles    does the same as EachImage but also refines the neighbors to the cells that fulfill
 *                            the EachImage criterion
 *
 *     Warning:     set imageinfo[i].ShouldNotRefine = 0 and imageinfo[i].uniform_mag = unchecked;
 *
 */
int ImageFinding::refine_grid_on_imageISOP(Lens *lens,Source *source,GridHndl grid
                                           ,ImageInfo *imageinfo,int *Nimages
                                           ,int Nsources,int NimageMax
                                           ,const PosType res_target,const PosType res_min
                                           ,ExitCriterion criterion,bool divide_images
                                           ){
  
	//printf("entering refine_grid\n");
  
  if((*Nimages) < 1) return 0;
  
  int number_of_refined; /* Ngrid_block must be odd */
  PosType total_area,r;
  Point *i_points,*current;
  unsigned long Ncells,Nold,j,i;
  bool reborder=false,redivide=false;
  
  for(i=0,total_area=0;i<(*Nimages);++i) total_area += imageinfo[i].area;
  for(i=0,Nold=0;i<(*Nimages);++i) Nold += imageinfo[i].imagekist->Nunits();
  if(total_area == 0.0) return 0;
  
  Kist<Point> * nearest=new Kist<Point>;
  std::vector<Point *> points_to_refine;
  
  if(verbose) printf(" entering refine_on_image(), number of points %li\n",grid->getNumberOfPoints());
  
  number_of_refined=0;
  Boo outcome;
  for(i=0,Ncells=0;i<(*Nimages);++i){
        
	  assert(imageinfo[i].area >= 0.0);
	  if(imageinfo[i].ShouldNotRefine == 0 && imageinfo[i].uniform_mag == no){  // If image was not refined on the last round do not revisit unless the images are re-divided.
		  imageinfo[i].ShouldNotRefine = 1;
      
		  imageinfo[i].imagekist->MoveToTop();
		  for(j = 0 ; j < imageinfo[i].imagekist->Nunits() ; ++j,imageinfo[i].imagekist->Down() ){
        
        current = imageinfo[i].imagekist->getCurrent();
        
        if(current->flag == false){
          ImageFinding::IntegrateFluxInCell(current,*source,1.0e-2,outcome);
          if(outcome == TRUE && current->gridsize < res_min){
            current->flag = true;
          }else if(
            ImageFinding::RefinePoint2(current,grid->i_tree
                                      ,imageinfo[i].area,total_area,1.0,criterion,res_target,nearest)
                   ){
          
            ++Ncells;
          
            // adjust contribution to image flux from point that will be refined.
            assert(imageinfo[i].area > 0.0);
          
            // Determine if image point to be refined is a inner border point
            if(reborder == false){
              if(nearest->Nunits() == 0 ) grid->i_tree->FindAllBoxNeighborsKist(current,nearest);
              nearest->MoveToTop();
              do{
                if(nearest->getCurrent()->in_image != TRUE){
                  reborder = true;
                  break;
                }
              }while(nearest->Down());
            }
          
            imageinfo[i].area -=  pow(current->gridsize,2)
            *(current->surface_brightness);
          
            assert(imageinfo[i].area >= 0);
          
            assert(current->leaf->child1 == NULL);
            assert(current->leaf->child2 == NULL);
          
            
            points_to_refine.push_back(current);
            
            imageinfo[i].area +=  pow(current->gridsize/grid->getNgrid_block(),2)
                    *(current->surface_brightness);
          
            imageinfo[i].ShouldNotRefine = 0;   // mark to continue refinement on next round
          
          }
        }
      }
      
        
      i_points = grid->RefineLeaves(lens,points_to_refine);
      ImageFinding::check_sb_add(source,&imageinfo[i],i_points,1.0,Nold,number_of_refined);
      points_to_refine.clear();
        
      
      // If inner border point was refined re-do the borders
      if(reborder){
        if(imageinfo[i].outerborder->MoveToTop()){
          do{
            imageinfo[i].outerborder->getCurrent()->in_image = FALSE;
            //if(imageinfo[i].outerborder->getCurrent()->surface_brightness > 0)
            //  point = imageinfo[i].outerborder->getCurrent();
          }while(imageinfo[i].outerborder->Down());
        }
        findborders4(grid->i_tree,&(imageinfo[i]));
        
        reborder = false;
      }
      
      assert(imageinfo[i].area >= 0.0);
      // Do the same for the outerborder.
      
      imageinfo[i].outerborder->MoveToTop();
      for(j = 0 ; j < imageinfo[i].outerborder->Nunits() ; ++j,imageinfo[i].outerborder->Down() ){
        
          // TODO: This was taken out and i'm not sure if it was needed.
          //assert(getCurrentKist(imageinfo[i].outerborder)->surface_brightness == 0);
        
        current = imageinfo[i].outerborder->getCurrent();
        if(current->flag == false){
          ImageFinding::IntegrateFluxInCell(current,*source,1.0e-2,outcome);
          if(outcome == TRUE && current->gridsize < res_min){
            current->flag = true;
          }else if(
             ImageFinding::RefinePoint2(current,grid->i_tree
                                      ,imageinfo[i].area,total_area,1.0,criterion,res_target,nearest)
           ){
          
            if(current->in_image != MAYBE){
              current->in_image = MAYBE;
            
              reborder = true;  // Since the border has been refined, re do the borders.
              ++Ncells;
            
              assert(current->leaf->child1 == NULL);
              assert(current->leaf->child2 == NULL);

              points_to_refine.push_back(current);
            
              imageinfo[i].ShouldNotRefine = 0;   // mark for another look next time
            
            }else{
              /// there is an overlap in images so the images will need to be re-divided
              redivide = true;
            }
          }
        }
		  } // loop through outer border
      
		  assert(imageinfo[i].area >= 0.0);
      
      i_points = grid->RefineLeaves(lens,points_to_refine);
      ImageFinding::check_sb_add(source,&imageinfo[i],i_points,1.0,Nold,number_of_refined);
      points_to_refine.clear();
		  
      
		  // if the outerborder was refined, but no image overlap has been detected
		  //   recalculate the border
		  if(reborder){
        
			  if(imageinfo[i].outerborder->MoveToTop()){
				  do{
					  imageinfo[i].outerborder->getCurrent()->in_image = FALSE;
					  //if(imageinfo[i].outerborder->getCurrent()->surface_brightness > 0) point = imageinfo[i].outerborder->getCurrent();
				  }while(imageinfo[i].outerborder->Down());
			  }
			  findborders4(grid->i_tree,&(imageinfo[i]));
			  assert(imageinfo[i].outerborder->Nunits() > 0);
        
			  // re-set markers to MAYBE so overlaps can be detected
			  if(imageinfo[i].outerborder->MoveToTop()){
				  do{
					  imageinfo[i].outerborder->getCurrent()->in_image = MAYBE;
				  }while(imageinfo[i].outerborder->Down());
			  }
        
			  reborder = false;
		  }
	  } // if previously refined
  } // loop through images
  
  delete nearest;
  
  // reset the flag in outer borders
  for(i=0;i<(*Nimages);++i){
    
	  if(imageinfo[i].outerborder->MoveToTop()){
		  do{
			  imageinfo[i].outerborder->getCurrent()->in_image = FALSE;
		  }while(imageinfo[i].outerborder->Down());
	  }
  }
  
  if( number_of_refined == 0 || (divide_images && redivide) ){
		// Put all the images together.
		imageinfo->imagekist->MoveToBottom();
		for(i=1 ; i < *Nimages ; ++i){
			imageinfo[i].imagekist->MoveToTop();
			while(imageinfo[i].imagekist->Nunits() > 0)
        imageinfo->imagekist->InsertAfterCurrent(imageinfo[i].imagekist->TakeOutCurrent());
      imageinfo->imagekist->Down();
		}
    
		// divide up images
		divide_images_kist(grid->i_tree,imageinfo,Nimages,NimageMax);
    
		// find borders again and include them in the image
		for(i=0 ; i < (*Nimages) ; ++i){
      
			findborders4(grid->i_tree,&(imageinfo[i]));
      
			//assert(imageinfo[i].outerborder->Nunits() > 0);
      
			// recalculate flux that has been changed by divide_images_kist()
			imageinfo[i].ShouldNotRefine = 0;
			imageinfo[i].uniform_mag = unchecked;

      imageinfo[i].area = 0.0;
      imageinfo[i].imagekist->MoveToTop();
      do{
        imageinfo[i].area += pow(imageinfo[i].imagekist->getCurrent()->gridsize,2)
          * (imageinfo[i].imagekist->getCurrent()->surface_brightness);
      }while(imageinfo[i].imagekist->Down());
			
		}
  }
  
  return number_of_refined;
}


/** \brief Integrate the flux within a cell using the isoparametric expansion.
 *
 *   If outcome returns as FALSE or MAYBE the surface brightness will be returned without 
 *   integrating.
 *
 *   The result is returned in the point->surface_brightness in surface brightness
 *   units.
 *
 *   The point must have already been linked and entered into a tree.
 *
 *   If 
 */
void ImageFinding::IntegrateFluxInCell(
                                       Point *point      /// center point of cell to be integrated
                                       ,Source &source   /// source to be integrated
                                       ,float tolerance  /// tolerance to which isop expantion is checked
                                       ,Boo &outcome     /// FALSE if not 8 neighbors, MAYBE if isop expansion not valid, TRUE is source is integarted successfully
                                       ){

  point->surface_brightness = source.SurfaceBrightness(point->x);

  // Check that
  if(point->leaf->neighbors.size() != 8){
    outcome = FALSE;
    return;
  }
  
  std::vector<Point *> neighbors(point->leaf->neighbors.size());
  int i = 0;
  for(std::list<Branch *>::iterator it = point->leaf->neighbors.begin();
      it != point->leaf->neighbors.end() ; ++it,++i){
		if((*it)->npoints == 0){
      outcome = FALSE;
      return;
    }
		neighbors[i] = (*it)->points;
	}
  
  // sort neighbors into counterclockwise order from bottom left
  int j=0,k=7;
  for(i=0;i<8;++i){
    if(neighbors[i]->x[1] < point->x[1]){
      std::swap(neighbors[i],neighbors[j]);
      ++j;
    }
    if(neighbors[i]->x[1] > point->x[1]){
      std::swap(neighbors[i],neighbors[k]);
      --k;
    }
  }
  assert(j == 2);
  assert(k == 5);
  j=0;
  k=2;
  for(i=0;i<3;++i){
    if(neighbors[i]->x[0] < point->x[0]){
      std::swap(neighbors[i],neighbors[j]);
      ++j;
    }
    if(neighbors[i]->x[0] > point->x[0]){
      std::swap(neighbors[i],neighbors[k]);
      --k;
    }
  }
  if(neighbors[3]->x[0] < neighbors[4]->x[0] ) std::swap(neighbors[4],neighbors[3]);
  std::swap(neighbors[4],neighbors[7]);
  j=4;
  k=6;
  for(i=4;i<7;++i){
    if(neighbors[i]->x[0] > point->x[0]){
      std::swap(neighbors[i],neighbors[j]);
      ++j;
    }
    if(neighbors[i]->x[0] < point->x[0]){
      std::swap(neighbors[i],neighbors[k]);
      --k;
    }
  }
  
  PosType x[8],y[8],flux=0.0;
  for(i=0;i<8;++i){
    x[i] = neighbors[i]->image->x[0];
    y[i] = neighbors[i]->image->x[1];
  }
  
  double scale = sqrt( (x[0]-x[4])*(x[0]-x[4]) + (y[0]-y[4])*(y[0]-y[4]) );
  // see if isop expansion is valid
  if( fabs(point->image->x[0] - SLsimLib::isop(x,0,0)) > tolerance*scale){outcome = MAYBE; return;}
  if( fabs(point->image->x[1] - SLsimLib::isop(y,0,0)) > tolerance*scale){outcome = MAYBE; return;}
  
  // goto each quadrant and integrate the function.
  {
    x[0] = neighbors[0]->image->x[0];
    y[0] = neighbors[0]->image->x[1];
    interpfrom2Points(neighbors[0],neighbors[1],&x[1],&y[1]);
    x[2] = neighbors[1]->image->x[0];
    y[2] = neighbors[1]->image->x[1];
    interpfrom2Points(neighbors[1],point,&x[3],&y[3]);
    x[4] = point->x[0];
    y[4] = point->x[1];
    interpfrom2Points(neighbors[7],point,&x[5],&y[5]);
    x[6] = neighbors[7]->image->x[0];
    y[6] = neighbors[7]->image->x[1];
    interpfrom2Points(neighbors[7],neighbors[0],&x[7],&y[7]);
    flux +=  SLsimLib::isop_render(source,x,y,0,0,1,1);
  }
  {
    x[0] = neighbors[1]->image->x[0];
    y[0] = neighbors[1]->image->x[1];
    interpfrom2Points(neighbors[1],neighbors[2],&x[1],&y[1]);
    x[2] = neighbors[2]->image->x[0];
    y[2] = neighbors[2]->image->x[1];
    interpfrom2Points(neighbors[2],neighbors[3],&x[3],&y[3]);
    x[4] = neighbors[3]->image->x[0];
    y[4] = neighbors[3]->image->x[1];
    interpfrom2Points(neighbors[3],point,&x[5],&y[5]);
    x[6] = point->x[0];
    y[6] = point->x[1];
    interpfrom2Points(neighbors[1],point,&x[7],&y[7]);
    flux +=  SLsimLib::isop_render(source,x,y,-1,0,0,1);
  }
  {
    x[0] = point->x[0];
    y[0] = point->x[1];
    interpfrom2Points(neighbors[3],point,&x[1],&y[1]);
    x[2] = neighbors[3]->image->x[0];
    y[2] = neighbors[3]->image->x[1];
    interpfrom2Points(neighbors[3],neighbors[4],&x[3],&y[3]);
    x[4] = neighbors[4]->image->x[0];
    y[4] = neighbors[4]->image->x[1];
    interpfrom2Points(neighbors[4],neighbors[5],&x[5],&y[5]);
    x[6] = neighbors[5]->image->x[0];
    y[6] = neighbors[5]->image->x[1];
    interpfrom2Points(neighbors[5],point,&x[7],&y[7]);
    flux +=  SLsimLib::isop_render(source,x,y,-1,-1,0,0);
  }
  {
    x[0] = neighbors[7]->image->x[0];
    y[0] = neighbors[7]->image->x[1];
    interpfrom2Points(neighbors[7],point,&x[1],&y[1]);
    x[2] = point->x[0];
    y[2] = point->x[1];
    interpfrom2Points(neighbors[5],point,&x[3],&y[3]);
    x[4] = neighbors[5]->image->x[0];
    y[4] = neighbors[5]->image->x[1];
    interpfrom2Points(neighbors[5],neighbors[6],&x[5],&y[5]);
    x[6] = neighbors[6]->image->x[0];
    y[6] = neighbors[6]->image->x[1];
    interpfrom2Points(neighbors[6],neighbors[7],&x[7],&y[7]);
    flux +=  SLsimLib::isop_render(source,x,y,0,-1,1,0);
  }
  
  point->surface_brightness = flux/point->gridsize/point->gridsize;
  outcome = TRUE;
}

/** \brief Finds the source position of a point that lies half way between points p1 and p2 on the image plane by third order interpolation.
 */
void ImageFinding::interpfrom2Points(Point const * p1,Point const * p2,PosType *x,PosType *y){
  PosType dx = p1->x[0] - p2->x[0],dy = p1->x[1] - p2->x[1];
  PosType l = sqrt(dx*dx +dy*dy);
  dx /=l;
  dy /=l;
  
  // The sign convention agrees with GLAMER paper II here which should be correct.
  double A1p1 = (1 - p1->kappa - p1->gamma[0])*dx - (p1->gamma[1] - p1->gamma[2])*dy;
  double A2p1 = (1 - p1->kappa + p1->gamma[0])*dy - (p1->gamma[1] + p1->gamma[2])*dx;
  double A1p2 = (1 - p2->kappa - p2->gamma[0])*dx - (p2->gamma[1] - p2->gamma[2])*dy;
  double A2p2 = (1 - p2->kappa + p2->gamma[0])*dy - (p2->gamma[1] + p2->gamma[2])*dx;
  
  /*
   x = p1->image->x[0] + (p2->image->x[0] - p1->image->x[0])*s*s*(3-2*s)
   + l*A1p1*s*(1-s)*(1-s) + l*A1p2*s*s*(1-s);
   y = p1->image->x[1] + (p2->image->x[1] - p1->image->x[1])*s*s*(3-2*s)
   + l*A2p1*s*(1-s)*(1-s) + l*A2p2*s*s*(1-s);
   */
  
  *x = (p2->image->x[0] + p1->image->x[0])/2
  + l*( A1p1 + A1p2 )/8;
  *y = (p2->image->x[1] + p1->image->x[1])/2
  + l*( A2p1 + A2p2 )/8;
  
  return;
}

