/*
 * map_images.c
 *
 *  Created on: Oct 6, 2010
 *      Author: bmetcalf
 */
#include "slsimlib.h"
#include "isop.h"
#include "map_images.h"
#include <thread>

const bool verbose = false;
const PosType FracResTarget = 3.0e-5;

/** \ingroup ImageFinding
 *  \brief Find images and refine them based on their surface brightness distribution.
 *
 *  Uses ImageFinding::find_images_kist() to initially find and refine images and then using 
 *  IF_routines::IntegrateFluxInCell().  If cells cannot be integrated refiments are madebased on a surface brightness
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
    ,std::vector<ImageInfo> &imageinfo   /// information on each image
		,PosType rmax            /// Maximum size of source on souce plane.  The entire source must be within this distance from
		                        ///  source->getX()[]
		,PosType res_min        /// requred resolution of image, typically the pixel size of the final image
		,PosType initial_size    /// Initial size of source for telescoping, 0 to start from the initial grid size.
		                        /// If < 0 no telescoping is used and only the already existing points are used to
		                        /// to initiate the image finding.
		,ExitCriterion criterion  /// see data type
		,bool divide_images    /// if true will divide images and apply the exit criterion to them separately.
    ,bool int_on           /// if true the flux in each cell is integrated, if false the surface brightness at the center point of the cell is used
    ,bool verbose
                                  ){

  if(imageinfo.size() < 5) imageinfo.resize(5);
	assert(lens);
	assert(grid->s_tree);
	assert(grid->i_tree);
	assert(imageinfo[0].imagekist);

	unsigned long Nimagepoints,Ntmp;
	PosType tmp,area_tot,flux;
	long i,j;
	int Nsources;

	assert(res_min > 0);
	assert(rmax > 0);

  if(source->getRadius() <= 0.0){ERROR_MESSAGE(); printf("ERROR: find_images, point source must have a resolution target\n"); exit(1);}

	if(initial_size == 0) initial_size = 3*(grid->i_tree->getTop()->boundary_p2[0]
                                          - grid->i_tree->getTop()->boundary_p1[0])/grid->getInitNgrid();


  if(verbose) std::cout << "number of grid points before ImageFinding::find_images_kist: "
    << grid->getNumberOfPoints() << std::endl;
    
  ImageFinding::find_images_kist(lens,source->getX(),source->getRadius(),grid,Nimages
                                 ,imageinfo,&Nimagepoints,0.0,true,0,verbose);
  
  //assert(*Nimages == 1);
  
// TODO: Seems to miss some images.  Do not know why.
// TODO: Need some why of refining when the expansion is not valid or the original cells from above are too large.
  
  if(verbose) std::cout << "number of grid points after ImageFinding::find_images_kist: "
      << grid->getNumberOfPoints() << std::endl;
  Nsources = 1;

	if(verbose) printf("total number of points after telescope: %li\n",grid->getNumberOfPoints());

	// Set all in_image flags to false.  This should not be necessary.  !!!!!!
	grid->ClearAllMarks();
	tmp = grid->RefreshSurfaceBrightnesses(source);

	if(tmp == 0.0){  // no flux was found
	  imageinfo[0].imagekist->Empty();
	  *Nimages = 0;
	  for(i=0; i < imageinfo.size() ; ++i) imageinfo[i].area = 0.0;
	  return;
	}
  
  area_tot = 0;
  for(i=0;i<*Nimages;++i){

    // take out points with no flux
    if(verbose) printf("before taking out zero surface brightness points: %li\n",imageinfo[i].imagekist->Nunits());
    Ntmp = imageinfo[i].imagekist->Nunits();
    PosType gridmax=0;
    size_t ii;
    for(ii=0,imageinfo[i].imagekist->MoveToTop(); ii < Ntmp ; ++ii ){

      flux = imageinfo[i].imagekist->getCurrent()->surface_brightness * pow(imageinfo[i].imagekist->getCurrent()->gridsize,2);

      area_tot += flux;

      //if(imageinfo[i].imagekist->getCurrent()->surface_brightness > 0){

        imageinfo[i].imagekist->getCurrent()->in_image = YES;
        imageinfo[i].imagekist->getCurrent()->image->in_image = YES;
        imageinfo[i].imagekist->Down();
        gridmax=MAX(imageinfo[i].imagekist->getCurrent()->gridsize,gridmax);
      /*}else{
        imageinfo[i].imagekist->getCurrent()->in_image = NO;
        imageinfo[i].imagekist->getCurrent()->image->in_image = NO;
        if(imageinfo[i].imagekist->AtTop()) go = false; else go = true;
        imageinfo[i].imagekist->TakeOutCurrent();
        if(go) imageinfo[i].imagekist->Down();
      }*/

    }
    if(imageinfo[i].gridrange[0] > gridmax) imageinfo[i].gridrange[0] = gridmax;
  
    // If cell is not smaller than res_min refine it and add the new points to the image
    if(imageinfo[i].gridrange[0] > res_min){
      Kist<Point>::iterator it = imageinfo[i].imagekist->BottomIt();
      std::vector<Point*> points_to_refine;
      Point *new_points;

      while(imageinfo[i].gridrange[0] > res_min){
        for(it = imageinfo[i].imagekist->BottomIt()
            ;!(it.atend());++it){
          if((*it).gridsize > res_min){
            points_to_refine.push_back(&(*it));
          }
        }
        new_points = grid->RefineLeaves(lens,points_to_refine);
        for(size_t ii=0;ii<new_points->head;++ii){
          imageinfo[i].imagekist->InsertBeforeCurrent(&new_points[ii]);
        }
        if(new_points->head > 1) imageinfo[i].gridrange[0] /= 3;
        if(imageinfo[i].gridrange[0] < imageinfo[i].gridrange[2])
          imageinfo[i].gridrange[1] = imageinfo[i].gridrange[2] = imageinfo[i].gridrange[0];
      }
    }
    
    if(verbose) printf("after taking out zero surface brightness points: %li\n"
                     ,imageinfo[i].imagekist->Nunits());
    
  }
	if(area_tot == 0){
		*Nimages = 0;
		for(i=0; i < imageinfo.size() ; ++i) imageinfo[i].area = 0.0;
		return;
	}

	/////////////////////////////////////////////
	// link image points lists for each image
	// and calculate surface brightness at each point
	/////////////////////////////////////////////

  long count=0;//,count_tot=0,count_moreNeighbors=0,count_isopfail=0,count_res=0;
  
	// ****** calculate surface brightnesses and flux of each image   ******
	for(i=0,area_tot=0.0; i < *Nimages ; ++i){

    size_t counti = 0;
		imageinfo[i].ShouldNotRefine = 0;

		//findborders4(grid->i_tree,&(imageinfo[i]));

		imageinfo[i].area = 0.0;
    
    // integrate cells in parallel
    if(int_on && imageinfo[i].imagekist->Nunits() > 0){
      int nthreads = N_THREADS < imageinfo[i].imagekist->Nunits() ? N_THREADS : imageinfo[i].imagekist->Nunits();
      std::thread thread[N_THREADS];
      Kist<Point>::iterator its[2];
      int block_size = (int)(imageinfo[i].imagekist->Nunits()/nthreads);
      PosType area[N_THREADS];
      size_t counts[N_THREADS];
      
      //std::cout << "block_size " << block_size << " nthreads " << nthreads << " points in image " << imageinfo[i].imagekist->Nunits() << std::endl;
      int jj;
      its[0] = its[1] = imageinfo[i].imagekist->BottomIt();
      for(jj=0;jj< block_size
          && its[1] !=  imageinfo[i].imagekist->TopIt();++jj) ++its[1];      
      
      nthreads = 0;
      //std::cout << "Start thread " << nthreads << "  = " << jj << "  " << (its[1] == imageinfo[i].imagekist->getTopIt()) << std::endl;
      thread[nthreads] = std::thread(IF_routines::IntegrateCellsParallel,its[0],its[1],source
                                     ,&area[nthreads],&counts[nthreads]);
      ++nthreads;
      while(its[1] != imageinfo[i].imagekist->TopIt()){
        ++its[1];
        its[0] = its[1];
        
        for(jj=0;jj < block_size
            && its[1] !=  imageinfo[i].imagekist->TopIt();++jj) ++its[1];
        //std::cout << "Start thread " << nthreads << " jj = " << jj << "  " << (its[1] == imageinfo[i].imagekist->getTopIt()) << std::endl;

        thread[nthreads] = std::thread(IF_routines::IntegrateCellsParallel,its[0],its[1],source
                                         ,&area[nthreads],&counts[nthreads]);
        ++nthreads;
      }
      assert(nthreads <= N_THREADS);
        
      for(int ii=0;ii<nthreads;++ii) thread[ii].join();
      for(int ii=0;ii<nthreads;++ii){ imageinfo[i].area += area[ii]; counti += counts[ii];}
      assert(counti == imageinfo[i].imagekist->Nunits());
      count += counti;
    }
    
    // Warning: this will overwrite the areas
    // divide up images
    //if(divide_images) divide_images_kist(grid->i_tree,imageinfo,Nimages,NimageMax);
    //else *Nimages = 1;

    if(verbose) printf("number of images after first division is %i\n",*Nimages);

    /***
		do{
           current = imageinfo[i].imagekist->getCurrent();
           ImageFinding::IF_routines::IntegrateFluxInCell(current,*source,1.0e-2,outcome);
           if(outcome != YES || current->gridsize > res_min ){
              ++count;
              imageinfo[i].imagekist->getCurrent()->flag = false;
               if(outcome == NO) ++count_moreNeighbors;
               if(outcome == MAYBE) ++count_isopfail;
               if(current->gridsize > res_min) ++count_res;
           }else{
              imageinfo[i].imagekist->getCurrent()->flag = true;
           }
         
           imageinfo[i].area += pow(imageinfo[i].imagekist->getCurrent()->gridsize,2)
               *(imageinfo[i].imagekist->getCurrent()->surface_brightness);

            ++count_tot;
		}while(imageinfo[i].imagekist->Down());
*/
		area_tot += imageinfo[i].area;
		//printf("   %i area = %e\n",i,imageinfo[i].area);
		assert(imageinfo[i].area >= 0);
	}
  /*
    std::cout << "number of IF_routines::IntegrateFluxInCell failures "<< count << "  "
    << count*100.0/count_tot << "%" << std::endl;
    std::cout << "    "<< count_moreNeighbors*100.0/count_tot << "% not 8 neighbors" << std::endl;
    std::cout << "    "<< count_isopfail*100.0/count_tot << "% isop expansion fails" << std::endl;
    std::cout << "    "<< count_res*100.0/count_tot << "% pixel too large" << std::endl;
    */
  if(count*0){

	/********************************************************
	 ******* refine images based on flux in each pixel ******
	 *******************************************************/
    i=0;
    while( 0 < ImageFinding::IF_routines::refine_grid_on_imageISOP(lens,source,grid,imageinfo,Nimages
                                                      ,Nsources,FracResTarget
                                                      ,res_min,res_min*res_min*1.0e-3
                                                      ,criterion,divide_images) ) ++i;
      
  }
    
    std::cout << "  map_imagesISOP: " << i << " refinement steps were taken " << std::endl;

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

			imageinfo[i].imagekist->getCurrent()->in_image = NO;  // re-set marks
			imageinfo[i].imagekist->getCurrent()->image->in_image = NO;  // re-set marks
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
int ImageFinding::IF_routines::refine_grid_on_imageISOP(Lens *lens,Source *source,GridHndl grid
                                           ,std::vector<ImageInfo> &imageinfo,int *Nimages
                                           ,int Nsources
                                           ,PosType res_target
                                           ,PosType res_min
                                           ,PosType res_source_area
                                           ,ExitCriterion criterion,bool divide_images
                                           ){
  
	//printf("entering refine_grid\n");
  
  if((*Nimages) < 1) return 0;
  
  int number_of_refined; /* Ngrid_block must be odd */
  PosType total_area;
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
    
	  if(imageinfo[i].ShouldNotRefine == 0 ){  // If image was not refined on the last round do not revisit unless the images are re-divided.
		  imageinfo[i].ShouldNotRefine = 1;
      
		  imageinfo[i].imagekist->MoveToTop();
		  for(j = 0 ; j < imageinfo[i].imagekist->Nunits() ; ++j,imageinfo[i].imagekist->Down() ){
        
        current = imageinfo[i].imagekist->getCurrent();
        
        if(current->flag == false){
          ImageFinding::IF_routines::IntegrateFluxInCell(current,*source,1.0e-2,outcome);
          if(outcome == YES && current->gridsize < res_min){
            current->flag = true;
          }else if( (outcome == MAYBE && current->image->leaf->area() > res_source_area ) || current->gridsize > res_min
                   || ImageFinding::IF_routines::RefinePoint2(current,grid->i_tree
                                      ,imageinfo[i].area,total_area,1.0,criterion,res_target,nearest)
                   ){
          
          ++Ncells;
          imageinfo[i].ShouldNotRefine = 0;   // mark to continue refinement on next round
              
          
          // adjust contribution to image flux from point that will be refined.
          assert(imageinfo[i].area > 0.0);
          
          // Determine if image point to be refined is a inner border point
          if(reborder == false){
            if(nearest->Nunits() == 0 ) grid->i_tree->FindAllBoxNeighborsKist(current,nearest);
            nearest->MoveToTop();
            do{
              if(nearest->getCurrent()->in_image != YES){
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
          
        }
        }
      }
      
        
      i_points = grid->RefineLeaves(lens,points_to_refine);
      ImageFinding::IF_routines::check_sb_add(source,&imageinfo[i],i_points,1.0,Nold,number_of_refined);
      points_to_refine.clear();
        
      
      // If inner border point was refined re-do the borders
      if(reborder){
        if(imageinfo[i].outerborder->MoveToTop()){
          do{
            imageinfo[i].outerborder->getCurrent()->in_image = NO;
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
      for(j = 0 ; j < imageinfo[i].outerborder->Nunits()*0 ; ++j,imageinfo[i].outerborder->Down() ){
        
          // TODO: This was taken out and i'm not sure if it was needed.
          //assert(getCurrentKist(imageinfo[i].outerborder)->surface_brightness == 0);
        
        current = imageinfo[i].outerborder->getCurrent();
        if(current->flag == false){
          ImageFinding::IF_routines::IntegrateFluxInCell(current,*source,1.0e-2,outcome);
          if(outcome == YES && current->gridsize < res_min){
            current->flag = true;
          }else if(
             ImageFinding::IF_routines::RefinePoint2(current,grid->i_tree
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
      ImageFinding::IF_routines::check_sb_add(source,&imageinfo[i],i_points,1.0,Nold,number_of_refined);
      points_to_refine.clear();
		  
      
		  // if the outerborder was refined, but no image overlap has been detected
		  //   recalculate the border
		  if(reborder){
        
			  if(imageinfo[i].outerborder->MoveToTop()){
				  do{
					  imageinfo[i].outerborder->getCurrent()->in_image = NO;
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
      }
  } // loop through images
  
  delete nearest;
  
  // reset the flag in outer borders
  for(i=0;i<(*Nimages);++i){
    
	  if(imageinfo[i].outerborder->MoveToTop()){
		  do{
			  imageinfo[i].outerborder->getCurrent()->in_image = NO;
		  }while(imageinfo[i].outerborder->Down());
	  }
  }
  
  if( number_of_refined == 0 || (divide_images && redivide) ){
		// Put all the images together.
		imageinfo[0].imagekist->MoveToBottom();
		for(i=1 ; i < *Nimages ; ++i){
			imageinfo[i].imagekist->MoveToTop();
			while(imageinfo[i].imagekist->Nunits() > 0)
        imageinfo[0].imagekist->InsertAfterCurrent(imageinfo[i].imagekist->TakeOutCurrent());
      imageinfo[0].imagekist->Down();
		}
    
		// divide up images
		divide_images_kist(grid->i_tree,imageinfo,Nimages);
    
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

/** \brief Integrate the flux within a series of cells using the isoparametric expansion 
 *    in a thread safe way.
 *
 *    The integration is done with ImageFinding::IF_routines::IntegrateFluxInCell().
 *    Points in the Kist from it1 to it2 inclusive are integrated. 
 *    source.SurfaceBrightness() must be thread safe.
 *    
 *    If outcome from ImageFinding::IF_routines::IntegrateFluxInCell() returns YES the 
 *    flag of that point is set to true.  Otherwise it is false.
 */

void ImageFinding::IF_routines::IntegrateCellsParallel(
      Kist<Point>::iterator it1   /// iterator pointing to first point in block
     ,Kist<Point>::iterator it2   /// iterator pointing to last point in block
     ,Source *source
     ,PosType *area              /// returns total flux in these cells
     ,size_t *count              /// returns total number of cells that were integrated
                                          ){
  
  Boo outcome;
  *area = 0.0;
  int i=0;
  *count = 0;
  for(Kist<Point>::iterator it = it1; it != it2 ; ++it,++i){
    if((*it).flag == false){
      ImageFinding::IF_routines::IntegrateFluxInCell(&(*it),*source,1.0e-2,outcome);
      ++(*count);
      *area += (*it).surface_brightness*(*it).gridsize*(*it).gridsize;
      if(outcome == YES) (*it).flag = true;
      else (*it).flag = false;
    }
  }
  if((*it2).flag == false){
    ImageFinding::IF_routines::IntegrateFluxInCell(&(*it2),*source,1.0e-2,outcome);
    ++(*count);
    *area += (*it2).surface_brightness*(*it2).gridsize*(*it2).gridsize;
    if(outcome == YES) (*it2).flag = true;
    else (*it2).flag = false;
  }
  return;
}

/** \brief Integrate the flux within a cell using the isoparametric expansion.
 *
 *   If outcome returns as NO or MAYBE the surface brightness will be returned without 
 *   integrating.
 *
 *   The result is returned in the point->surface_brightness in surface brightness
 *   units.
 *
 *   The point must have already been linked and entered into a tree.
 *
 *   If 
 */


void ImageFinding::IF_routines::IntegrateFluxInCell(
                                       Point *point      /// center point of cell to be integrated
                                       ,Source &source   /// source to be integrated
                                       ,float tolerance  /// tolerance to which isop expantion is checked
                                       ,Boo &outcome     /// NO if not 8 neighbors, MAYBE if isop expansion not valid, YES is source is integarted successfully
                                       ){

  point->surface_brightness = source.SurfaceBrightness(point->x);//*point->gridsize*point->gridsize;

  // Check that
  if(point->leaf->neighbors.size() != 8){
    outcome = NO;
    return;
  }
  
  std::vector<Point *> neighbors(8);
  int i = 0;
  for(std::list<Branch *>::iterator it = point->leaf->neighbors.begin();
      it != point->leaf->neighbors.end() ; ++it,++i){
		if((*it)->npoints == 0){
      outcome = NO;
      return;
    }
		neighbors[i] = (*it)->points;
	}
  
  // sort neighbors into counterclockwise order from bottom left
  
  std::sort(neighbors.begin(),neighbors.end(),Point::orderY);
  std::sort(neighbors.begin(),neighbors.begin() + 2,Point::orderX);
  std::sort(neighbors.begin() + 3,neighbors.end(),Point::orderXrev);
  if(neighbors[3]->x[1] > neighbors[4]->x[1] ) std::swap(neighbors[4],neighbors[3]);
  if(neighbors[6]->x[1] < neighbors[7]->x[1] ) std::swap(neighbors[6],neighbors[7]);
  
  //for(int ii=0;ii<8;++ii) std::cout << neighbors[ii]->x[0] - point->x[0] << "  " <<  neighbors[ii]->x[1] - point->x[1] << std::endl;
   
  PosType x[8],y[8],flux=0.0;
  for(i=0;i<8;++i){
    x[i] = neighbors[i]->image->x[0];
    y[i] = neighbors[i]->image->x[1];
  }
    
  // goto each quadrant and integrate the function.
  {
    x[0] = neighbors[0]->image->x[0];
    y[0] = neighbors[0]->image->x[1];
    IF_routines::interpfrom2Points(neighbors[0],neighbors[1],&x[1],&y[1]);
    x[2] = neighbors[1]->image->x[0];
    y[2] = neighbors[1]->image->x[1];
    IF_routines::interpfrom2Points(neighbors[1],point,&x[3],&y[3]);
    x[4] = point->x[0];
    y[4] = point->x[1];
    IF_routines::interpfrom2Points(neighbors[7],point,&x[5],&y[5]);
    x[6] = neighbors[7]->image->x[0];
    y[6] = neighbors[7]->image->x[1];
    IF_routines::interpfrom2Points(neighbors[7],neighbors[0],&x[7],&y[7]);
    flux +=  ISOP::isop_render(source,x,y,0,0,1,1);
  }
  {
    x[0] = neighbors[1]->image->x[0];
    y[0] = neighbors[1]->image->x[1];
    IF_routines::interpfrom2Points(neighbors[1],neighbors[2],&x[1],&y[1]);
    x[2] = neighbors[2]->image->x[0];
    y[2] = neighbors[2]->image->x[1];
    IF_routines::interpfrom2Points(neighbors[2],neighbors[3],&x[3],&y[3]);
    x[4] = neighbors[3]->image->x[0];
    y[4] = neighbors[3]->image->x[1];
    IF_routines::interpfrom2Points(neighbors[3],point,&x[5],&y[5]);
    x[6] = point->x[0];
    y[6] = point->x[1];
    IF_routines::interpfrom2Points(neighbors[1],point,&x[7],&y[7]);
    flux +=  ISOP::isop_render(source,x,y,-1,0,0,1);
  }
  {
    x[0] = point->x[0];
    y[0] = point->x[1];
    IF_routines::interpfrom2Points(neighbors[3],point,&x[1],&y[1]);
    x[2] = neighbors[3]->image->x[0];
    y[2] = neighbors[3]->image->x[1];
    IF_routines::interpfrom2Points(neighbors[3],neighbors[4],&x[3],&y[3]);
    x[4] = neighbors[4]->image->x[0];
    y[4] = neighbors[4]->image->x[1];
    IF_routines::interpfrom2Points(neighbors[4],neighbors[5],&x[5],&y[5]);
    x[6] = neighbors[5]->image->x[0];
    y[6] = neighbors[5]->image->x[1];
    IF_routines::interpfrom2Points(neighbors[5],point,&x[7],&y[7]);
    flux +=  ISOP::isop_render(source,x,y,-1,-1,0,0);
  }
  {
    x[0] = neighbors[7]->image->x[0];
    y[0] = neighbors[7]->image->x[1];
    IF_routines::interpfrom2Points(neighbors[7],point,&x[1],&y[1]);
    x[2] = point->x[0];
    y[2] = point->x[1];
    IF_routines::interpfrom2Points(neighbors[5],point,&x[3],&y[3]);
    x[4] = neighbors[5]->image->x[0];
    y[4] = neighbors[5]->image->x[1];
    IF_routines::interpfrom2Points(neighbors[5],neighbors[6],&x[5],&y[5]);
    x[6] = neighbors[6]->image->x[0];
    y[6] = neighbors[6]->image->x[1];
    IF_routines::interpfrom2Points(neighbors[6],neighbors[7],&x[7],&y[7]);
    flux +=  ISOP::isop_render(source,x,y,0,-1,1,0);
  }
  
  point->surface_brightness = flux;///point->gridsize/point->gridsize;
  
  double scale = sqrt( (x[0]-x[4])*(x[0]-x[4]) + (y[0]-y[4])*(y[0]-y[4]) );
  // see if isop expansion is valid
  if( fabs(point->image->x[0] - ISOP::isop(x,0,0)) > tolerance*scale){outcome = MAYBE; return;}
  if( fabs(point->image->x[1] - ISOP::isop(y,0,0)) > tolerance*scale){outcome = MAYBE; return;}

  outcome = YES;
}

/** \brief Finds the source position of a point that lies half way between points p1 and p2 on the image plane by third order interpolation.
 */
void ImageFinding::IF_routines::interpfrom2Points(Point const * p1,Point const * p2,PosType *x,PosType *y){
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

void partion(Point **neighbors,int n,PosType divide,int dimension){
  int j=0,k=n-1;
  while(j<k){
    if(neighbors[j]->x[dimension] < divide){
      ++j;
    }else{
      std::swap(neighbors[j],neighbors[k]);
      --k;
    }
  }
}
