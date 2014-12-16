 /*
 * map_images.c
 *
 *  Created on: Oct 6, 2010
 *      Author: bmetcalf
 */
#include "slsimlib.h"

//const float mumin = 0.3;  // actually the sqrt of the minimum magnification
const PosType FracResTarget = 3.0e-5;
//const float FracResTarget = 1.0e-4;
const PosType target_all = 1.0e-3;
//const int MinPoints = 100;  // Minimum number of points per image
const float tol_UniformMag = 1.0e-3;
const bool verbose = false;

/** \ingroup ImageFinding
 *  \brief Find images and refine them based on their surface brightness distribution.
 *
 *  Uses ImageFinding::find_images_kist() to initially find and refine images and then uses a surface brightness
 *  based criterion to refine the most important parts of the lens.
 *
 *  map_images is intended for mapping images of sources more complicated than simple circles.
 *
 */
void ImageFinding::map_images(
	    /// model
		LensHndl lens
		,Source *source
		,GridHndl grid          /// Tree of grid points
		,int *Nimages           /// number of images found
    ,std::vector<ImageInfo> &imageinfo   /// information on each image
		,PosType xmax            /// Maximum size of source on image plane.  The entire source must be within this distance from
		                        ///  source->getX()[]
		,PosType xmin            /// The smallest scale of the source measured on the lens plane.  The more accurate these
			                    /// 2 parameters are the less likely it is that an image will be missed.
		,PosType initial_size    /// Initial size of source for telescoping, 0 to start from the initial grid size.
		                        /// If < 0 no telescoping is used and only the already existing points are used to
		                        /// to initiate the image finding.
		,ExitCriterion criterion  /// see data type
		//,bool kappa_off         /// turns off calculation of surface density, shear, magnification and time delay
		,bool FindCenter        /// if the center of the source is not known this can be set to true and it will attempt to
		                        /// find the center, find the size of the source and determine if there is more than one source
			                    /// The entire source must be within xmax of source->getX() because this is the only
			                    /// region that will be scanned.
		,bool divide_images    /// if true will divide images and apply the exit criterion to them separately.
		){

	assert(lens);
	assert(grid->s_tree);
	assert(grid->i_tree);
      if(imageinfo.size() < 2) imageinfo.resize(2);
	assert(imageinfo[0].imagekist);
      
     

	unsigned long Nimagepoints,Ntmp;
	PosType tmp,area_tot,flux;
	static PosType oldy[2],oldr=0;
	//short moved;
	long i,j;
	//Point *i_points,*s_points;
	//Point *point;
	//time_t to;
	//ListHndl tmp_border_pointer;
	static int oldNimages=0;
	bool go;
	PosType center[2],y[2],sb,xx[2],source_flux = 0;
	int Nsources;
  PosType maxflux;

	assert(xmin > 0);
	assert(xmax > 0);

	if(initial_size == 0) initial_size = 3*(grid->i_tree->top->boundary_p2[0] - grid->i_tree->top->boundary_p1[0])/grid->getInitNgrid();

	unsigned long tmp_count = 0;

  std::vector<ImageInfo> sourceinfo(imageinfo.size());

	if(FindCenter){
		Point *newpoint;

		//*** find source and size by scanning source plane

		// Ntmp dictates the smallest source component that would be expected to be found
		Ntmp = 80*(unsigned long)(xmax/xmin + 1);
		xx[0] = xx[1] = 0.0;
		for(i=0;i<Ntmp*Ntmp;++i){
			Utilities::PositionFromIndex(i,y,Ntmp,2*xmax,xx);
			sb = source->SurfaceBrightness(y);
			//printf(" sb = %e  y = %e %e\n",sb,y[0],y[1]);
			if(sb > 0.0){

				newpoint = NewPointArray(1);
				newpoint->x[0] = y[0] + source->getX()[0];
				newpoint->x[1] = y[1] + source->getX()[1];
				newpoint->image = newpoint;
				newpoint->surface_brightness = sb;
				newpoint->gridsize = 2*xmax/(Ntmp-1);

				sourceinfo[0].imagekist->InsertAfterCurrent(newpoint);

				++tmp_count;
				source_flux += sb*pow(newpoint->gridsize,2);
			}
		}

		if(sourceinfo[0].imagekist->Nunits() == 0){
			Nsources = 0;
			*Nimages = 0;
			for(i=0; i < sourceinfo.size() ; ++i) imageinfo[i].area = 0.0;

			sourceinfo.clear();
			return;
		}else{
			DirtyFoF(sourceinfo,&Nsources,0);
		}
		if(verbose) printf(" number of source points found %li in %i sources\n",sourceinfo[0].imagekist->Nunits(),Nsources);

		// brakes source into four if there is only one
		//if(Nsources == 1 && sourceinfo->imagekist->Nunits() > 100) DirtyDivider(sourceinfo,&Nsources,NimageMax,(sourceinfo->imagekist->Nunits()/4+1));

		if(verbose) printf(" number of source points found %li in %i sources divided\n",sourceinfo[0].imagekist->Nunits(),Nsources);

		// split source into groups

		// calculate centroids and sizes of sources
		PosType *rs = new PosType[Nsources+1];
		for(i=0;i<Nsources; ++i){
			sourceinfo[i].centroid[0] = sourceinfo[i].centroid[1] = 0;
			sourceinfo[i].imagekist->MoveToTop();
			do{
				sourceinfo[i].centroid[0] += sourceinfo[i].imagekist->getCurrent()->x[0];
				sourceinfo[i].centroid[1] += sourceinfo[i].imagekist->getCurrent()->x[1];
			}while(sourceinfo[i].imagekist->Down());
			sourceinfo[i].centroid[0] /= sourceinfo[i].imagekist->Nunits();
			sourceinfo[i].centroid[1] /= sourceinfo[i].imagekist->Nunits();

			rs[i] = 0;
			sourceinfo[i].area = 0;  // actually the source flux
			sourceinfo[i].imagekist->MoveToTop();
			do{
				rs[i] = MAX(rs[i],pow(sourceinfo[i].imagekist->getCurrent()->x[0] - sourceinfo[i].centroid[0],2)
						+ pow(sourceinfo[i].imagekist->getCurrent()->x[1] - sourceinfo[i].centroid[1],2) );
				//xx[0] = getCurrentKist(sourceinfo[i].imagekist)->x[0] - source->getX()[0];
				//xx[1] = getCurrentKist(sourceinfo[i].imagekist)->x[1] - source->getX()[1];
				sourceinfo[i].area += source->SurfaceBrightness(sourceinfo[i].imagekist->getCurrent()->x)*pow(2*xmax/(Ntmp-1),2);

			}while(sourceinfo[i].imagekist->Down());

			rs[i] = (rs[i] == 0) ? xmin : sqrt(rs[i]);

			// filling factor or holiness of source
			if(verbose) printf("holiness of source %e\n",sourceinfo[i].imagekist->Nunits()*pow(2*xmax/(Ntmp-1)/rs[i],2)/pi);
			if(verbose) printf("     dx = %e %e rs = %e Npoints = %li\n",(source->getX()[0]-sourceinfo[i].centroid[0])/xmin
			                                               ,(source->getX()[1]-sourceinfo[i].centroid[1])/xmin
			                                               ,rs[i],sourceinfo[i].imagekist->Nunits());

		}

		if(sourceinfo[0].imagekist->Nunits() == 0) rs[0] = xmin;

		// **** what if undetected or one point in image
		// refine both grids to each source
		if(sourceinfo[0].imagekist->Nunits() > 0){
			center[0] = sourceinfo[0].centroid[0];
			center[1] = sourceinfo[0].centroid[1];
		}else{
			center[0] = source->getX()[0];
			center[1] = source->getX()[1];
		}
		ImageFinding::find_images_kist(lens,center,rs[0],grid,Nimages,imageinfo,&Nimagepoints,initial_size,true,0,false);
		for(i=1;i<Nsources;++i){
			if(sourceinfo[i].imagekist->Nunits() > 0){
				center[0] = sourceinfo[i].centroid[0];
				center[1] = sourceinfo[i].centroid[1];
			}else{
				center[0] = source->getX()[0];
				center[1] = source->getX()[1];
			}
			ImageFinding::find_images_kist(lens,center,rs[i],grid,Nimages,imageinfo,&Nimagepoints,xmax,true,0,false);
		}

		// free array of sourceinfo's'
		delete[] rs;
	}else{
		if(verbose) std::cout << "number of grid points before ImageFinding::find_images_kist: "<< grid->getNumberOfPoints() << std::endl;
		//ImageFinding::find_images_kist(lens,source->getX(),xmin,grid,Nimages,imageinfo,NimageMax,&Nimagepoints
		//		,0,true,0,false,true);
		ImageFinding::find_images_kist(lens,source->getX(),xmin,grid,Nimages,imageinfo,&Nimagepoints
				,0,false,0,verbose);
		if(verbose) std::cout << "number of grid points after ImageFinding::find_images_kist: "<< grid->getNumberOfPoints() << std::endl;
		Nsources = 1;
		sourceinfo[0].centroid[0] = source->getX()[0];
		sourceinfo[0].centroid[1] = source->getX()[1];
		sourceinfo[0].area = source->getTotalFlux();
	}

  
/*
	if(FindCenter){
		Point *points;
		TreeHndl tree;
		PosType **centers,*rs;

		//set up grid for scanning for source
		Ntmp = 2*(unsigned long)(source->source_r_out/xmin + 1);
		Ntmp = (Ntmp > 10000) ? 10000 : Ntmp;
		Ntmp = 2*prevpower(Ntmp);
		points = NewPointArray(Ntmp*Ntmp);
		xx[0] = xx[1] = 0.0;
		xygridpoints(points,2*source->source_r_out,xx,Ntmp,0);
		tree = BuildTree(points,Ntmp*Ntmp);

		// find points that have some surface brightness
		EmptyKist(imageinfo->imagekist);
		for(i=0;i<Ntmp*Ntmp;++i){
			points[i].image = &points[i];
			if(source->SurfaceBrightness(points[i].x) > 0.0){
				InsertAfterCurrentKist(imageinfo->imagekist,&points[i]);
			}
		}

		if(imageinfo->imagekist->Nunits() == 0){
			FreePointArray(points);
			*Nimages = 0;

			return ;
		}

		if(imageinfo->imagekist->Nunits() > 1){
			divide_images_kist(tree,imageinfo,&Nsources,NimageMax);
		}else{
			Nsources = 1;
		}

		printf("Nsources = %i\n",Nsources);
		centers = dmatrix(0,Nsources-1,0,1);
		rs = (PosType *) malloc(Nsources*sizeof(PosType));

		for(i=0;i<Nsources;++i){
			centers[i][0] = centers[i][1] = 0.0;
			MoveToTopKist(imageinfo[i].imagekist);
			do{
				centers[i][0] += getCurrentKist(imageinfo[i].imagekist)->x[0];
				centers[i][1] += getCurrentKist(imageinfo[i].imagekist)->x[1];
				printf("      id = %li x = %e %e sb = %e\n",getCurrentKist(imageinfo->imagekist)->id,getCurrentKist(imageinfo[i].imagekist)->x[0]
				          ,getCurrentKist(imageinfo[i].imagekist)->x[1],source->SurfaceBrightness(getCurrentKist(imageinfo[i].imagekist)->x));

			}while(MoveDownKist(imageinfo[i].imagekist));
			centers[i][0] /= imageinfo[i].imagekist->Nunits();
			centers[i][1] /= imageinfo[i].imagekist->Nunits();

			rs[i] = 0.0;
			MoveToTopKist(imageinfo[i].imagekist);
			do{
				rs[i] = MAX(rs[i],sqrt( pow(centers[i][0] - getCurrentKist(imageinfo[i].imagekist)->x[0],2)
						              + pow(centers[i][1] - getCurrentKist(imageinfo[i].imagekist)->x[1],2) ) );
			}while(MoveDownKist(imageinfo[i].imagekist));

			printf(" sb = %e\n",source->SurfaceBrightness(centers[i]));
			centers[i][0] += source->getX()[0];
			centers[i][1] += source->getX()[1];

			if(rs[i] == 0.0) rs[i] = 2*source->source_r_out/Ntmp;
		}

		ImageFinding::find_images_kist(centers[0],rs[0],grid,Nimages,imageinfo,NimageMax,&Nimagepoints,initial_size,true,0,false,true);
		for(i=1;i<Nsources;++i){
			ImageFinding::find_images_kist(centers[i],rs[i],grid,Nimages,imageinfo,NimageMax,&Nimagepoints,source->source_r_out,true,0,false,true);
		}

		free_dmatrix(centers,0,Nsources-1,0,1);
		free(rs);
		freeTree(tree);
	}else{

		ImageFinding::find_images_kist(source->getX(),xmin,grid,Nimages
			  ,imageinfo,NimageMax,&Nimagepoints,initial_size,true,0,false,true);
	}
*/

	//if(oldr==0) oldr=source->source_r_out;
	//if((oldy[0]==source->getX()[0])*(oldy[1]==source->getX()[1])*(oldr > source->source_r_out)) initial_size=oldr;


	if(source->getRadius() <= 0.0){ERROR_MESSAGE(); printf("ERROR: find_images, point source must have a resolution target\n"); exit(1);}

	// do an initial refinement to find all images and refine grid
	// the xmin is used as a characteristic small size for the source

	//if(FindCenter == false ){
/*
		// refine grid down to a small size
		Ntmp = NumberOfPoints(grid);
		if(initial_size >= 0) ImageFinding::find_images_kist(center,r_source,grid,Nimages
			  ,imageinfo,NimageMax,&Nimagepoints,initial_size,true,0,false,true);
		Ntmp = NumberOfPoints(grid) - Ntmp;
		if(Ntmp > 0) tmp = RefreshSurfaceBrightnesses(grid,lens);  // Calculate surface brightnesses for points added.
*/
/*	}else{
	// Find center and size of source based on points already found and telescope/refine to that source

		printf("total number of points before telescope: %li\n",NumberOfPoints(grid));

		PosType xsMax[2],xsMin[2],center[2],sbmax,ssize,r_source,rtemp,y[2];
		bool detected = false;

		xsMax[0] = xsMax[1] = -1.0e100;
		xsMin[0] = xsMin[1] = 1.0e100;
		sbmax = 0;
		ssize = 0.0;
		center[0] = source->getX()[0];
		center[1] = source->getX()[1];
		r_source = xmin;

		grid->ClearAllMarks();
	    rtemp = r_source*pow(Ngrid_block,(int)(log(initial_size/sqrt(pi)/fabs(r_source*mumin))/log(Ngrid_block) ) + 1);

		for( ; rtemp >= r_source*Ngrid_block ;rtemp /= Ngrid_block ){

			do{
				image_finder_kist(center,rtemp,grid
						  ,Nimages,imageinfo,NimageMax,&Nimagepoints,-1,0);
				detected = false;

				MoveToTopKist(imageinfo->imagekist);
				do{

					point = getCurrentKist(imageinfo->imagekist);

					//find surface brightnesses
					y[0] = point->image->x[0] - source->getX()[0];
					y[1] = point->image->x[1] - source->getX()[1];
					point->surface_brightness = source->SurfaceBrightness(y);
					point->image->surface_brightness  = point->surface_brightness;

					// Find new center and source size
					if(point->surface_brightness > 0){
						if(!detected){
							center[0] = center[1] = 0.0;
							tmp = 0.0;
							detected = true;
						}

						if(point->surface_brightness > 0.0){
							sbmax = point->surface_brightness;
							center[0] += point->image->x[0]*point->surface_brightness;
							center[1] += point->image->x[1]*point->surface_brightness;
							tmp += point->surface_brightness;
						}

						xsMax[0] = MAX(xsMax[0],point->image->x[0]);
						xsMin[0] = MIN(xsMin[0],point->image->x[0]);
						xsMax[1] = MAX(xsMax[1],point->image->x[1]);
						xsMin[1] = MIN(xsMin[1],point->image->x[1]);
					}

				}while(MoveDownKist(imageinfo->imagekist));

				if(detected){
					center[0] /= tmp;
					center[1] /= tmp;

					ssize = MIN(xsMax[0]-xsMin[0],xsMax[1]-xsMin[1])/2;
					if(ssize > 0) r_source = ssize;

					printf("     ssize = %e rtemp = %e r_source = %e center = %e %e Nimages = %i\n",ssize,rtemp,r_source,center[0],center[1],*Nimages);
				}

			}while(ImageFinding::refine_grid_kist(grid,imageinfo,*Nimages,rtemp*mumin/Ngrid_block,2,kappa_off,true,dummy_pnt));

			printf("      total number of points while telescoping: %li\n",NumberOfPoints(grid));

			while( (Ntmp = PrunePointsOutside(grid,rtemp*Ngrid_block,center,rtemp,rtemp*Ngrid_block)) )
				printf("Number of points pruned = %li\n",Ntmp);
		}

	}*/

	if(verbose) printf("total number of points after telescope: %li\n",grid->getNumberOfPoints());

	// Set all in_image flags to false.  This should not be necessary.  !!!!!!
	grid->ClearAllMarks();

	//freeKist(subkist);

	tmp = grid->RefreshSurfaceBrightnesses(source);

	if(tmp == 0.0){  // no flux was found
	  imageinfo[0].imagekist->Empty();
	  *Nimages = 0;
	  for(i=0; i < imageinfo.size() ; ++i) imageinfo[i].area = 0.0;
	  return;
	}
	assert(tmp > 0.0 || imageinfo[0].imagekist->Nunits() == 0);

	/*/********** test lines **********************
	PointsWithinKist_iter(grid->s_tree,source->getX(),0,source->source_r_out,imageinfo->imagekist);
	Ntmp = imageinfo->imagekist->Nunits();
	for(i=0,MoveToTopKist(imageinfo->imagekist); i < Ntmp ; ++i ){

		if(getCurrentKist(imageinfo->imagekist)->surface_brightness > 0){

			getCurrentKist(imageinfo->imagekist)->in_image = YES;
			MoveDownKist(imageinfo->imagekist);
		}else{
			getCurrentKist(imageinfo->imagekist)->in_image = NO;
			if(AtTopKist(imageinfo->imagekist)) go = false; else go = true;
			TakeOutCurrentKist(imageinfo->imagekist);
			if(go) MoveDownKist(imageinfo->imagekist);
		}
	}
	divide_images_kist(grid->s_tree,imageinfo,Nimages,NimageMax);
	printf("number of sources %i\n",*Nimages);
	/ *******************************************************/

	//PointsWithinKist(grid->s_tree,source->getX(),source->source_r_out,imageinfo->imagekist,0);
	grid->s_tree->PointsWithinKist_iter(source->getX(),0,source->getRadius(),imageinfo[0].imagekist);


	// move from source plane to image plane
	imageinfo[0].imagekist->TranformPlanes();

	// take out points with no flux
	if(verbose) printf("before taking out zero surface brightness points: %li\n",imageinfo[0].imagekist->Nunits());
	Ntmp = imageinfo[0].imagekist->Nunits();
	for(i=0,imageinfo[0].imagekist->MoveToTop(),area_tot=0,maxflux=0.0; i < Ntmp ; ++i ){

		flux = imageinfo[0].imagekist->getCurrent()->surface_brightness * pow(imageinfo[0].imagekist->getCurrent()->gridsize,2);
		maxflux = MAX(flux,maxflux);
		area_tot += flux;

		if(imageinfo[0].imagekist->getCurrent()->surface_brightness > 0){

			imageinfo[0].imagekist->getCurrent()->in_image = YES;
			imageinfo[0].imagekist->getCurrent()->image->in_image = YES;
			imageinfo[0].imagekist->Down();
		}else{
			imageinfo[0].imagekist->getCurrent()->in_image = NO;
			imageinfo[0].imagekist->getCurrent()->image->in_image = NO;
			if(imageinfo[0].imagekist->AtTop()) go = false; else go = true;
			imageinfo[0].imagekist->TakeOutCurrent();
			if(go) imageinfo[0].imagekist->Down();
		}

	}
	if(maxflux == 0 ) maxflux = 1.0;
	if(verbose) printf("after taking out zero surface brightness points: %li\n",imageinfo[0].imagekist->Nunits());

	if(imageinfo[0].imagekist->Nunits() == 0){
		*Nimages = 0;
		for(i=0; i < imageinfo.size() ; ++i) imageinfo[i].area = 0.0;
		return;
	}

	// divide up images
	if(divide_images) divide_images_kist(grid->i_tree,imageinfo,Nimages);
	else *Nimages = 1;
	if(verbose) printf("number of images after first division is %i\n",*Nimages);

	/////////////////////////////////////////////
	// link image points lists for each image
	// and calculate surface brightness at each point
	/////////////////////////////////////////////

	// ****** calculate surface brightnesses and flux of each image   ******
	for(i=0,area_tot=0.0; i < *Nimages ; ++i){

		imageinfo[i].ShouldNotRefine = 0;
		imageinfo[i].uniform_mag = unchecked;

		findborders4(grid->i_tree,&(imageinfo[i]));

		//assert(imageinfo[i].outerborder->Nunits() > 0);
		/*/  ***** test lines ****
		if(criterion != EachImage){
			MoveToTopKist(imageinfo[i].outerborder);
			do{
				if(getCurrentKist(imageinfo[i].outerborder)->surface_brightness > 0) point = getCurrentKist(imageinfo[i].outerborder);
				assert(getCurrentKist(imageinfo[i].outerborder)->surface_brightness == 0);

				if(getCurrentKist(imageinfo[i].outerborder)->id == 754939){
					ERROR_MESSAGE();  printf(" particle in border of image = %i  sb=%e  in_image=%i\n",i
							,getCurrentKist(imageinfo[i].outerborder)->surface_brightness
							,getCurrentKist(imageinfo[i].outerborder)->in_image);
				}

			}while(MoveDownKist(imageinfo[i].outerborder));

			if(getCurrentKist(imageinfo[i].imagekist)->id == 754939){
				ERROR_MESSAGE();  printf(" particle in image = %i  sb=%e  in_image=%i\n",i
						,getCurrentKist(imageinfo[i].imagekist)->surface_brightness
						,getCurrentKist(imageinfo[i].imagekist)->in_image);
			}
		}*/

		imageinfo[i].imagekist->MoveToTop();
		imageinfo[i].area = 0.0;
		do{
			imageinfo[i].area += pow(imageinfo[i].imagekist->getCurrent()->gridsize,2)
                             *(imageinfo[i].imagekist->getCurrent()->surface_brightness/maxflux);

		}while(imageinfo[i].imagekist->Down());

		area_tot += imageinfo[i].area;
		//printf("   %i area = %e\n",i,imageinfo[i].area);
		assert(imageinfo[i].area >= 0);
	}

	/*** need to make sure
	1) all images uniform_mag == unchecked initially
	2) uniform_mag get reset after process
	3) magnification is calculated properly.  Will need flux from source.  Careful of the normalization!
	4) What if there are more than one sources?
	***/
	/********************************************************
	 ******* refine images based on flux in each pixel ******
	 *******************************************************/
	i=0;
	while( ImageFinding::refine_grid_on_image(lens,source,grid,maxflux,imageinfo,Nimages
                                            ,sourceinfo,Nsources,FracResTarget
                                            ,criterion,divide_images) > 0 ) ++i;

	//printf("i=%i Nold=%li\n",i,Nold);
	//printf("%li\n",imagelist->Npoints);

	oldy[0] = source->getX()[0];
	oldy[1] = source->getX()[1];
	oldr = source->getRadius();

	oldNimages=*Nimages;

	// find image centroid
	for(i=0;i<*Nimages;++i){
		// reset these flags
		imageinfo[i].ShouldNotRefine = 0;
		imageinfo[i].uniform_mag = unchecked;
/*
		// Find borders to images
		for(j = 0,MoveToTopKist(imageinfo[i].imagekist) ; j < imageinfo[i].imagekist->Nunits() ; ++j,MoveDownKist(imageinfo[i].imagekist) ){
			getCurrentKist(imageinfo[i].imagekist)->in_image = YES;  // re-set marks
			getCurrentKist(imageinfo[i].imagekist)->image->in_image = YES;  // re-set marks
		}
		findborders4(grid->i_tree,&(imageinfo[i]));
*/
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

/** \ingroup ImageFinding
 *  \brief Find the images without any additional grid refinement or ray shooting.
 *
 *   The images and surface brighnesses are found.  This will be much faster than map_images(), 
 *   but will miss images that are not resolved on the current grid.
 *   
 *   The surface brightnesses on all of the grid are refreshed.
 *   The images centroids and gridranges are maintained.
 **/

void ImageFinding::map_images_fixedgrid(
                Source *source
                ,GridHndl grid          /// Tree of grid points
                ,int *Nimages           /// number of images found
                ,std::vector<ImageInfo> &imageinfo   /// information on each image
                ,PosType xmax            /// Maximum size of source on image plane.  The entire source must be within this distance from
                                        ///  source->getX()[].  Decreasing it will make the code run faster.  Making xmax much bigger than
                                        /// the grid boundaries will check all points for surface brightness.
                ,bool divide_images     /// if true will divide images.
                ,bool find_borders      /// if true will find the inner and outer borders of each image
  ){
  
  if(imageinfo.size() < 2) imageinfo.resize(2);

	assert(grid->s_tree);
	assert(grid->i_tree);
	assert(imageinfo[0].imagekist);
  
  imageinfo[0].area = 0.0;
  imageinfo[0].area_error = 0.0;

  imageinfo[0].innerborder->Empty();
  imageinfo[0].outerborder->Empty();
  
  grid->RefreshSurfaceBrightnesses(source);
  grid->ClearAllMarks();  // TODO: might be nice to elliminate the need for this

  PosType x[2];
  
  grid->s_tree->PointsWithinKist(source->getX(), xmax, imageinfo[0].imagekist, 0);
  
  bool move;
  // Assign surface brightnesses and remove points from image without flux
  for(imageinfo[0].imagekist->MoveToTop();!(imageinfo[0].imagekist->OffBottom());){
    
    x[0] = imageinfo[0].imagekist->getCurrent()->x[0] - source->getX()[0];
    x[1] = imageinfo[0].imagekist->getCurrent()->x[1] - source->getX()[1];
    
    imageinfo[0].imagekist->getCurrent()->surface_brightness
    = imageinfo[0].imagekist->getCurrent()->image->surface_brightness
      //= source->SurfaceBrightness( x );
    = source->SurfaceBrightness(imageinfo[0].imagekist->getCurrent()->x);
    
    imageinfo[0].area += imageinfo[0].imagekist->getCurrent()->surface_brightness;
    
    if(imageinfo[0].imagekist->getCurrent()->surface_brightness == 0.0){
      imageinfo[0].imagekist->getCurrent()->in_image = NO;  // re-set marks
			imageinfo[0].imagekist->getCurrent()->image->in_image = NO;  // re-set marks
      move = imageinfo[0].imagekist->AtTop();
      imageinfo[0].imagekist->TakeOutCurrent();
      if(!move) imageinfo[0].imagekist->Down();
    }else{
      imageinfo[0].imagekist->getCurrent()->in_image = YES;  // re-set marks
			imageinfo[0].imagekist->getCurrent()->image->in_image = YES;  // re-set marks
      imageinfo[0].imagekist->Down();
    }
    
  }
  
  if(imageinfo[0].imagekist->Nunits() == 0){
    *Nimages = 0;
    return;
  }
  
  *Nimages = 1;
  imageinfo[0].imagekist->TranformPlanes();
  
  int i;

  // divide up images
  if(divide_images) divide_images_kist(grid->i_tree,imageinfo,Nimages);
 
  // calculate grid range for images
  for(i=0;i<*Nimages;++i){
    
    imageinfo[i].gridrange[0] = imageinfo[0].gridrange[1] = 0.0;
    imageinfo[i].gridrange[2] = 1.0e99; // minimum grid size in image
    
    for(imageinfo[i].imagekist->MoveToTop();!(imageinfo[i].imagekist->OffBottom());imageinfo[i].imagekist->Down()){
      if(imageinfo[i].imagekist->getCurrent()->gridsize < imageinfo[i].gridrange[2])
        imageinfo[i].gridrange[2] = imageinfo[i].imagekist->getCurrent()->gridsize;
      if(imageinfo[i].imagekist->getCurrent()->gridsize > imageinfo[i].gridrange[0])
        imageinfo[i].gridrange[0] = imageinfo[i].imagekist->getCurrent()->gridsize;
    }
  }

  if(find_borders){
    for(i=0;i<*Nimages;++i) findborders4(grid->i_tree,&(imageinfo[i]));
  }

  // find image centroids
  for(i=0;i<*Nimages;++i){
		imageinfo[i].imagekist->MoveToTop();
		imageinfo[i].centroid[0] = 0.0;
		imageinfo[i].centroid[1] = 0.0;
    imageinfo[i].area = 0.0;
		for(long j = 0 ; j < imageinfo[i].imagekist->Nunits() ; ++j,imageinfo[i].imagekist->Down() ){
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
		}
		if(imageinfo[i].getNimagePoints() > 0 ){
			imageinfo[i].centroid[0] /= imageinfo[i].area;
			imageinfo[i].centroid[1] /= imageinfo[i].area;
		}
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
int ImageFinding::refine_grid_on_image(Lens *lens,Source *source,GridHndl grid,PosType maxflux
                                       ,std::vector<ImageInfo> &imageinfo,int *Nimages
                                       ,std::vector<ImageInfo> &sourceinfo
                                       ,int Nsources,const PosType res_target
                                       ,ExitCriterion criterion,bool divide_images,bool batch){

	//printf("entering refine_grid\n");

  if((*Nimages) < 1) return 0;

  int k,number_of_refined; /* Ngrid_block must be odd */
  PosType total_area,r,rmin;
  Point *i_points;
  unsigned long Ncells,Nold,j,i;
  bool reborder=false,redivide=false;

  for(i=0,total_area=0;i<(*Nimages);++i) total_area += imageinfo[i].area;
  for(i=0,Nold=0;i<(*Nimages);++i) Nold += imageinfo[i].imagekist->Nunits();
  if(total_area == 0.0) return 0;

  Kist<Point> * nearest=new Kist<Point>;
  std::vector<Point *> points_to_refine;

  if(verbose) printf(" entering refine_on_image(), number of points %li\n",grid->getNumberOfPoints());

  number_of_refined=0;
  for(i=0,Ncells=0;i<(*Nimages);++i){

	  // Constant magnification cut
	  if(imageinfo[i].uniform_mag == unchecked){
		  ImageFinding::UniformMagCheck(&imageinfo[i]);
		  if(imageinfo[i].uniform_mag == yes){
			  if(Nsources > 1){
				  // Find the closest source
				  for(j=0,rmin = 1.0e100;j<Nsources;++j){
					  r = pow(sourceinfo[j].centroid[0] - imageinfo[i].imagekist->getCurrent()->image->x[0],2)
						  + pow(sourceinfo[j].centroid[1] - imageinfo[i].imagekist->getCurrent()->image->x[1],2);
					  if(rmin > r){
						  rmin = r;
						  k = j;
					  }
				  }
			  }else{ k = 0; }

			  assert(sourceinfo[k].area > 0);
			  assert(fabs(imageinfo[i].imagekist->getCurrent()->invmag) > 0.0);
			  imageinfo[i].area = sourceinfo[k].area/fabs(imageinfo[i].imagekist->getCurrent()->invmag)/maxflux;
			  if(verbose) printf("magnification is uniform for image %lu\n",i);
		  }
	  }

	  assert(imageinfo[i].area >= 0.0);
	  if(imageinfo[i].ShouldNotRefine == 0 && imageinfo[i].uniform_mag == no){  // If image was not refined on the last round do not revisit unless the images are re-divided.
		  imageinfo[i].ShouldNotRefine = 1;

		  imageinfo[i].imagekist->MoveToTop();
		  for(j = 0 ; j < imageinfo[i].imagekist->Nunits() ; ++j,imageinfo[i].imagekist->Down() ){

			  if(
					  ImageFinding::RefinePoint2(imageinfo[i].imagekist->getCurrent(),grid->i_tree
					  	,imageinfo[i].area,total_area,maxflux,criterion,res_target,nearest)
           //RefinePoint_sb(getCurrentKist(imageinfo[i].imagekist),grid->i_tree
           //		,imageinfo[i].area,total_area,2*source->getSBlimit(),nearest)
           //RefinePoint_smallsize(getCurrentKist(imageinfo[i].imagekist),grid->i_tree
           //		,imageinfo[i].area,total_area,source->,nearest)

			  ){

				  ++Ncells;

				  // adjust contribution to image flux from point that will be refined.
				  assert(imageinfo[i].area > 0.0);

				  // Determine if image point to be refined is a inner border point
				  if(reborder == false){
					  if(nearest->Nunits() == 0 ) grid->i_tree->FindAllBoxNeighborsKist(imageinfo[i].imagekist->getCurrent(),nearest);
					  nearest->MoveToTop();
					  do{
						  if(nearest->getCurrent()->in_image != YES){
							  reborder = true;
							  break;
						  }
					  }while(nearest->Down());
				  }

				  imageinfo[i].area -=  pow(imageinfo[i].imagekist->getCurrent()->gridsize,2)
						  *(imageinfo[i].imagekist->getCurrent()->surface_brightness/maxflux);

				  assert(imageinfo[i].area >= 0);

				  assert(imageinfo[i].imagekist->getCurrent()->leaf->child1 == NULL);
				  assert(imageinfo[i].imagekist->getCurrent()->leaf->child2 == NULL);

				  if(batch){
					  points_to_refine.push_back(imageinfo[i].imagekist->getCurrent());

					  imageinfo[i].area +=  pow(imageinfo[i].imagekist->getCurrent()->gridsize/grid->getNgrid_block(),2)
							  *(imageinfo[i].imagekist->getCurrent()->surface_brightness/maxflux);
				  }else{
					  i_points = grid->RefineLeaf(lens,imageinfo[i].imagekist->getCurrent());
					  ImageFinding::check_sb_add(source,&imageinfo[i],i_points,maxflux,Nold,number_of_refined);

					  /*if(i_points != NULL){
						  // link new points into image kist and calculate surface brightnesses
						  for(k=0;k < i_points->head;++k){

							  // put point into image imageinfo[i].imagekist

							  //y[0] = i_points[k].image->x[0];// - source->getX()[0];
							  //y[1] = i_points[k].image->x[1];// - source->getX()[1];
							  i_points[k].surface_brightness = source->SurfaceBrightness(i_points[k].image->x);
							  i_points[k].image->surface_brightness  = i_points[k].surface_brightness;

							  // if new point has flux add to image
							  if(i_points[k].surface_brightness > 0.0){
								  InsertAfterCurrentKist(imageinfo[i].imagekist,&(i_points[k]));
								  MoveDownKist(imageinfo[i].imagekist);

								  i_points[k].in_image = YES;
								  i_points[k].image->in_image = YES;

								  imageinfo[i].area += pow(i_points[k].gridsize,2)*(i_points[k].surface_brightness/maxflux);

							  }else{
								  i_points[k].in_image =  i_points[k].image->in_image = NO;
							  }

							  ++number_of_refined;
						  }
						  Nold += i_points->head;
					  }*/

					  imageinfo[i].area +=  pow(imageinfo[i].imagekist->getCurrent()->gridsize,2)
							  *(imageinfo[i].imagekist->getCurrent()->surface_brightness/maxflux);
				  }

				  imageinfo[i].ShouldNotRefine = 0;   // mark to continue refinement on next round

			  }

		  }

		  if(batch){
			  i_points = grid->RefineLeaves(lens,points_to_refine);
			  ImageFinding::check_sb_add(source,&imageinfo[i],i_points,maxflux,Nold,number_of_refined);
			  points_to_refine.clear();
		  }

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
		  for(j = 0 ; j < imageinfo[i].outerborder->Nunits() ; ++j,imageinfo[i].outerborder->Down() ){

			  // TODO: This was taken out and i'm not sure if it was needed.
			  //assert(getCurrentKist(imageinfo[i].outerborder)->surface_brightness == 0);

			  if(
					  ImageFinding::RefinePoint2(imageinfo[i].outerborder->getCurrent(),grid->i_tree
					      ,imageinfo[i].area,total_area,maxflux,criterion,res_target,nearest)
					  //RefinePoint_sb(getCurrentKist(imageinfo[i].outerborder),grid->i_tree
					  //    ,imageinfo[i].area,total_area,2*source->getSBlimit(),nearest)
			  ){

				  if(imageinfo[i].outerborder->getCurrent()->in_image != MAYBE){
					  imageinfo[i].outerborder->getCurrent()->in_image = MAYBE;

					  reborder = true;  // Since the border has been refined, re do the borders.
					  ++Ncells;

					  assert(imageinfo[i].outerborder->getCurrent()->leaf->child1 == NULL);
					  assert(imageinfo[i].outerborder->getCurrent()->leaf->child2 == NULL);
					  if(batch){
						  points_to_refine.push_back(imageinfo[i].outerborder->getCurrent());
					  }else{
						  i_points = grid->RefineLeaf(lens,imageinfo[i].outerborder->getCurrent());
						  ImageFinding::check_sb_add(source,&imageinfo[i],i_points,maxflux,Nold,number_of_refined);
					  }
					  /*if(i_points != NULL){
						  // link new points into image kist and calculate surface brightnesses
						  for(k=0;k < i_points->head ;++k){

							  // put point into image imageinfo[i].outerborder
							  //y[0] = i_points[k].image->x[0] - source->getX()[0];
							  //y[1] = i_points[k].image->x[1] - source->getX()[1];
							  i_points[k].surface_brightness = source->SurfaceBrightness(i_points[k].image->x);
							  i_points[k].image->surface_brightness  = i_points[k].surface_brightness;

							  // if new point has flux add to image
							  if(i_points[k].surface_brightness > 0.0){
								  InsertAfterCurrentKist(imageinfo[i].imagekist,&(i_points[k]));
								  MoveDownKist(imageinfo[i].imagekist);

								  i_points[k].in_image = YES;
								  i_points[k].image->in_image = YES;
								  imageinfo[i].area += pow(i_points[k].gridsize,2)*(i_points[k].surface_brightness/maxflux);

							  }else{
								  i_points[k].in_image =  i_points[k].image->in_image = NO;
							  }

							  ++number_of_refined;
						  }
						  Nold += i_points->head;
					  }*/

					  imageinfo[i].ShouldNotRefine = 0;   // mark for another look next time

				  }else{
					/// there is an overlap in images so the images will need to be re-divided
					  redivide = true;
				  }

			  }
		  } // loop through outer border

		  assert(imageinfo[i].area >= 0.0);

		  if(batch){
			  i_points = grid->RefineLeaves(lens,points_to_refine);
			  ImageFinding::check_sb_add(source,&imageinfo[i],i_points,maxflux,Nold,number_of_refined);
			  points_to_refine.clear();
		  }

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
	  } // if previously refined
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

			ImageFinding::UniformMagCheck(&imageinfo[i]);
			if(imageinfo[i].uniform_mag == yes){
				if(Nsources > 1){
					// Find the closest source
					for(j=0,rmin = 1.0e200;j<Nsources;++j){
						r = pow(sourceinfo[j].centroid[0] - imageinfo[i].imagekist->getCurrent()->image->x[0],2)
									  + pow(sourceinfo[j].centroid[1] - imageinfo[i].imagekist->getCurrent()->image->x[1],2);
						if(rmin > r){
							rmin = r;
							k = j;
						}
					}
				}else{ k=0; }

				assert(sourceinfo[k].area > 0);
				assert(fabs(imageinfo[i].imagekist->getCurrent()->invmag) > 0);
				imageinfo[i].area = sourceinfo[k].area/fabs(imageinfo[i].imagekist->getCurrent()->invmag)/maxflux;

			}else{
				imageinfo[i].area = 0.0;
				imageinfo[i].imagekist->MoveToTop();
				do{
					imageinfo[i].area += pow(imageinfo[i].imagekist->getCurrent()->gridsize,2)
		    				 * (imageinfo[i].imagekist->getCurrent()->surface_brightness/maxflux);
				}while(imageinfo[i].imagekist->Down());
			}
		}
  }

  return number_of_refined;
}
// Assign surface brightness of new points and add the ones that are nonzero to the image
void ImageFinding::check_sb_add(Source *source,ImageInfo *imageinfo,Point *i_points,PosType maxflux,unsigned long &Nold,int &number_of_refined){
	if(i_points != NULL){
		// link new points into image kist and calculate surface brightnesses
		for(unsigned long k=0;k < i_points->head;++k){

			// put point into image imageinfo[i].imagekist

			//y[0] = i_points[k].image->x[0];// - source->getX()[0];
			//y[1] = i_points[k].image->x[1];// - source->getX()[1];
			i_points[k].surface_brightness = source->SurfaceBrightness(i_points[k].image->x);
			i_points[k].image->surface_brightness  = i_points[k].surface_brightness;

			// if new point has flux add to image
			if(i_points[k].surface_brightness > 0.0){
				imageinfo[0].imagekist->InsertAfterCurrent(&(i_points[k]));
				imageinfo[0].imagekist->Down();

				i_points[k].in_image = YES;
				i_points[k].image->in_image = YES;

				imageinfo[0].area += pow(i_points[k].gridsize,2)*(i_points[k].surface_brightness/maxflux);

			}else{
				i_points[k].in_image =  i_points[k].image->in_image = NO;
			}

			++number_of_refined;
		}
		Nold += i_points->head;
	}
}

bool ImageFinding::RefinePoint2(Point *point,TreeHndl i_tree,PosType image_area,PosType total_area,PosType maxflux
		,ExitCriterion criterion,PosType res_target,Kist<Point> * nearest){

	PosType borderSB = 0,error = 0,maxdiff_sb;

	nearest->Empty();
	// Prevent cell from getting so small that precision error prevents everything from working
	if(point->gridsize <= pow(10.,1 - DBL_DIG)) return false;  // this shouldn't be necessary every time

	if( image_area < target_all*res_target*total_area ) return false;

	if(criterion == FillHoles)
		if(pow(point->gridsize,2)*(point->surface_brightness/maxflux) > res_target*image_area) return true;

	// This stops the refinement when the cell is so small that the flux becomes zero because of underflow error.
	//if(point->surface_brightness > 0 && pow(point->gridsize,2)*point->surface_brightness == 0) return false;

	i_tree->FindAllBoxNeighborsKist(point,nearest);
	//assert(nearest->Nunits() < 300);

	maxdiff_sb = 0.0;
	nearest->MoveToTop();
	do{
		borderSB += nearest->getCurrent()->surface_brightness;
		maxdiff_sb = MAX(fabs(nearest->getCurrent()->surface_brightness - point->surface_brightness),maxdiff_sb);
	}while(nearest->Down());
	borderSB /= nearest->Nunits();

	error = pow(point->gridsize,2)*(fabs(borderSB - point->surface_brightness)/6/maxflux);
	//error = pow(point->gridsize,2)*maxdiff/6/maxflux;

	if( ( criterion == EachImage || criterion == FillHoles ) && error > res_target*image_area ) return true;
	if( criterion == TotalArea && error > res_target*total_area ) return true;

	return false;
}
// refinement criterion based on difference in surface brightness
bool ImageFinding::RefinePoint_sb(Point *point,TreeHndl i_tree,PosType image_area,PosType total_area
                    ,PosType sb_limit,PosType maxflux,Kist<Point> * nearest){
  
  PosType smallsize = 1.0e-7;
	// Prevent cell from getting so small that precision error prevents everything from working
	if(point->gridsize <= pow(10.,1 - DBL_DIG)) return false;  // this shouldn't be necessary every time
  
	if( image_area < 1.0e-3*total_area ) return false;
	if( point->gridsize*point->gridsize*point->surface_brightness/maxflux > 1.0e-4*image_area)
		return false;
  
  if(point->gridsize*point->gridsize*point->invmag > smallsize*smallsize) return true;
  else return false;
  
  //*****************************************************************
	nearest->Empty();
  
	i_tree->FindAllBoxNeighborsKist(point,nearest);
	nearest->MoveToTop();
	do{
		if(sb_limit < fabs(nearest->getCurrent()->surface_brightness - point->surface_brightness) )
			return true;
	}while(nearest->Down());
  
	return false;
}
// refinement criterion based on difference in surface brightness
bool ImageFinding::RefinePoint_smallsize(Point *point,TreeHndl i_tree,PosType image_area,PosType total_area
                    ,PosType smallsize,PosType maxflux,Kist<Point> * nearest){
  
	// Prevent cell from getting so small that precision error prevents everything from working
	if(point->gridsize <= pow(10.,1 - DBL_DIG)) return false;  // this shouldn't be necessary every time
  
	if( image_area < 1.0e-3*total_area ) return false;
	if( point->gridsize*point->gridsize*point->surface_brightness/maxflux > 1.0e-4*image_area)
		return false;
  
  if(point->gridsize*point->gridsize*point->invmag > smallsize*smallsize) return true;
  
  return false;
}

/*
 * Requires that all neighbors to point have a flux below a
 *
bool RefinePoint(Point *point,TreeHndl i_tree,PosType image_area,PosType total_area
		,ExitCriterion criterion,PosType res_target,Kist<Point> * nearest){

	PosType flux,tmp;

	flux = pow(point->gridsize,2)*point->surface_brightness;

	if(criterion == TotalArea && flux > res_target*total_area){
		return true;
	}

	if(criterion == EachImage && (flux > res_target*image_area
		&& flux > target_all*res_target*total_area) ){
		return true;
	}

	if(criterion == Resolution &&
			point->gridsize > res_target ){
		return true;
	}

	if(criterion == FillHoles ){

		if(res_target*total_area < image_area){
			if(flux > res_target*image_area
				&& flux > target_all*res_target*total_area
			){
				return true;
			}else{  // refine the borders by refining the current if any of its neighbors have enough flux

				FindAllBoxNeighborsKist(i_tree,point,nearest);

				MoveToTopKist(nearest);
				do{
					tmp = getCurrentKist(nearest)->surface_brightness*pow(getCurrentKist(nearest)->gridsize,2);


					if( tmp > res_target*image_area
							&& tmp > target_all*res_target*total_area){

						return true;
					}
				}while( MoveDownKist(nearest) );

			assert(nearest->Nunits() < 300);
			}
		}
	}
	return false;
}
*/
/** \ingroup mageFindingL2
 * \brief Checks to see if the image has a nearly uniform magnification across it and thus can be considered linearly
 * distorted.
 */
void ImageFinding::UniformMagCheck(ImageInfo *imageinfo){

	// find minimum and maximum magnification on border
	if(imageinfo->imagekist->Nunits() > 10 && imageinfo->uniform_mag == unchecked){

		PosType magmin,magmax;

		imageinfo->imagekist->MoveToBottom();
		magmin = magmax = 1.0/imageinfo->imagekist->getCurrent()->invmag;
		imageinfo->imagekist->MoveToTop();
		do{
			magmin = MIN(magmin,1.0/imageinfo->imagekist->getCurrent()->invmag);
			magmax = MAX(magmax,1.0/imageinfo->imagekist->getCurrent()->invmag);
		}while( imageinfo->imagekist->Down() && (magmax-magmin) < tol_UniformMag*fabs(magmax) );

		if((magmax-magmin) < tol_UniformMag*fabs(magmax)){
			imageinfo->uniform_mag = yes;
		}else{
			imageinfo->uniform_mag = no;
		}
	}
}


