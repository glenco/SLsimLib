
#include "slsimlib.h"

static const int NpointsRequired = 100;  // number of points required to be within an image
//static const int Ngrid_block = 3;       // each cell is divided into Ngrid_block^2 subcells
//static const float mumin = 0.5;  // actually the sqrt of the minimum magnification
//static const float mumin = 0.45;  // actually the sqrt of the minimum magnification
//static const float mumin = 0.1;
static const float mumin = 0.3;


static const float FracResTarget = 4.0e-4;
//static const float FracResTarget = 1.0e-4;
static const float telescope_high = 1.0e-3;
static const float telescope_low = 0.01;
extern const double initialgridsize;

/** \ingroup ImageFinding
 *
 * \brief  Finds images given a source position and size.
 *
 * find_image_kist returns finite refined images in images[0...*Nimages-1].imagekist
 *  It starts with a large source and reduces down to the right size refining at each step.
 *  It should not miss any image larger than ~ munin*r_source linear size, but seems to do
 *  much better than that.  It does nothing with the surface brightnesses.
 *
 * The routine can follow three different strategies for refining each image controlled by edge_refinement.
 *
 * edge_refinement
 *   - 0 does not do edge refinement, Every pixel in every image is refined until the criterion is met.
 *  The image(s) are found again after each refinement which can make it slower.
 *   - 1 uses refine_edge().  After an initial refinement of all the pixels in the image(s) the code switches
 *	 to refining only the edges of the images.  The images are found after each refinement.
 *   - 2 uses refine_edge2() Same as 1, but the images are not found after each refinement.  This can make the
 *	 routine run much faster, but has the disadvantage that the number of images will not change during the final
 *	 stage of refinement.  This is the setting generally recommended.
 *
 *  kappa_off - This is used to turn off the calculation of surface density, shear, magnification and time delay.
 *  When finding a finite sized source these quantities are generally not required and slow down the routine.
 *
 *
 */

void find_images_kist(
		LensHndl lens,          /// contains the lens/es and source/sources
		double *y_source        /// position of source center
		,double r_source        /// radius of source
		,GridHndl grid          /// grid provided to routine
		,int *Nimages           /// number of images found
		,ImageInfo *imageinfo   /// information on each image
		,const int NimageMax    /// maximum number of images allowed
		,unsigned long *Nimagepoints  /// number of points in final images
		,double initial_size    /// Initial size of source for telescoping, 0 to start from the initial grid size.
		,bool splitimages       /// TRUE each image is refined to target accuracy, otherwise all images are treated as one
		,short edge_refinement  /// see comment
		,bool verbose           /// verbose
		,bool kappa_off         /// turns off calculation of surface density, shear, magnification and time delay
		){


	if(  grid->s_tree->top->boundary_p1[0] > (y_source[0] + r_source)
	  || grid->s_tree->top->boundary_p2[0] < (y_source[0] - r_source)
	  || grid->s_tree->top->boundary_p1[1] > (y_source[1] + r_source)
	  || grid->s_tree->top->boundary_p2[1] < (y_source[1] - r_source)
	){
		// source is not within initialized grid
		*Nimages = 0;
		std::cout << "source not within initialized grid" << std::endl;
		ERROR_MESSAGE();
		return;
	}

	int Nsizes;
	double rtemp,tmp,maxgridsize;
	static double oldy[2],oldr=0;
	short moved,flag;
	int i,j,k;
	//Point *i_points,*s_points,*point;
	time_t to,t1,t2,t3,now;
	//Kist<Point> * tmp_border_kist;
	bool image_overlap;
	static int oldNimages=0;
	static unsigned long Npoints_old = 0;
	//Point **dummy_pnt = NULL;
	//unsigned long Ntmp;
	//Point *point,*closestpoint;


	int Ngrid_block = grid->getNgrid_block();

	if(r_source==0.0){ERROR_MESSAGE(); printf("ERROR: find_images, point source must have a resolution target\n"); exit(1);}

	if(verbose) printf("initialgridsize=%e\n",initialgridsize);
	if(initial_size==0) initial_size=initialgridsize;

	if(oldr==0){ oldr=r_source; Npoints_old = grid->i_tree->pointlist->Npoints;}
	if((Npoints_old <= grid->i_tree->pointlist->Npoints )* // if grid has not been refreshed
			(oldy[0]==y_source[0])*(oldy[1]==y_source[1])* // and source not moved
			(oldr > r_source)  // and source size has gotten smaller
	){
		Nsizes=(int)( log(oldr/r_source/mumin)/log(Ngrid_block) ); // round up
	    rtemp = r_source*pow(1.0*Ngrid_block,Nsizes);
	}else{
		Nsizes=(int)(log(initial_size/fabs(r_source*mumin))/log(Ngrid_block) ) + 1 ; // round up
	    rtemp = r_source*pow(1.0*Ngrid_block,Nsizes);
	}

	Npoints_old = grid->i_tree->pointlist->Npoints;


	// starting with a larger source size make sure all the grid sizes are small enough to find it
	Kist<Point> * subkist = new Kist<Point>;//,pointkist = NewKist();

	if(verbose) printf("entering find_image\n");
	time(&to);


    if(verbose) printf("Ntemp=%li\n",Nsizes);

    grid->ClearAllMarks();

    //////////////////////////////////////////
    // telescope source size down to target
    //////////////////////////////////////////

    int Nimages_old = 1;
    for(i=0
				//for(rtemp = fabs(r_source/mumin)*pow(Ngrid_block,Nsizes),Nold=0
		//		;rtemp >= 0.99*Ngrid_block*fabs(r_source)
    		;rtemp >= r_source
       		;rtemp /= Ngrid_block,++i ){

    	time(&t3);

 			/************* method that seporates images ****************
			// mark image points in tree
			PointsWithinKist(grid->s_tree,y_source,rtemp,subkist,1);
			moved = image_finder_kist(lens,y_source,rtemp,grid
				  ,Nimages,imageinfo,NimageMax,Nimagepoints,0,0);
			// unmark image points in tree
			PointsWithinKist(grid->s_tree,y_source,rtemp,subkist,-1);
			//***********************************************************/


			time(&t1);
			time(&t2);
			if(verbose)
				printf("\n   new source size = %e    Nimages = %i  telescoping rsource = %e\n",rtemp,*Nimages,r_source);

			j=0;
			//while(refine_grid_kist(lens,grid,imageinfo,*Nimages,telescope_res,3,kappa_off,NULL)){
			//while(refine_grid_kist(lens,grid,imageinfo,*Nimages,0.05/Ngrid_block/Ngrid_block,3,kappa_off,NULL)){
			//while(refine_grid_kist(lens,grid,imageinfo,*Nimages,0.05/Ngrid_block/Ngrid_block,1,kappa_off,NULL)){
			do{
				time(&t1);
				if(verbose) std::cout << "    refined images" << std::endl;

				//************* method that does not separate images ****************
				moved = image_finder_kist(lens,y_source,rtemp,grid
						,Nimages,imageinfo,NimageMax,Nimagepoints,-1,0);
				/************* method that separates images ****************/

				/************* method that separates images ****************
				// mark image points in tree
				PointsWithinKist(grid->s_tree,y_source,rtemp,subkist,1);
				moved = image_finder_kist(lens,y_source,rtemp,grid
						,Nimages,imageinfo,NimageMax,Nimagepoints,0,0);
				// unmark image points in tree
				PointsWithinKist(grid->s_tree,y_source,rtemp,subkist,-1);
				//***********************************************************/

				assert(*Nimages > 0);

				if(verbose){
					printf("      refound images after refinement\n        Nimagepoints=%li  Nimages = %i\n"
						,*Nimagepoints,*Nimages);
					for(int k=0; k < *Nimages; ++k){
							std::cout << "   " << imageinfo[k].area << " " << imageinfo[k].area_error << " " << imageinfo[k].area/pi/rtemp/rtemp
									<< " " << imageinfo[k].getNimagePoints() << std::endl;
					}
				}

				++j;
			}while(refine_grid_kist(lens,grid,imageinfo,*Nimages,rtemp*mumin/Ngrid_block,2,kappa_off,NULL));
			//}

			time(&t1);
			if(verbose)	printf("      time in refine grid %f sec\n",difftime(t1,t2));

			time(&now);
			if(verbose) printf("    time for one source size %f sec\n",difftime(now,t3));
    }

	time(&now);
	if(verbose) printf(" time for source size reduction %f sec\n",difftime(now,to));
	time(&to);

	////////////////////////////////////////////////////////////////////////////////*/
	//////////////////////////////////////////////////////////////////////////////////
	// target source size has been reached, do full image decomposition and refinement
	/////////////////////////////////////////////////////////////////////////////////

    grid->ClearAllMarks();

	i=0;

	if(splitimages) flag = 1; else flag = 0;
	//flag = 1; // Changed so that each image always has at least 100 points.

	time(&now);

	assert(*Nimages > 0);
	// do an initial uniform refinement to make sure there are enough point in
	//  the images
	i=0;
	do{
		time(&t3);
		if(verbose)
			printf("     time in image refinement %f min\n           points in grid=%li\n"
					,fabs(difftime(t3,now)/60.),grid->i_tree->pointlist->Npoints);

		// mark image points in tree
		grid->s_tree->PointsWithinKist(y_source,r_source,subkist,1);

		//moved=image_finder_kist(lens,y_source,fabs(r_source),grid
		//		,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);
		moved=image_finder_kist(lens,y_source,fabs(r_source),grid
				,Nimages,imageinfo,NimageMax,Nimagepoints,0,0);

		//if(*Nimages < 1) printf("  Nimages=%i i=%i\n",*Nimages,i);

		time(&now);
		if(verbose){
			printf("\n    i=%i\n     time in finding images %f min\n          Nimages=%i   Nimagepoints=%li\n"
					,i,difftime(now,t3)/60.,*Nimages,*Nimagepoints);
			printf("     image   # of points    error in area\n");
			for(j=0;j<*Nimages;++j) printf("       %i        %li         %e\n",j,imageinfo[j].imagekist->Nunits(),imageinfo[j].area_error);
		}
		if(i > 20 && *Nimagepoints == 100){
			// case where no image is found at any size
			*Nimages = 0;
			*Nimagepoints = 0;
			return ;
		}
		++i;
	}while( refine_grid_kist(lens,grid,imageinfo,*Nimages,1.0/NpointsRequired,flag,kappa_off,NULL)
			|| moved );
	assert(*Nimages > 0);

	// find points that are truly in the image and not just neighbors
	moved=image_finder_kist(lens,y_source,fabs(r_source),grid
			,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);
	assert(*Nimages > 0);


	// remove images with no points in them
	for(j=0;j<*Nimages;++j){
		if(imageinfo[j].imagekist->Nunits() < 1){
			ERROR_MESSAGE();
			for(k=j+1;k<*Nimages;++k) SwapImages(&imageinfo[k-1],&imageinfo[k]);
			//printf("image %i has no points\n",j);
			--*Nimages;
			--j;
		}
	}
	assert(*Nimages > 0);

	/////////////////////////////////////////////
	// second stage of refinement -
	// depends of choice of edge_refinement
	/////////////////////////////////////////////
	time(&now);

	if(splitimages) flag = 0; else flag = 2;

	k=i;
	if(edge_refinement==0){   // uniform refinement over image
		do{
			// mark image points in tree
			grid->s_tree->PointsWithinKist(y_source,r_source,subkist,1);

			moved=image_finder_kist(lens,y_source,fabs(r_source),grid
					,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);
			++i;
		}while( refine_grid_kist(lens,grid,imageinfo,*Nimages,FracResTarget,0,kappa_off,NULL)
				|| moved );

	}else if(edge_refinement==1){    // edge refinement with image finding at each step
		do{
			// mark image points in tree
			grid->s_tree->PointsWithinKist(y_source,r_source,subkist,1);

			moved=image_finder_kist(lens,y_source,fabs(r_source),grid
					,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);

			//for(i = 0; i < *Nimages; ++i) PrintImageInfo(&(imageinfo[i]));
			//printf("\n");

			++i;
		}while( refine_edges(lens,grid,imageinfo,*Nimages,FracResTarget,flag,kappa_off)
				|| moved );

	}else if(edge_refinement==2){  // edge refinement with no image finding at each step
		++i;
		while(refine_edges2(lens,y_source,r_source,grid
				,imageinfo,&image_overlap,*Nimages,FracResTarget,flag,kappa_off)){
			// if an overlap is detected find the images again

			if(image_overlap) moved=image_finder_kist(lens,y_source,fabs(r_source),grid
					,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);
			++i;
		}
	}
	// unmark image points so new source can be used
	grid->s_tree->PointsWithinKist(y_source,r_source,subkist,-1);

	if(verbose) printf("finished edge refinement i=%i\n",i);

	assert(*Nimages > 0);

	time(&t3);
	if(verbose) printf("     time in image refinement %f min\n",difftime(t3,now)/60.);

	// if point source take only closest image point
	if(r_source <= 0){
		ERROR_MESSAGE();
		exit(1);
	}

	time(&now);
	if(verbose) printf("time in find_images %f min\n",difftime(now,to)/60.);

	oldy[0]=y_source[0];
	oldy[1]=y_source[1];
	oldr=r_source;

	delete subkist;

	for(i=*Nimages;i<oldNimages;i++){  // save some space
		EmptyKist(imageinfo[i].innerborder);
		EmptyKist(imageinfo[i].outerborder);		EmptyKist(imageinfo[i].imagekist);
	}
	oldNimages=*Nimages;

	// remove images without points
	for(j=0;j<*Nimages;++j){
		if(imageinfo[j].imagekist->Nunits() < 1){
			assert(imageinfo[j].area == 0);
			assert(*Nimages < NimageMax);
			ERROR_MESSAGE();
			for(k=j+1;k<*Nimages;++k) SwapImages(&imageinfo[k-1],&imageinfo[k]);
			//printf("image %i has no points\n",j);
			--*Nimages;
			--j;
		}
	}
	assert(*Nimages > 0);

	// calculate the centroid of the images assuming uniform surface brightness
	for(i=0;i<*Nimages;++i){
		tmp=0.0;
		imageinfo[i].centroid[0] = 0.0;
		imageinfo[i].centroid[1] = 0.0;
		MoveToTopKist(imageinfo[i].imagekist);
		do{
			tmp += pow(getCurrentKist(imageinfo[i].imagekist)->gridsize,2);
			imageinfo[i].centroid[0] += getCurrentKist(imageinfo[i].imagekist)->x[0]
			                                         *pow(getCurrentKist(imageinfo[i].imagekist)->gridsize,2);
			imageinfo[i].centroid[1] += getCurrentKist(imageinfo[i].imagekist)->x[1]
			                                         *pow(getCurrentKist(imageinfo[i].imagekist)->gridsize,2);

		}while(MoveDownKist(imageinfo[i].imagekist));
		if(imageinfo[i].imagekist->Nunits() > 0 ){
			imageinfo[i].centroid[0] /= tmp;
			imageinfo[i].centroid[1] /= tmp;
		}

		// redefine error so that it is based on the smallest grid cell on the border of the image
		if(imageinfo[i].outerborder->Nunits() > 0 ) imageinfo[i].area_error = imageinfo[i].gridrange[2]/imageinfo[i].area;
	}

    grid->ClearAllMarks();

	//freeKist(pointkist);


	return;
}

/** \ingroup ImageFinding
 *
 * \brief  Finds images given a source position and size.
 *
 * find_image_kist returns finite refined images in images[0...*Nimages-1].imagekist
 *  It starts with a large source and reduces down to the right size refining at each step.
 *  It should not miss any image larger than ~ munin*r_source linear size, but seems to do
 *  much better than that.  It does nothing with the surface brightnesses.
 *
 * The routine can follow three different strategies for refining each image controlled by edge_refinement.
 *
 * edge_refinement
 *   - 0 does not do edge refinement, Every pixel in every image is refined until the criterion is met.
 *  The image(s) are found again after each refinement which can make it slower.
 *   - 1 uses refine_edge().  After an initial refinement of all the pixels in the image(s) the code switches
 *	 to refining only the edges of the images.  The images are found after each refinement.
 *   - 2 uses refine_edge2() Same as 1, but the images are not found after each refinement.  This can make the
 *	 routine run much faster, but has the disadvantage that the number of images will not change during the final
 *	 stage of refinement.  This is the setting generally recommended.
 *
 *  kappa_off - This is used to turn off the calculation of surface density, shear, magnification and time delay.
 *  When finding a finite sized source these quantities are generally not required and slow down the routine.
 *
 *
 */

void find_images_microlens(
		BaseAnaLens *lens,          /// contains the lens/es and source/sources
		double *y_source        /// position of source center
		,double r_source        /// radius of source
		,GridHndl grid          /// grid provided to routine
		,int *Nimages           /// number of images found
		,ImageInfo *imageinfo   /// information on each image
		,const int NimageMax    /// maximum number of images allowed
		,unsigned long *Nimagepoints  /// number of points in final images
		,double initial_size    /// Initial size of source for telescoping, 0 to start from the initial grid size.
        ,double mu_min
		,bool splitimages       /// TRUE each image is refined to target accuracy, otherwise all images are treated as one
		,short edge_refinement  /// see comment
		,bool verbose           /// verbose
		,bool kappa_off         /// turns off calculation of surface density, shear, magnification and time delay
		){


	if(  grid->s_tree->top->boundary_p1[0] > (y_source[0] + r_source)
	  || grid->s_tree->top->boundary_p2[0] < (y_source[0] - r_source)
	  || grid->s_tree->top->boundary_p1[1] > (y_source[1] + r_source)
	  || grid->s_tree->top->boundary_p2[1] < (y_source[1] - r_source)
	){
		// source is not within initialized grid
		*Nimages = 0;
		std::cout << "source not within initialized grid" << std::endl;
		ERROR_MESSAGE();
		return;
	}

	int Nsizes;
	double rtemp,tmp;
	static double oldy[2],oldr=0;
	short moved,flag;
	int i,j,k;
	//Point *i_points,*s_points,*point;
	time_t to,t1,t2,t3,now;
	//Kist<Point> * tmp_border_kist;
	bool image_overlap;
	static int oldNimages=0;
	static unsigned long Npoints_old = 0;

	bool time_on = false;

	double **xstars = lens->stars_xp;
	int Ngrid_block = grid->getNgrid_block();

	if(r_source==0.0){ERROR_MESSAGE(); printf("ERROR: find_images, point source must have a resolution target\n"); exit(1);}

	if(verbose) printf("initialgridsize=%e\n",initialgridsize);
	if(initial_size==0) initial_size=initialgridsize;

	if(oldr==0){ oldr=r_source; Npoints_old = grid->i_tree->pointlist->Npoints;}
	if((Npoints_old <= grid->i_tree->pointlist->Npoints )* // if grid has not been refreshed
			(oldy[0]==y_source[0])*(oldy[1]==y_source[1])* // and source not moved
			(oldr > r_source)  // and source size has gotten smaller
	){
		Nsizes=(int)( log(oldr/r_source/mumin)/log(Ngrid_block) ); // round up
	    rtemp = r_source*pow(1.0*Ngrid_block,Nsizes);
	}else{
		Nsizes=(int)(log(initial_size/fabs(r_source*mumin))/log(Ngrid_block) ) + 1 ; // round up
	    rtemp = r_source*pow(1.0*Ngrid_block,Nsizes);
	}

	Npoints_old = grid->i_tree->pointlist->Npoints;


	// starting with a larger source size make sure all the grid sizes are small enough to find it
	Kist<Point> * subkist = new Kist<Point>;//,pointkist = NewKist();

	if(verbose) printf("entering find_image\n");
	time(&to);

    /**** TODO test line **********************
	bool map_on = true;
    PixelMap map(2000,(grid->i_tree->top->boundary_p2[0]-grid->i_tree->top->boundary_p1[0])/2,grid->i_tree->top->center);
    char chrstr[100];
    std::string output = "image_proportional";
    int Nmaps=0;
    // ************************************/

    if(verbose) printf("Ntemp=%li\n",Nsizes);

    grid->ClearAllMarks();

    //////////////////////////////////////////
    // telescope source size down to target
    //////////////////////////////////////////
	//if(!( (oldy[0]==y_source[0])*(oldy[1]==y_source[1])*(oldr < r_source) )){

    //unsigned long minN = 0;
    //float telescope_factor = 1.0/Ngrid_block;
    //float telescope_factor = 0.5;
    float telescope_factor = 0.66;
    double time_in_refine = 0,time_in_find = 0;
    ImageInfo *critcurve = new ImageInfo[NimageMax];
    int Ncrits;
    bool dummybool;

    for(i=0
				//for(rtemp = fabs(r_source/mumin)*pow(Ngrid_block,Nsizes),Nold=0
		//		;rtemp >= 0.99*Ngrid_block*fabs(r_source)
    		;rtemp >= r_source
       		;rtemp *= telescope_factor,++i ){

    	time(&t3);

    	/************* method that seporates images ****************
    	// mark image points in tree
    	grid->s_tree->PointsWithinKist(y_source,rtemp,subkist,1);
    	moved = image_finder_kist(lens,y_source,rtemp,grid
    			,Nimages,imageinfo,NimageMax,Nimagepoints,0,0);
    	// unmark image points in tree
    	grid->s_tree->PointsWithinKist(y_source,rtemp,subkist,-1);
    	// ***********************************************************/

    	//imageinfo->imagekist->Print();

/*    	minN = imageinfo[0].getNimagePoints();
    	for(int k=1; k < *Nimages; ++k){
    		minN = minN < imageinfo[k].getNimagePoints() ? minN : imageinfo[k].getNimagePoints();
    		//imageinfo[k].area/pi/r_source/r_source;
    	}

    	if(minN < 5  && i > 0){
    		// the size of the images jumped too quickly when the size changed so it needs to
    		//  be done again at a higher resolution

    		rtemp /= telescope_factor;  //  turn back
    		--i;
    		++nstep;
    		telescope_factor = (1+telescope_factor)/2;

    		if(verbose){
    			printf("      repeating telescoping source scale \n"
    					,difftime(t2,t1),*Nimagepoints,*Nimages);
    			for(int k=0; k < *Nimages; ++k){
    				std::cout << "   " << imageinfo[k].area << " " << imageinfo[k].area_error << " " << imageinfo[k].area/pi/rtemp/rtemp
    						<< " " << imageinfo[k].getNimagePoints() << std::endl;
    			}
    		}

    		continue;
    	}else{
    		telescope_factor = 1.0/Ngrid_block;
    		nstep = 1;
    	}
*/

    	time(&t1);
    	time(&t2);
    	if(verbose)
				printf("\n   new source size = %e    Nimages = %i  telescoping rsource = %e\n",rtemp,*Nimages,r_source);
/*
   		for(int k=0; k < *Nimages; ++k){
            if( 5.0e-3 > imageinfo[k].area/pi/r_source/r_source ) imageinfo[k].ShouldNotRefine = true;
            else imageinfo[k].ShouldNotRefine = false;
        }
*/

 		moved = image_finder_kist(lens,y_source,rtemp,grid
    				,Nimages,imageinfo,NimageMax,Nimagepoints,-1,0);

    	j=0;
        time_in_refine = time_in_find = 0;
        time(&now);
    	//while(refine_grid_kist(lens,grid,imageinfo,*Nimages,telescope_res,3,kappa_off)){
       	//while(refine_grid_kist(lens,grid,imageinfo,*Nimages,mu_min,0,kappa_off)){
       	 while(refine_grid_kist(lens,grid,imageinfo,*Nimages,mu_min,0,false)){
       	    	//while( refine_grid_kist(lens,grid,imageinfo,*Nimages,mu_min*telescope_factor*telescope_factor,3,kappa_off) ){
    	//do{

    		time(&t1);
            time_in_refine += difftime(t1, now);

    		if(verbose) std::cout << "    refined images" << std::endl;

            /*for(int i=0;i<lens->stars_N;++i){
                 grid->zoom(lens,xstars[i],rtemp*0.1/Ngrid_block,kappa_off);//,tmp);
            }*/

    		// refine critical curves
    		//find_crit(lens,grid,critcurve,NimageMax,&Ncrits,rtemp*0.01,&dummybool,false,false,verbose);
            refine_crit_in_image(lens,grid,r_source,y_source,rtemp*0.01);

//    		std::cout << "    Ncrits = " << Ncrits << " with " << critcurve->imagekist->Nunits() << " points." << std::endl;

    		//************* method that does not separate images ****************
    		moved = image_finder_kist(lens,y_source,rtemp,grid
    				,Nimages,imageinfo,NimageMax,Nimagepoints,-1,0);
    		//imageinfo->imagekist->Print();
    		// *************  ****************/


       		/************* method that separates images ****************
      		for(int k=0; k < *Nimages; ++k) imageinfo[k].ShouldNotRefine = false;
      		// mark image points in tree
    		//grid->s_tree->PointsWithinKist(y_source,rtemp,subkist,1);
    		moved = image_finder_kist(lens,y_source,rtemp,grid
    				,Nimages,imageinfo,NimageMax,Nimagepoints,-1,0);
            divide_images_kist(grid->i_tree,imageinfo,Nimages,NimageMax);
            for(int k=0; k < *Nimages; ++k) imageinfo[k].outerborder->Empty();
    		// unmark image points in tree
    		grid->s_tree->PointsWithinKist(y_source,rtemp,subkist,-1);

    		for(int k=0; k < *Nimages; ++k) if( 1.0e-2 > imageinfo[k].area/pi/r_source/r_source ) imageinfo[k].ShouldNotRefine = true;

    		/***********************************************************/
            /**** TODO test line **********************
    		if(map_on){
    			map.AddImages(imageinfo, *Nimages, true);
    			snprintf(chrstr,100,"%i",Nmaps++);
    			map.printFITS(output+chrstr+".fits");
    			map.Clean();
    		}
            // ************************************/

            time(&now);
            time_in_find += difftime(now , t1);

            if(time_on) printf("    time in finding images %f sec\n",difftime(now,t1));
            assert(*Nimages > 0);

    		if(verbose){
    			printf("      refound images after refinement\n        Nimagepoints=%li  Nimages = %i\n"
						,*Nimagepoints,*Nimages);
    			for(int k=0; k < *Nimages; ++k){
    				std::cout << "   " << imageinfo[k].area << " " << imageinfo[k].area_error << " " << imageinfo[k].area/pi/rtemp/rtemp
    						<< " " << imageinfo[k].area/pi/r_source/r_source << " " << imageinfo[k].getNimagePoints() << std::endl;
    			}
    		}

    		++j;
    	//}while(refine_grid_kist(lens,grid,imageinfo,*Nimages,rtemp*mumin/Ngrid_block,2,kappa_off,NULL));
    	}

  		//for(int k=0; k < *Nimages; ++k) imageinfo[k].ShouldNotRefine = false;
        
    	time(&now);
    	if(time_on) printf("    time for one source size %f sec\n",difftime(now,t3));
    	if(time_on) printf("        time for refinement %f sec and image finding %f sec\n",time_in_refine,time_in_find);
    }

	time(&now);
	if(time_on) printf(" time for source size reduction %f sec\n",difftime(now,to));
	time(&to);

	// refine critical curves
	//find_crit(lens,grid,critcurve,NimageMax,&Ncrits,r_source*0.01,&dummybool,false,false,verbose);
    refine_crit_in_image(lens,grid,r_source,y_source,r_source*0.01);

	time(&now);
	if(time_on) printf(" time for refine critical curves %f sec\n",difftime(now,to));
	time(&to);

	/**** Refine grid around each star so that no images are missed
	if(lens->AreStarsImaplated()){
		if(verbose) std::cout << "  zooming in around " << lens->stars_N << " stars" << std::endl;
		//Branch *tmp = grid->i_tree->current;
		for(i=0;i<lens->stars_N;++i){
			//if(verbose) std::cout << "  xstars[i] " << i << " " << xstars[i][0] << " " << xstars[i][1] << std::endl;
			grid->zoom(lens,xstars[i],r_source*0.1,kappa_off);//,tmp);
		}
		if(verbose) std::cout << "  out of zoom" << std::endl;
	}
	time(&now);
	if(time_on) printf(" time for zooming on stars %f sec\n",difftime(now,to));
	time(&to);
	//**********************************************************************/

	////////////////////////////////////////////////////////////////////////////////*/
	//////////////////////////////////////////////////////////////////////////////////
	// target source size has been reached, do full image decomposition and refinement
	/////////////////////////////////////////////////////////////////////////////////

    grid->ClearAllMarks();

	i=0;

	//if(splitimages) flag = 1; else flag = 0;
	flag = 1; // Changed so that each image always has at least 100 points.

	time(&now);

	assert(*Nimages > 0);
	// do an initial uniform refinement to make sure there are enough point in
	//  the images
	i=0;
    double area_tot = 0;
    int count =0;
	do{
		time(&t3);
		if(verbose)
			printf("     time in image refinement %f min\n           points in grid=%li\n"
					,fabs(difftime(t3,now)/60.),grid->i_tree->pointlist->Npoints);

		// mark image points in tree
		grid->s_tree->PointsWithinKist(y_source,r_source,subkist,1);

        for(int k=0;k<*Nimages;++k) imageinfo[k].ShouldNotRefine = false;
		//moved=image_finder_kist(lens,y_source,fabs(r_source),grid
		//		,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);
		moved=image_finder_kist(lens,y_source,fabs(r_source),grid
				,Nimages,imageinfo,NimageMax,Nimagepoints,0,0);

        area_tot = 0;
        count = 0;
        for(int k=0;k < *Nimages;++k){
            //std::cout << "area_tot = " << area_tot << "  area " << imageinfo[k].area << std::endl;
            area_tot += imageinfo[k].area;
        }
        for(int k=0;k < *Nimages;++k){
            if(imageinfo[k].area < area_tot*1.0e-3){
                ++count;
                imageinfo[k].ShouldNotRefine = true;
            }else{
                imageinfo[k].ShouldNotRefine = false;
            }
        }
		//if(*Nimages < 1) printf("  Nimages=%i i=%i\n",*Nimages,i);

		time(&now);
		if(verbose){
			printf("\n    i=%i\n     time in finding images %f min\n          Nimages=%i   Nimagepoints=%li\n"
					,i,difftime(now,t3)/60.,*Nimages,*Nimagepoints);
			printf("     image   # of points    error in area\n");
			for(j=0;j<*Nimages;++j) printf("       %i        %li         %e\n",j,imageinfo[j].imagekist->Nunits(),imageinfo[j].area_error);
		}
		if(i > 20 && *Nimagepoints == 100){
			// case where no image is found at any size
			*Nimages = 0;
			*Nimagepoints = 0;
			return ;
		}
		++i;
	}while( refine_grid_kist(lens,grid,imageinfo,*Nimages,0.1,3,kappa_off,NULL)
			|| moved );
    for(int k=0;k<*Nimages;++k) imageinfo[k].ShouldNotRefine = false;

	time(&now);
	if(time_on) printf(" time for uniform refinement %f sec\n",difftime(now,to));
	time(&to);

	assert(*Nimages > 0);

	// find points that are truly in the image and not just neighbors
	moved=image_finder_kist(lens,y_source,fabs(r_source),grid
			,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);
	assert(*Nimages > 0);


	// remove images with no points in them
	for(j=0;j<*Nimages;++j){
		if(imageinfo[j].imagekist->Nunits() < 1){
			ERROR_MESSAGE();
			for(k=j+1;k<*Nimages;++k) SwapImages(&imageinfo[k-1],&imageinfo[k]);
			//printf("image %i has no points\n",j);
			--*Nimages;
			--j;
		}
	}
	assert(*Nimages > 0);

	/////////////////////////////////////////////
	// second stage of refinement -
	// depends of choice of edge_refinement
	/////////////////////////////////////////////
	time(&now);

	if(splitimages) flag = 0; else flag = 2;

	k=i;
	if(edge_refinement==0){   // uniform refinement over image
		do{
			// mark image points in tree
			grid->s_tree->PointsWithinKist(y_source,r_source,subkist,1);

			moved=image_finder_kist(lens,y_source,fabs(r_source),grid
					,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);
			++i;
		}while( refine_grid_kist(lens,grid,imageinfo,*Nimages,FracResTarget,0,kappa_off,NULL)
				|| moved );

	}else if(edge_refinement==1){    // edge refinement with image finding at each step
		do{
			// mark image points in tree
			grid->s_tree->PointsWithinKist(y_source,r_source,subkist,1);

			moved=image_finder_kist(lens,y_source,fabs(r_source),grid
					,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);

			//for(i = 0; i < *Nimages; ++i) PrintImageInfo(&(imageinfo[i]));
			//printf("\n");

			++i;
		}while( refine_edges(lens,grid,imageinfo,*Nimages,FracResTarget,flag,kappa_off)
				|| moved );

	}else if(edge_refinement==2){  // edge refinement with no image finding at each step
		++i;
        area_tot = 0;
        count=0;
        for(int k=0;k<*Nimages;++k) area_tot += imageinfo[k].area;
        for(int k=0;k<*Nimages;++k){            
            if(imageinfo[k].area < area_tot*1.0e-3){++count; imageinfo[k].ShouldNotRefine = true;}
            else imageinfo[k].ShouldNotRefine = false;
        }
		while(refine_edges2(lens,y_source,r_source,grid
				,imageinfo,&image_overlap,*Nimages,FracResTarget,flag,kappa_off)){
			// if an overlap is detected find the images again

			if(image_overlap) moved=image_finder_kist(lens,y_source,fabs(r_source),grid
					,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);
			++i;
		}
        for(int k=0;k<*Nimages;++k) imageinfo[k].ShouldNotRefine = false;

	}
	// unmark image points so new source can be used
	grid->s_tree->PointsWithinKist(y_source,r_source,subkist,-1);

	if(verbose) printf("finished edge refinement i=%i\n",i);

	assert(*Nimages > 0);

	time(&t3);
	if(time_on) printf("     time in edge refinement %f sec\n",difftime(t3,now));

	// if point source take only closest image point
	if(r_source <= 0){
		ERROR_MESSAGE();
		exit(1);
	}

	time(&now);
	if(time_on) printf("time in find_images %f min\n",difftime(now,to)/60.);

	oldy[0]=y_source[0];
	oldy[1]=y_source[1];
	oldr=r_source;

	delete subkist;

	for(i=*Nimages;i<oldNimages;i++){  // save some space
		EmptyKist(imageinfo[i].innerborder);
		EmptyKist(imageinfo[i].outerborder);		EmptyKist(imageinfo[i].imagekist);
	}
	oldNimages=*Nimages;

	// remove images without points
	for(j=0;j<*Nimages;++j){
		if(imageinfo[j].imagekist->Nunits() < 1){
			assert(imageinfo[j].area == 0);
			assert(*Nimages < NimageMax);
			ERROR_MESSAGE();
			for(k=j+1;k<*Nimages;++k) SwapImages(&imageinfo[k-1],&imageinfo[k]);
			//printf("image %i has no points\n",j);
			--*Nimages;
			--j;
		}
	}
	assert(*Nimages > 0);

	// calculate the centroid of the images assuming uniform surface brightness
	for(i=0;i<*Nimages;++i){
		tmp=0.0;
		imageinfo[i].centroid[0] = 0.0;
		imageinfo[i].centroid[1] = 0.0;
		MoveToTopKist(imageinfo[i].imagekist);
		do{
			tmp += pow(getCurrentKist(imageinfo[i].imagekist)->gridsize,2);
			imageinfo[i].centroid[0] += getCurrentKist(imageinfo[i].imagekist)->x[0]
			                                         *pow(getCurrentKist(imageinfo[i].imagekist)->gridsize,2);
			imageinfo[i].centroid[1] += getCurrentKist(imageinfo[i].imagekist)->x[1]
			                                         *pow(getCurrentKist(imageinfo[i].imagekist)->gridsize,2);

		}while(MoveDownKist(imageinfo[i].imagekist));
		if(imageinfo[i].imagekist->Nunits() > 0 ){
			imageinfo[i].centroid[0] /= tmp;
			imageinfo[i].centroid[1] /= tmp;
		}

		// redefine error so that it is based on the smallest grid cell on the border of the image
		if(imageinfo[i].outerborder->Nunits() > 0 ) imageinfo[i].area_error = imageinfo[i].gridrange[2]/imageinfo[i].area;
	}

    grid->ClearAllMarks();

	//freeKist(pointkist);


	return;
}


/** \ingroup ImageFindingL2
 *
 * \brief finds images for a given source position and size
 *
 * image points are put into imageinfo[].imagekist
 *      imageinfo[].points and imageinfo[].Npoints are not changed
 * imagemarker array must be large enough to include all images *
 * splitparities=  0 don't split attached negative and positive parity images
 *              =  1 do split parities  NOTE: this is now obsolete
 *              = -1 doesn't slit into images at all
 *	                 , also does not find borders or change in_image markers
 * true_images = 1 gives just the points that are in the image
 *             = 0 if there are not enough points in images this will include close points to be refined
 * side-effects :  Will make in_image = true for all image points if splitparities == 0
 */
short image_finder_kist(LensHndl lens, double *y_source,double r_source,GridHndl grid
		,int *Nimages,ImageInfo *imageinfo,const int NimageMax,unsigned long *Nimagepoints
		,short splitparities,short true_images){

  unsigned long i,Nsource_points=0;
  static long count=0,Nold_images;
  static double oldy[2];
  short moved;
  double r;

  //if(count==0) oldy[0]=oldy[1]=0;
  TreeHndl i_tree = grid->i_tree,s_tree = grid->s_tree;

  if(splitparities==1){ ERROR_MESSAGE(); printf("ERROR: image_finger, option splitparaties==1 is obsolete\n"); exit(1); }
  if(r_source <= 0.0){
	  ERROR_MESSAGE();
	  printf("ERROR: cannot do point source right now \n");
	  exit(1);
  }
//  printf("in image_finder\n");
  if( (oldy[0] != y_source[0]) ||  (oldy[1] != y_source[1]) ){
    oldy[0] = y_source[0];
    oldy[1] = y_source[1];
    moved=1;
  }else moved=0;

  ++count;

  grid->ClearAllMarks();
  assert(imageinfo->imagekist);

  // if source has moved make sure low res grids are included in first source list
  if( moved && (initialgridsize > mumin*r_source) && !true_images){
	  // if new source position use larger image to make sure new images are found on the coarser grid
	  s_tree->PointsWithinKist(y_source,initialgridsize/mumin,imageinfo->imagekist,0);

	  MoveToTopKist(imageinfo->imagekist);
	  for(i=0;i<imageinfo->imagekist->Nunits();++i){
		  r=sqrt( pow(getCurrentKist(imageinfo->imagekist)->x[0]-y_source[0],2)
				  + pow(getCurrentKist(imageinfo->imagekist)->x[1]-y_source[1],2) );
		  if(r > r_source && getCurrentKist(imageinfo->imagekist)->image->gridsize < r_source*mumin){
			  TakeOutCurrentKist(imageinfo->imagekist);
			  --i;
		  }else{
			  MoveDownKist(imageinfo->imagekist);
		  }
	  }
  }else{
	  // if source hasn't moved just take points within image
	  s_tree->PointsWithinKist(y_source,r_source,imageinfo->imagekist,0);

	  //	do something about non-circular sources
	  // the points outside of non-circular are removed from sourcelist
	  // in_source(y_source,sourcelist);

	  if(imageinfo->imagekist->Nunits() < 1  && true_images ){  // no points in the source
		  *Nimages=0;
		  Nimagepoints=0;
		  return moved;
	  }
  }

  Nsource_points = imageinfo->imagekist->Nunits();

  // if there are not enough points in source find nearest ones
  if(!true_images && imageinfo->imagekist->Nunits() < NpointsRequired){
	  if(imageinfo->imagekist->Nunits() == 0){
		  s_tree->NearestNeighborKist(y_source,NpointsRequired,imageinfo->imagekist);
	  }else{  // add nearest points to already found image points

		  bool redundant = false;

		  s_tree->NearestNeighborKist(y_source,NpointsRequired,imageinfo[1].imagekist);

		  i=0;
		  while(imageinfo[1].imagekist->Nunits() > 0){
			  imageinfo->imagekist->MoveToTop();
			  imageinfo->imagekist->JumpDown(i);
			  redundant = false;
			  do{
				  if(imageinfo->imagekist->getCurrent() == imageinfo[1].imagekist->getCurrent()){
					  imageinfo[1].imagekist->TakeOutCurrent();
					  redundant = true;
					  break;
				  }
			  }while(imageinfo->imagekist->Down());
			  if( !redundant ){
				  imageinfo->imagekist->MoveToTop();
				  imageinfo->imagekist->InsertBeforeCurrent(imageinfo[1].imagekist->TakeOutCurrent());
				  ++i;
			  }

		  }
	  }
  }

  // mark all image points
  if(imageinfo->imagekist->Nunits() > 0 && splitparities == 0){
	  MoveToTopKist(imageinfo->imagekist);
	  do{
		  getCurrentKist(imageinfo->imagekist)->in_image = TRUE;
		  getCurrentKist(imageinfo->imagekist)->image->in_image = TRUE;

	  }while(MoveDownKist(imageinfo->imagekist));
  }

  *Nimagepoints = imageinfo->imagekist->Nunits();

  // transform from source plane to image points
  TranformPlanesKist(imageinfo->imagekist);

  // At this point all the image points are in imageinfo->imagekist and not divided up into separate images

  // split images
  if( splitparities == 0 && imageinfo->imagekist->Nunits() < i_tree->pointlist->Npoints){
	  divide_images_kist(i_tree,imageinfo,Nimages,NimageMax);
  }else{
	  *Nimages = 1;

	  imageinfo->area = 0.0;
	  MoveToTopKist(imageinfo->imagekist);
	  do{
		  imageinfo->area += pow(getCurrentKist(imageinfo->imagekist)->gridsize,2);
	  }while(MoveDownKist(imageinfo->imagekist));
  }

  // don't copy information into array
  //for(i=0;i<*Nimages;++i) imageinfo[i].Npoints = 0;  // to make sure points array is not read beyond length

  // find borders
  if( splitparities == 0 ) for(i=0;i<*Nimages;++i) findborders4(i_tree,&imageinfo[i]);

  for(i=0;i<*Nimages;++i){

	  if(splitparities == -1){
		  // avoid finding image borders, but need to set grid range for each image
		  EmptyKist(imageinfo[i].innerborder);
		  EmptyKist(imageinfo[i].outerborder);

		  imageinfo[i].gridrange[2] = 1.0e99; // minimum grid size in image
		  imageinfo[i].gridrange[0] = 0.0;    // maximum grid size in outerborder
		  imageinfo[i].gridrange[1] = 0.0;    // maximum grid size in image

		  //for(j=0;j<imageinfo[i].imagekist->Nunits();++j){

		  MoveToTopKist(imageinfo[i].imagekist);
		  do{
		    if(imageinfo[i].gridrange[1] < getCurrentKist(imageinfo[i].imagekist)->gridsize)
		      imageinfo[i].gridrange[1] = getCurrentKist(imageinfo[i].imagekist)->gridsize;
		    if(imageinfo[i].gridrange[2] > getCurrentKist(imageinfo[i].imagekist)->gridsize)
		      imageinfo[i].gridrange[2] = getCurrentKist(imageinfo[i].imagekist)->gridsize;
		  }while(MoveDownKist(imageinfo[i].imagekist));

		  imageinfo[i].gridrange[0] = imageinfo[i].gridrange[1];
	  }

	  // find area of images
	  //findarea(&imageinfo[i]);  // ****this is now done in divide_images

	  assert(imageinfo[i].area >= 0.0);
	  if(Nsource_points < NpointsRequired || moved) imageinfo[i].area_error=1.0;
	  else imageinfo[i].area_error = pow(imageinfo[i].gridrange[1],2)/imageinfo[i].area;
 }

  // Empty the border lists of old images to save mem
  for(i=*Nimages;i<Nold_images;++i){
	  EmptyKist(imageinfo[i].innerborder);
	  EmptyKist(imageinfo[i].outerborder);
  }

	//for(i=0;i<*Nimages;++i) printf("  image %i  Npoints = %li Ninner = %li Noutter = %li  area = %e\n",i
  //,imageinfo[i].Npoints,imageinfo[i].innerborder->Nunits,imageinfo[i].outerborder->Nunits,imageinfo[i].area);

   return moved;
}


/** \ingroup ImageFindingL2
 *
 *\brief Refines every point in the given image and its outer border that satisfies the refinement criterion.
 *
 *
 * Same as refine_grid with additions
 *     - uses imageinfo->imagekist instead of imageinfo->points[]
 *     - can export rayshooting
 *
 *     - Does not affect imageinfo->points in any way
 *     - Allows for refinements to be done on the borders of the initial grid region
 *          The initial grid region is never expanded.
 *
 * The borders of the image must be found previously.
 *
 * criterion = 0 stops refining when error in total area reaches res_target
 * 	         = 1 stops refining when each image reaches error limit or is smaller than res_target, (imageinfo[i].area_error > res_target)*(imageinfo[i].area > 1.0e-2*res_target*total_area)
 *           = 2 stops refining when grid resolution is smaller than res_target in all images
 *           = 3 stops refining when each image reaches error limit or imageinfo[i].ShouldNotRefine == true
 *
 * Returns the number of points that were added to the grids.
 * 
 */
int refine_grid_kist(
	LensHndl lens            /// the lens model
	,GridHndl grid           /// the grid
	,ImageInfo *imageinfo    /// images
	,unsigned long Nimages   /// Number of images to refine
	,double res_target       /// meaning depends on criterion, see general notes
	,short criterion         /// see general notes
	,bool kappa_off          /// true = no kappa, gamma and dt are calculated
	,Kist<Point> * newpointskist  /// returns a Kist of the points that were added to the grid on this pass, if == NULL will not be added
	,bool batch              /// True, passes all points to rayshooter at once, False shoots rays each cell at a time and new points are in memory blocks of 8 or smaller
	){

  if(newpointskist)
	newpointskist->Empty();
	//printf("entering refine_grid\n");

  if(Nimages < 1) return 0;

  TreeHndl i_tree = grid->i_tree,s_tree = grid->s_tree;
  int Ngrid_block = grid->getNgrid_block();

  int i,j,k,number_of_refined,count;
  double rmax,total_area;
  Point *point;
  short pass=0;
  long Ncells=0,Ncells_o=0;
  Point *i_points;

  std::vector<Point *> points_to_refine;

  total_area=0;
  for(i=0;i<Nimages;++i) total_area += imageinfo[i].area;

  number_of_refined = Ncells = 0;
  for(i=0;i<Nimages;++i){
    count=0;

    if(criterion == 0) pass = imageinfo[i].area*imageinfo[i].area_error/total_area > res_target;
    if(criterion == 1) pass = (imageinfo[i].area_error > res_target)*(imageinfo[i].area > 1.0e-5*total_area);
    if(criterion == 2) pass = imageinfo[i].gridrange[1] > res_target;
    if(criterion == 3) pass = (imageinfo[i].area_error > res_target)*(!(imageinfo[i].ShouldNotRefine));

    // make sure no border point has a lower res than any image point

    if( pass || imageinfo[i].gridrange[0]>1.01*imageinfo[i].gridrange[1]){

    	/* largest gridsize in image */
    	rmax=MAX(imageinfo[i].gridrange[1],imageinfo[i].gridrange[0]);

      // count number of grid cells to be refined
      // cells in image
    	MoveToTopKist(imageinfo[i].imagekist);
    	do{
    		if( getCurrentKist(imageinfo[i].imagekist)->gridsize > 1.01*rmax/Ngrid_block) ++Ncells;
    	}while( MoveDownKist(imageinfo[i].imagekist) );

    	//printf("   initial image point refinements count = %li\n",Ncells);
       	// cells on border
    	if(imageinfo[i].outerborder->Nunits() > 0){
    		MoveToTopKist(imageinfo[i].outerborder);
    		do{
    			if( getCurrentKist(imageinfo[i].outerborder)->gridsize > 1.01*rmax/Ngrid_block){
    				// border point is marked to prevent refining more than once
    				//   it will be unmarked by the end of refine grid
    				getCurrentKist(imageinfo[i].outerborder)->in_image = TRUE;
    				++Ncells;
    			}
    		}while( MoveDownKist(imageinfo[i].outerborder) );
    		//printf("   initial outer border refinements count = %li\n",Ncells);
    	}
    }
  }

  Ncells_o = Ncells;
  Ncells = 0;

  for(i=0,Ncells=0;i<Nimages;++i){
	  count=0;

	  if(criterion == 0) pass = imageinfo[i].area*imageinfo[i].area_error/total_area > res_target;
	  if(criterion == 1) pass = (imageinfo[i].area_error > res_target)*(imageinfo[i].area > 1.0e-2*res_target*total_area);
	  if(criterion == 2) pass = imageinfo[i].gridrange[1] > res_target;
	  if(criterion == 3) pass = (imageinfo[i].area_error > res_target)*(!(imageinfo[i].ShouldNotRefine));

	  // make sure no border point has a lower res than any image point

	  if( pass || imageinfo[i].gridrange[0]>1.01*imageinfo[i].gridrange[1]){

		  rmax=MAX(imageinfo[i].gridrange[1],imageinfo[i].gridrange[0]);

		  // loop through points in ith image
		  MoveToTopKist(imageinfo[i].imagekist);
		  for(j=0 ; j<imageinfo[i].imagekist->Nunits() ; ++j,MoveDownKist(imageinfo[i].imagekist) ){

			  if( getCurrentKist(imageinfo[i].imagekist)->gridsize > 1.01*rmax/Ngrid_block){  /* only refine largest grid size in image*/

				  //assert(getCurrentKist(imageinfo[i].imagekist)->leaf->child1 == NULL);
				  //assert(getCurrentKist(imageinfo[i].imagekist)->leaf->child2 == NULL);
				  //assert(getCurrentKist(imageinfo[i].imagekist)->image->leaf->child1 == NULL);
				  //assert(getCurrentKist(imageinfo[i].imagekist)->image->leaf->child2 == NULL);

				  if(batch){
					  points_to_refine.push_back(getCurrentKist(imageinfo[i].imagekist));
				  }else{
					  i_points = grid->RefineLeaf(lens,getCurrentKist(imageinfo[i].imagekist),kappa_off);
					  if(newpointskist && i_points != NULL){
						  for(k=0; k < i_points->head ; ++k) newpointskist->InsertAfterCurrent(&i_points[k]);
						  ++count;
					  }
				  }

				  ++Ncells;
			  }
		  }

		  //printf("   actual image point refinements count = %li\n",Ncells);
		  // * loop through outer border of ith image *

		  MoveToTopKist(imageinfo[i].outerborder);
		  for(j=0;j<imageinfo[i].outerborder->Nunits();++j,MoveDownKist(imageinfo[i].outerborder)){
			  if( getCurrentKist(imageinfo[i].outerborder)->gridsize > 1.01*rmax/Ngrid_block){ // only refine largest grid size in image

				  point = getCurrentKist(imageinfo[i].outerborder);
				  //assert(point->gridsize > 0);

				  if(point->in_image){ // point has not been refined yet as border of another image

					  assert(point->leaf->child1 == NULL);
					  assert(point->leaf->child2 == NULL);
					  assert(point->image->leaf->child1 == NULL);
					  assert(point->image->leaf->child2 == NULL);

					  if(batch){
						  points_to_refine.push_back(point);
					  }else{
						  i_points = grid->RefineLeaf(lens,point,kappa_off);
						  if(newpointskist && i_points != NULL){
							  for(k=0;k < i_points->head; ++k) newpointskist->InsertAfterCurrent(&i_points[k]);
							  ++count;;
						  }
					  }

					  ++Ncells;
					  point->in_image = FALSE;  // unmark so that it wouldn't double refine
				  }
			  }
		  }
    }

    if(count > 0) ++number_of_refined;
  } // end of image loop

  if(batch){
	  //for(i=0;i<points_to_refine.size();++i) assert(points_to_refine[i]->leaf->child1 == NULL && points_to_refine[i]->leaf->child2 == NULL);
	  i_points = grid->RefineLeaves(lens,points_to_refine,kappa_off);
	  if(newpointskist && i_points != NULL) for(k=0;k < i_points->head; ++k) newpointskist->InsertAfterCurrent(&i_points[k]);
	  if(i_points == NULL) number_of_refined = 0;
	  else number_of_refined = i_points->head;

	  points_to_refine.clear();
  }

  return number_of_refined;
}


/** \ingroup ImageFindingL2
 *
 * \brief Finds inner and outer borders of an image using bordering box method.
 *
 *   uses the in_image markers
 *   uses imaginfo->imagekist instead of imageinfo->points
 *
 *   In the case of the entire grid being within the image, the innerborders
 *   is the border points of the grid and the outerborder contains no points.
 *
 *   Note:  Markers in_image must be preset to true for all image points and
 *   false for non-image points.
 */
void findborders4(TreeHndl i_tree,ImageInfo *imageinfo){
	int i;
	unsigned long m,j;
	bool addinner;
	bool allin = false;


	if(i_tree->pointlist->Npoints == imageinfo->imagekist->Nunits()) allin = true;    // all points on the grid are in the image

	//printf("beginning findborders2\n");
	//checkTree(i_tree);

	//point=(Point *)malloc(sizeof(Point));
	EmptyKist(imageinfo->innerborder);
	EmptyKist(imageinfo->outerborder);

	imageinfo->gridrange[2] = 1.0e99; // minimum grid size in image
	imageinfo->gridrange[0] = 0.0;   // maximum grid size in outerborder
	imageinfo->gridrange[1] = 0.0;   // maximum grid size in image

	if(imageinfo->imagekist->Nunits() < 1) return;

	Kist<Point> * neighborkist = new Kist<Point>;

	Kist<Point> * imagekist = imageinfo->imagekist; assert(imageinfo->imagekist);

	MoveToTopKist(imagekist);
	for(j=0;j<imagekist->Nunits();++j,MoveDownKist(imagekist)){

		if(imageinfo->gridrange[1] < getCurrentKist(imagekist)->gridsize)
			imageinfo->gridrange[1] = getCurrentKist(imagekist)->gridsize;
		if(imageinfo->gridrange[2] > getCurrentKist(imagekist)->gridsize)
			imageinfo->gridrange[2] = getCurrentKist(imagekist)->gridsize;

		addinner=false;

		i_tree->FindAllBoxNeighborsKist(getCurrentKist(imagekist),neighborkist);

		if( allin && neighborkist->Nunits() < 4){
			InsertAfterCurrentKist(imageinfo->innerborder,getCurrentKist(imagekist));
		}else{

			MoveToTopKist(neighborkist);
			for(i=0;i<neighborkist->Nunits();++i){

				if( getCurrentKist(neighborkist)->in_image != TRUE){  // point is a neighbor
					addinner=true;
					/*
					// check if point is already in list
					MoveToTopKist(imageinfo->outerborder);
					for(m=0;m<imageinfo->outerborder->Nunits();++m){
						if( getCurrentKist(imageinfo->outerborder) == getCurrentKist(neighborkist) ) break;
						MoveDownKist(imageinfo->outerborder);
					}
					
					if(m == imageinfo->outerborder->Nunits()){
						// add point to outerborder
						InsertAfterCurrentKist(imageinfo->outerborder,getCurrentKist(neighborkist));
						MoveDownKist(imageinfo->outerborder);
					}
					*/
					
					if(getCurrentKist(neighborkist)->in_image == FALSE){  // if point is not yet in outerborder
						// add point to outerborder
						getCurrentKist(neighborkist)->in_image = MAYBE;
						InsertAfterCurrentKist(imageinfo->outerborder,getCurrentKist(neighborkist));
						MoveDownKist(imageinfo->outerborder);
					}
				}
				MoveDownKist(neighborkist);
			}

			if(addinner){
				// add point to innerborderkist
				InsertAfterCurrentKist(imageinfo->innerborder,getCurrentKist(imagekist));
				MoveDownKist(imageinfo->innerborder);
			}
		}

	}

		// mark outer borders back to in_image=FALSE
		
	if(!allin  && imageinfo->outerborder->Nunits() > 0){
		MoveToTopKist(imageinfo->outerborder);
		do{
			getCurrentKist(imageinfo->outerborder)->in_image = FALSE;

			if(imageinfo->gridrange[0] < getCurrentKist(imageinfo->outerborder)->gridsize)
				imageinfo->gridrange[0] = getCurrentKist(imageinfo->outerborder)->gridsize;
		}while(MoveDownKist(imageinfo->outerborder));
	}

	delete neighborkist;

	return;
}

void SwapImages(ImageInfo *image1,ImageInfo *image2){
	unsigned long Npoints,i;
	double tmp;
	Kist<Point> * list;

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

	list = image1->imagekist;
	image1->imagekist = image2->imagekist;
	image2->imagekist = list;

	list = image1->innerborder;
	image1->innerborder = image2->innerborder;
	image2->innerborder = list;

	list = image1->outerborder;
	image1->outerborder = image2->outerborder;
	image2->outerborder = list;

	return ;
}

