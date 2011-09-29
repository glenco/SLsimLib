#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <nrutil.h>
#include <Tree.h>
#include <KistDriver.h>
#include <divide_images.h>

static const int NpointsRequired = 50;  // number of points required to be within an image
static const int Ngrid_block = 3;       // each cell is divided into Ngrid_block^2 subcells

//static const float mumin = 0.5;  // actually the sqrt of the minimum magnification
static const float mumin = 0.45;  // actually the sqrt of the minimum magnification
//static const float mumin = 0.5;
static const float FracResTarget = 4.0e-4;
//static const float FracResTarget = 1.0e-4;

Point *pointg;
double ysourceg[2],magsigng;
extern const double initialgridsize;

void find_images_kist(double *y_source,double r_source
		  ,TreeHndl s_tree,TreeHndl i_tree,int *Nimages
		  ,ImageInfo *imageinfo,const int NimageMax,unsigned long *Nimagepoints
		  ,double initial_size,Boolean splitimages,short edge_refinement
		  ,Boolean verbose,Boolean kappa_off){

	/*
	 * find_image returns finite refined images in images[].imagekist
	 *  it starts with a large source and reduces down to the right size refining at each step
	 *  should not miss any image larger than ~ munin*r_source linear size
	 *
	 *  if splitimage==TRUE each image is refined to target accuracy, otherwise all images are treated as one
	 *
	 *  edge_refinement = 0 does not do edge refinement
	 *                  = 1 uses refine_edge() which keeps the images fully up to date
	 *                  = 2 uses refine_edge2() which is faster but does not keep edges
	 *  kappa_off = True - turns off calculation of kappa and gamma
	 *              False - turns on calculation of kappa and gamma
	 *
	 */

	long Nsizes;
	double rtemp,tmp;
	static double oldy[2],oldr=0;
	short moved,flag;
	int i,j,k,m;
	//Point *i_points,*s_points,*point;
	time_t to,t1,t2,t3,now;
	//KistHndl tmp_border_kist;
	Boolean image_overlap;
	static int oldNimages=0;
	static unsigned long Npoints_old = 0;
	Point **dummy_pnt = NULL;
	//unsigned long Ntmp;
	//Point *point,*closestpoint;

	if(r_source==0.0){ERROR_MESSAGE(); printf("ERROR: find_images, point source must have a resolution target\n"); exit(1);}

	if(verbose) printf("initialgridsize=%e\n",initialgridsize);
	if(initial_size==0) initial_size=initialgridsize;

	if(oldr==0){ oldr=r_source; Npoints_old = i_tree->pointlist->Npoints;}
	if((Npoints_old <= i_tree->pointlist->Npoints )* // if grid has not been refreshed
			(oldy[0]==y_source[0])*(oldy[1]==y_source[1])* // and source not moved
			(oldr > r_source)  // and source size has gotten smaller
	){
		Nsizes=(int)( log(oldr/r_source)/log(Ngrid_block) ); // round up
	    rtemp = r_source*pow(Ngrid_block,Nsizes);
	}else{
		Nsizes=(int)(log(initial_size/sqrt(pi)/fabs(r_source*mumin))/log(Ngrid_block) ) + 1 ; // round up
	    rtemp = r_source*pow(Ngrid_block,Nsizes);
	}

	Npoints_old = i_tree->pointlist->Npoints;


	// starting with a larger source size make sure all the grid sizes are small enough to find it
	KistHndl subkist = NewKist();//,pointkist = NewKist();

	if(verbose) printf("entering find_image\n");
	time(&to);


    if(verbose) printf("Ntemp=%li\n",Nsizes);

    ClearAllMarks(i_tree);

    //////////////////////////////////////////
    // telescope source size down to target
    //////////////////////////////////////////
	//if(!( (oldy[0]==y_source[0])*(oldy[1]==y_source[1])*(oldr < r_source) )){

		for(
				//for(rtemp = fabs(r_source/mumin)*pow(Ngrid_block,Nsizes),Nold=0
		//		;rtemp >= 0.99*Ngrid_block*fabs(r_source)
				;rtemp >= r_source
				;rtemp /= Ngrid_block ){

			time(&t1);
			time(&t3);
			if(verbose)
			printf("\n   new source size = %e    telescoping rsource = %e\n",rtemp,r_source);

			do{
				time(&t1);
				if(verbose) printf("      time in refine grid %f sec\n",difftime(t1,t2));

				moved = image_finder_kist(y_source,rtemp,s_tree,i_tree
						  ,Nimages,imageinfo,NimageMax,Nimagepoints,-1,0);

//				moved=image_finder_kist(y_source,rtemp,s_tree,i_tree
//						,Nimages,imageinfo,NimageMax,Nimagepoints,0,0);

				//printf("%li\n",imageinfo->Npoints);
				//for(i=0;i<imageinfo->Npoints;++i) printf("%e %e\n",imageinfo->points[i].x[0]
				//                                                  ,imageinfo->points[i].x[1]);
				//printf("rtemp=%e\n",rtemp);
				//exit(0);
				//printf("     %li %li Nimages=%i \n",imageinfo->Npoints
				//		,i_tree->pointlist->Npoints,*Nimages);

				time(&t2);
				if(verbose)	printf("      time in image_finder %f sec\n        Nimagepoints=%li\n"
						,difftime(t2,t1),*Nimagepoints);

				//printf("  telescoping Nimages=%li\n",*Nimages);

			}while(refine_grid_kist(i_tree,s_tree,imageinfo,*Nimages,rtemp*mumin/Ngrid_block,2,kappa_off,True,dummy_pnt));

			// testing lines //////////////////////////////////////////
			/*
//			moved=image_finder_kist(y_source,rtemp,s_tree,i_tree,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);

			for(MoveToTopKist(pointkist),i=0;i<pointkist->Nunits;++i,MoveDownKist(pointkist)){
				point = getCurrentKist(pointkist);
				printf("     %i closest point dr/rtemp = %e\n",i,sqrt(pow(point->image->x[0]-y_source[0],2) + pow(point->image->x[1]-y_source[1],2))/rtemp);
			}
			EmptyKist(pointkist);


			printf("        Nimages=%i during telescoping\n",*Nimages);
			for(i=0;i<*Nimages;++i){

				imageinfo[i].centroid[0] = 0.0;
				imageinfo[i].centroid[1] = 0.0;
				tmp = 0.0;
				tmp2 = 1.0e99;
				MoveToTopKist(imageinfo[i].imagekist);
				do{

					point = getCurrentKist(imageinfo[i].imagekist);
					tmp += pow(getCurrentKist(imageinfo[i].imagekist)->gridsize,2);
					imageinfo[i].centroid[0] += pow(getCurrentKist(imageinfo[i].imagekist)->gridsize,2)*getCurrentKist(imageinfo[i].imagekist)->x[0];
					imageinfo[i].centroid[1] += pow(getCurrentKist(imageinfo[i].imagekist)->gridsize,2)*getCurrentKist(imageinfo[i].imagekist)->x[1];

					if(tmp2 > pow(point->image->x[0]-y_source[0],2) + pow(point->image->x[1]-y_source[1],2)){
						tmp2 = pow(point->image->x[0]-y_source[0],2) + pow(point->image->x[1]-y_source[1],2);
						closestpoint = point;
					}
				}while(MoveDownKist(imageinfo[i].imagekist));

				InsertAfterCurrentKist(pointkist,closestpoint);
				MoveDownKist(pointkist);

				imageinfo[i].centroid[0] /= tmp;
				imageinfo[i].centroid[1] /= tmp;
				printf("         %i telescoping %e %e   %e +/- %e  %li gridrang = %e %e %e\n",i,imageinfo[i].centroid[0],imageinfo[i].centroid[1]
				                          ,imageinfo[i].area,imageinfo[i].area_error
						,imageinfo[i].imagekist->Nunits,imageinfo[i].gridrange[0],imageinfo[i].gridrange[1],imageinfo[i].gridrange[2]);
			}

			///////////////////////////////////////////////////////////////*/

/*			// Erase marks for next round
			if(imageinfo->imagekist->Nunits > 0){
				MoveToTopKist(imageinfo->imagekist);
				do{
					getCurrentKist(imageinfo->imagekist)->in_image = False;
					getCurrentKist(imageinfo->imagekist)->image->in_image = False;
				}while(MoveDownKist(imageinfo->imagekist));
			}
*/
			time(&t1);
			if(verbose)	printf("      time in refine grid %f sec\n",difftime(t1,t2));

			time(&now);
			if(verbose) printf("    time for one source size %f sec\n",difftime(now,t3));
		}
//	}//else printf("No telescoping!\n");

	time(&now);
	if(verbose) printf(" time for source size reduction %f sec\n",difftime(now,to));
	time(&to);

	/********    tests to be removed ************ //////
	moved=image_finder_kist(y_source,r_source,s_tree,i_tree,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);

	printf("  Nimages=%i after telescoping\n",*Nimages);
	for(i=0,Ntmp=0;i<*Nimages;++i){

		Ntmp += imageinfo[i].imagekist->Nunits;
		imageinfo[i].centroid[0] = 0.0;
		imageinfo[i].centroid[1] = 0.0;
		tmp = 0.0;
		MoveToTopKist(imageinfo[i].imagekist);
		do{
			tmp += pow(getCurrentKist(imageinfo[i].imagekist)->gridsize,2);
			imageinfo[i].centroid[0] += pow(getCurrentKist(imageinfo[i].imagekist)->gridsize,2)*getCurrentKist(imageinfo[i].imagekist)->x[0];
			imageinfo[i].centroid[1] += pow(getCurrentKist(imageinfo[i].imagekist)->gridsize,2)*getCurrentKist(imageinfo[i].imagekist)->x[1];
		}while(MoveDownKist(imageinfo[i].imagekist));
		imageinfo[i].centroid[0] /= tmp;
		imageinfo[i].centroid[1] /= tmp;
		printf("    %e %e   %e\n",imageinfo[i].centroid[0],imageinfo[i].centroid[1],imageinfo[i].area);
	}

	tmp_border_kist = NewKist();
	PointsWithinKist(s_tree,y_source,r_source,tmp_border_kist,0);

	printf("     points within = %li in images = %li   Nimagepoints = %li\n",tmp_border_kist->Nunits,Ntmp,*Nimagepoints);
	freeKist(tmp_border_kist);
	//exit(0);



	////////////////////////////////////////////////////////////////////////////////*/
	//////////////////////////////////////////////////////////////////////////////////
	// target source size has been reached, do full image decomposition and refinement
	/////////////////////////////////////////////////////////////////////////////////

	ClearAllMarks(i_tree);

	i=0;

	if(splitimages) flag = 1; else flag = 0;

	flag = 1; // Changed so that each image always has at least 100 points.

	time(&now);

	// do an initial uniform refinement to make sure there are enough point in
	//  the images
	do{
		time(&t3);
		if(verbose)
			printf("     time in image refinement %f min\n           points in grid=%li\n"
					,fabs(difftime(t3,now)/60.),i_tree->pointlist->Npoints);

		// mark image points in tree
		PointsWithinKist(s_tree,y_source,r_source,subkist,1);

		moved=image_finder_kist(y_source,fabs(r_source),s_tree,i_tree
				,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);

		if(*Nimages < 1) printf("  Nimages=%i i=%i\n",*Nimages,i);

		time(&now);
		if(verbose){
			printf("\n    i=%i\n     time in finding images %f min\n          Nimages=%i   Nimagepoints=%li\n"
					,i,difftime(now,t3)/60.,*Nimages,*Nimagepoints);
			printf("     image   # of points    error in area\n");
			for(j=0;j<*Nimages;++j) printf("       %i        %li         %e\n",j,imageinfo[j].imagekist->Nunits,imageinfo[j].area_error);
		}
		++i;
	}while( refine_grid_kist(i_tree,s_tree,imageinfo,*Nimages,1.0e-1,flag,kappa_off,True,dummy_pnt)
			|| moved );

	// remove images with no points in them
	for(j=0;j<*Nimages;++j){
		if(imageinfo[j].imagekist->Nunits < 1){
			ERROR_MESSAGE();
			for(k=j+1;k<*Nimages;++k) SwapImages(&imageinfo[k-1],&imageinfo[k]);
			//printf("image %i has no points\n",j);
			--*Nimages;
			--j;
		}
	}
	/////////////////////////////////////////////
	// refine just image edges to high resolution
	/////////////////////////////////////////////
	time(&now);

	if(splitimages) flag = 0; else flag = 2;

	k=i;
	if(edge_refinement==0){
		do{
			// mark image points in tree
			PointsWithinKist(s_tree,y_source,r_source,subkist,1);

			moved=image_finder_kist(y_source,fabs(r_source),s_tree,i_tree
					,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);
			++i;
		}while( refine_grid_kist(i_tree,s_tree,imageinfo,*Nimages,FracResTarget,0,kappa_off,True,dummy_pnt)
				|| moved );

	}else if(edge_refinement==1){
		do{
			// mark image points in tree
			PointsWithinKist(s_tree,y_source,r_source,subkist,1);

			moved=image_finder_kist(y_source,fabs(r_source),s_tree,i_tree
					,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);

			//for(i = 0; i < *Nimages; ++i) PrintImageInfo(&(imageinfo[i]));
			//printf("\n");

			++i;
		}while( refine_edges(i_tree,s_tree,imageinfo,*Nimages,FracResTarget,flag,kappa_off)
				|| moved );

	}else if(edge_refinement==2){
		++i;
		while(refine_edges2(y_source,r_source,i_tree,s_tree
				,imageinfo,&image_overlap,*Nimages,FracResTarget,flag,kappa_off)){
			// if an overlap is detected find the images again

			if(image_overlap) moved=image_finder_kist(y_source,fabs(r_source),s_tree,i_tree
					,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);
			++i;
		}
	}
	// unmark image points so new source can be used
	PointsWithinKist(s_tree,y_source,r_source,subkist,-1);

	if(verbose) printf("finished edge refinement i=%i\n",i);


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

	freeKist(subkist);
	for(i=*Nimages;i<oldNimages;i++){  // save some space
		EmptyKist(imageinfo[i].innerborder);
		EmptyKist(imageinfo[i].outerborder);		EmptyKist(imageinfo[i].imagekist);
	}
	oldNimages=*Nimages;

	// remove images without points
	for(j=0;j<*Nimages;++j){
		if(imageinfo[j].imagekist->Nunits < 1){
			assert(imageinfo[j].area == 0);
			assert(*Nimages < NimageMax);
			ERROR_MESSAGE();
			for(k=j+1;k<*Nimages;++k) SwapImages(&imageinfo[k-1],&imageinfo[k]);
			//printf("image %i has no points\n",j);
			--*Nimages;
			--j;
		}
	}

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
		if(imageinfo[i].imagekist->Nunits > 0 ){
			imageinfo[i].centroid[0] /= tmp;
			imageinfo[i].centroid[1] /= tmp;
		}

		// redefine error so that it is based on the smallest grid cell on the border of the image
		if(imageinfo[i].outerborder->Nunits > 0 ) imageinfo[i].area_error = imageinfo[i].gridrange[2]/imageinfo[i].area;
	}

	ClearAllMarks(i_tree);

	//freeKist(pointkist);
	//for(i = 0; i < *Nimages; ++i) PrintImageInfo(&(imageinfo[i]));

	return;
}


short image_finder_kist(double *y_source,double r_source,TreeHndl s_tree,TreeHndl i_tree
		,int *Nimages,ImageInfo *imageinfo,const int NimageMax,unsigned long *Nimagepoints
		,short splitparities,short true_images){
  /* finds images for a given source position and size
   * image points are put into imageinfo[].imagekist
   *      imageinfo[].points and imageinfo[].Npoints are not changed
  * imagemarker array must be large enough to include all images *
  * splitparities=  0 don't split attached negative and positive parity images
  *              =  1 do split parities  NOTE: this is now obsolete
  *              = -1 doesn't slit into images at all
  *	                 , also does not find borders or change in_image markers
  * true_images = 1 gives just the points that are in the image
  *             = 0 if there are not enough points in images this will include close points to be refined
  */
  unsigned long i,Nsource_points=0;
  static long count=0,Nold_images;
  static double oldy[2];
  short moved;
  double r;

  //if(count==0) oldy[0]=oldy[1]=0;

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

  ClearAllMarks(i_tree);
  assert(imageinfo->imagekist);

  //ListHndl list = NewList();

  // if source has moved make sure low res grids are included in first source list
  if( moved && (initialgridsize > mumin*r_source) && !true_images){
	  // if new source position use larger image to make sure new images are found on the coarser grid
	  PointsWithinKist(s_tree,y_source,initialgridsize/mumin,imageinfo->imagekist,0);

	  MoveToTopKist(imageinfo->imagekist);
	  for(i=0;i<imageinfo->imagekist->Nunits;++i){
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
	  PointsWithinKist(s_tree,y_source,r_source,imageinfo->imagekist,0);

	  //	do something about non-circular sources
	  // the points outside of non-circular are removed from sourcelist
	  // in_source(y_source,sourcelist);

	  if(imageinfo->imagekist->Nunits < 1  && true_images ){  // no points in the source
		  //EmptyKist(imageinfo->imagekist);
		  imageinfo->Npoints=0;
		  *Nimages=0;
		  Nimagepoints=0;
		  return moved;
	  }
  }

  Nsource_points = imageinfo->imagekist->Nunits;

  // if there are not enough points in source find nearest ones
  if(!true_images && imageinfo->imagekist->Nunits < NpointsRequired)
	  NearestNeighborKist(s_tree,y_source,NpointsRequired,imageinfo->imagekist);


  // mark all image points
  if(imageinfo->imagekist->Nunits > 0){
	  MoveToTopKist(imageinfo->imagekist);
	  do{
		  getCurrentKist(imageinfo->imagekist)->in_image = True;
		  getCurrentKist(imageinfo->imagekist)->image->in_image = True;
	  }while(MoveDownKist(imageinfo->imagekist));
  }

  *Nimagepoints = imageinfo->imagekist->Nunits;

  // transform from source plane to image points
  TranformPlanesKist(imageinfo->imagekist);

  // At this point all the image points are in imageinfo->imagekist and not divided up into separate images

  // split images
  if( splitparities == 0 && imageinfo->imagekist->Nunits < i_tree->pointlist->Npoints){
	  divide_images_kist(i_tree,imageinfo,Nimages,NimageMax);
  }else{
	  *Nimages = 1;
	  imageinfo->Npoints = 0;
	  imageinfo->area = 0.0;
	  MoveToTopKist(imageinfo->imagekist);
	  do{
		  imageinfo->area += pow(getCurrentKist(imageinfo->imagekist)->gridsize,2);
	  }while(MoveDownKist(imageinfo->imagekist));
  }

  // don't copy information into array
  for(i=0;i<*Nimages;++i) imageinfo[i].Npoints = 0;  // to make sure points array is not read beyond length

  // find borders
  if( splitparities == 0 ) for(i=0;i<*Nimages;++i) findborders4(i_tree,&imageinfo[i]);

  assert(*Nimages < NimageMax-1);

  //#pragma omp parallel for firstprivate(i_tree)

  for(i=0;i<*Nimages;++i){

	  if(splitparities == -1){
		  // avoid finding image borders, but need to set grid range for each image
		  EmptyKist(imageinfo[i].innerborder);
		  EmptyKist(imageinfo[i].outerborder);

		  imageinfo[i].gridrange[2] = 1.0e99; // minimum grid size in image
		  imageinfo[i].gridrange[0] = 0.0;    // maximum grid size in outerborder
		  imageinfo[i].gridrange[1] = 0.0;    // maximum grid size in image

		  //for(j=0;j<imageinfo[i].imagekist->Nunits;++j){

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




int refine_grid_kist(TreeHndl i_tree,TreeHndl s_tree,ImageInfo *imageinfo,unsigned long Nimages,double res_target
		,short criterion,Boolean kappa_off,Boolean shootrays,Point **point_pnt){

	/*
	 * Same as refine_grid with additions
	 *     - uses imageinfo->imagekist instead of imageinfo->points[]
	 *     - can export rayshooting
	 *         when shootrays == False no rayshooting is done and *i_points is the array
	 *             of point that are being added, if there are more than one image this
	 *             will be a pointer only to the new points in the first image.  Source
	 *             points are NOT added to s_tree.
	 *
	 *     - Does not affect imageinfo->points in any way
	 *     - Allows for refinements to be done on the borders of the initial grid region
	 *          The initial grid region is never expanded.
	 *
	 * criterion = 0 stops refining when error in total area reaches res_target
	 * 	         = 1 stops refining when each image reaches error limit or is smaller than res_target
	 *           = 2 stops refining when grid resolution is smaller than res_target in all images
	 */

	//printf("entering refine_grid\n");

  if(Nimages < 1) return 0;

  int i,j,number_of_refined,count,kk;
  double rmax,total_area;
  Point *s_points,*point;
  short pass=0;
  long Ncells=0,Ncells_o=0;
  Point *i_points;
  unsigned long Nmarker = 0,Nout = 0;

  for(i=0,total_area=0;i<Nimages;++i) total_area += imageinfo[i].area;

  number_of_refined=0;
  for(i=0,Ncells=0;i<Nimages;++i){
    count=0;

    if(criterion == 0) pass = imageinfo[i].area*imageinfo[i].area_error/total_area > res_target;
    if(criterion == 1) pass = (imageinfo[i].area_error > res_target)*(imageinfo[i].area > 1.0e-5*total_area);
    if(criterion == 2) pass = imageinfo[i].gridrange[1] > res_target;

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
    	if(imageinfo[i].outerborder->Nunits > 0){
    		MoveToTopKist(imageinfo[i].outerborder);
    		do{
    			if( getCurrentKist(imageinfo[i].outerborder)->gridsize > 1.01*rmax/Ngrid_block){
    				// border point is marked to prevent refining more than once
    				//   it will be unmarked by the end of refine grid
    				getCurrentKist(imageinfo[i].outerborder)->in_image = True;
    				++Ncells;
    			}
    		}while( MoveDownKist(imageinfo[i].outerborder) );
    		//printf("   initial outer border refinements count = %li\n",Ncells);
    	}
    }
  }

  Ncells_o = Ncells;

  i_points = NewPointArray((Ngrid_block*Ngrid_block-1)*Ncells,True);
  Nmarker = 0;
  Ncells = 0;

  for(i=0,Ncells=0;i<Nimages;++i){
	  count=0;

	  if(criterion == 0) pass = imageinfo[i].area*imageinfo[i].area_error/total_area > res_target;
	  if(criterion == 1) pass = (imageinfo[i].area_error > res_target)*(imageinfo[i].area > 1.0e-5*total_area);
	  if(criterion == 2) pass = imageinfo[i].gridrange[1] > res_target;

	  // make sure no border point has a lower res than any image point

	  //printf("imageinfo[%i].area_error=%e N=%i\n",i,imageinfo[i].area_error,imageinfo[i].Npoints);
	  if( pass || imageinfo[i].gridrange[0]>1.01*imageinfo[i].gridrange[1]){

		  rmax=MAX(imageinfo[i].gridrange[1],imageinfo[i].gridrange[0]);

       // loop through points in ith image
		//for(j=0,Ncells=0;j<imageinfo[i].Npoints;++j){
		  MoveToTopKist(imageinfo[i].imagekist);
		  for(j=0 ; j<imageinfo[i].imagekist->Nunits ; ++j,MoveDownKist(imageinfo[i].imagekist) ){

			//if( imageinfo[i].points[j].gridsize > 1.01*rmax/Ngrid_block){  /* only refine largest grid size in image*/
			  if( getCurrentKist(imageinfo[i].imagekist)->gridsize > 1.01*rmax/Ngrid_block){  /* only refine largest grid size in image*/

				  // get real point in tree
				  //point = imageinfo[i].points[j].image->image;
				  point = getCurrentKist(imageinfo[i].imagekist);
				  assert(point->gridsize > 0);

				  //imageinfo[i].points[j].gridsize /= Ngrid_block;
				  ++count;

				  xygridpoints(&i_points[Nmarker],point->gridsize*(Ngrid_block-1)/Ngrid_block,point->x,Ngrid_block,1);


				  // check if new points are outside of initial grid region
				  Nout = 0;
				  if( (point->x[0] == i_tree->top->boundery_p1[0]) || (point->x[0] == i_tree->top->boundery_p2[0])
			       || (point->x[1] == i_tree->top->boundery_p1[1]) || (point->x[1] == i_tree->top->boundery_p2[1]) ){

					  // remove the points that are outside initial grid
					  for(kk=0,Nout=0;kk < (Ngrid_block*Ngrid_block-1);++kk){
						  if(!inbox(i_points[Nmarker + kk - Nout].x,i_tree->top->boundery_p1,i_tree->top->boundery_p2)){

							  SwapPointsInArray(&i_points[Nmarker + kk - Nout],&i_points[Nmarker + Ngrid_block*Ngrid_block - 2 - Nout]);
							  ++Nout;

							  //printf("  point taken out\n");
							  //ERROR_MESSAGE();
						  }
					  }
					  assert(Nout > 0);
				  }

				  assert(Nout >= 0);
				  Nmarker += (Ngrid_block*Ngrid_block-1) - Nout;
				  ++Ncells;

				  point->gridsize /= Ngrid_block;
				  point->image->gridsize /= Ngrid_block;

			  }
		  }

		  //printf("   actual image point refinements count = %li\n",Ncells);
		  // * loop through outer border of ith image *

		  MoveToTopKist(imageinfo[i].outerborder);
		  for(j=0;j<imageinfo[i].outerborder->Nunits;++j,MoveDownKist(imageinfo[i].outerborder)){
			  if( getCurrentKist(imageinfo[i].outerborder)->gridsize > 1.01*rmax/Ngrid_block){ // only refine largest grid size in image

				  point = getCurrentKist(imageinfo[i].outerborder);
				  assert(point->gridsize > 0);

				  if(point->in_image){ // point has not been refined yet as border of another image
					  ++count;

					  xygridpoints(&i_points[Nmarker],point->gridsize*(Ngrid_block-1)/Ngrid_block,point->x,Ngrid_block,1);

					  // check if new points are outside of initial grid region
					  Nout = 0;
					  if( (point->x[0] == i_tree->top->boundery_p1[0]) || (point->x[0] == i_tree->top->boundery_p2[0])
				       || (point->x[1] == i_tree->top->boundery_p1[1]) || (point->x[1] == i_tree->top->boundery_p2[1]) ){

						  // remove the points that are outside initial grid
						  for(kk=0,Nout=0;kk<(Ngrid_block*Ngrid_block-1);++kk){
							  if( !inbox(i_points[Nmarker + kk - Nout].x,i_tree->top->boundery_p1,i_tree->top->boundery_p2) ){

								  SwapPointsInArray(&i_points[Nmarker + kk - Nout],&i_points[Nmarker + Ngrid_block*Ngrid_block - 2 - Nout]);
								  ++Nout;

								  //printf("  point taken out\n");
								  //ERROR_MESSAGE();
							  }
						  }
						  assert(Nout > 0);

					  }

					  assert(Nout >= 0);
					  Nmarker += (Ngrid_block*Ngrid_block-1) - Nout;

					  ++Ncells;
					  point->gridsize /= Ngrid_block;
					  point->image->gridsize /= Ngrid_block;
					  point->in_image = False;  // unmak so that it wouldn't double refine
				  }
			  }
		  }


      //printf("should be Ncells=%i Ngrid_block=%i   %i\n",Ncells,Ngrid_block,
      //		(Ngrid_block*Ngrid_block-1)*Ncells);

    }

    if(count > 0) ++number_of_refined;
  } // end of image loop

  assert( Nmarker <= (Ngrid_block*Ngrid_block-1)*Ncells_o );
  if( Nmarker != (Ngrid_block*Ngrid_block-1)*Ncells_o ) i_points = AddPointToArray(i_points,Nmarker,(Ngrid_block*Ngrid_block-1)*Ncells_o);

  s_points=LinkToSourcePoints(i_points,Nmarker);

  if(shootrays){
	  rayshooterInternal(Nmarker,i_points,kappa_off);
  }else{
	  assert(Nimages == 1);
	  // The new points could be put into a kist if there where more than one image
      // for(j=0;j<(Ngrid_block*Ngrid_block-1)*Ncells;++j){
	  //	  i_points[i].image->x[0] = i_points[i].x[0];
	      //	  i_points[i].image->x[1] = i_points[i].x[1];
      //InsertAfterCurrentKist(newkist,i_points[j]);
      //  }
  }

  // add points to trees
  AddPointsToTree(i_tree,i_points,Nmarker);
  if(shootrays) AddPointsToTree(s_tree,s_points,Nmarker);

  if(!shootrays) *point_pnt = i_points;
  return number_of_refined;
}


void findborders4(TreeHndl i_tree,ImageInfo *imageinfo){
	/* finds inner and outer borders of an image using
	 * bordering box method
	 *   uses the in_image markers
	 *   uses imaginfo->imagekist instead of imageinfo->points
	 *
	 *   In the case of the entire grid being within the image, the innerborders
	 *   is the border points of the grid and the outerborder contains no points.
	 *
	 *   Note:  markers in_image must be set
	 */
	int i;
	unsigned long m,j;
	Boolean addinner;
	Boolean allin = False;


	if(i_tree->pointlist->Npoints == imageinfo->imagekist->Nunits) allin = True;    // all points on the grid are in the image

	//printf("beginning findborders2\n");
	//checkTree(i_tree);

	//point=(Point *)malloc(sizeof(Point));
	EmptyKist(imageinfo->innerborder);
	EmptyKist(imageinfo->outerborder);

	imageinfo->gridrange[2] = 1.0e99; /* minimum grid size in image */
	imageinfo->gridrange[0] = 0.0; /* maximum grid size in outerborder */
	imageinfo->gridrange[1] = 0.0;      /* maximum grid size in image */

	if(imageinfo->imagekist->Nunits < 1) return;

	KistHndl neighborkist = NewKist(),imagekist;

	imagekist = imageinfo->imagekist; assert(imageinfo->imagekist);

	MoveToTopKist(imagekist);
	for(j=0;j<imagekist->Nunits;++j,MoveDownKist(imagekist)){

		if(imageinfo->gridrange[1] < getCurrentKist(imagekist)->gridsize)
			imageinfo->gridrange[1] = getCurrentKist(imagekist)->gridsize;
		if(imageinfo->gridrange[2] > getCurrentKist(imagekist)->gridsize)
			imageinfo->gridrange[2] = getCurrentKist(imagekist)->gridsize;

		addinner=False;

		FindAllBoxNeighborsKist(i_tree,getCurrentKist(imagekist),neighborkist);

		if( allin && neighborkist->Nunits < 4){
			InsertAfterCurrentKist(imageinfo->innerborder,getCurrentKist(imagekist));
		}else{

			MoveToTopKist(neighborkist);
			for(i=0;i<neighborkist->Nunits;++i){

				if( getCurrentKist(neighborkist)->in_image == False){  // point is a neighbor
					addinner=True;
					// check if point is already in list
					MoveToTopKist(imageinfo->outerborder);
					for(m=0;m<imageinfo->outerborder->Nunits;++m){
						if( getCurrentKist(imageinfo->outerborder) == getCurrentKist(neighborkist) ) break;
						MoveDownKist(imageinfo->outerborder);
					}
					if(m==imageinfo->outerborder->Nunits){
						// add point to outerborder
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

	if(!allin){
		MoveToTopKist(imageinfo->outerborder);
		do{
			if(imageinfo->gridrange[0] < getCurrentKist(imageinfo->outerborder)->gridsize)
				imageinfo->gridrange[0] = getCurrentKist(imageinfo->outerborder)->gridsize;
		}while(MoveDownKist(imageinfo->outerborder));
	}

	freeKist(neighborkist);

	return;
}

