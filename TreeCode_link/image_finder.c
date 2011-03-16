#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <nrutil.h>
#include "Tree.h"
#include "KistDriver.h"
#include "divide_images.h"

static const int NpointsRequired = 50;  // number of points required to be within an image
static const int Ngrid_block = 3;       // each cell is divided into Ngrid_block^2 subcells

static const float mumin = 0.3;  // actually the sqrt of the minimum magnification
//static const float mumin = 0.2;
static const float FracResTarget = 4.0e-4;
//static const float FracResTarget 1.0e-3

Point *pointg;
double ysourceg[2],magsigng;
static double initialgridsize=0;

void find_images(double *y_source,double r_source
		  ,TreeHndl s_tree,TreeHndl i_tree,int *Nimages
		  ,ImageInfo *imageinfo,const int NimageMax,unsigned long *Nimagepoints
		  ,double initial_size,Boolean splitimages,short edge_refinement
		  ,Boolean verbose,Boolean kappa_off){

	/*
	 * find_image returns finite refined images
	 *  it starts with a large source and reduces down to the right size refining at each step
	 *  should not miss any image larger than ~ munin*r_source linear size
	 *
	 *  if splitimage==TRUE each image is refined to target accuracy, otherwise all images are treated as one
	 *  edge_refinement = 0 does not do edge refinement
	 *                  = 1 uses refine_edge() which keeps the images fully up to date
	 *                  = 2 uses refine_edge2() which is faster but does not keep edges
	 *  kappa_off = True - turns off calculation of kappa and gamma
	 *              False - turns on calculation of kappa and gamma
	 *
	 *  fromthetop = True -  The telescoping steps will start from the largest gird size and
	 *  	                 got down to make sure all images are found
	 *             = False - If the source has not moved and source size is smaller than
	 *                       the last the telescaoping will begin with the last source size
	 *                       ,otherwise the same as true
	 */

	long Nsizes,Nold;
	double rtemp,tmp;
	static double oldy[2],oldr=0;
	short moved;
	int i,j,k,m;
	//Point *i_points,*s_points,*point;
	time_t to,t1,t2,t3,now;
	KistHndl tmp_border_kist;
	Boolean image_overlap;
	static int oldNimages=0;
	static unsigned long Npoints_old = 0;

	if(oldr==0){ oldr=r_source; Npoints_old = i_tree->pointlist->Npoints;}
	if((Npoints_old <= i_tree->pointlist->Npoints )* // if grid has not been refreshed
			(oldy[0]==y_source[0])*(oldy[1]==y_source[1])* // and source not moved
			(oldr > r_source)  // and source size has gotten smaller
			) initial_size=oldr;

	Npoints_old = i_tree->pointlist->Npoints;

	if(r_source==0.0){ERROR_MESSAGE(); printf("ERROR: find_images, point source must have a resolution target\n"); exit(1);}

	// starting with a larger source size make sure all the grid sizes are small enough to find it
	KistHndl subkist = NewKist();
	if(verbose) printf("entering find_image\n");
	time(&to);

	if(verbose) printf("initialgridsize=%e\n",initialgridsize);
	if(initial_size==0) initial_size=initialgridsize;
    Nsizes=(int)(log(initial_size/sqrt(pi)/fabs(r_source*mumin))/log(Ngrid_block) ) + 1.0; // round up
    if(verbose) printf("Ntemp=%li\n",Nsizes);

    ClearAllMarkes(i_tree);

    //////////////////////////////////////////
    // telescope source size down to target
    //////////////////////////////////////////
	if(!( (oldy[0]==y_source[0])*(oldy[1]==y_source[1])*(oldr < r_source) )){

		for(rtemp = fabs(r_source)*pow(Ngrid_block,Nsizes),Nold=0
//				for(rtemp = fabs(r_source/mumin)*pow(Ngrid_block,Nsizes),Nold=0
				;rtemp >= 0.99*Ngrid_block*fabs(r_source)
//				;rtemp > 0.99*fabs(r_source)
				;rtemp /= Ngrid_block ){

			time(&t1);
			time(&t3);
			if(verbose)
				printf("\n   new source size = %e    telescoping rsource=%e\n",rtemp,r_source);

			do{
				time(&t1);
				if(verbose) printf("      time in refine grid %f sec\n",difftime(t1,t2));

				moved=image_finder(y_source,rtemp,s_tree,i_tree
						,Nimages,imageinfo,NimageMax,Nimagepoints,-1,0);

				//printf("%li\n",imageinfo->Npoints);
				//for(i=0;i<imageinfo->Npoints;++i) printf("%e %e\n",imageinfo->points[i].x[0]
				//                                                  ,imageinfo->points[i].x[1]);
				//printf("rtemp=%e\n",rtemp);
				//exit(0);
				//printf("     %li %li images=%i \n",imageinfo->Npoints
				//		,i_tree->pointlist->Npoints,*Nimages);
				time(&t2);
				if(verbose)	printf("      time in image_finder %f sec\n        Nimagepoints=%li\n"
						,difftime(t2,t1),*Nimagepoints);

			}while(refine_grid(i_tree,s_tree,imageinfo,*Nimages,rtemp*mumin/Ngrid_block,2,kappa_off));

			time(&t1);
			if(verbose)	printf("      time in refine grid %f sec\n",difftime(t1,t2));

			time(&now);
			if(verbose) printf("    time for one source size %f sec\n",difftime(now,t3));
		}
	}
	time(&now);
	if(verbose) printf(" time for source size reduction %f sec\n",difftime(now,to));
	time(&to);
	/*
	 moved=image_finder(y_source,fabs(r_source),s_tree,i_tree
			,Nimages,imageinfo,Nimagepoints,0,1);

	printf("  Nimages=%li\n",*Nimages);
	exit(0);
	 */
	//////////////////////////////////////////////////////////////////////////////////
	// target source size has been reached, do full image decomposition and refinement
	/////////////////////////////////////////////////////////////////////////////////
	i=0;
	if(splitimages){

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

			moved=image_finder(y_source,fabs(r_source),s_tree,i_tree
					,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);

			if(*Nimages < 1) printf("  Nimages=%i i=%i\n",*Nimages,i);

			time(&now);
			if(verbose){
				printf("\n    i=%i\n     time in finding images %f min\n          Nimages=%i   Nimagepoints=%li\n"
						,i,difftime(now,t3)/60.,*Nimages,*Nimagepoints);
				printf("     image   # of points    error in area\n");
				for(j=0;j<*Nimages;++j) printf("       %i        %li         %e\n",j,imageinfo[j].Npoints,imageinfo[j].area_error);
			}
			++i;
		}while( refine_grid(i_tree,s_tree,imageinfo,*Nimages,1.0e-2,1,kappa_off)
				|| moved );

		time(&now);

		// remove images with no points in them
		for(j=0;j<*Nimages;++j){
			if(imageinfo[j].Npoints < 1){
				ERROR_MESSAGE();
				for(k=j+1;k<*Nimages;++k){
					imageinfo[k-1].Nencircled=imageinfo[k].Nencircled;
					imageinfo[k-1].Npoints=imageinfo[k].Npoints;
					imageinfo[k-1].area=imageinfo[k].area;
					imageinfo[k-1].area_error=imageinfo[k].area_error;
					imageinfo[k-1].points=imageinfo[k].points;
					for(m=0;m<3;++m) imageinfo[k-1].gridrange[m]=imageinfo[k].gridrange[m];
					tmp_border_kist=imageinfo[k-1].innerborder;
					imageinfo[k-1].innerborder=imageinfo[k].innerborder;
					imageinfo[k].innerborder=tmp_border_kist;
					tmp_border_kist=imageinfo[k-1].outerborder;
					imageinfo[k-1].outerborder=imageinfo[k].outerborder;
					imageinfo[k].outerborder=tmp_border_kist;
				}
				//printf("image %i has no points\n",j);
				--*Nimages;
				--j;
			}
		}
		if(verbose) printf("starting edge refinement\n");

		/////////////////////////////////////////////
		// refine just image edges to high resolution
		/////////////////////////////////////////////
		k=i;
		if(edge_refinement==0){
			do{
				// mark image points in tree
				PointsWithinKist(s_tree,y_source,r_source,subkist,1);

				moved=image_finder(y_source,fabs(r_source),s_tree,i_tree
						,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);
				++i;
			}while( refine_grid(i_tree,s_tree,imageinfo,*Nimages,FracResTarget
					,0,kappa_off)
					|| moved );

		}else if(edge_refinement==1){
			do{
				// mark image points in tree
				PointsWithinKist(s_tree,y_source,r_source,subkist,1);

				moved=image_finder(y_source,fabs(r_source),s_tree,i_tree
						,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);

				//for(i = 0; i < *Nimages; ++i) PrintImageInfo(&(imageinfo[i]));
				//printf("\n");

				++i;
			}while( refine_edges(i_tree,s_tree,imageinfo,*Nimages,FracResTarget,0,kappa_off)
					|| moved );

		}else if(edge_refinement==2){
			++i;
			while(refine_edges2(y_source,r_source,i_tree,s_tree
					,imageinfo,&image_overlap,*Nimages,FracResTarget,0,kappa_off)){
				// if an overlap is detected find the images again

				if(image_overlap) moved=image_finder(y_source,fabs(r_source),s_tree,i_tree
						,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);
				++i;
			}
		}
		// unmark image points so new source can be used
		PointsWithinKist(s_tree,y_source,r_source,subkist,-1);

		if(verbose) printf("finished edge refinement i=%i\n",i);

	}else{ // splitimages == False

		// refine grid to fractional error of total image
		do{
			time(&t3);
			if(verbose && i>0) printf("     time in image refinement %f min\n",difftime(t3,now)/60.);
			// printf("finding images\n    linkinglength=%e\n",linkinglength);

			// mark image points in tree
			PointsWithinKist(s_tree,y_source,r_source,subkist,1);

			moved=image_finder(y_source,fabs(r_source),s_tree,i_tree
					,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);

			if(*Nimages < 1) printf("  Nimages-%i i=%i\n",*Nimages,i);

			time(&now);
			if(verbose) printf("    i=%i\n     time in finding images %f min\n",i,difftime(now,t3)/60.);
			++i;
		}while( refine_grid(i_tree,s_tree,imageinfo,*Nimages,1.0e-2,0,kappa_off)
				|| moved );


		// remove images with no points in them
		for(j=0;j<*Nimages;++j){
			if(imageinfo[j].Npoints < 1){
				ERROR_MESSAGE();
				for(k=j+1;k<*Nimages;++k){
					imageinfo[k-1].Nencircled=imageinfo[k].Nencircled;
					imageinfo[k-1].Npoints=imageinfo[k].Npoints;
					imageinfo[k-1].area=imageinfo[k].area;
					imageinfo[k-1].area_error=imageinfo[k].area_error;
					imageinfo[k-1].points=imageinfo[k].points;
					for(m=0;m<3;++m) imageinfo[k-1].gridrange[m]=imageinfo[k].gridrange[m];
					tmp_border_kist=imageinfo[k-1].innerborder;
					imageinfo[k-1].innerborder=imageinfo[k].innerborder;
					imageinfo[k].innerborder=tmp_border_kist;
					tmp_border_kist=imageinfo[k-1].outerborder;
					imageinfo[k-1].outerborder=imageinfo[k].outerborder;
					imageinfo[k].outerborder=tmp_border_kist;
				}
				//printf("image %i has no points\n",j);
				--*Nimages;
				--j;
			}
		}
		/////////////////////////////////////////////
		// refine just image edges to high resolution
		/////////////////////////////////////////////
		time(&now);

		if(edge_refinement==0){
			do{
				// mark image points in tree
				PointsWithinKist(s_tree,y_source,r_source,subkist,1);

				moved=image_finder(y_source,fabs(r_source),s_tree,i_tree
						,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);
				++i;
			}while( refine_grid(i_tree,s_tree,imageinfo,*Nimages,FracResTarget,0,kappa_off)
					|| moved );

		}else if(edge_refinement==1){
			do{
				// mark image points in tree
				PointsWithinKist(s_tree,y_source,r_source,subkist,1);

				moved=image_finder(y_source,fabs(r_source),s_tree,i_tree
						,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);
				++i;
			}while( refine_edges(i_tree,s_tree,imageinfo,*Nimages,FracResTarget,2,kappa_off)
					|| moved );

		}else if(edge_refinement==2){

			while(refine_edges2(y_source,r_source,i_tree,s_tree
					,imageinfo,&image_overlap,*Nimages,FracResTarget,2,kappa_off)){

			if(image_overlap) moved=image_finder(y_source,fabs(r_source),s_tree,i_tree
					,Nimages,imageinfo,NimageMax,Nimagepoints,0,1);
			++i;
			}
		}
		// unmark image points so new source can be used
		PointsWithinKist(s_tree,y_source,r_source,subkist,-1);

	}
	time(&t3);
	if(verbose) printf("     time in image refinement %f min\n",difftime(t3,now)/60.);

	// if point source take only closest image point
	if(r_source <= 0){

		double r,rmin;

		if(*Nimages > 1){
			for(i=0,*Nimagepoints=0;i<*Nimages;++i){

				// find closest point in each image
				for(j=0,rmin=1.0e99;j<imageinfo[i].Npoints;++j){
					r=sqrt( pow(y_source[0]-imageinfo[i].points[j].image->x[0],2)
							+ pow(y_source[1]-imageinfo[i].points[j].image->x[1],2) );
						  if(r<rmin){
							  rmin=r;
							  PointCopyData(&imageinfo[0].points[*Nimagepoints],&imageinfo[i].points[j]);
						  }
				}
				imageinfo[*Nimagepoints].points=&imageinfo[0].points[*Nimagepoints];
				imageinfo[*Nimagepoints].Npoints=1;
				++*Nimagepoints;
			}
		}
	}

	time(&now);
	if(verbose) printf("time in find_images %f min\n",difftime(now,to)/60.);

	oldy[0]=y_source[0];
	oldy[1]=y_source[1];
	oldr=r_source;

	freeKist(subkist);
	for(i=*Nimages;i<oldNimages;i++){
		EmptyKist(imageinfo[i].innerborder);
		EmptyKist(imageinfo[i].outerborder);
	}
	oldNimages=*Nimages;

	// remove images without points
	for(j=0;j<*Nimages;++j){
		if(imageinfo[j].Npoints < 1){
			assert(imageinfo[j].area == 0);
			for(k=j+1;k<*Nimages;++k){
				imageinfo[k-1].Nencircled=imageinfo[k].Nencircled;
				imageinfo[k-1].Npoints=imageinfo[k].Npoints;
				imageinfo[k-1].area=imageinfo[k].area;
				imageinfo[k-1].area_error=imageinfo[k].area_error;
				imageinfo[k-1].points=imageinfo[k].points;
				for(m=0;m<3;++m) imageinfo[k-1].gridrange[m]=imageinfo[k].gridrange[m];
				tmp_border_kist=imageinfo[k-1].innerborder;
				imageinfo[k-1].innerborder=imageinfo[k].innerborder;
				imageinfo[k].innerborder=tmp_border_kist;
				tmp_border_kist=imageinfo[k-1].outerborder;
				imageinfo[k-1].outerborder=imageinfo[k].outerborder;
				imageinfo[k].outerborder=tmp_border_kist;
			}
			//printf("image %i has no points\n",j);
			--*Nimages;
			--j;
		}
	}
	for(i=0;i<*Nimages;++i){
		tmp=0.0;
		imageinfo[i].centroid[0] = 0.0;
		imageinfo[i].centroid[1] = 0.0;
		for(j=0;j<imageinfo[i].Npoints;++j){
			tmp += pow(imageinfo[i].points[j].gridsize,2);
			imageinfo[i].centroid[0] += imageinfo[i].points[j].x[0]*pow(imageinfo[i].points[j].gridsize,2);
			imageinfo[i].centroid[1] += imageinfo[i].points[j].x[1]*pow(imageinfo[i].points[j].gridsize,2);
		}
		if(imageinfo[i].Npoints >0 ){
			imageinfo[i].centroid[0] /= tmp;
			imageinfo[i].centroid[1] /= tmp;
		}
	}

	   ClearAllMarkes(i_tree);

	   //for(i = 0; i < *Nimages; ++i) PrintImageInfo(&(imageinfo[i]));

	return;
}

short image_finder(double *y_source,double r_source,TreeHndl s_tree,TreeHndl i_tree
		,int *Nimages,ImageInfo *imageinfo,const int NimageMax,unsigned long *Nimagepoints
		,short splitparities,short true_images){
  /* finds images for a given source position and size
  * remember to free output pointer if image_finder is used multiple times
  * imagemarker array must be large enough to include all images *
  * splitparities=  0 don't split attached negative and positive parity images
  *              =  1 do split parities  Note: this is now obsolete
  *              = -1 doesn't slit into images at all
  *	                 , also does not find borders or change in_image markers
  * true_images = 1 gives just the points that are in the image
  *             = 0 if there are not enough points in images this will include close points to be refined
	*/
  unsigned long i,j,Nsource_points=0;
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

  // if(moved) printf("***** new source ******\n");
  // printf("In image_finder\n");

  ++count;

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
			  EmptyKist(imageinfo->imagekist);
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

  *Nimagepoints = imageinfo->imagekist->Nunits;

	  // free points if this not the first time
  if(count>1 && imageinfo->points != NULL){
  //if(count>1 ){
	  free(imageinfo[0].points);
	  Nold_images = *Nimages;
  }

  // transform from source plane to image points
  TranformPlanesKist(imageinfo->imagekist);

  //  x[] is not allocated
  imageinfo->points = NewPointArray(imageinfo->imagekist->Nunits,False);

  // split images
  if( splitparities == 0 ){
	  divide_images2(i_tree,imageinfo,Nimages,NimageMax);
  }else{
	  *Nimages = 1;
	  imageinfo->Npoints = imageinfo->imagekist->Nunits;
	  imageinfo->area = 0.0;
	  MoveToTopKist(imageinfo->imagekist);
	  do{
		  imageinfo->area += pow(getCurrentKist(imageinfo->imagekist)->gridsize,2);
	  }while(MoveDownKist(imageinfo->imagekist));
  }

/*  with divide_images
  // copy information into array
  MoveToTopKist(imageinfo->imagekist);
  for(i=0;i < imageinfo->imagekist->Nunits;++i){
	  //PointCopyData(&(imageinfo->points[i]),getCurrentKist(imageinfo->imagekist));
	  PointCopyData(&(imageinfo->points[i]),TakeOutCurrentKist(imageinfo->imagekist));
	  //MoveDownKist(imageinfo->imagekist);
  }
*/

  // for divide_images2
  // copy information into array
  for(i=0,j=0;i<*Nimages;++i){
	  imageinfo[i].Npoints = imageinfo[i].imagekist->Nunits;

	  MoveToTopKist(imageinfo[i].imagekist);
	  PointCopyData(&(imageinfo->points[j]),TakeOutCurrentKist(imageinfo[i].imagekist));
	  imageinfo[i].points = &(imageinfo->points[j]);
	  ++j;
	  while(imageinfo[i].imagekist->Nunits > 0){
		  PointCopyData(&(imageinfo->points[j]),TakeOutCurrentKist(imageinfo[i].imagekist));
		  ++j;
	  }
  }

  //imageinfo->Npoints = imageinfo->imagekist->Nunits; // ***
/*
  // set first point in each image in the array
  for(i=1,j=0 ; i < *Nimages ; ++i){
	  //printf("  image %i  Npoints = %li \n",i,imageinfo[i].Npoints);
	  j += imageinfo[i-1].Npoints;
	  imageinfo[i].points = &(imageinfo->points[j]);
  }
*/
  //printf("Nimages =%i splitparities == %i r_source = %e\n",*Nimages,splitparities,r_source);
  /*if(splitparities == -1 && r_source > 0.0){

	  *Nimages=1;

  }else{ // *************** split into separate images ********************

	  //printf("Nimages =%i\n",*Nimages);
	  //split_images(i_tree,imageinfo,NimageMax,Nimages,True);
	  split_images2(i_tree,imageinfo,NimageMax,Nimages);    // brute force, but reliable
	  //split_images3(i_tree,imageinfo,NimageMax,Nimages,True);  // split by edge method

	  //printf("Nimages =%i\n\n",*Nimages);
  }*/

  if( splitparities == 0 ) for(i=0;i<*Nimages;++i) findborders3(i_tree,&imageinfo[i]);

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

		  for(j=0;j<imageinfo[i].Npoints;++j){

		    if(imageinfo[i].gridrange[1] < imageinfo[i].points[j].gridsize)
		      imageinfo[i].gridrange[1] = imageinfo[i].points[j].gridsize;
		    if(imageinfo[i].gridrange[2] > imageinfo[i].points[j].gridsize)
		      imageinfo[i].gridrange[2] = imageinfo[i].points[j].gridsize;
		  }
		  imageinfo[i].gridrange[0] = imageinfo[i].gridrange[1];
	  }

	  // find area of images
	  //findarea(&imageinfo[i]);  // ****this is now done in divide_images

	  assert(imageinfo[i].area >= 0.0);
	  if(Nsource_points < NpointsRequired || moved) imageinfo[i].area_error=1.0;
 }

  for(i=*Nimages;i<Nold_images;++i){
	  EmptyKist(imageinfo[i].innerborder);
	  EmptyKist(imageinfo[i].outerborder);
  }

	//for(i=0;i<*Nimages;++i) printf("  image %i  Npoints = %li Ninner = %li Noutter = %li  area = %e\n",i
  //,imageinfo[i].Npoints,imageinfo[i].innerborder->Nunits,imageinfo[i].outerborder->Nunits,imageinfo[i].area);

   return moved;
}


int refine_grid(TreeHndl i_tree,TreeHndl s_tree,ImageInfo *imageinfo
		,unsigned long Nimages,double res_target,short criterion,Boolean kappa_off){

	/* criterion = 0 stops refining when error in total area reaches res_target
	 * 	         = 1 stops refining when each image reaches error limit or is smaller than res_target
	 *           = 2 stops refining when grid resolution is smaller than res_target in all images
	 */

	//printf("entering refine_grid\n");

  if(Nimages < 1) return 0;

  int i,j,number_of_refined,count; /* Ngrid_block must be odd */
  double rmax,total_area;
  Point *i_points,*s_points,*point;
  short pass=0;
  long Ncells,Ncells_o;

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
    	for(j=0;j<imageinfo[i].outerborder->Nunits;++j){
    		if( getCurrentKist(imageinfo[i].outerborder)->gridsize > 1.01*rmax/Ngrid_block ){
    			// border point is marked to prevent refining more than once
    			//   it will be unmarked by the end of refine grid
    			getCurrentKist(imageinfo[i].outerborder)->in_image = True;
    			++Ncells;
    		}
    		MoveDownKist(imageinfo[i].outerborder);
    	}

        i_points=NewPointArray((Ngrid_block*Ngrid_block-1)*Ncells,True);

		Ncells_o=Ncells;
    	//printf("should be Ncells=%i Ngrid_block=%i\n",Ncells,Ngrid_block);

       // loop through points in ith image
    	for(j=0,Ncells=0;j<imageinfo[i].Npoints;++j){

    		if( imageinfo[i].points[j].gridsize > 1.01*rmax/Ngrid_block){  /* only refine largest grid size in image*/

    			// get real point in tree
    			point = imageinfo[i].points[j].image->image;
    			assert(point->gridsize > 0);

    			imageinfo[i].points[j].gridsize /= Ngrid_block;
    			++count;

    			xygridpoints(&i_points[(Ngrid_block*Ngrid_block-1)*Ncells]
    			                       ,point->gridsize*(Ngrid_block-1)/Ngrid_block
    			                       ,imageinfo[i].points[j].x,Ngrid_block,1);

    			++Ncells;

    			point->gridsize /= Ngrid_block;
    			point->image->gridsize /= Ngrid_block;

    		}
    	}

      /* loop through outer border of ith image */

      MoveToTopKist(imageinfo[i].outerborder);
      for(j=0;j<imageinfo[i].outerborder->Nunits;++j){
    	  if( getCurrentKist(imageinfo[i].outerborder)->gridsize > 1.01*rmax/Ngrid_block){/* only refine largest grid size in image*/

    		  point = getCurrentKist(imageinfo[i].outerborder);
    		  assert(point->gridsize > 0);

    		  if(point->in_image){ /* point has not been refined yet */
    			  ++count;

    			  xygridpoints(&i_points[(Ngrid_block*Ngrid_block-1)*Ncells]
    			                         ,point->gridsize*(Ngrid_block-1)/Ngrid_block
    			                         ,point->x,Ngrid_block,1);

    			  if( inbox(i_points[(Ngrid_block*Ngrid_block-1)*Ncells ].x
    					  ,i_tree->top->boundery_p1,i_tree->top->boundery_p2) == 0 ){
    				  ERROR_MESSAGE();
    			  }
    			  ++Ncells;
       			  point->gridsize /= Ngrid_block;
       			  point->image->gridsize /= Ngrid_block;
       			  point->in_image = False;  // unmak so that it wouldn't bouble refine
    		  }
    		  //imageinfo[i].outerborderlist->current->gridsize /= Ngrid_block;/**/
    	  }
    	  MoveDownKist(imageinfo[i].outerborder);
      }

      if(Ncells != Ncells_o){
       	  i_points=AddPointToArray(i_points,(Ngrid_block*Ngrid_block-1)*Ncells
       			  ,(Ngrid_block*Ngrid_block-1)*Ncells_o);
      }

      s_points=LinkToSourcePoints(i_points,(Ngrid_block*Ngrid_block-1)*Ncells);

      // lens the added points
      //printf("refine_grid\n");
      rayshooterInternal((Ngrid_block*Ngrid_block-1)*Ncells,i_points,i_tree,kappa_off);

      //printf("should be Ncells=%i Ngrid_block=%i   %i\n",Ncells,Ngrid_block,
      //		(Ngrid_block*Ngrid_block-1)*Ncells);

      // add points to trees
      //printf("refine grid\n");
		AddPointsToTree(i_tree,i_points,Ncells*(Ngrid_block*Ngrid_block-1));
		//printf("   s-plane\n");
		AddPointsToTree(s_tree,s_points,Ncells*(Ngrid_block*Ngrid_block-1));
    }

    if(count > 0) ++number_of_refined;
  } /* end of image loop */

   return number_of_refined;
}

long refine_edges(TreeHndl i_tree,TreeHndl s_tree,ImageInfo *imageinfo
		,unsigned long Nimages,double res_target,short criterion,Boolean kappa_off){
	/*	refines only inner and outer edges of all images
	 * criterion = 0 stops refining when each image reaches error limit
	 *           = 1 stops refining when grid resolution is smaller than res_target in all images
	 *           = 2 stop when area of a cell reaches res_target * area of all images
	 */
	 //printf("entering refine_edges\n");

	if(Nimages < 1) return 0;

	long i,j,Ncells=0,Ncells_o=0,count=0;
	Point *i_points,*s_points,*point;
	double area_total=0;

	// count border points
	if( criterion == 2 ) for(i=0,area_total=0.0;i<Nimages;++i) area_total += imageinfo[i].area;

	for(i=0,Ncells=0;i<Nimages;++i){

		MoveToTopKist(imageinfo[i].outerborder);
		for(j=0;j<imageinfo[i].outerborder->Nunits;++j){
			if( ( criterion==0 && pow(getCurrentKist(imageinfo[i].outerborder)->gridsize,2)/imageinfo[i].area > res_target )
							|| ( criterion==1 && getCurrentKist(imageinfo[i].outerborder)->gridsize > res_target)
							|| ( criterion==2 && pow(getCurrentKist(imageinfo[i].outerborder)->gridsize,2)/area_total > res_target) ){
				++Ncells;
				getCurrentKist(imageinfo[i].outerborder)->in_image = True;  // Temporarily mark point so they are not double refined
			}
			//if( criterion==0 && pow( getCurrentKist(imageinfo[i].outerborder)->gridsize,2)/imageinfo[i].area > res_target) ++Ncells;
			//if( criterion==1 && getCurrentKist(imageinfo[i].outerborder)->gridsize > res_target) ++Ncells;
			//if( criterion==2 && pow(getCurrentKist(imageinfo[i].outerborder)->gridsize,2)/area_total > res_target) ++Ncells;
			MoveDownKist(imageinfo[i].outerborder);
		}
		MoveToTopKist(imageinfo[i].innerborder);
		for(j=0;j<imageinfo[i].innerborder->Nunits;++j){
			if( criterion==0 && pow(getCurrentKist(imageinfo[i].innerborder)->gridsize,2)/imageinfo[i].area > res_target ) ++Ncells;
			if( criterion==1 && getCurrentKist(imageinfo[i].innerborder)->gridsize > res_target) ++Ncells;
			if( criterion==2 && pow(getCurrentKist(imageinfo[i].innerborder)->gridsize,2)/area_total > res_target ) ++Ncells;
			MoveDownKist(imageinfo[i].innerborder);
		}

		//printf("       %i Nouter=%i Ninner=%i\n",i,imageinfo[i].outerborder->Npoints,imageinfo[i].innerborderkist->Npoints);
	}

	if(Ncells==0) return 0;

	i_points=NewPointArray((Ngrid_block*Ngrid_block-1)*Ncells,True);
	Ncells_o=Ncells;

	for(i=0,Ncells=0;i<Nimages;++i){
			/* loop through outer border of ith image */

		MoveToTopKist(imageinfo[i].outerborder);
		for(j=0;j<imageinfo[i].outerborder->Nunits;++j){

			//if( imageinfo[i].outerborder->current->gridsize > res_target){
			if( ( criterion==0 && pow(getCurrentKist(imageinfo[i].outerborder)->gridsize,2)/imageinfo[i].area > res_target )
				|| ( criterion==1 && getCurrentKist(imageinfo[i].outerborder)->gridsize > res_target)
				|| ( criterion==2 && pow(getCurrentKist(imageinfo[i].outerborder)->gridsize,2)/area_total > res_target) ){

				//point=imageinfo[i].outerborderlist->current->leaf->points;
				point = getCurrentKist(imageinfo[i].outerborder);
    			assert(point->gridsize > 0);

				if(point->in_image){ // point has not been refined yet
					++count;

					xygridpoints(&i_points[(Ngrid_block*Ngrid_block-1)*Ncells]
						              ,point->gridsize*(Ngrid_block-1)/Ngrid_block
					                  ,point->x,Ngrid_block,1);

					point->gridsize /= Ngrid_block;
					point->image->gridsize /= Ngrid_block;

					point->in_image = False;

					++Ncells;
				}
				//imageinfo[i].outerborderlist->current->gridsize /= Ngrid_block;/**/
			}
			MoveDownKist(imageinfo[i].outerborder);
		}

		MoveToTopKist(imageinfo[i].innerborder);
		for(j=0;j<imageinfo[i].innerborder->Nunits;++j){

			//if( getCurrentKist(imageinfo[i].innerborderkist)->gridsize > res_target){
			if( ( criterion==0 && pow(getCurrentKist(imageinfo[i].innerborder)->gridsize,2)/imageinfo[i].area > res_target )
					|| (criterion==1 && getCurrentKist(imageinfo[i].innerborder)->gridsize > res_target)
					|| (criterion==2 && pow(getCurrentKist(imageinfo[i].innerborder)->gridsize,2)/area_total > res_target) ){

				point = getCurrentKist(imageinfo[i].innerborder);
    			assert(point->gridsize > 0);

				//if( getCurrentKist(imageinfo[i].innerborderkist)->gridsize == point->gridsize){ /* point has not been refined yet */
    			++count;

    			xygridpoints(&i_points[(Ngrid_block*Ngrid_block-1)*Ncells]
    			                       ,point->gridsize*(Ngrid_block-1)/Ngrid_block
    			                       ,point->x,Ngrid_block,1);

				point->gridsize /= Ngrid_block;
    			point->image->gridsize /= Ngrid_block;

       			++Ncells;
				//}
				//getCurrentKist(imageinfo[i].innerborderkist)->gridsize /= Ngrid_block;
			}
			MoveDownKist(imageinfo[i].innerborder);
		}

	}

    if(Ncells != Ncells_o){
     	  i_points=AddPointToArray(i_points,(Ngrid_block*Ngrid_block-1)*Ncells
     			  ,(Ngrid_block*Ngrid_block-1)*Ncells_o);
    }

	/* lens the added points */
	s_points=LinkToSourcePoints(i_points,(Ngrid_block*Ngrid_block-1)*Ncells);
	//printf("refine_edges\n");
	rayshooterInternal((Ngrid_block*Ngrid_block-1)*Ncells,i_points,i_tree,kappa_off);

	/* add points to trees */
	//printf("edges1\n");
	AddPointsToTree(i_tree,i_points,Ncells*(Ngrid_block*Ngrid_block-1));
	//printf("   s-plane\n");
	AddPointsToTree(s_tree,s_points,Ncells*(Ngrid_block*Ngrid_block-1));

	return count;
}

long refine_edges2(double *y_source,double r_source,TreeHndl i_tree,TreeHndl s_tree
		,ImageInfo *imageinfo,Boolean *image_overlap,unsigned long Nimages,double res_target
		,short criterion,Boolean kappa_off){
	/*	refines only inner and outer edges of all images
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
	 *       *image_overlap = True if the images merge, but not guaranteed when they separate.
	 *
	 * criterion = 0 stops refining when each image reaches error limit or is smaller than res_target
	 *           = 1 stops refining when grid resolution is smaller than res_target in all images
	 *           = 2 stop when area of a cell reaches res_target * area of all images
	 */
	 //printf("entering refine_edges2\n");

	if(Nimages < 1) return 0;

	long i,j,k,n,Ncells=0,Ncells_o=0,count=0,Npoints=0;
	Point *i_points,*s_points,*point;
	KistHndl neighborkist = NewKist();
	double tmp_area=0,area_total=0;
	Boolean addinner;

	// count border points
	if( criterion==2) for(i=0,area_total = 0.0;i<Nimages;++i) area_total += imageinfo[i].area;

	*image_overlap=False;

	// loop through outer border of ith image

	for(i=0;i<Nimages;++i){

		// count border points that needs to be refined
		MoveToTopKist(imageinfo[i].outerborder);
		for(j=0,Ncells=0;j<imageinfo[i].outerborder->Nunits;++j){
			if( ( criterion==0 && pow( getCurrentKist(imageinfo[i].outerborder)->gridsize,2)/imageinfo[i].area > res_target )
					|| ( criterion==1 && getCurrentKist(imageinfo[i].outerborder)->gridsize > res_target)
					|| ( criterion==2 && pow( getCurrentKist(imageinfo[i].outerborder)->gridsize,2)/area_total > res_target)){
				++Ncells;
				getCurrentKist(imageinfo[i].outerborder)->in_image = True; // temporary mark
			}
			/*if( criterion==0 && pow( getCurrentKist(imageinfo[i].outerborder)->gridsize,2)/imageinfo[i].area > res_target ) ++Ncells;
			if( criterion==1 && getCurrentKist(imageinfo[i].outerborder)->gridsize > res_target) ++Ncells;
			if( criterion==2 && pow( getCurrentKist(imageinfo[i].outerborder)->gridsize,2)/area_total > res_target) ++Ncells;*/
			MoveDownKist(imageinfo[i].outerborder);
		}
		MoveToTopKist(imageinfo[i].innerborder);
		for(j=0;j<imageinfo[i].innerborder->Nunits;++j){
			if( criterion==0 && pow(getCurrentKist(imageinfo[i].innerborder)->gridsize,2)/imageinfo[i].area > res_target ) ++Ncells;
			if( criterion==1 && getCurrentKist(imageinfo[i].innerborder)->gridsize > res_target) ++Ncells;
			if( criterion==2 && pow(getCurrentKist(imageinfo[i].innerborder)->gridsize,2)/area_total > res_target ) ++Ncells;
			MoveDownKist(imageinfo[i].innerborder);
		}
		//printf("Ncells=%i total \n",Ncells);
		//printf("       %i Nouter=%i Ninner=%i\n",i,imageinfo[i].outerborder->Npoints,imageinfo[i].innerborderkist->Nunits);

		EmptyKist(neighborkist);

		if(Ncells>0){

			i_points=NewPointArray((Ngrid_block*Ngrid_block-1)*Ncells,True);
			Ncells_o=Ncells;

			MoveToTopKist(imageinfo[i].outerborder);
			for(j=0,Ncells=0;j<imageinfo[i].outerborder->Nunits;++j){

			//if( imageinfo[i].outerborder->current->gridsize > res_target){
				if( ( criterion==0 && pow( getCurrentKist(imageinfo[i].outerborder)->gridsize,2)/imageinfo[i].area > res_target )
						|| ( criterion==1 && getCurrentKist(imageinfo[i].outerborder)->gridsize > res_target)
						|| ( criterion==2 && pow( getCurrentKist(imageinfo[i].outerborder)->gridsize,2)/area_total > res_target)){

					point = getCurrentKist(imageinfo[i].outerborder);
	    			assert(point->gridsize > 0);

	                // point has not been refined yet
					if(point->in_image){
						++count;

						xygridpoints(&i_points[(Ngrid_block*Ngrid_block-1)*Ncells]
					                       ,point->gridsize*(Ngrid_block-1)/Ngrid_block
					                       ,point->x,Ngrid_block,1);

						++Ncells;
						point->gridsize /= Ngrid_block;
						point->image->gridsize /= Ngrid_block;

						point->in_image = False;
					}else *image_overlap=True;

				}
				MoveDownKist(imageinfo[i].outerborder);
			}
			//printf("Ncells=%i in outer border the second time \n",Ncells);

			MoveToTopKist(imageinfo[i].innerborder);
			tmp_area=imageinfo[i].area;
			for(j=0;j<imageinfo[i].innerborder->Nunits;++j){

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

						xygridpoints(&i_points[(Ngrid_block*Ngrid_block-1)*Ncells]
					                       ,point->gridsize*(Ngrid_block-1)/Ngrid_block
					                       ,point->x,Ngrid_block,1);

						++Ncells;
						point->gridsize /= Ngrid_block;
						point->image->gridsize /= Ngrid_block;

						// subtract area of new points from image area
						imageinfo[i].area += pow(point->gridsize,2)*(1-Ngrid_block*Ngrid_block);

					//}else *image_overlap = True;

					if(imageinfo[i].gridrange[2] > point->gridsize)
						imageinfo[i].gridrange[2] = point->gridsize;

					//getCurrentKist(imageinfo[i].innerborderkist)->gridsize /= Ngrid_block;
				}
				MoveDownKist(imageinfo[i].innerborder);
			}

		//printf("image %i Ncells=%i after splitting process\n",i,Ncells);
		//printf("image %i area = %e  inner %i outer %i\n",i,imageinfo[i].area,imageinfo[i].innerborderkist->Nunits,
		//		imageinfo[i].outerborder->Npoints);
		//if(tmp_area==imageinfo[i].area){printf("no subtraction of area\n"); exit(0);}

			if(Ncells != Ncells_o){
				i_points=AddPointToArray(i_points,(Ngrid_block*Ngrid_block-1)*Ncells
						,(Ngrid_block*Ngrid_block-1)*Ncells_o);
			}

			//if(Ncells == 0){
			// lens the added points
			s_points=LinkToSourcePoints(i_points,(Ngrid_block*Ngrid_block-1)*Ncells);
			rayshooterInternal((Ngrid_block*Ngrid_block-1)*Ncells,i_points,i_tree,kappa_off);

			// add points to trees
			AddPointsToTree(i_tree,i_points,Ncells*(Ngrid_block*Ngrid_block-1));
			AddPointsToTree(s_tree,s_points,Ncells*(Ngrid_block*Ngrid_block-1));

			//  Update image borders and area

			if(Ncells) imageinfo[i].gridrange[0]/=Ngrid_block; /* maximum grid size in outerborder */

			//  sort new points into in and out of image
			//    and add them to inner and outer borders
			for(j=0;j<(Ngrid_block*Ngrid_block-1)*Ncells;++j){
				if( sqrt(pow(i_points[j].image->x[0]-y_source[0],2)
						+ pow(i_points[j].image->x[1]-y_source[1],2)) < r_source){

					// mark points
					i_points[j].in_image = True;
					i_points[j].image->in_image = True;

					InsertAfterCurrentKist(imageinfo[i].innerborder,&(i_points[j]));
					MoveDownKist(imageinfo[i].innerborder);
					//PointCopyData(imageinfo[i].innerborder->current,&(i_points[j]));

					imageinfo[i].area += pow(i_points[j].gridsize,2);
				}else{

					// un-mark points
					i_points[j].in_image = False;
					i_points[j].image->in_image = False;

				}
			}

			/*
			 *     Weed out of the borders the points that are no longer on the border
			 */

			//if(Ncells){
			EmptyKist(imageinfo[i].outerborder);

			Npoints=imageinfo[i].innerborder->Nunits;
			MoveToTopKist(imageinfo[i].innerborder);
			for(j=0;j<Npoints;++j){
				addinner=False;

				// update leaf pointer of inner border point if necessary
				//  *** don't think this is necessary anymore
				if(getCurrentKist(imageinfo[i].innerborder)->leaf->npoints > 1){
					i_tree->current = getCurrentKist(imageinfo[i].innerborder)->leaf;
					_FindBox(i_tree,getCurrentKist(imageinfo[i].innerborder)->x);
					getCurrentKist(imageinfo[i].innerborder)->leaf = i_tree->current;
				}

				FindAllBoxNeighborsKist(i_tree,getCurrentKist(imageinfo[i].innerborder),neighborkist);

				MoveToTopKist(neighborkist);
				for(k=0;k < neighborkist->Nunits ;++k){
					if(getCurrentKist(neighborkist)->in_image == False){
						addinner=True;

						MoveToTopKist(imageinfo[i].outerborder);
						for(n=0;n<imageinfo[i].outerborder->Nunits;++n){
							if(getCurrentKist(imageinfo[i].outerborder) == getCurrentKist(neighborkist)) break;
								MoveDownKist(imageinfo[i].outerborder);
							}
							// add to outer border
						if(n==imageinfo[i].outerborder->Nunits){
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
	}  // end loop through images

	freeKist(neighborkist);

	return count;
}

void findborders2(TreeHndl i_tree,ImageInfo *imageinfo){
	/* finds inner and outer borders of an image using
	 * bordering box method
	 */
	int i;
	unsigned long l,m,j;
	Boolean addinner;

	//printf("beginning findborders2\n");
	//checkTree(i_tree);

	//point=(Point *)malloc(sizeof(Point));
	EmptyKist(imageinfo->innerborder);
	EmptyKist(imageinfo->outerborder);

	imageinfo->gridrange[2] = 1.0e99; /* minimum grid size in image */
	imageinfo->gridrange[0] = 0.0; /* maximum grid size in outerborder */
	imageinfo->gridrange[1] = 0.0;      /* maximum grid size in image */

	if(imageinfo->Npoints < 1) return;

	KistHndl neighborkist = NewKist();

	for(j=0;j<imageinfo->Npoints;++j){

		if(imageinfo->gridrange[1] < imageinfo->points[j].gridsize)
			imageinfo->gridrange[1] = imageinfo->points[j].gridsize;
		if(imageinfo->gridrange[2] > imageinfo->points[j].gridsize)
			imageinfo->gridrange[2] = imageinfo->points[j].gridsize;

		addinner=False;

		FindAllBoxNeighborsKist(i_tree,&(imageinfo->points[j]),neighborkist);

		MoveToTopKist(neighborkist);
		for(i=0;i<neighborkist->Nunits;++i){

			for(l=0;l<imageinfo->Npoints;++l) if( getCurrentKist(neighborkist)->id
					== imageinfo->points[l].id) break;

			if(l==imageinfo->Npoints){  // point is a neighbor
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

						addinner=True;
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
	freeKist(neighborkist);

	return;
}

void findborders3(TreeHndl i_tree,ImageInfo *imageinfo){
	/* finds inner and outer borders of an image using
	 * bordering box method
	 *   same as findborder2() but it uses the in_image markers
	 *   which makes it faster
	 *   Note:  markers in_image must be set
	 */
	int i;
	unsigned long m,j;
	Boolean addinner;

	//printf("beginning findborders2\n");
	//checkTree(i_tree);

	//point=(Point *)malloc(sizeof(Point));
	EmptyKist(imageinfo->innerborder);
	EmptyKist(imageinfo->outerborder);

	imageinfo->gridrange[2] = 1.0e99; /* minimum grid size in image */
	imageinfo->gridrange[0] = 0.0; /* maximum grid size in outerborder */
	imageinfo->gridrange[1] = 0.0;      /* maximum grid size in image */

	if(imageinfo->Npoints < 1) return;

	KistHndl neighborkist = NewKist();

	for(j=0;j<imageinfo->Npoints;++j){

		if(imageinfo->gridrange[1] < imageinfo->points[j].gridsize)
			imageinfo->gridrange[1] = imageinfo->points[j].gridsize;
		if(imageinfo->gridrange[2] > imageinfo->points[j].gridsize)
			imageinfo->gridrange[2] = imageinfo->points[j].gridsize;

		addinner=False;

		FindAllBoxNeighborsKist(i_tree,&(imageinfo->points[j]),neighborkist);

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
			InsertAfterCurrentKist(imageinfo->innerborder,imageinfo->points[j].image->image);  // need to put in real point
			MoveDownKist(imageinfo->innerborder);
		}
	}

	freeKist(neighborkist);

	return;
}
void findborders4(TreeHndl i_tree,ImageInfo *imageinfo){
	/* finds inner and outer borders of an image using
	 * bordering box method
	 *   uses the in_image markers
	 *   uses imaginfo->imagekist instead of
	 *   which makes it faster
	 *   Note:  markers in_image must be set
	 */
	int i;
	unsigned long m,j;
	Boolean addinner;

	//printf("beginning findborders2\n");
	//checkTree(i_tree);

	//point=(Point *)malloc(sizeof(Point));
	EmptyKist(imageinfo->innerborder);
	EmptyKist(imageinfo->outerborder);

	imageinfo->gridrange[2] = 1.0e99; /* minimum grid size in image */
	imageinfo->gridrange[0] = 0.0; /* maximum grid size in outerborder */
	imageinfo->gridrange[1] = 0.0;      /* maximum grid size in image */

	if(imageinfo->Npoints < 1) return;

	KistHndl neighborkist = NewKist(),imagekist;

	imagekist = imageinfo->imagekist;

	MoveToTopKist(imagekist);
	for(j=0;j<imagekist->Nunits;++j,MoveDownKist(imagekist)){

		if(imageinfo->gridrange[1] < getCurrentKist(imagekist)->gridsize)
			imageinfo->gridrange[1] = getCurrentKist(imagekist)->gridsize;
		if(imageinfo->gridrange[2] > getCurrentKist(imagekist)->gridsize)
			imageinfo->gridrange[2] = getCurrentKist(imagekist)->gridsize;

		addinner=False;

		FindAllBoxNeighborsKist(i_tree,&(imageinfo->points[j]),neighborkist);

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

	freeKist(neighborkist);

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

      if( (2*(i/Ngrid_1d)/(Ngrid_1d-1) == 1) && (i%Ngrid_1d== Ngrid_1d/2+1) ) j=1;
      i_points[i-j].id=id;
      ++id;
      PositionFromIndex(i,i_points[i-j].x,Ngrid_1d,range,center);
      //i_points[i-j].x[0] = center[0] + range*( 1.0*(i%Ngrid_1d)/(Ngrid_1d-1) - 0.5 );
      //i_points[i-j].x[1] = center[1] + range*( 1.0*(i/Ngrid_1d)/(Ngrid_1d-1) - 0.5 );
      i_points[i-j].gridsize=range/(Ngrid_1d-1);
    }

  }else{
    /*i_points=NewPointArray(Ngrid_1d*Ngrid_1d);*/
    for(i=0;i<Ngrid_1d*Ngrid_1d;++i){
      i_points[i].id=id;
      ++id;
      PositionFromIndex(i,i_points[i].x,Ngrid_1d,range,center);
      //i_points[i].x[0] = center[0] + range*( 1.0*(i%Ngrid_1d)/(Ngrid_1d-1) - 0.5 );
      //i_points[i].x[1] = center[1] + range*( 1.0*(i/Ngrid_1d)/(Ngrid_1d-1) - 0.5 );
      i_points[i].gridsize=range/(Ngrid_1d-1);
    }
  }

  return;
}

void combineCloseImages(double linkinglength,ImageInfo *imageinfo,int *Nimages
		,int *NewNimages){
	unsigned long i,j,k;

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

void SwapImages(ImageInfo *image1,ImageInfo *image2){
	Point *point;
	unsigned long Npoints,i;
	double tmp;
	KistHndl list;

	point = image1->points;
	image1->points = image2->points;
	image2->points = point;

	Npoints = image1->Npoints;
	image1->Npoints = image2->Npoints;
	image2->Npoints = Npoints;
	Npoints = image1->Nencircled;
	image1->Nencircled = image2->Nencircled;
	image2->Nencircled = Npoints;

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

	list = image1->innerborder;
	image1->innerborder = image2->innerborder;
	image2->innerborder = list;

	list = image1->outerborder;
	image1->outerborder = image2->outerborder;
	image2->outerborder = list;

	return ;
}
