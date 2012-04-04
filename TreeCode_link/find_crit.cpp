/*
 * find_crit.c
 *
 *  Created on: Sep 8, 2009
 *      Author: R.B. Metcalf
 */

#include <slsimlib.h>

#define NMAXCRITS 1000

/** \ingroup ImageFinding
  *
  * \brief Finds critical curves and caustics.
  *
  * OUTPUT: each critical curve is in a array of IamgeInfo's
  *         result.parity = 1 tangential caustic, 2 radial, 0 not enough points to determine
  *  the inner out outer boundaries of the result are the estimated critical curves
  *
  * This routine needs to be updated.  It still uses List instead of kist and perhaps some
  * older and slower image splitting and ordering.
  */
ImageInfo *find_crit(
	    LensHndl lens
		,GridHndl grid             /// The grid.  It must be initialized.
		,int *Ncrits              /// The number of critical curves found.
		,double resolution        /// The target resolution that the critical curve is mapped on the image plane.
		,bool *orderingsuccess    /// true if ordering was successful.
		,bool ordercurve          /// Order the curve so that it can be drawn or used to find the winding number.
		,bool verbose
		){

  Point *minpoint;
  ImageInfo *critexport;
  //unsigned long j,k,m,jold;
  unsigned long Npoints,i=0,j;
  short refinements;
  //short spur,closed;
  double maxgridsize,mingridsize,x[2];
  ListHndl negpointlist;

  negpointlist=NewList();
  minpoint=NewPoint(x,0);
  minpoint->invmag=1.0e99;
  /*point=NmewPoint(x,0);*/

  OldImageInfo *critcurve = new OldImageInfo[NMAXCRITS];


  for(;;){

	  // find list of points with negative magnification
	  EmptyList(negpointlist);
	  MoveToTopList(grid->i_tree->pointlist);
	  for(i=0,minpoint->kappa=0;i<grid->i_tree->pointlist->Npoints;++i){
		  if(grid->i_tree->pointlist->current->invmag < 0){

			  InsertAfterCurrent(negpointlist,grid->i_tree->pointlist->current->x,grid->i_tree->pointlist->current->id
					  ,grid->i_tree->pointlist->current->image);
			  MoveDownList(negpointlist);
			  PointCopyData(negpointlist->current,grid->i_tree->pointlist->current);
		  }

		  // record point of maximum kappa
		  if(grid->i_tree->pointlist->current->kappa > minpoint->kappa) PointCopyData(minpoint,grid->i_tree->pointlist->current);
		  MoveDownList(grid->i_tree->pointlist);
	  }

	  Npoints=negpointlist->Npoints;
	  critcurve[0].Npoints=Npoints;

	  //std::cout << "Npoints " << Npoints << "\n";

	  if(Npoints == 0){
		  if(minpoint->gridsize <= resolution){  // no caustic found at this resolution
			  *Ncrits=0;
			  delete[] critcurve;
			  critexport = new ImageInfo[1];
			  return critexport;
		  }

		  /* if there is no negative magnification points use maximum mag point */
		  ++Npoints;
		  //critcurve[0].points=(Point *) malloc(sizeof(Point));
		  //critcurve[0].points->head=1;
		  critcurve[0].points=NewPointArray(1,false);
		  critcurve[0].points->in_image = FALSE;
		  critcurve[0].Npoints=1;
		  PointCopyData(critcurve[0].points,minpoint);
	  }else{
		  /*critcurve[0].points=(Point *) malloc(Npoints*sizeof(Point));
		  critcurve[0].Npoints=Npoints;
		  critcurve[0].points->head=Npoints;
		   */

		  critcurve[0].points=NewPointArray(Npoints,false);

		  MoveToTopList(negpointlist);
		  for(i=0;i<negpointlist->Npoints;++i){
			  PointCopyData(&(critcurve[0].points[i]),negpointlist->current);
    	  //std::printf("     negpointlist = %e %e \n       critcurve[0].point[%i]=%e %e\n",negpointlist->current->x[0]
    	 // 	       ,negpointlist->current->x[1],i,critcurve[0].points[i].x[0],critcurve[0].points[i].x[1]);
			  MoveDownList(negpointlist);
		  }
	  }

    if(verbose) std::printf("find_crit, going into findborders 1\n");
    findborders2(grid->i_tree,&critcurve[0]);
    if(verbose) std::printf("find_crit, came out of findborders 1\n");

    /* make inner border the image */
    MoveToTopKist(critcurve[0].innerborder);
    for(i=0,maxgridsize=0.0,mingridsize=1.0e99;i<critcurve[0].innerborder->Nunits();++i){
    	PointCopyData(&(critcurve[0].points[i]),getCurrentKist(critcurve[0].innerborder));
    	if(critcurve[0].points[i].gridsize > maxgridsize) maxgridsize=critcurve[0].points[i].gridsize;
    	if(critcurve[0].points[i].gridsize < mingridsize) mingridsize=critcurve[0].points[i].gridsize;
    	MoveDownKist(critcurve[0].innerborder);
    }
    critcurve[0].Npoints=critcurve[0].innerborder->Nunits();

    // find borders again to properly define outer border
    //std::printf("going into findborders 2\n");
	findborders2(grid->i_tree,critcurve);
	//std::printf("came out of findborders 2\n");

	if(verbose) std::printf("find_crit, going into refine_grid\n");
     //std::printf("  Npoints=%i\n",critcurve->Npoints);
	refinements=refine_grid(lens,grid->i_tree,grid->s_tree,critcurve,1,resolution,2,false);
    if(verbose) std::printf("find_crit, came out of refine_grid\n");

    if(verbose) cout << "Npoints " << critcurve[0].Npoints << endl;

     if(refinements==0){
      break;
    }else free(critcurve[0].points);
  }

  if(verbose) std::printf("find_crit, number of caustic points: %li\n",critcurve[0].Npoints);

// order points in curve
  if(ordercurve) split_order_curve4(critcurve,NMAXCRITS,Ncrits);
  else if(critcurve->Npoints > 0) *Ncrits=1;
  if(critcurve->Npoints == 0) *Ncrits=0;

/*
//   print out the critical curves and caustics
  std::printf("Ncrits=%i\n",*Ncrits);
  for(j=0;j<*Ncrits;++j){
	std::printf("%li\n",critcurve[j].Npoints);
	for(i=0;i<critcurve[j].Npoints;++i)
		std::printf("%e %e\n",critcurve[j].points[i].x[0]
	                    ,critcurve[j].points[i].x[1]);
  }
  std::printf("Ncrits=%i\n",*Ncrits);
  for(j=0;j<*Ncrits;++j){
	std::printf("%li\n",critcurve[j].Npoints);
	for(i=0;i<critcurve[j].Npoints;++i)
		std::printf("%e %e\n",critcurve[j].points[i].image->x[0]
	                    ,critcurve[j].points[i].image->x[1]);
  }
  exit(0);
	*/

  if(*Ncrits==0 && critcurve->Npoints > 0 ){
	  *Ncrits=1;
	  *orderingsuccess=false;
  }else{ *orderingsuccess=true;}

  if(*Ncrits > NMAXCRITS){ERROR_MESSAGE(); std::printf("ERROR: in find_crit, too many critical curves Ncrits=%i > NMAXCRITS gridsize=%e\n"
			       ,*Ncrits,critcurve[0].points[0].gridsize); exit(1);}

  /* find area of critical curves */
  x[0]=x[1]=0.0;
  for(i=0;i<*Ncrits;++i){
	  if(critcurve[i].Npoints > 5){
			   windings(x,critcurve[i].points,critcurve[i].Npoints,&(critcurve[i].area),0);
			   //std::printf("critarea = %e\n",critcurve[i].area);
	  }else critcurve[i].area=0;
  }

  EmptyList(negpointlist);
  free(negpointlist);
  free(minpoint);

  // resize crit array so it doesn't use more mem than necessary
  critexport= new ImageInfo[*Ncrits];
  for(i=0;i<*Ncrits;++i){
	  critexport[i].ShouldNotRefine = critcurve[i].ShouldNotRefine;
	  //critexport[i].Npoints=critcurve[i].Npoints;
	  critexport[i].area = critcurve[i].area;
	  critexport[i].area_error = critcurve[i].area_error;
	  //critexport[i].points=critcurve[i].points;

	  for(j=0;j<critcurve[i].Npoints;++j){
		  critexport[i].imagekist->InsertAfterCurrent(&(critcurve[i].points[j]));
	  }
  }

  //freeImageInfo(critcurve,NMAXCRITS);
  delete[] critcurve;
  return critexport;
}

