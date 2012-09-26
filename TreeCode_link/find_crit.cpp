/*
 * find_crit.c
 *
 *  Created on: Sep 8, 2009
 *      Author: R.B. Metcalf
 */

#include <slsimlib.h>

#define NMAXCRITS 1000

using namespace std;

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
  *
  * The critical curve is found by refining the edges of regions of negative magnification.
  * If there are no regions of negative magnification in the original grid the grid is refined
  * around the point of highest kappa.  If there are other points of high kappa that are not of
  * interest one should be careful that the region of interest is not missed.
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

  ImageInfo *critexport;
  //unsigned long j,k,m,jold;
  unsigned long Npoints,i=0,j;
  short refinements;
  //short spur,closed;
  double maxgridsize,mingridsize,x[2];
  ListHndl negpointlist = NewList();
  Point *minpoint = NewPoint(x,0);


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
	  critcurve->Npoints=Npoints;

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
		  critcurve->points=NewPointArray(1,false);
		  critcurve->points->in_image = FALSE;
		  critcurve->Npoints=1;
		  PointCopyData(critcurve->points,minpoint);
	  }else{

		  critcurve->points=NewPointArray(Npoints,false);

		  MoveToTopList(negpointlist);
		  for(i=0;i<negpointlist->Npoints;++i){
			  PointCopyData(&(critcurve->points[i]),negpointlist->current);
 			  MoveDownList(negpointlist);
		  }
	  }

	  EmptyList(negpointlist);

	  if(verbose) std::printf("find_crit, going into findborders 1\n");
	  findborders2(grid->i_tree,critcurve);
	  if(verbose) std::printf("find_crit, came out of findborders 1\n");


	  /* make inner border the image */
	  MoveToTopKist(critcurve->innerborder);
	  for(i=0,maxgridsize=0.0,mingridsize=1.0e99;i<critcurve->innerborder->Nunits();++i){
		  PointCopyData(&(critcurve->points[i]),getCurrentKist(critcurve->innerborder));
		  if(critcurve->points[i].gridsize > maxgridsize) maxgridsize=critcurve->points[i].gridsize;
		  if(critcurve->points[i].gridsize < mingridsize) mingridsize=critcurve->points[i].gridsize;

		  MoveDownKist(critcurve->innerborder);
	  }
	  critcurve->Npoints=critcurve->innerborder->Nunits();

	  // find borders again to properly define outer border
	  //std::printf("going into findborders 2\n");
	  findborders2(grid->i_tree,critcurve);
	  //std::printf("came out of findborders 2\n");

	  if(verbose) std::printf("find_crit, going into refine_grid\n");
     //std::printf("  Npoints=%i\n",critcurve->Npoints);
	  refinements=refine_grid(lens,grid,critcurve,1,resolution,2,false);
	  if(verbose) std::printf("find_crit, came out of refine_grid\n");

	  if(verbose) cout << "Npoints " << critcurve->Npoints << endl;

	  if(refinements==0){
		  break;
	  }else free(critcurve->points);
  }

  if(verbose) std::printf("find_crit, number of caustic points: %li\n",critcurve->Npoints);

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
			       ,*Ncrits,critcurve->points[0].gridsize); exit(1);}

  /* find area of critical curves */
  x[0]=x[1]=0.0;
  for(i=0;i<*Ncrits;++i){
	  if(critcurve[i].Npoints > 5){
			   windings(x,critcurve[i].points,critcurve[i].Npoints,&(critcurve[i].area),0);
			   //std::printf("critarea = %e\n",critcurve[i].area);
	  }else critcurve[i].area=0;
  }

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

/** \ingroup ImageFinding
  *
  * \brief Finds critical curves and caustics.
  *
  * OUTPUT: each critical curve is in a array of IamgeInfo's
  *         result.parity = 1 tangential caustic, 2 radial, 0 not enough points to determine
  *  the inner out outer boundaries of the result are the estimated critical curves
  *
  * The critical curve is found by refining the edges of regions of negative magnification.
  * If there are no regions of negative magnification in the original grid the grid is refined
  * around the point of highest kappa.  If there are other points of high kappa that are not of
  * interest one should be careful that the region of interest is not missed.
  */

void find_crit_kist(
	    LensHndl lens             /// The lens model.
		,GridHndl grid            /// The grid.  It must be initialized.
		,ImageInfo *critcurve     /// Structure to hold critical curve.  Must be pre-allocated with maxNcrit elements. Stored in critcurve[i].imagekist.
		,int maxNcrits            /// Maximum number of critical curves.
		,int *Ncrits              /// The number of critical curves found.
		,double resolution        /// The target resolution that the critical curve is mapped on the image plane.
		,bool *orderingsuccess    /// true if ordering was successful.
		,bool ordercurve          /// Order the curve so that it can be drawn or used to find the winding number.
		,bool verbose
		){

  //ImageInfo *critexport;
  //unsigned long j,k,m,jold;
  unsigned long i=0;
  short refinements;
  //short spur,closed;
  double maxgridsize,mingridsize,x[2];
  ImageInfo negimage;
  Point *minpoint = NewPoint(x,0);
  Kist newpoint_kist;

  minpoint->invmag=1.0e99;
  /*point=NmewPoint(x,0);*/

  // find list of points with negative magnification
  negimage.imagekist->Empty();
  MoveToTopList(grid->i_tree->pointlist);
  minpoint = grid->i_tree->pointlist->current;
  for(i=0;i<grid->i_tree->pointlist->Npoints;++i){
	  if(grid->i_tree->pointlist->current->invmag < 0){

		  //InsertAfterCurrent(negpointlist,grid->i_tree->pointlist->current->x,grid->i_tree->pointlist->current->id
		  //		  ,grid->i_tree->pointlist->current->image);
		  //MoveDownList(negpointlist);
		  //PointCopyData(negpointlist->current,grid->i_tree->pointlist->current);
		  negimage.imagekist->InsertAfterCurrent(grid->i_tree->pointlist->current);
		  negimage.imagekist->Down();
	  }

	  // record point of maximum kappa
	  if(grid->i_tree->pointlist->current->kappa > minpoint->kappa) minpoint = grid->i_tree->pointlist->current;
	  MoveDownList(grid->i_tree->pointlist);
  }

  if(negimage.imagekist ->Nunits() == 0){
	  if(minpoint->gridsize <= resolution){  // no caustic found at this resolution
		  *Ncrits=0;
		  *orderingsuccess = false;
		  return;
	  }

	  /* if there is no negative magnification points use maximum mag point */
	  negimage.imagekist->InsertAfterCurrent(minpoint);
  }

  for(;;){

	  //EmptyList(negpointlist);
	  //negpointkist->Empty();

	  negimage.imagekist->MoveToTop();
	  do{negimage.imagekist ->getCurrent()->in_image = TRUE;} while(negimage.imagekist->Down());

	  if(verbose) std::printf("find_crit, going into findborders 1\n");
	  //findborders2(grid->i_tree,critcurve);
	  findborders4(grid->i_tree,&negimage);
	  if(verbose) std::printf("find_crit, came out of findborders 1\n");

	  // unmark image points
	  negimage.imagekist->MoveToTop();
	  do{negimage.imagekist->getCurrent()->in_image = FALSE;} while(negimage.imagekist->Down());

	  // make inner border of the image
	  critcurve->imagekist->Empty();
	  MoveToTopKist(negimage.innerborder);
	  for(i=0,maxgridsize=0.0,mingridsize=1.0e99;i<negimage.innerborder->Nunits();++i){
		  //PointCopyData(&(critcurve->points[i]),getCurrentKist(negpoints.innerborder));
		  //if(critcurve->points[i].gridsize > maxgridsize) maxgridsize=critcurve->points[i].gridsize;
		  //if(critcurve->points[i].gridsize < mingridsize) mingridsize=critcurve->points[i].gridsize;

		  if(getCurrentKist(negimage.innerborder)->gridsize > maxgridsize) maxgridsize = getCurrentKist(negimage.innerborder)->gridsize;
		  if(getCurrentKist(negimage.innerborder)->gridsize < mingridsize) mingridsize = getCurrentKist(negimage.innerborder)->gridsize;

		  critcurve->imagekist->InsertAfterCurrent(getCurrentKist(negimage.innerborder));
		  critcurve->imagekist->Down();
		  critcurve->imagekist->getCurrent()->in_image = TRUE;

		  MoveDownKist(negimage.innerborder);
	  }
	  findborders4(grid->i_tree,critcurve);
	  //std::printf("came out of findborders 2\n");

	  if(verbose) std::printf("find_crit, going into refine_grid\n");
     //std::printf("  Npoints=%i\n",critcurve->Npoints);
	  //refinements=refine_grid(lens,grid->i_tree,grid->s_tree,critcurve,1,resolution,2,false);
	  refinements=refine_grid_kist(lens,grid,critcurve,1,resolution,2,false,&newpoint_kist);
	  if(verbose) std::printf("find_crit, came out of refine_grid\n");

	  if(verbose) cout << "Npoints " << critcurve->imagekist->Nunits() << endl;

	  /* to limit the size of the crit curve *
	  if(critcurve->imagekist->Nunits() > 1000)
		  break;
	  // */

	  if(refinements==0) break;
	  //}else free(critcurve->points);

	  // add new negative points to negpoints
	  newpoint_kist.MoveToTop();
	  negimage.imagekist->MoveToBottom();
	  do{
		  if(newpoint_kist.getCurrent()->invmag < 0)
			  negimage.imagekist->InsertAfterCurrent(newpoint_kist.getCurrent());
	  }while(newpoint_kist.Down());
  }
  newpoint_kist.Empty();
  critcurve->imagekist->MoveToTop();
  do{critcurve->imagekist->getCurrent()->in_image = FALSE;}while(critcurve->imagekist->Down());

  if(verbose) std::printf("find_crit, number of caustic points: %li\n",critcurve->imagekist->Nunits());

  *orderingsuccess = true;
  if(ordercurve){

	  unsigned long Npoints = critcurve->imagekist->Nunits();
	  unsigned long NewNumber;
	  divide_images_kist(grid->i_tree,critcurve,Ncrits,maxNcrits);

	  if(*Ncrits > maxNcrits){ERROR_MESSAGE(); std::printf("ERROR: in find_crit, too many critical curves Ncrits=%i > maxNcrits gridsize=%e\n"
				       ,*Ncrits,critcurve->imagekist->getCurrent()->gridsize); exit(1);}

	  // order points in curve

	  x[0]=x[1]=0.0;
	  Point *tmp_points = NewPointArray(Npoints,false);
	  unsigned long ii;
	  for(i=0;i<*Ncrits;++i){
		  // return in_image value to false
		  critcurve[i].imagekist->MoveToTop();
		  do{critcurve[i].imagekist->getCurrent()->in_image = FALSE;}while(critcurve[i].imagekist->Down());

		  critcurve[i].imagekist->MoveToTop();
		  //copy points into a point array for compatibility with curve ordering routines
		  for(ii=0; ii < critcurve[i].imagekist->Nunits() ; ++ii, critcurve[i].imagekist->Down())
			  PointCopyData(&tmp_points[ii],getCurrentKist(critcurve[i].imagekist));

		  // order the curve
		  NewNumber = order_curve4(tmp_points,critcurve[i].imagekist->Nunits());
		  if(NewNumber != critcurve[i].imagekist->Nunits() ) *orderingsuccess = false;

		  // find area of critical curves
		  if(critcurve[i].imagekist->Nunits() > 5){
				   windings(x,tmp_points,critcurve[i].imagekist->Nunits(),&(critcurve[i].area),0);
		  }else critcurve[i].area=0;

		  // resort points in imagekist
		  for(ii=0;ii<NewNumber;++ii){
			  critcurve[i].imagekist->MoveToTop();
			  do{
				  if(tmp_points[ii].id == critcurve[i].imagekist->getCurrent()->id)
					  critcurve[i].imagekist->MoveCurrentToBottom();
			  }while(critcurve[i].imagekist->Down());
		  }
		  // remove points that were not linked up into a closed curve
		  critcurve[i].imagekist->MoveToTop();
		  critcurve[i].imagekist->JumpDown(NewNumber);
		  while(critcurve[i].imagekist->Nunits() > NewNumber){
			  critcurve[i].imagekist->TakeOutCurrent();
			  critcurve[i].imagekist->Down();
		  }

	  }
	  FreePointArray(tmp_points,false);
	  //split_order_curve4(critcurve,maxNcrits,Ncrits);
  }else if(critcurve->imagekist->Nunits() > 0) *Ncrits=1;
  if(critcurve->imagekist->Nunits() == 0) *Ncrits=0;

  return ;
}

