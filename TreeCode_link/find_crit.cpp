/*
 * find_crit.c
 *
 *  Created on: Sep 8, 2009
 *      Author: R.B. Metcalf
 */

#include "slsimlib.h"

#define NMAXCRITS 1000

using namespace std;
/** \ingroup ImageFinding
  *
  * \brief Finds critical curves and caustics.
  *
  * OUTPUT: each critical curve is in a array of ImageInfo's
  *         result.parity = 1 tangential caustic, 2 radial, 0 not enough points to determine
  *  the inner out outer boundaries of the result are the estimated critical curves
  *
  * The critical curve is found by refining the edges of regions of negative magnification.
  * If there are no regions of negative magnification in the original grid the grid is refined
  * around the point of highest kappa.  If there are other points of high kappa that are not of
  * interest one should be careful that the region of interest is not missed.
  */

void ImageFinding::find_crit(
	    LensHndl lens             /// The lens model.
		,GridHndl grid            /// The grid.  It must be initialized.
		,ImageInfo *critcurve     /// Structure to hold critical curve.  Must be pre-allocated with maxNcrit elements. Stored in critcurve[i].imagekist.
		,int maxNcrits            /// Maximum number of critical curves.
		,int *Ncrits              /// The number of critical curves found.
		,PosType resolution        /// The target resolution that the critical curve is mapped on the image plane.
		,bool *orderingsuccess    /// true if ordering was successful.
		,bool ordercurve          /// Order the curve so that it can be drawn or used to find the winding number.
    ,bool dividecurves        /// Divide the critical curves into separate curves by whether they are attached
    ,PosType invmag_min        /// finds regions with 1/magnification < invmag_min, set to zero for caustics               
		,bool verbose
		){

  long i=0;
  short refinements;
  //short spur,closed;
  PosType maxgridsize,mingridsize,x[2];
  ImageInfo negimage;
  Kist<Point> newpoint_kist;

  ImageInfo *pseudocurve = new ImageInfo[maxNcrits];
  bool pseuodcaustic = false;
  PosType pseudolimit = -0.01;

  // find kist of points with negative magnification
  negimage.imagekist->Empty();
  MoveToTopList(grid->i_tree->pointlist);
  Point *minpoint = grid->i_tree->pointlist->current;
    
  for(i=0;i<grid->i_tree->pointlist->Npoints;++i){
	  if(grid->i_tree->pointlist->current->invmag < invmag_min){
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

	  // if there is no negative magnification points use maximum mag point
	  negimage.imagekist->InsertAfterCurrent(minpoint);
  }

  for(;;){

	  negimage.imagekist->MoveToTop();
	  do{negimage.imagekist ->getCurrent()->in_image = YES;} while(negimage.imagekist->Down());

	  if(verbose) std::printf("find_crit, going into findborders 1\n");
	  //findborders2(grid->i_tree,critcurve);
	  findborders4(grid->i_tree,&negimage);
	  if(verbose) std::printf("find_crit, came out of findborders 1\n");

	  // unmark image points
	  negimage.imagekist->MoveToTop();
	  do{negimage.imagekist->getCurrent()->in_image = NO;} while(negimage.imagekist->Down());

	  // make inner border of the image
	  critcurve->imagekist->Empty();
	  negimage.innerborder->MoveToTop();
	  for(i=0,maxgridsize=0.0,mingridsize=1.0e99;i<negimage.innerborder->Nunits();++i){

		  if(negimage.innerborder->getCurrent()->gridsize > maxgridsize) maxgridsize = negimage.innerborder->getCurrent()->gridsize;
		  if(negimage.innerborder->getCurrent()->gridsize < mingridsize) mingridsize = negimage.innerborder->getCurrent()->gridsize;

		  critcurve->imagekist->InsertAfterCurrent(negimage.innerborder->getCurrent());
		  critcurve->imagekist->Down();
		  critcurve->imagekist->getCurrent()->in_image = YES;

		  negimage.innerborder->Down();
	  }
	  findborders4(grid->i_tree,critcurve);
	  //std::printf("came out of findborders 2\n");

	  if(verbose) std::printf("find_crit, going into refine_grid\n");
     //std::printf("  Npoints=%i\n",critcurve->Npoints);
	  //refinements=refine_grid(lens,grid->i_tree,grid->s_tree,critcurve,1,resolution,2,false);
	  refinements=ImageFinding::refine_grid_kist(lens,grid,critcurve,1,resolution,2,&newpoint_kist,true);
	  if(verbose) std::printf("find_crit, came out of refine_grid\n");

	  if(verbose) cout << "Npoints " << critcurve->imagekist->Nunits() << endl;

	  if(refinements==0) break;
	  //}else free(critcurve->points);

	  // add new negative points to negpoints
	  newpoint_kist.MoveToTop();
	  negimage.imagekist->MoveToBottom();
	  do{
		  if(newpoint_kist.getCurrent()->invmag < invmag_min)
			  negimage.imagekist->InsertAfterCurrent(newpoint_kist.getCurrent());
	  }while(newpoint_kist.Down());
  }

  newpoint_kist.Empty();
  if(critcurve->imagekist->Nunits()){
	  critcurve->imagekist->MoveToTop();
	  do{critcurve->imagekist->getCurrent()->in_image = NO;}while(critcurve->imagekist->Down());
  }

  if(verbose) std::printf("find_crit, number of caustic points: %li\n",critcurve->imagekist->Nunits());

  unsigned long Npoints = critcurve->imagekist->Nunits();

  if(dividecurves) divide_images_kist(grid->i_tree,critcurve,Ncrits,maxNcrits);
  else *Ncrits = 1;

  if(*Ncrits > maxNcrits){ERROR_MESSAGE(); std::printf("ERROR: in find_crit, too many critical curves Ncrits=%i > maxNcrits gridsize=%e\n"
			       ,*Ncrits,critcurve->imagekist->getCurrent()->gridsize); exit(1);}


  if(pseuodcaustic){
  	  // Find points within each critical curve that have invmag < pseudolimit
  	  // If there are none use the minimum invmag value point.
  	  Point *minmupoint;
  	  PosType mumin = 0.0;
  	  Kist<Point> *newpoints = new Kist<Point>;

  	  pseudocurve->imagekist->copy(negimage.imagekist);
  	  divide_images_kist(grid->i_tree,pseudocurve,Ncrits,maxNcrits);

  	  for(int i =1;i<*Ncrits;++i){
  		  mumin = 0.0;
  		  pseudocurve[i].imagekist->MoveToBottom();
  		  do{
  			  if(pseudocurve[i].imagekist->getCurrent()->invmag < mumin){
  				  minmupoint = pseudocurve[i].imagekist->getCurrent();
  				  mumin = pseudocurve[i].imagekist->getCurrent()->invmag;
  			  }
  			  if(pseudocurve[i].imagekist->getCurrent()->invmag > pseudolimit){
  				  pseudocurve[i].imagekist->getCurrent()->in_image = NO;
  				  pseudocurve[i].imagekist->TakeOutCurrent();
  			  }
  		  }while(pseudocurve[i].imagekist->Up());

  		  // in case one before top was taken out
  		  if(pseudocurve[i].imagekist->getCurrent()->invmag > pseudolimit){
  			  pseudocurve[i].imagekist->getCurrent()->in_image = NO;
  			  pseudocurve[i].imagekist->TakeOutCurrent();
  		  }

  		  if(pseudocurve[i].imagekist->Nunits() == 0){
  			  pseudocurve[i].imagekist->InsertAfterCurrent(minmupoint);
  			  pseudocurve[i].ShouldNotRefine = 0;  // marks that region has not been found
  		  }else{
  			  pseudocurve[i].ShouldNotRefine = 1;
  		  }


  		  while( pseudocurve[i].imagekist->Nunits() < 100 &&
  				  refine_edges(lens,grid,&pseudocurve[i],1,0.01*resolution/sqrt(fabs(pseudolimit)),1,newpoints)
  				  ){
  			  // update region
  			  if(pseudocurve[i].ShouldNotRefine == 0){
  				  mumin = pseudocurve[i].imagekist->getCurrent()->invmag;
  				  minmupoint = pseudocurve[i].imagekist->getCurrent();
  				  pseudocurve[i].imagekist->getCurrent()->in_image = NO;
  				  pseudocurve[i].imagekist->TakeOutCurrent();
  			  }

  			  newpoints->MoveToTop();
  			  do{
  				  if(newpoints->getCurrent()->invmag < mumin){
  					  mumin = newpoints->getCurrent()->invmag;
  					  minmupoint = newpoints->getCurrent();
  				  }
  				  if(newpoints->getCurrent()->invmag < 0.0)
  					  negimage.imagekist->InsertAfterCurrent(newpoints->getCurrent());
  				  if(newpoints->getCurrent()->invmag < pseudolimit){
  					  newpoints->getCurrent()->in_image = YES;
  					  pseudocurve[i].imagekist->InsertAfterCurrent(newpoints->getCurrent());
  				  }
  			  }while(newpoints->Down());

  			  if(pseudocurve[i].ShouldNotRefine == 0){

  				  if(pseudocurve[i].imagekist->Nunits() == 0){
  					  minmupoint->in_image = YES;
  					  pseudocurve[i].imagekist->InsertAfterCurrent(minmupoint);
  				  }else{
  					  pseudocurve[i].ShouldNotRefine = 1;
  				  }
  			  }
  			  findborders4(grid->i_tree,&pseudocurve[i]);
  		  }

  		  pseudocurve[i].imagekist->MoveToTop();
  		  do{
  			pseudocurve[i].imagekist->getCurrent()->in_image = NO;
  		  }while(pseudocurve[i].imagekist->Down());
  	  }
  }

  *orderingsuccess = true;
  if(ordercurve){

	  unsigned long NewNumber;
	  unsigned long ii;

	  // order points in curve

	  x[0]=x[1]=0.0;
	  Point *tmp_points = NewPointArray(Npoints);
	  for(i=0;i<*Ncrits;++i){
		  // return in_image value to false
		  critcurve[i].imagekist->MoveToTop();
		  do{critcurve[i].imagekist->getCurrent()->in_image = NO;}while(critcurve[i].imagekist->Down());

		  critcurve[i].imagekist->MoveToTop();
		  //copy points into a point array for compatibility with curve ordering routines
		  for(ii=0; ii < critcurve[i].imagekist->Nunits() ; ++ii, critcurve[i].imagekist->Down())
			  PointCopyData(&tmp_points[ii],critcurve[i].imagekist->getCurrent());

		  // order the curve
		  NewNumber = Utilities::order_curve4(tmp_points,critcurve[i].imagekist->Nunits());
		  if(NewNumber != critcurve[i].imagekist->Nunits() ){
        *orderingsuccess = false;
        std::cout << "find_crit(), i = " << i << " ordering curve NewNumber = " << NewNumber << " oldnumber = " << critcurve[i].imagekist->Nunits() << std::endl;
      }
		  // find area of critical curves
		  if(critcurve[i].imagekist->Nunits() > 5){
				   Utilities::windings(x,tmp_points,NewNumber,&(critcurve[i].area),0);
		  }else critcurve[i].area=0;

		  // re-sort points in imagekist
		  for(ii=0;ii<NewNumber;++ii){
			  critcurve[i].imagekist->MoveToTop();
			  do{
				  if(tmp_points[ii].id == critcurve[i].imagekist->getCurrent()->id){
					  critcurve[i].imagekist->MoveCurrentToBottom();
            break;
          }
			  }while(critcurve[i].imagekist->Down());
		  }
      
		  // remove points that were not linked up into a closed curve
		  critcurve[i].imagekist->MoveToTop();
		  //critcurve[i].imagekist->JumpDown(NewNumber);
		  while(critcurve[i].imagekist->Nunits() > NewNumber){
			  critcurve[i].imagekist->TakeOutCurrent();
		  }

	  }
	  FreePointArray(tmp_points,false);
	  //split_order_curve4(critcurve,maxNcrits,Ncrits);
  }//else if(critcurve->imagekist->Nunits() > 0) *Ncrits=1;
  if(critcurve->imagekist->Nunits() == 0) *Ncrits=0;

  for(i=0;i<*Ncrits;++i){
	  critcurve[i].centroid[0] = 0;
	  critcurve[i].centroid[1] = 0;
	  if(critcurve[i].imagekist->Nunits() > 1){
		  critcurve[i].imagekist->MoveToTop();
		  do{
			  critcurve[i].centroid[0] += critcurve[i].imagekist->getCurrent()->x[0];
			  critcurve[i].centroid[1] += critcurve[i].imagekist->getCurrent()->x[1];
		  }while(critcurve[i].imagekist->Down());
		  critcurve[i].centroid[0] /= critcurve[i].imagekist->Nunits();
		  critcurve[i].centroid[1] /= critcurve[i].imagekist->Nunits();
	  }else{
		  // take out curves with no points
		  critcurve[i].copy(critcurve[*Ncrits-1]);
		  *Ncrits -= 1;
		  --i;
	  }
  }

  /*/******* TODO: test line *****************
  char chrstr[100];
  std::string output = "test_crit";
  PosType tmp_range=0,tmp,xx[2];
  for(i=0; i< MIN(*Ncrits,50) ; ++i){
    
    if(critcurve[i].getNimagePoints() > 10){
      tmp_range=0;
      for(critcurve[i].imagekist->MoveToTop();!(critcurve[i].imagekist->OffBottom()); critcurve[i].imagekist->Down()){
        tmp = pow(critcurve[i].centroid[0] - critcurve[i].imagekist->getCurrent()->x[0],2)
            + pow(critcurve[i].centroid[1] - critcurve[i].imagekist->getCurrent()->x[1],2);
        if(tmp > tmp_range) tmp_range = tmp;
      }
    
      tmp_range = sqrt(tmp_range)*2;
    
      snprintf(chrstr,100,"%i",i);
      PixelMap map(critcurve[i].centroid,500,tmp_range/500);
    
      //map.AddImages(&critcurve[i],1,0);
      map.AddCurve(&critcurve[i],1);
     
      map.printFITS(output + chrstr + ".0.fits");
      
      critcurve[i].imagekist->TranformPlanes();

      tmp_range=0;
      xx[0] = critcurve[i].imagekist->getCurrent()->x[0];
      xx[1] = critcurve[i].imagekist->getCurrent()->x[1];
      for(critcurve[i].imagekist->MoveToTop();!(critcurve[i].imagekist->OffBottom()); critcurve[i].imagekist->Down()){
        tmp = pow(xx[0] - critcurve[i].imagekist->getCurrent()->x[0],2)
        + pow(xx[1] - critcurve[i].imagekist->getCurrent()->x[1],2);
        if(tmp > tmp_range) tmp_range = tmp;
      }
      
      tmp_range = sqrt(tmp_range)*2;

      PixelMap map2(xx,500,tmp_range/500);
     
      map2.AddCurve(&critcurve[i],1);
    
      map2.printFITS(output + chrstr + ".1.fits");
      critcurve[i].imagekist->TranformPlanes();

    }
  }
  //************************************/

  return ;
}
void ImageFinding::find_crit2(
	    LensHndl lens             /// The lens model.
		,GridHndl grid            /// The grid.  It must be initialized.
		,ImageInfo *critcurve     /// Structure to hold critical curve.  Must be pre-allocated with maxNcrit elements. Stored in critcurve[i].imagekist.
		,int maxNcrits            /// Maximum number of critical curves.
		,int *Ncrits              /// The number of critical curves found.
		,PosType resolution        /// The target resolution that the critical curve is mapped on the image plane.
		,bool *orderingsuccess    /// true if ordering was successful.
		,bool ordercurve          /// Order the curve so that it can be drawn or used to find the winding number.
    ,bool dividecurves        /// Divide the critical curves into seporate curves by whether they are attached
    ,PosType invmag_min        /// finds regions with 1/magnification < invmag_min
		,bool verbose
		){

  long i=0;
  long refinements;
  //short spur,closed;
  PosType maxgridsize,mingridsize,x[2];
  Kist<Point> newpoint_kist,neighborkist;

  // find kist of points with negative magnification
  critcurve->imagekist->Empty();
  MoveToTopList(grid->i_tree->pointlist);
  Point *minpoint = grid->i_tree->pointlist->current;

  for(i=0;i<grid->i_tree->pointlist->Npoints;++i){
	  if(grid->i_tree->pointlist->current->invmag < invmag_min){
		  critcurve->imagekist->InsertAfterCurrent(grid->i_tree->pointlist->current);
		  critcurve->imagekist->Down();
	  }

	  // record point of maximum kappa
	  if(grid->i_tree->pointlist->current->kappa > minpoint->kappa) minpoint = grid->i_tree->pointlist->current;
	  MoveDownList(grid->i_tree->pointlist);
  }
  bool maxpoint = false;

  if(critcurve->imagekist->Nunits() == 0){
	  if(minpoint->gridsize <= resolution){  // no caustic found at this resolution
		  *Ncrits=0;
		  *orderingsuccess = false;
		  return;
	  }

	  // if there is no negative magnification points use maximum mag point
	  critcurve->imagekist->InsertAfterCurrent(minpoint);
	  maxpoint =true;
  }

  critcurve->imagekist->SetInImage(YES);

  findborders4(grid->i_tree,critcurve);

  for(int k=0;;++k){

	  refinements=refine_edges(lens,grid,critcurve,1,resolution/2,1,&newpoint_kist);
	  //refinements=refine_edges(lens,grid,critcurve,1,1.0e-3,0,&newpoint_kist);
	  //refinements=refine_grid_kist(lens,grid,critcurve,1,resolution,2,&newpoint_kist);

	  if(!refinements) break;
	  critcurve->outerborder->SetInImage(MAYBE);

	  // add new points to negative region
	  newpoint_kist.MoveToTop();
	  critcurve->imagekist->MoveToBottom();
	  do{
		  if(newpoint_kist.getCurrent()->invmag < invmag_min){
			  newpoint_kist.getCurrent()->in_image = YES;
			  critcurve->imagekist->InsertAfterCurrent(newpoint_kist.getCurrent());

			  // It is possible that gridrange[] will not be maintained
			  if(critcurve->gridrange[1] < newpoint_kist.getCurrent()->gridsize)
				  critcurve->gridrange[1] = newpoint_kist.getCurrent()->gridsize;

			  if(critcurve->gridrange[2] > newpoint_kist.getCurrent()->gridsize)
				  critcurve->gridrange[2] = newpoint_kist.getCurrent()->gridsize;

		  }else{
			  newpoint_kist.getCurrent()->in_image = NO;
		  }
	  }while(newpoint_kist.Down());

	  if(maxpoint){
		  if(critcurve->imagekist->Nunits() > 1){
			  // take out old max point
			  critcurve->imagekist->MoveToTop();
			  do{
				  if(critcurve->imagekist->getCurrent()->invmag > 0){
					  critcurve->imagekist->getCurrent()->in_image = NO;
					  critcurve->imagekist->TakeOutCurrent();
					  break;
				  }
			  }while(critcurve->imagekist->Down());
			  maxpoint = false;
		  }else{
			  // update maximum kappa point if no negative magnification points have been found
			  newpoint_kist.MoveToTop();
			  do{
				  if(newpoint_kist.getCurrent()->kappa
						  > critcurve->imagekist->getCurrent()->kappa ){
					  critcurve->imagekist->getCurrent()->in_image = NO;
					  critcurve->imagekist->TakeOutCurrent();
					  newpoint_kist.getCurrent()->in_image = YES;
					  critcurve->imagekist->InsertAfterCurrent(newpoint_kist.getCurrent());
				  }
			  }while(newpoint_kist.Down());
		  }
	  }

	  // check which points in inner border are still in the border
	  bool ininner;
	  unsigned long Ntmp = critcurve->innerborder->Nunits();
	  critcurve->innerborder->MoveToTop();
	  for(unsigned long j=0;j < Ntmp;++j){

		  ininner=false;

		  grid->i_tree->FindAllBoxNeighborsKist(critcurve->innerborder->getCurrent(),&neighborkist);

		  neighborkist.MoveToTop();
		  do{

			  if( neighborkist.getCurrent()->in_image != YES){  // point is a neighbor
				  ininner=true;

				  if(neighborkist.getCurrent()->in_image == NO){  // if point is not yet in outerborder
					  // add point to outerborder
					  neighborkist.getCurrent()->in_image = MAYBE;
					  critcurve->outerborder->InsertAfterCurrent(neighborkist.getCurrent());
					  critcurve->outerborder->Down();
				  }
			  }

		  }while(neighborkist.Down());

		  if(!ininner){
			  bool tmp = critcurve->innerborder->AtTop();
			  critcurve->innerborder->TakeOutCurrent();
			  if(!tmp) critcurve->innerborder->Down();
		  }else{
			  critcurve->innerborder->Down();
		  }
	  }

	  // Take out outer border points that are no longer in outer border
	  critcurve->gridrange[0] = 0.0;
	  Ntmp = critcurve->outerborder->Nunits();
	  critcurve->outerborder->MoveToTop();
	  bool tmpbool;
	  for(unsigned long j=0;j<Ntmp;++j){

		  tmpbool = true;
		  assert(critcurve->outerborder->getCurrent()->in_image == MAYBE);
		  grid->i_tree->FindAllBoxNeighborsKist(critcurve->outerborder->getCurrent(),&neighborkist);
		  neighborkist.MoveToTop();
		  do{
			  if(neighborkist.getCurrent()->in_image == YES){
				  //critcurve->outerborder->getCurrent()->in_image = NO;
				  if(critcurve->outerborder->getCurrent()->gridsize
						  > critcurve->gridrange[0]) critcurve->gridrange[0]
						           = critcurve->outerborder->getCurrent()->gridsize;
				  tmpbool = false;
				  break;
			  }
		  }while(neighborkist.Down());

		  if(tmpbool){  // no neighbor in image was found
			  bool tmp = critcurve->outerborder->AtTop();
			  critcurve->outerborder->getCurrent()->in_image = NO;
			  critcurve->outerborder->TakeOutCurrent();
			  if(!tmp) critcurve->outerborder->Down();
		  }else{
			  critcurve->outerborder->Down();
		  }
	  }

	  //critcurve->outerborder->SetInImage(MAYBE);
	  // sort new points into inner and outer borders
	  newpoint_kist.MoveToTop();
	  do{

		  if(newpoint_kist.getCurrent()->in_image != MAYBE){
			  grid->i_tree->FindAllBoxNeighborsKist(newpoint_kist.getCurrent(),&neighborkist);

			  tmpbool = true;
			  neighborkist.MoveToTop();
			  do{
				  if( newpoint_kist.getCurrent()->in_image == YES){
					  if(neighborkist.getCurrent()->in_image != YES){
						  if(tmpbool){
							  critcurve->innerborder->InsertAfterCurrent(newpoint_kist.getCurrent());
							  tmpbool = false;
						  }
						  if(neighborkist.getCurrent()->in_image == NO){
							  neighborkist.getCurrent()->in_image = MAYBE;
							  critcurve->outerborder->InsertAfterCurrent(neighborkist.getCurrent());
						  }
					  }
				  }else{
					  if(neighborkist.getCurrent()->in_image == YES){
						  newpoint_kist.getCurrent()->in_image = MAYBE;
						  critcurve->outerborder->InsertAfterCurrent(newpoint_kist.getCurrent());
						  break;
					  }
				  }
			  }while(neighborkist.Down());
		  }
	  }while(newpoint_kist.Down());

	  critcurve->outerborder->SetInImage(NO);
  }

  if(maxpoint){
 	  *Ncrits = 0;
 	  assert(critcurve->imagekist->Nunits() == 1);
 	  critcurve->imagekist->getCurrent()->in_image = NO;
 	  critcurve->imagekist->Empty();
 	  critcurve->outerborder->Empty();
 	  critcurve->innerborder->Empty();
 	  return;
   }

  // make inner border the image
  critcurve->imagekist->SetInImage(NO);
  critcurve->imagekist->Empty();
  critcurve->imagekist->copy(critcurve->innerborder);

  unsigned long Npoints = critcurve->imagekist->Nunits();
  if(dividecurves) divide_images_kist(grid->i_tree,critcurve,Ncrits,maxNcrits);
  else *Ncrits = 1;

  if(*Ncrits > maxNcrits){ERROR_MESSAGE(); std::printf("ERROR: in find_crit, too many critical curves Ncrits=%i > maxNcrits gridsize=%e\n"
				       ,*Ncrits,critcurve->imagekist->getCurrent()->gridsize); exit(1);}


  for(i=0;i<*Ncrits;++i) critcurve[i].imagekist->SetInImage(NO);


  *orderingsuccess = true;
  if(ordercurve && dividecurves){

	  //unsigned long NewNumber;

	  // order points in curve

	  x[0]=x[1]=0.0;
	  //Point *tmp_points = NewPointArray(Npoints);
	  //unsigned long ii;
	  for(i=0;i<*Ncrits;++i){

		  //critcurve[i].imagekist->MoveToTop();
		  //copy points into a point array for compatibility with curve ordering routines
		  //for(ii=0; ii < critcurve[i].imagekist->Nunits() ; ++ii, critcurve[i].imagekist->Down())
			//  PointCopyData(&tmp_points[ii],getCurrentKist(critcurve[i].imagekist));

		  // order the curve
//		  NewNumber = Utilities::order_curve4(tmp_points,critcurve[i].imagekist->Nunits());
//		  NewNumber = Utilities::order_curve4(critcurve[i].imagekist);
//		  NewNumber = Utilities::order_curve5(critcurve[i].imagekist);
      
      //Utilities::ordered_convexhull(critcurve[i].imagekist);
      Utilities::ordered_shrink_wrap(critcurve[i].imagekist);
      
/*		  if(NewNumber != critcurve[i].imagekist->Nunits() ){
        std::cout << "find_crit2(), i = " << i << " ordering curve NewNumber = " << NewNumber << " oldnumber = " << critcurve[i].imagekist->Nunits() << std::endl;
        *orderingsuccess = false;
      }
*/

		  // resort points in imagekist
		  /*for(ii=0;ii<NewNumber;++ii){
			  critcurve[i].imagekist->MoveToTop();
			  do{
				  if(tmp_points[ii].id == critcurve[i].imagekist->getCurrent()->id)
					  critcurve[i].imagekist->MoveCurrentToBottom();
			  }while(critcurve[i].imagekist->Down());
		  }*/
      
		  // remove points that were not linked up into a closed curve
      //std::cout << "This needs to be uncommented" << std::endl;
		  //critcurve[i].imagekist->MoveToBottom();
		  //while(critcurve[i].imagekist->Nunits() > NewNumber) critcurve[i].imagekist->TakeOutCurrent();

      // find area of critical curves
		  if(critcurve[i].imagekist->Nunits() > 5){
        Utilities::windings(critcurve[i].centroid,critcurve[i].imagekist,&(critcurve[i].area),0);
		  }else critcurve[i].area=0;

	  }
	  //FreePointArray(tmp_points,false);
	  //split_order_curve4(critcurve,maxNcrits,Ncrits);
  }//else if(critcurve->imagekist->Nunits() > 0) *Ncrits=1;

  if(critcurve->imagekist->Nunits() == 0) *Ncrits=0;

  // Find centroid
  for(i=0;i<*Ncrits;++i){
    critcurve[i].centroid[0] = 0;
    critcurve[i].centroid[1] = 0;
    
    if(critcurve[i].imagekist->Nunits() >0){
      critcurve[i].imagekist->MoveToTop();
      do{
        critcurve[i].centroid[0] += critcurve[i].imagekist->getCurrent()->x[0];
        critcurve[i].centroid[1] += critcurve[i].imagekist->getCurrent()->x[1];
      }while(critcurve[i].imagekist->Down());
      critcurve[i].centroid[0] /= critcurve[i].imagekist->Nunits();
      critcurve[i].centroid[1] /= critcurve[i].imagekist->Nunits();
    }else{
      // take out curves with no points
      critcurve[i].copy(critcurve[*Ncrits-1]);
      *Ncrits -= 1;
      --i;
    }
  }
  

  /******* TODO: test line *****************
  char chrstr[100];
  std::string output = "test_crit";
  PosType tmp_range=0,tmp,xx[2];
  for(i=0; i< MIN(*Ncrits,50) ; ++i){
    
    if(critcurve[i].getNimagePoints() > 10){
      tmp_range=0;
      for(critcurve[i].imagekist->MoveToTop();!(critcurve[i].imagekist->OffBottom()); critcurve[i].imagekist->Down()){
        tmp = pow(critcurve[i].centroid[0] - critcurve[i].imagekist->getCurrent()->x[0],2)
        + pow(critcurve[i].centroid[1] - critcurve[i].imagekist->getCurrent()->x[1],2);
        if(tmp > tmp_range) tmp_range = tmp;
      }
      
      tmp_range = sqrt(tmp_range)*2;
      
      snprintf(chrstr,100,"%i",i);
      PixelMap map(critcurve[i].centroid,500,tmp_range/500);
      
      //map.AddImages(&critcurve[i],1,0);
      map.AddCurve(&critcurve[i],1);
      
        //Utilities::order_curve5(critcurve[i].imagekist);
      //map.AddCurve(&critcurve[i],1);
      
      map.printFITS(output + chrstr + ".0.fits");
      
      critcurve[i].imagekist->TranformPlanes();
      
      tmp_range=0;
      xx[0] = critcurve[i].imagekist->getCurrent()->x[0];
      xx[1] = critcurve[i].imagekist->getCurrent()->x[1];
      for(critcurve[i].imagekist->MoveToTop();!(critcurve[i].imagekist->OffBottom()); critcurve[i].imagekist->Down()){
        tmp = pow(xx[0] - critcurve[i].imagekist->getCurrent()->x[0],2)
        + pow(xx[1] - critcurve[i].imagekist->getCurrent()->x[1],2);
        if(tmp > tmp_range) tmp_range = tmp;
      }
      
      tmp_range = sqrt(tmp_range)*2;
      
      PixelMap map2(xx,500,tmp_range/500);
      
      map2.AddCurve(&critcurve[i],1);
      
      map2.printFITS(output + chrstr + ".1.fits");
      critcurve[i].imagekist->TranformPlanes();
      
    }
  }
  //************************************/
  
  return ;
}

/**
 *  This is a stripped down version of find_crit() for use in find_image_microlens() that
 *  refines the critical lines that are within the image.
 */
void refine_crit_in_image(
               LensHndl lens             /// The lens model.
               ,GridHndl grid            /// The grid.  It must be initialized.
               ,PosType r_source
               ,PosType x_source[]
               ,PosType resolution        /// The target resolution that the critical curve is mapped on the image plane.
               ){
    
    unsigned long i=0;
    short refinements;
    //short spur,closed;
    PosType maxgridsize,mingridsize,x[2];
    ImageInfo negimage,critcurve;
    Kist<Point> newpoint_kist;
        
    // find kist of points with negative magnification
    negimage.imagekist->Empty();
    MoveToTopList(grid->i_tree->pointlist);    
    for(i=0;i<grid->i_tree->pointlist->Npoints;++i){
        x[0] = grid->i_tree->pointlist->current->image->x[0] - x_source[0];
        x[1] = grid->i_tree->pointlist->current->image->x[1] - x_source[1];
        
        if(grid->i_tree->pointlist->current->invmag < 0 && r_source*r_source > (x[0]*x[0] + x[1]*x[1]) ){
            negimage.imagekist->InsertAfterCurrent(grid->i_tree->pointlist->current);
            negimage.imagekist->Down();
        }
        MoveDownList(grid->i_tree->pointlist);
    }
    
    if(negimage.imagekist->Nunits() == 0) return;
    
    for(;;){
        
        negimage.imagekist->MoveToTop();
        do{negimage.imagekist ->getCurrent()->in_image = YES;} while(negimage.imagekist->Down());
        
        findborders4(grid->i_tree,&negimage);
        
        // unmark image points
        negimage.imagekist->MoveToTop();
        do{negimage.imagekist->getCurrent()->in_image = NO;} while(negimage.imagekist->Down());
        
        // make inner border of the image
        critcurve.imagekist->Empty();
        negimage.innerborder->MoveToTop();
        for(i=0,maxgridsize=0.0,mingridsize=1.0e99;i<negimage.innerborder->Nunits();++i){
            
            if(negimage.innerborder->getCurrent()->gridsize > maxgridsize) maxgridsize
              = negimage.innerborder->getCurrent()->gridsize;
            if(negimage.innerborder->getCurrent()->gridsize < mingridsize) mingridsize
              = negimage.innerborder->getCurrent()->gridsize;
            
            critcurve.imagekist->InsertAfterCurrent(negimage.innerborder->getCurrent());
            critcurve.imagekist->Down();
            critcurve.imagekist->getCurrent()->in_image = YES;
            
            negimage.innerborder->Down();
        }
        findborders4(grid->i_tree,&critcurve);
        
        refinements=ImageFinding::refine_grid_kist(lens,grid,&critcurve,1,resolution,2,&newpoint_kist);
        
        if(refinements==0) break;
        //}else free(critcurve->points);
        
        // add new negative points to negpoints
        newpoint_kist.MoveToTop();
        negimage.imagekist->MoveToBottom();
        do{
            x[0] = grid->i_tree->pointlist->current->image->x[0] - x_source[0];
            x[1] = grid->i_tree->pointlist->current->image->x[1] - x_source[1];
            
            if(newpoint_kist.getCurrent()->image->invmag < 0 && r_source*r_source > (x[0]*x[0] + x[1]*x[1]) )
                negimage.imagekist->InsertAfterCurrent(newpoint_kist.getCurrent());
        }while(newpoint_kist.Down());
    }
    
    newpoint_kist.Empty();
    if(critcurve.imagekist->Nunits()){
        critcurve.imagekist->MoveToTop();
        do{critcurve.imagekist->getCurrent()->in_image = NO;}while(critcurve.imagekist->Down());
    }
    
    //std::cout << " number of points in critical curves: " << critcurve.imagekist->Nunits() << std::endl;
    return ;
}



