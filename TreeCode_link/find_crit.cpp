/*
 * find_crit.c
 *
 *  Created on: Sep 8, 2009
 *      Author: R.B. Metcalf
 */

#include "slsimlib.h"
#include "map_images.h"

#define NMAXCRITS 1000

using namespace std;
/** \ingroup ImageFinding
 *
 * \brief Finds critical curves and caustics.
 *
 *
 * The critical curve is found by refining the edges of regions of negative magnification.
 * If there are no regions of negative magnification in the original grid the grid is refined
 * around the point of highest kappa.  If there are other points of high kappa that are not of
 * interest one should be careful that the region of interest is not missed.
 */

void ImageFinding::find_crit(
                             LensHndl lens             /// The lens model.
                             ,GridHndl grid            /// The grid.  It must be initialized.
                             ,std::vector<CriticalCurve> &crtcurve     /// Structure to hold critical curve.
                             ,int *Ncrits              /// The number of critical curves found.
                             ,PosType resolution        /// The target resolution that the critical curve is mapped on the image plane.
                             ,PosType invmag_min        /// finds regions with 1/magnification < invmag_min, set to zero for caustics
                             ,bool verbose
                             ){
  
  long i=0;
  short refinements;
  //short spur,closed;
  //PosType maxgridsize,mingridsize;
  std::vector<ImageInfo> negimage(1);
  Kist<Point> newpoint_kist;
  
  
  bool pseuodcaustic = true;
  
  // find kist of points with 1/magnification less than invmag_min
  negimage[0].imagekist->Empty();
  PointList::iterator i_tree_pointlist_current(grid->i_tree->pointlist->Top());
  Point *minpoint = *i_tree_pointlist_current;
  
  for(i=0;i<grid->i_tree->pointlist->size();++i){
    if((*i_tree_pointlist_current)->invmag < invmag_min){
      negimage[0].imagekist->InsertAfterCurrent(*i_tree_pointlist_current);
      negimage[0].imagekist->Down();
    }
    
    // record point of maximum kappa
    if((*i_tree_pointlist_current)->kappa > minpoint->kappa) minpoint = *i_tree_pointlist_current;
    --i_tree_pointlist_current;
  }
  
  if(negimage[0].imagekist->Nunits() == 0){
    if(minpoint->gridsize <= resolution){  // no caustic found at this resolution
      *Ncrits=0;
      return;
    }
    
    // if there is no negative magnification points use maximum mag point
    negimage[0].imagekist->InsertAfterCurrent(minpoint);
  }
  
  divide_images_kist(grid->i_tree,negimage,Ncrits);

  std::vector<ImageInfo> critcurve(*Ncrits);     /// Structure to hold critical curve.  Must be pre-

  //assert(negimage[0].imagekist->CheckInImage(NO));
  for(int ii=0;ii<negimage.size();++ii){
    negimage[ii].imagekist->SetInImage(YES);
    // refine borders until target resolution is reached
    for(;;){
      
      findborders4(grid->i_tree,&negimage[ii]);
      //negimage[ii].imagekist->SetInImage(NO);

      if(negimage[ii].innerborder->Nunits() > 2000) break;
      
      if(verbose) std::printf("find_crit, going into refine_grid\n");

      refinements=ImageFinding::IF_routines::refine_edges(lens,grid,&negimage[ii]
                        ,1,resolution,1,&newpoint_kist,true);

      if(verbose) std::printf("find_crit, came out of refine_grid\n");
      
      if(verbose) cout << "find_crit, came out of refine_grid, Npoints " << critcurve[ii].imagekist->Nunits() << endl;
      
      if(refinements==0) break;
      
      // add new negative points to negpoints
      newpoint_kist.MoveToTop();
      negimage[ii].imagekist->MoveToBottom();
      do{
        if(newpoint_kist.getCurrent()->invmag < invmag_min){
          negimage[ii].imagekist->InsertAfterCurrent(newpoint_kist.getCurrent());
          newpoint_kist.getCurrent()->in_image = YES;
        }
      }while(newpoint_kist.Down());
    }
    
    critcurve[ii].imagekist->copy(negimage[ii].innerborder);
    // set the old regions back to yes incase there are overlaping regions
    newpoint_kist.Empty();
    //negimage[ii].imagekist->SetInImage(YES);
  }
  
  for(int ii=0;ii<negimage.size();++ii){
    negimage[ii].imagekist->SetInImage(NO);
  }
  
  // gather all the curves together and re-divide them to avoid overlaps
  for(int ii=1;ii<negimage.size();++ii) critcurve[0].imagekist->add(critcurve[ii].imagekist);
  if(verbose) ;
  std::printf("find_crit, number of caustic points: %li\n",critcurve[0].imagekist->Nunits());
  divide_images_kist(grid->i_tree,critcurve,Ncrits);
  for(int ii=0;ii<*Ncrits;++ii) critcurve[ii].imagekist->SetInImage(NO);
  
  if(critcurve[0].imagekist->Nunits() == 0) *Ncrits=0;
  
  // ****  Convert the imagekist into a CriticalCurve structure
  
  {
    crtcurve.resize(*Ncrits);
    
    Kist<Point> neighbors;
    
    for(size_t ii=0;ii<*Ncrits;++ii){
      
      // classify critical curve
      
      grid->i_tree->FindAllBoxNeighborsKist(critcurve[ii].imagekist->getCurrent(),&neighbors);
      Kist<Point>::iterator it = neighbors.TopIt();
      while((*it)->invmag < 0 && !it.atend() ) --it;
      if( 1 < ( (*it)->kappa - sqrt( (*it)->gamma[0]*(*it)->gamma[0] + (*it)->gamma[1]*(*it)->gamma[1]) ) ) crtcurve[ii].type = radial;
      else crtcurve[ii].type = tangential;
      
      std::vector<Point *> hull = critcurve[ii].imagekist->copytovector();
      hull = Utilities::concave_hull(hull,10);
      
      crtcurve[ii].critical_curve.resize(hull.size());
      crtcurve[ii].caustic_curve_intersecting.resize(hull.size());
      crtcurve[ii].critical_center[0] = 0;
      crtcurve[ii].critical_center[1] = 0;
      
      for(size_t jj=0;jj<hull.size();++jj){
        crtcurve[ii].critical_curve[jj] = *hull[jj];
        crtcurve[ii].caustic_curve_intersecting[jj] = *(hull[jj]->image);
        crtcurve[ii].critical_center += *hull[jj];
      }
      
      crtcurve[ii].critical_center /= hull.size();
      
      Utilities::windings(crtcurve[ii].critical_center.x,hull.data(),hull.size(),&(crtcurve[ii].critical_area));
      
      
      critcurve[ii].imagekist->TranformPlanes();
      hull = critcurve[ii].imagekist->copytovector();
      hull = Utilities::concave_hull(hull,5);
      
      crtcurve[ii].caustic_curve_outline.resize(hull.size());
      crtcurve[ii].caustic_center[0] = 0;
      crtcurve[ii].caustic_center[1] = 0;
      
      for(size_t jj=0;jj<hull.size();++jj){
        crtcurve[ii].caustic_curve_outline[jj] = *hull[jj];
        crtcurve[ii].caustic_center += *hull[jj];
      }
      
      crtcurve[ii].caustic_center[0] /= hull.size();
      crtcurve[ii].caustic_center[1] /= hull.size();
      
      Utilities::windings(crtcurve[ii].caustic_center.x,hull.data(),hull.size(),&(crtcurve[ii].caustic_area));
      
      crtcurve[ii].caustic_intersections = Utilities::Geometry::intersect(crtcurve[ii].caustic_curve_intersecting);
      
    }
  }
  
  if(pseuodcaustic && negimage[0].imagekist->Nunits() > 1){
    
    std::vector<ImageInfo> pseudocurve(negimage.size());
    //std::vector<CriticalCurve> psecurve;
    const PosType pseudolimit = -100.0;
    int Npseudo = -1;
    Point *current;
    
    
    // Find points within each critical curve that have invmag < pseudolimit
    // If there are none use the minimum invmag value point.
    Point *minmupoint;
    PosType mumin = 0.0;
    Kist<Point> *newpoints = new Kist<Point>;
    Kist<Point> paritypoints;
    bool found;
    
    //pseudocurve[0].imagekist->copy(negimage.imagekist);
    //divide_images_kist(grid->i_tree,pseudocurve,&Npseudo);
    
    for(int ii = 0;ii<negimage.size();++ii){
      mumin = 0.0;
      found = false;
      
      //negimage[ii].imagekist->SetInImage(YES);
      //findborders4(grid->i_tree,&negimage[i]);

      // find if negative region has a radial caustic border that was already detected
      negimage[ii].outerborder->MoveToTop();
      do{
        current = negimage[ii].outerborder->getCurrent();
        if(1 < ( current->kappa - sqrt( current->gamma[0]*current->gamma[0]
                                       + current->gamma[1]*current->gamma[1]) )){
          // radial caustic must already have been found
          found = true;
          break;
        }
      }while(negimage[ii].outerborder->Down());
      
      if(found) continue;
      
      pseudocurve[++Npseudo].imagekist->copy(negimage[ii].imagekist);

      std::cout << " pseudocurve size " << pseudocurve[Npseudo].imagekist->Nunits() << std::endl;
      
      /// remove all but the points below tmp_pseudolimit
      pseudocurve[Npseudo].imagekist->MoveToTop();
      do{
        if(pseudocurve[Npseudo].imagekist->getCurrent()->invmag < mumin){
          minmupoint = pseudocurve[Npseudo].imagekist->getCurrent();
          mumin = pseudocurve[Npseudo].imagekist->getCurrent()->invmag;
        }
        if(pseudocurve[Npseudo].imagekist->getCurrent()->invmag > pseudolimit){
          pseudocurve[Npseudo].imagekist->getCurrent()->in_image = NO;
          pseudocurve[Npseudo].imagekist->TakeOutCurrent();
        }
      }while(pseudocurve[Npseudo].imagekist->Down());
      
      
      // in case one before top was taken out
      if(pseudocurve[Npseudo].imagekist->getCurrent()->invmag > pseudolimit){
        pseudocurve[Npseudo].imagekist->getCurrent()->in_image = NO;
        pseudocurve[Npseudo].imagekist->TakeOutCurrent();
      }
      
      if(pseudocurve[Npseudo].imagekist->Nunits() == 0){
        pseudocurve[Npseudo].imagekist->InsertAfterCurrent(minmupoint);
        pseudocurve[Npseudo].ShouldNotRefine = 0;  // marks that region has not been found
      }else{
        pseudocurve[Npseudo].ShouldNotRefine = 1;
      }
      
      std::cout << "mumin = " << mumin << " pseudocurve size " << pseudocurve[Npseudo].imagekist->Nunits() << std::endl;
      
      
      pseudocurve[Npseudo].imagekist->SetInImage(YES);
      findborders4(grid->i_tree,&pseudocurve[Npseudo]);
      
      while(
            paritypoints.Nunits() == 0 &&
            IF_routines::refine_edges(lens,grid,&pseudocurve[Npseudo],1,0.1*resolution/sqrt(fabs(pseudolimit)),1,newpoints)
            ){
        // update region
        if(pseudocurve[Npseudo].ShouldNotRefine == 0){
          mumin = pseudocurve[Npseudo].imagekist->getCurrent()->invmag;
          minmupoint = pseudocurve[Npseudo].imagekist->getCurrent();
          pseudocurve[Npseudo].imagekist->getCurrent()->in_image = NO;
          pseudocurve[Npseudo].imagekist->TakeOutCurrent();
        }
        
        newpoints->MoveToTop();
        do{
          current = newpoints->getCurrent();
          if(current->invmag < mumin){
            mumin = current->invmag;
            minmupoint = current;
          }
          if(current->invmag < invmag_min)
            negimage[ii].imagekist->InsertAfterCurrent(newpoints->getCurrent());
          
          if(current->invmag < pseudolimit){
            current->in_image = YES;
            pseudocurve[Npseudo].imagekist->InsertAfterCurrent(current);
          }
          
          if( 1 < ( current->kappa - sqrt( current->gamma[0]*current->gamma[0]
                                          + current->gamma[1]*current->gamma[1]) ) ){
            paritypoints.InsertAfterCurrent(current);
          }
          
        }while(newpoints->Down());
        
        std::cout << "mumin = " << mumin << std::endl;
        
        if(pseudocurve[Npseudo].ShouldNotRefine == 0){
          
          if(pseudocurve[Npseudo].imagekist->Nunits() == 0){
            minmupoint->in_image = YES;
            pseudocurve[Npseudo].imagekist->InsertAfterCurrent(minmupoint);
          }else{
            pseudocurve[Npseudo].ShouldNotRefine = 1;
          }
        }
        findborders4(grid->i_tree,&pseudocurve[Npseudo]);
      }  // refinement loop
      
      pseudocurve[Npseudo].imagekist->SetInImage(NO);
      
      if(paritypoints.Nunits() > 0 ){
        // pseudo-caustic detected
        pseudocurve[Npseudo].imagekist->copy(paritypoints);
        
        pseudocurve[Npseudo].imagekist->SetInImage(YES);
        findborders4(grid->i_tree,&pseudocurve[Npseudo]);
        while(
              IF_routines::refine_edges(lens,grid,&pseudocurve[Npseudo],1,0.1*resolution,1,newpoints)
              ){
          
          newpoints->MoveToTop();
          do{
            current = newpoints->getCurrent();
            
            if( 1 < ( current->kappa - sqrt( current->gamma[0]*current->gamma[0]
                                            + current->gamma[1]*current->gamma[1]) ) ){
              pseudocurve[Npseudo].imagekist->InsertAfterCurrent(current);
              current->in_image = YES;
            }
            
            if( current->invmag < invmag_min )
              negimage[ii].imagekist->InsertAfterCurrent(current);
            

          }while(newpoints->Down());
          
          findborders4(grid->i_tree,&pseudocurve[Npseudo]);
        }  // refinement loop
        
        pseudocurve[Npseudo].imagekist->SetInImage(NO);
      
        paritypoints.Empty();

      }else if(mumin > pseudolimit){
        // neither a region with a magnification below pseudolimit or a radiual caustic were found
        --Npseudo;
      }
    }
    ++Npseudo;
    pseudocurve.resize(Npseudo);

    int Nc = crtcurve.size();
    crtcurve.resize(Npseudo+Nc);
    
    for(size_t ii=Nc,i=0;ii<crtcurve.size();++ii,++i){

      Point *current = pseudocurve[i].imagekist->getCurrent();
      bool rad = ( 1 < ( current->kappa - sqrt( current->gamma[0]*current->gamma[0]
                                      + current->gamma[1]*current->gamma[1]) ) );

      std::vector<Point *> hull = pseudocurve[i].innerborder->copytovector();
      hull = Utilities::concave_hull(hull,10);
      
      if(rad){
        crtcurve[ii].critical_curve.resize(hull.size());
        crtcurve[ii].caustic_curve_intersecting.resize(hull.size());
      }else{
        crtcurve[ii].critical_curve.clear();
        crtcurve[ii].caustic_curve_intersecting.clear();
      }
      crtcurve[ii].critical_center[0] = 0;
      crtcurve[ii].critical_center[1] = 0;
      
      for(size_t jj=0;jj<hull.size();++jj){
        if(rad){
          crtcurve[ii].critical_curve[jj] = *hull[jj];
          crtcurve[ii].caustic_curve_intersecting[jj] = *(hull[jj]->image);
        }
        crtcurve[ii].critical_center += *hull[jj];
      }
      
      crtcurve[ii].critical_center /= hull.size();
      
      if(rad) Utilities::windings(crtcurve[ii].critical_center.x,hull.data(),hull.size(),&(crtcurve[ii].critical_area));
      else crtcurve[ii].critical_area = 0.0;
      
      
      pseudocurve[i].imagekist->TranformPlanes();
      hull = pseudocurve[i].imagekist->copytovector();
      hull = Utilities::concave_hull(hull,10);
      
      crtcurve[ii].caustic_curve_outline.resize(hull.size());
      crtcurve[ii].caustic_center[0] = 0;
      crtcurve[ii].caustic_center[1] = 0;
      
      for(size_t jj=0;jj<hull.size();++jj){
        crtcurve[ii].caustic_curve_outline[jj] = *hull[jj];
        crtcurve[ii].caustic_center += *hull[jj];
      }
      
      crtcurve[ii].caustic_center[0] /= hull.size();
      crtcurve[ii].caustic_center[1] /= hull.size();
      
      Utilities::windings(crtcurve[ii].caustic_center.x,hull.data(),hull.size(),&(crtcurve[ii].caustic_area));
      
      crtcurve[ii].type = radial;
      
    }
    
    /***** test lines *******
    if(Npseudo >= 0){
      PosType rmax,rmin,rave;
      psecurve[0].CausticRadius(rmax,rmin,rave);
      std::cout << "caustic " << rmax << " " << rmin << " " << rave << std::endl;
      PixelMap map(psecurve[0].critical_center.x,1000,rmax/500);
      map.AddCurve(psecurve[0].critical_curve,1.0);
      map.AddCurve(psecurve[0].caustic_curve_outline,2.0);
      map.printFITS("!test_pseudo.fits");
      
      psecurve[0].CriticalRadius(rmax,rmin,rave);
      std::cout << "critical " << rmax << " " << rmin << " " << rave << std::endl;
      
    }/**/
  }
  
  for(int ii=0;ii<negimage.size();++ii)
    negimage[ii].imagekist->SetInImage(NO);
  
  *Ncrits = crtcurve.size();
  return ;
}
/*
 
 Uses max point if negative region is not found
 
 The refinement is done in find_crit2() is done while keeping the
 whole regions and then the borders are stripped off.  This is the opposite
 order from find_crit().
 
 Unlike find_crit() there is no pseuodcaustic option.
 */
void ImageFinding::find_crit2(
                              LensHndl lens             /// The lens model.
                              ,GridHndl grid            /// The grid.  It must be initialized.
                              ,std::vector<CriticalCurve> &crtcurve     /// Structure to hold critical curve.  Must be pre-allocated with maxNcrit elements. Stored in critcurve[i].imagekist.
                              ,int *Ncrits              /// The number of critical curves found.
                              ,PosType resolution        /// The target resolution that the critical curve is mapped on the image plane.
                              ,bool *orderingsuccess    /// true if ordering was successful.
                              ,bool ordercurve          /// Order the curve so that it can be drawn or used to find the winding number.
                              ,bool dividecurves        /// Divide the critical curves into seporate curves by whether they are attached
                              ,PosType invmag_min        /// finds regions with 1/magnification < invmag_min
                              ,bool verbose
                              ){
  
  
  std::cerr << "There are known bugs in ImageFinding::find_crit2() that we are trying to remove.  Use ImageFinding::find_crit()."
  << std::endl;
  throw std::runtime_error("Under construction");
  
  long i=0;
  long refinements;
  int Nregions=0;
  //short spur,closed;
  Kist<Point> newpoint_kist,neighborkist;
  
  std::vector<ImageInfo> critcurve(10);
  
  // find kist of points with 1/magnification less than invmag_min
  critcurve[0].imagekist->Empty();
  PointList::iterator i_tree_pointlist_current(grid->i_tree->pointlist->Top());
  Point *minpoint = *i_tree_pointlist_current;
  
  for(i=0;i<grid->i_tree->pointlist->size();++i){
    if((*i_tree_pointlist_current)->invmag < invmag_min){
      critcurve[0].imagekist->InsertAfterCurrent(*i_tree_pointlist_current);
      critcurve[0].imagekist->Down();
    }
    
    // record point of maximum kappa
    if((*i_tree_pointlist_current)->kappa > minpoint->kappa) minpoint = *i_tree_pointlist_current;
    --i_tree_pointlist_current;
  }
  bool maxpoint = false;
  
  if(critcurve[0].imagekist->Nunits() == 0){ // no point with small enough invmag
    
    if(minpoint->gridsize <= resolution){  // no caustic found at this resolution
      *Ncrits=0;
      *orderingsuccess = false;
      return;
    }
    
    // if there is no negative magnification points use maximum mag point
    critcurve[0].imagekist->InsertAfterCurrent(minpoint);
    maxpoint =true;
  }
  
  // divide into regions that are widely seporated
  if(critcurve[0].imagekist->Nunits() >1 ) divide_images_kist(grid->i_tree,critcurve,&Nregions);
  for(int ii=0;ii<Nregions;++ii){
    critcurve[ii].imagekist->SetInImage(YES);
    findborders4(grid->i_tree,&critcurve[ii]);
  }
  
  // loop through seporated regions, this could be done in parrellel
  for(int ii=0;ii<Nregions;++ii){
    
    // refine borders until target resolution is found
    for(int k=0;;++k){
      
      refinements=IF_routines::refine_edges(lens,grid,&critcurve[ii],1,resolution/2,1,&newpoint_kist);
      //refinements=refine_edges(lens,grid,&critcurve[ii],1,1.0e-3,0,&newpoint_kist);
      //refinements=IF_routines::refine_grid_kist(lens,grid,&critcurve[ii],1,resolution,2,&newpoint_kist);
      
      if(!refinements) break;
      critcurve[ii].outerborder->SetInImage(MAYBE);
      
      // add new points to negative region
      newpoint_kist.MoveToTop();
      critcurve[ii].imagekist->MoveToBottom();
      do{
        if(newpoint_kist.getCurrent()->invmag < invmag_min){
          newpoint_kist.getCurrent()->in_image = YES;
          critcurve[ii].imagekist->InsertAfterCurrent(newpoint_kist.getCurrent());
          
          // It is possible that gridrange[] will not be maintained
          if(critcurve[ii].gridrange[1] < newpoint_kist.getCurrent()->gridsize)
            critcurve[ii].gridrange[1] = newpoint_kist.getCurrent()->gridsize;
          
          if(critcurve[ii].gridrange[2] > newpoint_kist.getCurrent()->gridsize)
            critcurve[ii].gridrange[2] = newpoint_kist.getCurrent()->gridsize;
          
        }else{
          newpoint_kist.getCurrent()->in_image = NO;
        }
      }while(newpoint_kist.Down());
      
      if(maxpoint){
        if(critcurve[ii].imagekist->Nunits() > 1){
          // take out old max point
          critcurve[ii].imagekist->MoveToTop();
          do{
            if(critcurve[ii].imagekist->getCurrent()->invmag > 0){
              critcurve[ii].imagekist->getCurrent()->in_image = NO;
              critcurve[ii].imagekist->TakeOutCurrent();
              break;
            }
          }while(critcurve[ii].imagekist->Down());
          maxpoint = false;
        }else{
          // update maximum kappa point if no negative magnification points have been found
          newpoint_kist.MoveToTop();
          do{
            if(newpoint_kist.getCurrent()->kappa
               > critcurve[ii].imagekist->getCurrent()->kappa ){
              critcurve[ii].imagekist->getCurrent()->in_image = NO;
              critcurve[ii].imagekist->TakeOutCurrent();
              newpoint_kist.getCurrent()->in_image = YES;
              critcurve[ii].imagekist->InsertAfterCurrent(newpoint_kist.getCurrent());
            }
          }while(newpoint_kist.Down());
        }
      }
      
      // check which points in inner border are still in the border
      bool ininner;
      unsigned long Ntmp = critcurve[0].innerborder->Nunits();
      critcurve[0].innerborder->MoveToTop();
      for(unsigned long j=0;j < Ntmp;++j){
        
        ininner=false;
        
        grid->i_tree->FindAllBoxNeighborsKist(critcurve[ii].innerborder->getCurrent(),&neighborkist);
        
        neighborkist.MoveToTop();
        do{
          
          if( neighborkist.getCurrent()->in_image != YES){  // point is a neighbor
            ininner=true;
            
            if(neighborkist.getCurrent()->in_image == NO){  // if point is not yet in outerborder
              // add point to outerborder
              neighborkist.getCurrent()->in_image = MAYBE;
              critcurve[ii].outerborder->InsertAfterCurrent(neighborkist.getCurrent());
              critcurve[ii].outerborder->Down();
            }
          }
          
        }while(neighborkist.Down());
        
        if(!ininner){
          bool tmp = critcurve[ii].innerborder->AtTop();
          critcurve[ii].innerborder->TakeOutCurrent();
          if(!tmp) critcurve[ii].innerborder->Down();
        }else{
          critcurve[ii].innerborder->Down();
        }
      }
      
      // Take out outer border points that are no longer in outer border
      critcurve[ii].gridrange[0] = 0.0;
      Ntmp = critcurve[ii].outerborder->Nunits();
      critcurve[ii].outerborder->MoveToTop();
      bool tmpbool;
      for(unsigned long j=0;j<Ntmp;++j){
        
        tmpbool = true;
        assert(critcurve[ii].outerborder->getCurrent()->in_image == MAYBE);
        grid->i_tree->FindAllBoxNeighborsKist(critcurve[ii].outerborder->getCurrent(),&neighborkist);
        neighborkist.MoveToTop();
        do{
          if(neighborkist.getCurrent()->in_image == YES){
            if(critcurve[ii].outerborder->getCurrent()->gridsize
               > critcurve[ii].gridrange[0])
              critcurve[ii].gridrange[0] = critcurve[ii].outerborder->getCurrent()->gridsize;
            tmpbool = false;
            break;
          }
        }while(neighborkist.Down());
        
        if(tmpbool){  // no neighbor in image was found
          bool tmp = critcurve[ii].outerborder->AtTop();
          critcurve[ii].outerborder->getCurrent()->in_image = NO;
          critcurve[ii].outerborder->TakeOutCurrent();
          if(!tmp) critcurve[ii].outerborder->Down();
        }else{
          critcurve[ii].outerborder->Down();
        }
      }
      
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
                  critcurve[ii].innerborder->InsertAfterCurrent(newpoint_kist.getCurrent());
                  tmpbool = false;
                }
                if(neighborkist.getCurrent()->in_image == NO){
                  neighborkist.getCurrent()->in_image = MAYBE;
                  critcurve[ii].outerborder->InsertAfterCurrent(neighborkist.getCurrent());
                }
              }
            }else{
              if(neighborkist.getCurrent()->in_image == YES){
                newpoint_kist.getCurrent()->in_image = MAYBE;
                critcurve[ii].outerborder->InsertAfterCurrent(newpoint_kist.getCurrent());
                break;
              }
            }
          }while(neighborkist.Down());
        }
      }while(newpoint_kist.Down());
      
      critcurve[ii].outerborder->SetInImage(NO);
    }
  }
  
  //  borders are no refined and critcurve[ii].imagekist contains all of the region
  
  if(maxpoint){
 	  *Ncrits = 0;
 	  assert(critcurve[0].imagekist->Nunits() == 1);
 	  critcurve[0].imagekist->getCurrent()->in_image = NO;
 	  critcurve[0].imagekist->Empty();
 	  critcurve[0].outerborder->Empty();
 	  critcurve[0].innerborder->Empty();
 	  return;
  }
  
  
  // make inner border of all regions the image of region 0
  //  This is done so that regions that have grown together during refinement can be joined
  critcurve[0].imagekist->SetInImage(NO);
  critcurve[0].imagekist->Empty();
  critcurve[0].imagekist->copy(critcurve[0].innerborder);
  
  for(int ii=1;ii<Nregions;++ii){
    // make inner borders the image
    critcurve[ii].imagekist->SetInImage(NO);
    critcurve[ii].imagekist->Empty();
    critcurve[0].imagekist->add(critcurve[ii].innerborder);
  }
  
  // Now critcurve[0].imagekist constains all the inner borders
  
  //size_t Npoints = critcurve[0].imagekist->Nunits();
  if(dividecurves) divide_images_kist(grid->i_tree,critcurve,Ncrits);
  else *Ncrits = 1;
  
  for(i=0;i<*Ncrits;++i) critcurve[i].imagekist->SetInImage(NO);
  
  *orderingsuccess = true;
  
  if(critcurve[0].imagekist->Nunits() == 0) *Ncrits=0;
  
  for(i=0;i<*Ncrits;++i){
    if(critcurve[i].imagekist->Nunits() == 0){      // take out curves with no points
      critcurve[i].copy(critcurve[*Ncrits-1],true);
      *Ncrits -= 1;
      --i;
    }
  }
  
  // ****  Convert the imagekist into a CriticalCurve structure
  
  
  crtcurve.resize(*Ncrits);
  
  for(size_t ii=0;ii<*Ncrits;++ii){
    
    
    std::vector<Point *> hull = critcurve[ii].imagekist->copytovector();
    if(ordercurve) hull = Utilities::concave_hull(hull,10);
    
    crtcurve[ii].critical_curve.resize(hull.size());
    crtcurve[ii].caustic_curve_intersecting.resize(hull.size());
    crtcurve[ii].critical_center[0] = 0;
    crtcurve[ii].critical_center[1] = 0;
    
    for(size_t jj=0;jj<hull.size();++jj){
      crtcurve[ii].critical_curve[jj] = *hull[jj];
      crtcurve[ii].caustic_curve_intersecting[jj] = *(hull[jj]->image);
      crtcurve[ii].critical_center += *hull[jj];
    }
    
    crtcurve[ii].critical_center /= hull.size();
    
    Utilities::windings(crtcurve[ii].critical_center.x,hull.data(),hull.size(),&(crtcurve[ii].critical_area));
    
    critcurve[ii].imagekist->TranformPlanes();
    hull = critcurve[ii].imagekist->copytovector();
    if(ordercurve) hull = Utilities::concave_hull(hull,5);
    
    crtcurve[ii].caustic_curve_outline.resize(hull.size());
    crtcurve[ii].caustic_center[0] = 0;
    crtcurve[ii].caustic_center[1] = 0;
    
    for(size_t jj=0;jj<hull.size();++jj){
      crtcurve[ii].caustic_curve_outline[jj] = *hull[jj];
      crtcurve[ii].caustic_center += *hull[jj];
    }
    
    crtcurve[ii].caustic_center[0] /= hull.size();
    crtcurve[ii].caustic_center[1] /= hull.size();
    
    Utilities::windings(crtcurve[ii].caustic_center.x,hull.data(),hull.size(),&(crtcurve[ii].caustic_area));
    
    
  }
  
  
  return ;
}

/**
 *  This is a stripped down version of find_crit() for use in find_image_microlens() that
 *  refines the critical lines that are within the image.
 */
void ImageFinding::IF_routines::refine_crit_in_image(
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
  PointList::iterator i_tree_pointlist_current(grid->i_tree->pointlist->Top());
  for(i=0;i<grid->i_tree->pointlist->size();++i){
    x[0] = (*i_tree_pointlist_current)->image->x[0] - x_source[0];
    x[1] = (*i_tree_pointlist_current)->image->x[1] - x_source[1];
    
    if( (*i_tree_pointlist_current)->invmag < 0 && r_source*r_source > (x[0]*x[0] + x[1]*x[1]) ){
      negimage.imagekist->InsertAfterCurrent(*i_tree_pointlist_current);
      negimage.imagekist->Down();
    }
    --i_tree_pointlist_current;
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
    
    refinements=ImageFinding::IF_routines::refine_grid_kist(lens,grid,&critcurve,1,resolution,2,&newpoint_kist);
    
    if(refinements==0) break;
    //}else free(critcurve->points);
    
    // add new negative points to negpoints
    newpoint_kist.MoveToTop();
    negimage.imagekist->MoveToBottom();
    do{
      x[0] = (*i_tree_pointlist_current)->image->x[0] - x_source[0];
      x[1] = (*i_tree_pointlist_current)->image->x[1] - x_source[1];
      
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






/** \ingroup ImageFinding
 *
 * \brief Finds iso kappa contours.
 *
 *
 */
/* void ImageFinding::find_contour(
 LensHndl lens             /// The lens model.
 ,GridHndl grid            /// The grid.  It must be initialized.
 ,std::vector<CriticalCurve> &crtcurve     /// Structure to hold iso kappa contour.
 ,int *Ncrits              /// The number of critical curves found.
 ,PosType isokappa        /// finds regions with 1/magnification < invmag_min, set to zero for caustics
 ,bool verbose
 ){
 
 long i=0;
 short refinements;
 PosType maxgridsize,mingridsize;
 ImageInfo contour;
 Kist<Point> newpoint_kist;
 
 std::vector<ImageInfo> critcurve(10);     /// Structure to hold critical curve.  Must be pre-
 
 
 std::vector<ImageInfo> pseudocurve(critcurve.size());
 bool pseuodcaustic = false;
 PosType pseudolimit = -0.01;
 PosType kappalimit = 0.1;
 
 // find kist of points with negative magnification
 contour.imagekist->Empty();
 
 PointList::iterator i_tree_pointlist_current(grid->i_tree->pointlist->Top());
 Point *minpoint = *i_tree_pointlist_current;
 
 for(i=0;i<grid->i_tree->pointlist->size();++i){
 if ((*i_tree_pointlist_current)->kappa>isokappa){
 std::cout << (*i_tree_pointlist_current)->kappa << std::endl;
 contour.imagekist->InsertAfterCurrent(*i_tree_pointlist_current);
 contour.imagekist->Down();
 }
 
 // record point of maximum kappa
 if((*i_tree_pointlist_current)->kappa > minpoint->kappa) minpoint = *i_tree_pointlist_current;
 --i_tree_pointlist_current;
 }
 
 findborders4(grid->i_tree,&critcurve[ii])
 
 
 // ****  Convert the kinagekist into a CriticalCurve structure
 
 
 crtcurve.resize(*Ncrits);
 
 for(size_t ii=0;ii<*Ncrits;++ii){
 
 
 std::vector<Point *> hull = critcurve[ii].imagekist->copytovector();
 if(ordercurve) hull = Utilities::concave_hull(hull,10);
 
 crtcurve[ii].critical_curve.resize(hull.size());
 crtcurve[ii].critical_center[0] = 0;
 crtcurve[ii].critical_center[1] = 0;
 
 for(size_t jj=0;jj<hull.size();++jj){
 crtcurve[ii].critical_curve[jj] = *hull[jj];
 crtcurve[ii].critical_center += *hull[jj];
 }
 
 crtcurve[ii].critical_center /= hull.size();
 
 Utilities::windings(crtcurve[ii].critical_center.x,hull.data(),hull.size(),&(crtcurve[ii].critical_area));
 
 critcurve[ii].imagekist->TranformPlanes();
 hull = critcurve[ii].imagekist->copytovector();
 if(ordercurve) hull = Utilities::concave_hull(hull,5);
 
 crtcurve[ii].caustic_curve_outline.resize(hull.size());
 crtcurve[ii].caustic_center[0] = 0;
 crtcurve[ii].caustic_center[1] = 0;
 
 for(size_t jj=0;jj<hull.size();++jj){
 crtcurve[ii].caustic_curve_outline[jj] = *hull[jj];
 crtcurve[ii].caustic_center += *hull[jj];
 }
 
 crtcurve[ii].caustic_center[0] /= hull.size();
 crtcurve[ii].caustic_center[1] /= hull.size();
 
 Utilities::windings(crtcurve[ii].caustic_center.x,hull.data(),hull.size(),&(crtcurve[ii].caustic_area));
 
 
 }
 
 
 return ;
 }
 */

void ImageFinding::find_contour(
                                LensHndl lens             /// The lens model.
                                ,GridHndl grid            /// The grid.  It must be initialized.
                                ,std::vector<CriticalCurve> &crtcurve     /// Structure to hold critical curve.  Must be pre-allocated with maxNcrit elements. Stored in critcurve[i].imagekist.
                                ,int *Ncrits              /// The number of critical curves found.
                                ,PosType resolution        /// The target resolution that the critical curve is mapped on the image plane.
                                ,bool *orderingsuccess    /// true if ordering was successful.
                                ,bool ordercurve          /// Order the curve so that it can be drawn or used to find the winding number.
                                ,bool dividecurves        /// Divide the critical curves into seporate curves by whether they are attached
                                ,double contour_value    /// value at which the contour is wanted
                                ,LensingVariable contour_type  /// KAPPA, INVMAG or DT
                                ,bool verbose
                                ){
  
  std::vector<ImageInfo> critcurve(10);
  
  long i=0;
  long refinements;
  int Nregions=0;
  //short spur,closed;
  Kist<Point> newpoint_kist,neighborkist;
  
  // find kist of points with negative magnification
  critcurve[0].imagekist->Empty();
  PointList::iterator i_tree_pointlist_current(grid->i_tree->pointlist->Top());
  Point *minpoint = *i_tree_pointlist_current;
  
  KappaType value,maxval=0;
  
  for(i=0;i<grid->i_tree->pointlist->size();++i){
    
    switch (contour_type) {
      case KAPPA:
        value = (*i_tree_pointlist_current)->kappa;
        maxval = minpoint->kappa;
        break;
      case INVMAG:
        value = (*i_tree_pointlist_current)->invmag;
        maxval = minpoint->invmag;
        break;
      case DT:
        value = (*i_tree_pointlist_current)->dt;
        maxval = minpoint->dt;
        break;
      default:
        break;
    }
    
    if(value > contour_value){
      critcurve[0].imagekist->InsertAfterCurrent(*i_tree_pointlist_current);
      critcurve[0].imagekist->Down();
    }
    
    // record point of maximum kappa
    if(value > maxval) minpoint = *i_tree_pointlist_current;
    --i_tree_pointlist_current;
  }
  bool maxpoint = false;
  
  if(critcurve[0].imagekist->Nunits() == 0){
    if(minpoint->gridsize <= resolution){  // no caustic found at this resolution
      *Ncrits=0;
      *orderingsuccess = false;
      return;
    }
    
    // if there is no negative magnification points use maximum mag point
    critcurve[0].imagekist->InsertAfterCurrent(minpoint);
    maxpoint =true;
  }
  
  // divide into regions that are widely seporated
  if(critcurve[0].imagekist->Nunits() >1 ) divide_images_kist(grid->i_tree,critcurve,&Nregions);
  for(int ii=0;ii<Nregions;++ii){
    critcurve[ii].imagekist->SetInImage(YES);
    findborders4(grid->i_tree,&critcurve[ii]);
  }
  
  // loop through seporated regions, this could be done in parrellel
  for(int ii=0;ii<Nregions;++ii){
    
    for(int k=0;;++k){
      
      refinements=IF_routines::refine_edges(lens,grid,&critcurve[ii],1,resolution/2,1,&newpoint_kist);
      //refinements=refine_edges(lens,grid,&critcurve[ii],1,1.0e-3,0,&newpoint_kist);
      //refinements=IF_routines::refine_grid_kist(lens,grid,&critcurve[ii],1,resolution,2,&newpoint_kist);
      
      if(!refinements) break;
      critcurve[ii].outerborder->SetInImage(MAYBE);
      
      // add new points to negative region
      newpoint_kist.MoveToTop();
      critcurve[ii].imagekist->MoveToBottom();
      do{
        switch (contour_type) {
          case KAPPA:
            value = newpoint_kist.getCurrent()->kappa;
            break;
          case INVMAG:
            value = newpoint_kist.getCurrent()->invmag;
            break;
          case DT:
            value = newpoint_kist.getCurrent()->dt;
            break;
          default:
            break;
        }
        
        if(value > contour_value){
          newpoint_kist.getCurrent()->in_image = YES;
          critcurve[ii].imagekist->InsertAfterCurrent(newpoint_kist.getCurrent());
          
          // It is possible that gridrange[] will not be maintained
          if(critcurve[ii].gridrange[1] < newpoint_kist.getCurrent()->gridsize)
            critcurve[ii].gridrange[1] = newpoint_kist.getCurrent()->gridsize;
          
          if(critcurve[ii].gridrange[2] > newpoint_kist.getCurrent()->gridsize)
            critcurve[ii].gridrange[2] = newpoint_kist.getCurrent()->gridsize;
          
        }else{
          newpoint_kist.getCurrent()->in_image = NO;
        }
      }while(newpoint_kist.Down());
      
      // check which points in inner border are still in the border
      bool ininner;
      unsigned long Ntmp = critcurve[0].innerborder->Nunits();
      critcurve[0].innerborder->MoveToTop();
      for(unsigned long j=0;j < Ntmp;++j){
        
        ininner=false;
        
        grid->i_tree->FindAllBoxNeighborsKist(critcurve[ii].innerborder->getCurrent(),&neighborkist);
        
        neighborkist.MoveToTop();
        do{
          
          if( neighborkist.getCurrent()->in_image != YES){  // point is a neighbor
            ininner=true;
            
            if(neighborkist.getCurrent()->in_image == NO){  // if point is not yet in outerborder
              // add point to outerborder
              neighborkist.getCurrent()->in_image = MAYBE;
              critcurve[ii].outerborder->InsertAfterCurrent(neighborkist.getCurrent());
              critcurve[ii].outerborder->Down();
            }
          }
          
        }while(neighborkist.Down());
        
        if(!ininner){
          bool tmp = critcurve[ii].innerborder->AtTop();
          critcurve[ii].innerborder->TakeOutCurrent();
          if(!tmp) critcurve[ii].innerborder->Down();
        }else{
          critcurve[ii].innerborder->Down();
        }
      }
      
      // Take out outer border points that are no longer in outer border
      critcurve[ii].gridrange[0] = 0.0;
      Ntmp = critcurve[ii].outerborder->Nunits();
      critcurve[ii].outerborder->MoveToTop();
      bool tmpbool;
      for(unsigned long j=0;j<Ntmp;++j){
        
        tmpbool = true;
        assert(critcurve[ii].outerborder->getCurrent()->in_image == MAYBE);
        grid->i_tree->FindAllBoxNeighborsKist(critcurve[ii].outerborder->getCurrent(),&neighborkist);
        neighborkist.MoveToTop();
        do{
          if(neighborkist.getCurrent()->in_image == YES){
            if(critcurve[ii].outerborder->getCurrent()->gridsize
               > critcurve[ii].gridrange[0])
              critcurve[ii].gridrange[0] = critcurve[ii].outerborder->getCurrent()->gridsize;
            tmpbool = false;
            break;
          }
        }while(neighborkist.Down());
        
        if(tmpbool){  // no neighbor in image was found
          bool tmp = critcurve[ii].outerborder->AtTop();
          critcurve[ii].outerborder->getCurrent()->in_image = NO;
          critcurve[ii].outerborder->TakeOutCurrent();
          if(!tmp) critcurve[ii].outerborder->Down();
        }else{
          critcurve[ii].outerborder->Down();
        }
      }
      
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
                  critcurve[ii].innerborder->InsertAfterCurrent(newpoint_kist.getCurrent());
                  tmpbool = false;
                }
                if(neighborkist.getCurrent()->in_image == NO){
                  neighborkist.getCurrent()->in_image = MAYBE;
                  critcurve[ii].outerborder->InsertAfterCurrent(neighborkist.getCurrent());
                }
              }
            }else{
              if(neighborkist.getCurrent()->in_image == YES){
                newpoint_kist.getCurrent()->in_image = MAYBE;
                critcurve[ii].outerborder->InsertAfterCurrent(newpoint_kist.getCurrent());
                break;
              }
            }
          }while(neighborkist.Down());
        }
      }while(newpoint_kist.Down());
      
      critcurve[ii].outerborder->SetInImage(NO);
    }
  }
  
  // make inner border of all regions the image of region 0
  //  This is done so that regions that have grown together during refinement can be joined
  critcurve[0].imagekist->SetInImage(NO);
  critcurve[0].imagekist->Empty();
  critcurve[0].imagekist->copy(critcurve[0].innerborder);
  
  for(int ii=0;ii<Nregions;++ii){
    // make inner borders the image
    critcurve[ii].imagekist->SetInImage(NO);
    critcurve[ii].imagekist->Empty();
    critcurve[0].imagekist->add(critcurve[ii].innerborder);
    
  }
  
  //size_t Npoints = critcurve[0].imagekist->Nunits();
  if(dividecurves) divide_images_kist(grid->i_tree,critcurve,Ncrits);
  else *Ncrits = 1;
  
  for(i=0;i<*Ncrits;++i) critcurve[i].imagekist->SetInImage(NO);
  
  *orderingsuccess = true;
  
  if(critcurve[0].imagekist->Nunits() == 0) *Ncrits=0;
  
  for(i=0;i<*Ncrits;++i){
    if(critcurve[i].imagekist->Nunits() == 0){      // take out curves with no points
      critcurve[i].copy(critcurve[*Ncrits-1],true);
      *Ncrits -= 1;
      --i;
    }
  }
  
  // ****  Convert the kinagekist into a CriticalCurve structure
  
  
  crtcurve.resize(*Ncrits);
  
  
  for(size_t ii=0;ii<*Ncrits;++ii){
    std::vector<Point *> hull = critcurve[ii].imagekist->copytovector();
    if(ordercurve) hull = Utilities::concave_hull(hull,10);
    
    crtcurve[ii].critical_curve.resize(hull.size());
    crtcurve[ii].critical_center[0] = 0;
    crtcurve[ii].critical_center[1] = 0;
    
    for(size_t jj=0;jj<hull.size();++jj){
      crtcurve[ii].critical_curve[jj] = *hull[jj];
      crtcurve[ii].critical_center += *hull[jj];
      //std::cout << crtcurve[ii].critical_curve[jj].x[0] << " " << crtcurve[ii].critical_curve[jj].x[1]<< std::endl;
      //std::cout << crtcurve[ii].critical_center.x[0] << " " << crtcurve[ii].critical_center.x[1]<< std::endl;
      
      /*tmp_dist=sqrt((crtcurve[ii].critical_curve[jj].x[0]-crtcurve[ii].critical_center.x[0])*(crtcurve[ii].critical_curve[jj].x[0]-crtcurve[ii].critical_center.x[0])+(crtcurve[ii].critical_curve[jj].x[1]-crtcurve[ii].critical_center.x[1])*(crtcurve[ii].critical_curve[jj].x[1]-crtcurve[ii].critical_center.x[1]));
       
       if (nearest_contour_pnt>tmp_dist && tmp_dist!=0){nearest_contour_pnt=tmp_dist;}
       if (mostdistant_contour_pnt<tmp_dist){mostdistant_contour_pnt=tmp_dist;}*/
      
    }
    
    crtcurve[ii].critical_center /= hull.size();
    
    
    Utilities::contour_ellipse(crtcurve[ii].critical_curve ,Utilities::contour_center(crtcurve[ii].critical_curve, hull.size()) , hull.size(), crtcurve[ii].ellipse_curve ,&(crtcurve[ii].contour_ell), &(crtcurve[ii].ellipse_area));
    
    
    Utilities::windings(crtcurve[ii].critical_center.x,hull.data(),hull.size(),&(crtcurve[ii].critical_area));
    
    critcurve[ii].imagekist->TranformPlanes();
    hull = critcurve[ii].imagekist->copytovector();
    if(ordercurve) hull = Utilities::concave_hull(hull,5);
    
    crtcurve[ii].caustic_curve_outline.resize(hull.size());
    crtcurve[ii].caustic_center[0] = 0;
    crtcurve[ii].caustic_center[1] = 0;
    
    for(size_t jj=0;jj<hull.size();++jj){
      crtcurve[ii].caustic_curve_outline[jj] = *hull[jj];
      crtcurve[ii].caustic_center += *hull[jj];
    }
    
    crtcurve[ii].caustic_center[0] /= hull.size();
    crtcurve[ii].caustic_center[1] /= hull.size();
    
    Utilities::windings(crtcurve[ii].caustic_center.x,hull.data(),hull.size(),&(crtcurve[ii].caustic_area));
    
    
    
  }
  
  
  return ;
}


