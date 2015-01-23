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
                             ,std::vector<CriticalCurve> &crtcurve     /// Structure to hold critical curve.
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
  PosType maxgridsize,mingridsize;
  ImageInfo negimage;
  Kist<Point> newpoint_kist;
  
  std::vector<ImageInfo> critcurve(10);     /// Structure to hold critical curve.  Must be pre-
  
  
  std::vector<ImageInfo> pseudocurve(critcurve.size());
  bool pseuodcaustic = false;
  PosType pseudolimit = -0.01;
  
  // find kist of points with negative magnification
  negimage.imagekist->Empty();

  PointList::iterator i_tree_pointlist_current(grid->i_tree->pointlist->Top());
  Point *minpoint = *i_tree_pointlist_current;
  
  for(i=0;i<grid->i_tree->pointlist->size();++i){
    if((*i_tree_pointlist_current)->invmag < invmag_min){
      negimage.imagekist->InsertAfterCurrent(*i_tree_pointlist_current);
      negimage.imagekist->Down();
    }
    
    // record point of maximum kappa
    if((*i_tree_pointlist_current)->kappa > minpoint->kappa) minpoint = *i_tree_pointlist_current;
    --i_tree_pointlist_current;
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
    critcurve[0].imagekist->Empty();
    negimage.innerborder->MoveToTop();
    for(i=0,maxgridsize=0.0,mingridsize=1.0e99;i<negimage.innerborder->Nunits();++i){
      
      if(negimage.innerborder->getCurrent()->gridsize > maxgridsize) maxgridsize = negimage.innerborder->getCurrent()->gridsize;
      if(negimage.innerborder->getCurrent()->gridsize < mingridsize) mingridsize = negimage.innerborder->getCurrent()->gridsize;
      
      critcurve[0].imagekist->InsertAfterCurrent(negimage.innerborder->getCurrent());
      critcurve[0].imagekist->Down();
      critcurve[0].imagekist->getCurrent()->in_image = YES;
      
      negimage.innerborder->Down();
    }
    findborders4(grid->i_tree,critcurve.data());
    //std::printf("came out of findborders 2\n");
    
    if(verbose) std::printf("find_crit, going into refine_grid\n");
    //std::printf("  Npoints=%i\n",critcurve->Npoints);
    //refinements=refine_grid(lens,grid->i_tree,grid->s_tree,critcurve,1,resolution,2,false);
    refinements=ImageFinding::IF_routines::refine_grid_kist(lens,grid,critcurve.data(),1,resolution,2,&newpoint_kist,true);
    if(verbose) std::printf("find_crit, came out of refine_grid\n");
    
    if(verbose) cout << "Npoints " << critcurve[0].imagekist->Nunits() << endl;
    
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
  if(critcurve[0].imagekist->Nunits()){
    critcurve[0].imagekist->MoveToTop();
    do{critcurve[0].imagekist->getCurrent()->in_image = NO;}while(critcurve[0].imagekist->Down());
  }
  
  if(verbose) std::printf("find_crit, number of caustic points: %li\n",critcurve[0].imagekist->Nunits());
  
  if(dividecurves) divide_images_kist(grid->i_tree,critcurve,Ncrits);
  else *Ncrits = 1;
  
  if(pseuodcaustic){
    // Find points within each critical curve that have invmag < pseudolimit
    // If there are none use the minimum invmag value point.
    Point *minmupoint;
    PosType mumin = 0.0;
    Kist<Point> *newpoints = new Kist<Point>;
    
    pseudocurve[0].imagekist->copy(negimage.imagekist);
    divide_images_kist(grid->i_tree,pseudocurve,Ncrits);
    
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
              IF_routines::refine_edges(lens,grid,&pseudocurve[i],1,0.01*resolution/sqrt(fabs(pseudolimit)),1,newpoints)
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
 
  if(critcurve[0].imagekist->Nunits() == 0) *Ncrits=0;
  
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
      critcurve[i].copy(critcurve[*Ncrits-1],true);
      *Ncrits -= 1;
      --i;
    }
  }
  // ****  Convert the imagekist into a CriticalCurve structure
  
  {
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
      
      crtcurve[ii].caustic_curve.resize(hull.size());
      crtcurve[ii].caustic_center[0] = 0;
      crtcurve[ii].caustic_center[1] = 0;
      
      for(size_t jj=0;jj<hull.size();++jj){
        crtcurve[ii].caustic_curve[jj] = *hull[jj];
        crtcurve[ii].caustic_center += *hull[jj];
      }
      
      crtcurve[ii].caustic_center[0] /= hull.size();
      crtcurve[ii].caustic_center[1] /= hull.size();
      
      Utilities::windings(crtcurve[ii].caustic_center.x,hull.data(),hull.size(),&(crtcurve[ii].caustic_area));
      
      
    }
  }
  
  return ;
}
/*
 
 Uses max point if negative region is not found
 
 Does not have the pseuodcaustic option
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
      }
      
      crtcurve[ii].critical_center /= hull.size();
      
      Utilities::windings(crtcurve[ii].critical_center.x,hull.data(),hull.size(),&(crtcurve[ii].critical_area));
      
      critcurve[ii].imagekist->TranformPlanes();
      hull = critcurve[ii].imagekist->copytovector();
       if(ordercurve) hull = Utilities::concave_hull(hull,5);
      
      crtcurve[ii].caustic_curve.resize(hull.size());
      crtcurve[ii].caustic_center[0] = 0;
      crtcurve[ii].caustic_center[1] = 0;
      
      for(size_t jj=0;jj<hull.size();++jj){
        crtcurve[ii].caustic_curve[jj] = *hull[jj];
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
    
    crtcurve[ii].caustic_curve.resize(hull.size());
    crtcurve[ii].caustic_center[0] = 0;
    crtcurve[ii].caustic_center[1] = 0;
    
    for(size_t jj=0;jj<hull.size();++jj){
      crtcurve[ii].caustic_curve[jj] = *hull[jj];
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
                              ,double kappa_val        /// finds regions with 1/magnification < invmag_min
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
  
  for(i=0;i<grid->i_tree->pointlist->size();++i){
    if((*i_tree_pointlist_current)->kappa > kappa_val){
      critcurve[0].imagekist->InsertAfterCurrent(*i_tree_pointlist_current);
      critcurve[0].imagekist->Down();
    }
    
    // record point of maximum kappa
    if((*i_tree_pointlist_current)->kappa > minpoint->kappa) minpoint = *i_tree_pointlist_current;
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
        if(newpoint_kist.getCurrent()->kappa > kappa_val){
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
            if(critcurve[ii].imagekist->getCurrent()->kappa > 0){
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
    
    crtcurve[ii].caustic_curve.resize(hull.size());
    crtcurve[ii].caustic_center[0] = 0;
    crtcurve[ii].caustic_center[1] = 0;
    
    for(size_t jj=0;jj<hull.size();++jj){
      crtcurve[ii].caustic_curve[jj] = *hull[jj];
      crtcurve[ii].caustic_center += *hull[jj];
    }
    
    crtcurve[ii].caustic_center[0] /= hull.size();
    crtcurve[ii].caustic_center[1] /= hull.size();
    
    Utilities::windings(crtcurve[ii].caustic_center.x,hull.data(),hull.size(),&(crtcurve[ii].caustic_area));
    
    
    
  }
  
  
  return ;
}


