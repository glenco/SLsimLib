/*
 * find_crit.c
 *
 *  Created on: Sep 8, 2009
 *      Author: R.B. Metcalf
 */

#include "slsimlib.h"
#include "map_images.h"
#include "concave_hull.h"
#include "gridmap.h"


#define NMAXCRITS 1000

using namespace std;
/** Finding
 *
 * \brief Finds critical curves and caustics.
 *
 *
 * The critical curve is found by refining the edges of regions of negative magnification.
 * If there are no regions of negative magnification in the original grid the grid is refined
 * around the point of highest kappa.  If there are other points of high kappa that are not of
 * interest one should be careful that the region of interest is not missed. For this reason
 * critical curves that are smaller than the grid resolution (could be un-uniform)
 * are not guaranteed to be found.
 *
 * After the borders of the negative regions are found the code looks for radial and pseudo
 * caustics within the island.  It usually finds at least one, but if there are more than one
 * per island some might be missed.
 *
 * All the critical curve / caustic pairs are classified as radial, tangential or pseudo.
 * small enough radial critical curve could be miss classified as a pseudo caustic.
 */

void ImageFinding::find_crit(
                             LensHndl lens             /// The lens model.
                             ,GridHndl grid            /// The grid.  It must be initialized.
                             ,std::vector<CriticalCurve> &crtcurve     /// Structure to hold critical curve.
                             ,int *Ncrits              /// The number of critical curves found.
                             ,PosType resolution        /// The target resolution that the critical curve is mapped on the image plane.
                             ,PosType invmag_min        /// finds regions with 1/magnification < invmag_min, set to zero for caustics
                             ,bool verbose
                             ,bool TEST
                             ){
  
  long i=0;
  short refinements;
  //short spur,closed;
  //PosType maxgridsize,mingridsize;
  std::vector<ImageInfo> negimage(1);
  Kist<Point> newpoint_kist;
  bool usingminpoint = false;
  
  if(verbose) std::cout << "****  find_crit() ****" << std::endl;
  
  bool pseuodcaustic = true;
  
  // find kist of points with 1/magnification less than invmag_min
  negimage[0].imagekist->Empty();
  PointList::iterator i_tree_pointlist_it;
  i_tree_pointlist_it.current = (grid->i_tree->pointlist.Top());
  Point *minpoint = *i_tree_pointlist_it;
  
  for(i=0;i<grid->i_tree->pointlist.size();++i){
    if((*i_tree_pointlist_it)->invmag() < invmag_min){
      negimage[0].imagekist->InsertAfterCurrent(*i_tree_pointlist_it);
      negimage[0].imagekist->Down();
    }else{
      (*i_tree_pointlist_it)->in_image = NO;
    }
    
    // record point of maximum kappa
    if((*i_tree_pointlist_it)->kappa() > minpoint->kappa()) minpoint = *i_tree_pointlist_it;
    --i_tree_pointlist_it;
  }
  
  // case where all the points have negative magnification
  if(negimage[0].imagekist->Nunits() == grid->i_tree->pointlist.size()){
    crtcurve.resize(0);
    *Ncrits=0;
    
    return;
  }
  
  if(negimage[0].imagekist->Nunits() == 0){
    if(minpoint->gridsize <= resolution){  // no caustic found at this resolution
      *Ncrits=0;
      if(verbose) std::cout << "********* find_crit() out **************" << std::endl;
      return;
    }
    
    // if there is no negative magnification points use maximum mag point
    negimage[0].imagekist->InsertAfterCurrent(minpoint);
    usingminpoint = true;
  }
  
  /******* test *****************
   {
   Point_2d c = grid->getInitCenter();
   grid->writeFits(c.x, grid->getInitNgrid() ,grid->getInitRange()/grid->getInitNgrid(),INVMAG,"!infind_crit");
   PixelMap map(c.x, grid->getInitNgrid() ,grid->getInitRange()/grid->getInitNgrid());
   map.AddImages(negimage,1,0);
   map.printFITS("!infind_crit_cit");
   }
   *******************************/
  
  
  divide_images_kist(grid->i_tree,negimage,Ncrits);
  /******* test *****************
   {
   Point_2d c = grid->getInitCenter();
   PixelMap map(c.x, grid->getInitNgrid() ,grid->getInitRange()/grid->getInitNgrid());
   map.AddImages(negimage,*Ncrits,0);
   map.printFITS("!infind_crit_cit2");
   }
   *******************************/
  
  std::vector<ImageInfo> critcurve(*Ncrits);     /// Structure to hold critical curve.  Must be pre-
  negimage.resize(*Ncrits);
  std::vector<bool> touches_edge(*Ncrits);
  bool tmpbool;
  if(verbose) std::cout << *Ncrits << " negative islands found." << std::endl;
  
  
  /******* test *****************
   
   Point_2d c = grid->getInitCenter();
   PixelMap map(c.x, grid->getInitNgrid() ,grid->getInitRange()/grid->getInitNgrid());
   **********************************/
  
  //assert(negimage[0].imagekist->CheckInImage(NO));
  for(int ii=0;ii<negimage.size();++ii){
    negimage[ii].imagekist->SetInImage(YES);
    // refine borders until target resolution is reached
    if(verbose) std::cout << "  refining island " << ii << std::endl;
    for(;;){
      
      findborders4(grid->i_tree,&negimage[ii],tmpbool);
      touches_edge[ii] = tmpbool;
      
      /******* test *****************
       map.AddCurve(negimage[ii].outerborder,1);
       // *******************************/
      
      //negimage[ii].imagekist->SetInImage(NO);
      //if(negimage[ii].innerborder->Nunits() > 2000) break;
      
      refinements=ImageFinding::IF_routines::refine_edges(
                                lens,grid,&negimage[ii]
                                ,1,resolution,1,&newpoint_kist,true);
      
      if(refinements==0) break;
      
      if(verbose) cout << "      adding " << newpoint_kist.Nunits() << " points to grid" << endl;
      // add new negative points to negpoints
      newpoint_kist.MoveToTop();
      negimage[ii].imagekist->MoveToBottom();
      if(newpoint_kist.Nunits()>0){
        do{
          if(newpoint_kist.getCurrent()->invmag() < invmag_min){
            if(usingminpoint){
              negimage[ii].imagekist->TakeOutCurrent()->in_image = NO;
              usingminpoint = false;
            }
            negimage[ii].imagekist->InsertAfterCurrent(newpoint_kist.getCurrent());
            newpoint_kist.getCurrent()->in_image = YES;
          }
          
          // update minpoint
          if(usingminpoint && newpoint_kist.getCurrent()->kappa() > minpoint->kappa()) minpoint = newpoint_kist.getCurrent();
          
        }while(newpoint_kist.Down());
      }
      // if no negative island has been found update negimage to minpoint
      if(usingminpoint && minpoint != negimage[ii].imagekist->getCurrent()){
        negimage[ii].imagekist->TakeOutCurrent()->in_image = NO;
        negimage[ii].imagekist->InsertAfterCurrent(minpoint);
      }
    }
  
    critcurve[ii].imagekist->copy(negimage[ii].innerborder);
    // set the old regions back to yes incase there are overlaping regions
    newpoint_kist.Empty();
    negimage[ii].imagekist->SetInImage(YES);
  }
  
  /******* test *****************
   map.printFITS("!infind_crit_outer");
   *******************************/
  
  for(int ii=0;ii<negimage.size();++ii){
    negimage[ii].imagekist->SetInImage(NO);
  }
  
  /******* test *****************
   map.Clean();
   map.AddImages(critcurve,*Ncrits,0);
   map.printFITS("!infind_crit_critcurve");
   map.Clean();
   *******************************/
  
  // gather all the curves together and re-divide them to avoid overlaps
  for(int ii=1;ii<negimage.size();++ii) critcurve[0].imagekist->add(critcurve[ii].imagekist);
  if(verbose)
  std::printf("find_crit, number of caustic points: %li\n",critcurve[0].imagekist->Nunits());
  divide_images_kist(grid->i_tree,critcurve,Ncrits);
  for(int ii=0;ii<*Ncrits;++ii) critcurve[ii].imagekist->SetInImage(NO);
  if(verbose) std::cout << *Ncrits << " borders found." << std::endl;
  
  if(critcurve[0].imagekist->Nunits() == 0) *Ncrits=0;
  
  // ****  Convert the imagekist into a CriticalCurve structure
  
  {
    crtcurve.resize(*Ncrits);
    
    Kist<Point> neighbors;
    size_t ii = 0;
    for(size_t jj=0;jj<*Ncrits;++jj){
      
      if(critcurve[jj].imagekist->Nunits() < 1) continue;
      // classify critical curve
      
      {
        bool catagolized=false;
        critcurve[jj].imagekist->MoveToTop();
        while(!catagolized && !(critcurve[jj].imagekist->AtBottom()) ){
          grid->i_tree->FindAllBoxNeighborsKist(critcurve[jj].imagekist->getCurrent(),&neighbors);
          Kist<Point>::iterator it = neighbors.TopIt();
          while(!it.atend() && (*it).invmag() < 0 ) --it;
          if(it.atend()){ //case where point is on the edge of the field
            critcurve[jj].imagekist->Down();
          }else{
            if( (*it).inverted() ){
              crtcurve[ii].type = CritType::radial;
            }else{
              crtcurve[ii].type = CritType::tangential;
            }
            catagolized=true;
          }
        }
      }
      crtcurve[ii].touches_edge = touches_edge[jj];
      
      /************ test line ****************
       std::cout << "neighbors" << std::endl;
       for(it = neighbors.TopIt(); !it.atend() ; --it){
       std::cout << (*it)->invmag << " " << 1 - ( (*it)->kappa() - sqrt( (*it)->gamma1()*(*it)->gamma1() + (*it)->gamma2()*(*it)->gamma2() ) ) << std::endl;
       }
       ***************************************/
      
      std::vector<Point> points(critcurve[jj].imagekist->Nunits());
      auto iter = critcurve[jj].imagekist->begin();
      for(Point &p : points){
        p = *iter;
        ++iter;
      }
      
      std::vector<Point> hull;
      
      /******* test *****************
       it = critcurve[jj].imagekist->TopIt();
       for(auto pp : hull){
       assert(pp->x[0] == (*it)->x[0]);
       assert(pp->x[1] == (*it)->x[1]);
       --it;
       }
       // *******************************/
      {
        int k=10;
        hull = Utilities::concaveK(points,k);
      }
      
      crtcurve[ii].critcurve.resize(hull.size());
      crtcurve[ii].caustic_curve_intersecting.resize(hull.size());
      crtcurve[ii].critical_center[0] = 0;
      crtcurve[ii].critical_center[1] = 0;

      size_t kk=0;
      for(auto &p : hull){
        crtcurve[ii].critcurve[kk] = p;
        crtcurve[ii].caustic_curve_intersecting[kk] = *(p.image);
        crtcurve[ii].critical_center[0] += p[0];
        crtcurve[ii].critical_center[1] += p[1];
        ++kk;
      }
      crtcurve[ii].critical_center /= crtcurve[ii].critcurve.size();
      
      Utilities::windings(crtcurve[ii].critical_center,crtcurve[ii].critcurve,&(crtcurve[ii].critical_area));
      
      //***************** move to source plane ************/
      
      std::vector<Point_2d> &short_cac = crtcurve[ii].caustic_curve_outline;

      short_cac.resize(points.size());
      
      kk=0;
      PosType scale=0,tmp;
      for(Point &p : points){
        short_cac[kk++] = *(p.image);
        tmp =  p.leaf->area();
        if(scale < tmp) scale = tmp;
      }
      
      //**** size scale ???
      {
        int k = 10;
        short_cac = Utilities::concaveK(short_cac,k);
      }
      
      assert(short_cac.size() > 0);

      // center of caustic
      crtcurve[ii].caustic_center[0] = 0;
      crtcurve[ii].caustic_center[1] = 0;
      for(auto p  : short_cac){
        crtcurve[ii].caustic_center[0] += p[0];
        crtcurve[ii].caustic_center[1] += p[1];
      }
      crtcurve[ii].caustic_center[0] /= short_cac.size();
      crtcurve[ii].caustic_center[1] /= short_cac.size();
      
      
      Utilities::windings(crtcurve[ii].caustic_center,short_cac,&(crtcurve[ii].caustic_area));
      
      crtcurve[ii].caustic_intersections = Utilities::Geometry::intersect(crtcurve[ii].caustic_curve_intersecting);
      
      // take out infinitesimal cases
      if(crtcurve[ii].type == CritType::tangential && crtcurve[ii].critical_area == 0.0) continue;
      if(crtcurve[ii].type != CritType::tangential && crtcurve[ii].caustic_area == 0.0) continue;

      ++ii;
    }
    /******* test *****************
     map.printFITS("!infind_crit_hulled");
     map.Clean();
     // *******************************/
    
    
    *Ncrits = ii;
    crtcurve.resize(*Ncrits);
  }
  
  if(TEST){
    
    //*********************  test lines ****************************
    // This tests that every radial or pseudo critical line is near at
    // least one negative mag point
    Kist<Point> nkist;
    for(auto &crit : crtcurve){
      
      // check that all non-tangent critical line points have a neighbor
      //  with negative magnification
      if(crit.type != CritType::tangential){
        Point *pointp = nullptr;
        if(crit.type == CritType::radial){
          for(RAY &p : crit.critcurve){
            pointp = grid->i_tree->FindBoxPoint(p.x.x);
            grid->i_tree->FindAllBoxNeighborsKist(pointp,&nkist);
            bool good = false;
            for(auto &np : nkist){
              if(np.invmag() < 0){ good = true; break;}
            }
            if(!good){
              
              grid->i_tree->Test();
              grid->s_tree->Test();
              
              std::cout << "Caustic point without negative neighbor"
              << std::endl;
              std::cout << "invmag " << pointp->invmag() << std::endl;
              std::cout << "inverted ? " << pointp->inverted() << std::endl;
              std::cout << "id " << pointp->id << std::endl;
              std::cout << "neighbors: " << std::endl;
              if(pointp->leaf->boundary_p1[0] == grid->i_tree->getTop()->boundary_p1[0])
                std::cout << "At left boundary of grid." << std::endl;
              if(pointp->leaf->boundary_p1[1] == grid->i_tree->getTop()->boundary_p1[1])
                std::cout << "At bottom boundary of grid." << std::endl;
              if(pointp->leaf->boundary_p2[0] == grid->i_tree->getTop()->boundary_p2[0])
                std::cout << "At right boundary of grid." << std::endl;
              if(pointp->leaf->boundary_p2[1] == grid->i_tree->getTop()->boundary_p2[1])
                std::cout << "At top boundary of grid." << std::endl;
              for(auto &np : nkist){
                std::cout << np.id << "    inverted ? " << np.inverted() << std::endl;
              }
              std::cout << " # of points in crit curve: " << crit.critcurve.size() << std::endl;
            }
            assert(good);
          }
        }else if(crit.type == CritType::pseudo){
          pointp = grid->i_tree->FindBoxPoint(crit.critical_center.x);
          grid->i_tree->FindAllBoxNeighborsKist(pointp,&nkist);
          bool good = false;
          for(auto &np : nkist){
            if(np.invmag() < 0){ good = true; break;}
          }
          assert(good);
        }
      }
      
      
      // check for no tangential orphans
      if(crit.type != CritType::tangential){
        
        int count = 0;
        PosType rmax,rmin,rave,r;
        PosType r_closest = grid->i_tree->getTop()->boundary_p2[0]
        - grid->i_tree->getTop()->boundary_p1[0];
        
        PosType crmax,crmin,crave;
        
        CriticalCurve *crit_closest = nullptr;
        
        for(auto &critt : crtcurve){
          if(critt.type == CritType::tangential){
            critt.CriticalRadius(rmax,rmin,rave);
            assert(rmax >= rmin);
            assert(rave >= rmin);
            assert(rmax >= rave);
            r = (critt.critical_center - crit.critical_center).length();
            if( rmax > r) ++count;
            if( r < r_closest){
              r_closest = r;
              crit_closest = &critt;
              crmax = rmax;
              crmin = rmin;
              crave = rave;
            }
          }
        }
        
        Point *pointp;
        if(count == 0){
          std::cout << "Radial or pseudo caustic without tangential partner, plots have been made named ophan*.fits"
          << std::endl;
          pointp = grid->i_tree->FindBoxPoint(crit.critical_center.x);
          grid->i_tree->FindAllBoxNeighborsKist(pointp,&nkist);
          pointp->Print();
          
          std::cout << "closest tangential caustic is " << r_closest << " radians away and has a radius of (" << crmax <<","<<crave<<","<<crmin<<")" << std::endl;
          
          /// make some figures
          Point_2d p1,p2;
          crit.CritRange(p1,p2);
          PosType range = 2.3*r_closest;
          PixelMap map(crit.critical_center.x,1000,range/1000);
          map.AddCurve(crit.critcurve,1.0);
          map.printFITS("!orphin_pseudo.fits");
          
          grid->writeFits(crit.critical_center.x,1000,range/1000,LensingVariable::INVMAG,"!orphin_pseudo");
          map.Clean();
          
          for(auto &critt : crtcurve){
            map.AddCurve(critt.critcurve,1.0);
          }
          
          map.printFITS("!orphin_pseudo_all.fits");
          
        }
        if(count > 1){
          std::cout << "Radial or pseudo caustic has " << count << " tangential critical line within rmax"
          << std::endl;
        }
        
        assert(count > 0);
        assert(count < 2);
        
        if(crit.type == CritType::tangential){
          
          int count = 0;
          PosType rmax,rmin,rave,r;
          PosType r_closest = grid->i_tree->getTop()->boundary_p2[0]
          - grid->i_tree->getTop()->boundary_p1[0];
          
          PosType crmax,crmin,crave;
          
          CriticalCurve *crit_closest = nullptr;
          
          for(auto &critt : crtcurve){
            if(critt.type != CritType::tangential){
              crit.CriticalRadius(rmax,rmin,rave);
              assert(rmax >= rmin);
              assert(rave >= rmin);
              assert(rmax >= rave);
              r = (critt.critical_center - crit.critical_center).length();
              if( rmax > r) ++count;
              if( r < r_closest){
                r_closest = r;
                crit_closest = &critt;
                crmax = rmax;
                crmin = rmin;
                crave = rave;
              }
            }
          }
          
          Point *pointp;
          if(count == 0){
            std::cout << " Tangential caustic without radial or pseudo partner, plots have been made named ophan*.fits"
            << std::endl;
            pointp = grid->i_tree->FindBoxPoint(crit.critical_center.x);
            grid->i_tree->FindAllBoxNeighborsKist(pointp,&nkist);
            pointp->Print();
            
            std::cout << "closest radial or pseudo caustic is " << r_closest << " radians away and tangential radius is (" << crmax <<","<<crave<<","<<crmin<<")" << std::endl;
            
            /// make some figures
            Point_2d p1,p2;
            crit.CritRange(p1,p2);
            PosType range = 2.3*r_closest;
            PixelMap map(crit.critical_center.x,1000,range/1000);
            map.AddCurve(crit.critcurve,1.0);
            map.printFITS("!orphin_pseudo.fits");
            
            grid->writeFits(crit.critical_center.x,1000,range/1000,LensingVariable::INVMAG,"!orphin_pseudo");
            map.Clean();
            
            for(auto &critt : crtcurve){
              map.AddCurve(critt.critcurve,1.0);
            }
            
            map.printFITS("!orphin_pseudo_all.fits");
            
          }
          
          assert(count > 0);
          
        }
        
      }
      
      
    }
    
    std::cout << "No orphan critical curves where found and every critical curve point has a negative magnification neighbour before pseudos are found." << std::endl;
    
    //**************************************************************/
    
  }
  if(pseuodcaustic && negimage[0].imagekist->Nunits() > 1){
    
    
    // regroup the negative islands
    for(int ii=1;ii<negimage.size();++ii) negimage[0] += negimage[ii];
    
    
    // find all negative island again *******************************************
    negimage.resize(1);
    int tmp;
    divide_images_kist(grid->i_tree,negimage,&tmp);
    negimage.resize(tmp);
    for(int ii=0;ii<negimage.size();++ii){
      findborders4(grid->i_tree,&negimage[ii],tmpbool);
      touches_edge[ii]=tmpbool;
    }
    if(verbose) std::cout << " found " << tmp << " negative islands in re-sorting." << std::endl;
    //******************************************************************************
    
    std::vector<ImageInfo> pseudocurve(negimage.size());
    std::vector<CritType> types(negimage.size());
    const PosType pseudolimit = -100.0;
    int Npseudo = 0;
    Point *current;
    
    // Find points within each critical curve that have invmag < pseudolimit
    // If there are none use the minimum invmag value point.
    bool found;
    Kist<Point> paritypoints;
    
    int nfound = 0,Nnotdefined = 0;
    for(int ii = 0;ii<negimage.size();++ii){
      found = false;
      paritypoints.Empty();
      
      // find if negative region has a radial caustic border that was already detected
      negimage[ii].outerborder->MoveToTop();
      do{
        current = negimage[ii].outerborder->getCurrent();
        //if(1 < ( current->kappa() - sqrt( current->gamma1()*current->gamma1()
        //                               + current->gamma2()*current->gamma2()) )){
        
        if( current->inverted() ){
          // radial caustic must already have been found
          found = true;
          break;
        }
      }while(negimage[ii].outerborder->Down());
      
      if(found){
        ++nfound;
        if(verbose) std::cout << "Radial caustic already found in negative island " << ii << std::endl;
        continue;
      }
      types[Npseudo] = ImageFinding::find_pseudo(pseudocurve[Npseudo],negimage[ii]
                                                 ,pseudolimit,lens,grid
                                                 ,resolution,paritypoints,false);
      if(types[Npseudo] != CritType::ND ) ++Npseudo;
      else ++Nnotdefined;
      
    }
    
    if(verbose) std::cout << "  " << nfound << " radial caustics found with tangential caustics " << std::endl;
    if(verbose) std::cout << "  " << Npseudo << " additional radial or pseudo caustics found after further refinement " << std::endl;
    if(verbose) std::cout << "  " << Nnotdefined << " not defined caustics " << std::endl;
    
    
    
    pseudocurve.resize(Npseudo);
    
    int Nc = crtcurve.size();
    crtcurve.resize(Npseudo+Nc);
    
    size_t ii,i;
    // convert to CriticalCurve structure
    for(ii=Nc-1,i=0;i<Npseudo;++i){
      
      //Point *current = pseudocurve[i].imagekist->getCurrent();
      
      if(types[i] == CritType::ND) continue;
      ++ii;
      
      crtcurve[ii].type = types[i];
      
      //std::vector<Point *> hull = pseudocurve[i].innerborder->copytovector();
      //std::vector<Point *> hull = pseudocurve[i].outerborder->copytovector();
      
      
      std::vector<Point> points(pseudocurve[i].outerborder->Nunits());
      auto iter = pseudocurve[i].outerborder->begin();
      PosType scale = 0,tmp;
      for(Point &p : points){
        p = *iter;
        tmp = p.leaf->area();
        if(scale < tmp ) scale = tmp;
        ++iter;
      }
      
      if(verbose) std::cout << " doing concave hull with " << points.size() << " points..." << std::endl;
      
      std::vector<Point> hull;
      
      {
        int k=10;
        hull = Utilities::concaveK(points,k);
      }
      
      crtcurve[ii].critcurve.resize(hull.size());
      crtcurve[ii].caustic_curve_intersecting.resize(hull.size());
      crtcurve[ii].critical_center[0] = 0;
      crtcurve[ii].critical_center[1] = 0;
      size_t kk=0;
      for(auto &p : hull){
        crtcurve[ii].critcurve[kk] = p;
        crtcurve[ii].caustic_curve_intersecting[kk] = *(p.image);
        crtcurve[ii].critical_center[0] += p[0];
        crtcurve[ii].critical_center[1] += p[1];
        ++kk;
      }
      crtcurve[ii].critical_center /= crtcurve[ii].critcurve.size();
      
      Utilities::windings(crtcurve[ii].critical_center,crtcurve[ii].critcurve,&(crtcurve[ii].critical_area));
      
      //***************** move to source plane ************/
      
      std::vector<Point_2d> &short_cac = crtcurve[ii].caustic_curve_outline;
      
      short_cac.resize(points.size());
      
      kk=0;
      scale = 0;
      for(Point &p : points){
        short_cac[kk++] = *(p.image);
        tmp = p.image->leaf->area();
        if(scale < tmp ) scale = tmp;
      }
    
      {
//        int k=10;
//        short_cac = Utilities::concaveK<Point_2d>(short_cac,k);
//        *** try deintersection instead
        
        Utilities::RemoveIntersections(short_cac);
      }
      
      assert(short_cac.size() > 0);
      
      // center of caustic
      crtcurve[ii].caustic_center[0] = 0;
      crtcurve[ii].caustic_center[1] = 0;
      for(auto p  : short_cac){
        crtcurve[ii].caustic_center[0] += p[0];
        crtcurve[ii].caustic_center[1] += p[1];
      }
      crtcurve[ii].caustic_center[0] /= short_cac.size();
      crtcurve[ii].caustic_center[1] /= short_cac.size();
      
      
      if(crtcurve[ii].type != CritType::pseudo)
        Utilities::windings(crtcurve[ii].caustic_center,short_cac,&(crtcurve[ii].caustic_area));
      else crtcurve[ii].critical_area = 0.0;

      
      crtcurve[ii].caustic_intersections = Utilities::Geometry::intersect(crtcurve[ii].caustic_curve_intersecting);
      
      
      /*if (crtcurve[ii].type == tangential){
        hull = Utilities::concave_hull(hull,10);
      }else{
        hull = Utilities::convex_hull(hull);
      }
      
      //hull = Utilities::concave_hull(hull,10);
      // hull = Utilities::convex_hull(hull);
      assert(hull.size() <= pseudocurve[i].outerborder->Nunits());
      
      if(crtcurve[ii].type != CritType::pseudo){
        crtcurve[ii].critical_curve.resize(hull.size());
        crtcurve[ii].caustic_curve_intersecting.resize(hull.size());
      }else{
        crtcurve[ii].critical_curve.clear();
        crtcurve[ii].caustic_curve_intersecting.clear();
      }
      crtcurve[ii].critical_center[0] = 0;
      crtcurve[ii].critical_center[1] = 0;
      
      for(size_t jj=0;jj<hull.size();++jj){
        if(crtcurve[ii].type != CritType::pseudo){
          crtcurve[ii].critical_curve[jj] = *hull[jj];
          crtcurve[ii].caustic_curve_intersecting[jj] = *(hull[jj]->image);
        }
       }
      
      if(hull.size()){
        Point_2d po;
        po = *hull[0];
        for(size_t jj=0;jj<hull.size();++jj){
          crtcurve[ii].critical_center -= po - *hull[jj];
        }
        
        crtcurve[ii].critical_center /= hull.size();
        crtcurve[ii].critical_center += po;
      }
      
      if(crtcurve[ii].type != CritType::pseudo)
        Utilities::windings(crtcurve[ii].critical_center.x,hull.data(),hull.size(),&(crtcurve[ii].critical_area));
      else crtcurve[ii].critical_area = 0.0;
      */
      // caustic
      /*pseudocurve[i].imagekist->TranformPlanes();
      hull = pseudocurve[i].imagekist->copytovector();
      if(verbose) std::cout << " doing concave hull with " << hull.size() << " points..." << std::endl;
      if (crtcurve[ii].type == tangential){
        hull = Utilities::concave_hull(hull,10);
      }else{
        hull = Utilities::convex_hull(hull);
      }*/
      
      
       /*crtcurve[ii].caustic_curve_outline.resize(hull.size());
       size_t i = 0;
       for(Point* &p : hull){
       crtcurve[ii].caustic_curve_outline[i++] = *(p->image);
       }
       
       std::vector<Point_2d> &short_caus_curve = crtcurve[ii].caustic_curve_outline;
       
       crtcurve[ii].caustic_intersections =
       Utilities::RemoveIntersections(short_caus_curve);
*/
      
      //hull = Utilities::concave_hull(hull,10);
      // hull = Utilities::convex_hull(hull);
      
      //crtcurve[ii].caustic_curve_outline.resize(hull.size());
      /*crtcurve[ii].caustic_center[0] = 0;
      crtcurve[ii].caustic_center[1] = 0;
      // center of caustic
      for(auto p  : short_caus_curve){
        crtcurve[ii].caustic_center[0] += p[0];
        crtcurve[ii].caustic_center[1] += p[1];
      }*/
      /*for(size_t jj=0;jj<hull.size();++jj){
        crtcurve[ii].caustic_curve_outline[jj] = *hull[jj];
        crtcurve[ii].caustic_center[0] += hull[jj]->x[0];
        crtcurve[ii].caustic_center[1] += hull[jj]->x[1];
      }*/
      /*
       crtcurve[ii].caustic_center[0] /= hull.size();
      crtcurve[ii].caustic_center[1] /= hull.size();
      */
      //Utilities::windings(crtcurve[ii].caustic_center.x,hull.data(),hull.size(),&(crtcurve[ii].caustic_area));
      
      //Utilities::windings(crtcurve[ii].caustic_center,short_caus_curve,&(crtcurve[ii].caustic_area));

      
      //Utilities::windings(crtcurve[ii].caustic_center.x,hull.data(),hull.size(),&(crtcurve[ii].caustic_area));
    
      // take out infinitesimal cases
      if(crtcurve[ii].type == CritType::tangential && crtcurve[ii].critical_area == 0.0) --ii;
      if(crtcurve[ii].type != CritType::tangential && crtcurve[ii].caustic_area == 0.0) --ii;
    }
    
    // remove cases that were ND type
    crtcurve.resize(ii+1);
  }
  
  
  for(int ii=0;ii<negimage.size();++ii)
    negimage[ii].imagekist->SetInImage(NO);
  
  *Ncrits = crtcurve.size();
  
  if(TEST){
    
    //*********************  test lines ****************************
    // This tests that every every radial or pseudo critical line is near at
    // least one negative mag point
    Kist<Point> nkist;
    for(auto &crit : crtcurve){
      
      // check that all non-tangent critical line points have a neighbor
      //  with negative magnification
      if(crit.type != CritType::tangential){
        Point *pointp = nullptr;
        if(crit.type == CritType::radial){
          for(RAY &p : crit.critcurve){
            pointp = grid->i_tree->FindBoxPoint(p.x.x);
            grid->i_tree->FindAllBoxNeighborsKist(pointp,&nkist);
            bool good = false;
            for(auto &np : nkist){
              if(np.invmag() < 0){ good = true; break;}
            }
            if(!good){
              
              grid->i_tree->Test();
              grid->s_tree->Test();
              
              std::cout << "Critical point without negative neighbor"
              << std::endl;
              std::cout << "invmag " << pointp->invmag() << std::endl;
              std::cout << "inverted ? " << pointp->inverted() << std::endl;
              std::cout << "id " << pointp->id << std::endl;
              std::cout << "neighbors: " << std::endl;
              if(pointp->leaf->boundary_p1[0] == grid->i_tree->getTop()->boundary_p1[0])
                std::cout << "At left boundary of grid." << std::endl;
              if(pointp->leaf->boundary_p1[1] == grid->i_tree->getTop()->boundary_p1[1])
                std::cout << "At bottom boundary of grid." << std::endl;
              if(pointp->leaf->boundary_p2[0] == grid->i_tree->getTop()->boundary_p2[0])
                std::cout << "At right boundary of grid." << std::endl;
              if(pointp->leaf->boundary_p2[1] == grid->i_tree->getTop()->boundary_p2[1])
                std::cout << "At top boundary of grid." << std::endl;
              for(auto &np : nkist){
                std::cout << np.id << "    inverted ? " << np.inverted() <<
                " invmag " << np.invmag() << std::endl;
                np.Print();
              }
              std::cout << " # of points in crit curve: " << crit.critcurve.size()
              << " type: " << to_string(crit.type)
              << std::endl;
            }
            //assert(good);
          }
        }else if(crit.type == CritType::pseudo){
          pointp = grid->i_tree->FindBoxPoint(crit.critical_center.x);
          grid->i_tree->FindAllBoxNeighborsKist(pointp,&nkist);
          bool good = false;
          for(auto &np : nkist){
            if(np.invmag() < 0){ good = true; break;}
          }
          //assert(good);
        }
      }
      
      
      // check for no tangential orphans
      if(crit.type != CritType::tangential){
        
        int count = 0;
        PosType rmax,rmin,rave,r;
        PosType r_closest = grid->i_tree->getTop()->boundary_p2[0]
        - grid->i_tree->getTop()->boundary_p1[0];
        
        PosType crmax,crmin,crave;
        
        CriticalCurve *crit_closest = nullptr;
        
        for(auto &critt : crtcurve){
          if(critt.type == CritType::tangential){
            critt.CriticalRadius(rmax,rmin,rave);
            assert(rmax >= rmin);
            assert(rave >= rmin);
            assert(rmax >= rave);
            r = (critt.critical_center - crit.critical_center).length();
            if( rmax > r) ++count;
            if( r < r_closest){
              r_closest = r;
              crit_closest = &critt;
              crmax = rmax;
              crmin = rmin;
              crave = rave;
            }
          }
        }
        
        Point *pointp = nullptr;
        if(count == 0){
          std::cout << to_string(crit.type) << " caustic without tangential partner, plots have been made named ophan*.fits"
          << std::endl;
          
          
          pointp = grid->i_tree->FindBoxPoint(crit.critical_center.x);
          grid->i_tree->FindAllBoxNeighborsKist(pointp,&nkist);
          pointp->Print();
          
          std::cout << "closest tangential caustic is " << r_closest << " radians away and has a radius of (" << crmax <<","<<crave<<","<<crmin<<")" << std::endl;
          
          /// make some figures
          Point_2d p1,p2;
          crit.CritRange(p1,p2);
          PosType range = 2.3*r_closest;
          PixelMap map(crit.critical_center.x,1000,range/1000);
          map.AddCurve(crit.critcurve,1.0);
          map.printFITS("!orphin_pseudo.fits");
          
          grid->writeFits(crit.critical_center.x,1000,range/1000,LensingVariable::INVMAG,"!orphin_pseudo");
          map.Clean();
          
          for(auto &critt : crtcurve){
            map.AddCurve(critt.critcurve,1.0);
          }
          
          map.printFITS("!orphin_pseudo_all.fits");
          
        }
        if(count > 1){
          std::cout << to_string(crit.type) << " caustic has " << count << " tangential critical line within rmax"
          << std::endl;
        }
        
        assert(count > 0 || crit_closest->critcurve.size() < 3);  // no partners
        assert(count < 2);  // more than one partner
        
        if(crit.type == CritType::tangential){
          
          int count = 0;
          PosType rmax,rmin,rave,r;
          PosType r_closest = grid->i_tree->getTop()->boundary_p2[0]
          - grid->i_tree->getTop()->boundary_p1[0];
          
          PosType crmax,crmin,crave;
          
          CriticalCurve *crit_closest = nullptr;
          
          for(auto &critt : crtcurve){
            if(critt.type != CritType::tangential){
              crit.CriticalRadius(rmax,rmin,rave);
              assert(rmax >= rmin);
              assert(rave >= rmin);
              assert(rmax >= rave);
              r = (critt.critical_center - crit.critical_center).length();
              if( rmax > r) ++count;
              if( r < r_closest){
                r_closest = r;
                crit_closest = &critt;
                crmax = rmax;
                crmin = rmin;
                crave = rave;
              }
            }
          }
          
          Point *pointp;
          if(count == 0){
            std::cout << " Tangential caustic without radial or pseudo partner, plots have been made named ophan*.fits"
            << std::endl;
            pointp = grid->i_tree->FindBoxPoint(crit.critical_center.x);
            grid->i_tree->FindAllBoxNeighborsKist(pointp,&nkist);
            pointp->Print();
            
            std::cout << "closest radial or pseudo caustic is " << r_closest << " radians away and tangential radius is (" << crmax <<","<<crave<<","<<crmin<<")" << std::endl;
            
            /// make some figures
            Point_2d p1,p2;
            crit.CritRange(p1,p2);
            PosType range = 2.3*r_closest;
            PixelMap map(crit.critical_center.x,1000,range/1000);
            map.AddCurve(crit.critcurve,1.0);
            map.printFITS("!orphin_pseudo.fits");
            
            grid->writeFits(crit.critical_center.x,1000,range/1000,LensingVariable::INVMAG,"!orphin_pseudo");
            map.Clean();
            
            for(auto &critt : crtcurve){
              map.AddCurve(critt.critcurve,1.0);
            }
            
            map.printFITS("!orphin_pseudo_all.fits");
            
          }
          
          assert(count > 0);
          
        }
        
      }
      
      
    }
    
    std::cout << "No orphan critical curves where found and every critical curve point has a negative magnification neighbour." << std::endl;
    
    //**************************************************************/
    
  }

  if(verbose) std::cout << "********* find_crit() out **************" << std::endl;
  
  return;
}

void ImageFinding::find_crit(
                             Lens &lens             /// The lens model.
                             ,GridMap &gridmap            /// The grid.  It must be initialized.
                             ,std::vector<CriticalCurve> &critcurves     /// Structure to hold critical curve.
                             ,bool verbose
                             ){
  
  double z_source = lens.getSourceZ();
  
  std::vector<std::vector<Point_2d> > crits;
  std::vector<bool> hits_boundary;
  std::vector<CritType> crit_type;
  
  gridmap.find_crit(crits,hits_boundary,crit_type);
  
  int Ncrit=crits.size();
  critcurves.resize(Ncrit);
  for(size_t i=0 ; i < Ncrit ; ++i){
    critcurves[i].type = crit_type[i];
    critcurves[i].touches_edge = hits_boundary[i];
 
    size_t Npoints = crits[i].size();
    
    critcurves[i].critcurve.resize(Npoints);
    critcurves[i].critical_center *= 0;
    for(size_t j=0 ; j < Npoints ; ++j){
      critcurves[i].critcurve[j].x = crits[i][j];
      critcurves[i].critical_center += crits[i][j];
      critcurves[i].critcurve[j].z = z_source;
    }
    critcurves[i].critical_center /= Npoints;
 
    Utilities::windings(critcurves[i].critical_center,critcurves[i].critcurve,&(critcurves[i].critical_area));
   
    
    lens.rayshooterInternal(Npoints,critcurves[i].critcurve.data());
    
    critcurves[i].caustic_center *= 0;
    critcurves[i].caustic_curve_intersecting.resize(Npoints);
    for(size_t j=0 ; j < Npoints ; ++j){
      critcurves[i].caustic_curve_intersecting[j] = critcurves[i].critcurve[j].y;
      critcurves[i].caustic_center += critcurves[i].critcurve[j].y;
    }
    critcurves[i].caustic_center /= Npoints;
    
    critcurves[i].caustic_curve_outline = Utilities::Geometry::MagicHull(critcurves[i].caustic_curve_intersecting);

    Utilities::windings(critcurves[i].caustic_center,critcurves[i].caustic_curve_outline,&(critcurves[i].caustic_area));
    critcurves[i].caustic_intersections = -1;
    
  }

}

/*  This function is not meant for an external user.  It is only used by
 find_crit(). paritypoints must be empty on first entry.
 */
CritType ImageFinding::find_pseudo(ImageInfo &pseudocurve,ImageInfo &negimage
                                   ,PosType pseudolimit,LensHndl lens,GridHndl grid
                                   ,PosType resolution,Kist<Point> &paritypoints,bool TEST){
  
  Kist<Point> newpoints;
  Point *current;
  
  // case where radial caustic has been detected, but not refined
  if(paritypoints.Nunits() > 0 ){
    // pseudo-caustic detected
    pseudocurve.imagekist->copy(paritypoints);
    
    pseudocurve.imagekist->SetInImage(YES);
    bool touches_edge;
    findborders4(grid->i_tree,&pseudocurve,touches_edge);
    if(TEST){
      for(auto &p : *(pseudocurve.imagekist) ){
        assert( p.inverted());
        assert( p.invmag() > 0 );
      }

      size_t i=0;
      for(auto &p : *(pseudocurve.outerborder) ){
        if( p.invmag() > 0 ){
          p.Print();
          std::cout << "invmag recalculated = " << (1-p.kappa())*(1-p.kappa()) - (p.gamma1()*p.gamma1()) - (p.gamma2()*p.gamma2()) + p.gamma3()*p.gamma3() <<
          std::endl;
        }
        ++i;
        assert( p.invmag() < 0 );
      }
      for(auto &p : *(pseudocurve.innerborder) ){
        assert( p.invmag() > 0 );
      }
   }

    while(pseudocurve.innerborder->Nunits() < 2000 &&
          IF_routines::refine_edges(lens,grid,&pseudocurve,1,0.1*resolution,1,&newpoints)
          ){
      
      newpoints.MoveToTop();
      do{
        current = newpoints.getCurrent();
        
        if( current->inverted() ){
          pseudocurve.imagekist->InsertAfterCurrent(current);
          current->in_image = YES;
        }else current->in_image = NO;
        
        if( current->invmag() < 0 )
          negimage.imagekist->InsertAfterCurrent(current);
        
        
      }while(newpoints.Down());
      
      findborders4(grid->i_tree,&pseudocurve,touches_edge);
      
      if(TEST){
        for(auto &p : *(pseudocurve.outerborder) ){
          assert( p.invmag() < 0 );
        }
        for(auto &p : *(pseudocurve.innerborder) ){
          assert( p.invmag() > 0 );
        }
      }
    }  // refinement loop
    
    pseudocurve.imagekist->SetInImage(NO);
    
    paritypoints.Empty();
    pseudocurve.ShouldNotRefine = 0;

    return CritType::radial;
  }
  
  
  // case where no radial caustic has been detected yet
  pseudocurve.imagekist->copy(negimage.imagekist);
  Point *minmupoint = pseudocurve.imagekist->getCurrent();
  PosType mumin = minmupoint->invmag();
  
  //std::cout << " pseudocurve size " << pseudocurve.imagekist->Nunits() << std::endl;
  
  /// remove all but the points below tmp_pseudolimit
  pseudocurve.imagekist->MoveToTop();
  do{
    if(pseudocurve.imagekist->getCurrent()->invmag() < mumin){
      minmupoint = pseudocurve.imagekist->getCurrent();
      mumin = pseudocurve.imagekist->getCurrent()->invmag();
    }
    if(pseudocurve.imagekist->getCurrent()->invmag() > pseudolimit){
      pseudocurve.imagekist->getCurrent()->in_image = NO;
      pseudocurve.imagekist->TakeOutCurrent();
    }
  }while(pseudocurve.imagekist->Down());
  
  
  // in case one before top was taken out
  if( pseudocurve.imagekist->Nunits() > 0 && pseudocurve.imagekist->getCurrent()->invmag() > pseudolimit){
    pseudocurve.imagekist->getCurrent()->in_image = NO;
    pseudocurve.imagekist->TakeOutCurrent();
  }
  
  if(pseudocurve.imagekist->Nunits() == 0){
    pseudocurve.imagekist->InsertAfterCurrent(minmupoint);
    pseudocurve.ShouldNotRefine = 0;  // marks that region has not been found
  }else{
    pseudocurve.ShouldNotRefine = 1;
  }
  
  //std::cout << "mumin = " << mumin << " pseudocurve size " << pseudocurve.imagekist->Nunits() << std::endl;
  
  pseudocurve.imagekist->SetInImage(YES);
  bool touches_edge;
  findborders4(grid->i_tree,&pseudocurve,touches_edge);
  
  while(
        paritypoints.Nunits() == 0 && pseudocurve.imagekist->Nunits() < 200 &&
        //IF_routines::refine_edges(lens,grid,&pseudocurve,1,0.1*resolution/sqrt(fabs(pseudolimit)),1,&newpoints)  &&
        IF_routines::refine_grid_kist(lens,grid,&pseudocurve,1, 0.1*resolution/sqrt(fabs(pseudolimit)),2,&newpoints)
        ){
    // update region
    if(pseudocurve.ShouldNotRefine == 0){
      mumin = pseudocurve.imagekist->getCurrent()->invmag();
      minmupoint = pseudocurve.imagekist->getCurrent();
      pseudocurve.imagekist->getCurrent()->in_image = NO;
      pseudocurve.imagekist->TakeOutCurrent();
    }
    
    newpoints.MoveToTop();
    do{
      current = newpoints.getCurrent();
      if(current->invmag() < mumin){
        mumin = current->invmag();
        minmupoint = current;
      }
      if(current->invmag() < 0)
        negimage.imagekist->InsertAfterCurrent(newpoints.getCurrent());
      
      if(current->invmag() < pseudolimit){
        current->in_image = YES;
        pseudocurve.imagekist->InsertAfterCurrent(current);
      }
      
      if( current->inverted() ){
        paritypoints.InsertAfterCurrent(current);
      }
      
    }while(newpoints.Down());
    
    //std::cout << "mumin = " << mumin << std::endl;
    
    if(pseudocurve.ShouldNotRefine == 0){
      
      if(pseudocurve.imagekist->Nunits() == 0){
        minmupoint->in_image = YES;
        pseudocurve.imagekist->InsertAfterCurrent(minmupoint);
      }else{
        pseudocurve.ShouldNotRefine = 1;
      }
    }
    findborders4(grid->i_tree,&pseudocurve,touches_edge);
  }  // refinement loop
  
  pseudocurve.imagekist->SetInImage(NO);
  
  // radial caustic was detected
  if(paritypoints.Nunits() > 0){
    return find_pseudo(pseudocurve,negimage,pseudolimit,lens,grid,resolution,paritypoints,TEST);
  }
  // neither a region with a magnification below pseudolimit or a radiual caustic were found
  //  now minimize the largest eigenvalue and see if it is negative
  if(mumin > pseudolimit){
    PosType eigmin,tmp;
    Point *pmin;
    
    pseudocurve.imagekist->Empty();
    
    // find point with minimum largest Eigenvalue
    negimage.imagekist->MoveCurrentToTop();
    current = negimage.imagekist->getCurrent();
    eigmin = 1 - current->kappa() + sqrt( current->gamma1()*current->gamma1()
                                      + current->gamma2()*current->gamma2()) ;
    //eigmin = 1-current->kappa() + sqrt( (1-current->kappa())*(1-current->kappa())
    //                                 + current->invmag()) ;
    assert(eigmin == eigmin);
    
    if(eigmin < 0) paritypoints.InsertAfterCurrent(current);
    
    pmin = negimage.imagekist->getCurrent();
    while(negimage.imagekist->Down()){
      current = negimage.imagekist->getCurrent();
      tmp = 1 - current->kappa() + sqrt( current->gamma1()*current->gamma1()
                                     + current->gamma2()*current->gamma2()) ;
     //tmp = 1-current->kappa() + sqrt( (1-current->kappa())*(1-current->kappa())
     //                               + current->invmag()) ;
      
      assert(tmp == tmp);
      
      if(eigmin > tmp ){
        eigmin = tmp;
        pmin = current;
      }
      
      if(tmp < 0) paritypoints.InsertAfterCurrent(current);
    }

    pseudocurve.imagekist->InsertAfterCurrent(pmin);
    pmin->in_image = YES;
    bool touches_edge;
    findborders4(grid->i_tree,&pseudocurve,touches_edge);
    
    while(eigmin >= 0 &&
          IF_routines::refine_edges(lens,grid,&pseudocurve,1,0.01*resolution,1,&newpoints)){
      
      newpoints.MoveToTop();
      do{
        current = newpoints.getCurrent();
        tmp = 1 - current->kappa() + sqrt( current->gamma1()*current->gamma1()
                                       + current->gamma2()*current->gamma2()) ;
        //tmp = 1-current->kappa() + sqrt( (1-current->kappa())*(1-current->kappa())
        //                              + current->invmag()) ;
        assert(tmp == tmp);
       if(eigmin > tmp ){
          eigmin = tmp;
          pmin->in_image = NO;
          pmin = current;
          pmin->in_image = YES;
          pseudocurve.imagekist->Empty();
          pseudocurve.imagekist->InsertAfterCurrent(pmin);
        }
        if(current->in_image < 0){
          negimage.imagekist->InsertAfterCurrent(current);
        }
        if(tmp < 0) paritypoints.InsertAfterCurrent(current);
        
      }while(newpoints.Down());
      
      findborders4(grid->i_tree,&pseudocurve,touches_edge);
    }
    
    pseudocurve.ShouldNotRefine = 0;
    if( eigmin < 0){
      assert(paritypoints.Nunits() > 0);
      // a radial caustic has been detected, repeat
      return find_pseudo(pseudocurve,negimage,pseudolimit,lens,grid,resolution,paritypoints);
    }else{
      // everything has failed
      assert(paritypoints.Nunits() == 0);
      pseudocurve.Empty();
      return CritType::ND;
    }
  }
  
  pseudocurve.ShouldNotRefine = 0;
  // pseudo-caustic was found
  return CritType::pseudo;
}


/*
 
 Uses max point if negative region is not found
 
 The refinement is done in find_crit2() is done while keeping the
 whole regions and then the borders are stripped off.  This is the opposite
 order from find_crit().
 
 Unlike find_crit() there is no pseuodcaustic option.
 *
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
 PointList::iterator i_tree_pointlist_it(grid->i_tree->pointlist->Top());
 Point *minpoint = *i_tree_pointlist_it;
 
 for(i=0;i<grid->i_tree->pointlist->size();++i){
 if((*i_tree_pointlist_it)->invmag() < invmag_min){
 critcurve[0].imagekist->InsertAfterCurrent(*i_tree_pointlist_it);
 critcurve[0].imagekist->Down();
 }
 
 // record point of maximum kappa
 if((*i_tree_pointlist_it)->kappa() > minpoint->kappa()) minpoint = *i_tree_pointlist_it;
 --i_tree_pointlist_it;
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
 if(newpoint_kist.getCurrent()->invmag() < invmag_min){
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
 if(critcurve[ii].imagekist->getCurrent()->invmag() > 0){
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
 if(newpoint_kist.getCurrent()->kappa()
 > critcurve[ii].imagekist->getCurrent()->kappa() ){
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
 */
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
  bool touches_edge;
  
  // find kist of points with negative magnification
  negimage.imagekist->Empty();
  PointList::iterator i_tree_pointlist_it;
  i_tree_pointlist_it.current = (grid->i_tree->pointlist.Top());
  for(i=0;i<grid->i_tree->pointlist.size();++i){
    x[0] = (*i_tree_pointlist_it)->image->x[0] - x_source[0];
    x[1] = (*i_tree_pointlist_it)->image->x[1] - x_source[1];
    
    if( (*i_tree_pointlist_it)->invmag() < 0 && r_source*r_source > (x[0]*x[0] + x[1]*x[1]) ){
      negimage.imagekist->InsertAfterCurrent(*i_tree_pointlist_it);
      negimage.imagekist->Down();
    }
    --i_tree_pointlist_it;
  }
  
  if(negimage.imagekist->Nunits() == 0) return;
  
  for(;;){
    
    negimage.imagekist->MoveToTop();
    do{negimage.imagekist ->getCurrent()->in_image = YES;} while(negimage.imagekist->Down());
    
    findborders4(grid->i_tree,&negimage,touches_edge);
    
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
    findborders4(grid->i_tree,&critcurve,touches_edge);
    
    refinements=ImageFinding::IF_routines::refine_grid_kist(lens,grid,&critcurve,1,resolution,2,&newpoint_kist);
    
    if(refinements==0) break;
    //}else free(critcurve->points);
    
    // add new negative points to negpoints
    newpoint_kist.MoveToTop();
    negimage.imagekist->MoveToBottom();
    do{
      x[0] = (*i_tree_pointlist_it)->image->x[0] - x_source[0];
      x[1] = (*i_tree_pointlist_it)->image->x[1] - x_source[1];
      
      if(newpoint_kist.getCurrent()->image->invmag() < 0 && r_source*r_source > (x[0]*x[0] + x[1]*x[1]) )
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






/** Finding
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
 
 PointList::iterator i_tree_pointlist_it(grid->i_tree->pointlist->Top());
 Point *minpoint = *i_tree_pointlist_it;
 
 for(i=0;i<grid->i_tree->pointlist->size();++i){
 if ((*i_tree_pointlist_it)->kappa()>isokappa){
 std::cout << (*i_tree_pointlist_it)->kappa() << std::endl;
 contour.imagekist->InsertAfterCurrent(*i_tree_pointlist_it);
 contour.imagekist->Down();
 }
 
 // record point of maximum kappa
 if((*i_tree_pointlist_it)->kappa() > minpoint->kappa()) minpoint = *i_tree_pointlist_it;
 --i_tree_pointlist_it;
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
                                ,LensingVariable contour_type  /// KAPPA, INVMAG or DELAYT
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
  PointList::iterator i_tree_pointlist_it;
  i_tree_pointlist_it.current = (grid->i_tree->pointlist.Top());
  Point *minpoint = *i_tree_pointlist_it;
  
  KappaType value,maxval=0;
  
  for(i=0;i<grid->i_tree->pointlist.size();++i){
    
    switch (contour_type) {
      case LensingVariable::KAPPA:
        value = (*i_tree_pointlist_it)->kappa();
        maxval = minpoint->kappa();
        break;
      case LensingVariable::INVMAG:
        value = (*i_tree_pointlist_it)->invmag();
        maxval = minpoint->invmag();
        break;
      case LensingVariable::DELAYT:
        value = (*i_tree_pointlist_it)->dt;
        maxval = minpoint->dt;
        break;
      default:
        break;
    }
    
    if(value > contour_value){
      critcurve[0].imagekist->InsertAfterCurrent(*i_tree_pointlist_it);
      critcurve[0].imagekist->Down();
    }
    
    // record point of maximum kappa
    if(value > maxval) minpoint = *i_tree_pointlist_it;
    --i_tree_pointlist_it;
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
  bool touches_edge;
  for(int ii=0;ii<Nregions;++ii){
    critcurve[ii].imagekist->SetInImage(YES);
    findborders4(grid->i_tree,&critcurve[ii],touches_edge);
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
          case LensingVariable::KAPPA:
            value = newpoint_kist.getCurrent()->kappa();
            break;
          case LensingVariable::INVMAG:
            value = newpoint_kist.getCurrent()->invmag();
            break;
          case LensingVariable::DELAYT:
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
    
    crtcurve[ii].critcurve.resize(hull.size());
    crtcurve[ii].critical_center[0] = 0;
    crtcurve[ii].critical_center[1] = 0;
    
    for(size_t jj=0;jj<hull.size();++jj){
      crtcurve[ii].critcurve[jj] = *hull[jj];
      crtcurve[ii].critical_center -= crtcurve[ii].critcurve[0].x - *hull[jj];
    }
    
    crtcurve[ii].critical_center /= hull.size();
    crtcurve[ii].critical_center += crtcurve[ii].critcurve[0].x;
    
    Utilities::contour_ellipse(crtcurve[ii].critcurve ,Utilities::contour_center(crtcurve[ii].critcurve, hull.size()) , hull.size(), crtcurve[ii].ellipse_curve ,&(crtcurve[ii].contour_ell), &(crtcurve[ii].ellipse_area));
    
    
    Utilities::windings(crtcurve[ii].critical_center.x,hull.data(),hull.size(),&(crtcurve[ii].critical_area));
    
    critcurve[ii].imagekist->TranformPlanes();
    hull = critcurve[ii].imagekist->copytovector();
    if(ordercurve) hull = Utilities::concave_hull(hull,10);
    
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

std::string to_string(CritType crit){
  std::string s;
  switch (crit) {
    case CritType::ND:
      s = "NotDefined";
      break;
    case CritType::radial:
      s = "radial";
      break;
    case CritType::tangential:
      s = "tangential";
      break;
    case CritType::pseudo:
      s = "pseudo";
      break;
    default:
      break;
  }
  
  return s;
}


