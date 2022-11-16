//
//  causticdata.h
//  SLsimLib
//
//  Created by bmetcalf on 5/6/13.
//
//

#ifndef __SLsimLib__causticdata__
#define __SLsimLib__causticdata__

#include <iostream>
#include "standard.h"
#include "utilities_slsim.h"
#include "simpleTree.h"
#include "simpleTreeVec.h"
#include "grid_maintenance.h"


struct CausticSummary{
  CausticSummary(){}
  CausticSummary(const CausticSummary &tmp){
    redshift = tmp.redshift;
    crit_center = tmp.crit_center;
    crit_radius[0] = tmp.crit_radius[0];
    crit_radius[1] = tmp.crit_radius[1];
    crit_radius[2] = tmp.crit_radius[2];
    crit_area = tmp.crit_area;
    crit_type = tmp.crit_type;
    
    caustic_center = tmp.caustic_center;
    caustic_radius[0] = tmp.caustic_radius[0];
    caustic_radius[1] = tmp.caustic_radius[1];
    caustic_radius[2] = tmp.caustic_radius[2];
    caustic_area = tmp.caustic_area;
    
  };
  
  CausticSummary &operator=(CausticSummary &tmp){
    if(&tmp==this) return *this;
    
    redshift = tmp.redshift;
    crit_center = tmp.crit_center;
    crit_radius[0] = tmp.crit_radius[0];
    crit_radius[1] = tmp.crit_radius[1];
    crit_radius[2] = tmp.crit_radius[2];
    crit_area = tmp.crit_area;
    crit_type = tmp.crit_type;
    
    caustic_center = tmp.caustic_center;
    caustic_radius[0] = tmp.caustic_radius[0];
    caustic_radius[1] = tmp.caustic_radius[1];
    caustic_radius[2] = tmp.caustic_radius[2];
    caustic_area = tmp.caustic_area;
    
    return *this;
  }
  /// redshift of source plane
  double redshift;
  /// center of critical line in radians
  Point_2d crit_center;
  /// largest radius, average radius and smallest radius of critical curve
  double crit_radius[3];
  /// area of critical curve in radian^2
  double crit_area;
  /// caustic type : tangential or radial
  CritType crit_type;
  /// center of caustic curve in radians
  Point_2d caustic_center;
  /// largest radius, average radius and smallest radius of caustic curve
  double caustic_radius[3];
  /// area of caustic curve in radian^2
  double caustic_area;
};


/** \brief A class for holding, printing and reading the information about the caustics and critical curves in an image.
 *
 *  Information about individual caustics/critical curves are stored in a CausticStructure class that can be obtaine
 *  with the [] operator.
 */
class CausticDataStore{
  
  /// A class used within CausticDataStore class to store information on a caustic and critical curve pair
  
public:
  CausticDataStore(std::string filename,bool verbose = false);
  /// creates an empty object 
  CausticDataStore(std::vector<ImageFinding::CriticalCurve> &crticurve_vec);
  CausticDataStore(const CausticDataStore &input);
  ~CausticDataStore();

 
  /// add some critical curves to the store
  void addcrits(std::vector<ImageFinding::CriticalCurve> &crticurve_vec);
  
  double getTotalCritArea(){return totalcritarea;}
  double getTotalCausticArea(){return totalcausticarea;}
  
  /// print to a file,  Produces a file that can be later read with the constructor.
  void printfile(std::string filename,std::string paramfile,double fieldofview,double minscale);

  /// change the number of cautics in the object.
  void resize(size_t size){ data.resize(size);}

  size_t numberOfCaustics(){return data.size();}
  
  /// returns in no particular order
  CausticSummary & operator[](size_t index){return data[index];}
  
  /// returns the CausticStructure with n'th largest cirtical curve area
  CausticSummary & CritAreaOrder(size_t n){return data[crit_area_index[n]];}
  /// returns the CausticStructure with n'th largest caustic curve area
  CausticSummary & CausticAreaOrder(size_t n){return data[caus_area_index[n]];}
  /// returns the CausticStructure with n'th largest critical curve radius
  CausticSummary & CritRadiusOrder(size_t n){return data[caus_area_index[n]];}
  
  /// returns the index of the n'th largest critical curve area
  size_t getNthIndexCritArea(size_t n){return crit_area_index[n];}
  /// returns the index of the n'th largest caustic curve area
  size_t getNthIndexCaustArea(size_t n){return caus_area_index[n];}
  /// returns the index of the n'th largest critical curve radius
  size_t getNthIndexCritRadius(size_t n){return crit_radius_index[n];}
 
  
  /// sort caustics by size of critical curve radius from largest to smallest
  //void SortByCritSize();
  /// sort caustics by size of critical curve area from largest to smallest
  //void SortByCritArea();
  /// sort caustics by size of caustic curve area from largest to smallest
  //void SortByCausticArea();
  
  /// initialize for selection with CausticDataStore::RandomLens()
  size_t init_for_random(
        short type              /// select according to: (1) critical curve area, (2) caustic curve area
        ,double limit = 0.0     /// minimum accepted area
  );
  
  /**
   \brief Returns the index for a caustic that is chosen randomly in proportion to area of 
   its critical curve or caustic depending on how init_for_random() was set.
  */
#ifdef ENABLE_CLANG
  int RandomLens(Utilities::RandomNumbers &ran){
    if(cummulative_area.size() == 0 ) throw std::runtime_error("CausticDataStore::RandomLens - CausticDataStore::init_for_random() must be set before using this!   You can also not reorder the caustics without reinitializing.");
      return Utilities::locate<double>(cummulative_area,ran()*cummulative_area.back());
  }
#endif

  int RandomLens(Utilities::RandomNumbers_NR &ran);
  bool findNearestCrit(PosType *x,long &index);
  bool findNearestCrit(PosType *x,long &index,CritType type);
  
  std::vector<CausticSummary> data;
private:
  void readfile(std::string filename,bool verbose);

  int sort_type = -1;
  double totalcritarea = 0.0;
  double totalcausticarea = 0.0;
  std::vector<size_t> crit_area_index;
  std::vector<size_t> crit_radius_index;
  std::vector<size_t> caus_area_index;

  int ncolumns;
  std::vector<double> cummulative_area;

  void SetSearchTree();
  //TreeSimple *searchtree;
  TreeSimpleVec<CausticSummary> *searchtreevec;
  //PosType **xp;
  size_t Nxp;
  
  void constructIndexes();
  
};

std::ostream &operator<<(std::ostream &os, CausticSummary const &caust);

#endif /* defined(__SLsimLib__causticdata__) */
