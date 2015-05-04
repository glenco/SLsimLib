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


struct CausticStructure{
  CausticStructure(){}
  CausticStructure(const CausticStructure &tmp){
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
  /// redshift of source plane
  double redshift;
  /// center of critical line in radians
  Point_2d crit_center;
  /// average radius, smallest radius and largest radius of critical curve
  double crit_radius[3];
  /// area of critical curve in radian^2
  double crit_area;
  /// caustic type : tangential or radial
  CritType crit_type;
  /// center of caustic curve in radians
  Point_2d caustic_center;
  /// average radius, smallest radius and largest radius of caustic curve
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
  CausticDataStore(std::string filename);
  /// creates an empty object 
  CausticDataStore(std::vector<ImageFinding::CriticalCurve> &crticurve_vec);
  CausticDataStore(const CausticDataStore &input);
  ~CausticDataStore();

 
  /// add some critical curves to the store
  void addcrits(std::vector<ImageFinding::CriticalCurve> &crticurve_vec);
  
  /// print to a file,  Produces a file that can be later read with the constructor.
  void printfile(std::string filename,std::string paramfile,double fieldofview,double minscale);

  /// change the number of cautics in the object.
  void resize(size_t size){ data.resize(size);}

  size_t numberOfCaustics(){return data.size();}
  
  CausticStructure & operator[](size_t index){return data[index];}
  
  void readfile(std::string filename);
  /// sort caustics by size of critical curve radius from largest to smallest
  void SortByCritSize();
  /// sort caustics by size of critical curve area from largest to smallest
  void SortByCritArea();
  /// sort caustics by size of caustic curve area from largest to smallest
  void SortByCausticArea();
  
  /// initialize for selection with CausticDataStore::RandomLens()
  size_t init_for_random(
        short type              /// select according to: (1) critical curve area, (2) caustic curve area
        ,double limit = 0.0     /// minimum accepted area
  ){
    size_t i;
    if(type == 1){
      SortByCritArea();
      cummulative_area.resize(data.size());
      cummulative_area[0] = 0;
      for(i=0;(i<data.size()-1) && (data[i].crit_area > limit) ;++i) cummulative_area[i+1] = data[i].crit_area + cummulative_area[i];
      cummulative_area.resize(i);
    }
    if(type == 2){
      SortByCausticArea();
      cummulative_area.resize(data.size());
      cummulative_area[0] = 0;
      for(i=0;(i<data.size()-1) && (data[i].caustic_area > limit) ;++i) cummulative_area[i+1] = data[i].caustic_area + cummulative_area[i];
      cummulative_area.resize(i);
    }
      
      return cummulative_area.size();
  }
  
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

  int RandomLens(Utilities::RandomNumbers_NR &ran){
    if(cummulative_area.size() == 0 ) throw std::runtime_error("CausticDataStore::RandomLens - CausticDataStore::init_for_random() must be set before using this!  You can also not reorder the caustics without reinitializing");
    return Utilities::locate<double>(cummulative_area,ran()*cummulative_area.back());
  }

  bool findNearestCrit(PosType x[2],long &index);
  bool findNearestCrit(PosType x[2],long &index,CritType type,bool &found_type);
  
private:
  
  int ncolumns;
//  std::vector<CausticStructure> data;
  std::vector<CausticStructure> data;
  std::vector<double> cummulative_area;

  void SetSearchTree();
  //TreeSimple *searchtree;
  TreeSimpleVec<CausticStructure> *searchtreevec;
  //PosType **xp;
  size_t Nxp;
  
  
};

std::ostream &operator<<(std::ostream &os, CausticStructure const &caust);

#endif /* defined(__SLsimLib__causticdata__) */
