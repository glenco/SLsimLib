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
/// A class used within CausticData class to store information on a caustic and critical curve pair
struct CausticStructure{
  CausticStructure(){}
  CausticStructure(const CausticStructure &tmp){
    redshift = tmp.redshift;
    crit_center[0] = tmp.crit_center[0];
    crit_center[1] = tmp.crit_center[1];
    crit_radius[0] = tmp.crit_radius[0];
    crit_radius[1] = tmp.crit_radius[1];
    crit_radius[2] = tmp.crit_radius[2];
    crit_area = tmp.crit_area;

    caustic_center[0] = tmp.caustic_center[0];
    caustic_center[1] = tmp.caustic_center[1];
    caustic_radius[0] = tmp.caustic_radius[0];
    caustic_radius[1] = tmp.caustic_radius[1];
    caustic_radius[2] = tmp.caustic_radius[2];
    caustic_area = tmp.caustic_area;

  };
  /// redshift of source plane
  double redshift;
  /// center of critical line in radians
  double crit_center[2];
  /// average radius, smallest radius and largest radius of critical curve
  double crit_radius[3];
  /// area of critical curve in radian^2
  double crit_area;
  
  /// center of caustic curve in radians
  double caustic_center[2];
  /// average radius, smallest radius and largest radius of caustic curve
  double caustic_radius[3];
  /// area of caustic curve in radian^2
  double caustic_area;
};
/** \brief A class for holding, printing and reading the information about the caustics and critical curves in an image.
 */
class CausticData{
  
public:
  CausticData(std::string filename);
  /// creates an empty object 
  CausticData(size_t size);
  ~CausticData();

  /// print to a file,  Produces a file that can be later read with the constructor.
  void printfile(std::string filename,std::string paramfile,double fieldofview,double minscale);

  /// change the number of cautics in the object.
  void resize(size_t size){ data.resize(size);}

  size_t numberOfCaustics(){return data.size();}
  

  CausticStructure & operator[](size_t index){return data[index];}
  void readfile(std::string filename);
  void SortByCritSize();
private:
  
  int ncolumns;
  std::vector<CausticStructure> data;
};

bool comparcritsize(const CausticStructure &caust1,const CausticStructure &caust2);

#endif /* defined(__SLsimLib__causticdata__) */
