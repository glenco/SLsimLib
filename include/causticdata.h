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
  void resize(size_t size);

  size_t numberOfCaustics(){return data.size();}
  
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
    double redshift;        /// redshift of source plane
    double crit_center[2];  /// center of critical line in radians
    double crit_radius[3];  /// average radius, smallest radius and largest radius of critical curve
    double crit_area;       /// area of critical curve in radian^2
    
    double caustic_center[2];
    double caustic_radius[3];
    double caustic_area;
  };

  CausticStructure & operator[](size_t index){return data[index];}
private:
  
  int ncolumns;
  void readfile(std::string filename);
  std::vector<CausticStructure> data;
};

#endif /* defined(__SLsimLib__causticdata__) */
