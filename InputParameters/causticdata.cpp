//
//  causticdata.cpp
//  SLsimLib
//
//  Created by bmetcalf on 5/6/13.
//
//

#include "causticdata.h"

CausticData::CausticData(std::string filename)
:ncolumns(13),Nxp(0)
{
  readfile(filename);
}

CausticData::CausticData(size_t size)
:ncolumns(13),Nxp(0)
{
  data.resize(size);
}

CausticData::CausticData(const CausticData &input)
:ncolumns(13),Nxp(0)
{
  data = input.data;
}


CausticData::~CausticData(){
  
  if(Nxp){
    //delete searchtree;
    //Utilities::free_PosTypeMatrix(xp,Nxp, 2);
    delete searchtreevec;
  }
  data.clear();
}

/// Read in data from a caustic catalog file
void CausticData::readfile(std::string filename){

  data.clear();

  std::ifstream file_in(filename.c_str());
  std::string myline;
  std::string space = " ";
	double mydouble;

	std::string strg;
	std::string f=",";
	std::stringstream buffer;
  
  if(!file_in){
    std::cout << "Can't open file " << filename << std::endl;
    ERROR_MESSAGE();
    throw std::runtime_error(" Cannot open file.");
  }
  
  std::cout << "Reading caustic information from " << filename << std::endl;
  size_t i=0;
  while(file_in.peek() == '#'){
    file_in.ignore(10000,'\n');
    ++i;
  }
  std::cout << "skipped "<< i << " comment lines in " << filename << std::endl;
  
  size_t pos;
  CausticStructure tmp_data;
  // read in data
  while(getline(file_in,myline)){
    
		if(myline[0] == '#'){
      std::cout << "skipped line " << i << std::endl;
			continue;
    }
	/*
	 if(causticdata[i].crit_radius[0] > 0 ){
	 // *** need to calculate area and max min sizes
	 catalog_caustic << z_sources
	 << " " << causticdata[i].crit_center[0] << " " << causticdata[i].crit_center[1]
	 << " " << causticdata[i].crit_radius[0] << " " << causticdata[i].crit_radius[2]
	 << " "<< causticdata[i].crit_radius[1]  << " " << imageinfo[i].area
	 << " " << causticdata[i].caustic_center[0] << " " << causticdata[i].caustic_center[1]
	 << " " << causticdata[i].caustic_radius[0] << " " << causticdata[i].caustic_radius[2]   << " "
	 << causticdata[i].caustic_radius[1]
	 << " " << causticdata[i].caustic_area
	 << std::endl;
	 }*/

		for(int l=0;l<ncolumns; l++){
			pos = myline.find(space);
			strg.assign(myline,0,pos);
			buffer << strg;
      switch (l) {
        case 0:
          buffer >> mydouble;
          tmp_data.redshift = mydouble;
          break;
        case 1:
          buffer >> mydouble;
          tmp_data.crit_center[0] = mydouble;
          break;
        case 2:
          buffer >> mydouble;
          tmp_data.crit_center[1] = mydouble;
          break;
        case 3:
          buffer >> mydouble;
          tmp_data.crit_radius[0] = mydouble;
          break;
        case 4:
          buffer >> mydouble;
          tmp_data.crit_radius[2] = mydouble;
          break;
        case 5:
          buffer >> mydouble;
          tmp_data.crit_radius[1] = mydouble;
          break;
        case 6:
          buffer >> mydouble;
          tmp_data.crit_area = mydouble;
          break;

        case 7:
          buffer >> mydouble;
          tmp_data.caustic_center[0] = mydouble;
          break;
        case 8:
          buffer >> mydouble;
          tmp_data.caustic_center[1] = mydouble;
          break;
        case 9:
          buffer >> mydouble;
          tmp_data.caustic_radius[0] = mydouble;
          break;
        case 10:
          buffer >> mydouble;
          tmp_data.caustic_radius[2] = mydouble;
          break;
        case 11:
          buffer >> mydouble;
          tmp_data.caustic_radius[1] = mydouble;
          break;
        case 12:
          buffer >> mydouble;
          tmp_data.caustic_area = mydouble;
          break;
          
        default:
          break;
      }
      
      myline.erase(0,pos+1);
      pos= myline.find_first_not_of(space);
      myline.erase(0,pos);
      
      strg.clear();
      buffer.clear();
      buffer.str(std::string());
    }
    myline.clear();
    data.push_back(tmp_data);
  }
  
  /*for(i=0; i < 20; ++i){
    std::cout << data[i].redshift << " " << data[i].crit_center[0] << " " << data[i].crit_center[1]
    << " " << data[i].crit_radius[0] << " " << data[i].crit_radius[1]<< " " << data[i].crit_radius[2] << std::endl;
  }*/

  SetSearchTree();
 }

  /// create tree for searching. Assumes xp is not allocated yet!
void CausticData::SetSearchTree(){
  
  if(Nxp){
    //delete searchtree;
    delete searchtreevec;
    //Utilities::free_PosTypeMatrix(xp,Nxp, 2);
  }
 /* xp = Utilities::PosTypeMatrix(data.size(), 2);
  for(size_t i=0;i<data.size();++i){
    xp[i][0] = data[i].crit_center[0];
    xp[i][1] = data[i].crit_center[1];
  }*/
  //searchtree = new TreeSimple(xp,data.size(),1);
  searchtreevec = new TreeSimpleVec<CausticStructure>(data.data(),data.size(),1,2,true,CausticData::critCenter);
  
  Nxp = data.size();
}

/// Finds the nearest critical curve to the point x[].  If that point is within the largest radius of the critical curve it returns true.
bool CausticData::findNearestCrit(PosType x[2],size_t &index){
  
  if(data.size() == 0){
    index = 0;
    return false;
  }
  if(Nxp == 0 && data.size() != 0){
    SetSearchTree();
  }

  float radius;
  size_t tmp_index;
  //searchtree->NearestNeighbors(x,1,&radius,&tmp_index);
  searchtreevec->NearestNeighbors(x,1,&radius,&index);
   
  return (data[index].crit_radius[2] > radius);
}

void CausticData::printfile(std::string filename,std::string paramfile,double fieldofview,double minscale){
  std::ofstream catalog_caustic(filename.c_str());
  
  catalog_caustic << "# column 1 redshift of source plane" << std::endl;
  catalog_caustic << "# column 2 critical curve center x position in degrees" << std::endl;
  catalog_caustic << "# column 3 critical curve center y position in degrees" << std::endl;
  catalog_caustic << "# column 4 critical curve average radius" << std::endl;
  catalog_caustic << "# column 5 critical curve max radius" << std::endl;
  catalog_caustic << "# column 6 critical curve min radius" << std::endl;
  catalog_caustic << "# column 7 critical curve area" << std::endl;
  catalog_caustic << "# column 8 caustic center x position in degrees" << std::endl;
  catalog_caustic << "# column 9 caustic center y position in degrees" << std::endl;
  catalog_caustic << "# column 10 caustic average radius" << std::endl;
  catalog_caustic << "# column 11 caustic max radius" << std::endl;
  catalog_caustic << "# column 12 caustic min radius" << std::endl;
  catalog_caustic << "# column 13 caustic area" << std::endl;
  catalog_caustic << "# parameter file: " << paramfile << std::endl;
  
  catalog_caustic << "# " << " all critical lines above a scale of " << 180*60*60*minscale/pi << " arcsec,  field of view: " << fieldofview << " square degrees" << std::endl;
  catalog_caustic << "# " << " number of caustics found: " << data.size() << std::endl;


  for(size_t i = 0; i < data.size(); ++i){
    if(data[i].crit_radius[0] != 0){
      catalog_caustic << data[i].redshift
      << " " << data[i].crit_center[0] << " " << data[i].crit_center[1]
      << " " << data[i].crit_radius[0] << " " << data[i].crit_radius[2] << " " << data[i].crit_radius[1]
      << " " << data[i].crit_area
      << " " << data[i].caustic_center[0] << " " << data[i].caustic_center[1]
      << " " << data[i].caustic_radius[0] << " " << data[i].caustic_radius[2] << " " << data[i].caustic_radius[1]
      << " " << data[i].caustic_area
      << std::endl;
    }
  }
}

void CausticData::SortByCritSize(){
  cummulative_area.resize(0);  // this stops RandomLens from being miss used.
  std::sort(data.begin(),data.end(),CausticData::comparcritsize);
  SetSearchTree();
}

void CausticData::SortByCritArea(){
  cummulative_area.resize(0);  // this stops RandomLens from being miss used.
  std::sort(data.begin(),data.end(),CausticData::comparcritarea);
  SetSearchTree();
}

void CausticData::SortByCausticArea(){
  cummulative_area.resize(0);  // this stops RandomLens from being miss used.
  std::sort(data.begin(),data.end(),CausticData::comparcausticarea);
  SetSearchTree();
}

