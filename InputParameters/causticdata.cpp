//
//  causticdata.cpp
//  SLsimLib
//
//  Created by bmetcalf on 5/6/13.
//
//

#include "causticdata.h"

CausticDataStore::CausticDataStore(std::string filename)
:ncolumns(14),Nxp(0)
{
  readfile(filename);
}

CausticDataStore::CausticDataStore(std::vector<ImageFinding::CriticalCurve> &critcurves_vec)
:ncolumns(14),Nxp(0)
{
  PosType rmax,rmin,rave;
  data.resize(critcurves_vec.size());
  for(size_t ii=0;ii<critcurves_vec.size();++ii){
   
    data[ii].redshift = critcurves_vec[ii].z_source;
    data[ii].crit_center[0] = critcurves_vec[ii].critical_center[0];
    data[ii].crit_center[1] = critcurves_vec[ii].critical_center[1];
    critcurves_vec[ii].CriticalRadius(rmax,rmin,rave);
    data[ii].crit_radius[0] = rmax;
    data[ii].crit_radius[2] = rmin;
    data[ii].crit_radius[1] = rave;
    data[ii].crit_area = critcurves_vec[ii].critical_area;
    data[ii].crit_type = critcurves_vec[ii].type;
    
    data[ii].caustic_center[0] = critcurves_vec[ii].caustic_center[0];
    data[ii].caustic_center[1] = critcurves_vec[ii].caustic_center[1];
    critcurves_vec[ii].CausticRadius(rmax,rmin,rave);
    data[ii].caustic_radius[0] = rmax;
    data[ii].caustic_radius[2] = rmin;
    data[ii].caustic_radius[1] = rave;
    data[ii].caustic_area = critcurves_vec[ii].caustic_area;
  }
}

void CausticDataStore::addcrits(std::vector<ImageFinding::CriticalCurve> &critcurves_vec){
  
  PosType rmax,rmin,rave;
  size_t ii = data.size();
  data.resize(ii + critcurves_vec.size());
  for(size_t jj = 0 ;ii < data.size();++ii,++jj){
    
    data[ii].redshift = critcurves_vec[jj].z_source;
    data[ii].crit_center[0] = critcurves_vec[jj].critical_center[0];
    data[ii].crit_center[1] = critcurves_vec[jj].critical_center[1];
    critcurves_vec[jj].CriticalRadius(rmax,rmin,rave);
    data[ii].crit_radius[0] = rmax;
    data[ii].crit_radius[2] = rmin;
    data[ii].crit_radius[1] = rave;
    data[ii].crit_area = critcurves_vec[jj].critical_area;
    data[ii].crit_type = critcurves_vec[jj].type;
    
    data[ii].caustic_center[0] = critcurves_vec[jj].caustic_center[0];
    data[ii].caustic_center[1] = critcurves_vec[jj].caustic_center[1];
    critcurves_vec[jj].CausticRadius(rmax,rmin,rave);
    data[ii].caustic_radius[0] = rmax;
    data[ii].caustic_radius[2] = rmin;
    data[ii].caustic_radius[1] = rave;
    data[ii].caustic_area = critcurves_vec[jj].caustic_area;
  }
  // rebuild search tree
  if(critcurves_vec.size() > 0){
    SetSearchTree();
  }
}


CausticDataStore::CausticDataStore(const CausticDataStore &input)
:ncolumns(14),Nxp(0)
{
  data = input.data;
}


CausticDataStore::~CausticDataStore(){
  
  if(Nxp){
    //delete searchtree;
    //Utilities::free_PosTypeMatrix(xp,Nxp, 2);
    delete searchtreevec;
  }
  data.clear();
}

/// Read in data from a caustic catalog file
void CausticDataStore::readfile(std::string filename){

  data.clear();

  std::ifstream file_in(filename.c_str());
  std::string myline;
  std::string space = " ";
	double mydouble;
  int myint;
  

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
			pos = myline.find("|");
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
          tmp_data.crit_radius[1] = mydouble;
          break;
        case 4:
          buffer >> mydouble;
          tmp_data.crit_radius[0] = mydouble;
          break;
        case 5:
          buffer >> mydouble;
          tmp_data.crit_radius[2] = mydouble;
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
          tmp_data.caustic_radius[1] = mydouble;
          break;
        case 10:
          buffer >> mydouble;
          tmp_data.caustic_radius[0] = mydouble;
          break;
        case 11:
          buffer >> mydouble;
          tmp_data.caustic_radius[2] = mydouble;
          break;
        case 12:
          buffer >> mydouble;
          tmp_data.caustic_area = mydouble;
          break;
        case 13:
          buffer >> myint;
          tmp_data.crit_type = static_cast<CritType>(myint);
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

  if(data.size() > 1) SetSearchTree();
 }

  /// create tree for searching. Assumes xp is not allocated yet!
void CausticDataStore::SetSearchTree(){
  
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
  //searchtreevec = new TreeSimpleVec<CausticStructure>(data.data(),data.size(),1,2,true,crit_center);
  
  searchtreevec = new TreeSimpleVec<CausticStructure>(data.data(),data.size(),1,2,true,[](CausticStructure &c){return c.crit_center.x;});
  Nxp = data.size();
}

/// Finds the nearest critical curve to the point x[].  If that point is within the largest radius of the critical curve it returns true.  If there are no caustics index = -1
bool CausticDataStore::findNearestCrit(PosType *x,long &index){
  
  if(data.size() == 0){
    index = -1;
    return false;
  }
  if(Nxp != data.size() ){
    SetSearchTree();
  }

  float radius;
  size_t tmp_index;
  searchtreevec->NearestNeighbors(x,1,&radius,&tmp_index);
  index = tmp_index;


//searchtreevec->NearestNeighbors(x,1,&radius,&index);
  
  /*********************** test result ***************
  PosType rmin = 1.0e60,r;
  size_t t_index=0;
  
  for(size_t i=0;i<data.size();++i){

      r = (data[i].crit_center[0]-x[0])*(data[i].crit_center[0]-x[0])
      + (data[i].crit_center[1]-x[1])*(data[i].crit_center[1]-x[1]);
      if(rmin > r){
        rmin = r;
        t_index = i;
      }
    
  }
  assert(index == t_index);*/

  return (data[index].crit_radius[2] > radius);

}

/** \brief Finds the nearest critical curve to the point x[] of type 'type'.
 There are two bool types: found returns true if any neighbour was found, found_type is true if a neighbour with correct CritType "type" was found. If there are no caustics of CritType "type" index will be set to -1.
 */
bool CausticDataStore::findNearestCrit(PosType *x,long &index,CritType type){
  
  if(data.size() == 0){
    index = -1;
    return false;
  }
  if(Nxp != data.size() ){
    SetSearchTree();
  }
  
  std::vector<float> radius(10);
  std::vector<IndexType> neighbors(10);
  bool found_type = false;
  int i;

  PosType rmin;
  for(int N=2 ; !found_type && N < data.size(); N += 10){
    
    if(N > data.size()) N = data.size();
    
    if(N > radius.size()){
      radius.resize(N+10);
      neighbors.resize(N+10);
    }
    searchtreevec->NearestNeighbors(x,N,radius.data(),neighbors.data());
    rmin=1E12; //radius[0];
    for(i=0;i<N;++i){
      if(data[neighbors[i]].crit_type == type){
        index = neighbors[i];
        found_type = true;
        rmin = radius[i];
        break;
      }
      if(found_type){break;}
    }
  };
  
  /*PosType rmin = 1.0e60,r;
  size_t t_index;

  for(size_t i=0;i<data.size();++i){
    if(data[i].crit_type == type){
      r = (data[i].crit_center[0]-x[0])*(data[i].crit_center[0]-x[0])
      + (data[i].crit_center[1]-x[1])*(data[i].crit_center[1]-x[1]);
      if(rmin > r){
        rmin = r;
        t_index = i;
      }
    }
  }
  assert(index == t_index);*/
  
  return (data[index].crit_radius[2] > radius[i]);
}

void CausticDataStore::printfile(std::string filename,std::string paramfile,double fieldofview,double minscale){
  std::ofstream catalog_caustic(filename.c_str());
  
  catalog_caustic << "# column 1 redshift of source plane" << std::endl;
  catalog_caustic << "# column 2 critical curve center x position in radians" << std::endl;
  catalog_caustic << "# column 3 critical curve center y position in radians" << std::endl;
  catalog_caustic << "# column 4 critical curve average radius" << std::endl;
  catalog_caustic << "# column 5 critical curve max radius" << std::endl;
  catalog_caustic << "# column 6 critical curve min radius" << std::endl;
  catalog_caustic << "# column 7 critical curve area" << std::endl;
  catalog_caustic << "# column 8 caustic center x position in radians" << std::endl;
  catalog_caustic << "# column 9 caustic center y position in radians" << std::endl;
  catalog_caustic << "# column 10 caustic average radius" << std::endl;
  catalog_caustic << "# column 11 caustic max radius" << std::endl;
  catalog_caustic << "# column 12 caustic min radius" << std::endl;
  catalog_caustic << "# column 13 caustic area" << std::endl;
  catalog_caustic << "# column 14 caustic type" << std::endl;
  catalog_caustic << "# parameter file: " << paramfile << std::endl;
  
  catalog_caustic << "# " << " all critical lines above a scale of " << 180*60*60*minscale/pi << " arcsec,  field of view: " << fieldofview << " square degrees" << std::endl;
  catalog_caustic << "# " << " number of caustics found: " << data.size() << std::endl;


  for(size_t i = 0; i < data.size(); ++i){
      catalog_caustic << data[i].redshift
      << " | " << data[i].crit_center[0] << " | " << data[i].crit_center[1]
      << " | " << data[i].crit_radius[1] << " | " << data[i].crit_radius[0]
      << " | " << data[i].crit_radius[2] << " | " << data[i].crit_area
      << " | " << data[i].caustic_center[0] << " | " << data[i].caustic_center[1]
      << " | " << data[i].caustic_radius[1] << " | " << data[i].caustic_radius[0]
      << " | " << data[i].caustic_radius[2]
      << " | " << data[i].caustic_area << " | " << data[i].crit_type
      << std::endl;
  }
}

/*ImageFinding::log_caustic_data(std::string filename,std::string paramfile,double fieldofview,double minscalestd::vector<ImageFinding::CriticalCurve> critcurves){
  
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

  for(size_t i = 0; i < critcurve.size(); ++i){
    if(critcurve[i].crit_radius[0] != 0){
      catalog_caustic << critcurve[i].redshift
      << " " << critcurve[i].crit_center[0] << " " << critcurve[i].crit_center[1]
      << " " << critcurve[i].crit_radius[0] << " " << critcurve[i].crit_radius[2] << " " << critcurve[i].crit_radius[1]
      << " " << critcurve[i].crit_area
      << " " << critcurve[i].caustic_center[0] << " " << critcurve[i].caustic_center[1]
      << " " << critcurve[i].caustic_radius[0] << " " << critcurve[i].caustic_radius[2] << " " << critcurve[i].caustic_radius[1]
      << " " << critcurve[i].caustic_area
      << std::endl;
    }
  }

}*/

void CausticDataStore::SortByCritSize(){
  cummulative_area.resize(0);  // this stops RandomLens from being miss used.
  std::sort(data.begin(),data.end(),[](const CausticStructure &c1,const CausticStructure &c2){
              return (c1.crit_radius[0] > c2.crit_radius[0]);});

  if(data.size() > 1) SetSearchTree();
}

void CausticDataStore::SortByCritArea(){
  cummulative_area.resize(0);  // this stops RandomLens from being miss used.
  //std::sort(data.begin(),data.end(),CausticDataStore::comparcritarea);

  std::sort(data.begin(),data.end(),[](const CausticStructure &c1,const CausticStructure &c2){
    return (c1.crit_area > c2.crit_area);});
  if(data.size() > 1) SetSearchTree();
}

void CausticDataStore::SortByCausticArea(){
  cummulative_area.resize(0);  // this stops RandomLens from being miss used.
  //std::sort(data.begin(),data.end(),CausticDataStore::comparcausticarea);
  std::sort(data.begin(),data.end(),[](const CausticStructure &c1,const CausticStructure &c2){
    return (c1.caustic_area > c2.caustic_area);});
  if(data.size() > 1) SetSearchTree();
}


std::ostream &operator<<(std::ostream &os, CausticStructure const &caust) {
  return os << " crit center: " << caust.crit_center << " caustic center: "
  << caust.caustic_center << " type : " << caust.crit_type;
}


