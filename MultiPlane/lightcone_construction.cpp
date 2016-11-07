//
//  lightcone_construction.cpp
//  GLAMER
//
//  Created by Ben Metcalf on 26/10/16.
//
//

#include "lightcone_construction.h"

LightCone::LightCone(
                     double angular_radius    /// angular radius of lightcone in radians
):
r_theta(angular_radius)
{
  sin_theta_sqrt = pow(sin(r_theta),2);
};

/** \brief Read in a lightcone data file and render them into NFW halos
 *
 *   filename file should be the output of LightCone::WriteLightCone() to
 *   be in the right format.
 */
void LightCone::ReadLightConeNFW(
                                 std::string filename   /// name of file to read
                                 ,COSMOLOGY &cosmo      /// cosmology for changing distance to redshift
                                 ,std::vector<LensHalo* > &lensVec  /// output LensHalos
){
  
  Utilities::delete_container(lensVec);
  
  std::ifstream file(filename.c_str());
  std::string myline;
  
  while(file.peek() == '#') file.ignore(10000,'\n');
  
  const int ncolumns = 7;
  void *addr[ncolumns];
  size_t haloid;
  double mass,Rvir,Rscale;
  Utilities::Geometry::SphericalPoint sph_point;
  Point_3d tmp_point;
  
  addr[0] = &haloid;
  addr[1] = &(tmp_point.x[0]);
  addr[2] = &(tmp_point.x[1]);
  addr[3] = &(tmp_point.x[2]);
  addr[4] = &mass;
  addr[5] = &Rvir;
  addr[6] = &Rscale;
  
  size_t myint;
  std::string strg;
  std::string delim=",";
  double myPosType;
  double z;
  
  // skip first line
  std::getline(file,myline);
  
  std::stringstream buffer;
  
  // read a line of data
  while(std::getline(file,myline)){
    
    if(myline[0] == '#')
      break;
    for(int l=0;l<ncolumns; l++){
      int pos = myline.find(delim);
      strg.assign(myline,0,pos);
      buffer << strg;
      if(l == 0){
        buffer >> myint;
        *((size_t *)addr[l]) = myint;
      }else{
        buffer >> myPosType;
        *((PosType *)addr[l]) = myPosType;
      }
      myline.erase(0,pos+1);
      strg.clear();
      buffer.clear();
      buffer.str(std::string());
    }
    
    // convert to shereical coordinates and redshift
    sph_point = tmp_point;
    z = cosmo.invCoorDist(sph_point.r);
    
    std::cout << tmp_point << std::endl << sph_point << std::endl;
    
    //***** read in data from light cone file
    LensHaloNFW * halo = new LensHaloNFW(mass,Rvir,z,Rvir/Rscale,1,0,0);
    halo->extendRadiius(2);
    lensVec.push_back( halo );
    lensVec.back()->setTheta(sph_point.phi,sph_point.theta);
  }
  
  file.close();
}

void LightCone::ReadBoxRockStar(
                  std::string filename    /// input file name
                 ,Point_3d xo     /// observers location within the box in comoving Mpc
                 ,Point_3d V      /// direction of view, doesn't need to be normalized
                 ,double rlow     /// lowest radial distance accepted
                 ,double rhigh    /// highest radial distance accepted
                 ,std::vector<DataRockStar> &conehalos  /// output list of halos, added to and not cleared
                 ,bool periodic_boundaries  /// if false the cone will intersect the box at most once, otherwise when the cone can be extended through an unending series of repeated boxes
                 ,bool allow_subhalos   /// if false only host halos are used
){
  
  std::cout <<" Opening " << filename << std::endl;
  std::ifstream file(filename.c_str());
  if(!file){
    std::cout << "Can't open file " << filename << std::endl;
    ERROR_MESSAGE();
    throw std::runtime_error(" Cannot open file.");
  }
  conehalos.clear();
  
  std::string myline;
  
  const size_t blocksize = 100000;
  std::vector<DataRockStar> boxhalos(1);
  boxhalos.reserve(blocksize);
  
  //const int ncolumns = 11;
  const int ncolumns = 42;
  
  void *addr[ncolumns];
  DataRockStar halo;
  double tmp;
  double parent_id;
  
  addr[0] = &(halo.id);
  addr[1] = &tmp;
  addr[2] = &(halo.mass);
  addr[3] = &tmp;
  addr[4] = &tmp;
  addr[5] = &(halo.Rvir);
  addr[6] = &(halo.Rscale);
  addr[7] = &tmp;
  addr[8] = &(halo.x[0]);
  addr[9] = &(halo.x[1]);
  addr[10] = &(halo.x[2]);
  for(int i=11 ; i < 41 ;++i) addr[i] = &tmp;
  addr[41] = &parent_id;
  
  unsigned int mysize_t;
  int myint;
  std::string strg;
  std::string delim=" ";
  double mydouble;
  double h,scale_factor;
  double BoxLength =0;
  
  std::stringstream buffer;
  int i_block = 0;
  
  do{
    boxhalos.clear();
    
    while ( boxhalos.size() < blocksize && getline(file,myline)) {
      if(myline[0] == '#'){
        int pos;
        
        pos = myline.find("a = ");
        if(pos != -1){
          myline.erase(0,pos+4);
          buffer << myline;
          buffer >> scale_factor;
          buffer.clear();
        }
        pos = myline.find("h = ");
        if(pos != -1){
          myline.erase(0,pos+4);
          buffer << myline;
          buffer >> h;
          buffer.clear();
        }
        pos = myline.find("Box size: ");
        if(pos != -1){
          myline.erase(0,pos+10);
          myline.erase(9,1000);
          buffer << myline;
          buffer >> BoxLength;
          buffer.clear();
          BoxLength /= h;
        }
        
        continue;
      }
      
      for(int l=0;l<ncolumns; l++){
        int pos = myline.find(delim);
        strg.assign(myline,0,pos);
        buffer << strg;
        
        std::cout << l << "  " << strg << std::endl;
        if(l == 0){
          buffer >> mysize_t;
          *((unsigned int *)addr[l]) = mysize_t;
        }else if(l == 1 ){
          buffer >> myint;
          *((int *)addr[l]) = myint;
        }else{
          buffer >> mydouble;
          *((double *)addr[l]) = mydouble;
        }
        myline.erase(0,pos+1);
        strg.clear();
        buffer.clear();
        buffer.str(std::string());
      }
      
      halo.x /= h;
      halo.mass /= h;
      halo.Rscale *= scale_factor/h;
      halo.Rvir *= scale_factor/h;
      
      if(parent_id == -1 || allow_subhalos )boxhalos.push_back(halo);
      
      //std::cout << boxhalos.back().id << " " << boxhalos.back().x[1] << std::endl;
    }
    
    select(xo,V,BoxLength,rlow,rhigh,boxhalos,conehalos,periodic_boundaries);
    
    i_block += boxhalos.size();
    
    if(boxhalos.size() > 0) std::cout << i_block << " halos from " << filename << ", total in cone currently:"
    << conehalos.size() << "...." << std::endl;
  }while(boxhalos.size() > 0);
  std::cout << "done" << std::endl;
  
  file.close();
}

void LightCone::WriteLightCone(std::string filename,std::vector<DataRockStar> &vec){
  
  std::ofstream file(filename);
  file << "ID,x,y,z,mass,Rvir,Rscale" << std::endl;
  
  for(auto h : vec){
    file << h.id << "," << h.x[0] << "," << h.x[1] << "," << h.x[2] << "," << h.mass
    << "," << h.Rvir << "," << h.Rscale << std::endl;
  }
  
  file.close();
}



