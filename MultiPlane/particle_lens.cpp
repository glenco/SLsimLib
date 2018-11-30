#include <fstream>
#include <mutex>
#include <thread>
#include "slsimlib.h"
#include "particle_types.h"
#include "particle_halo.h"
#include "simpleTree.h"
#include "utilities_slsim.h"
#include "gadget.hh"

#ifdef ENABLE_HDF5
#include "H5Cpp.h"
#endif

#ifdef ENABLE_FITS
#include <CCfits/CCfits>
using namespace CCfits;
#endif

// *************************************************************************************
// ******************** methods for MakeParticleLenses *********************************
// *************************************************************************************
// *************************************************************************************

MakeParticleLenses::MakeParticleLenses(
                    const std::string &filename  /// path / root name of gadget-2 snapshot
                   ,SimFileFormat format
                   ,int Nsmooth   /// number of nearest neighbors used for smoothing
                   ,bool recenter /// recenter so that the LenHalos are centered on the center of mass of all the particles
                   ):filename(filename),Nsmooth(Nsmooth)
{
  
  
  
  if(format == glmb ){
    nparticles.resize(6,0);
    readSizesB(filename,data,Nsmooth,nparticles,z_original);
  }else{
    
    std::string sizefile = filename + "_S"
    + std::to_string(Nsmooth) + ".glmb";
    
    if(access( sizefile.c_str(), F_OK ) != -1){
      nparticles.resize(6,0);
      readSizesB(sizefile,data,Nsmooth,nparticles,z_original);
    }else{
      // processed file does not exist read the gadget file and find sizes
      switch (format) {
        case gadget2:
          readGadget2();
          break;
        case csv3:
          readCSV(3);
          break;
        case csv4:
          readCSV(4);
          break;
        case csv5:
          readCSV(5);
          break;
        case csv6:
          readCSV(6);
        default:
          break;
      }
      
      // write size file so that next time we wont have to do this
      writeSizesB(sizefile,data,Nsmooth,nparticles,z_original);
    }
  }
  
  // find center of mass
  cm[0]=cm[1]=cm[2]=0;
  double m=0;
  for(auto p : data){
    for(int i=0 ; i < 3 ; ++i){
      cm[i] += p[i]*p.Mass;
    
      if(bbox_ll[i] > p[i] ) bbox_ll[i] = p[i];
      if(bbox_ur[i] < p[i] ) bbox_ur[i] = p[i];
    }
    m += p.Mass;
  }
  cm /= m;

  // subtract center of mass
  if(recenter){
     for(auto &p : data){
      p[0] -= cm[0];
      p[1] -= cm[1];
      p[2] -= cm[2];
    }
    
    for(int i=0 ; i < 3 ; ++i){
      bbox_ll[i] -= cm[i];
      bbox_ur[i] -= cm[i];
    }
    cm[0]=cm[1]=cm[2]=0;
  }
}

MakeParticleLenses::MakeParticleLenses(const std::string &filename  /// path / name of glmb file
                   ,bool recenter /// recenter so that the LenHalos are centered on the center of mass
                   ):filename(filename)
{
  
  
  nparticles.resize(6,0);
  readSizesB(filename,data,Nsmooth,nparticles,z_original);
  
  // find center of mass
  cm[0]=cm[1]=cm[2]=0;
  double m=0;
  for(auto p : data){
    for(int i=0 ; i < 3 ; ++i){
      cm[i] += p[i]*p.Mass;
      
      if(bbox_ll[i] > p[i] ) bbox_ll[i] = p[i];
      if(bbox_ur[i] < p[i] ) bbox_ur[i] = p[i];
    }
    m += p.Mass;
  }
  cm /= m;
  
  // subtract center of mass
  if(recenter){
    for(auto &p : data){
      p[0] -= cm[0];
      p[1] -= cm[1];
      p[2] -= cm[2];
    }
    
    for(int i=0 ; i < 3 ; ++i){
      bbox_ll[i] -= cm[i];
      bbox_ur[i] -= cm[i];
    }
    cm[0]=cm[1]=cm[2]=0;
  }
}

void MakeParticleLenses::CreateHalos(const COSMOLOGY &cosmo,double redshift){

  // put into proper units
  float h = cosmo.gethubble();
  
  double length_unit = h;
  double mass_unit = h;

  for(auto &p : data){
    p[0] *= length_unit;
    p[1] *= length_unit;
    p[2] *= length_unit;
    p.Size *= length_unit;
    p.Mass *= mass_unit;
  }
  
  // remove halos if they already exist
  while(halos.size() > 0){
    delete halos.back();
    halos.pop_back();
  }
  
  // create halos
  ParticleType<float> *pp;
  size_t skip = 0;
  Point_2d theta_rotate;
  for(int i = 0 ; i < 6 ; ++i){  //loop through types
    if(nparticles[i] > 0){
      
      pp = data.data() + skip;  // pointer to first particle of type
      /* halos.emplace_back(pp
       ,nparticles[i]
       ,z_original
       ,cosmo
       ,theta_rotate
       ,false
       ,0);*/
      halos.push_back(new LensHaloParticles<ParticleType<float> >(
                                                                  pp
                                                                  ,nparticles[i]
                                                                  ,redshift
                                                                  ,cosmo
                                                                  ,theta_rotate
                                                                  ,false
                                                                  ,0)
                      );
      /*/
       /*LensHaloParticles<ParticleType<float> > halo(pp
       ,nparticles[i]
       ,z_original
       ,cosmo
       ,theta_rotate
       ,false
       ,0);
       */
    }
    skip += nparticles[i];
  }
};

bool MakeParticleLenses::readCSV(int columns_used){
  
  std::string delimiter = ",";
  
  
  std::ifstream file(filename);
  std::string line = "";
  size_t ntot = 0;
  while (getline(file, line) && line[0] == '#');
  std::vector<std::string> vec;
  Utilities::splitstring(line,vec,delimiter);
  const int ncolumns = vec.size();
  
  if(ncolumns < columns_used){
    std::cerr << "file " << filename <<" does not have enough columns." << std::endl;
    throw std::invalid_argument("bad file");
  }
  
  ++ntot;
  while (getline(file, line)) ntot++;

  std::cout << "counted " << ntot << " entries in CSV file "
  << filename << std::endl;
  std::cout << " attempting to read them...";
  data.resize(ntot);
  
  //*** be able to read different types of csv files
  
  ntot = 0;
  file.clear();
  file.seekg(0);  // return to begining
  while (getline(file, line) && line[0] == '#');
  nparticles = {0,0,0,0,0,0};

  // Iterate through each line and split the content using delimeter
  if(columns_used == 3){
    do{
      std::vector<std::string> vec;
      Utilities::splitstring(line,vec,delimiter);
    
      data[ntot][0] =  stof(vec[0]);
      data[ntot][1] =  stof(vec[1]);
      data[ntot][2] =  stof(vec[2]);
      data[ntot].Mass = 1.0;
      ++ntot;
    }while( getline(file, line) );
    nparticles = {ntot,0,0,0,0,0};
    masses = {0,0,0,0,0,0};
  }else if(columns_used == 4){
    do{
      std::vector<std::string> vec;
      Utilities::splitstring(line,vec,delimiter);
      
      data[ntot][0] =  stof(vec[0]);
      data[ntot][1] =  stof(vec[1]);
      data[ntot][2] =  stof(vec[2]);
      data[ntot].Mass = stof(vec[3]);
      data[ntot].type = 0;
      ++ntot;
    }while( getline(file, line) );
    nparticles = {ntot,0,0,0,0,0};
    masses = {0,0,0,0,0,0};
  }else if(columns_used == 5){
    std::cout << "Using the particle sizes from " << filename << std::endl;
    do{
      std::vector<std::string> vec;
      Utilities::splitstring(line,vec,delimiter);
    
      data[ntot][0] =  stof(vec[0]);
      data[ntot][1] =  stof(vec[1]);
      data[ntot][2] =  stof(vec[2]);
      data[ntot].Mass = stof(vec[3]);
      data[ntot].Size = stof(vec[4]);
      data[ntot].type = 0;
      ++ntot;
    }while( getline(file, line) );
    nparticles = {ntot,0,0,0,0,0};
    masses = {0,0,0,0,0,0};
  }else if(columns_used == 6){
    std::cout << "Using the particle sizes from " << filename << std::endl;
    std::cout << "Using the particle type information from " << filename << std::endl;

    do{
      std::vector<std::string> vec;
      Utilities::splitstring(line,vec,delimiter);
      
      data[ntot][0] =  stof(vec[0]);
      data[ntot][1] =  stof(vec[1]);
      data[ntot][2] =  stof(vec[2]);
      data[ntot].Mass = stof(vec[3]);
      data[ntot].Size = stof(vec[4]);
      data[ntot].type = stoi(vec[5]);
      ++ntot;
      ++nparticles[data[ntot].type];
    }while( getline(file, line) );
    
    int ntypes = 0;
    for(size_t n : nparticles){
      if(n > 0) ++ntypes;
    }
    
    if(ntypes > 1){
      // sort by type
      std::sort(data.begin(),data.end(),[](ParticleType<float> &a1,ParticleType<float> &a2){return a1.type < a2.type;});
    }
    masses = {0,0,0,0,0,0};
  }

  // Close the File
  file.close();
  
  std::cout << ntot << " particle positions read from file " << filename << std::endl;
  
  if(columns_used < 5){
    ParticleType<float> *pp;
    size_t skip = 0;
    for(int i = 0 ; i < 6 ; ++i){  //loop through types
      if(nparticles[i] > 0){
        pp = data.data() + skip;  // pointer to first particle of type
        size_t N = nparticles[i];
      
        LensHaloParticles<ParticleType<float> >::calculate_smoothing(Nsmooth,pp,N);
      }
      skip += nparticles[i];
    }
  }
  
  //LensHaloParticles<ParticleType<float> >::calculate_smoothing(Nsmooth,data.data(),data.size());
  
  return true;
};



bool MakeParticleLenses::readGadget2(){
  
  GadgetFile<ParticleType<float> > gadget_file(filename,data);
   z_original = gadget_file.redshift;
   
   for(int n=0 ; n < gadget_file.numfiles ; ++n){
     gadget_file.openFile();
     gadget_file.readBlock("POS");
     gadget_file.readBlock("MASS");
     gadget_file.closeFile();
   }
  
  double a = 1.0e-3/(1 + z_original);
  // **** convert to physical Mpc/h and Msun/h
  for(auto &p : data){
    p[0] *= a;
    p[1] *= a;
    p[2] *= a;
    p.Mass *= 1.0e10;
  }
  
  // sort by type
   std::sort(data.begin(),data.end(),[](ParticleType<float> &a1,ParticleType<float> &a2){return a1.type < a2.type;});
   
   ParticleType<float> *pp;
   size_t skip = 0;
   for(int i = 0 ; i < 6 ; ++i){  //loop through types
   
     nparticles.push_back(gadget_file.npart[i]);
     masses.push_back(gadget_file.masstab[i]);
   
     if(gadget_file.npart[i] > 0){
       pp = data.data() + skip;  // pointer to first particle of type
       size_t N = gadget_file.npart[i];
   
       LensHaloParticles<ParticleType<float> >::calculate_smoothing(Nsmooth,pp,N);
     }
     skip += gadget_file.npart[i];
   }
  
  return true;
};

/*
#ifdef ENABLE_HDF5
bool MakeParticleLenses::readHDF5(){
  
  H5::H5File file(filename.c_str(), H5F_ACC_RDONLY );
  
  std::vector<std::string> sets = {"MASS","TYPE"};
  
  for(auto set : sets){
  
    H5::DataSet dataset = file.openDataSet(set.c_str());
    H5T_class_t type_class = dataset.getTypeClass();

    // Get class of datatype and print message if it's an integer.
      if( type_class == H5T_INTEGER )
      {
        cout << "Data set has INTEGER type" << endl;
        //Get the integer datatype
        H5::IntType intype = dataset.getIntType();
        // Get order of datatype and print message if it's a little endian.
        H5std_string order_string;
        H5T_order_t order = intype.getOrder( order_string );
        cout << order_string << endl;
      
        // Get size of the data element stored in file and print it.
        size_t size = intype.getSize();
        cout << "Data size is " << size << endl;
      }else if(type_class == H5T_FLOAT ){
        cout << "Data set has FLOAT type" << endl;
        //Get the integer datatype
        H5::FloatType intype = dataset.getFloatType();
        // Get order of datatype and print message if it's a little endian.
        H5std_string order_string;
        H5T_order_t order = intype.getOrder( order_string );
        cout << order_string << endl;
      
        // Get size of the data element stored in file and print it.
        size_t size = intype.getSize();
        cout << "Data size is " << size << endl;
      }
  
    // Get dataspace of the dataset.
   
    H5::DataSpace dataspace = dataset.getSpace();
  
   // Get the number of dimensions in the dataspace.
    int rank = dataspace.getSimpleExtentNdims();
  }
  return true;
};

#endif
*/
// remove particles that are beyond radius (Mpc/h) of center
void MakeParticleLenses::radialCut(Point_2d center,double radius){
  
  double radius2 = radius*radius;
  double r2;
  auto end = data.end();
  auto it = data.begin();
  size_t ntot = data.size();
  
  while(it != end){
    r2 = ( (*it)[0]-center[0] )*( (*it)[0]-center[0] )
    + ( (*it)[1]-center[1] )*( (*it)[1]-center[1] )
    +( (*it)[2]-center[2] )*( (*it)[2]-center[2] );
    if(r2 > radius2){
      --nparticles[(*it).type];
      --end;
      --ntot;
      std::swap(*it,*end);
    }else ++it;
  }
  
  data.resize(ntot);
  /// resort by type
  int ntypes = 0;
  for(size_t n : nparticles){
    if(n > 0) ++ntypes;
  }
  
  if(ntypes > 1){
    // sort by type
    std::sort(data.begin(),data.end(),[](ParticleType<float> &a1,ParticleType<float> &a2){return a1.type < a2.type;});
  }
}

// remove particles that are beyond cylindrical radius (Mpc/h) of center
void MakeParticleLenses::cylindricalCut(Point_2d center,double radius){
  
  double radius2 = radius*radius;
  double r2;
  auto end = data.end();
  auto it = data.begin();
  size_t ntot = data.size();
  
  while(it != end){
    r2 = ( (*it)[0]-center[0] )*( (*it)[0]-center[0] )
    + ( (*it)[1]-center[1] )*( (*it)[1]-center[1] );
    if(r2 > radius2){
      --nparticles[(*it).type];
      --end;
      --ntot;
      std::swap(*it,*end);
    }else ++it;
  }
  
  data.resize(ntot);
  /// resort by type
  int ntypes = 0;
  for(size_t n : nparticles){
    if(n > 0) ++ntypes;
  }
  
  if(ntypes > 1){
    // sort by type
    std::sort(data.begin(),data.end(),[](ParticleType<float> &a1,ParticleType<float> &a2){return a1.type < a2.type;});
  }
}

