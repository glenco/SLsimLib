#include <fstream>
#include <mutex>
#include <thread>
#include "slsimlib.h"
#include "particle_types.h"
#include "particle_halo.h"
#include "simpleTree.h"
#include "utilities_slsim.h"
#include "gadget.hh"


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
                   ,bool recenter /// recenter so that the LenHalos are centered on the center of mass
                   ):filename(filename),Nsmooth(Nsmooth)
{
  
  
  
  if(format == glamb ){
    nparticles.resize(6,0);
    readSizesB(filename,data,Nsmooth,nparticles,z_original);
  }else{
    
    std::string sizefile = filename + "_S"
    + std::to_string(Nsmooth) + ".glamb";
    
    if(access( sizefile.c_str(), F_OK ) != -1){
      nparticles.resize(6,0);
      readSizesB(sizefile,data,Nsmooth,nparticles,z_original);
    }else{
      // processed file does not exist read the gadget file and find sizes
      switch (format) {
        case gadget2:
          readGadget2();
          break;
        case csv:
          readCSV();
        default:
          break;
      }
      
      // write size file so that next time we wont have to do this
      writeSizesB(sizefile,data,Nsmooth,nparticles,z_original);
    }
  }
  
  // find center of mass
  if(recenter){
    double m=0;
    Point_3d cm(0,0,0);
    for(auto p : data){
      cm[0] += p[0]*p.Mass;
      cm[1] += p[1]*p.Mass;
      cm[2] += p[2]*p.Mass;
      m += p.Mass;
    }
    cm /= m;
    for(auto &p : data){
      p[0] -= cm[0];
      p[1] -= cm[1];
      p[2] -= cm[2];
    }
  }
}

MakeParticleLenses::MakeParticleLenses(const std::string &filename  /// path / name of galmb file
                   ,bool recenter /// recenter so that the LenHalos are centered on the center of mass
                   ):filename(filename)
{
  
  
  nparticles.resize(6,0);
  readSizesB(filename,data,Nsmooth,nparticles,z_original);
  
  
  // find center of mass
  if(recenter){
    double m=0;
    Point_3d cm(0,0,0);
    for(auto p : data){
      cm[0] += p[0]*p.Mass;
      cm[1] += p[1]*p.Mass;
      cm[2] += p[2]*p.Mass;
      m += p.Mass;
    }
    cm /= m;
    for(auto &p : data){
      p[0] -= cm[0];
      p[1] -= cm[1];
      p[2] -= cm[2];
    }
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

bool MakeParticleLenses::readCSV(){
  
  std::ifstream file(filename);
  std::string line = "";
  size_t ntot = 0;
  //while (getline(file, line) && ntot < 1000) ntot++; // ????
  while (getline(file, line) && line[0] == '#');
  ++ntot;
  while (getline(file, line)) ntot++;

  std::cout << "counted " << ntot << " entries in CSV file "
  << filename << std::endl;
  std::cout << " attempting to read them...";
  data.resize(ntot);
  
  std::string delimiter = ",";
  
  //*** be able to read different types of csv files
  
  ntot = 0;
  file.clear();
  file.seekg(0);  // return to begining
  //while (getline(file, line) && line[0] == '#');

  // Iterate through each line and split the content using delimeter
  while (getline(file, line)){
    std::vector<std::string> vec;
    Utilities::splitstring(line,vec,delimiter);
    
    data[ntot][0] =  stof(vec[0]);
    data[ntot][1] =  stof(vec[1]);
    data[ntot][2] =  stof(vec[2]);
    data[ntot].Mass = stof(vec[3]);
    data[ntot].type = 0;
    ++ntot;
  }
     
  // Close the File
  file.close();
  
  std::cout << ntot << " particle positions read from file " << filename << std::endl;
  
  nparticles = {ntot,0,0,0,0,0};
  masses = {0,0,0,0,0,0};
  LensHaloParticles<ParticleType<float> >::calculate_smoothing(Nsmooth,data.data(),data.size());
  
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
  
  // ????? **** convert to physical Mpc/h and Msun/h
  // ???? *** can we store sizes in gadget blocks
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


