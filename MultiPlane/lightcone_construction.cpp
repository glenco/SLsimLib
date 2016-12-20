//
//  lightcone_construction.cpp
//  GLAMER
//
//  Created by Ben Metcalf on 26/10/16.
//
//

#include <mutex>
#include <thread>

#include <boost/iostreams/device/mapped_file.hpp> // for mmap
#include <boost/iostreams/stream.hpp>             // for stream

#include "particle_halo.h"
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
                                 ,PosType &theta_max   /// largest angular seporation from theta=0,phi=0
){
  
  Utilities::delete_container(lensVec);
  
  std::ifstream file(filename.c_str());
  if(!file){
    std::cout << "Can't open file " << filename << std::endl;
    ERROR_MESSAGE();
    throw std::runtime_error(" Cannot open file.");
  }

  
  // make lookup table for sigm(m)
  Utilities::LogLookUpTable<double> tophat(
      std::bind(&COSMOLOGY::TopHatVarianceM,cosmo,std::placeholders::_1,0)
                                           ,1.0e6,1.0e16,1000);
  
  std::string myline;
  
  while(file.peek() == '#') file.ignore(10000,'\n');
  
  const int ncolumns = 7;
  void *addr[ncolumns];
  size_t haloid;
  double mass,Rvir,Rscale;
  Utilities::Geometry::SphericalPoint sph_point,sph_center(0,0,0);
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
  
  //double mass_range[2],z_range[2],rvir_range[2],theta_range[2];
  //mass_range[0] = z_range[0] = rvir_range[0] = theta_range[0] = HUGE_VALF;
  //mass_range[1] = z_range[1] = rvir_range[1] = theta_range[1] = 0.0;
  
  Utilities::Range<double> mass_range(HUGE_VALF,0),z_range(HUGE_VALF,0),rvir_range(HUGE_VALF,0),theta_range(HUGE_VALF,0);
  
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
    
    //std::cout << tmp_point << std::endl << sph_point << std::endl;
    
    
    double nu;
    tophat(mass,nu);  // look up sigm(M)^2
    nu = cosmo.delta_c()/sqrt(nu)/cosmo.Dgrowth(z);

    //mass *= 0.82*(1+0.63*exp(-nu/3.52));
    
    //***** read in data from light cone file
    LensHaloNFW * halo = new LensHaloNFW(mass,Rvir/1.0e3,z,Rvir/Rscale,1,0,0);
    
    // Formula for splash redius from More,Diemer & Kravtsov, 2014
    //halo->extendRadius(0.81*(1+0.97*exp(-nu/2.44)));
    
    //halo->LensHalo::set_Rsize(Rvir*0.81*(1+0.97*exp(-nu/2.44)));
    
    lensVec.push_back( halo );
    lensVec.back()->setTheta(sph_point.phi,sph_point.theta);
    
    theta_range.update(sph_center.AngleSeporation(sph_point));
    mass_range.update(mass);
    z_range.update(z);
    rvir_range.update(Rvir);
  }
  theta_max = theta_range.max();
  
  file.close();
  std::cout << mass_range.min() << " <= mass <= " << mass_range.max() << std::endl;
  std::cout << z_range.min() << " <= z <= " << z_range.max() << std::endl;
  std::cout << rvir_range.min() << " <= Rvir <= " << rvir_range.max() << std::endl;
}

/** \brief Read in a lightcone data file and render them into LensHaloParticles
*
*   filename file should be the output of LightCone::WriteLightCone() to
*   be in the right format.
*/
using Utilities::Geometry::SphericalPoint;

void LightCone::ReadLightConeParticles(
                                       std::string filename   /// name of file to read
                                       ,COSMOLOGY &cosmo      /// cosmology for changing distance to redshift
                                       ,std::vector<LensHalo* > &lensVec  /// output LensHalos
                                       ,std::vector<float> &particle_masses
                                       ,std::vector<float> &particle_sizes
                                       ,int Nplanes
                                       ){
  bool smooth = true;
  SphericalPoint center;
  
  Utilities::delete_container(lensVec);
  
  std::ifstream file(filename.c_str());
  if(!file){
    std::cout << "Can't open file " << filename << std::endl;
    ERROR_MESSAGE();
    throw std::runtime_error(" Cannot open file.");
  }
  
  std::string myline;
  
  // skip first line
  std::getline(file,myline);
  std::string delim =",";
  size_t i=0;
  
  //boost::iostreams::mapped_file_source mmap(filename);
  // boost::iostreams::stream<boost::iostreams::mapped_file_source> is(mmap, std::ios::binary);
  
  size_t Npoints = 100;
  int Nsmooth = 64;
  
  //double z_range[2],theta_range[2];
  //z_range[0] = theta_range[0] = HUGE_VALF;
  //z_range[1] = theta_range[1] = 0.0;
  
  std::vector<SphericalPoint> particles;
  particles.reserve(1000000);
  Point_3d point;
  Utilities::Range<double> theta_range(HUGE_VALF,-1);

  while (getline(file,myline) && i < 100) {
    
    std::cout << myline << std::endl;
    
    int pos1 = 0;
    int pos = myline.find(delim);
    
    point[0] = std::stod(myline.substr(pos1,pos));
    
    //std::cout << xp[i][0] << std::endl;
    
    pos1 = pos + 1;
    pos = myline.find(delim,pos1);
    
    point[1] = std::stod(myline.substr(pos1,pos));
    
    //std::cout << xp[i][1] << std::endl;
    point[2] =std::stod(myline.substr(pos+1));
    
    std::cout << point << std::endl;
    
    particles.push_back(point);
    theta_range.update(center.AngleSeporation(particles.back()));
    
    ++i;
  }
  
  file.close();
  
  // sort by distance
  std::sort(particles.begin(),particles.end()
            ,[](SphericalPoint &p1,SphericalPoint &p2){return p1.r < p2.r;});
  assert(particles[0].r <= particles.back().r);
  
  std::vector<float> sizes = {0.0};
  //PosType **xp = Utilities::PosTypeMatrix(100,3);

// find smoothing
  //sizefile = simfile + "." + std::to_string(Nsmooth) + "sizes";
/*  if(smooth){
//  if(!readSizesFile(sizefile,Nsmooth,min_size)){
    // calculate sizes
    sizes.resize(Npoints);
    //calculate_smoothing(Nsmooth);
    
    std::cout << "Calculating smoothing of particles ..." << std::endl
    << Nsmooth << " neighbors.  If there are a lot of particles this could take a while." << std::endl;
    
    LensHaloParticles::find_smoothing(xp,Npoints,sizes,Nsmooth);
    
    // save result to a file for future use
    //writeSizes(sizefile,Nsmooth);
  }
*/
  
  std::vector<double> planes(Nplanes);
  std::vector<size_t> pnumbers(Nplanes);
  
  double vol = pi*theta_range.max()*theta_range.max()*(pow(particles.back().r,3) - pow(particles[0].r,3))/3;
  double psize = pow(vol/particles.size(),1.0/3.0);
  PosType **xp;
  
  // break up into slices redshift bins
  for(int i = 0 ; i < Nplanes ; ++i ) planes[i] = (2*i+1)*(particles.back().r - particles[0].r)/Nplanes/2 + particles[0].r;
  std::vector<SphericalPoint>::iterator it1 = particles.begin(),it2;
  for(int i = 0 ; i < Nplanes ; ++i ){
    it2 = std::lower_bound(it1,particles.end(),(i+1)*particles.back().r/Nplanes
                           ,[](SphericalPoint &p,double v){return p.r < v;});
    
    assert(it2-it1 > 0);
    
    pnumbers[i] = it2 - it1;
    
    double z = cosmo.invCoorDist(planes[i]);
    double Dl = planes[i]/(1+z);
    double sigma_back = particle_mass*pnumbers[i]/pi/Dl/Dl/theta_range.max()/theta_range.max();
    
    // project onto planes
    xp = Utilities::PosTypeMatrix(it2-it1,2);
    for(auto it = it1 ; it != it2 + 1 ; ++it){
      xp[i][0] = (*it).theta*Dl;
      xp[i][1] = (*it).phi*Dl;
    }
    // volume of shell
    std::vector<float> sizes(pnumbers[i],psize/(1+z));
    std::vector<float> masses = {particle_mass};
    lensVec.push_back(new LensHaloParticles(xp,sizes,masses,z,cosmo,false,sigma_back));
    
    it1 = it2 + 1;
  }
}


void LightCone::WriteLightCone(std::string filename,std::vector<DataRockStar> &vec){
  
  std::ofstream file(filename);
  file << "ID,x,y,z,mass,Rvir(kpc),Rscale(kpc),Vmax(km/s)" << std::endl;
  
  for(auto h : vec){
    file << h.id << "," << h.x[0] << "," << h.x[1] << "," << h.x[2] << "," << h.mass
    << "," << h.Rvir << "," << h.Rscale << "," << h.Vmax << std::endl;
  }
  
  file.close();
}

void LightCone::WriteLightCone(std::string filename,std::vector<Point_3d> &vec){
  
  std::ofstream file(filename);
  file << "x,y,z" << std::endl;
  
  for(auto h : vec){
    file << h[0] << "," << h[1] << "," << h[2] << std::endl;
  }
  
  file.close();
}

void MultiLightCone::ReadBoxRockStar(std::string filename
                                     ,double rlow,double rhigh
                                     ,std::vector<std::vector<LightCone::DataRockStar> > &conehalos
                                     ,bool periodic_boundaries
                                     ,bool allow_subhalos){
  
  
  std::cout <<" Opening " << filename << std::endl;
  std::ifstream file(filename.c_str());
  if(!file){
    std::cout << "Can't open file " << filename << std::endl;
    ERROR_MESSAGE();
    throw std::runtime_error(" Cannot open file.");
  }
  if(conehalos.size() != cones.size() ){
    if(conehalos.size() ==0 ){
      conehalos.resize(cones.size());
    }else{
      std::cerr << " conhalos vector passed into MultiLightCone::ReadBoxRockStar does does not mathch the size expected." << std::endl;
      throw std::invalid_argument("conehalos wrong size.");
    }
  }
  
  for(auto c : conehalos ) std::cout << c.size() << " halos already in cone." << std::endl;
  
  if(rlow > rhigh) std::swap(rlow,rhigh);
  if(rlow == rhigh) return;
  if(rlow < 0) rlow = 0.0;
  
  std::string myline;
  
  const size_t blocksize = 500000;
  std::vector<LightCone::DataRockStar> boxhalos;
  boxhalos.reserve(blocksize);
  
  //const int ncolumns = 11;
  const int ncolumns = 42;
  
  void *addr[ncolumns];
  LightCone::DataRockStar halo;
  double tmp;
  double parent_id;
  double m200c,m200b,R200,Mvir_all,rs_klypin;
  for(int i=0 ; i < ncolumns ;++i) addr[i] = &tmp;
  
  addr[0] = &(halo.id);
  
  addr[2] = &(halo.mass);
  addr[3] = &(halo.Vmax);

  addr[5] = &(halo.Rvir);
  addr[6] = &(halo.Rscale);

  addr[8] = &(halo.x[0]);
  addr[9] = &(halo.x[1]);
  addr[10] = &(halo.x[2]);
  
  addr[18] = &rs_klypin;
  addr[19] = &Mvir_all;
  addr[20] = &m200b;
  addr[21] = &m200c;

  addr[41] = &parent_id;

  unsigned int mysize_t;
  int myint;
  std::string strg;
  std::string delim=" ";
  double mydouble;
  double h,scale_factor,Omega;
  double BoxLength =0;
  double totalmass = 0.0;
  
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
        pos = myline.find("Om = ");
        if(pos != -1){
          myline.erase(0,pos+5);
          buffer << myline.substr(0,6);
          buffer >> Omega;
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
        
        //std::cout << l << "  " << strg << std::endl;
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
      
      //std::cout << halo.id << " " << halo.x << std::endl;
      
      halo.x /= h;
      halo.mass /= h;
      halo.Rscale *= scale_factor/h;
      halo.Rvir *= scale_factor/h;
      
      totalmass += halo.mass;
      
      if(parent_id == -1 || allow_subhalos )boxhalos.push_back(halo);
      
    }
    
    if(boxhalos.size() > 0){
      int nthreads = Utilities::GetNThreads();
      int chunk_size;
      do{
        chunk_size =  boxhalos.size()/nthreads;
        if(chunk_size == 0) nthreads /= 2;
      }while(chunk_size == 0);
      
      int remainder =  boxhalos.size()%chunk_size;
      
      
      for(int j=0; j< cones.size() ; ++j ){
        
        std::vector<std::thread> thr(nthreads);
        Utilities::LockableContainer<std::vector<LightCone::DataRockStar> > halos;
        halos.swap(conehalos[j]);
        
        assert(nthreads*chunk_size + remainder == boxhalos.size() );
        for(int ii =0; ii< nthreads ; ++ii){
          
          //std::cout << ii*chunk_size << " " << n << std::endl;
          
          thr[ii] = std::thread(&LightCone::select<LightCone::DataRockStar>,cones[j]
                                ,xos[j],vs[j],BoxLength,rlow,rhigh
                                ,boxhalos.data() + ii*chunk_size
                                ,boxhalos.data() + (ii+1)*chunk_size + (ii==nthreads-1)*remainder
                                //,std::ref(conehalos[j])
                                ,std::ref( halos )
                                ,periodic_boundaries );
        }
        for(int ii = 0; ii < nthreads ;++ii){ if(thr[ii].joinable() ) thr[ii].join();}
        halos.swap(conehalos[j]);
      }
      
      
      /*for(int i =0 ; i < cones.size() ;++i){
       cones[i].select<LightCone::DataRockStar>(xos[i],vs[i],BoxLength,rlow,rhigh
       ,boxhalos.data()
       ,boxhalos.data() + boxhalos.size()
       ,conehalos[i],periodic_boundaries);
       }*/
      i_block += boxhalos.size();
      
      if(boxhalos.size() > 0 && i_block % 10000000 == 0 ){
        std::cout << i_block << " halos from " << filename << ", total in cones currently: "
        << std::endl;
        for(auto c : conehalos) std::cout << "    " << c.size() << "...." << std::endl;
      }
    }
  }while(boxhalos.size() > 0);
  std::cout << "done" << std::endl;
  std::cout << "Total mass in halos: " << totalmass << " Msun." << std::endl
  << "Mass density in halos: " << totalmass/BoxLength/BoxLength/BoxLength
  << " Msun/Mpc^3 comoving." << std::endl;
  
  
  
  file.close();
}

void MultiLightCone::ReadBoxXYZ(std::string filename
                                ,double rlow,double rhigh
                                ,std::vector<std::vector<Point_3d> > &conehalos
                                ,double hubble
                                ,double BoxLength
                                ,bool periodic_boundaries
                                ){
  
  
  std::cout <<" Opening " << filename << std::endl;
  std::ifstream file(filename.c_str());
  if(!file){
    std::cout << "Can't open file " << filename << std::endl;
    ERROR_MESSAGE();
    throw std::runtime_error(" Cannot open file.");
  }
  if(conehalos.size() != cones.size() ){
    if(conehalos.size() ==0 ){
      conehalos.resize(cones.size());
    }else{
      std::cerr << " conhalos vector passed into MultiLightCone::ReadBoxXYZ does does not mathch the size expected." << std::endl;
      throw std::invalid_argument("conehalos wrong size.");
    }
  }
  
  for(auto c : conehalos ) std::cout << c.size() << " halos already in cone." << std::endl;
  
  if(rlow > rhigh) std::swap(rlow,rhigh);
  if(rlow == rhigh) return;
  if(rlow < 0) rlow = 0.0;
  
  std::string myline;
  
  const size_t blocksize = 500000;
  std::vector<Point_3d> boxhalos;
  boxhalos.reserve(blocksize);
  
  const int ncolumns = 3;
  
  Point_3d halo;
  std::string strg;
  std::string delim = " ";
  
  std::stringstream buffer;
  int i_block = 0;
  
  do{
    boxhalos.clear();
    
    while ( boxhalos.size() < blocksize
           && getline(file,myline)) {
      
      int pos = myline.find_first_not_of(delim);
      myline.erase(0,pos);
      
      for(int l=0;l<ncolumns; l++){
        pos = myline.find(delim);
        strg.assign(myline,0,pos);
        buffer << strg;
        
        //std::cout << l << "  " << strg << std::endl;
        buffer >> halo[l];
        //std::cout << halo << std::endl;
        
        myline.erase(0,pos+1);
        pos = myline.find_first_not_of(delim);
        myline.erase(0,pos);
        
        strg.clear();
        buffer.clear();
        buffer.str(std::string());
      }
      
      halo /= hubble;
      
      boxhalos.push_back(halo);
    }
    
    if(boxhalos.size() > 0){
      int nthreads = Utilities::GetNThreads();
      int chunk_size;
      do{
        chunk_size =  boxhalos.size()/nthreads;
        if(chunk_size == 0) nthreads /= 2;
      }while(chunk_size == 0);
      
      int remainder =  boxhalos.size()%chunk_size;
      
      
      for(int j=0; j< cones.size() ; ++j ){
        
        std::vector<std::thread> thr(nthreads);
        Utilities::LockableContainer<std::vector<Point_3d> > halos;
        halos.swap(conehalos[j]);
        
        assert(nthreads*chunk_size + remainder == boxhalos.size() );
        for(int ii =0; ii< nthreads ; ++ii){
          
          //std::cout << ii*chunk_size << " " << n << std::endl;
          
          thr[ii] = std::thread(&LightCone::select<Point_3d>,cones[j]
                                ,xos[j],vs[j],BoxLength,rlow,rhigh
                                ,boxhalos.data() + ii*chunk_size
                                ,boxhalos.data() + (ii+1)*chunk_size + (ii==nthreads-1)*remainder
                                //,std::ref(conehalos[j])
                                ,std::ref( halos )
                                ,periodic_boundaries );
        }
        for(int ii = 0; ii < nthreads ;++ii){ if(thr[ii].joinable() ) thr[ii].join();}
        halos.swap(conehalos[j]);
      }
      
      i_block += boxhalos.size();
      
      if(boxhalos.size() > 0 && i_block % 10000000 == 0 ){
        std::cout << i_block << " halos from " << filename << ", total in cones currently: "
        << std::endl;
        for(auto c : conehalos) std::cout << "    " << c.size() << "...." << std::endl;
      }
    }
  }while(boxhalos.size() > 0);
  
  std::cout << "done" << std::endl;
  
  file.close();
}


