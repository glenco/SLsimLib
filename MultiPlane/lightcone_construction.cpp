//
//  lightcone_construction.cpp
//  GLAMER
//
//  Created by Ben Metcalf on 26/10/16.
//
//

//#include <mutex>
//#include <thread>

#include <boost/iostreams/device/mapped_file.hpp> // for mmap
#include <boost/iostreams/stream.hpp>             // for stream
#include <deque>

#include "particle_halo.h"
#include "MOKAlens.h"
#include "lightcone_construction.h"
#include "utilities_slsim.h"
#include "profiles.h"

namespace LightCones{
  
  
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
  void ReadLightConeNFW(
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
  
  /** \brief Reads particle data in from a csv file and puts them on multiple planes using tree force solver.
   
   The data file must be in csv format with a header.  The first three columns must be labeled
   "x,y,z". These coordinates should be in comoving Mpc (no h factor).  The file can have additional columns "mass" and "size" to give each particle a different size or mass.  If both are present "mass" must come after "z".  There can be no other columns.  "size" values should be in physical Mpc.
   
   If there is no "size" size column the input "particle_size" will be assigned to all particles.  If "angulare_size" is set to true, "particle_size" is interpreted as an
   angular size and thus will be a different physical size on different planes.  It is
   ignored when a "size" column is present.
   */
  void ReadLightConeParticles(
                              std::string filename   /// name of file to read
                              ,COSMOLOGY &cosmo      /// cosmology for changing distance to redshift
                              ,std::vector<LensHaloParticles *> &lensVec  /// output LensHalos
                              ,int Nplanes           /// number of planes to create
                              ,float particle_mass   /// these will be used only if they are not in the input file, solar masses
                              ,float particle_size   /// these will be used only if they are not in the input file, physical Mpc or radians see angular_size parameter.
                              ,bool angular_sizes  /// if yes partical_size will be interpreted as an angular size in radians
                              ,bool verbose
                              ){
    std::ifstream file(filename.c_str());
    if(!file){
      std::cerr << "Can't open file " << filename << std::endl;
      ERROR_MESSAGE();
      throw std::runtime_error(" Cannot open file.");
    }
    
    SphericalPoint center;
    
    bool multimass;
    bool multisize;
    
    std::deque<float> particle_masses;
    std::deque<float> particle_sizes;
    
    Utilities::delete_container(lensVec);
    //lensVec.clear();
    
    std::string myline;
    
    // skip first line
    std::string delim =",";
    size_t i=0;
    
    //boost::iostreams::mapped_file_source mmap(filename);
    // boost::iostreams::stream<boost::iostreams::mapped_file_source> is(mmap, std::ios::binary);
    
    //double z_range[2],theta_range[2];
    //z_range[0] = theta_range[0] = HUGE_VALF;
    //z_range[1] = theta_range[1] = 0.0;
    
    std::deque<SphericalPoint> particles;
    //particles.reserve(1000000);
    Point_3d point;
    Utilities::Range<double> theta_range(HUGE_VALF,-1);
    
    /// header line
    getline(file,myline);
    if(verbose) std::cout << myline << std::endl;
    if(myline[0] != 'x'){
      std::cerr << "File " << filename << " is in an unexpected format." << std::endl;
      throw std::runtime_error(" Cannot open file.");
    }
    
    
    int pos = myline.find("mass");
    if(pos == std::string::npos){
      particle_masses.push_back(particle_mass);
      multimass = false;
    }else{
      multimass = true;
    }
    
    int pos1 = pos;
    pos = myline.find("size");
    if(pos1 > pos && pos1 != std::string::npos){
      std::cout << myline << std::endl;
      std::cerr << "File " << filename << " is in an unexpected format. Particle mass must come before particles size."<< std::endl;
      throw std::runtime_error(" Cannot open file.");
    }
    if(pos == myline.npos){
      particle_sizes.push_back(particle_size);
      multisize = false;
    }else{
      multisize = true;
    }
    
    while (getline(file,myline) ) {
      
      if(verbose) std::cout << myline << std::endl;
      
      int pos1 = 0;
      int pos = myline.find(delim);
      
      point[0] = std::stod(myline.substr(pos1,pos));
      
      //std::cout << xp[i][0] << std::endl;
      
      pos1 = pos + 1;
      pos = myline.find(delim,pos1);
      
      point[1] = std::stod(myline.substr(pos1,pos));
      
      pos1 = pos + 1;
      pos = myline.find(delim,pos1);
      //std::cout << xp[i][1] << std::endl;
      point[2] =std::stod(myline.substr(pos1,pos));
      //std::cout << point << std::endl;
      
      particles.push_back(point);
      theta_range.update(center.AngleSeporation(particles.back()));
      
      //std::cout << myline.length() << std::endl;
      
      if(multimass){
        if(pos != std::string::npos){
          pos1 = pos + 1;
          pos = myline.find(delim,pos1);
          //std::cout << xp[i][1] << std::endl;
          particle_masses.push_back(std::stof(myline.substr(pos1,pos)));
        }else{
          particle_masses.push_back(0.0);
        }
        std::cout << particle_masses.back() << std::endl;
      }
      if(multisize){
        if(pos != std::string::npos){
          pos1 = pos + 1;
          pos = myline.find(delim,pos1);
          //std::cout << xp[i][1] << std::endl;
          particle_sizes.push_back(std::stof(myline.substr(pos1,pos)));
        }else{
          particle_sizes.push_back(0.0);
        }
        std::cout << particle_sizes.back() << std::endl;
      }
      
      ++i;
    }
    
    size_t Npoints = particles.size();
    
    file.close();
    
    // sort by distance
    if(!multisize && !multimass){
      std::sort(particles.begin(),particles.end()
                ,[](SphericalPoint &p1,SphericalPoint &p2){return p1.r < p2.r;});
      
    }else{
      std::vector<size_t> index(particles.size());
      
      std::sort(index.begin(),index.end()
                ,[&particles](size_t i,size_t j){return particles[i].r < particles[j].r;});
      
      {
        Utilities::apply_permutation_in_place(particles,index);
      }
      
      assert(particles[0].r < particles.back().r);
      
      if(multimass){
        assert(Npoints == particle_masses.size());
        Utilities::apply_permutation_in_place(particle_masses,index);
      }
      
      if(multisize){
        Utilities::apply_permutation_in_place(particle_sizes,index);
      }
      
    }
    
    assert(particles[0].r <= particles.back().r);
    std::cout << particles.size() << " particles read in from cone."
    << std::endl;
    
    std::vector<double> D_planes(Nplanes);
    
    //double vol = pi*theta_range.max()*theta_range.max()*(pow(particles.back().r,3) - pow(particles[0].r,3))/3;
    //double psize = pow(vol/particles.size(),1.0/3.0);
    
    // break up into slices redshift bins
    for(int i = 0 ; i < Nplanes ; ++i ) D_planes[i] = (2*i+1)*(particles.back().r - particles[0].r)/Nplanes/2 + particles[0].r;
    
    //Utilities::RandomNumbers_NR ran(1234);
    PosType **xp;
    std::deque<SphericalPoint>::iterator it2;
    size_t i1 = 0;
    for(int i = 0 ; i < Nplanes ; ++i ){
      it2 = std::lower_bound(particles.begin(),particles.end()
                             ,(i+1)*particles.back().r*1.001/Nplanes
                             ,[](SphericalPoint &p,double v){return p.r < v;});
      
      double z = cosmo.invCoorDist(D_planes[i]);
      double Dl = D_planes[i]/(1+z);
      double total_mass = 0.0;
      
      // project onto planes
      size_t Np_on_plane = it2 - particles.begin();
      
      if(verbose) std::cout << "Number of particles on plane " << i << " : " << Np_on_plane
        << " z = " << z << std::endl;
      
      xp = Utilities::PosTypeMatrix(Np_on_plane,2);
      size_t j = i1;
      for(size_t ii = 0 ; ii < Np_on_plane ; ++ii,++j){
        //for(auto it = it1 ; it != it2 + 1 ; ++it){
        //xp[i][0] = (*it).theta*Dl;
        //xp[i][1] = (*it).phi*Dl;
        xp[ii][0] = particles.front().theta*Dl;
        xp[ii][1] = particles.front().phi*Dl;
        
        /// test lines for randomizing positions
        //double tmp = theta_range.max()*sqrt(ran()),tmp_t = 2*pi*ran();
        //xp[i][0] = tmp*cos(tmp_t)*Dl;
        //xp[i][1] = tmp*sin(tmp_t)*Dl;
        
        total_mass += particle_masses[j*multimass];
        particles.pop_front();
      }
      
      double sigma_back = total_mass/pi/Dl/Dl/theta_range.max()/theta_range.max();
      
      // transfer the sizes and masses to a vector
      
      std::vector<float> sizes(Np_on_plane*multisize + !multisize);
      if(multisize){
        for(size_t i=0 ; i < Np_on_plane ; ++i){
          sizes[i] = particle_sizes.front();
          particle_sizes.pop_front();
        }
      }else{
        sizes[0] = particle_size;
        if(angular_sizes) sizes[0] *= Dl;
      }
      
      std::vector<float> masses(Np_on_plane*multimass + !multimass);
      if(multisize){
        for(size_t i=0 ; i < Np_on_plane ; ++i){
          masses[i] = particle_masses.front();
          particle_masses.pop_front();
        }
      }else{
        masses[0] = particle_mass;
      }
      
      lensVec.push_back(new LensHaloParticles(xp,Np_on_plane,sizes,masses,z,cosmo,sigma_back));
      
      i1 += Np_on_plane;
    }
    
    assert(i1 == Npoints);
  }
  
  /** \brief Reads particle data in from a csv file and puts them on multiple planes using FFT force solver.
   
   The data file must be in csv format with a header.  The first three columns must be labeled
   "x,y,z". These coordinates should be in comoving Mpc (no h factor).  The file can have additional columns "mass" and "size" to give each particle a different mass.  In the current implementation "size" is not used.  If both are present "mass" must come after "z".  There can be no other columns.
   
   If "angulare_resolution" is the angular resolution of the mass maps which, in the current implementation, is the same on all planes.
   */
  void ReadLightConeParticles(
                              std::string filename   /// name of file to read
                              ,COSMOLOGY &cosmo      /// cosmology for changing distance to redshift
                              ,std::vector<LensHaloMassMap *> &lensVec  /// output LensHalos
                              ,int Nplanes           /// number of planes to create
                              ,float particle_mass   /// these will be used only if they are not in the input file, solar masses
                              ,float angular_resolution  /// angular resolution of mass maps in redians
                              ,bool verbose
                              ){
    std::ifstream file(filename.c_str());
    if(!file){
      std::cerr << "Can't open file " << filename << std::endl;
      ERROR_MESSAGE();
      throw std::runtime_error(" Cannot open file.");
    }
    
    SphericalPoint center;
    
    bool multimass;
    
    std::deque<float> particle_masses;
    
    Utilities::delete_container(lensVec);
    //lensVec.clear();
    
    std::string myline;
    
    // skip first line
    std::string delim =",";
    size_t i=0;
    
    //boost::iostreams::mapped_file_source mmap(filename);
    // boost::iostreams::stream<boost::iostreams::mapped_file_source> is(mmap, std::ios::binary);
    
    //double z_range[2],theta_range[2];
    //z_range[0] = theta_range[0] = HUGE_VALF;
    //z_range[1] = theta_range[1] = 0.0;
    
    std::deque<SphericalPoint> particles;
    //particles.reserve(1000000);
    Point_3d point;
    Utilities::Range<double> theta_range(HUGE_VALF,-1);
    
    /// header line
    getline(file,myline);
    if(verbose) std::cout << myline << std::endl;
    if(myline[0] != 'x'){
      std::cerr << "File " << filename << " is in an unexpected format." << std::endl;
      throw std::runtime_error(" Cannot open file.");
    }
    
    
    int pos = myline.find("mass");
    if(pos == std::string::npos){
      particle_masses.push_back(particle_mass);
      multimass = false;
    }else{
      multimass = true;
    }
    
    int pos1 = pos;
    pos = myline.find("size");
    if(pos1 > pos && pos1 != std::string::npos){
      std::cout << myline << std::endl;
      std::cerr << "File " << filename << " is in an unexpected format. Particle mass must come before particles size."<< std::endl;
      throw std::runtime_error(" Cannot open file.");
    }
    
    while (getline(file,myline) ) {
      
      if(verbose) std::cout << myline << std::endl;
      
      int pos1 = 0;
      int pos = myline.find(delim);
      
      point[0] = std::stod(myline.substr(pos1,pos));
      
      //std::cout << xp[i][0] << std::endl;
      
      pos1 = pos + 1;
      pos = myline.find(delim,pos1);
      
      point[1] = std::stod(myline.substr(pos1,pos));
      
      pos1 = pos + 1;
      pos = myline.find(delim,pos1);
      //std::cout << xp[i][1] << std::endl;
      point[2] =std::stod(myline.substr(pos1,pos));
      //std::cout << point << std::endl;
      
      particles.push_back(point);
      theta_range.update(center.AngleSeporation(particles.back()));
      
      //std::cout << myline.length() << std::endl;
      
      if(multimass){
        if(pos != std::string::npos){
          pos1 = pos + 1;
          pos = myline.find(delim,pos1);
          //std::cout << xp[i][1] << std::endl;
          particle_masses.push_back(std::stof(myline.substr(pos1,pos)));
        }else{
          particle_masses.push_back(0.0);
        }
        std::cout << particle_masses.back() << std::endl;
      }
      
      ++i;
    }
    
    size_t Npoints = particles.size();
    
    file.close();
    
    // sort by distance
    if(!multimass){
      std::sort(particles.begin(),particles.end()
                ,[](SphericalPoint &p1,SphericalPoint &p2){return p1.r < p2.r;});
      
    }else{
      std::vector<size_t> index(particles.size());
      
      std::sort(index.begin(),index.end()
                ,[&particles](size_t i,size_t j){return particles[i].r < particles[j].r;});
      
      {
        Utilities::apply_permutation_in_place(particles,index);
      }
      
      assert(particles[0].r < particles.back().r);
      
      if(multimass){
        assert(Npoints == particle_masses.size());
        Utilities::apply_permutation_in_place(particle_masses,index);
      }
      
    }
    
    assert(particles[0].r <= particles.back().r);
    std::cout << particles.size() << " particles read in from cone."
    << std::endl;
    
    std::vector<double> D_planes(Nplanes);
    
    //double vol = pi*theta_range.max()*theta_range.max()*(pow(particles.back().r,3) - pow(particles[0].r,3))/3;
    //double psize = pow(vol/particles.size(),1.0/3.0);
    
    // break up into slices redshift bins
    for(int i = 0 ; i < Nplanes ; ++i ) D_planes[i] = (2*i+1)*(particles.back().r - particles[0].r)/Nplanes/2 + particles[0].r;
    
    //Utilities::RandomNumbers_NR ran(1234);
    //PosType **xp;
    std::deque<SphericalPoint>::iterator it2;
    size_t i1 = 0;
    
    Point_2d center_cart;
    
    PixelMap massmap(center_cart.x,(size_t)(sqrt(2)*theta_range.max()/angular_resolution)
                     ,(size_t)(sqrt(2)*theta_range.max()/angular_resolution),angular_resolution);
    //PixelMap massmap(center_cart.x,(size_t)(2*theta_range.max()/angular_resolution)
    //                 ,(size_t)(2*theta_range.max()/angular_resolution),angular_resolution);
    
    size_t count = 0;
    for(int i = 0 ; i < Nplanes ; ++i ){
      
      it2 = std::lower_bound(particles.begin(),particles.end()
                             ,(i+1)*particles.back().r*1.001/Nplanes
                             ,[](SphericalPoint &p,double v){return p.r < v;});
      
      double z = cosmo.invCoorDist(D_planes[i]);
      //double Dl = D_planes[i]/(1+z);
      double total_mass = 0.0;
      
      // project onto planes
      size_t Np_on_plane = it2 - particles.begin();
      
      if(verbose) std::cout << "Number of particles on plane " << i << " : " << Np_on_plane
        << " z = " << z << std::endl;
      
      //xp = Utilities::PosTypeMatrix(Np_on_plane,2);
      size_t j = i1;
      long index = -1;
      for(int ii = 0 ; ii < Np_on_plane ; ++ii,++j){
        //for(auto it = it1 ; it != it2 + 1 ; ++it){
        //xp[i][0] = (*it).theta*Dl;
        //xp[i][1] = (*it).phi*Dl;
        //xp[i][0] = particles.front().theta*Dl;
        //xp[i][1] = particles.front().phi*Dl;
        
        index = massmap.find_index(particles.front().theta,particles.front().phi);
        
        if(index != -1){
          massmap[index] += particle_masses[j*multimass];
          total_mass += particle_masses[j*multimass];
          
          //std::cout << " index " << index << "  " << particles.front().theta
          //  << " " << particles.front().phi << std::endl;
          ++count;
        }
        /// test lines for randomizing positions
        //double tmp = theta_range.max()*sqrt(ran()),tmp_t = 2*pi*ran();
        //xp[i][0] = tmp*cos(tmp_t)*Dl;
        //xp[i][1] = tmp*sin(tmp_t)*Dl;
        
        particles.pop_front();
      }
      
      //double sigma_back = total_mass/pi/Dl/Dl/theta_range.max()/theta_range.max();
      
      //std::cout << " total mass " << massmap.ave()*massmap.size() << std::endl;
      //massmap.printFITS("!test_map.fits");
      lensVec.push_back(new LensHaloMassMap(massmap,1.0,z,1.0,true,cosmo));
      
      /*{
       size_t tmp = lensVec.back()->map->nx*lensVec.back()->map->ny;
       double temp =0;
       for(size_t i=0;i<tmp ; ++i){
       std::cout <<
       lensVec.back()->map->convergence[i] << "," <<
       lensVec.back()->map->gamma1[i] << "," <<
       std::endl;
       temp += lensVec.back()->map->convergence[i];
       }
       std::cout << " total " << temp << std::endl;
       }*/
      
      i1 += Np_on_plane;
      
      massmap.Clean();
    }
    
    //std::cout << "fraction of particles on map : "  << count*1.0/Npoints << std::endl;
    
    assert(i1 == Npoints);
  }
  
  void WriteLightCone(std::string filename,std::vector<DataRockStar> &vec){
    
    std::ofstream file(filename);
    file << "ID,x,y,z,mass,Rvir(kpc),Rscale(kpc),Vmax(km/s)" << std::endl;
    
    for(auto h : vec){
      file << h.id << "," << h.x[0] << "," << h.x[1] << "," << h.x[2] << "," << h.mass
      << "," << h.Rvir << "," << h.Rscale << "," << h.Vmax << std::endl;
    }
    
    file.close();
  }
  
  void WriteLightCone(std::string filename,std::vector<Point_3d> &vec){
    
    std::ofstream file(filename);
    file << "x,y,z" << std::endl;
    
    for(auto h : vec){
      file << h[0] << "," << h[1] << "," << h[2] << std::endl;
    }
    
    file.close();
  }
  
  void MultiLightCone::ReadBoxRockStar(std::string filename
                                       ,double rlow,double rhigh
                                       ,std::vector<std::vector<DataRockStar> > &conehalos
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
    std::vector<DataRockStar> boxhalos;
    boxhalos.reserve(blocksize);
    
    //const int ncolumns = 11;
    const int ncolumns = 42;
    
    void *addr[ncolumns];
    DataRockStar halo;
    double tmp;
    double parent_id;
    double m200c,m200b,Mvir_all,rs_klypin;
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
          Utilities::LockableContainer<std::vector<DataRockStar> > halos;
          halos.swap(conehalos[j]);
          
          assert(nthreads*chunk_size + remainder == boxhalos.size() );
          for(int ii =0; ii< nthreads ; ++ii){
            
            //std::cout << ii*chunk_size << " " << n << std::endl;
            
            thr[ii] = std::thread(&LightCone::select<DataRockStar>,cones[j]
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
  
  void ASCII_XV::fastplanes_parallel(
                                     Point_3d *begin
                                     ,Point_3d *end
                                     ,const COSMOLOGY &cosmo
                                     ,std::vector<std::vector<Point_3d> > &boxes
                                     ,std::vector<Point_3d> &observers
                                     ,std::vector<Quaternion> &rotationQs
                                     ,std::vector<double> &dsources
                                     ,std::vector<std::vector<PixelMap> > &maps
                                     ,double dmin
                                     ,double dmax
                                     ,double BoxLength
                                     ){
    // loop lines / read
    int Ncones = maps.size();
    int Nmaps = maps[0].size();
    
    
    for(Point_3d *phalo = begin ; phalo != end ; ++phalo){
      
      //*phalo /= cosmo.gethubble();
      
      // loop cones
      for(int icone=0;icone<Ncones;++icone){
        
        // loop through repitions of box ??? this could be done better
        for(auto dn : boxes[icone]){
          
          //for(dn[0] = min_box[icone][0] ; dn[0] <= max_box[icone][0] ; ++dn[0]){
          //for(dn[1] = min_box[icone][1] ; dn[1] <= max_box[icone][1] ; ++dn[1]){
          //for(dn[2] = min_box[icone][2] ; dn[2] <= max_box[icone][2] ; ++dn[2]){
          
          Point_3d x = *phalo - observers[icone] + dn*BoxLength;
          double r = x.length();
          if( r > dmin && r < dmax ){
            // rotate to cone frame - direction[i] is the x-axis
            x = rotationQs[icone].Rotate(x);
            SphericalPoint sp(x);
            
            // find pixel
            long image_index = maps[icone][0].find_index(sp.theta,sp.phi);
            
            if(image_index != -1){
              for(int isource = 0 ; isource < Nmaps ; ++isource){
                if(dsources[isource] > sp.r  ){
                  //std::lock_guard<std::mutex> lock(moo);
                  // add mass or distribute mass to pixels
                  maps[icone][isource][image_index] += (dsources[isource] - sp.r)/sp.r;  // this is assuming flat ???
                }
              }
            }
          }
        }//}}
      }
      
    }
    
  }
  
  
  void ASCII_XM::fastplanes_parallel(
                                     LightCones::DatumXM *begin
                                     ,LightCones::DatumXM *end
                                     ,const COSMOLOGY &cosmo
                                     ,std::vector<std::vector<Point_3d> > &boxes
                                     ,std::vector<Point_3d> &observers
                                     ,std::vector<Quaternion> &rotationQs
                                     ,std::vector<double> &dsources
                                     ,std::vector<std::vector<PixelMap> > &maps
                                     ,double dmin
                                     ,double dmax
                                     ,double BoxLength
                                     ){
    // loop lines / read
    
    int Ncones = maps.size();
    int Nmaps = maps[0].size();
    //double hubble = cosmo.gethubble();
    
    for(auto *phalo = begin ; phalo != end ; ++phalo){
      
      //phalo->x /= hubble;
      
      // loop cones
      for(int icone=0;icone<Ncones;++icone){
        
        // loop through repitions of box ??? this could be done better
        //Point_3d dn;
        for(auto dn : boxes[icone]){
          
          //for(dn[0] = min_box[icone][0] ; dn[0] <= max_box[icone][0] ; ++dn[0]){
          //for(dn[1] = min_box[icone][1] ; dn[1] <= max_box[icone][1] ; ++dn[1]){
          //for(dn[2] = min_box[icone][2] ; dn[2] <= max_box[icone][2] ; ++dn[2]){
          
          Point_3d x = phalo->x - observers[icone] + dn*BoxLength;
          double r = x.length();
          if( r > dmin && r < dmax ){
            // rotate to cone frame - direction[i] is the x-axis
            x = rotationQs[icone].Rotate(x);
            SphericalPoint sp(x);
            
            // find pixel
            long image_index = maps[icone][0].find_index(sp.theta,sp.phi);
            
            if(image_index != -1){
              for(int isource = 0 ; isource < Nmaps ; ++isource){
                if(dsources[isource] > sp.r  ){
                  //std::lock_guard<std::mutex> lock(moo);
                  // add mass or distribute mass to pixels
                  maps[icone][isource][image_index]
                  += phalo->mass*(dsources[isource] - sp.r)/sp.r;  // this is assuming flat ???
                }
              }
            }
          }
        }//}}
      }
    }
  }
  
  
  
  void ASCII_XMR::fastplanes_parallel(
                                      LightCones::DatumXMR *begin
                                      ,LightCones::DatumXMR *end
                                      ,const COSMOLOGY &cosmo
                                      ,std::vector<std::vector<Point_3d> > &boxes
                                      ,std::vector<Point_3d> &observers
                                      ,std::vector<Quaternion> &rotationQs
                                      ,std::vector<double> &dsources
                                      ,std::vector<std::vector<PixelMap> > &maps
                                      ,double dmin
                                      ,double dmax
                                      ,double BoxLength
  //,std::mutex &moo
  ){
    // loop lines / read
    
    int Ncones = maps.size();
    int Nmaps = maps[0].size();
    const double resolution = maps[0][0].getResolution();
    Point_2d p1 = maps[0][0].getLowerLeft();
    const size_t Nx = maps[0][0].getNx();
    const size_t Ny = maps[0][0].getNy();
    const size_t Nxm1 = Nx-1;
    const size_t Nym1 = Ny-1;
    const double half_range =maps[0][0].getRangeX()/2;
    const double hubble = cosmo.gethubble();
    const double pix_area = resolution*resolution;

    
    for(auto *phalo = begin ; phalo != end ; ++phalo){
      
      // loop cones
      for(int icone=0;icone<Ncones;++icone){
        
        for(auto dn : boxes[icone]){
          // loop through repitions of box ??? this could be done better
          
          Point_3d x = phalo->x - observers[icone] + dn*BoxLength;
          double r = x.length();
          if( r > dmin && r < dmax){
            // rotate to cone frame - direction[i] is the x-axis
            x = rotationQs[icone].Rotate(x);
            if(x[0] < 0 ) continue;
            
            SphericalPoint sp(x);
            
            double ang_radius = 2*phalo->r/sp.r;
            if(fabs(sp.theta) - ang_radius > half_range || fabs(sp.phi) - ang_radius > half_range) continue;
            ang_radius /= resolution;
            
            double dx = (sp.theta - p1[0])/resolution ;
            double dy = (sp.phi   - p1[1])/resolution ;
            //*****************************************
            long jjmin = (long)MAX(dy-ang_radius,0);
            long jjmax = (long)MIN(dy+ang_radius,Nym1);
            if(jjmin > jjmax) continue;
            long iimin = (long)MAX(dx-ang_radius,0);
            long iimax = (long)MIN(dx+ang_radius,Nxm1);
            if(iimin > iimax) continue;
            
            ang_radius /= 2;
            double ascale = 1.0/(cosmo.invCoorDist(sp.r/hubble)+1);
            
            
            if(ang_radius < 2){
              double m_part = phalo->mass/ascale/sp.r;

              size_t index = (long)( (iimin + iimax)/2 + 0.5 ) + Nx*(size_t)( (jjmin + jjmax)/2 + 0.5 );
              
              for(int isource = 0 ; isource < Nmaps ; ++isource){
                if(dsources[isource] > sp.r  ){
                  // add mass or distribute mass to pixels
                  maps[icone][isource][index]
                  += m_part*(dsources[isource] - sp.r);  // this is assuming flat ???
                }
              }
            }else{

              double m_halo = phalo->mass*pix_area/ascale
              /phalo->r/phalo->r;

              //double area = ang_radius*ang_radius;
              
              for(long jj = jjmin ; jj <= jjmax ; ++jj){
                size_t index = iimin + Nx*jj;
                for(long ii = iimin ; ii <= iimax ; ++ii,++index){
                  
                  double q = sqrt( (ii - dx)*(ii - dx) + (jj - dy)*(jj - dy) )/ang_radius;
                  //double m = Profiles::Bspline<2>(q) * m_Dl/area;
                  double m = Profiles::Bspline<2>(q)*m_halo;

                    for(int isource = 0 ; isource < Nmaps ; ++isource){
                      if(dsources[isource] > sp.r  ){
                        // add mass or distribute mass to pixels
                        maps[icone][isource][index]
                        += m*(dsources[isource] - sp.r);  // this is assuming flat ???
                      }
                    }
                  
                }
              }
            }
            
            /*****************************************
             
             if(dx + size > 0 && dx - size < Nx &&
             dy + size > 0 && dy - size < Ny){
             long image_index = (long)(dx) + Nx*(long)(dy);
             
             if(dx < 0 || dy < 0 ){
             dx = 0;
             }else{
             long ix = (long)(dx);
             long iy = (long)(dy);
             
             //long image_index = ix + Nx*iy;
             
             dx -= ix;
             dy -= iy;
             }
             
             if( dx < size || (1 - dx) < size ||
             dy < size || (1 - dy) < size ){  /// halo needs to be resolved
             
             // subdivide halo
             
             Point_3d p;
             for(size_t ii = 0 ; ii < Nsub ; ++ii){;
             pgen.draw(p);
             SphericalPoint sp(p*phalo->r + x);
             subindex[ii] = maps[icone][0].find_index(sp.theta,sp.phi);
             }
             
             for(int isource = 0 ; isource < Nmaps ; ++isource){
             if(dsources[isource] > sp.r  ){
             
             std::lock_guard<std::mutex> lock(moo);
             // add mass or distribute mass to pixels
             for(long ii : subindex){
             if(ii != -1) maps[icone][isource][ii]
             += phalo->mass*(dsources[isource] - sp.r)/sp.r/Nsub;  // this is assuming flat ???
             }
             }
             }
             
             }else{
             
             for(int isource = 0 ; isource < Nmaps ; ++isource){
             if(dsources[isource] > sp.r  ){
             std::lock_guard<std::mutex> lock(moo);
             // add mass or distribute mass to pixels
             maps[icone][isource][image_index]
             += phalo->mass*(dsources[isource] - sp.r)/sp.r;  // this is assuming flat ???
             }
             }
             }
             
             }*/
          }
        }
      }
    }
  }
  
  
  void ASCII_XMRRT::fastplanes_parallel(
                                        LightCones::DatumXMRmRs *begin
                                        ,LightCones::DatumXMRmRs *end
                                        ,const COSMOLOGY &cosmo
                                        ,std::vector<std::vector<Point_3d> > &boxes
                                        ,std::vector<Point_3d> &observers
                                        ,std::vector<Quaternion> &rotationQs
                                        ,std::vector<double> &dsources
                                        ,std::vector<std::vector<PixelMap> > &maps
                                        ,double dmin
                                        ,double dmax
                                        ,double BoxLength
                                        ){
    // loop lines / read
    
    int Ncones = maps.size();
    int Nmaps = maps[0].size();
    const double resolution = maps[0][0].getResolution();
    Point_2d p1 = maps[0][0].getLowerLeft();
    const size_t Nx = maps[0][0].getNx();
    const size_t Ny = maps[0][0].getNy();
    const size_t Nxm1 = Nx-1;
    const size_t Nym1 = Ny-1;
    const double half_range =maps[0][0].getRangeX()/2;
    const double hubble = cosmo.gethubble();
    const double pix_area = resolution*resolution;

    Profiles::NFW2D profile;
    
    for(auto *phalo = begin ; phalo != end ; ++phalo){
      
      // loop cones
      for(int icone=0;icone<Ncones;++icone){
        
        for(auto dn : boxes[icone]){
          // loop through repetitions of box ??? this could be done better
          
          Point_3d x = phalo->x - observers[icone] + dn*BoxLength;
          double r = x.length();
          if( r > dmin && r < dmax){
            // rotate to cone frame - direction[i] is the x-axis
            x = rotationQs[icone].Rotate(x);
            if(x[0] < 0 ) continue;
            
            SphericalPoint sp(x);
            
            double size = phalo->r_max/sp.r;
            
            if(fabs(sp.theta) - size > half_range || fabs(sp.phi) - size > half_range) continue;
            
            size /= resolution;
            
            double dx = (sp.theta - p1[0])/resolution ;
            double dy = (sp.phi   - p1[1])/resolution ;
            
            //*****************************************
            
            
            long jjmin = (long)MAX(dy-size,0);
            long jjmax = (long)MIN(dy+size,Nym1);
            if(jjmin > jjmax) continue;
            long iimin = (long)MAX(dx-size,0);
            long iimax = (long)MIN(dx+size,Nxm1);
            if(iimin > iimax) continue;
            
            // this is m / Dl with no hubble constant
            double ascale = 1.0/(cosmo.invCoorDist(sp.r/hubble)+1);
            

            if(size < 2){
              double m_part = phalo->mass/ascale/sp.r;
             
              size_t index = (long)( (iimin + iimax)/2 + 0.5 ) + Nx*(size_t)( (jjmin + jjmax)/2 + 0.5 );
              
              for(int isource = 0 ; isource < Nmaps ; ++isource){
                if(dsources[isource] > sp.r  ){
                  // add mass or distribute mass to pixels
                  maps[icone][isource][index]
                  += m_part*(dsources[isource] - sp.r);  // this is assuming flat ???
                }
              }
            }else{
              double m_halo = phalo->mass*pix_area/ascale
              /phalo->r_scale/phalo->r_scale;

              float con = phalo->r_max/phalo->r_scale;
              
              // change to pixel scale size
              size /= con;
              
              //Profiles::TNFW2D profile(con);
              
              for(long jj = jjmin ; jj <= jjmax ; ++jj){
                size_t index = iimin + Nx*jj;
                for(long ii = iimin ; ii <= iimax ; ++ii,++index){
                  
                  double q = sqrt( (ii - dx)*(ii - dx) + (jj - dy)*(jj - dy) )/size;
                  double m = profile(q,con)*m_halo;
                  //std::cout << m << " q " << q << " " << std::endl;
                  
                    for(int isource = 0 ; isource < Nmaps ; ++isource){
                      if(dsources[isource] > sp.r  ){
                        // add mass or distribute mass to pixels
                        maps[icone][isource][index]
                        += m*(dsources[isource] - sp.r);  // this is assuming flat ???
                      }
                    }
                }
              }
              
              
              
            }
          }
        }
      }
    }
  }
  
  void random_observers(std::vector<Point_3d> &observers
                        ,std::vector<Point_3d> &directions
                        ,int Ncones
                        ,double BoxLength
                        ,double cone_opening_radius
                        ,Utilities::RandomNumbers_NR &ran
                        
                        ){
    if(cone_opening_radius > pi/4){
      throw std::invalid_argument("cone too wide");
    }
    observers.resize(Ncones);
    directions.resize(Ncones);
    double costheta = cos(cone_opening_radius);
    
    
    for(auto &x : observers){
      x[0] = ran()*BoxLength;
      x[1] = ran()*BoxLength;
      x[2] = ran()*BoxLength;
    }
    for(auto &v : directions){
      double theta;
      do{
        theta = 2*pi*ran();
        v[0] = 2*ran() - 1;
        v[1] = sqrt(1-v[0]*v[0])*cos(theta);
        v[2] = sqrt(1-v[0]*v[0])*sin(theta);
        // prevent the coordinate axes from being inside the cone
      }while( fabs(v[0]) > costheta || fabs(v[1]) > costheta || fabs(v[2]) > costheta );
    }
  }
  
}

