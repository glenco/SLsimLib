//
//  lightcone_construction.h
//  GLAMER
//
//  Created by Ben Metcalf on 26/10/16.
//
//

#ifndef __GLAMER__lightcone_construction__
#define __GLAMER__lightcone_construction__

#include <stdio.h>
#include "point.h"
#include "lens_halos.h"
#include "geometry.h"


namespace LightCones{
  
  struct DataRockStar{
    unsigned int id;
    Point_3d x;
    double mass;
    double Rvir;     // in comoving units
    double Rscale;   // in comoving units
    double Vmax;
  };
  
  struct DatumXM{
    Point_3d x;
    double mass;
  };
  struct DatumXMRmRs{
    Point_3d x;
    double mass;
    double r_max;
    double r_scale;
  };
  struct DatumXMR{
    Point_3d x;
    double mass;
    double r;
  };
  

  
  void ReadLightConeNFW(std::string filename,COSMOLOGY &cosmo
                        ,std::vector<LensHalo* > &lensVec
                        ,PosType &theta_max);
  
  void ReadLightConeParticles(std::string filename,COSMOLOGY &cosmo
                              ,std::vector<LensHaloParticles *> &lensVec
                              ,int Nplanes
                              ,float particle_mass
                              ,float particle_size
                              ,bool angular_sizes = false
                              ,bool verbose = false);
  
  void ReadLightConeParticles(std::string filename,COSMOLOGY &cosmo
                              ,std::vector<LensHaloMassMap *> &lensVec
                              ,int Nplanes
                              ,float particle_mass
                              ,float angular_resolution
                              ,bool verbose = false);
  
  /*** \brief Class for constructing light cones form simulation boxes.
   *
   *
   */
  class LightCone{
  public:
    LightCone(double angular_radius);
    
    void WriteLightCone(std::string filename,std::vector<DataRockStar> &vec);
    void WriteLightCone(std::string filename,std::vector<Point_3d> &vec);
    
    friend class MultiLightCone;
    
    template <typename T>
    void select(Point_3d xo,Point_3d v,double Length,double rlow,double rhigh
                ,T* begin
                ,T* end
                ,Utilities::LockableContainer<std::vector<T> > &incone
                ,bool periodic_boundaries = true);
    
  private:
    
    double r_theta; // angular radius of light cone (radians)
    double sin_theta_sqrt;
  };
  /// select the halos from the box that are within the light cone
  template <typename T>
  void LightCone::select(Point_3d xo,Point_3d v,double Length,double rlow,double rhigh
                         ,T* begin
                         ,T* end
                         ,Utilities::LockableContainer<std::vector<T> > &incone
                         ,bool periodic_boundaries){
    
    if(begin == end) return;
    
    Point_3d dx,x;
    double r2,xp;
    Point_3d n;
    double rhigh2 = rhigh*rhigh;
    double rlow2 = rlow*rlow;
    T b;
    
    v /= v.length();
    
    // standard coordinates so cone is always centered on the x-axis
    Point_3d y_axis,z_axis;
    
    y_axis[0] = -v[1];
    y_axis[1] =  v[0];
    y_axis[2] =  0;
    
    y_axis /= y_axis.length();
    z_axis = v.cross(y_axis);
    
    //std::cout << v*y_axis << " " << v*z_axis << " " << y_axis*z_axis << std::endl;
    
    // we can do better than this !!!!
    int n0[2] = {(int)((rhigh + xo[0])/Length),(int)((xo[0] - rhigh)/Length) - 1 };
    int n1[2] = {(int)((rhigh + xo[1])/Length),(int)((xo[1] - rhigh)/Length) - 1 };
    int n2[2] = {(int)((rhigh + xo[2])/Length),(int)((xo[2] - rhigh)/Length) - 1 };
    if(!periodic_boundaries){
      n0[0] = n0[1] = 0;
      n1[0] = n1[1] = 0;
      n2[0] = n2[1] = 0;
    }
    for(auto a = begin ; a != end ; ++a){
      dx[0] = (*a).x[0] - xo[0];
      dx[1] = (*a).x[1] - xo[1];
      dx[2] = (*a).x[2] - xo[2];
      
      int count = 0;
      
      for(n[0] = n0[1]; n[0] <= n0[0] ; ++n[0]){
        for(n[1] = n1[1]; n[1] <= n1[0] ; ++n[1]){
          for(n[2] = n2[1]; n[2] <= n2[0] ; ++n[2]){
            
            x = dx + n*Length;
            xp = x*v;
            if(xp > 0){
              r2 = x.length_sqr();
              
              if(r2 > rlow2 && r2 < rhigh2){
                if( r2*sin_theta_sqrt > r2 - xp*xp){
                  //std::cout << n << " | " << xo << " | "
                  //          << v << " |  " << x << std::endl;
                  
                  b = *a;
                  //b.x = x;
                  // rotate to standard reference frame
                  b.x[0] = xp;
                  b.x[1] = y_axis*x;
                  b.x[2] = z_axis*x;
                  
                  //std::cout << "b = " << b.x << std::endl;
                  
                  incone.push_back(b);
                  ++count;
                }
              }
            }
          }
        }
      }
      if(count > 1) std::cout << " Warning: Same halo appears " << count << " times in light-cone." << std::endl;
    }
  }
  
  class MultiLightCone{
  public:
    MultiLightCone(
                   double angular_radius
                   ,const std::vector<Point_3d> &observers    /// postion of observers within the simulation box
                   ,const std::vector<Point_3d> &directions   /// direction of light cones
    ):
    xos(observers),vs(directions)
    {
      assert(observers.size() == directions.size());
      for(int i=0;i<observers.size(); ++i) cones.push_back(angular_radius);
    }
    
    /** Read the points in from a snapshot of the halos created by RockStar.
     The output will be the halos in the cone in coordinates where the x-axis is
     the center of the cone.  Changing the observer, xo, and direction of view, V,
     effectively translates and rotates the box, respectively.  Only halos with radii
     between rlow and rhigh are added.
     */
    void ReadBoxRockStar(std::string filename
                         ,double rlow,double rhigh
                         ,std::vector<std::vector<LightCones::DataRockStar> > &conehalos
                         ,bool periodic_boundaries = true
                         ,bool allow_subhalos = false);
    
    void ReadBoxXYZ(std::string filename
                    ,double rlow,double rhigh
                    ,std::vector<std::vector<Point_3d> > &conehalos
                    ,double hubble
                    ,double BoxLength
                    ,bool periodic_boundaries = true
                    );
    
    void ReadBoxToMaps(std::string filename
                       ,double rlow,double rhigh
                       ,std::vector<std::vector<PixelMap>  > &maps
                       ,double hubble
                       ,double BoxLength
                       ,bool periodic_boundaries = true
                       );
    
    
  private:
    
    std::vector<Point_3d> xos;
    std::vector<Point_3d> vs;
    std::vector<LightCone> cones;
  };
  
  /** \brief
   
   The Born approximation is used.
   
   *
  template<typename T>
  void FastLightCones(
                      const COSMOLOGY &cosmo
                      ,const std::vector<double> &zsources    // vector of source redshifts desired
                      ,std::vector<std::vector<PixelMap> > &maps
                      ,double range
                      ,double angular_resolution
                      ,std::vector<Point_3d> &observers    /// position of observers within the simulation box
                      ,std::vector<Point_3d> &directions   /// direction of light cones
                      ,const std::vector<std::string> &snap_filenames
                      ,const std::vector<float> &snap_redshifts
                      ,double BoxLength
                      ,double particle_mass
                      ,bool verbose = false
                      ,bool addtocone = false  /// if false the maps will be cleared and new maps made, if false particles are added to the existing maps
  );/**/

  using Utilities::Geometry::Quaternion;
  using Utilities::Geometry::SphericalPoint;
  
  /*************************************************************************************************
   templated functions for projecting mass onto planes
   *************************************************************************************************/
  
  template<typename T>
  void _fastplanes_parallel_(T *begin,T *end
                             ,const COSMOLOGY &cosmo
                             ,std::vector<Point_3d> &max_box
                             ,std::vector<Point_3d> &min_box
                             ,std::vector<Point_3d> &observers
                             ,std::vector<Quaternion> &rotationQs
                             ,std::vector<double> &dsources
                             ,std::vector<std::vector<PixelMap> > &maps
                             ,double dmin
                             ,double dmax
                             ,double BoxLength
                             ,std::mutex &moo
                             ){
    // loop lines / read
    
    int Ncones = maps.size();
    int Nmaps = maps[0].size();
    
    for(Point_3d *phalo = begin ; phalo != end ; ++phalo){
      
      *phalo /= cosmo.gethubble();
      
      // loop cones
      for(int icone=0;icone<Ncones;++icone){
        
        // loop through repitions of box ??? this could be done better
        Point_3d dn;
        for(dn[0] = min_box[icone][0] ; dn[0] <= max_box[icone][0] ; ++dn[0]){
          for(dn[1] = min_box[icone][1] ; dn[1] <= max_box[icone][1] ; ++dn[1]){
            for(dn[2] = min_box[icone][2] ; dn[2] <= max_box[icone][2] ; ++dn[2]){
              
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
                      std::lock_guard<std::mutex> lock(moo);
                      // add mass or distribute mass to pixels
                      maps[icone][isource][image_index] += (dsources[isource] - sp.r)/sp.r;  // this is assuming flat ???
                    }
                  }
                }
              }
            }}}
      }
      
    }
    
  }
  //template<typename T>

  

  /*************************************************************************************************
   templated functions for reading a block from a file
   *************************************************************************************************/
   

  template <typename T>
  size_t scan_block(size_t blocksize,std::vector<T> &points,FILE *pFile){
    double tmpf;
    // read in a block of points
    size_t i=0;
    
    points.resize(blocksize);
    while(i < blocksize &&
          fscanf(pFile,"%lf %lf %lf %lf %lf %lf"
                 ,&points[i][0],&points[i][1],&points[i][2],&tmpf,&tmpf,&tmpf) != EOF)
      ++i;
    points.resize(i);
    
    return i;
  }
  
  struct ASCII_XV{
    std::vector<Point_3d> points;
    
    size_t scan_block(size_t blocksize,FILE *pFile){
      double tmpf;
      // read in a block of points
      size_t i=0;
      
      points.resize(blocksize);
      while(i < blocksize &&
            fscanf(pFile,"%lf %lf %lf %lf %lf %lf"
                   ,&points[i][0],&points[i][1],&points[i][2],&tmpf,&tmpf,&tmpf) != EOF)
        ++i;
      points.resize(i);
      
      return i;
    }
    
    static void fastplanes_parallel(
                         Point_3d *begin
                         ,Point_3d *end
                         ,const COSMOLOGY &cosmo
                         ,std::vector<Point_3d> &max_box
                         ,std::vector<Point_3d> &min_box
                         ,std::vector<Point_3d> &observers
                         ,std::vector<Quaternion> &rotationQs
                         ,std::vector<double> &dsources
                         ,std::vector<std::vector<PixelMap> > &maps
                         ,double dmin
                         ,double dmax
                         ,double BoxLength
                         ,std::mutex &moo
                         );
    
  };
  struct ASCII_XM{
    std::vector<DatumXM> points;
    size_t scan_block(size_t blocksize,FILE *pFile){
      
      // read in a block of points
      size_t i=0;
      
      points.resize(blocksize);
      while(i < blocksize &&
            fscanf(pFile,"%lf %lf %lf %lf"
                   ,&points[i].x[0],&points[i].x[1],&points[i].x[2],&points[i].mass) != EOF)
        ++i;
      points.resize(i);
      
      return i;
    }
    static void fastplanes_parallel(
                                    DatumXM *begin
                                    ,DatumXM *end
                                    ,const COSMOLOGY &cosmo
                                    ,std::vector<Point_3d> &max_box
                                    ,std::vector<Point_3d> &min_box
                                    ,std::vector<Point_3d> &observers
                                    ,std::vector<Quaternion> &rotationQs
                                    ,std::vector<double> &dsources
                                    ,std::vector<std::vector<PixelMap> > &maps
                                    ,double dmin
                                    ,double dmax
                                    ,double BoxLength
                                    ,std::mutex &moo
                                    );
    

  };
  struct ASCII_XMR{
    std::vector<DatumXMR> points;

    size_t scan_block(size_t blocksize,FILE *pFile){
      
      // read in a block of points
      size_t i=0;
      
      points.resize(blocksize);
      while(i < blocksize &&
            fscanf(pFile,"%lf %lf %lf %lf %lf"
                   ,&points[i].x[0],&points[i].x[1],&points[i].x[2]
                   ,&points[i].mass,&points[i].r) != EOF)
        ++i;
      points.resize(i);
      
      return i;
    }
    static void fastplanes_parallel(
                                    DatumXMR *begin
                                    ,DatumXMR *end
                                    ,const COSMOLOGY &cosmo
                                    ,std::vector<Point_3d> &max_box
                                    ,std::vector<Point_3d> &min_box
                                    ,std::vector<Point_3d> &observers
                                    ,std::vector<Quaternion> &rotationQs
                                    ,std::vector<double> &dsources
                                    ,std::vector<std::vector<PixelMap> > &maps
                                    ,double dmin
                                    ,double dmax
                                    ,double BoxLength
                                    ,std::mutex &moo
                                    );
    


  };
  struct ASCII_XMRmRs{
    std::vector<DatumXMRmRs> points;
    size_t scan_block(size_t blocksize,FILE *pFile){
      
      // read in a block of points
      size_t i=0;
      
      points.resize(blocksize);
      while(i < blocksize &&
            fscanf(pFile,"%lf %lf %lf %lf %lf %lf"
                   ,&points[i].x[0],&points[i].x[1],&points[i].x[2]
                   ,&points[i].mass,&points[i].r_max,&points[i].r_scale) != EOF)
        ++i;
      points.resize(i);
      
      return i;
    }
    
    static void fastplanes_parallel(
                                    DatumXMRmRs *begin
                                    ,DatumXMRmRs *end
                                    ,const COSMOLOGY &cosmo
                                    ,std::vector<Point_3d> &max_box
                                    ,std::vector<Point_3d> &min_box
                                    ,std::vector<Point_3d> &observers
                                    ,std::vector<Quaternion> &rotationQs
                                    ,std::vector<double> &dsources
                                    ,std::vector<std::vector<PixelMap> > &maps
                                    ,double dmin
                                    ,double dmax
                                    ,double BoxLength
                                    ,std::mutex &moo
                                    );
    

  };

  
  /*************************************************************************************************
   *************************************************************************************************/

  /** \brief class for generating positions in proportion to mass in an NFW profiles
   */
  class NFWgenerator{
  public:
    NFWgenerator(Utilities::RandomNumbers_NR &ran_in,double max_cons);
    /// returns a vector of points drawn from a spherical halo
    void drawSpherical(std::vector<Point_3d> &points  /// output points
                       ,double cons                   /// concentration
                       ,double Rvir                   /// maximum elliptical radius
    );
    ///  STILL UNDER CONSTRUCTION returns a vector of points drawn from a triaxial halo,
    void drawTriAxial(std::vector<Point_3d> &points  /// output points
                      ,double cons                   /// concentration
                      ,double Rvir                   /// maximum elliptical radius
                      ,double f1                     /// axis ratio 1 to 3
                      ,double f2                     /// axis ratio 2 to 3
                      ,SphericalPoint v              /// direction of axis 3
    );
  private:
    Utilities::RandomNumbers_NR &ran;
    double dx;
    std::vector<double> X;
    std::vector<double> F;
    const int N = 1000;
  };
  
  struct DkappaDz{
    DkappaDz(const COSMOLOGY &cos,double zsource):cosmo(cos),zs(zsource){
      rho = cosmo.getOmega_matter()*cosmo.rho_crit(0);
    };
    
    double operator()(double z){
      double x = 1+z;
      return cosmo.drdz(x)*rho*x*x*x/cosmo.SigmaCrit(z,zs)/cosmo.getHubble();
    }
    
    const COSMOLOGY &cosmo;
    double zs;
    double rho;
  };

  using Utilities::Geometry::Quaternion;
  using Utilities::Geometry::SphericalPoint;
  
  template <typename T>
  void FastLightCones(
                      const COSMOLOGY &cosmo
                      ,const std::vector<double> &zsources
                      ,std::vector<std::vector<PixelMap> > &maps
                      ,double range
                      ,double angular_resolution
                      ,std::vector<Point_3d> &observers
                      ,std::vector<Point_3d> &directions
                      ,const std::vector<std::string> &snap_filenames
                      ,const std::vector<float> &snap_redshifts
                      ,double BoxLength
                      ,double particle_mass
                      ,bool verbose = false
                      ,bool addtocone = false
                      ){
    
    assert(cosmo.getOmega_matter() + cosmo.getOmega_lambda() == 1.0);
    // set coordinate distances for planes
    
    int Nmaps = zsources.size();    // number of source planes per cone
    int Ncones = observers.size();  // number of cones
    
    if(directions.size() != observers.size()){
      std::cerr << "Size of direction and observers must match." << std::endl;
      throw std::invalid_argument("");
    }
    
    double zs_max=0;
    for(auto z: zsources) zs_max = (z > zs_max) ? z : zs_max;
    
    Point_2d center;
    
    if(!addtocone){
      // allocate mamory for all maps
      maps.clear();
      maps.resize(Ncones);
      for(auto &map_v : maps){  // loop through cones
        map_v.reserve(Nmaps);
        for(int i=0;i<Nmaps;++i){
          map_v.emplace_back(center.x
                             ,(size_t)(range/angular_resolution)
                             ,(size_t)(range/angular_resolution)
                             ,angular_resolution);
        }
      }
    }else{
      /// check the sizes
      if(maps.size() != Ncones){
        std::cerr << "FastLightCones: You must have enough maps allocated or use addtocones = false " << std::endl;
        throw std::invalid_argument("Too few maps.");
      }
      for(auto &map_v : maps){  // loop through cones
        if(map_v.size() != Nmaps){
          std::cerr << "FastLightCones: You must have enough maps allocated or use addtocones = false " << std::endl;
          throw std::invalid_argument("Too few maps");
        }
      }
    }
    
    // find unique redshifts
    std::vector<double> z_unique(1,snap_redshifts[0]);
    for(auto z : snap_redshifts) if(z != z_unique.back() ) z_unique.push_back(z);
    std::sort(z_unique.begin(),z_unique.end());
    
    // find redshift ranges for each snapshot
    // shift the redshifts to between the snapshots
    std::vector<double> abins(z_unique.size() + 1);
    abins[0] = 1.0;
    for(int i=1;i<z_unique.size();++i){
      abins[i] = ( 1/(1+z_unique[i]) + 1/(1+z_unique[i-1]) )/2 ;
    }
    abins.back() = 1.0/( 1 + std::max(zs_max,z_unique.back() ) );
    
    std::vector<double> dbins(abins.size());
    for(int i=0 ; i<abins.size() ; ++i) dbins[i] = cosmo.coorDist(1.0/abins[i] - 1);
    
    std::vector<double> dsources(zsources.size());
    for(int i=0 ; i<Nmaps ; ++i) dsources[i] = cosmo.coorDist(zsources[i]);
    
    // make rotation Quaturnions to the observer frames
    std::vector<Quaternion> rotationQs(Ncones);
    for(int i = 0 ; i<Ncones ; ++i ){
      directions[i].unitize();
      SphericalPoint sp(directions[i]);
      rotationQs[i] = Quaternion::q_y_rotation(-sp.theta)*Quaternion::q_z_rotation(-sp.phi);
    }
    
    //const int blocksize = 1000000;  ????
    const int blocksize = 100000;
    
    //std::vector<T> points(blocksize);
    
    T unit;
    
    // loop through files
    for(int i_file=0 ; i_file < snap_filenames.size() ; ++i_file){
      
      // read header ??
      
      // find r range using snap_redshifts
      int i;
      for(i=0;i<z_unique.size();++i) if(snap_redshifts[i_file] == z_unique[i]) break;
      double dmin = dbins[i],dmax = dbins[i+1];
      
      // find the box range for each cone
      std::vector<Point_3d> max_box(Ncones),min_box(Ncones);
      for(int icone=0;icone<Ncones;++icone){
        for(int i=0;i<3;++i){
          
          double cos1,cos2;
          {
            double theta1,theta2;
            theta1 = acos(fabs(directions[icone][i])) - range/sqrt(2.);
            theta2 = pi - acos(fabs(directions[icone][i])) - range/sqrt(2.);
            
            if(theta1 > 0.0){
              cos1 = cos(theta1);
            }else cos1 = 1.0;
            
            if(theta2 > pi/2){
              cos2 = 0.0;
            }else if(theta2 < 0.0){
              cos2 = 1.0;
            }else{
              cos2 = cos(theta2);
            }
          }
          
          if(directions[icone][i] < 0.0) std::swap(cos1,cos2);
          
          max_box[icone][i] = (int)((observers[icone][i] + cos1*dmax)/BoxLength);
          double d2 = ((observers[icone][i] - cos2*dmax)/BoxLength);
          if(d2 > 0) min_box[icone][i] = 0;
          else min_box[icone][i] = (int)(d2) - 1;
          
        }
      }
      
      //open file
      if(verbose) std::cout <<" Opening " << snap_filenames[i_file] << std::endl;
      FILE *pFile = fopen(snap_filenames[i_file].c_str(),"r");
      if(pFile == nullptr){
        std::cout << "Can't open file " << snap_filenames[i_file] << std::endl;
        ERROR_MESSAGE();
        throw std::runtime_error(" Cannot open file.");
      }
      //points.resize(blocksize);
      size_t Nlines = 0;
      std::mutex clmoo;
      while(!feof(pFile)){  // loop through blocks
        
        //scan_block<T>(blocksize,points,pFile);
        //scan_block<Point_3d>(blocksize,points,pFile);
        
        long Nbatch = unit.scan_block(blocksize,pFile);
        
        if(Nbatch > 0){  // multi-thread the sorting into cones and projection onto planes
          int nthreads = Utilities::GetNThreads();
          int chunk_size;
          do{
            chunk_size =  unit.points.size()/nthreads;
            if(chunk_size == 0) nthreads /= 2;
          }while(chunk_size == 0);
          if(nthreads == 0) nthreads = 1;
          
          int remainder =  unit.points.size()%chunk_size;
          
          assert(nthreads*chunk_size + remainder == unit.points.size() );
          
          std::vector<std::thread> thr(nthreads);
          for(int ii =0; ii< nthreads ; ++ii){
            
            //std::cout << ii*chunk_size << " " << n << std::endl;
            
            thr[ii] = std::thread(unit.fastplanes_parallel
                                  ,unit.points.data() + ii*chunk_size
                                  ,unit.points.data() + (ii+1)*chunk_size + (ii==nthreads-1)*remainder
                                  ,cosmo,std::ref(max_box),std::ref(min_box)
                                  ,std::ref(observers),std::ref(rotationQs)
                                  ,std::ref(dsources),std::ref(maps)
                                  ,dmin,dmax,BoxLength,std::ref(clmoo));
          }
          for(int ii = 0; ii < nthreads ;++ii){ if(thr[ii].joinable() ) thr[ii].join();}
        }
        
        Nlines += blocksize;
        if(Nlines%1000000 == 0) std::cout << "=" << std::flush ;
        std::cout << std::endl;
      }
      fclose(pFile);
    }
    
    
    size_t Npixels = maps[0][0].size();
    for(int isource=0;isource<Nmaps;++isource){
      // renormalize map
      double norm = 4*pi*particle_mass*Grav/dsources[isource]/angular_resolution/angular_resolution;
      
      // calculate expected average kappa
      DkappaDz dkappadz(cosmo,zsources[isource]);
      double avekappa = Utilities::nintegrate<DkappaDz,double>(dkappadz,0,zsources[isource],1.0e-6);
      
      for(int icone=0 ; icone<Ncones ; ++icone){
        //for(auto &cone_maps: maps){
        maps[icone][isource] *= norm;
        // subtract average
        //double ave = cone_maps[i].ave();
        for(size_t ii=0 ; ii<Npixels ; ++ii)
          maps[icone][isource][ii] -= avekappa;
      }
    }
  }
  
}
#endif /* defined(__GLAMER__lightcone_construction__) */
