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

/**  \brief The LightCones namespace is for classes and functions related to making light cones and weak lensing maps from 3D snapshots of particles or halos.
 */
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
  
  
  using Utilities::Geometry::Quaternion;
  using Utilities::Geometry::SphericalPoint;
  
  /*************************************************************************************************
   templated functions for projecting mass onto planes
   *************************************************************************************************/
  
  
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
  
  /*************************************************************************************************
   structures to implement templated FastLightCones<>()
   *************************************************************************************************/
  
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
                                    ,std::vector<std::vector<Point_3d> > &boxes
                                    ,std::vector<Point_3d> &observers
                                    ,std::vector<Quaternion> &rotationQs
                                    ,std::vector<double> &dsources
                                    ,std::vector<std::vector<PixelMap> > &maps
                                    ,double dmin
                                    ,double dmax
                                    ,double BoxLength
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
                                    ,std::vector<std::vector<Point_3d> > &boxes
                                    ,std::vector<Point_3d> &observers
                                    ,std::vector<Quaternion> &rotationQs
                                    ,std::vector<double> &dsources
                                    ,std::vector<std::vector<PixelMap> > &maps
                                    ,double dmin
                                    ,double dmax
                                    ,double BoxLength
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
        if(points[i].mass > 0) ++i;
      points.resize(i);
      
      return i;
    }
    static void fastplanes_parallel(
                                    DatumXMR *begin
                                    ,DatumXMR *end
                                    ,const COSMOLOGY &cosmo
                                    ,std::vector<std::vector<Point_3d> > &boxes
                                    ,std::vector<Point_3d> &observers
                                    ,std::vector<Quaternion> &rotationQs
                                    ,std::vector<double> &dsources
                                    ,std::vector<std::vector<PixelMap> > &maps
                                    ,double dmin
                                    ,double dmax
                                    ,double BoxLength
                                    );
    
  };
  
  
  struct ASCII_XMRRT{
    
    std::vector<DatumXMRmRs> points;
    size_t scan_block(size_t blocksize,FILE *pFile){
      
      // read in a block of points
      size_t i=0;
      int tmp;
      
      points.resize(blocksize);
      while(i < blocksize &&
            fscanf(pFile,"%lf %lf %lf %lf %lf %lf %i"
                   ,&points[i].x[0],&points[i].x[1],&points[i].x[2]
                   ,&points[i].mass,&points[i].r_max,&points[i].r_scale,&tmp) != EOF)
        ++i;
      points.resize(i);
      
      for(auto &h: points){
        h.r_max /= 1.0e3;
        h.r_scale /= 1.0e3;
      }
      return i;
    }
    
    static void fastplanes_parallel(
                                    DatumXMRmRs *begin
                                    ,DatumXMRmRs *end
                                    ,const COSMOLOGY &cosmo
                                    ,std::vector<std::vector<Point_3d> > &boxes
                                    ,std::vector<Point_3d> &observers
                                    ,std::vector<Quaternion> &rotationQs
                                    ,std::vector<double> &dsources
                                    ,std::vector<std::vector<PixelMap> > &maps
                                    ,double dmin
                                    ,double dmax
                                    ,double BoxLength
                                    );
    
    
  };
  
  struct ASCII_XMRRT12:public ASCII_XMRRT{
    size_t scan_block(size_t blocksize,FILE *pFile){
      
      // read in a block of points
      size_t i=0;
      int tmp;
      
      points.resize(blocksize);
      while(i < blocksize &&
            fscanf(pFile,"%lf %lf %lf %lf %lf %lf %i"
                   ,&points[i].x[0],&points[i].x[1],&points[i].x[2]
                   ,&points[i].mass,&points[i].r_max,&points[i].r_scale,&tmp) != EOF){
              if(tmp == 1 || tmp == 2) ++i;
            }
      points.resize(i);
      
      for(auto &h: points){
        h.r_max /= 1.0e3;
        h.r_scale /= 1.0e3;
      }
      return i;
    }
    
  };
  
  /**************************************************************************/
  /**************************************************************************/
  /**************************************************************************/
  
  void random_observers(std::vector<Point_3d> &observers
                        ,std::vector<Point_3d> &directions
                        ,int Ncones
                        ,double BoxLength
                        ,double cone_opening_radius
                        ,Utilities::RandomNumbers_NR &ran
                        );  
  
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
  
  /** \brief This function goes directly from snapshots to lensing maps with Born approximation and linear propogations.
   
   <p>
   The template parameter allows this function to use different input data formats and
   different ways of distributing the mass into grid cells.  The current options are
   
   LightCones::ASCII_XV   -- for 6 column ASCII file with position and velocity
   
   LightCones::ASCII_XM   -- for 5 column ASCII file with position and mass
   
   LightCones::ASCII_XMR   -- for 6 column ASCII file with position, mass and size
   
   LightCones::ASCII_XMRRT  -- for 7 column ASCII file with position, mass,Rmax,Rscale and an integer denoting type, the halos are rendered as truncated NFWs
   
   LightCone::ASCII_XMRRT12 -- same as LightCones::ASCII_XMRRT but only takes entries with the 7th column equal to 1 or 2
   <\p>
   
   */
  
  template <typename T>
  void FastLightCones(
                      const COSMOLOGY &cosmo
                      ,const std::vector<double> &zsources   /// vector of source redshifts
                      ,std::vector<std::vector<PixelMap> > &maps  /// output kappa maps, do not need to be allocated initially
                      ,double range                           /// angular range in radians of maps
                      ,double angular_resolution              /// angular resolution in radians of maps
                      ,std::vector<Point_3d> &observers     /// position of observers within the simulation box coordinates 0 < x < Lengths
                      ,std::vector<Point_3d> &directions      /// direction of light cones
                      ,const std::vector<std::string> &snap_filenames  /// names of the data files
                      ,const std::vector<float> &snap_redshifts  /// the redshift of the data files
                      ,double BoxLength
                      ,double mass_units
                      ,double subtract_mean = true /// subtracts expected mean kappa from map, otherwise the average density of the universe is not subtracted
                      ,bool verbose = false
                      ,bool addtocone = false  /// if false the maps will be cleared and new maps made, if true particles are added to the existing maps
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
    
    /// make a seporate set of maps for each thread
    std::vector<std::vector<std::vector<PixelMap> > > map_pack(Utilities::GetNThreads());
    for(auto &maps : map_pack){
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
    
    const int blocksize = 1000000;
    T unit;
    
    // loop through files
    for(int i_file=0 ; i_file < snap_filenames.size() ; ++i_file){
      
      // read header ??
      
      // find r range using snap_redshifts
      int i;
      for(i=0;i<z_unique.size();++i) if(snap_redshifts[i_file] == z_unique[i]) break;
      double dmin = dbins[i],dmax = dbins[i+1];
      
      // find the box range for each cone
      std::vector<std::vector<Point_3d> > boxes(Ncones);
      for(int icone=0;icone<Ncones;++icone){
        Point_3d max_box,min_box;
        for(int i=0;i<3;++i){
          
          {  // dumb box range the contains total sphere of radius dmax
            max_box[i] = (int)( (observers[icone][i] + dmax)/BoxLength );
            if(observers[icone][i] < dmax){
              min_box[i] = (int)( (observers[icone][i] - dmax)/BoxLength ) - 1;
            }else{
              min_box[i] = 0;
            }
          }
        }
        
        Point_3d n;
        Utilities::Geometry::Cone cone(observers[icone],directions[icone],range/sqrt(2));
        for(n[0] = min_box[0] ; n[0] <= max_box[0] ; ++n[0]){
          for(n[1] = min_box[1] ; n[1] <= max_box[1] ; ++n[1]){
            for(n[2] = min_box[2] ; n[2] <= max_box[2] ; ++n[2]){
              
              Point_3d p1(BoxLength*n[0],BoxLength*n[1],BoxLength*n[2]);
              Point_3d p2(BoxLength*(n[0]+1),BoxLength*(n[1]+1),BoxLength*(n[2]+1));
              
              //if( cone.intersect_box(p1,p2) ) boxes[icone].push_back(n);
              if( cone.intersect_box(p1,p2) ){
                
                // require that at least one corner is outside the sphere R=dmin around observer
                p1 = p1 - observers[icone];
                p2 = p2 - observers[icone];
                if(p1.length() > dmin || p2.length() > dmin ||
                   (p1 + Point_3d(BoxLength,0,0)).length() > dmin ||
                   (p1 + Point_3d(0,BoxLength,0)).length() > dmin ||
                   (p1 + Point_3d(0,0,BoxLength)).length() > dmin ||
                   (p1 + Point_3d(0,BoxLength,BoxLength)).length() > dmin ||
                   (p1 + Point_3d(BoxLength,0,BoxLength)).length() > dmin ||
                   (p1 + Point_3d(BoxLength,BoxLength,0)).length() > dmin
                   )  boxes[icone].push_back(n);
              }
              
            }
          }
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
      
      size_t Nlines = 0,Nblocks=0;
      while(!feof(pFile)){  // loop through blocks
        
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
            
            thr[ii] = std::thread(&unit.fastplanes_parallel
                                  ,unit.points.data() + ii*chunk_size
                                  ,unit.points.data() + (ii+1)*chunk_size + (ii==nthreads-1)*remainder
                                  ,cosmo,std::ref(boxes)
                                  ,std::ref(observers),std::ref(rotationQs)
                                  ,std::ref(dsources),std::ref(map_pack[ii])
                                  ,dmin,dmax,BoxLength);
          }
          for(int ii = 0; ii < nthreads ;++ii){ if(thr[ii].joinable() ) thr[ii].join();}
        }
        
        Nlines += blocksize;
        ++Nblocks;
        std::cout << "=" << std::flush ;
        if(Nblocks%10 == 0) std::cout << " " << std::flush;
        if(Nblocks%50 == 0) std::cout << std::endl;
        if(Nblocks%100 == 0) std::cout << std::endl;
        //std::cout << std::endl;
        
      }
      fclose(pFile);
      std::cout << std::endl;
    }
    
    // add maps from different threads together
    for(auto &tmaps : map_pack){
        for(int isource=0;isource<Nmaps;++isource){
          for(int icone=0 ; icone<Ncones ; ++icone){
            maps[icone][isource] += tmaps[icone][isource];
        }
      }
    }

    
    size_t Npixels = maps[0][0].size();
    //double norm = 4*pi*mass_units*Grav*cosmo.gethubble()*cosmo.gethubble()/angular_resolution/angular_resolution;
    double norm = 4*pi*mass_units*Grav/angular_resolution/angular_resolution;

    // normalize maps
    for(int isource=0;isource<Nmaps;++isource){
      // renormalize map
      // ??? need to put factors of hubble parameters in ?
      double norm2 = norm/dsources[isource];
      
      for(int icone=0 ; icone<Ncones ; ++icone){
        //for(auto &cone_maps: maps){
        maps[icone][isource] *= norm2;
      }
    }
    
    // subtract average kappa
    
    if(subtract_mean){
      for(int isource=0;isource<Nmaps;++isource){
        // calculate expected average kappa
        DkappaDz dkappadz(cosmo,zsources[isource]);
        double avekappa = Utilities::nintegrate<DkappaDz,double>(dkappadz,0,zsources[isource],1.0e-6);
      
        for(int icone=0 ; icone<Ncones ; ++icone){
          for(size_t ii=0 ; ii<Npixels ; ++ii)
            maps[icone][isource][ii] -= avekappa;
        }
      }
    }
  }
  
  
}
#endif /* defined(__GLAMER__lightcone_construction__) */
