//
//  OTreeNB.h
// 3d octal tree
//  GLAMER
//
//  Created by Ben Metcalf on 5/2025.
//

#ifndef OTreeNB_h
#define OTreeNB_h

#include "utilities.h"
#include "cosmo.h"
#include "point.h"

#ifndef PI
#define PI 3.141593
#endif

/** type for particle positions and boundaries etc **/

#ifndef IndexType_declare
#define IndexType_declare
typedef unsigned long IndexType;
#endif


/** \brief
 * OTreeNB: Tree structure used for force calculation with particles (i.e. stars, Nbody or halos).
 *
 * The tree also contains pointers to the list of positions, sizes and masses of the particles.
 * Also flags for the number of dimensions the tree is defined in (2 or 3), and if multiple
 * masses and sizes should be used.
 * 
 * The template parameter PType is the structure 
 * used to store the particle positions.  
 * It must operators [0], [1], [2] and operator - defined.
 * Point_3d<float> for example
 */
template <typename PType = Point_3d<PosType>,typename DType = PosType >
struct OTreeNB
{

  ///  classes
  /** \brief Box representing a branch in a tree.  It has four children.  Used in OTreeNB which is used in TreeQuad.
 */
struct Branch
{
  Branch()
  {
    prev = nullptr;
    child = nullptr;
    brother = nullptr;
    particles = nullptr;
    nparticles = 0;
    level = -1;
    boxsize = 0.0;
    root_inv_density = 0.0;
    ++Nbranches;
  };

  int Nbranches;  ///< number of branches in tree, used for debugging
  //Branch(Branch *parent) : prev(parent)
  //{
  //  child = nullptr;
  //  brother = nullptr;
  //  particles = nullptr;
  //  nparticles = 0;
  //  level = parent->level + 1;
  //};

  ~Branch() {
    //std::cout << " deleting branch " << Nbranches << " "; 
    --Nbranches; 
    if(children != nullptr) delete[] children;
  };//{std::cout << "Branch destructor called" << std::endl; };

  /// array of particles in Branch
  IndexType *particles;
  IndexType nparticles;

  /// bottom, left, back corner of box
  PType boundary_p1;
  /// top, right, front corner of box
  PType boundary_p2;
  ///< center of box
  PType center;  

  PosType boxsize;  // one dimensional size of box
  PosType root_inv_density;  // cube root of number density in branch
  /// level in tree
  int level;
  //unsigned long number;

  // std::vector<Branch> children;
  Branch *children = nullptr;

  Branch *child;
  /// father of branch
  Branch *prev;
  /// Either child2 of father is branch is child1 and child2 exists or the brother of the father.
  /// Used for iterative tree walk.
  Branch *brother;

  // projected quantities   
  // center of mass  
  PosType xcm[3] = {0.0, 0.0, 0.0}; 
  // quadropole moment of branch in particle mass units
  PosType quad[3] = {0.0, 0.0, 0.0};
  PosType mass = 0;
};

  /**
   \brief A iterator class that allows for movement through the tree without changing
       anything in the tree itself.
 */
  class iterator
  {

  private:
    Branch *current;

  public:
    /// Sets the top or root to the top of "tree".
    iterator(OTreeNB<PType> *tree) { current = tree->top; }
    /// Sets the root to the input branch so that this will be a subtree in branch is not the real root.
    iterator(Branch *branch) { current = branch; }

    iterator& operator=(iterator &it)
    {
      current = it.current;
      return *this;
    }

    iterator& operator=(const Branch *branch)
    {
      current = branch;
      return *this;
    }

    /// Returns a pointer to the current Branch.
    Branch *operator*() { return current; }

    bool up()
    {
      if (current->prev == nullptr)
      {
        current = nullptr;
        return false; // at top
      }
      current = current->prev;
      return true;
    }

    /// Move to child
    bool down(short child)
    {
      assert(child >= 0 && child < 8);
      current = current->child;
    }

    bool atTop() const
    {
      if (current->prev == nullptr)
        return true;
      return false;
    }

    bool atLeaf() const
    {
      for (int i = 0; i < 8; ++i)
      {
        if (current->child != nullptr)
          return false;
      }
      return true;
    }

    /**  walk the tree, returns false if the next step would be to move to the 
     brother of the end iterator

     end = tree.begin() if it is going to the end of the tree
     allowDescent = true if it is allowed to go down to the child
                   , false if it should only go to the next brother
     */
    bool walk(bool allowDescent,iterator end)
    {
      if (current == nullptr)
        return false; // at top

      // if there is a valid chaild, go down
      if (allowDescent && current->child != nullptr){
        current = current->child;
        return true;
      }

      if (current->brother != end.current->brother)
        {
          current = current->brother;
          return true;
        }
        current = end.current;
        return false;  // at end
      }

  };
  

  OTreeNB(PType *xp,IndexType nparticles):
  xxp(xp)
  {
    index.resize(nparticles);
    for(unsigned long i=0;i<nparticles;++i) index[i] = i;
    Nbranches = 1;
    top.particles = index.data();
    //std::cout << top.particles[10] << std::endl;
    top.nparticles = nparticles;
     
    top.level = 0;
    top.boundary_p1[0] = top.boundary_p2[0] = xp[0][0];
    top.boundary_p1[1] = top.boundary_p2[1] = xp[0][1];
    top.boundary_p1[2] = top.boundary_p2[2] = xp[0][2];
  
    for(IndexType i = 1 ; i < nparticles ; ++i)
    {
      for(int j=0;j<3;++j){
        top.boundary_p1[j] = top.boundary_p1[j] < xp[i][j] ? top.boundary_p1[j] : xp[i][j];
        top.boundary_p2[j] = top.boundary_p2[j] > xp[i][j] ? top.boundary_p2[j] : xp[i][j];

        //top.boundary_p1[j] = min<PType>(top.boundary_p1[j], xp[i][j]);
        //top.boundary_p2[j] = max<PType>(top.boundary_p2[j], xp[i][j]);
      }
    }


    // make it into a cube
    PType length = top.boundary_p2 - top.boundary_p1;
    PosType max_length = length[0];
    int dim = 0;
    for(int j=1;j<3;++j){
      if(length[j] > max_length){
        max_length = length[j];
        dim = j;
      }
    }
    for(int j=0;j<3;++j){
      if(j != dim){
        top.boundary_p2[j] = top.boundary_p1[j] + max_length;
      }
    }

    top.center = (top.boundary_p1 + top.boundary_p2) * 0.5;
    top.boxsize =  max_length;
    top.root_inv_density = max_length/pow(top.nparticles,1.0/3.);
  }
  ~OTreeNB(){ };  // all the Branches are deleted in the destructor of Branch


  // returns number of branches in tree
  unsigned long size() { return Nbranches; };

  // build the tree down to bucket size, does not require the particles to have sizes
  void build(int Nbucket){ 
    auto it = begin();
    if((*it)->child != nullptr){
      throw std::runtime_error("OTreeNB Error: calling build() on already built tree");
    }
    span8(*it);
    assert(total_branches ==  Nbranches);
    while(it.walk(true,begin())){
     
      //std::cout << "level : " << (*it)->level << std::endl;
      //std::cout << (*it)->boundary_p1 << std::endl;
      //std::cout << (*it)->boundary_p2 << std::endl;
      //std::cout << "N : " << (*it)->nparticles << std::endl << std::endl;

      if((*it)->nparticles > Nbucket) span8(*it);
    }
    calcMoments();  // assumes all particles are same mass
  }

  // build the tree down to bucket size, does not require the particles to have sizes
  void build_size(int Nbucket){ 
    auto it = begin();
    if((*it)->child != nullptr){
      throw std::runtime_error("OTreeNB Error: calling build() on already built tree");
    }
    span8(*it);
    while(it.walk(true,begin())){
     
      // find largest particle in branch
      double max_size = 0.0;
      for(IndexType i=0;i<(*it)->nparticles;++i){
        max_size = max(max_size, xxp[(*it)->particles[i]].size);
      }
      //std::cout << "level : " << (*it)->level << std::endl;
      //std::cout << (*it)->boundary_p1 << std::endl;
      //std::cout << (*it)->boundary_p2 << std::endl;
      //std::cout << "N : " << (*it)->nparticles << std::endl << std::endl;
      if((*it)->nparticles > Nbucket && max_size < (*it)->boxsize ) span8(*it);
    }
  }

  // calculate the moments assuming the particles are equal mass and symmetric
  void calcMoments(){
    auto it = begin();
    while(it.walk(true,begin())){

      Branch *branch = *it;
      if(branch->nparticles == 0) continue;
      
      // calculate center of mass
      branch->xcm[0] = 0.0;
      branch->xcm[1] = 0.0;
      branch->xcm[2] = 0.0;
      for(IndexType i=0;i<branch->nparticles;++i){
        PType &x = xxp[branch->particles[i]];
        branch->xcm[0] += x[0];
        branch->xcm[1] += x[1];
        branch->xcm[2] += x[2];
      }
      branch->xcm[0] /= branch->nparticles;
      branch->xcm[1] /= branch->nparticles;
      branch->xcm[2] /= branch->nparticles;
      branch->mass = branch->nparticles;

      /** calculate quadropole moment of branch ?? needs to have mass
      Point_2d dxcm;
      branch->quad[0]=branch->quad[1]=branch->quad[2]=0;
      for(IndexType i=0;i<branch->nparticles;++i){
        PType &x = xxp[branch->particles[i]];
        dxcm[0] = x[0] - branch->xcm[0];
        dxcm[1] = x[1] - branch->xcm[1];

        double r2 = dxcm[0]*dxcm[0] + dxcm[1]*dxcm[1];
        branch->quad[0] += (r2-2*dxcm[0]*dxcm[0]);
        branch->quad[1] += (r2-2*dxcm[1]*dxcm[1]);
        branch->quad[2] += -2*dxcm[0]*dxcm[1];
      }*/
    }
  }

  // calculate the moments assuming the particles have different masses
  void calcMoments(DType *masses){
    moments_set = true;
    auto it = begin();
    while(it.walk(true,begin())){

      Branch *branch = *it;
      if(branch->nparticles == 0) continue;
      
      // calculate center of mass
      branch->xcm[0] = 0.0;
      branch->xcm[1] = 0.0;
      branch->xcm[2] = 0.0;
      branch->mass = 0;

      for(IndexType i=0;i<branch->nparticles;++i){
        IndexType j = branch->particles[i];
        PType &x = xxp[j];
        double m = masses[j];
        branch->xcm[0] += x[0]*m;
        branch->xcm[1] += x[1]*m;
        branch->xcm[2] += x[2]*m;
        branch->mass += m;
      }
      branch->xcm[0] /= branch->mass;
      branch->xcm[1] /= branch->mass;
      branch->xcm[2] /= branch->mass;

      /** calculate quadropole moment of branch
      Point_2d dxcm;
      branch->quad[0]=branch->quad[1]=branch->quad[2]=0;
      for(IndexType i=0;i<branch->nparticles;++i){
        IndexType j = branch->particles[i];
        PType &x = xxp[j];
        
        dxcm[0] = x[0] - branch->xcm[0];
        dxcm[1] = x[1] - branch->xcm[1];

        double r2 = dxcm[0]*dxcm[0] + dxcm[1]*dxcm[1];
        branch->quad[0] += (r2-2*dxcm[0]*dxcm[0]) * masses[j];
        branch->quad[1] += (r2-2*dxcm[1]*dxcm[1]) * masses[j];
        branch->quad[2] += -2*dxcm[0]*dxcm[1] * masses[i];
      }*/
    }
  }

  void force2D(const PosType *ray
      ,PosType particle_mass
      ,PosType smooth_factor
      ,PosType theta2   // opening angle in radians
      ,PosType inv_area // compensating negative mass area 
      ,PosType *alpha   // not zeroed
      ,KappaType *kappa // not zeroed
      ,KappaType *gamma // not zeroed
      ,KappaType *phi   // not zeroed
    ) {

      //alpha[0]=alpha[1]=gamma[0]=gamma[1]=gamma[2]=0.0;
      //*kappa=*phi=0.0;

      PosType xcm[2],r2cm;
      auto it = begin();
      bool decend = true;
      while(it.walk(decend,begin())){

        xcm[0] = (*it)->xcm[0] - ray[0];
        xcm[1] = (*it)->xcm[1] - ray[1];
        r2cm = xcm[0]*xcm[0] + xcm[1]*xcm[1];

        if( r2cm*theta2 > ((*it)->boxsize)*((*it)->boxsize) ){

          if(r2cm*inv_area < 1){
            //????? use moments
            double mass = (*it)->mass * particle_mass;
            double prefac = mass/r2cm/PI;
            double tmp = -( prefac - mass*inv_area);
            
            alpha[0] -= tmp*xcm[0];
            alpha[1] -= tmp*xcm[1];
            
            tmp = -2.0*prefac/r2cm;
            
            gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
            gamma[1] += xcm[0]*xcm[1]*tmp;
            
            *kappa -= mass*inv_area;
            *phi += (prefac*log(r2cm) - mass*inv_area)*r2cm*0.5;
          }
          decend = false;
        }else{
        
          if((*it)->child == nullptr){ // leaf
          
            for(IndexType i=0;i<(*it)->nparticles;++i){

              PType &x = xxp[(*it)->particles[i]];
              xcm[0] = x[0] - ray[0];
              xcm[1] = x[1] - ray[1];
              r2cm = xcm[0]*xcm[0] + xcm[1]*xcm[1];
              if(r2cm*inv_area < 1){
                double prefac = particle_mass /r2cm/PI;
                double tmp = -( prefac - particle_mass*inv_area);
             
                alpha[0] -= tmp*xcm[0];
                alpha[1] -= tmp*xcm[1];
            
                tmp = -2.0*prefac/r2cm;
            
                gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
                gamma[1] += xcm[0]*xcm[1]*tmp;
            
                *kappa -= particle_mass*inv_area;
                *phi += (prefac*log(r2cm)- particle_mass*inv_area)*r2cm*0.5;

                double scale = smooth_factor*(*it)->root_inv_density;
                if(r2cm < 4*scale*scale){
                  // cubic B-spline profile
                  b_spline_profile(
                    xcm
                    ,sqrt(r2cm)
                    ,particle_mass
                    ,scale
                    ,alpha,kappa,gamma,phi
                  );
                }
              }
            }
            decend = false;
          }else{
            // decend
            decend = true;
          }

        }
    }
  }
  void force2D(const PosType *ray
      ,DType *masses
      ,PosType smooth_factor
      ,PosType theta2   // opening angle in radians
      ,PosType inv_area // compensating negative mass area 
      ,PosType *alpha   // not zeroed
      ,KappaType *kappa // not zeroed
      ,KappaType *gamma // not zeroed
      ,KappaType *phi   // not zeroed
    ) {

      assert( moments_set );
      //alpha[0]=alpha[1]=gamma[0]=gamma[1]=gamma[2]=0.0;
      //*kappa=*phi=0.0;

      PosType xcm[2],r2cm;
      auto it = begin();
      bool decend = true;
      while(it.walk(decend,begin())){

        xcm[0] = (*it)->xcm[0] - ray[0];
        xcm[1] = (*it)->xcm[1] - ray[1];
        r2cm = xcm[0]*xcm[0] + xcm[1]*xcm[1];

        if( r2cm*theta2 > ((*it)->boxsize)*((*it)->boxsize) ){

          if(r2cm*inv_area < 1){
            //????? use moments
            DType &mass = (*it)->mass;
            double prefac = mass/r2cm/PI;
            double tmp = -( prefac - mass*inv_area);
            
            alpha[0] -= tmp*xcm[0];
            alpha[1] -= tmp*xcm[1];
            
            tmp = -2.0*prefac/r2cm;
            
            gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
            gamma[1] += xcm[0]*xcm[1]*tmp;
            
            *kappa -= mass*inv_area;
            *phi += (prefac*log(r2cm) - mass*inv_area)*r2cm*0.5;
          }
          decend = false;
        }else{
        
          if((*it)->child == nullptr){ // leaf
          
            for(IndexType i=0;i<(*it)->nparticles;++i){
              IndexType j = (*it)->particles[i];
              PType &x = xxp[j];
              DType &mass = masses[j];

              xcm[0] = x[0] - ray[0];
              xcm[1] = x[1] - ray[1];
              r2cm = xcm[0]*xcm[0] + xcm[1]*xcm[1];
              if(r2cm*inv_area < 1){
                double prefac = mass /r2cm/PI;
                double tmp = -( prefac - mass * inv_area);
             
                alpha[0] -= tmp*xcm[0];
                alpha[1] -= tmp*xcm[1];
            
                tmp = -2.0*prefac/r2cm;
            
                gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
                gamma[1] += xcm[0]*xcm[1]*tmp;
            
                *kappa -= mass * inv_area;
                *phi += (prefac*log(r2cm)- mass * inv_area)*r2cm*0.5;

                double scale = smooth_factor*(*it)->root_inv_density;
                if(r2cm < 4*scale*scale){
                  // cubic B-spline profile
                  b_spline_profile(
                    xcm
                    ,sqrt(r2cm)
                    ,mass
                    ,scale
                    ,alpha,kappa,gamma,phi
                  );
                }
              }
            }
            decend = false;
          }else{
            // decend
            decend = true;
          }

        }
    }
  }
  
  void force2D_hole(const PosType *ray
      ,PosType particle_mass
      ,PosType smooth_factor
      ,PosType theta2   // opening angle in radians
      ,PosType inv_area // compensating negative mass area 
      ,PosType rmin     // minimum radius for hole
      ,PosType *xo      // center of hole in Mpc
      ,PosType *alpha   // not zeroed
      ,KappaType *kappa // not zeroed
      ,KappaType *gamma // not zeroed
      ,KappaType *phi   // not zeroed
    ) {

      //alpha[0]=alpha[1]=gamma[0]=gamma[1]=gamma[2]=0.0;
      //*kappa=*phi=0.0;

      PosType xcm[2],r2cm,dx[2];

      auto it = begin();
      bool decend = true;
      while(it.walk(decend,begin())){

        dx[0] = (*it)->center[0] - xo[0];
        dx[1] = (*it)->center[1] - xo[1];
        PosType d = sqrt(dx[0]*dx[0] + dx[1]*dx[1]); // distance to center of branch
        
        if(d < rmin - (*it)->boxsize) {
          // branch is inside hole, skip it
          decend = false;
          continue;
        }

        xcm[0] = (*it)->xcm[0] - ray[0];
        xcm[1] = (*it)->xcm[1] - ray[1];
        r2cm = xcm[0]*xcm[0] + xcm[1]*xcm[1];

        if(d > rmin + (*it)->boxsize){
          if( r2cm*theta2 > ((*it)->boxsize)*((*it)->boxsize) ){

            if(r2cm*inv_area < 1){
              //????? use moments
              double mass = (*it)->nparticles * particle_mass;
              double prefac = mass/r2cm/PI;
              double tmp = -( prefac - mass*inv_area);
            
              alpha[0] -= tmp*xcm[0];
              alpha[1] -= tmp*xcm[1];
            
              tmp = -2.0*prefac/r2cm;
            
              gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
              gamma[1] += xcm[0]*xcm[1]*tmp;
            
              *kappa -= mass*inv_area;
              *phi += (prefac*log(r2cm) - mass*inv_area)*r2cm*0.5;
            }
            decend = false;
          }else{
        
            if((*it)->child == nullptr){ // leaf
          
              for(IndexType i=0;i<(*it)->nparticles;++i){

                PType &x = xxp[(*it)->particles[i]];
                xcm[0] = x[0] - ray[0];
                xcm[1] = x[1] - ray[1];
                r2cm = xcm[0]*xcm[0] + xcm[1]*xcm[1];
                if(r2cm*inv_area < 1){
                  double prefac = particle_mass /r2cm/PI;
                  double tmp = -( prefac - particle_mass*inv_area);
             
                  alpha[0] -= tmp*xcm[0];
                  alpha[1] -= tmp*xcm[1];
            
                  tmp = -2.0*prefac/r2cm;
            
                  gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
                  gamma[1] += xcm[0]*xcm[1]*tmp;
            
                  *kappa -= particle_mass*inv_area;
                  *phi += (prefac*log(r2cm)- particle_mass*inv_area)*r2cm*0.5;

                  double scale = smooth_factor*(*it)->root_inv_density;
                  if(r2cm < 4*scale*scale){
                    // cubic B-spline profile
                    b_spline_profile(
                      xcm
                      ,sqrt(r2cm)
                      ,particle_mass
                      ,scale
                      ,alpha,kappa,gamma,phi
                    );
                  }
               }
             }
             decend = false;
            }else{
              // decend
              decend = true;
            }
          }
        }else{
          // possible particle in and outside of hole
            if((*it)->child == nullptr){ // leaf
          
              for(IndexType i=0;i<(*it)->nparticles;++i){

                PType &x = xxp[(*it)->particles[i]];
                xcm[0] = x[0] - ray[0];
                xcm[1] = x[1] - ray[1];
                r2cm = xcm[0]*xcm[0] + xcm[1]*xcm[1];
                
                dx[0] = x[0] - xo[0];
                dx[1] = x[1] - xo[1];

                if( dx[0]*dx[0] + dx[1]*dx[1] > rmin*rmin &&
                  r2cm*inv_area < 1){
                  double prefac = particle_mass /r2cm/PI;
                  double tmp = -( prefac - particle_mass*inv_area);
             
                  alpha[0] -= tmp*xcm[0];
                  alpha[1] -= tmp*xcm[1];
            
                  tmp = -2.0*prefac/r2cm;
            
                  gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
                  gamma[1] += xcm[0]*xcm[1]*tmp;
            
                  *kappa -= particle_mass*inv_area;
                  *phi += (prefac*log(r2cm)- particle_mass*inv_area)*r2cm*0.5;

                  double scale = smooth_factor*(*it)->root_inv_density;
                  if(r2cm < 4*scale*scale){
                    // cubic B-spline profile
                    b_spline_profile(
                      xcm
                      ,sqrt(r2cm)
                      ,particle_mass
                      ,scale
                      ,alpha,kappa,gamma,phi
                    );
                  }
                }
              }
              decend = false;
            }else{
              // decend
              decend = true;
            }
        }
      }
  }

    void force2D_hole(
      const PosType *ray
      ,PType *masses
      ,PosType smooth_factor
      ,PosType theta2   // opening angle in radians
      ,PosType inv_area // compensating negative mass area 
      ,PosType rmin     // minimum radius for hole
      ,PosType *xo      // center of hole in Mpc
      ,PosType *alpha   // not zeroed
      ,KappaType *kappa // not zeroed
      ,KappaType *gamma // not zeroed
      ,KappaType *phi   // not zeroed
    ) {

      assert( moments_set );
      //alpha[0]=alpha[1]=gamma[0]=gamma[1]=gamma[2]=0.0;
      //*kappa=*phi=0.0;

      PosType xcm[2],r2cm,dx[2];

      auto it = begin();
      bool decend = true;
      while(it.walk(decend,begin())){

        dx[0] = (*it)->center[0] - xo[0];
        dx[1] = (*it)->center[1] - xo[1];
        PosType d = sqrt(dx[0]*dx[0] + dx[1]*dx[1]); // distance to center of branch
        
        if(d < rmin - (*it)->boxsize) {
          // branch is inside hole, skip it
          decend = false;
          continue;
        }

        xcm[0] = (*it)->xcm[0] - ray[0];
        xcm[1] = (*it)->xcm[1] - ray[1];
        r2cm = xcm[0]*xcm[0] + xcm[1]*xcm[1];

        if(d > rmin + (*it)->boxsize){
          if( r2cm*theta2 > ((*it)->boxsize)*((*it)->boxsize) ){

            if(r2cm*inv_area < 1){
              //????? use moments
              double &mass = (*it)->mass;
              double prefac = mass/r2cm/PI;
              double tmp = -( prefac - mass*inv_area);
            
              alpha[0] -= tmp*xcm[0];
              alpha[1] -= tmp*xcm[1];
            
              tmp = -2.0*prefac/r2cm;
            
              gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
              gamma[1] += xcm[0]*xcm[1]*tmp;
            
              *kappa -= mass*inv_area;
              *phi += (prefac*log(r2cm) - mass*inv_area)*r2cm*0.5;
            }
            decend = false;
          }else{
        
            if((*it)->child == nullptr){ // leaf
          
              for(IndexType i=0;i<(*it)->nparticles;++i){
                IndexType j = (*it)->particles[i];
                PType &x = xxp[j];
                xcm[0] = x[0] - ray[0];
                xcm[1] = x[1] - ray[1];
                r2cm = xcm[0]*xcm[0] + xcm[1]*xcm[1];
                if(r2cm*inv_area < 1){
                  double &mass = masses[j];
                  double prefac = mass /r2cm/PI;
                  double tmp = -( prefac - mass * inv_area);
             
                  alpha[0] -= tmp*xcm[0];
                  alpha[1] -= tmp*xcm[1];
            
                  tmp = -2.0*prefac/r2cm;
            
                  gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
                  gamma[1] += xcm[0]*xcm[1]*tmp;
            
                  *kappa -= mass * inv_area;
                  *phi += (prefac*log(r2cm) - mass * inv_area)*r2cm*0.5;

                  double scale = smooth_factor*(*it)->root_inv_density;
                  if(r2cm < 4*scale*scale){
                    // cubic B-spline profile
                    b_spline_profile(
                      xcm
                      ,sqrt(r2cm)
                      ,mass
                      ,scale
                      ,alpha,kappa,gamma,phi
                    );
                  }
               }
             }
             decend = false;
            }else{
              // decend
              decend = true;
            }
          }
        }else{
          // possible particle in and outside of hole
            if((*it)->child == nullptr){ // leaf
          
              for(IndexType i=0;i<(*it)->nparticles;++i){
                IndexType j = (*it)->particles[i];
                PType &x = xxp[j];
                xcm[0] = x[0] - ray[0];
                xcm[1] = x[1] - ray[1];
                r2cm = xcm[0]*xcm[0] + xcm[1]*xcm[1];
                
                dx[0] = x[0] - xo[0];
                dx[1] = x[1] - xo[1];

                if( dx[0]*dx[0] + dx[1]*dx[1] > rmin*rmin &&
                  r2cm*inv_area < 1){
                  double &mass = masses[j];
                  double prefac = mass /r2cm/PI;
                  double tmp = -( prefac - mass * inv_area);
             
                  alpha[0] -= tmp*xcm[0];
                  alpha[1] -= tmp*xcm[1];
            
                  tmp = -2.0*prefac/r2cm;
            
                  gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
                  gamma[1] += xcm[0]*xcm[1]*tmp;
            
                  *kappa -= mass * inv_area;
                  *phi += (prefac*log(r2cm)- mass * inv_area)*r2cm*0.5;

                  double scale = smooth_factor*(*it)->root_inv_density;
                  if(r2cm < 4*scale*scale){
                    // cubic B-spline profile
                    b_spline_profile(
                      xcm
                      ,sqrt(r2cm)
                      ,mass
                      ,scale
                      ,alpha,kappa,gamma,phi
                    );
                  }
                }
              }
              decend = false;
            }else{
              // decend
              decend = true;
            }
        }
      }
  }
  size_t getUsedBranches() { return Nbranches; }
  size_t getTotalBranches() { return total_branches; }

  iterator begin()
  {
    return iterator(&top);
  }
 
  // maximum depth of tree
  int getDepth() const
  {
    return depth;
  }

  void getBoundingBox(
    PType &p1,  ///< bottom, left, back corner of box
    PType &p2   ///< top, right, front corner of box
  ) const {
    p1 = top.boundary_p1;
    p2 = top.boundary_p2;
  }
private:

  void span8(Branch *current);
 
  Branch top;
  /// Array of particle positions
  PType *xxp;
  std::vector<unsigned long> index;

  /// number of branches in tree
  unsigned long Nbranches;
  unsigned long total_branches = 1;
  int depth = 1; // maximum depth of tree, used for debugging
  bool moments_set = false;

  PosType phiintconst = (120*log(2.) - 631.)/840 + 19./70; // ????

  /* cubic B-spline kernel for particle profile
   
   The lensing quantities are added to and a point mass is subtracted
   */
  inline void b_spline_profile(
                               PosType *xcm       // vector in Mpc connecting ray to center of particle
                               ,PosType r         // distance from center in Mpc
                               ,PosType Mass      // mass in solar masses
                               ,PosType size      // size scale in Mpc
                               ,PosType *alpha    // deflection angle times Sigma_crit
                               ,KappaType *kappa  // surface density
                               ,KappaType *gamma  // shear times Sigma_crit
                               ,KappaType *phi
                               ) const {
    
    
    PosType q = r/size;
    PosType M,sigma;
    if(q > 2){
      sigma = 0;
      M = 1;
    }else{
      PosType q2=q*q,q3=q2*q,q4=q2*q2,q5=q4*q;
      
      sigma = (8 - 12*q + 6*q2 - q3)/4;
      if(q > 1){
        sigma *= 10/size/size/7/PI;
        M = (-1 + 20*q2*(1 - q + 3*q2/8 - q3/20) )/7;
        *phi += Mass*(-1232. + 1200*q2 - 800.*q3 + 225.*q4 - 24*q5 + 120*log(2./q) )/840/PI;
      }else{
        sigma = 10*( sigma - 1 + 3*q - 3*q2 + q3)/size/size/7/PI;
        M = 10*q2*(1 - 3*q/4 + 3*q3/10)/7;
        
        *phi += Mass*( phiintconst + 10*(q2/2 - 3*q4/4 + 3*q5/50)/7
                      )/PI;
      }
    }
    
    PosType alpha_r,gt;  // deflection * Sig_crit / Mass
    alpha_r = (M-1)/PI/r;
    gt = alpha_r/r - sigma;
    
    alpha[0] += Mass*alpha_r*xcm[0]/r;
    alpha[1] += Mass*alpha_r*xcm[1]/r;
    gamma[0] -= gt*Mass*(xcm[0]*xcm[0]-xcm[1]*xcm[1])/r/r;
    gamma[1] -= gt*Mass*2*xcm[0]*xcm[1]/r/r;
    *kappa += Mass*sigma;
    *phi -= Mass*log(r)/PI;
  }
  
};

// creates 8 children to current with links to their parent and brother
// and sets their center and boxsize
template <typename PType,typename DType>
void OTreeNB<PType,DType>::span8(Branch *current)
{
  if (current == nullptr)
  {
    ERROR_MESSAGE();
    std::cerr << "OTreeNB Error: calling spon() on empty tree" << std::endl;
    exit(1);
  }

  current->children = new Branch[8];
  Branch *children = current->children;

  for(int i=0; i<8; ++i)
  {
    children[i].prev = current;
    children[i].level = current->level + 1;
  }
  total_branches += 8;
  if(depth < current->level + 1) depth = current->level + 1;

  PType &center = current->center;
  //center[0] = (current->boundary_p1[0] + current->boundary_p2[0]) / 2.0;
  //center[1] = (current->boundary_p1[1] + current->boundary_p2[1]) / 2.0;
  //center[2] = (current->boundary_p1[2] + current->boundary_p2[2]) / 2.0;
  PosType boxsize = current->boxsize / 2.0;

  size_t m = 0;
  for (int i = -1; i < 2; i += 2)
  {
    for (int j = -1; j < 2; j += 2)
    {
      for (int k = -1; k < 2; k += 2)
      {

        children[m].boundary_p1[0] = ((1 - i) / 2) * current->boundary_p1[0] 
                                              + ((1 + i) / 2) * center[0];
        children[m].boundary_p1[1] = ((1 - j) / 2) * current->boundary_p1[1] 
                                              + ((1 + j) / 2) * center[1];
        children[m].boundary_p1[2] = ((1 - k) / 2) * current->boundary_p1[2] 
                                              + ((1 + k) / 2) * center[2];

        children[m].boundary_p2[0] = ((1 + i) / 2) * current->boundary_p2[0]
                                              + ((1 - i) / 2) * center[0];
        children[m].boundary_p2[1] = ((1 + j) / 2) * current->boundary_p2[1]
                                              + ((1 - j) / 2) * center[1];
        children[m].boundary_p2[2] = ((1 + k) / 2) * current->boundary_p2[2]
                                              + ((1 - k) / 2) * center[2];

        children[m].center = (children[m].boundary_p1 + children[m].boundary_p2) * 0.5;

        children[m].boxsize = boxsize;
        ++m;
      }
    }
  }

  // sort particles into children
  IndexType *p = current->particles;
  //std::cout << p[0] << " " << p[1] << " " << p[2] << std::endl;
  IndexType *p_end = p + current->nparticles;
  IndexType np = 0;
  for (int m = 0; m < 7; ++m)
  {
    children[m].particles = p;
    IndexType *pp = p;
    //std::cout << pp[0] << " " << pp[1] << " " << pp[2] << std::endl;
    //std::cout << xxp[pp[0]][1] << " " << xxp[*pp][1] << " " << pp[2] << std::endl;
  
    while (pp != p_end)
    {
      PType &x = xxp[*pp];
      if (x[0] >= children[m].boundary_p1[0] &&
          x[0] < children[m].boundary_p2[0] &&
          x[1] >= children[m].boundary_p1[1] &&
          x[1] < children[m].boundary_p2[1] &&
          x[2] >= children[m].boundary_p1[2] &&
          x[2] < children[m].boundary_p2[2])
      {
        std::swap(*p, *pp);
        children[m].nparticles++;
        ++p;
      }
      ++pp;
    }
    np += children[m].nparticles;
    
    children[m].root_inv_density = children[m].boxsize/pow(children[m].nparticles,1.0/3.); // reset density
  }
  children[7].particles = p;
  children[7].nparticles = current->nparticles - np;
  children[7].root_inv_density = children[7].boxsize/pow(children[7].nparticles,1.0/3.);
  // remove empty branches from the brotherhood
  int i=0;
  while(children[i].nparticles==0 && i<8) ++i;
  assert(i != 8);
  current->child = children+i;
  ++Nbranches;
  int j=i+1;
  while(j<8){
    if(children[j].nparticles > 0){
      children[i].brother = children + j;
      i=j;
      ++Nbranches;
    }
    ++j; 
  }
  children[i].brother = current->brother;
}

#endif /* OTreeNB_h */
