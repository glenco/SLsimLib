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

/** \brief Box representing a branch in a tree.  It has four children.  Used in OTreeNB which is used in TreeQuad.
 */
struct OBranchNB
{
  OBranchNB()
  {
    prev = nullptr;
    child = nullptr;
    brother = nullptr;
    particles = nullptr;
    nparticles = 0;
    level = -1;
    boxsize = 0.0;
    root_density = 0.0;
  };

  //OBranchNB(OBranchNB *parent) : prev(parent)
  //{
  //  child = nullptr;
  //  brother = nullptr;
  //  particles = nullptr;
  //  nparticles = 0;
  //  level = parent->level + 1;
  //};

  ~OBranchNB() {};

  /// array of particles in OBranchNB
  IndexType *particles;
  IndexType nparticles;

  /// bottom, left, back corner of box
  Point_3d<PosType> boundary_p1;
  /// top, right, front corner of box
  Point_3d<PosType> boundary_p2;

  PosType boxsize;  // one dimensional size of box
  PosType root_density;  // cube root of number density in branch
  /// level in tree
  int level;
  //unsigned long number;

  // std::vector<OBranchNB> children;
  std::unique_ptr<OBranchNB> children;

  OBranchNB *child;
  /// father of branch
  OBranchNB *prev;
  /// Either child2 of father is branch is child1 and child2 exists or the brother of the father.
  /// Used for iterative tree walk.
  OBranchNB *brother;

  // projected quantities   
  // center of mass  
  PosType xcm[3] = {0.0, 0.0, 0.0}; 
  // quadropole moment of branch in particle mass units
  PosType quad[3] = {0.0, 0.0, 0.0};
};

/** \brief
 * OTreeNB: Tree structure used for force calculation with particles (i.e. stars, Nbody or halos).
 *
 * The tree also contains pointers to the list of positions, sizes and masses of the particles.
 * Also flags for the number of dimensions the tree is defined in (2 or 3), and if multiple
 * masses and sizes should be used.
 */
template <typename PType = double *>
struct OTreeNB
{
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
        top.boundary_p1[j] = min(top.boundary_p1[j], xp[i][j]);
        top.boundary_p2[j] = max(top.boundary_p2[j], xp[i][j]);
      }
    }

    // make it into a cube
    Point_3d<PosType> length = top.boundary_p2 - top.boundary_p1;
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

    top.boxsize =  max_length;
    top.root_density = pow(top.nparticles,1.0/3.)/max_length;
  }
  ~OTreeNB(){};

  // returns number of branches in tree
  unsigned long size() { return Nbranches; };

  // build the tree down to bucket size, does not require the particles to have sizes
  void build(int Nbucket){ 
    auto it = begin();
    if((*it)->child != nullptr){
      throw std::runtime_error("OTreeNB Error: calling build() on already built tree");
    }
    span8(*it);
    while(it.walk(true,begin())){
     
      std::cout << "level : " << (*it)->level << std::endl;
      std::cout << (*it)->boundary_p1 << std::endl;
      std::cout << (*it)->boundary_p2 << std::endl;
      std::cout << "N : " << (*it)->nparticles << std::endl << std::endl;

      if((*it)->nparticles > Nbucket) span8(*it);
    }
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
  void calcMoments_point(double particle_mass){
    auto it = begin();
    while(it.walk(true,begin())){

      OBranchNB *branch = *it;
      if(branch->nparticles == 0) continue;
      
      // calculate center of mass
      branch->xcm[0] = 0.0; // reset center
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

      // calculate quadropole moment of branch
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
      }
      branch->quad[0] *= particle_mass;
      branch->quad[1] *= particle_mass;
      branch->quad[2] *= particle_mass;
    }
  }


  /**
   \brief A iterator class that allows for movement through the tree without changing
       anything in the tree itself.
 */
  class iterator
  {

  private:
    OBranchNB *current;

  public:
    /// Sets the top or root to the top of "tree".
    iterator(OTreeNB<PType> *tree) { current = tree->top; }
    /// Sets the root to the input branch so that this will be a subtree in branch is not the real root.
    iterator(OBranchNB *branch) { current = branch; }

    
    iterator& operator=(iterator &it)
    {
      current = it.current;
      return *this;
    }
  

    iterator& operator=(OBranchNB *branch)
    {
      current = branch;
      return *this;
    }

    /// Returns a pointer to the current Branch.
    OBranchNB *operator*() { return current; }

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

  size_t getNbranches() { return Nbranches; }

  iterator begin()
  {
    return iterator(&top);
  }

  void span8(OBranchNB *current);
  
  // maximum depth of tree
  int getDepth() const
  {
    return depth;
  }
  
private:
  OBranchNB top;
  /// Array of particle positions
  PType *xxp;
  std::vector<unsigned long> index;

  /// number of branches in tree
  unsigned long Nbranches;
  unsigned long total_branches = 1;
  int depth = 1; // maximum depth of tree, used for debugging
};

// creates 8 children to current with links to their parent and brother
// and sets their center and boxsize
template <typename PType>
void OTreeNB<PType>::span8(OBranchNB *current)
{
  if (current == nullptr)
  {
    ERROR_MESSAGE();
    std::cerr << "OTreeNB Error: calling spon() on empty tree" << std::endl;
    exit(1);
  }

  current->children.reset(new OBranchNB[8]);
  OBranchNB *children = current->children.get();

  for(int i=0; i<8; ++i)
  {
    children[i].prev = current;
    children[i].level = current->level + 1;
  }
  total_branches += 8;
  if(depth < current->level + 1) depth = current->level + 1;

  PosType center[3];
  center[0] = (current->boundary_p1[0] + current->boundary_p2[0]) / 2.0;
  center[1] = (current->boundary_p1[1] + current->boundary_p2[1]) / 2.0;
  center[2] = (current->boundary_p1[2] + current->boundary_p2[2]) / 2.0;
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
    
    children[m].root_density = pow(children[m].nparticles,1.0/3.)/children[m].boxsize; // reset density
  }
  children[7].particles = p;
  children[7].nparticles = current->nparticles - np;
  children[7].root_density = pow(children[7].nparticles,1.0/3.)/children[7].boxsize; // reset root_density
  
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
