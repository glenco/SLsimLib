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
  OBranchNB(OBranchNB *parent) : prev(parent)
  {
    static unsigned long n = 0;

    child = NULL;
    brother = NULL;
    particles = NULL;
    nparticles = 0;
    number = n;
    ++n;
  };

  ~OBranchNB() {};

  /// array of particles in OBranchNB
  IndexType *particles;
  IndexType nparticles;

  /// bottom, left, back corner of box
  PosType boundary_p1[3];
  /// top, right, front corner of box
  PosType boundary_p2[3];

  bool cubic = false;
  PosType boxsize;
  PosType boxsize2;
  /// level in tree
  int level;
  unsigned long number;

  // std::vector<OBranchNB> children;
  std::unique_ptr<int[]> children;

  OBranchNB *child;
  /// father of branch
  OBranchNB *prev;
  /// Either child2 of father is branch is child1 and child2 exists or the brother of the father.
  /// Used for iterative tree walk.
  OBranchNB *brother;

  /* projected quantities */
  /// quadropole moment of branch
  /// center of mass
  // PosType center[3];
  // PosType quad[3];
  // PosType mass;
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
  OTreeNB(PType *xp, IndexType *particles, IndexType nparticles, PosType boundary_p1[], PosType boundary_p2[]);
  ~OTreeNB();

  unsigned long Nbranches() { return Nbranches; };

  /**
   \brief A iterator class that allows for movement through the tree without changing
       anything in the tree itself.
 */
  template <typename PType>
  class iterator
  {

  private:
    OBranchNB *current;

  public:
    /// Sets the top or root to the top of "tree".
    iterator(OTreeNB<PType> *tree) { current = tree->top; }
    /// Sets the root to the input branch so that this will be a subtree in branch is not the real root.
    iterator(OBranchNB *branch) { current = branch; }

    iterator operator=(OBranchNB *branch)
    {
      current = branch;
      return *this;
    }

    /// Returns a pointer to the current Branch.
    OBranchNB *operator*() { return current; }

    bool up()
    {
      if (current->prev == NULL)
      {
        current = NULL;
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
      if (current->prev == NULL)
        return true;
      return false;
    }

    bool atLeaf() const
    {
      for (int i = 0; i < 8; ++i)
      {
        if (current->children[i] != NULL)
          return false;
      }
      return true;
    }

    // returns true if not at the end of the tree
    bool TreeWalkStep(bool allowDescent)
    {
      if (allowDescent)
      {
        if (current->children[i] != NULL)
        {
          down();
          return true;
        }
      }
      if (current->brother != NULL)
      {
        current = current->brother;
        return true;
      }
      return false;
    }
  };

  size_t getNbranches() { return nbranches; }

  iterator<PType> begin()
  {
    return iterator<PType>(&top);
  }

  void span8(OBranchNB *current);
  void span4(OBranchNB *current, int dim); // dim is the dimension that is not divided
  void span2(OBranchNB *current, int dim); // dim is the dimension that IS divided

private:
  OBranchNB top;
  /// Array of particle positions
  PType *xxp;

  /// number of branches in tree
  unsigned long Nbranches;
};

// creates 8 children to current with links to their parent and brother
// and sets their center and boxsize
template <typename PType>
void OTreeNB<PType>::span8(OBranchNB *current, IndexType *particles, IndexType nparticles)
{
  if (current == NULL)
  {
    ERROR_MESSAGE();
    std::cerr << "OTreeNB Error: calling spon() on empty tree" << std::endl;
    exit(1);
  }

  current->children.reset(new OBranchNB[8]);
  current->child = children.get();

  for (auto &child : current->children)
  {
    child->prev = current;
    child->level = current->level + 1;
    child->boxsize = current->boxsize / 2.0;
    child->boxsize2 = child->boxsize * child->boxsize;
  }
  for (int i = 0; i < 7; ++i)
  {
    current->children[i]->brother = current->children[i + 1];
  }
  current->children[7]->brother = current->brother;

  PosType half_box = current->boxsize / 2.0;

  PosType center[3];
  center[0] = (current->boundary_p1[0] + current->boundary_p2[0]) / 2.0;
  center[1] = (current->boundary_p1[1] + current->boundary_p2[1]) / 2.0;
  center[2] = (current->boundary_p1[2] + current->boundary_p2[2]) / 2.0;

  size_t m = 0;
  for (int i = -1; i < 2; i += 2)
  {
    for (int j = -1; j < 2; j += 2)
    {
      for (int k = -1; k < 2; k += 2)
      {

        current->children[m].boundary_p1[0] = ((1 - i) / 2) * current->boundary_p1[0] + ((1 + i) / 2) * center[0];

        current->children[m].boundary_p1[1] = ((1 - j) / 2) * current->boundary_p1[1] + ((1 + j) / 2) * center[1];

        current->children[m].boundary_p1[2] = ((1 - k) / 2) * current->boundary_p1[2] + ((1 + k) / 2) * center[2];

        current->children[m].boundary_p2[0] = ((1 + i) / 2) * current->boundary_p2[0] + ((1 - i) / 2) * center[0];

        current->children[m].boundary_p2[1] = ((1 + j) / 2) * current->boundary_p2[1] + ((1 - j) / 2) * center[1];

        current->children[m].boundary_p2[2] = ((1 + k) / 2) * current->boundary_p2[2] + ((1 - k) / 2) * center[2];
        ++m;
      }
    }
  }

  // sort particles into children
  IndexType *p = current->particles;
  IndexType *p_end = p + current->nparticles;
  IndexType np = 0;
  for (int m = 0; m < 7; ++m)
  {
    children[m].particles = p;
    children[m].nparticles = 0;
    if (m == 3)
    {
      children[m].nparticles = current->nparticles - np;
    }
    else
    {
      IndexType *pp = p;
      while (pp != p_end)
      {
        Ptype &x = xxp[*pp];
        if (x[0] >= current->children[m].boundary_p1[0] &&
            x[0] < current->children[m].boundary_p2[0] &&
            x[1] >= current->children[m].boundary_p1[1] &&
            x[1] < current->children[m].boundary_p2[1] &&
            x[2] >= current->children[m].boundary_p1[2] &&
            x[2] < current->children[m].boundary_p2[2])
        {
          std::swap(*p, *pp);
          children[m].nparticles++;
          ++p;
        }
        ++pp;
      }
      np += children[m].nparticles;
    }
  }
  Nbranches += 8;
}

template <typename PType>
void OTreeNB<PType>::span4(OBranchNB *current, int dim, )
{
  if (current == NULL)
  {
    ERROR_MESSAGE();
    std::cerr << "OTreeNB Error: calling spon() on empty tree" << std::endl;
    exit(1);
  }

  current->children.reset(new OBranchNB[4]);
  current->child = children.get();

  for (auto &child : current->children)
  {
    child->prev = current;
    child->level = current->level + 1;
    child->boxsize = current->boxsize / 2.0;
    child->boxsize2 = child->boxsize * child->boxsize;
    child->cubic = current->cubic;
  }
  for (int i = 0; i < 3; ++i)
  {
    current->children[i]->brother = current->children[i + 1];
  }
  current->children[3]->brother = current->brother;

  PosType center[3];
  center[0] = (current->boundary_p1[0] + current->boundary_p2[0]) / 2.0;
  center[1] = (current->boundary_p1[1] + current->boundary_p2[1]) / 2.0;
  center[2] = (current->boundary_p1[2] + current->boundary_p2[2]) / 2.0;

  size_t m = 0;
  for (int i = -1; i < 2; i += 2)
  {
    for (int j = -1; j < 2; j += 2)
    {
      if (dim == 0)
      {
        current->children[m].boundary_p1[0] = current->boundary_p1[0];

        current->children[m].boundary_p1[1] = ((1 - i) / 2) * current->boundary_p1[1] + ((1 + i) / 2) * center[1];

        current->children[m].boundary_p1[2] = ((1 - j) / 2) * current->boundary_p1[2] + ((1 + j) / 2) * center[2];

        current->children[m].boundary_p2[0] = current->boundary_p2[0];

        current->children[m].boundary_p2[1] = ((1 + i) / 2) * current->boundary_p2[1] + ((1 - i) / 2) * center[1];

        current->children[m].boundary_p2[2] = ((1 + j) / 2) * current->boundary_p2[2] + ((1 - j) / 2) * center[2];
      }
      else if (dim == 1)
      {
        current->children[m].boundary_p1[0] = ((1 - i) / 2) * current->boundary_p1[0] + ((1 + i) / 2) * center[0];

        current->children[m].boundary_p1[1] = current->boundary_p1[1];

        current->children[m].boundary_p1[2] = ((1 - j) / 2) * current->boundary_p1[2] + ((1 + j) / 2) * center[2];

        current->children[m].boundary_p2[0] = ((1 + i) / 2) * current->boundary_p2[0] + ((1 - i) / 2) * center[0];

        current->children[m].boundary_p2[1] = current->boundary_p2[1];

        current->children[m].boundary_p2[2] = ((1 + j) / 2) * current->boundary_p2[2] + ((1 - j) / 2) * center[2];
      }
      else
      {
        current->children[m].boundary_p1[0] = ((1 - i) / 2) * current->boundary_p1[0] + ((1 + i) / 2) * center[0];

        current->children[m].boundary_p1[1] = ((1 - j) / 2) * current->boundary_p1[1] + ((1 + j) / 2) * center[1];

        current->children[m].boundary_p1[2] = current->boundary_p1[2];

        current->children[m].boundary_p2[0] = ((1 + i) / 2) * current->boundary_p2[0] + ((1 - i) / 2) * center[0];

        current->children[m].boundary_p2[1] = ((1 + j) / 2) * current->boundary_p2[1] + ((1 - j) / 2) * center[1];

        current->children[m].boundary_p2[2] = current->boundary_p2[2];
      }
      ++m;
    }
  }

  // sort particles into children
  IndexType *p = current->particles;
  IndexType *p_end = p + current->nparticles;
  IndexType np = 0;
  for (int m = 0; m < 4; ++m)
  {
    children[m].particles = p;
    children[m].nparticles = 0;
    if (m == 3)
    {
      children[m].nparticles = current->nparticles - np;
    }
    else
    {
      IndexType *pp = p;
      while (pp != p_end)
      {
        Ptype &x = xxp[*pp];
        if (x[0] >= current->children[m].boundary_p1[0] &&
            x[0] < current->children[m].boundary_p2[0] &&
            x[1] >= current->children[m].boundary_p1[1] &&
            x[1] < current->children[m].boundary_p2[1] &&
            x[2] >= current->children[m].boundary_p1[2] &&
            x[2] < current->children[m].boundary_p2[2])
        {
          std::swap(*p, *pp);
          children[m].nparticles++;
          ++p;
        }
        ++pp;
      }
      np += children[m].nparticles;
    }
  }
  Nbranches += 4;
}

template <typename PType>
void OTreeNB<PType>::span2(OBranchNB *current, int dim)
{
  if (current == NULL)
  {
    ERROR_MESSAGE();
    std::cerr << "OTreeNB Error: calling spon() on empty tree" << std::endl;
    exit(1);
  }

  current->children.reset(new OBranchNB[2]);
  current->child = children.get();

  for (auto &child : current->children)
  {
    child->prev = current;
    child->level = current->level + 1;
  }

  current->children[0]->brother = current->children[1];
  current->children[1]->brother = current->brother;

  for (int m = 0; m < 2; ++m)
  {
    current->children[m].boundary_p1[0] = current->boundary_p1[0];
    current->children[m].boundary_p1[1] = current->boundary_p1[1];
    current->children[m].boundary_p1[2] = current->boundary_p1[2];

    current->children[m].boundary_p2[0] = current->boundary_p2[0];
    current->children[m].boundary_p2[1] = current->boundary_p2[1];
    current->children[m].boundary_p2[2] = current->boundary_p2[2];
  }

  PosType center = (current->boundary_p1[dim] + current->boundary_p2[dim]) / 2.0;
  current->children[0].boundary_p2[dim] = center;
  current->children[1].boundary_p1[dim] = center;

  // sort particles into children
  IndexType *p = current->particles;
  IndexType *p_end = p + current->nparticles;
  for (int m = 0; m < 2; ++m)
  {
    children[m].particles = p;
    children[m].nparticles = 0;
    if (m == 1)
    {
      children[m].nparticles = current->nparticles - current->children[0].nparticles;
    }
    else
    {
      IndexType *pp = p;
      while (pp != p_end)
      {
        Ptype &x = xxp[*pp];
        if (x[0] >= current->children[m].boundary_p1[0] &&
            x[0] < current->children[m].boundary_p2[0] &&
            x[1] >= current->children[m].boundary_p1[1] &&
            x[1] < current->children[m].boundary_p2[1] &&
            x[2] >= current->children[m].boundary_p1[2] &&
            x[2] < current->children[m].boundary_p2[2])
        {
          std::swap(*p, *pp);
          children[m].nparticles++;
          ++p;
        }
        ++pp;
      }
    }
  }
  Nbranches += 2;
}

#endif /* OTreeNB_h */
