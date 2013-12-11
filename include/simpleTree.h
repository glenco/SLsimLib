/**
 * simpleTree.h
 *
 *  Created on: Oct 14, 2011
 *      Author: bmetcalf
 */

#ifndef SIMP_TREE_H_
#define SIMP_TREE_H_

#include "standard.h"
#include "Tree.h"
#include "lens_halos.h"

/** \brief Box representing a branch in a tree.  It has four children.  Used in TreeNBStruct which is used in TreeForce.
 */
struct BranchNB{
  BranchNB(int Ndim);
  ~BranchNB();

	/// array of particles in BranchNB
  IndexType *particles;
  IndexType nparticles;
  /// the number of particles that aren't in children
  IndexType big_particle;
  /// Size of largest particle in branch
  PosType maxrsph;
  /// center of mass
  PosType *center;
  PosType mass;
  /// level in tree
  int level;
  unsigned long number;
  /// bottom, left, back corner of box
  PosType *boundary_p1;
  /// top, right, front corner of box
  PosType *boundary_p2;
  BranchNB *child1;
  BranchNB *child2;
  /// father of branch
  BranchNB *prev;
  /// Either child2 of father is branch is child1 and child2 exists or the brother of the father.
  /// Used for iterative tree walk.
  BranchNB *brother;

  /* projected quantities */
  /// quadropole moment of branch
  PosType quad[3];
  /// largest dimension of box
  PosType rmax;
  /// the critical distance below which a branch is opened in the
  PosType rcrit_angle;
                /* force calculation */
  PosType rcrit_part;
  //PosType cm[2]; /* projected center of mass */

};

/** \brief
 * TreeNB: Tree structure used for force calculation with particles (i.e. stars, Nbody or halos).
 *
 * The tree also contains pointers to the list of positions, sizes and masses of the particles.
 * Also flags for the number of dimensions the tree is defined in (2 or 3), and if multiple
 * masses and sizes should be used.
 */
typedef struct TreeNBStruct{
  BranchNB *top;
  BranchNB *current;
  /// number of branches in tree
  unsigned long Nbranches;
  /// Dimension of tree, 2 or 3.  This will dictate how the force is calculated.
  short Ndimensions;
  /// Array of particle positions
  PosType **xp;
} TreeNBStruct;

typedef struct TreeNBStruct * TreeNBHndl;
typedef int TreeNBElement;

/**
 * \brief A C++ class wrapper for the bianary treeNB used in the Nobody force calculation, but
 * also useful for general purpose searches.
 *
 * Most of the code in TreeNB.c and TreeDriverNB.c is duplicated here as private methods and
 * a few public ones.
 */
class TreeSimple {
public:
	TreeSimple(PosType **xp,IndexType Npoints,int bucket = 5,int dimensions = 2,bool median = true);
	virtual ~TreeSimple();

	/// \brief Finds the points within a circle around center and puts their index numbers in a list
	void PointsWithinCircle(PosType center[2],float radius,std::list<unsigned long> &neighborkist);
	/// \brief Finds the points within an ellipse around center and puts their index numbers in a list
	void PointsWithinEllipse(PosType center[2],float a_max,float a_min,float posangle,std::list<unsigned long> &neighborkist);
	/// \brief Finds the nearest N neighbors and puts their index numbers in an array, also returns the distance to the Nth neighbor for calculating smoothing
	void NearestNeighbors(PosType *ray,int Nneighbors,float *rsph,IndexType *neighbors);

protected:

	int Ndim,incell,incell2;
	TreeNBHndl tree;
	IndexType *index;
	IndexType Nparticles;
	bool median_cut;
	int Nbucket;
	PosType realray[2];
	PosType **xp;


	TreeNBHndl BuildTreeNB(PosType **xp,IndexType Nparticles,IndexType *particles,int Ndimensions,PosType theta);
	void _BuildTreeNB(TreeNBHndl tree,IndexType nparticles,IndexType *particles);

	void _PointsWithin(PosType *ray,float *rmax,std::list<unsigned long> &neighborkist);
	void _NearestNeighbors(PosType *ray,int Nneighbors,unsigned long *neighbors,PosType *rneighbors);

	BranchNB *NewBranchNB(IndexType *particles,IndexType nparticles
			  ,PosType boundary_p1[],PosType boundary_p2[]
			  ,PosType center[],int level,unsigned long branchNBnumber);
	void FreeBranchNB(BranchNB *branchNB);
	TreeNBHndl NewTreeNB(IndexType *particles,IndexType nparticles
			 ,PosType boundary_p1[],PosType boundary_p2[],
			     PosType center[],short Ndimensions);
	void freeTreeNB(TreeNBHndl tree);
	short emptyTreeNB(TreeNBHndl tree);
	void _freeTreeNB(TreeNBHndl tree,short child);
	bool isEmptyNB(TreeNBHndl tree);
	bool atTopNB(TreeNBHndl tree);
	bool noChildNB(TreeNBHndl tree);
	bool offEndNB(TreeNBHndl tree);
	void getCurrentNB(TreeNBHndl tree,IndexType *particles,IndexType *nparticles);
	unsigned long getNbranchesNB(TreeNBHndl tree);
	void moveTopNB(TreeNBHndl tree);
	void moveUpNB(TreeNBHndl tree);
	void moveToChildNB(TreeNBHndl tree,int child);
	void insertChildToCurrentNB(TreeNBHndl tree, IndexType *particles,IndexType nparticles
				  ,PosType boundary_p1[],PosType boundary_p2[]
				  ,PosType center[],int child);
	void attachChildToCurrentNB(TreeNBHndl tree,BranchNB &data,int child);
	bool TreeNBWalkStep(TreeNBHndl tree,bool allowDescent);

	inline bool atLeaf(){
		return (tree->current->child1 == NULL)*(tree->current->child2 == NULL);
	}
	inline bool inbox(const PosType* center,PosType *p1,PosType *p2){
	  return (center[0]>=p1[0])*(center[0]<=p2[0])*(center[1]>=p1[1])*(center[1]<=p2[1]);
	}
	
	/*TreeNBHndl rotate_simulation(PosType **xp,IndexType Nparticles,IndexType *particles
			,PosType **coord,PosType theta,float *rsph,float *mass
			,bool MultiRadius,bool MultiMass);
	TreeNBHndl rotate_project(PosType **xp,IndexType Nparticles,IndexType *particles
			,PosType **coord,PosType theta,float *rsph,float *mass
			,bool MultiRadius,bool MultiMass);
	TreeNBHndl spread_particles(PosType **xp,IndexType Nparticles,IndexType *particles
			,PosType theta,float *rsph,float *mass,bool MultiRadius,bool MultiMass);
	 void cuttoffscale(TreeNBHndl tree,PosType *theta);*/
};

//PosType **PosTypeMatrix(long nrl, long nrh, long ncl, long nch);
//void free_PosTypeMatrix(PosType **m, long nrl, long nrh, long ncl, long nch);

#endif /* SIMP_TREE_H_ */
