/**
 * simpleTree.h
 *
 *  Created on: Oct 14, 2011
 *      Author: bmetcalf
 */

#ifndef QUAD_TREE_H_
#define QUAD_TREE_H_

#include "simpleTree.h"

//short const treeNBdim = 2;

/** \brief Box representing a branch in a tree.  It has four children.  Used in QTreeNB which is used in TreeQuad.
 */
struct QBranchNB{
  QBranchNB();
  ~QBranchNB();

	/// array of particles in QBranchNB
  IndexType *particles;
  IndexType *big_particles;
  IndexType nparticles;
  /// the number of particles that aren't in children
  IndexType Nbig_particles;
  /// Size of largest particle in branch
  PosType maxrsph;
  /// center of mass
  PosType center[2];
  PosType mass;
  /// level in tree
  int level;
  unsigned long number;
  /// bottom, left, back corner of box
  PosType boundary_p1[2];
  /// top, right, front corner of box
  PosType boundary_p2[2];
  QBranchNB *child0;
  QBranchNB *child1;
  QBranchNB *child2;
  QBranchNB *child3;

  /// father of branch
  QBranchNB *prev;
  /// Either child2 of father is branch is child1 and child2 exists or the brother of the father.
  /// Used for iterative tree walk.
  QBranchNB *brother;

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
 * QTreeNB: Tree structure used for force calculation with particles (i.e. stars, Nbody or halos).
 *
 * The tree also contains pointers to the list of positions, sizes and masses of the particles.
 * Also flags for the number of dimensions the tree is defined in (2 or 3), and if multiple
 * masses and sizes should be used.
 */
struct QTreeNB{
	QTreeNB(PosType **xp,IndexType *particles,IndexType nparticles
			 ,PosType boundary_p1[],PosType boundary_p2[]);
	~QTreeNB();

	const bool isEmpty();
	const bool atTop();
	const bool noChild();
	const bool offEnd();
	const unsigned long getNbranches();
	const bool atLeaf(){
		return (current->child0 == NULL)*(current->child1 == NULL)
    *(current->child2 == NULL)*(current->child3 == NULL);
	}
	const bool atLeaf(QBranchNB *branch){
		return (branch->child0 == NULL)*(branch->child1 == NULL)
    *(branch->child2 == NULL)*(branch->child3 == NULL);
	}
  
  
	void getCurrent(IndexType *particles,IndexType *nparticles);
	void moveTop();
	void moveUp();
	void moveToChild(int child);
	bool WalkStep(bool allowDescent);
  
  short empty();
	void attachChildrenToCurrent(QBranchNB *branch0,QBranchNB *branch1
			,QBranchNB *branch2,QBranchNB *branch3);

	QBranchNB *top;
	QBranchNB *current;
	/// Array of particle positions
	PosType **xp;

  /**
   *  \brief A iterator class fore TreeStruct that allows for movement through the tree without changing
   *      anything in the tree itself.
   *
   *   This class should be able to preform all of the constant movements within the tree without causing
   *   any change to the tree.
   */
  class iterator{
    
  private:
    QBranchNB *current;
    QBranchNB *top;
    
  public:
    /// Sets the top or root to the top of "tree".
    iterator(QTreeNB * tree){current = top = tree->top;}
    /// Sets the root to the input branch so that this will be a subtree in branch is not the real root.
    iterator(QBranchNB *branch){current = top = branch;}
    
    /// Returns a pointer to the current Branch.
    QBranchNB *operator*(){return current;}
    
    void movetop(){current = top;}
    
    /// Same as up()
    bool operator++(){ return up();}
    
    /// Same as up()
    bool operator++(int){ return up();}
    
    bool up();
    /// Move to child
    bool down(short child);
    bool TreeWalkStep(bool allowDescent);
  };

private:
  /// number of branches in tree
	unsigned long Nbranches;
	void _freeQTree(short child);
};

typedef struct QTreeNB * QTreeNBHndl;
typedef int QTreeNBElement;

/**
 * \brief TreeQuad is a class for calculating the deflection, kappa and gamma by tree method.
 *
 * TreeQuad is evolved from TreeSimple and TreeForce.  It splits each cell into four equal area
 * subcells instead of being a binary tree like TreeSimple.  When the "particles" are given sizes
 * the tree is built in such a way the large particles are stored in branches that are no smaller
 * than their size.  In this way particles are stored on all levels of the tree and not just in the
 * leaves.  This improves efficiency when particles of a wide range of sizes overlap in 2D.
 *
  * The default value of theta = 0.1 generally gives better than 1% accuracy on alpha.
 * The shear and kappa is always more accurate than the deflection.
 *
 */

class TreeQuad {
public:
	TreeQuad(
			PosType **xpt
			,float *my_masses
			,float *my_sizes
			,IndexType Npoints
			,bool Multimass
			,bool Multisize
			,PosType my_sigma_background = 0
			,int bucket = 5
			,PosType theta_force = 0.1
			);
	TreeQuad(
			PosType **xpt
			,LensHaloHndl *my_halos
			,IndexType Npoints
			,PosType my_sigma_background = 0
			,int bucket = 5
			,PosType theta_force = 0.1
			);
	virtual ~TreeQuad();

	// virtual void force2D(PosType const *ray,PosType *alpha,KappaType *kappa,KappaType *gamma,bool no_kappa);


    // PHI BY Fabien :
    virtual void force2D(PosType const *ray,PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,bool no_kappa);

    // PHI BY Fabien :
    virtual void force2D_recur(PosType const *ray,PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,bool no_kappa);
    
    virtual void printParticlesInBranch(unsigned long number);

	virtual void printBranchs(int level = -1);

protected:

	PosType **xp;
	bool MultiMass;
	bool MultiRadius;
	float *masses;
	float *sizes;

	IndexType Nparticles;
	PosType sigma_background;
	int Nbucket;

	PosType force_theta;

	QTreeNBHndl tree;
	IndexType *index;

	bool haloON;
	LensHaloHndl *halos;

	PosType realray[2];
	int incell,incell2;

	QTreeNBHndl BuildQTreeNB(PosType **xp,IndexType Nparticles,IndexType *particles);
	void _BuildQTreeNB(IndexType nparticles,IndexType *particles);

	inline short WhichQuad(PosType *x,QBranchNB &branch);

	//inline bool atLeaf();
	inline bool inbox(PosType *ray,PosType *p1,PosType *p2){
	  return (ray[0]>=p1[0])*(ray[0]<=p2[0])*(ray[1]>=p1[1])*(ray[1]<=p2[1]);
	}
	int cutbox(PosType *ray,PosType *p1,PosType *p2,float rmax);
	
	void CalcMoments();
	void rotate_coordinates(PosType **coord);

	// Internal profiles for a Gaussian particle
	virtual inline PosType alpha_h(PosType r2s2,PosType sigma){
	  return (sigma > 0.0 ) ? ( exp(-0.5*r2s2) - 1.0 ) : -1.0;
	}
	virtual inline PosType kappa_h(PosType r2s2,PosType sigma){
	  return 0.5*r2s2*exp(-0.5*r2s2);
	}
	virtual inline PosType gamma_h(PosType r2s2,PosType sigma){
	  return (sigma > 0.0 ) ? (-2.0 + (2.0 + r2s2)*exp(-0.5*r2s2) ) : -2.0;
	}
	virtual inline PosType phi_h(PosType r2s2,PosType sigma){
		ERROR_MESSAGE();  // not yet written
		exit(1);
		return 0;
	}

	QTreeNBHndl rotate_simulation(PosType **xp,IndexType Nparticles,IndexType *particles
			,PosType **coord,PosType theta,float *rsph,float *mass
			,bool MultiRadius,bool MultiMass);
	QTreeNBHndl rotate_project(PosType **xp,IndexType Nparticles,IndexType *particles
			,PosType **coord,PosType theta,float *rsph,float *mass
			,bool MultiRadius,bool MultiMass);
	 void cuttoffscale(QTreeNBHndl tree,PosType *theta);

	 void walkTree_recur(QBranchNB *branch,PosType const *ray,PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,bool no_kappa);

   };

#endif /* QUAD_TREE_H_ */
