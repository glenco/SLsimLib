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

	short empty();
	bool isEmpty();
	bool atTop();
	bool noChild();
	bool offEnd();
	void getCurrent(IndexType *particles,IndexType *nparticles);
	unsigned long getNbranches();
	void moveTop();
	void moveUp();
	void moveToChild(int child);
	void attachChildrenToCurrent(QBranchNB *branch0,QBranchNB *branch1
			,QBranchNB *branch2,QBranchNB *branch3);
	bool WalkStep(bool allowDescent);
	inline bool atLeaf(){
		return (current->child0 == NULL)*(current->child1 == NULL)
				*(current->child2 == NULL)*(current->child3 == NULL);
	}
	inline bool atLeaf(QBranchNB *branch){
		return (branch->child0 == NULL)*(branch->child1 == NULL)
				*(branch->child2 == NULL)*(branch->child3 == NULL);
	}
	QBranchNB *top;
	QBranchNB *current;
	/// number of branches in tree
	unsigned long Nbranches;
	/// Array of particle positions
	PosType **xp;

private:
	void _freeQTree(short child);
};

typedef struct QTreeNB * QTreeNBHndl;
typedef int QTreeNBElement;

/**
 * \brief QuadTree is a class for calculating the deflection, kappa and gamma by tree method.
 *
 * QuadTree is evolved from SimpleTree and ForceTree.  It splits each cell into four equal area
 * subcells instead of being a binary tree like SimpleTree.  When the "particles" are given sizes
 * the tree is built in such a way the large particles are stored in branches that are no smaller
 * than their size.  In this way particles are stored on all levels of the tree and not just in the
 * leaves.  This improves efficiency when particles of a wide range of sizes overlap in 2D.
 *
  * The default value of theta = 0.1 generally gives better than 1% accuracy on alpha.
 * The shear and kappa is always more accurate than the deflection.
 *
 */

class QuadTree {
public:
	QuadTree(
			PosType **xpt
			,float *my_masses
			,float *my_sizes
			,IndexType Npoints
			,bool Multimass
			,bool Multisize
			,double my_sigma_background = 0
			,int bucket = 5
			,double theta_force = 0.1
			);
	virtual ~QuadTree();

	virtual void force2D(double *ray,double *alpha,KappaType *kappa,KappaType *gamma,bool no_kappa);
	virtual void force2D_recur(double *ray,double *alpha,KappaType *kappa,KappaType *gamma,bool no_kappa);
	virtual void printParticlesInBranch(unsigned long number);
	virtual void printBranchs(int level = -1);
	void set_force_theta(double ft){force_theta=ft;};

protected:

	QuadTree(
			PosType **xpt
			,HaloStructure *my_halo_params
			,IndexType Npoints
			,double my_sigma_background = 0
			,int bucket = 5
			,double theta_force = 0.1
			,bool NSIE_ON = false
			);

	PosType **xp;
	bool MultiMass;
	bool MultiRadius;
	float *masses;
	float *sizes;

	IndexType Nparticles;
	double sigma_background;
	int Nbucket;

	double force_theta;

	QTreeNBHndl tree;
	IndexType *index;

	bool haloON;
	bool NSIE_ON;
	HaloStructure *halo_params;

	PosType realray[2];
	int incell,incell2;

	QTreeNBHndl BuildQTreeNB(PosType **xp,IndexType Nparticles,IndexType *particles);
	void _BuildQTreeNB(IndexType nparticles,IndexType *particles);

	inline short WhichQuad(double *x,QBranchNB &branch);

	//inline bool atLeaf();
	inline bool inbox(PosType *ray,PosType *p1,PosType *p2){
	  return (ray[0]>=p1[0])*(ray[0]<=p2[0])*(ray[1]>=p1[1])*(ray[1]<=p2[1]);
	}
	int cutbox(PosType *ray,PosType *p1,PosType *p2,float rmax);
	
	void CalcMoments();
	void rotate_coordinates(double **coord);

	virtual void force_halo(double *alpha,KappaType *kappa,KappaType *gamma,double *xcm,HaloStructure& halo_params,bool no_kappa);

	// Internal profiles for a Gaussian particle
	virtual inline double alpha_h(double r2s2,double sigma){
	  return (sigma > 0.0 ) ? ( exp(-0.5*r2s2) - 1.0 ) : -1.0;
	}
	virtual inline double kappa_h(double r2s2,double sigma){
	  return 0.5*r2s2*exp(-0.5*r2s2);
	}
	virtual inline double gamma_h(double r2s2,double sigma){
	  return (sigma > 0.0 ) ? (-2.0 + (2.0 + r2s2)*exp(-0.5*r2s2) ) : -2.0;
	}
	virtual inline double phi_h(double r2s2,double sigma){
		ERROR_MESSAGE();  // not yet written
		exit(1);
		return 0;
	}

	QTreeNBHndl rotate_simulation(PosType **xp,IndexType Nparticles,IndexType *particles
			,double **coord,double theta,float *rsph,float *mass
			,bool MultiRadius,bool MultiMass);
	QTreeNBHndl rotate_project(PosType **xp,IndexType Nparticles,IndexType *particles
			,double **coord,double theta,float *rsph,float *mass
			,bool MultiRadius,bool MultiMass);
	 void cuttoffscale(QTreeNBHndl tree,double *theta);

	 void walkTree_recur(QBranchNB *branch,double *ray,double *alpha,KappaType *kappa,KappaType *gamma,bool no_kappa);
};

/** \ingroup DeflectionL2
 *
 * \brief A class for calculating the deflection, kappa and gamma caused by a collection of halos
 * with truncated power-law mass profiles.
 *
 * Derived from the QuadTree class.  The "particles" are replaced with spherical halos.
 *The truncation is in 2d not 3d. \f$ \Sigma \propto r^\beta \f$ so beta would usually be negative.
 *
 *
 * The default value of theta = 0.1 generally gives better than 1% accuracy on alpha.
 * The shear and kappa is always more accurate than the deflection.
 */
class QuadTreePowerLaw : public QuadTree{

public:
	QuadTreePowerLaw(float beta,PosType **xp,IndexType Npoints,HaloStructure *par_internals
			,double my_sigma_bk=0.0,int bucket = 5,PosType theta = 0.1);

	~QuadTreePowerLaw();

protected:
	float beta; // logorithmic slop of 2d mass profile

	// Override internal structure of halos
	inline double alpha_h(double x,double xmax){
	  if(x==0) x=1e-5;
		return -1.0*pow(x/xmax,beta+2);
	}
	inline double kappa_h(double x,double xmax){
	  if(x==0) x=1e-5;
		return 0.5*(beta+2)*pow(x/xmax,beta)*x*x/(xmax*xmax);
	}
	inline double gamma_h(double x,double xmax){
	  if(x==0) x=1e-5;
		return 0.5*beta*pow(x/xmax,beta+2);
	}
	inline double phi_h(double x,double xmax){
		ERROR_MESSAGE();
		std::cout << "time delay has not been fixed for PowerLaw profile yet." << std::endl;
		exit(1);
		return 0.0;
	}
};

/** \ingroup DeflectionL2
 *
 * \brief A class for calculating the deflection, kappa and gamma caused by a collection of NFW
 * halos.
 *
 * Derived from the QuadTree class.  The "particles" are replaced with spherical NFW halos.
 *
 * This class uses the true expressions for the NFW profile.  This is
 * time consuming and not usually necessary. See QuadTreePseudoNFW for a faster alternative.
 *
* The default value of theta = 0.1 generally gives better than 1% accuracy on alpha.
 * The shear and kappa is always more accurate than the deflection.
 *
 */
class QuadTreeNFW : public QuadTree{

public:
	QuadTreeNFW(PosType **xp,IndexType Npoints,HaloStructure *par_internals
			,double my_sigma_bk = 0.0,int bucket = 5,PosType theta = 0.1);
	~QuadTreeNFW();

protected:

	//static double *ft, *gt, *g2t;
	static double *ftable,*gtable,*g2table,*xtable;
	static long ob_count;


	// Override internal structure of halos
	inline double alpha_h(double x,double xmax){
		return -1.0*InterpolateFromTable(gtable,xtable,x)/InterpolateFromTable(gtable,xtable,xmax);
	}
	inline double kappa_h(double x,double xmax){
		return 0.5*x*x*InterpolateFromTable(ftable,xtable,x)/InterpolateFromTable(gtable,xtable,xmax);
	}
	inline double gamma_h(double x,double xmax){
		return -0.25*x*x*InterpolateFromTable(g2table,xtable,x)/InterpolateFromTable(gtable,xtable,xmax);
	}
	inline double phi_h(double x,double xmax){
		ERROR_MESSAGE();
		std::cout << "time delay has not been fixed for NFW profile yet." << std::endl;
		exit(1);
		return 0.0;
	}
	void make_tables();

	double gfunctionRmax(double rm,double x);
	double ffunctionRmax(double rm,double x);
	double g2functionRmax(double rm,double x);
};

/** \ingroup DeflectionL2
 *
 * \brief A class for calculating the deflection, kappa and gamma caused by a collection of
 * halos with a double power-law mass profile.
 *
 * Derived from the QuadTree class.  The "particles" are replaced with spherical halos
 * with \f$ \Sigma \propto 1/(1 + r/r_s )^\beta \f$ so beta would usually be positive.
 *
 * An NFW profile is approximated beta = 2 .
 *
 * The default value of theta = 0.1 generally gives better than 1% accuracy on alpha.
 * The shear and kappa is always more accurate than the deflection.
 */
class QuadTreePseudoNFW : public QuadTree{

public:
	QuadTreePseudoNFW(double beta,PosType **xp,IndexType Npoints,HaloStructure *par_internals
			,double my_sigma_bk = 0.0,int bucket = 5,PosType theta = 0.1);
	~QuadTreePseudoNFW();

protected:

	double beta;
	static double *mhattable,*xtable;
	static long ob_count;


	// Override internal structure of halos
	inline double alpha_h(double x,double xmax){
		return -1.0*InterpolateFromTable(mhattable,xtable,x)/InterpolateFromTable(mhattable,xtable,xmax);
	}
	inline double kappa_h(double x,double xmax){
		return 0.5*x*x/InterpolateFromTable(mhattable,xtable,xmax)/pow(1+x,beta);
	}
	inline double gamma_h(double x,double xmax){
		return (0.5*x*x/pow(1+x,beta) - InterpolateFromTable(mhattable,xtable,x))/InterpolateFromTable(mhattable,xtable,xmax);
	}
	inline double phi_h(double x,double xmax){
		ERROR_MESSAGE();
		std::cout << "time delay has not been fixed for PseudoNFW profile yet." << std::endl;
		exit(1);
		return 0.0;
	}
	void make_tables();
};
/** \ingroup DeflectionL2
 *
 * \brief A class for calculating the deflection, kappa and gamma caused by a collection of
 * halos that are each non-singular isothermal ellipsoids.
 *
 */
class QuadTreeNSIE : public QuadTree{
public:
	  QuadTreeNSIE(PosType **xp,IndexType Npoints,HaloStructure *par_internals
				,double my_sigma_bk = 0.0,int bucket = 5,PosType theta = 0.1);
	  ~QuadTreeNSIE();

	  void test_force_halo(HaloStructure &halo_params);
protected:

	  virtual void force_halo(double *alpha,KappaType *kappa,KappaType *gamma,double *xcm
	  		,HaloStructure& halo_params,bool no_kappa);
private:
	  /// not used
	  inline double alpha_h(double x,double xmax){assert(0); return 0.0;}
	  /// not used
	  inline double kappa_h(double x,double xmax){assert(0); return 0.0;}
	  /// not used
	  inline double gamma_h(double x,double xmax){assert(0); return 0.0;}
	  /// not used
	  inline double phi_h(double x,double){assert(0); return 0.0;}
};

/** \ingroup DeflectionL2
 *
 * \brief A class for calculating the deflection, kappa and gamma caused by a collection of
 * halos that are each non-singular isothermal ellipsoids.
 *
 */
class QuadTreeNFW_NSIE : public QuadTreeNFW{
public:
	  QuadTreeNFW_NSIE(PosType **xp,IndexType Npoints,HaloStructure *par_internals
				,double my_sigma_bk = 0.0,int bucket = 5,PosType theta = 0.1);
	  ~QuadTreeNFW_NSIE();

	  virtual void force2D(double *ray,double *alpha,KappaType *kappa,KappaType *gamma,bool no_kappa);
	  virtual void force2D_recur(double *ray,double *alpha,KappaType *kappa,KappaType *gamma,bool no_kappa);

protected:

	  QuadTreeNSIE* qtreensie;
};

#endif /* QUAD_TREE_H_ */
