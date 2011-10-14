 /*
 * Code Name:     tree.h                                       
 * Programmer:    Ben Metcalf
 * Discription:                                                     
 * Comments:                           
 */

#include <Tree.h>
#include <cosmo.h>

#ifndef pi
#define pi  3.141593
#endif

#ifndef treeNBdim
#define treeNBdim 2  // dimension of boxes in tree
#endif

/** type for particle positions and boundaries etc **/

#ifndef PosType_declare
#define PosType_declare
typedef double PosType;
#endif

#ifndef IndexType_declare
#define IndexType_declare
typedef unsigned long IndexType;
#endif

/***** Exported Types *****/

#ifndef bool_declare
#define bool_declare
typedef enum {false, true} bool;
#endif

#ifndef treeNBtypes_declare
#define treeNBtypes_declare
/** \brief
 * Branch of treeNB for tree in force calculation.
 *
 * Contains list of particle indices, boundaries of box, moments of mass distribution
 * , pointers to other branches, etc.
 */
typedef struct branchNBstruct{
	/// array of particles in BranchNB
  IndexType *particles;
  IndexType nparticles;
  /// the number of particles that aren't in children
  IndexType big_particle;
  /// Size of largest particle in branch
  PosType maxrsph;
  /// center of mass
  PosType center[treeNBdim];
  PosType mass;
  /// level in tree
  int level;
  unsigned long number;
  /// bottom, left, back corner of box
  PosType boundary_p1[treeNBdim];
  /// top, right, front corner of box
  PosType boundary_p2[treeNBdim];
  struct branchNBstruct *child1;
  struct branchNBstruct *child2;
  /// father of branch
  struct branchNBstruct *prev;
  /// Either child2 of father is branch is child1 and child2 exists or the brother of the father.
  /// Used for iterative tree walk.
  struct branchNBstruct *brother;

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

} BranchNB;

/** \brief
 * Obsolete - Version of branchNBstruct used to save the tree to disk.
 */
typedef struct branchNBstruct2{
	 // list of particles in BranchNB
  IndexType particles;
  IndexType nparticles;
  // the number of particles that aren't in children
  IndexType big_particle;
  PosType maxrsph;
  PosType center[treeNBdim];
  int level;
  PosType boundary_p1[treeNBdim];
  PosType boundary_p2[treeNBdim];
  int child1;
  int child2;
  int prev;
} RelativeBranchNB;

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
  /// true if particles have different masses.
  bool MultiMass;
  /// true if particles have different sizes.
  bool MultiRadius;

  /// Array of particle positions
  PosType **xp;
  /// Array of particle sizes
  float *rsph;
  /// Array of particle masses
  float *mass;
} TreeNBStruct;

/* typedef struct Point{ */
/*   Point *up; */
/*   Point *down; */
/*   unsigned long id; */
/* } Point; */

typedef struct TreeNBStruct * TreeNBHndl;
typedef int TreeNBElement;

/** \brief Structure for a lens consisting of simulation particles */
typedef struct simlens{
  double zlens;
  double zsource;

  int Nspecies;  /** number of particle species */
  double phi;    /** orientation of the simulation projection */
  double theta;    /** orientation of the simulation projection */
  double mass_units; /** mass units for particles */
  double theta_force;

  double r_source;

  int Nsph; /** mass units for particles */

  double interpolation_scale; // grid scale below which interpolation is used

  char treefilenames[50];
  char simfilename[50];

  /* private data */
  TreeNBHndl tree;  /** pointers to the tree structures of each species */
  //int *firstparticle;  // the first particle of each species
  double *masses;     /** masses of particle species */
  double Sigma_crit;  /** critical surface density */
  double **coord;    /** orientation for internal use */

  IndexType Nparticles;
  IndexType *particles;
  float *rsph;    /** SPH smoothing scale for each particle */
  PosType **xp;    /** array of particle positions */

  PosType **xp_2d;    /** array of projected particle positions */
  
} SimLens;

#endif

/***** Constructors/Destructors*****/

/** Creates a new TreeNB structure */
TreeNBHndl NewTreeNB(IndexType *particles,IndexType nparticles
		 ,PosType boundary_p1[],PosType boundary_p2[],
		     PosType center[],short Ndimensions);
void freeTreeNB(TreeNBHndl tree);
short emptyTreeNB(TreeNBHndl tree);
void _freeTreeNB(TreeNBHndl tree,short child);

BranchNB *NewBranchNB(IndexType *particles,IndexType nparticles
		  ,PosType boundary_p1[],PosType boundary_p2[]
		  ,PosType center[2],int level,unsigned long branchNBnumber);
void FreeBranchNB(BranchNB *branchNB);

/***** Access functions *****/

bool isEmptyNB(TreeNBHndl tree);
bool atTopNB(TreeNBHndl tree);
bool offEndNB(TreeNBHndl tree);
bool noChildNB(TreeNBHndl tree);

/*unsigned long *getCurrent(TreeNBHndl tree,IndexType *nparticles);*/
void getCurrentNB(TreeNBHndl tree,IndexType *particles,IndexType *nparticles);

unsigned long getNbranchNBes(TreeNBHndl tree);


/***** Manipulation procedures *****/

void moveTopNB(TreeNBHndl tree);

void moveUpNB(TreeNBHndl tree);

void moveToChildNB(TreeNBHndl tree,int child);

void insertChildToCurrentNB(TreeNBHndl tree, IndexType *particles,IndexType nparticles
			  ,PosType boundary_p1[],PosType boundary_p2[]
			  ,PosType center[],int child);

void attachChildToCurrentNB(TreeNBHndl tree,BranchNB data,int child);

bool TreeNBWalkStep(TreeNBHndl tree,bool allowDescent);

/***** Other operations *****/

void printTreeNB(TreeNBHndl tree,PosType **xp);

void printBranchNB(BranchNB *branchNB,PosType **xp,short Ndim);

void saveTreeNB(TreeNBHndl tree,IndexType *particles,float *rsph,char *filename);
void _saveTreeNB(TreeNBHndl tree,RelativeBranchNB *tree_arr,IndexType *particles);
TreeNBHndl readTreeNB(IndexType *particles,float *rsph,IndexType Nparticles,char *filename);
void _readTreeNB(TreeNBHndl tree,RelativeBranchNB *tree_arr,IndexType *particles
		,unsigned long current);
void saveSPHsmoothing(TreeNBHndl tree,IndexType *particles,float *rsph,char *filename);
void readSmoothingNB(float *rsph,char *filename);

/** routines in TreeNBDriver.c **/

TreeNBHndl BuildTreeNB(PosType **xp,float *rsph,float *mass,bool MultiRadius
		,bool MultiMass,IndexType Nparticles,IndexType *particles,int Ndimensions
		,double theta);
IndexType *NearestNeighborNB(TreeNBHndl tree,double *ray,int Nneighbors
				 ,float *rsph);
int inboxNB(double ray[],PosType *p1,PosType *p2);
int cutboxNB(double ray[],PosType *p1,PosType *p2,PosType rmax);
TreeNBHndl rotate_simulation(PosType **xp,IndexType Nparticles,IndexType *particles
		,double **coord,double theta,float *rsph,float *mass
		,bool MultiRadius,bool MultiMass);
TreeNBHndl rotate_project(PosType **xp,IndexType Nparticles,IndexType *particles
		,double **coord,double theta,float *rsph,float *mass
		,bool MultiRadius,bool MultiMass);
TreeNBHndl spread_particles(PosType **xp,IndexType Nparticles,IndexType *particles
		,double theta,float *rsph,float *mass,bool MultiRadius,bool MultiMass);
 void cuttoffscale(TreeNBHndl tree,double *theta);

/** routines in TreeNBForce **/

float *FindRSPH(TreeNBHndl tree,int Nsph);
void TreeNBForce2D(TreeNBHndl tree,double *ray
		 ,double *alpha,double *kappa,double *gamma,bool no_kappa);
void _TreeNBForce2D(TreeNBHndl tree,double *ray
		  ,double *alpha,double *kappa,double *gamma,bool no_kappa);
void _TreeNBParticleForce2D(TreeNBHndl tree,double *ray
		  ,double *alpha,double *kappa,double *gamma,bool no_kappa);
void TreeNBParticleForce2Diter(TreeNBHndl tree,double *ray
		  ,double *alpha,double *kappa,double *gamma,bool no_kappa
		  ,double (*alpha_internal)(double r,float rmax)
		  ,double (*kappa_internal)(double r,float rmax)
		  ,double (*gamma_internal)(double r,float rmax));

double alpha_o(double r,float sigma);
double kappa_o(double r,float sigma);
double gamma_o(double r,float sigma);
int OpenBox(TreeNBHndl tree,PosType r);

/** routines in readfiles **/

SimLens *readparams(char *filename,struct cosmology *cosmo);
void PrintSimLens(SimLens *lens);
void readpositions(SimLens *lens);
/*void rayshooter(double *ray,unsigned long Nrays,double *alpha,double *gamma
  ,double *kappa,double *invmag,char *paramfile);*/
/*void rayshooter(double *ray,unsigned long Nrays,double *alpha,double *gamma
		,double*kappa,double *invmag,char *paramfile);*/
/*void rayshooter(unsigned long Nrays,Point *points,char *paramfile);*/
void rayshooterNB(unsigned long Nrays,Point *points,TreeHndl i_tree,char *paramfile,bool no_kappa);
PosType **PosTypeMatrix(long nrl, long nrh, long ncl, long nch);
void free_PosTypeMatrix(PosType **m, long nrl, long nrh, long ncl, long nch);
