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

#ifndef treeNBtypes_declare
#define treeNBtypes_declare
/** \brief
 * Branch of treeNB for tree in force calculation.
 *
 * Contains list of particle indices, boundaries of box, moments of mass distribution
 * , pointers to other branches, etc.
 */


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


/* typedef struct Point{ */
/*   Point *up; */
/*   Point *down; */
/*   unsigned long id; */
/* } Point; */

/** \brief Structure for a lens consisting of simulation particles */
typedef struct SimLens{
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

/** routines in readfiles **/

SimLens *readparams(char *filename,struct COSMOLOGY *cosmo);
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
