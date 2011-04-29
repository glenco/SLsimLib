 /*
 * Code Name:     tree.h                                       
 * Programmer:    Ben Metcalf
 * Description:
 * Comments:  This version uses a linked list for the points so that
 the tree can be expended dynamically
 */
#include <point.h>
#include <Kist.h>
#include <List.h>

/***** Exported Types *****/

#ifndef error_message
#define ERROR_MESSAGE() printf("ERROR: file: %s line: %i\n",__FILE__,__LINE__)
#endif

#ifndef line_message
#define PRINT_LINE() printf("file: %s line: %i\n",__FILE__,__LINE__)
#endif

#ifndef Boolean_declare
#define Boolean_declare
typedef enum {False, True} Boolean;
#endif

#ifndef criterion_declare
#define criterion_declare
typedef enum{TotalArea,EachImage,Resolution,FillHoles} ExitCriterion;
#endif

#ifndef swap_declare
#define swap_declare
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#endif

/*typedef struct branchstruct{
   struct Point *points;
  unsigned long npoints;
  double center[2];
  int level;
  unsigned long number;
  double boundery_p1[2];
  double boundery_p2[2];
  struct branchstruct *child1;
  struct branchstruct *child2;
  struct branchstruct *prev;
} Branch;*/

/*
typedef struct branchstruct2{
  unsigned long points;        // list of points in Branch
  unsigned long npoints;
  double center[2];
  int level;
  double boundery_p1[2];
  double boundery_p2[2];
  int child1;
  int child2;
  int prev;
} RelativeBranch;
*/

#ifndef treetypes_declare
#define treetypes_declare

// Tree: Exported struct
typedef struct TreeStruct{
  Branch *top;
  Branch *current;
  unsigned long Nbranches;  /* number of barnches in tree */
  PointList *pointlist;
} TreeStruct;

typedef struct TreeStruct *TreeHndl;
typedef int TreeElement;

// this can be used for image or curve
typedef struct ImageInfo{
  Point *points;
  KistHndl imagekist;     // later addition, holds all points in image, will replace points eventually
  unsigned long Npoints;
  double gridrange[3];
  // gridrange[2] minimum grid size in image
  // gridrange[0] maximum grid size in outerborder
  // gridrange[1] maximum grid size in image

  double centroid[2];
  double area;
  double area_error;
  KistHndl innerborder;
  KistHndl outerborder;
  short Nencircled;
} ImageInfo;

typedef struct ImageInfoKist{
	int NimagesMax;
	int Nimages;
	unsigned long Npoints;
	double gridrange[3];

	// seporate kists for each image
	KistHndl *imagekist;     // later addition, holds all points in image, will replace points eventually
	KistHndl *innerborder;
	KistHndl *outerborder;
	double **centroid[2];
	double *area;
	double *area_error;
} ImageInfoKist;


typedef struct Grid{
	TreeHndl i_tree; // tree on image plane
	TreeHndl s_tree; // tree on source plane
	Boolean initialized;
} Grid;

typedef struct Grid *GridHndl;

#endif

/***** Constructors/Destructors*****/

TreeHndl NewTree(Point *xp,unsigned long npoints
		 ,double boundery_p1[2],double boundery_p2[2],
		 double center[2]);

/***** Access functions *****/

Boolean isEmpty(TreeHndl tree);
Boolean atTop(TreeHndl tree);
Boolean offEnd(TreeHndl tree);
Boolean noChild(TreeHndl tree);

/*unsigned long *getCurrent(TreeHndl tree,unsigned long *npoints);*/
void getCurrent(TreeHndl tree,Point *points,unsigned long *npoints);
unsigned long getNbranches(TreeHndl tree);

/***** Manipulation procedures *****/

void moveTop(TreeHndl tree);
void moveUp(TreeHndl tree);

Boolean moveToChild(TreeHndl tree,int child);

void insertChildToCurrent(TreeHndl tree, Point *points,unsigned long npoints
			  ,double boundery_p1[2],double boundery_p2[2]
			  ,double center[2],int child);

void attachChildToCurrent(TreeHndl tree,Branch data,int child);
void attachChildrenToCurrent(TreeHndl tree,Branch child1,Branch child2);
Boolean TreeWalkStep(TreeHndl tree,Boolean allowDescent);
double ClosestBorder(double *ray,double *p1,double *p2);
inline double FurthestBorder(double *ray,double *p1,double *p2);
void PointsInCurrent(TreeHndl tree,unsigned long *ids,double **x);

/***** Other operations *****/

void printTree(TreeHndl tree);
void printBranch(Branch *branch);

void saveTree(TreeHndl tree,char *filename);
TreeHndl readTree(char *filename);
short emptyTree(TreeHndl tree);
short freeTree(TreeHndl tree);
void _freeTree(TreeHndl tree,short child);
void _freeTree_iter(TreeHndl tree);
void checkTree(TreeHndl tree);

/** routines in TreeDriver.c **/

TreeHndl BuildTree(Point *xp,unsigned long Npoints);
void _BuildTree(TreeHndl tree);
void FillTree(TreeHndl tree,Point *xp,unsigned long Npoints);
Point *NearestNeighbor(TreeHndl tree,double *ray,int Nneighbors,ListHndl neighborlist
		,short direction);
inline int inbox(double ray[2],double *p1,double *p2);
Boolean boxinbox(Branch *branch1,Branch *branch2);
double BoxIntersection(Branch *branch1,Branch *branch2);
int cutbox(double ray[2],double *p1,double *p2,double rmax);
void FindBoxPoint(TreeHndl tree,double *ray,Point *point);
void _FindBox(TreeHndl tree,double *ray);
Boolean AreBoxNeighbors(Point *point1,Point *point2);


// Point arrays

void PrintPoint(Point *point);
Point *NewPointArray(unsigned long N,Boolean NewXs);
Point *AddPointToArray(Point *points,unsigned long N,unsigned long Nold);
void FreePointArray(Point *array);
void SwapPointsInArray(Point *p1,Point *p2);
void PointCopy(Point *pcopy,Point *pins);
void PointCopyData(Point *pcopy,Point *pins);

// image info routines

ImageInfo *NewImageInfo(int Nimages);
void freeImageInfo(ImageInfo *imageinfo,int Nimages);
void combineCloseImages(double linkinglength,ImageInfo *imageinfo,int *Nimages
		,int *NewNimages);
void SwapImages(ImageInfo *image1,ImageInfo *image2);
void PrintImages(ImageInfo *images,long Nimages);
void PrintImageInfo(ImageInfo *image);

// routines using tree

int AddPointsToTree(TreeHndl tree,Point *xpoint,unsigned long Nadd);
void PointsWithin(TreeHndl tree,double *ray,float rmax,ListHndl neighborlist,short markpoints);
void PointsWithin_iter(TreeHndl tree,double *ray,float rmax,ListHndl neighborlist,short markpoints);
void NeighborsOfNeighbors(ListHndl neighbors,ListHndl wholelist);
void FriendsOfFriends(TreeHndl tree,double *starting_point,float rlink
		      ,ListHndl neighborlist,Point *filter
		      ,unsigned long Nfilter,unsigned long *filter_place);
void _PointsWithin2(TreeHndl tree,double *ray,float *rmax,ListHndl neighborlist
		   ,Point *filter,unsigned long Nfilter
		   ,unsigned long *filter_place,short compliment);
void _PointsWithin(TreeHndl tree,double *ray,float *rmax,ListHndl neighborlist
		,short markpoints);
void ClearAllMarkes(TreeHndl tree);
void FindAllBoxNeighbors(TreeHndl tree,Point *point,ListHndl neighbors);

// in image_finder.c

void find_images(double *y_source,double r_source,TreeHndl s_tree,TreeHndl i_tree
		,int *Nimages,ImageInfo *imageinfo,const int NimageMax,unsigned long *Nimagepoints
		  ,double initial_size,Boolean splitimages,short edge_refinement
		  ,Boolean verbose,Boolean kappa_off);
short image_finder(double *y_source,double r_source,TreeHndl s_tree,TreeHndl i_tree
		,int *Nimages,ImageInfo *imageinfo,const int NimageMax,unsigned long *Nimagepoints
		,short splitparities,short true_images);
int refine_grid(TreeHndl i_tree,TreeHndl s_tree,ImageInfo *imageinfo
		,unsigned long Nimages,double res_target,short criterion,Boolean kappa_off);
int refine_grid2(TreeHndl i_tree,TreeHndl s_tree,ImageInfo *imageinfo
		,unsigned long Nimages,double res_target,short criterion,Boolean kappa_off,Boolean shootrays,Point **i_points);
long refine_edges(TreeHndl i_tree,TreeHndl s_tree,ImageInfo *imageinfo
		,unsigned long Nimages,double res_target,short criterion,Boolean kappa_off);
long refine_edges2(double *y_source,double r_source,TreeHndl i_tree,TreeHndl s_tree
		,ImageInfo *imageinfo,Boolean *image_overlap,unsigned long Nimages,double res_target
		,short criterion,Boolean kappa_off);
void xygridpoints(Point *points,double range,double *center,long Ngrid
		,short remove_center);
void initialize_grid(double center[],double range,long Ngrid,TreeHndl s_tree,TreeHndl i_tree);

// in find_crit.c
void findborders(TreeHndl i_tree,ImageInfo *imageinfo);
void findborders2(TreeHndl i_tree,ImageInfo *imageinfo);
void findborders3(TreeHndl i_tree,ImageInfo *imageinfo);
void findborders4(TreeHndl i_tree,ImageInfo *imageinfo);
ImageInfo *find_crit(TreeHndl s_tree,TreeHndl i_tree,int *Ncrits,double resolution
		,Boolean *orderingsuccess,Boolean ordercurve,Boolean verbose);

/* in double_sort.c */
void double_sort(unsigned long n, double *arr, unsigned long *brr);
void double_sort_points(unsigned long n, double *arr, Point *brr);
void quicksortPoints(Point *pointarray,double *arr,unsigned long N);
void quicksort(unsigned long *particles,double *arr,unsigned long N);
void quickPartition(double pivotvalue,unsigned long *pivotindex,unsigned long *particles
		,double *arr,unsigned long N);
void quickPartitionPoints(double pivotvalue,unsigned long *pivotindex
		,Point *pointsarray,double *arr,unsigned long N);

/* in utilities.c */
inline double MIN(double x,double y);
inline double MAX(double x,double y);
Point *LinkToSourcePoints(Point *i_points,unsigned long Npoints);
void log_polar_grid(Point *i_points,double rmax,double rmin,double *center,long Ngrid);
void findarea(ImageInfo *imageinfo);
int windings(double *x,Point *points,unsigned long Npoints,double *area,short image);
long IndexFromPosition(double *x,long Npixels,double range,double *center);
void PositionFromIndex(unsigned long i,double *x,long Npixels,double range,double *center);
inline float isLeft( Point *p0, Point *p1, double *x );

// in curve_routines.c
void nesting_curve(ImageInfo *curves,int Ncurves);
void split_order_curve(ImageInfo *curves,int Maxcurves,int *Ncurves);
void split_order_curve2(ImageInfo *curves,int Maxcurves,int *Ncurves);
void split_order_curve3(ImageInfo *curves,int Maxcurves,int *Ncurves);
void split_order_curve4(ImageInfo *curves,int Maxcurves,int *Ncurves);
void walkcurve(Point *points,long Npoints,long *j,long *end);
short backtrack(Point *points,long Npoints,long *j,long jold,long *end);
void split_images(TreeHndl i_tree,ImageInfo *images,int Maximages,int *Nimages,Boolean sortallpoints);
void split_images2(TreeHndl i_tree,ImageInfo *images,int Maximages
		,int *Nimages);
void split_images3(TreeHndl i_tree,ImageInfo *images,int Maximages
		,int *Nimages,Boolean sortallpoints);
void splitter(ImageInfo *images,int Maximages,int *Nimages);
void splitlist(ListHndl imagelist,ImageInfo *images,int *Nimages,int Maximages);

/* externally provided functions */
/*********************************/

/*  void rayshooterInternal(double *x,double *alpha,double *gamma,double *kappa,double *invmag);*/
void rayshooterInternal(unsigned long Npoints,Point *i_points,TreeHndl i_tree
		,Boolean kappa_off);
void in_source(double *y_source,ListHndl sourcelist);


