 /*
 * Code Name:     tree.h                                       
 * Programmer:    Ben Metcalf
 * Description:
 * Comments:  This version uses a linked list for the points so that
 the tree can be expended dynamically
 */

#ifndef treetypes_declare
#define treetypes_declare

#include "pointlist.h"
#include "Kist.h"
#include "image_info.h"

/***** Exported Types *****/

#ifndef criterion_declare
#define criterion_declare
typedef enum{TotalArea,EachImage,Resolution,FillHoles} ExitCriterion;
#endif

/** \brief Tree: Exported struct */
typedef struct TreeStruct{
  Branch *top;
  Branch *current;
  /// number of barnches in tree */
  unsigned long Nbranches;
  /// list of points
  PointList *pointlist;
  /// number of points allowed in leaves of tree
  int Nbucket;
  short median_cut;
} TreeStruct;

typedef struct TreeStruct *TreeHndl;
typedef int TreeElement;


/*  *** Constructor****/
TreeHndl NewTree(Point *xp,unsigned long npoints
		 ,double boundary_p1[2],double boundary_p2[2]
		 ,double center[2],int Nbucket);
short freeTree(TreeHndl tree);

/***** Access functions *****/

bool isEmpty(TreeHndl tree);
bool atTop(TreeHndl tree);
//inline bool atLeaf(TreeHndl tree);
inline bool atLeaf(TreeHndl tree){
  return( (tree->current->child1==NULL)*(tree->current->child2==NULL) );
};
bool BoxInCircle(double *ray,double radius,double *p1,double *p2);
bool offEnd(TreeHndl tree);
bool CurrentIsSquareTree(TreeHndl tree);
bool noChild(TreeHndl tree);

/*unsigned long *getCurrent(TreeHndl tree,unsigned long *npoints);*/
void getCurrent(TreeHndl tree,Point *points,unsigned long *npoints);
unsigned long getNbranches(TreeHndl tree);

/***** Manipulation procedures *****/

void moveTop(TreeHndl tree);
bool moveUp(TreeHndl tree);

bool moveToChild(TreeHndl tree,int child);

void insertChildToCurrent(TreeHndl tree, Point *points,unsigned long npoints
			  ,double boundary_p1[2],double boundary_p2[2]
			  ,double center[2],int child);

void attachChildToCurrent(TreeHndl tree,Branch data,int child);
void attachChildrenToCurrent(TreeHndl tree,Branch child1,Branch child2);
bool TreeWalkStep(TreeHndl tree,bool allowDescent);
double ClosestBorder(double *ray,double *p1,double *p2);

inline double MIN(double x,double y){
	return (x < y) ? x : y;
};
inline double MAX(double x,double y){
	return (x > y) ? x : y;
};

/** \ingroup Utill
 * \brief Separation squared between two positions in 2 dimensions.
 */
inline double sepSQR(double *xx,double *yy){
	return pow(xx[0]-yy[0],2) + pow(xx[1]-yy[1],2);
}

/*  returns the distance from ray[] to the furthest point on the
 *    border of the box,
 */
inline double FurthestBorder(double *ray,double *p1,double *p2){
  return sqrt( pow(MAX(ray[0]-p1[0],p2[0]-ray[0]),2) + pow(MAX(ray[1]-p1[1],p2[1]-ray[1]),2) );
};
void PointsInCurrent(TreeHndl tree,unsigned long *ids,double **x);

/***** Other operations *****/

void printTree(TreeHndl tree);
void printBranch(Branch *branch);

void saveTree(TreeHndl tree,char *filename);
TreeHndl readTree(char *filename);
void checkTree(TreeHndl tree);

/** routines in TreeDriver.c **/

Point *NearestNeighbor(TreeHndl tree,double *ray,int Nneighbors,ListHndl neighborlist
		,short direction);
//inline int inbox(double ray[2],double *p1,double *p2);
/* return 1 (0) if ray is (not) in the cube */
inline int inbox(double *ray,double *p1,double *p2){
  return (ray[0]>=p1[0])*(ray[0]<=p2[0])*(ray[1]>=p1[1])*(ray[1]<=p2[1]);
};
bool boxinbox(Branch *branch1,Branch *branch2);
double BoxIntersection(Branch *branch1,Branch *branch2);
void FindBoxPoint(TreeHndl tree,double *ray,Point *point);
void _FindBox(TreeHndl tree,double *ray);
bool AreBoxNeighbors(Point *point1,Point *point2);
bool CircleInBox(double *ray,double radius,double *p1,double *p2);
bool BoxInCircle(double *ray,double radius,double *p1,double *p2);
void _FindLeaf(TreeHndl tree,double *ray,unsigned long Nadd);

// Point arrays

void PrintPoint(Point *point);
Point *NewPointArray(unsigned long N,bool NewXs);
Point *AddPointToArray(Point *points,unsigned long N,unsigned long Nold);
void FreePointArray(Point *array,bool NewXs = true);
void SwapPointsInArray(Point *p1,Point *p2);
void SwapPointsInArrayData(Point *p1,Point *p2);
void PointCopy(Point *pcopy,Point *pins);
void PointCopyData(Point *pcopy,Point *pins);
void PointCopyAll(Point *pcopy,Point *pins);

// image info routines

//ImageInfo *NewImageInfo(int Nimages);
//void freeImageInfo(ImageInfo *imageinfo,int Nimages);
void combineCloseImages(double linkinglength,ImageInfo *imageinfo,int *Nimages
		,int *NewNimages,int NimageMax);
void SwapImages(ImageInfo *image1,ImageInfo *image2);
void SwapImages(OldImageInfo *image1,OldImageInfo *image2);
void PrintImages(ImageInfo *images,long Nimages);
void PrintImageInfo(ImageInfo *image);

// routines using tree

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
void FindAllBoxNeighbors(TreeHndl tree,Point *point,ListHndl neighbors);

// in image_finder.c

/*void find_images(double *y_source,double r_source,GridHndl grid
		,int *Nimages,ImageInfo *imageinfo,const int NimageMax,unsigned long *Nimagepoints
		,double initial_size,bool splitimages,short edge_refinement
		,bool verbose,bool kappa_off);
short image_finder(double *y_source,double r_source,TreeHndl s_tree,TreeHndl i_tree
		,int *Nimages,ImageInfo *imageinfo,const int NimageMax,unsigned long *Nimagepoints
		,short splitparities,short true_images);*/
void findborders2(TreeHndl i_tree,OldImageInfo *imageinfo);
void findborders3(TreeHndl i_tree,OldImageInfo *imageinfo);

// in image_finder_kist.c

void findborders4(TreeHndl i_tree,ImageInfo *imageinfo);

// in find_crit.c
void findborders(TreeHndl i_tree,ImageInfo *imageinfo);

/* in double_sort.c */
void double_sort(unsigned long n, double *arr, unsigned long *brr);
void double_sort_points(unsigned long n, double *arr, Point *brr);
void quicksortPoints(Point *pointarray,double *arr,unsigned long N);
void quicksort(unsigned long *particles,double *arr,unsigned long N);
void quickPartition(double pivotvalue,unsigned long *pivotindex,unsigned long *particles
		,double *arr,unsigned long N);
void quickPartitionPoints(double pivotvalue,unsigned long *pivotindex
		,Point *pointsarray,double *arr,unsigned long N);
int cutbox(PosType *ray,PosType *p1,PosType *p2,float rmax);

/* in utilities.c */

Point *LinkToSourcePoints(Point *i_points,unsigned long Npoints);
void log_polar_grid(Point *i_points,double rmax,double rmin,double *center,long Ngrid);
void findarea(ImageInfo *imageinfo);
int windings2(double *x,Point *points,unsigned long Npoints,double *area,short image);
void writeCurves(int m, ImageInfo *critical, int Ncrit, int index);

long IndexFromPosition(double *x,long Npixels,double range,double *center);
void PositionFromIndex(unsigned long i,double *x,long Npixels,double range,double *center);
int IndexFromPosition(double x,long Npixels,double range,double center);
//inline float isLeft( Point *p0, Point *p1, double *x );

// isLeft(): tests if a point is Left|On|Right of an infinite line.
// Input:three points P0, P1, and x
// Return: >0 for x left of the line through P0 and P1
//         =0 for x on the line
//         <0 for x right of the line
inline float isLeft( Point *p0, Point *p1, double *x ){

	return (p1->x[0] - p0->x[0])*(x[1] - p0->x[1])
			- (x[0] - p0->x[0])*(p1->x[1] - p0->x[1]);
};
unsigned long prevpower(unsigned long k);

// in curve_routines.c
void nesting_curve(OldImageInfo *curves,int Ncurves);
void split_order_curve(OldImageInfo *curves,int Maxcurves,int *Ncurves);
void split_order_curve2(OldImageInfo *curves,int Maxcurves,int *Ncurves);
void split_order_curve3(OldImageInfo *curves,int Maxcurves,int *Ncurves);
void split_order_curve4(OldImageInfo *curves,int Maxcurves,int *Ncurves);
unsigned long order_curve4(Point *curve,long Npoints);
bool order_curve4(KistHndl curve);
bool order_ExteriorBoundary(Point *curve,long Npoints,long *NewNpoints,double *area);
double findAreaOfCurve(TreeHndl tree,ImageInfo *curve,int NimageMax);
void walkcurve(Point *points,long Npoints,long *j,long *end);
void walkcurveRight(Point *points,long Npoints,long *j,long *end);
short backtrack(Point *points,long Npoints,long *j,long jold,long *end);
void split_images(TreeHndl i_tree,ImageInfo *images,int Maximages,int *Nimages,bool sortallpoints);
void split_images2(TreeHndl i_tree,ImageInfo *images,int Maximages,int *Nimages);
void split_images3(TreeHndl i_tree,ImageInfo *images,int Maximages,int *Nimages,bool sortallpoints);
void splitter(OldImageInfo *images,int Maximages,int *Nimages);
void splitlist(ListHndl imagelist,OldImageInfo *images,int *Nimages,int Maximages);
int windings(double *x,Point *points,unsigned long Npoints,double *area,short image);
int windings(double *x,KistHndl kist,double *area,short image);

/* externally provided functions */
/*********************************/

/*  void rayshooterInternal(double *x,double *alpha,double *gamma,double *kappa,double *invmag);*/
void rayshooterInternal(unsigned long Npoints,Point *i_points,bool kappa_off);
void in_source(double *y_source,ListHndl sourcelist);
bool tree_count_test(TreeHndl tree);
bool testLeafs(TreeHndl tree);

#endif
