 /*
 * Code Name:     tree.h                                       
 * Programmer:    Ben Metcalf
 * Description:
 * Comments:  This version uses a linked list for the points so that
 the tree can be expended dynamically
 */
#include <math.h>
#include <point.h>
#include <Kist.h>
#include <List1.h>
#include <lens.h>

/***** Exported Types *****/

#ifndef error_message
#define error_message
#define ERROR_MESSAGE() printf("ERROR: file: %s line: %i\n",__FILE__,__LINE__)
#endif

#ifndef line_message
#define line_message
#define PRINT_LINE() printf("file: %s line: %i\n",__FILE__,__LINE__)
#endif

#ifndef criterion_declare
#define criterion_declare
typedef enum{TotalArea,EachImage,Resolution,FillHoles} ExitCriterion;
#endif

#ifndef swap_declare
#define swap_declare
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#endif

#ifndef treetypes_declare
#define treetypes_declare

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
} TreeStruct;

typedef struct TreeStruct *TreeHndl;
typedef int TreeElement;

#include <grid_maintenance.h>

/** \brief Structure for storing information about images or curves */
typedef struct ImageInfo{

    /// Array of points in image,  SHOULD NOT BE USED IN FAVOR OF imagekist!  Still used by caustic finding routines.
  Point *points;
  /// Number of points in image, SHOULD NOT BE USED IN FAVOR OF imagekist->Nunits().  Still used by caustic finding routines.
  unsigned long Npoints;
  /// later addition, holds all points in image, will replace points eventually
  KistHndl imagekist;
  /// gridrange[2] minimum grid size in image, gridrange[0] maximum grid size in outerborder, gridrange[1] maximum grid size in image
  double gridrange[3];
  /// Centroid of image
  double centroid[2];
  /// area of image or, when using map_images(), the total brightness of the image
  double area;
  /// error on the estimate of area
  double area_error;
  /// the points on the inner border of the image
  KistHndl innerborder;
  /// the points on the outer border of the image, i.e. not in the image
  KistHndl outerborder;
  short Nencircled;

} ImageInfo;


//#include <point.h>
//#include <List.h>
//#include <Kist.h>
//#include <KistDriver.h>

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
int cutbox(double ray[2],double *p1,double *p2,double rmax);
void FindBoxPoint(TreeHndl tree,double *ray,Point *point);
void _FindBox(TreeHndl tree,double *ray);
bool AreBoxNeighbors(Point *point1,Point *point2);

// Point arrays

void PrintPoint(Point *point);
Point *NewPointArray(unsigned long N,bool NewXs);
Point *AddPointToArray(Point *points,unsigned long N,unsigned long Nold);
void FreePointArray(Point *array);
void SwapPointsInArray(Point *p1,Point *p2);
void PointCopy(Point *pcopy,Point *pins);
void PointCopyData(Point *pcopy,Point *pins);

// image info routines

ImageInfo *NewImageInfo(int Nimages);
void freeImageInfo(ImageInfo *imageinfo,int Nimages);
void combineCloseImages(double linkinglength,ImageInfo *imageinfo,int *Nimages
		,int *NewNimages,int NimageMax);
void SwapImages(ImageInfo *image1,ImageInfo *image2);
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
void ClearAllMarks(TreeHndl tree);
void FindAllBoxNeighbors(TreeHndl tree,Point *point,ListHndl neighbors);

// in image_finder.c

/*void find_images(double *y_source,double r_source,GridHndl grid
		,int *Nimages,ImageInfo *imageinfo,const int NimageMax,unsigned long *Nimagepoints
		,double initial_size,bool splitimages,short edge_refinement
		,bool verbose,bool kappa_off);
short image_finder(double *y_source,double r_source,TreeHndl s_tree,TreeHndl i_tree
		,int *Nimages,ImageInfo *imageinfo,const int NimageMax,unsigned long *Nimagepoints
		,short splitparities,short true_images);*/
int refine_grid(LensHndl lens,TreeHndl i_tree,TreeHndl s_tree,ImageInfo *imageinfo
		,unsigned long Nimages,double res_target,short criterion,bool kappa_off);
long refine_edges(LensHndl lens,TreeHndl i_tree,TreeHndl s_tree,ImageInfo *imageinfo
		,unsigned long Nimages,double res_target,short criterion,bool kappa_off);
long refine_edges2(LensHndl lens,double *y_source,double r_source,TreeHndl i_tree,TreeHndl s_tree
		,ImageInfo *imageinfo,bool *image_overlap,unsigned long Nimages,double res_target
		,short criterion,bool kappa_off);
void xygridpoints(Point *points,double range,double *center,long Ngrid
		,short remove_center);
void initialize_grid(double center[],double range,long Ngrid,TreeHndl s_tree,TreeHndl i_tree);
void findborders2(TreeHndl i_tree,ImageInfo *imageinfo);
void findborders3(TreeHndl i_tree,ImageInfo *imageinfo);

// in image_finder_kist.c

void find_images_kist(LensHndl lens,double *y_source,double r_source,GridHndl grid
		,int *Nimages,ImageInfo *imageinfo,const int NimageMax,unsigned long *Nimagepoints
		,double initial_size,bool splitimages,short edge_refinement
		,bool verbose,bool kappa_off);

short image_finder_kist(LensHndl lens,double *y_source,double r_source,TreeHndl s_tree,TreeHndl i_tree
		,int *Nimages,ImageInfo *imageinfo,const int NimageMax,unsigned long *Nimagepoints
		,short splitparities,short true_images);

int refine_grid_kist(LensHndl lens,TreeHndl i_tree,TreeHndl s_tree,ImageInfo *imageinfo
		,unsigned long Nimages,double res_target,short criterion,bool kappa_off,bool shootrays,Point **i_points);
void findborders4(TreeHndl i_tree,ImageInfo *imageinfo);

// in find_crit.c
void findborders(TreeHndl i_tree,ImageInfo *imageinfo);

ImageInfo *find_crit(LensHndl lens,GridHndl grid,int *Ncrits,double resolution,bool *orderingsuccess
		,bool ordercurve,bool verbose);

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

Point *LinkToSourcePoints(Point *i_points,unsigned long Npoints);
void log_polar_grid(Point *i_points,double rmax,double rmin,double *center,long Ngrid);
void findarea(ImageInfo *imageinfo);
int windings(double *x,Point *points,unsigned long Npoints,double *area,short image);
long IndexFromPosition(double *x,long Npixels,double range,double *center);
void PositionFromIndex(unsigned long i,double *x,long Npixels,double range,double *center);
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
void nesting_curve(ImageInfo *curves,int Ncurves);
void split_order_curve(ImageInfo *curves,int Maxcurves,int *Ncurves);
void split_order_curve2(ImageInfo *curves,int Maxcurves,int *Ncurves);
void split_order_curve3(ImageInfo *curves,int Maxcurves,int *Ncurves);
void split_order_curve4(ImageInfo *curves,int Maxcurves,int *Ncurves);
void walkcurve(Point *points,long Npoints,long *j,long *end);
short backtrack(Point *points,long Npoints,long *j,long jold,long *end);
void split_images(TreeHndl i_tree,ImageInfo *images,int Maximages,int *Nimages,bool sortallpoints);
void split_images2(TreeHndl i_tree,ImageInfo *images,int Maximages
		,int *Nimages);
void split_images3(TreeHndl i_tree,ImageInfo *images,int Maximages
		,int *Nimages,bool sortallpoints);
void splitter(ImageInfo *images,int Maximages,int *Nimages);
void splitlist(ListHndl imagelist,ImageInfo *images,int *Nimages,int Maximages);

/* externally provided functions */
/*********************************/

/*  void rayshooterInternal(double *x,double *alpha,double *gamma,double *kappa,double *invmag);*/
void rayshooterInternal(unsigned long Npoints,Point *i_points,bool kappa_off);
void in_source(double *y_source,ListHndl sourcelist);

#endif
