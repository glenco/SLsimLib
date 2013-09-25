 /*
 * Code Name:     tree.h                                       
 * Programmer:    Ben Metcalf
 * Description:
 * Comments:  This version uses a linked list for the points so that
 the tree can be expended dynamically
 */

#ifndef treetypes_declare
#define treetypes_declare

//#include "pointlist.h"
#include "point.h"
#include "Kist.h"
#include "image_info.h"

/***** Exported Types *****/

#ifndef criterion_declare
#define criterion_declare
typedef enum{TotalArea,EachImage,Resolution,FillHoles} ExitCriterion;
#endif

/** \brief Tree: Exported struct */
struct TreeStruct{
public:
	TreeStruct(Point *xp,unsigned long Npoints,short my_median_cut = 1,double buffer = 0.0);
	TreeStruct(Point *xp,unsigned long npoints
			 ,double boundary_p1[2],double boundary_p2[2]
			 ,double center[2],int Nbucket);

	 ~TreeStruct();

	 /// root branch
	 Branch *top;
	 Branch *current;
	 /// list of points
	 PointList *pointlist;

	 /*  *** Constructor****
	TreeHndl NewTree(Point *xp,unsigned long npoints
			 ,double boundary_p1[2],double boundary_p2[2]
			 ,double center[2],int Nbucket);
	short freeTree(TreeHndl tree);*/

  //void FindAllBoxNeighbors(Point *point,ListHndl neighbors);
  void FindAllBoxNeighborsKist(Point *point,Kist<Point> * neighbors);
  void PointsWithinEllipKist(double *ray,float rmax,float rmin,float posangle,Kist<Point> * neighborkist);
  double PointsWithinKist(double *ray,float rmax,Kist<Point> * neighborkist,short markpoints);
  void PointsWithinKist_iter(double *ray,float rmin,float rmax,Kist<Point> * neighborkist);
  Point *NearestNeighborKist(double *ray,int Nneighbors,Kist<Point> * neighborkist);

  void PointsInCurrent(unsigned long *ids,double **x);

  /***** Movement on tree *****/

  void moveTop(){ current = top;}
  bool moveUp();
  bool moveToChild(int child);
  bool TreeWalkStep(bool allowDescent);

  // Adding and removing to branches of tree
  void insertChildToCurrent(Branch *branch,int child);
  void attachChildrenToCurrent(Branch *child1,Branch *child2);
  Point *RemoveLeafFromTree(unsigned long *Npoints);

  // Higher level builds
  void FillTree(Point *xp,unsigned long Npoints);
  int AddPointsToTree(Point *xpoint,unsigned long Nadd);

  short emptyTree();
  short freeTree();
  void RebuildTreeFromList();

  /***** State of tree functions *****/
  bool isEmpty();
  bool atTop();
  bool atLeaf(){ return( (current->child1==NULL)*(current->child2==NULL) ); }
  bool offEnd();
  bool CurrentIsSquareBranch();
  bool noChild();

  void getCurrent(Point *points,unsigned long *npoints);
  unsigned long getNbranches();
  void printTree();
  void checkTree();

  void FindBoxPoint(double *ray,Point *point);

  void _FindLeaf(double *ray,unsigned long Nadd = 0);

private:

  /// number of barnches in tree */
  unsigned long Nbranches;
   /// number of points allowed in leaves of tree
  int Nbucket;
  short median_cut;

  void construct_root(Point *xp,unsigned long npoints
			 ,double boundary_p1[2],double boundary_p2[2]
			 ,double center[2],int Nbucket);


  //void _FindAllBoxNeighbors(Branch *leaf,ListHndl neighbors);
  void _FindAllBoxNeighborsKist(Branch *leaf,Kist<Point> * neighbors);
  void _FindAllBoxNeighborsKist_iter(Branch *leaf,Kist<Point> * neighbors);
  void _PointsWithinKist(double *ray,float *rmax,Kist<Point> * neighborkist
  		,short markpoints,double *maxgridsize);

  void _freeBranches(short child);
  void _freeBranches_iter();
  void _AddPoint();
  void _BuildTree();

  void _checkTree(unsigned long *count);
  void _FindBox(double *ray);

  // Should be obsolete
  Point *NearestNeighbor(double *ray,int Nneighbors,ListHndl neighborlist
  		,short direction);
  void _NearestNeighbor(double *ray,int Nneighbors,Point **neighborpoints,double *rneighbors,short *direction);


  // Are obsolete
  //void PointsWithin(double *ray,float rmax,ListHndl neighborlist,short markpoints);
  //void PointsWithin_iter(double *ray,float rmax,ListHndl neighborlist,short markpoints);
  //void FriendsOfFriends(double *starting_point,float rlink
  //		      ,ListHndl neighborlist,Point *filter
  //		      ,unsigned long Nfilter,unsigned long *filter_place);
  //void _PointsWithin2(double *ray,float *rmax,ListHndl neighborlist
  //		   ,Point *filter,unsigned long Nfilter
  //		   ,unsigned long *filter_place,short compliment);
  //void _PointsWithin(double *ray,float *rmax,ListHndl neighborlist
  //		,short markpoints);

};

typedef struct TreeStruct *TreeHndl;
typedef int TreeElement;

bool BoxInCircle(double *ray,double radius,double *p1,double *p2);
double ClosestBorder(double *ray,double *p1,double *p2);

inline double MIN(double x,double y){
	return (x < y) ? x : y;
};
inline double MAX(double x,double y){
	return (x > y) ? x : y;
};


/*  returns the distance from ray[] to the furthest point on the
 *    border of the box,
 */
inline double FurthestBorder(double *ray,double *p1,double *p2){
  return sqrt( pow(MAX(ray[0]-p1[0],p2[0]-ray[0]),2) + pow(MAX(ray[1]-p1[1],p2[1]-ray[1]),2) );
};

/***** Other operations *****/
void printBranch(Branch *branch);

void saveTree(TreeHndl tree,char *filename);
TreeHndl readTree(char *filename);

/** routines in TreeDriver.c **/

//inline int inbox(double ray[2],double *p1,double *p2);
/* return 1 (0) if ray is (not) in the cube */
inline int inbox(double *ray,double *p1,double *p2){
  return (ray[0]>=p1[0])*(ray[0]<=p2[0])*(ray[1]>=p1[1])*(ray[1]<=p2[1]);
};
bool boxinbox(Branch *branch1,Branch *branch2);
double BoxIntersection(Branch *branch1,Branch *branch2);
bool AreBoxNeighbors(Point *point1,Point *point2);
bool AreBoxNeighbors(Branch *branch1,Branch *branch2);
bool CircleInBox(double *ray,double radius,double *p1,double *p2);
bool BoxInCircle(double *ray,double radius,double *p1,double *p2);

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

Point *LinkToSourcePoints(Point *i_points,unsigned long Npoints);

/// \ingroup Util
namespace Utilities{
    ///Separation squared between two positions in 2 dimensions.
	inline double sepSQR(double *xx,double *yy){
		return pow(xx[0]-yy[0],2) + pow(xx[1]-yy[1],2);
	}
	void double_sort(unsigned long n, double *arr, unsigned long *brr);
	void double_sort_points(unsigned long n, double *arr, Point *brr);

	void quicksortPoints(Point *pointarray,double *arr,unsigned long N);
	void quicksort(unsigned long *particles,double *arr,unsigned long N);
	void quickPartition(double pivotvalue,unsigned long *pivotindex,unsigned long *particles
		,double *arr,unsigned long N);
	void quickPartitionPoints(double pivotvalue,unsigned long *pivotindex
		,Point *pointsarray,double *arr,unsigned long N);
	int cutbox(PosType *ray,PosType *p1,PosType *p2,float rmax);
	void log_polar_grid(Point *i_points,double rmax,double rmin,double *center,long Ngrid);
	void findarea(ImageInfo *imageinfo);
	int windings2(double *x,Point *points,unsigned long Npoints,double *area,short image);
	void writeCurves(int m, ImageInfo *critical, int Ncrit, int index);
  double cross(const Point *O, const Point *A, const Point *B);
  bool xorder(Point *p1,Point *p2);
  std::vector<Point *> convex_hull(std::vector<Point *> P);


	long IndexFromPosition(double *x,long Npixels,double range,double *center);
	void PositionFromIndex(unsigned long i,double *x,long Npixels,double range,double *center);
	long IndexFromPosition(double x,long Npixels,double range,double center);
  double TwoDInterpolator(double *x,int Npixels,double range,double *center,double *map,bool init=true);
  double TwoDInterpolator(double *map);

  /** \ingroup Utill
   * \brief Bilinear interpolation class for interpolating from a 2D uniform grid.
   *
   *  Out of bounds points return 0.  map is a i dimensional array representing a 2 dimensional map.
   *  
   *  Later calls can use interpolator(map) for the same point
   *  in the same coordinate system to save time in calculating the indexes.
   */

  template <typename T>
  class Interpolator{
  public:
    Interpolator(
                 double *x          /// position of point
                 ,int Npixels       /// Number of pixels in one dimension
                 ,double my_range   /// Range of map in same units as x[]
                 ,double *my_center /// Center of map in same units as x[]
                 ):
    N(Npixels),range(my_range),map_p(NULL)
    {
      center[0] = my_center[0];
      center[1] = my_center[1];
      initparams(x);
    };
        
    /**
     This constructor takes the map as a pointer to an array of values and stores it.
     The resulting object can then be used with the () operator as a function.
     Warning: Be sure to distroy the object before distroying map.
     */
    Interpolator(
                 int Npixels          /// Number of pixels in one dimension
                 ,double my_range     /// Range of map in same units as x[]
                 ,double *my_center   /// Center of map in same units as x[]
                 ,const T *map        /// One dimensional array of fundamental type
                 ):
    N(Npixels),range(my_range),map_p(map)
    {
      center[0] = my_center[0];
      center[1] = my_center[1];
    };
    
    
    /** 
     Does interpolation of map at point that object was constructed with or last called with.
     Can use and map type that has a [] operator that returns a double.
     */
    double interpolate(
                       T& map    /// map that supports the [] operator 
                       ){
      if(map.size() != N*N){
        ERROR_MESSAGE();
        std::cout << "ERROR: Interpolator:interpolator(double *,T&)" << std::endl;
      }
      if(index == -1) return 0;
      
      return (1.-fx)*(1.-fy)*map[index] + fx*(1.-fy)*map[index+1] + fx*fy*map[index+1+N]
          + (1.-fx)*fy*map[index+N];
    };
    /// reinitializes to a new position
    double interpolate(
                       double *x   /// position of point
                       ,T& map     /// map that supports the [] operator
                       ){
      if(map.size() != N*N){
        ERROR_MESSAGE();
        std::cout << "ERROR: Interpolator:interpolator(double *,T&)" << std::endl;
      }
      initparams(x);
      return interpolate(map);
    }
    
    /**
     Does interpolation of store map at point x. Only for use with the second constructor.
     */
    double operator ()(double *x){
      if(map_p == NULL){
        std::cout << "Didn't use the right constructor for Interpolator class" << std::endl;
        throw std::runtime_error("Did not use the right constructor.");
        return 0.0;
      }
      initparams(x);
      return (1-fx)*(1-fy)*map_p[index] + fx*(1-fy)*map_p[index+1] + fx*fy*map_p[index+1+N]
      + (1-fx)*fy*map_p[index+N];
    }
        
    void test(void){
      std::valarray<float> map;
      double tmp,x[2];
          
      map.resize(N*N);
      
      for(int i=0;i<N;++i){
        tmp = cos( 6.*i*2.*pi/(N-1) );
        for(int j=0;j<N;++j){
            map[i+N*j] = tmp;
        }
      }

      x[1] = 0*range/2.;
      for(int i=0;i<N;++i){
        x[0] = x[1] = center[0] + range*( 1.0*(i)/(N-1) - 0.5 ) + 0.5*range/(N-1);
        std::cout << i << "  " << map[i+N*i] << " " << interpolate(x,map) << std::endl;
      }
    }

  private:
    double range,center[2];
    const T *map_p;
    double fx, fy;
    long index;
    int N;

    void initparams(double *x){
      long ix,iy;
      // position in pixel coordinates
      fx = ((x[0] - center[0])/range + 0.5)*(N-1);
      fy = ((x[1] - center[1])/range + 0.5)*(N-1);
      //std::cout << "(  " << fx << " " << fy << "   ";
     
      if (fx < 0. || fx > N-1){index = -1; return;}
      else ix = (unsigned long)(fx);
      if(ix == N-1) ix = N-2;
      
      if (fy < 0. || fy > N-1){index = -1; return;}
      else iy = (unsigned long)(fy);
      if(iy == N-1) iy = N-2;
      
      index = ix + N*iy;

      //std::cout << "  " << ix << " " << iy << " " << index << "   ";
      
      /** bilinear interpolation */
      fx = center[0] + range*( 1.0*(ix)/(N-1) - 0.5 );
      fy = center[1] + range*( 1.0*(iy)/(N-1) - 0.5 );
      fx=(x[0] - fx)*(N-1)/range;
      fy=(x[1] - fy)*(N-1)/range;
      //std::cout  << fx << " " << fy << "   )";
     
      /*
      //printf("%i %i\n",i,j);
      if(j > nx-1 || i > ny-1) return 0.0;
      if(j==0 || i==0 ) return 0.0;
      
      // bilinear interpolation 
      t=(x-xtab[j])/(xtab[j+1]-xtab[j]);
      u=(y-ytab[i])/(ytab[i+1]-ytab[i]);
      
      return (1-t)*(1-u)*tab[i][j] + t*(1-u)*tab[i][j+1] + t*u*tab[i+1][j+1]
      + (1-t)*u*tab[i+1][j];
      */
    }
        
   };

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

	int windings(double *x,Point *points,unsigned long Npoints,double *area,short image = 0 );
	int windings(double *x,Kist<Point> * kist,double *area,short image = 0);
}
// in curve_routines.c
void nesting_curve(OldImageInfo *curves,int Ncurves);
void split_order_curve(OldImageInfo *curves,int Maxcurves,int *Ncurves);
void split_order_curve2(OldImageInfo *curves,int Maxcurves,int *Ncurves);
void split_order_curve3(OldImageInfo *curves,int Maxcurves,int *Ncurves);
void split_order_curve4(OldImageInfo *curves,int Maxcurves,int *Ncurves);
namespace Utilities{
	unsigned long order_curve4(Point *curve,long Npoints);
	unsigned long order_curve4(Kist<Point> * curve);
  unsigned long order_curve5(Kist<Point> * curve);
  void ordered_convexhull(Kist<Point> * curve);
}
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


/* externally provided functions */
/*********************************/

/*  void rayshooterInternal(double *x,double *alpha,double *gamma,double *kappa,double *invmag);*/
void rayshooterInternal(unsigned long Npoints,Point *i_points,bool kappa_off);
void in_source(double *y_source,ListHndl sourcelist);
bool tree_count_test(TreeHndl tree);
bool testLeafs(TreeHndl tree);

void NeighborsOfNeighbors(ListHndl neighbors,ListHndl wholelist);

#endif
