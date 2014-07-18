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
	TreeStruct(Point *xp,unsigned long Npoints,short my_median_cut = 1,PosType buffer = 0.0);
	TreeStruct(Point *xp,unsigned long npoints
			 ,PosType boundary_p1[2],PosType boundary_p2[2]
			 ,PosType center[2],int Nbucket);

	 ~TreeStruct();

	 /// root branch
	 Branch *top;
	 Branch *current;
	 /// list of points
	 PointList *pointlist;

	 /*  *** Constructor****
	TreeHndl NewTree(Point *xp,unsigned long npoints
			 ,PosType boundary_p1[2],PosType boundary_p2[2]
			 ,PosType center[2],int Nbucket);
	short freeTree(TreeHndl tree);*/

  //void FindAllBoxNeighbors(Point *point,ListHndl neighbors);
  void FindAllBoxNeighborsKist(Point *point,Kist<Point> * neighbors);
  void PointsWithinEllipKist(const PosType* center,float rmax,float rmin,float posangle,Kist<Point> * neighborkist);
  PosType PointsWithinKist(const PosType* center,PosType rmax,Kist<Point> * neighborkist,short markpoints);
  void PointsWithinKist_iter(const PosType* center,float rmin,float rmax,Kist<Point> * neighborkist);
  Point *NearestNeighborKist(const PosType* center,int Nneighbors,Kist<Point> * neighborkist);

  void PointsInCurrent(unsigned long *ids,PosType **x);

  /***** Movement on tree *****/

  void moveTop(){ current = top;}
  bool moveUp();
  bool moveToChild(int child);
  bool TreeWalkStep(bool allowDescent);
  bool Test();

  // Adding and removing to branches of tree
  void insertChildToCurrent(Branch *branch,int child);
  void attachChildrenToCurrent(Branch *child1,Branch *child2);
  Point *RemoveLeafFromTree(unsigned long *Npoints);

  // Higher level builds
  void FillTree(Point *xp,unsigned long Npoints);
  int AddPointsToTree(Point *xpoint,unsigned long Nadd);

  short emptyTree();
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

  void FindBoxPoint(const PosType* ray,Point *point);

  void _FindLeaf(const PosType* ray,unsigned long Nadd = 0);
  
  TreeStruct * spawn();

  /**
   *  \brief A iterator class fore TreeStruct that allows for movement through the tree without changing
   *      anything in the tree itself.
   *
   *   This class should be able to preform all of the constant movements within the tree without causing
   *   any change to the tree.
   */
  class iterator{
    
  private:
    Branch *current;
    Branch *top;
    
  public:
    /// Sets the top or root to the top of "tree".
    iterator(TreeStruct * tree){current = top = tree->top;}
    /// Sets the root to the input branch so that this will be a subtree in branch is not the real root.
    iterator(Branch *branch){current = top = branch;}
    
    /// Returns a pointer to the current Branch.
    Branch *operator*(){return current;}
    
    void movetop(){current = top;}
    
    /// Same as up()
    bool operator++(){ return up();}
    
    /// Same as up()
    bool operator++(int){ return up();}
    
    bool up();
    /// Move to brother if it exists
    bool brother();
    /// Move to child
    bool down(short child);
    bool TreeWalkStep(bool allowDescent);
  };

private:

  TreeStruct(){};

  Point **temp_points;
  std::vector<Point *> tmp_point;
  
  /// number of barnches in tree */
  unsigned long Nbranches;
   /// number of points allowed in leaves of tree
  int Nbucket;
  short median_cut;
  int incell;
  PosType realray[2];


  void construct_root(Point *xp,unsigned long npoints
			 ,PosType boundary_p1[2],PosType boundary_p2[2]
			 ,PosType center[2],int Nbucket);


  //void _FindAllBoxNeighbors(Branch *leaf,ListHndl neighbors);
  void _FindAllBoxNeighborsKist(Branch *leaf,Kist<Point> * neighbors);
  void _FindAllBoxNeighborsKist_iter(Branch *leaf,Kist<Point> * neighbors);
  void _PointsWithinKist(PosType *ray,PosType *rmax,Kist<Point> * neighborkist
  		,short markpoints,PosType *maxgridsize);

  void _freeBranches(short child);
  void _AddPoint();
  void _BuildTree();

  void _checkTree(unsigned long *count);
  void _freeBranches_iter();
  void _FindBox(const PosType* ray);

  // Should be obsolete
  //Point *NearestNeighbor(const PosType* center,int Nneighbors,ListHndl neighborlist
  //		,short direction);
  void _NearestNeighbor(PosType* ray,int Nneighbors,Point **neighborpoints,PosType *rneighbors,short *direction);


  // Are obsolete
  //void PointsWithin(PosType *ray,float rmax,ListHndl neighborlist,short markpoints);
  //void PointsWithin_iter(PosType *ray,float rmax,ListHndl neighborlist,short markpoints);
  //void FriendsOfFriends(PosType *starting_point,float rlink
  //		      ,ListHndl neighborlist,Point *filter
  //		      ,unsigned long Nfilter,unsigned long *filter_place);
  //void _PointsWithin2(PosType *ray,float *rmax,ListHndl neighborlist
  //		   ,Point *filter,unsigned long Nfilter
  //		   ,unsigned long *filter_place,short compliment);
  //void _PointsWithin(PosType *ray,float *rmax,ListHndl neighborlist
  //		,short markpoints);

};

typedef struct TreeStruct *TreeHndl;
typedef int TreeElement;

bool BoxInCircle(PosType *ray,PosType radius,PosType *p1,PosType *p2);
PosType ClosestBorder(PosType *ray,PosType *p1,PosType *p2);

inline PosType MIN(PosType x,PosType y){
	return (x < y) ? x : y;
};
inline PosType MAX(PosType x,PosType y){
	return (x > y) ? x : y;
};

template <class T>
inline bool BETWEEN(T x,T xmin,T xmax){
	return (x > xmin)*(x < xmax);
};

/*  returns the distance from ray[] to the furthest point on the
 *    border of the box,
 */
inline PosType FurthestBorder(const PosType* center,PosType *p1,PosType *p2){
  return sqrt( pow(MAX(center[0]-p1[0],p2[0]-center[0]),2) + pow(MAX(center[1]-p1[1],p2[1]-center[1]),2) );
};

/***** Other operations *****/
void printBranch(Branch *branch);

void saveTree(TreeHndl tree,char *filename);
TreeHndl readTree(char *filename);

/** routines in TreeDriver.c **/

//inline int inbox(PosType ray[2],PosType *p1,PosType *p2);
/* return 1 (0) if ray is (not) in the cube */
inline int inbox(const PosType* center,PosType *p1,PosType *p2){
  return (center[0]>=p1[0])*(center[0]<=p2[0])*(center[1]>=p1[1])*(center[1]<=p2[1]);
};
bool boxinbox(Branch *branch1,Branch *branch2);
PosType BoxIntersection(Branch *branch1,Branch *branch2);
bool AreBoxNeighbors(Point *point1,Point *point2);
bool AreBoxNeighbors(Branch *branch1,Branch *branch2);
bool CircleInBox(const PosType* center,PosType radius,PosType *p1,PosType *p2);
bool BoxInCircle(const PosType* center,PosType radius,PosType *p1,PosType *p2);

// Point arrays

void PrintPoint(Point *point);
Point *NewPointArray(unsigned long N);
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
void combineCloseImages(PosType linkinglength,ImageInfo *imageinfo,int *Nimages
		,int *NewNimages,int NimageMax);
void SwapImages(ImageInfo *image1,ImageInfo *image2);
void SwapImages(OldImageInfo *image1,OldImageInfo *image2);
void PrintImages(ImageInfo *images,long Nimages);
//void PrintImageInfo(ImageInfo *image);

// routines using tree

// in image_finder.c

/*void find_images(PosType *y_source,PosType r_source,GridHndl grid
		,int *Nimages,ImageInfo *imageinfo,const int NimageMax,unsigned long *Nimagepoints
		,PosType initial_size,bool splitimages,short edge_refinement
		,bool verbose);
short image_finder(PosType *y_source,PosType r_source,TreeHndl s_tree,TreeHndl i_tree
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
	inline PosType sepSQR(PosType *xx,PosType *yy){
		return pow(xx[0]-yy[0],2) + pow(xx[1]-yy[1],2);
	}
	void double_sort(unsigned long n, PosType *arr, unsigned long *brr);
	void double_sort_points(unsigned long n, PosType *arr, Point *brr);

	void quicksortPoints(Point *pointarray,PosType *arr,unsigned long N);
	void quicksort(unsigned long *particles,PosType *arr,unsigned long N);
	void quickPartition(PosType pivotvalue,unsigned long *pivotindex,unsigned long *particles
		,PosType *arr,unsigned long N);
	void quickPartitionPoints(PosType pivotvalue,unsigned long *pivotindex
		,Point *pointsarray,PosType *arr,unsigned long N);
	int cutbox(const PosType* center,PosType *p1,PosType *p2,float rmax);
	void log_polar_grid(Point *i_points,PosType rmax,PosType rmin,PosType *center,long Ngrid);
	void findarea(ImageInfo *imageinfo);
	int windings2(PosType *x,Point *points,unsigned long Npoints,PosType *area,short image);
	void writeCurves(int m, ImageInfo *critical, int Ncrit, int index);
  PosType cross(const Point *O, const Point *A, const Point *B);
  bool xorder(Point *p1,Point *p2);
  std::vector<Point *> convex_hull(std::vector<Point *> P);
  std::vector<Point *> shrink_rap(std::vector<Point *> P);
  std::vector<Point *> concave_hull(std::vector<Point *> P);



	long IndexFromPosition(PosType *x,long Npixels,PosType range,const PosType *center);
	void PositionFromIndex(unsigned long i,PosType *x,long Npixels,PosType range,PosType const *center);
	long IndexFromPosition(PosType x,long Npixels,PosType range,PosType center);
  PosType TwoDInterpolator(PosType *x,int Npixels,PosType range,PosType *center,PosType *map,bool init=true);
  PosType TwoDInterpolator(PosType *map);

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
                 PosType const *x          /// position of point
                 ,int Npixels       /// Number of pixels in one dimension
                 ,PosType my_range   /// Range of map in same units as x[]
                 ,PosType *my_center /// Center of map in same units as x[]
                 ):
    N(Npixels),range(my_range),map_p(NULL),Ny(Npixels),range_y(my_range)
    {
      center[0] = my_center[0];
      center[1] = my_center[1];

      initparams(x);
    };

    /**
     *  Constructor for case when region is a rectangle and not a square.
     *  Array must be indexed i = ix + iy * Nx 
     */
    Interpolator(
                 PosType const *x          /// position of point
                 ,int my_Nx       /// Number of pixels in x dimension
                 ,PosType my_range_x   /// Range of map in x in same units as x[]
                 ,int my_Ny       /// Number of pixels in y dimension
                 ,PosType my_range_y   /// Range of map in y in same units as x[]
                 ,PosType *my_center /// Center of map in same units as x[]
                 ):
    N(my_Nx),range(my_range_x),map_p(NULL),Ny(my_Ny),range_y(my_range_y)
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
                 ,PosType my_range     /// Range of map in same units as x[]
                 ,PosType *my_center   /// Center of map in same units as x[]
                 ,const T *map        /// One dimensional array of fundamental type
                 ):
    N(Npixels),range(my_range),map_p(map),Ny(Npixels),range_y(my_range)
    {
      center[0] = my_center[0];
      center[1] = my_center[1];
    };
    
    /** 
     Does interpolation of map at point that object was constructed with or last called with.
     Can use any map type that has a [] operator that returns a PosType.
     */
    PosType interpolate(
                       T& map    /// map that supports the [] operator 
                       ){
      if(map.size() != N*Ny){
        ERROR_MESSAGE();
        std::cout << "ERROR: Interpolator:interpolator(T&), wrong size map" << std::endl;
      }
      if(index == -1) return 0;
      
      return (1.-fx)*(1.-fy)*map[index] + fx*(1.-fy)*map[index+1] + fx*fy*map[index+1+N]
          + (1.-fx)*fy*map[index+N];
    };
    /// reinitializes to a new position
    PosType interpolate(
                       PosType *x   /// position of point
                       ,T& map     /// map that supports the [] operator
                       ){
      if(map.size() != N*Ny){
        ERROR_MESSAGE();
        std::cout << "ERROR: Interpolator:interpolator(PosType *,T&), wrong size map" << std::endl;
      }
      initparams(x);
      return interpolate(map);
    }
    
    /**
     Does interpolation of store map at point x. Only for use with the second constructor.
     */
    PosType operator ()(PosType *x){
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
      std::valarray<PosType> map;
      PosType tmp,x[2];
          
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
    PosType range,range_y,center[2];
    const T *map_p;
    PosType fx, fy;
    long index;
    int N,Ny;

    void initparams(PosType const *x){
      long ix,iy;
      // position in pixel coordinates
      fx = ((x[0] - center[0])/range + 0.5)*(N-1);
      fy = ((x[1] - center[1])/range_y + 0.5)*(Ny-1);
      //std::cout << "(  " << fx << " " << fy << "   ";
     
      if (fx < 0. || fx > N-1){index = -1; return;}
      else ix = (unsigned long)(fx);
      if(ix == N-1) ix = N-2;
      
      if (fy < 0. || fy > Ny-1){index = -1; return;}
      else iy = (unsigned long)(fy);
      if(iy == Ny-1) iy = Ny-2;
      
      index = ix + N*iy;

      //std::cout << "  " << ix << " " << iy << " " << index << "   ";
      
      /** bilinear interpolation */
      fx = center[0] + range*( 1.0*(ix)/(N-1) - 0.5 );
      fy = center[1] + range_y*( 1.0*(iy)/(Ny-1) - 0.5 );
      fx=(x[0] - fx)*(N-1)/range;
      fy=(x[1] - fy)*(Ny-1)/range_y;
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

	//inline float isLeft( Point *p0, Point *p1, PosType *x );

	// isLeft(): tests if a point is Left|On|Right of an infinite line.
	// Input:three points P0, P1, and x
	// Return: >0 for x left of the line through P0 and P1
//         =0 for x on the line
//         <0 for x right of the line
	inline float isLeft( Point *p0, Point *p1, PosType *x ){

		return (p1->x[0] - p0->x[0])*(x[1] - p0->x[1])
			- (x[0] - p0->x[0])*(p1->x[1] - p0->x[1]);
	};
	unsigned long prevpower(unsigned long k);

        int windings(PosType *x,Point *points,unsigned long Npoints,PosType *area,short image = 0 );
        int windings(PosType *x,Point **points,unsigned long Npoints,PosType *area,short image = 0 );
	int windings(PosType *x,Kist<Point> * kist,PosType *area,short image = 0);
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
  void ordered_concavehull(Kist<Point> * curve);
  PosType ConvexHullArea(Kist<Point> * curve);
}
bool order_ExteriorBoundary(Point *curve,long Npoints,long *NewNpoints,PosType *area);
PosType findAreaOfCurve(TreeHndl tree,ImageInfo *curve,int NimageMax);
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

//void rayshooterInternal(unsigned long Npoints,Point *i_points);
void in_source(PosType *y_source,ListHndl sourcelist);
bool tree_count_test(TreeHndl tree);
bool testLeafs(TreeHndl tree);

void NeighborsOfNeighbors(ListHndl neighbors,ListHndl wholelist);

#endif
