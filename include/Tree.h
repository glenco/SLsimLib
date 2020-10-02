/*
 * Code Name:     tree.h
 * Programmer:    Ben Metcalf
 * Description:
 * Comments:  This version uses a linked list for the points so that
 the tree can be expended dynamically
 */

#ifndef treetypes_declare
#define treetypes_declare

#include <mutex>
//#include "pointlist.h"
#include "point.h"
#include "Kist.h"
#include "image_info.h"
#include "utilities_slsim.h"
#include <future>


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
    iterator(iterator &p){
      top = p.top;
      current = p.current;
    }
    
    iterator &operator=(iterator &p){
      if(&p == this){
        return *this;
      }
      
      top = p.top;
      current = p.current;
      
      return *this;
    }
    
    iterator &operator=(Branch *branch){
      current = branch;
      return *this;
    }
    
    /// Returns a pointer to the current Branch.
    Branch *operator*(){return current;}
    
    void movetop(){current = top;}
    
    /// Same as up()
    //bool operator++(){ return up();}
    /// Same as up()
    //bool operator++(int){ return up();}
    
    bool up();
    /// Move to brother if it exists
    bool brother();
    /// Move to child
    bool down(short child);
    bool atLeaf(){ return( (current->child1==NULL)*(current->child2==NULL) ); }
    
    bool TreeWalkStep(bool allowDescent);
    
    bool noChild();
    bool offEnd();
    bool IsSquareBranch();
    bool atTop(){return current==top;}
  };
  
  struct Globals{
    bool incell;
    PosType ray[2];
    PosType realray[2];
    std::vector<Point *> tmp_point;
  };
  
  TreeStruct::iterator begin() const{
    iterator it(top);
    return it;
  };
  
  TreeStruct::iterator end() const{
    iterator it(top->brother);
    return it;
  };
  
  /// root branch
  Branch * getTop(){return top;}
  
  /// list of points
  PointList *pointlist;
  
  
  //void FindAllBoxNeighbors(Point *point,ListHndl neighbors);
  void FindAllBoxNeighborsKist(Point *point,Kist<Point> * neighbors) const;
  void PointsWithinEllipKist(const PosType* center,float rmax,float rmin,float posangle,Kist<Point> * neighborkist) const;
  PosType PointsWithinKist(const PosType* center,PosType rmax,Kist<Point> * neighborkist,short markpoints) const;
  void PointsWithinKist_iter(const PosType* center,float rmin,float rmax,Kist<Point> * neighborkist) const;
  Point *NearestNeighborKist(const PosType* center,int Nneighbors,Kist<Point> * neighborkist) const;
  
  
  bool Test();
  
  Point *RemoveLeafFromTree(TreeStruct::iterator &current,unsigned long *Npoints);
  
  // Higher level builds
  void FillTree(Point *xp,unsigned long Npoints);
  int AddPointsToTree(Point *xpoint,unsigned long Nadd);
  
  short emptyTree();
  void RebuildTreeFromList();
  
  /***** State of tree functions *****/
  bool isEmpty();
  
  //void getCurrent(Point *points,unsigned long *npoints);
  unsigned long getNbranches();
  void printTree(TreeStruct::iterator &current);
  void checkTree();
  
  Point * FindBoxPoint(const PosType* ray) const;
  
  TreeStruct * spawn(TreeStruct::iterator &current);
  
  void _FindLeaf(TreeStruct::iterator &current,const PosType* ray,unsigned long Nadd = 0) const;
  
  static std::mutex mutex;
private:
  
  // Adding and removing to branches of tree
  void insertChildToCurrent(Branch *current,Branch *branch,int child);
  void attachChildrenToCurrent(Branch *current,Branch *child1,Branch *child2);
  
  /// root branch
  Branch *top;
  
  TreeStruct(){};
  
  Point **temp_points;
  
  /// number of barnches in tree */
  unsigned long Nbranches;
  /// number of points allowed in leaves of tree
  int Nbucket;
  short median_cut;
  
  void construct_root(Point *xp,unsigned long npoints
                      ,PosType boundary_p1[2],PosType boundary_p2[2]
                      ,PosType center[2],int Nbucket);
  
  void _freeBranches(TreeStruct::iterator &current, short child);
  void _AddPoint(TreeStruct::iterator &current);
  void _BuildTree(TreeStruct::iterator &current);
  
  void _checkTree(TreeStruct::iterator &current, unsigned long *count);
  void _freeBranches_iter();
  
  // const functions
  void _FindAllBoxNeighborsKist(Branch *leaf,TreeStruct::iterator &treeit,Kist<Point> * neighbors) const;
  void _FindAllBoxNeighborsKist_iter(Branch *leaf,TreeStruct::iterator &treeit,Kist<Point> * neighbors) const;
  void _PointsWithinKist(TreeStruct::iterator &treeit,PosType *rmax,Kist<Point> * neighborkist
                         ,short markpoints,PosType *maxgridsize,TreeStruct::Globals &glabs) const;
  
  void _FindBox(TreeStruct::iterator &current,const PosType* ray) const;
  
  void _NearestNeighbor(TreeStruct::iterator &current,int Nneighbors,Point **neighborpoints,PosType *rneighbors,short *direction,TreeStruct::Globals &glabs) const;
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
inline T MIN(T x,T y){
  return (x < y) ? x : y;
};
template <class T>
inline T MAX(T x,T y){
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
bool BoxIntersectCircle(const PosType* center,PosType radius,PosType *p1,PosType *p2);

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

void findborders2(TreeHndl i_tree,OldImageInfo *imageinfo);
void findborders3(TreeHndl i_tree,OldImageInfo *imageinfo);

// in image_finder_kist.c

void findborders4(TreeHndl i_tree,ImageInfo *imageinfo);

// in find_crit.c
void findborders(TreeHndl i_tree,ImageInfo *imageinfo);

Point *LinkToSourcePoints(Point *i_points,unsigned long Npoints);

///
namespace Utilities{
///Separation squared between two positions in 2 dimensions.
inline PosType sepSQR(PosType *xx,PosType *yy){
  return pow(xx[0]-yy[0],2) + pow(xx[1]-yy[1],2);
}
void double_sort(unsigned long n, PosType *arr, unsigned long *brr);
void double_sort_points(unsigned long n, PosType *arr, Point *brr);

void quicksortPoints(Point *pointarray,PosType *arr,unsigned long N);
void quicksortPoints(Point *pointarray,double (*func)(Point &),unsigned long N);

// sort particles and arr according to order of arr
//void quicksort(unsigned long *particles,PosType *arr,unsigned long N);

template<typename D>
void quicksort(unsigned long *particles,D *arr,unsigned long N){
  
  std::vector<size_t> index(N);
  
  Utilities::sort_indexes(arr,index,N);
  
  Utilities::apply_permutation(particles,index);
  Utilities::apply_permutation(arr,index);
  
  //std::cout << arr[0] << " " << arr[1] << " " << arr[2] << std::endl;
  assert(arr[0] <= arr[N-1]);
  return;
//
//
//  D pivotvalue;
//  unsigned long pivotindex,newpivotindex,i;
//
//  if(N <= 1) return ;
//
//  // pick pivot as the median of the first, last and middle values
//  if ((arr[0] >= arr[N/2] && arr[0] <= arr[N-1])
//      || (arr[0] >= arr[N-1] && arr[0] <= arr[N/2])) pivotindex = 0;
//  else if ((arr[N/2] >= arr[0] && arr[N/2] <= arr[N-1])
//           || (arr[N/2] >= arr[N-1] && arr[N/2] <= arr[0])) pivotindex = N/2;
//  else pivotindex = N-1;
//  pivotvalue=arr[pivotindex];
//
//  // move pivet to end of array
//  std::swap(arr[pivotindex],arr[N-1]);
//  std::swap(particles[pivotindex],particles[N-1]);
//  newpivotindex=0;
//
//  // partition list and array
//  for(i=0;i<N;++i){
//    if(arr[i] <= pivotvalue){
//      std::swap(arr[newpivotindex],arr[i]);
//      std::swap(particles[newpivotindex],particles[i]);
//      ++newpivotindex;
//    }
//  }
//  if(newpivotindex != 0) --newpivotindex;
//
//  quicksort(particles,arr,newpivotindex);
//  quicksort(&particles[newpivotindex+1],&arr[newpivotindex+1],N-newpivotindex-1);
//
//  return ;
}

/*
 void quickPartition(PosType pivotvalue,unsigned long *pivotindex,unsigned long *particles
 ,PosType *arr,unsigned long N);
 
 void quickPartitionPoints(PosType pivotvalue,unsigned long *pivotindex
 ,Point *pointsarray,PosType *arr,unsigned long N);
 void quickPartitionPoints(PosType pivotvalue,unsigned long *pivotindex
 ,Point *pointarray,PosType (*func)(Point &p),unsigned long N);
 */
/*
 * Partitions arr[] and particles[] into those with x <= pivotvalue and those with
 *  x > pivotvalue  pivotindex is left at first array value with x > pivotvalue
 */
template<typename T>
void quickPartition(T pivotvalue,unsigned long *pivotindex
                    ,unsigned long *particles
                    ,T *arr,unsigned long N){
  unsigned long i;
  
  *pivotindex=0;
  
  for(i=0;i<N;++i){
    if(arr[i] <= pivotvalue){
      std::swap(arr[*pivotindex],arr[i]);
      std::swap(particles[*pivotindex],particles[i]);
      ++(*pivotindex);
    }
  }
  
  return ;
}
template<typename T>
void quickPartitionPoints(T pivotvalue,unsigned long *pivotindex
                          ,Point *pointarray,T *arr,unsigned long N){
  unsigned long i;
  
  *pivotindex=0;
  
  for(i=0;i<N;++i){
    if(arr[i] <= pivotvalue){
      std::swap(arr[*pivotindex],arr[i]);
      assert(*pivotindex < N);
      SwapPointsInArray(&pointarray[*pivotindex],&pointarray[i]);
      ++(*pivotindex);
    }
  }
  
  return ;
}

template<typename T>
void quickPartitionPoints(T pivotvalue,unsigned long *pivotindex
                          ,Point *pointarray,T (*func)(Point &p),unsigned long N){
  unsigned long i;
  
  *pivotindex=0;
  
  for(i=0;i<N;++i){
    if(func(pointarray[i]) <= pivotvalue){
      assert(*pivotindex < N);
      SwapPointsInArray(&pointarray[*pivotindex],&pointarray[i]);
      ++(*pivotindex);
    }
  }
  
  return ;
}

//}


void log_polar_grid(Point *i_points,PosType rmax,PosType rmin,PosType *center,long Ngrid);
void findarea(ImageInfo *imageinfo);
void writeCurves(int m, ImageInfo *critical, int Ncrit, int index);
PosType cross(const Point *O, const Point *A, const Point *B);
bool xorder(Point *p1,Point *p2);
bool yorder(Point *p1,Point *p2);
//

PosType cross(const Point *O, const Point *A, const Point *B);
PosType crossD(const double *O, const double *A, const double *B);
PosType crossD(Point_2d &O,Point_2d &A,Point_2d &B);

std::vector<Point *> convex_hull(std::vector<Point *> &P);
std::vector<double *> convex_hull(std::vector<double *> &P);
//void convex_hull(std::vector<Point_2d> &P,std::vector<Point_2d> &hull);

std::vector<Point *> concave_hull(std::vector<Point *> &P,int k,bool test=false);
std::vector<double *> concave_hull(std::vector<double *> &P,int k);


void contour_ellipse(std::vector<Point_2d> &P, Point_2d center, unsigned long Npoints ,std::vector<Point_2d> &C, double *ellipticity, double *ellipse_area) ;
Point_2d contour_center(std::vector<Point_2d> &P, unsigned long Npoints);

long IndexFromPosition(PosType *x,long Npixels,PosType range,const PosType *center);
void PositionFromIndex(unsigned long i,PosType *x,long Npixels,PosType range,PosType const *center);
long IndexFromPosition(PosType x,long Npixels,PosType range,PosType center);
PosType TwoDInterpolator(PosType *x,int Npixels,PosType range,PosType *center,PosType *map,bool init=true);
PosType TwoDInterpolator(PosType *map);

// return 1 (0) if box is (not) within rmax of ray
template<typename T>
int cutbox(const T* center,T *p1,T *p2,float rmax){
  /*  returns:  0 if whole box is outside rmax from ray[]
   *            1 if whole box is inside circle but ray is not in the box
   *            2 if ray[] is inside box
   *            3 if box intersects circle but ray[] is not inside box
   */
  short i,tick=0;
  PosType close[2],rtmp;
  PosType tmp1,tmp2;
  
  // find closest point on box borders to ray[]
  for(i=0;i<2;++i){
    if( center[i] < p1[i] ){
      close[i]=p1[i];
    }else if(center[i] > p2[i]){
      close[i]=p2[i];
    }else{
      close[i]=center[i];
      ++tick;
    }
  }
  
  if(tick==2) return 2;  // ray is inside box
  
  for(i=0,rtmp=0;i<2;++i) rtmp += pow(center[i] - close[i],2);
  
  if(rtmp>rmax*rmax) return 0;  // box is all outside circle
  
  // find farthest point on box border from ray[]
  for(i=0,rtmp=0;i<2;++i) rtmp += ((tmp1 = pow(center[i]-p1[i],2)) > (tmp2=pow(center[i]-p2[i],2))) ? tmp1 : tmp2;
  //for(i=0,rtmp=0;i<2;++i) rtmp += DMAX(pow(ray[i]-p1[i],2),pow(ray[i]-p2[i],2));
  
  if(rtmp<rmax*rmax) return 1;  // box is all inside circle
  
  return 3;  // box intersects circle
}

/// Multi-threaded quicksort.  The maximum number of threads used is 2^lev.  The last parameter should be left out when calling so that it takes the default value
template<int lev>
void quicksortPoints_multithread(Point *pointarray,PosType *arr,unsigned long N,int level=0){
  PosType pivotvalue;
  unsigned long pivotindex,newpivotindex,i;
  
  if(N <= 1) return ;
  
  // pick pivot as the median of the first, last and middle values
  if ((arr[0] >= arr[N/2] && arr[0] <= arr[N-1])
      || (arr[0] >= arr[N-1] && arr[0] <= arr[N/2])) pivotindex = 0;
  else if ((arr[N/2] >= arr[0] && arr[N/2] <= arr[N-1])
           || (arr[N/2] >= arr[N-1] && arr[N/2] <= arr[0])) pivotindex = N/2;
  else pivotindex = N-1;
  pivotvalue=arr[pivotindex];
  
  // move pivot to end of array
  std::swap(arr[pivotindex],arr[N-1]);
  assert(pivotindex < N);
  SwapPointsInArray(&pointarray[pivotindex],&pointarray[N-1]);
  newpivotindex=0;
  
  // partition list and array
  for(i=0;i<N;++i){
    if(arr[i] <= pivotvalue){
      std::swap(arr[newpivotindex],arr[i]);
      SwapPointsInArray(&pointarray[newpivotindex],&pointarray[i]);
      ++newpivotindex;
    }
  }
  if(newpivotindex != 0) --newpivotindex;
  
  if(level < lev && N > 500){
    auto thread1 = std::async(std::launch::async, [&] {
      return quicksortPoints_multithread<lev>(pointarray,arr,newpivotindex,level + 1); });
    
    quicksortPoints_multithread<lev>(&pointarray[newpivotindex+1],&arr[newpivotindex+1],N-newpivotindex-1,level + 1);
    
    //thread1.wait();
  }else{
    quicksortPoints_multithread<lev>(pointarray,arr,newpivotindex,level + 1);
    quicksortPoints_multithread<lev>(&pointarray[newpivotindex+1],&arr[newpivotindex+1],N-newpivotindex-1,level + 1);
  }
  return ;
}
/** \brief Multi-threaded quicksort.  The maximum number of threads used is 2^lev.  The function `func` takes a point and returns the value that is should be sorted by.  The last parameter should be left out when calling so that it takes the default value.
 
 This function is different from quicksort_multithread() in that it uses SwapPointsInArray() instead of std::swap() which is needed to make the image pointers follow the swap.
 */
template<int lev>
void quicksortPoints_multithread(Point *pointarray,double (*func)(Point &),unsigned long N,int level=0){
  
  PosType pivotvalue;
  unsigned long pivotindex,newpivotindex,i;
  
  if(N <= 1) return ;
  
  // pick pivot as the median of the first, last and middle values
  if ((func(pointarray[0]) >= func(pointarray[N/2]) && func(pointarray[0]) <= func(pointarray[N-1]))
      || (func(pointarray[0]) >= func(pointarray[N-1]) && func(pointarray[0]) <= func(pointarray[N/2]))) pivotindex = 0;
  else if ((func(pointarray[N/2]) >= func(pointarray[0]) && func(pointarray[N/2]) <= func(pointarray[N-1]))
           || (func(pointarray[N/2]) >= func(pointarray[N-1]) && func(pointarray[N/2]) <= func(pointarray[0]))) pivotindex = N/2;
  else pivotindex = N-1;
  pivotvalue=func(pointarray[pivotindex]);
  
  // move pivot to end of array
  assert(pivotindex < N);
  SwapPointsInArray(&pointarray[pivotindex],&pointarray[N-1]);
  newpivotindex=0;
  
  // partition list and array
  for(i=0;i<N;++i){
    if(func(pointarray[i]) <= pivotvalue){
      assert(newpivotindex < N);
      SwapPointsInArray(&pointarray[newpivotindex],&pointarray[i]);
      ++newpivotindex;
    }
  }
  if(newpivotindex != 0) --newpivotindex;
  
  if(level < lev && N > 500){
    auto thread1 = std::async(std::launch::async, [&] {
      return quicksortPoints_multithread<lev>(pointarray,func,newpivotindex,level + 1); });
    quicksortPoints_multithread<lev>(&pointarray[newpivotindex+1],func
                                     ,N-newpivotindex-1,level + 1);
    
    //thread1.wait();
  }else{
    quicksortPoints_multithread<lev>(pointarray,func,newpivotindex,level + 1);
    quicksortPoints_multithread<lev>(&pointarray[newpivotindex+1],func
                                     ,N-newpivotindex-1,level + 1);
  }
  return ;
}
/** \brief Multi-threaded quicksort.  The maximum number of threads used is 2^lev.  The function `func` takes a T type and returns the value that is should be sorted by.  The last parameter should be left out when calling so that it takes the default value.  std::swap() is used to swap elements of the array.
 */
template<typename T,int lev>
void quicksort_multithread(T *array,double (*func)(T &),unsigned long N,int level=0){
  
  PosType pivotvalue;
  unsigned long pivotindex,newpivotindex,i;
  
  if(N <= 1) return ;
  
  // pick pivot as the median of the first, last and middle values
  if ((func(array[0]) >= func(array[N/2]) && func(array[0]) <= func(array[N-1]))
      || (func(array[0]) >= func(array[N-1]) && func(array[0]) <= func(array[N/2]))) pivotindex = 0;
  else if ((func(array[N/2]) >= func(array[0]) && func(array[N/2]) <= func(array[N-1]))
           || (func(array[N/2]) >= func(array[N-1]) && func(array[N/2]) <= func(array[0]))) pivotindex = N/2;
  else pivotindex = N-1;
  pivotvalue=func(array[pivotindex]);
  
  // move pivot to end of array
  assert(pivotindex < N);
  SwapPointsInArray(&array[pivotindex],&array[N-1]);
  newpivotindex=0;
  
  // partition list and array
  for(i=0;i<N;++i){
    if(func(array[i]) <= pivotvalue){
      std::swap(array[newpivotindex],array[i]);
      ++newpivotindex;
    }
  }
  if(newpivotindex != 0 ) --newpivotindex;
  
  if(level < lev && N > 500){
    auto thread1 = std::async(std::launch::async, [&] {
      return quicksortPoints_multithread<lev>(array,func
                                              ,newpivotindex,level + 1); });
    quicksort_multithread<lev>(&array[newpivotindex+1],func
                               ,N-newpivotindex-1,level + 1);
    //thread1.wait();
  }else{
    quicksort_multithread<lev>(array,func,newpivotindex,level + 1);
    quicksort_multithread<lev>(&array[newpivotindex+1],func
                               ,N-newpivotindex-1,level + 1);
  }
  return ;
}


/**
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
  N(Npixels),range(my_range),map_p(NULL),range_y(my_range),Ny(Npixels)
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
  range(my_range_x),range_y(my_range_y),map_p(NULL),N(my_Nx),Ny(my_Ny)
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
      tmp = cos( 6.*i*2.*PI/(N-1) );
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
inline float isLeft(Point_2d &p0,Point_2d &p1,Point_2d &x ){
  return (p1[0] - p0[0])*(x[1] - p0[1]) - (x[0] - p0[0])*(p1[1] - p0[1]);
};
unsigned long prevpower(unsigned long k);

int windings(PosType *x,Point *points,unsigned long Npoints,PosType *area,short image = 0 );
int windings(PosType *x,Point **points,unsigned long Npoints,PosType *area,short image = 0 );
int windings(Point_2d &x,std::vector<Point_2d> &point,PosType *area);
int windings(PosType *x,Kist<Point> * kist,PosType *area,short image = 0);
int windings2(PosType *x,Point *points,unsigned long Npoints,PosType *area,short image);
/// returns 1 if it is in the curve and 0 if it is out.  Borders count as in.
int incurve(PosType x[],std::vector<Point *> curve);
/// returns 1 if it is in the curve and 0 if it is out.  Borders count as in.
int incurve(PosType x[],std::vector<Point_2d> curve);

unsigned long order_curve4(Point *curve,long Npoints);
unsigned long order_curve4(Kist<Point> * curve);
unsigned long order_curve5(Kist<Point> * curve);
void ordered_convexhull(Kist<Point> * curve);
void ordered_shrink_wrap(Kist<Point> * curve);
void ordered_concavehull(Kist<Point> * curve);
PosType ConvexHullArea(Kist<Point> * curve);
}
// in curve_routines.c
void nesting_curve(OldImageInfo *curves,int Ncurves);
void split_order_curve(OldImageInfo *curves,int Maxcurves,int *Ncurves);
void split_order_curve2(OldImageInfo *curves,int Maxcurves,int *Ncurves);
void split_order_curve3(OldImageInfo *curves,int Maxcurves,int *Ncurves);
void split_order_curve4(OldImageInfo *curves,int Maxcurves,int *Ncurves);

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
