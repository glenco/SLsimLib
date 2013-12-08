 /*
 * Code Name:     TreeKist.h
 * Programmer:    Ben Metcalf
 * Description:
 * Comments:  This version uses a linked list for the points so that
 the tree can be expended dynamically
 */

#ifndef treekist_declare
#define treekist_declare

//#include "pointlist.h"
#include "point.h"
#include "Kist.h"
#include "image_info.h"

/***** Exported Types *****/

#ifndef criterion_declare
#define criterion_declare
typedef enum{TotalArea,EachImage,Resolution,FillHoles} ExitCriterion;
#endif

/** \brief Tree data structure for points.  This is a rewrite of TreeStruct for 
 *  the moment that uses a Kist instead of a PointList and makes some other improvements.
 *  It is under construction.
 */
struct TreeKist{
public:
	TreeKist(Point *xp,unsigned long Npoints,short my_median_cut = 1,double buffer = 0.0);
	TreeKist(Point *xp,unsigned long npoints
			 ,double boundary_p1[2],double boundary_p2[2]
			 ,double center[2],int Nbucket);

	 ~TreeKist();

	 /// list of points
	 Kist<Point> pointkist;

    //void FindAllBoxNeighbors(Point *point,ListHndl neighbors);
    void FindAllBoxNeighborsKist(Point *point,Kist<Point> * neighbors);
    void PointsWithinEllipKist(const double* center,float rmax,float rmin,float posangle,Kist<Point> * neighborkist);
    double PointsWithinKist(const double* center,float rmax,Kist<Point> * neighborkist,short markpoints);
    void PointsWithinKist_iter(const double* center,float rmin,float rmax,Kist<Point> * neighborkist);
    Point *NearestNeighborKist(const double* center,int Nneighbors,Kist<Point> * neighborkist);

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
    void RebuildTreeFromList();

    /***** State of tree functions *****/
    bool isEmpty();
    bool atTop();
    bool atLeaf(){ return( (current->child1==NULL)*(current->child2==NULL) ); }
    bool offEnd();
    bool CurrentIsSquareBranch();
    bool noChild();

    Kist<Point>::iterator getCurrent(unsigned long *npoints);
    unsigned long getNbranches();
    void printTree();
    void checkTree();

    void FindBoxPoint(const double* ray,Point *point);

    void _FindLeaf(const double* ray,unsigned long Nadd = 0);
  
    // functions for creating subtrees and merging trees
    TreeKist * spawn();
    void merge(TreeKist *treekist1);

private:

  /// private constructor that does nothing but the default actions
  TreeKist(){};

    /// root branch
    Branch *top;
    Branch *current;

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
  void _AddPoint();
  void _BuildTree();

  void _checkTree(unsigned long *count);
  void _freeBranches_iter();
  void _FindBox(const double* ray);

  // Should be obsolete
  Point *NearestNeighbor(const double* center,int Nneighbors,ListHndl neighborlist
  		,short direction);
  void _NearestNeighbor(double* ray,int Nneighbors,Point **neighborpoints,double *rneighbors,short *direction);

};

#endif
