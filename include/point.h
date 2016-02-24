/*
 * point.h
 *
 *  Created on: Nov 15, 2010
 *      Author: bmetcalf
 *
 *      Defines Point and Branch type.
 *
 *      The Branch type needs to be defined here so that point.leaf can be defined.
 */


#ifndef pointtypes_declare
#define pointtypes_declare

#include <standard.h>
#include "Kist.h"

#ifndef pi
#define pi  3.141593
#endif

#ifndef error_message
#define error_message
#define ERROR_MESSAGE() std::cout << "ERROR: file: " << __FILE__ << " line: " << __LINE__ << std::endl;
#endif

#ifndef boo_declare
#define boo_declare
typedef enum {NO, YES, MAYBE} Boo;
#endif

/**  \brief Class for representing points or vectors in 2 dimensions.  Not that the dereferencing operator is overridden.
 
 */

struct Point_2d{
  Point_2d(){
    x[0]=x[1]=0.0;
  }
  Point_2d(double x1,double x2){
    x[0]=x1;
    x[1]=x2;
  }
  
  ~Point_2d(){};
  
  Point_2d(const Point_2d &p){
    x[0]=p.x[0];
    x[1]=p.x[1];
  }
  Point_2d & operator=(const Point_2d &p){
    if(this == &p) return *this;
    x[0]=p.x[0];
    x[1]=p.x[1];
    return *this;
  }
  
  bool operator==(const Point_2d &p){
    return (x[0] == p.x[0])*(x[1] == p.x[1]);
  }
  
  Point_2d  operator+(const Point_2d &p) const{
    Point_2d tmp;
    tmp.x[0] = x[0] + p.x[0];
    tmp.x[1] = x[1] + p.x[1];
    return tmp;
  }
  Point_2d  operator-(const Point_2d &p) const{
    Point_2d tmp;
    tmp.x[0] = x[0] - p.x[0];
    tmp.x[1] = x[1] - p.x[1];
    return tmp;
  }
  Point_2d & operator+=(const Point_2d &p){
    x[0]+=p.x[0];
    x[1]+=p.x[1];
    return *this;
  }
  Point_2d & operator/=(PosType value){
    x[0]/=value;
    x[1]/=value;
    return *this;
  }
  Point_2d & operator/(PosType value){
    x[0]/=value;
    x[1]/=value;
    return *this;
  }
  Point_2d & operator*=(PosType value){
    x[0]*=value;
    x[1]*=value;
    return *this;
  }
  /// scalar product
  PosType operator*(const Point_2d &p){
    return x[0]*p.x[0] + x[1]*p.x[1];
  }
  /// outer product
  PosType operator^(const Point_2d &p){
    return x[0]*p.x[1] - x[1]*p.x[0];
  }
  
  /// length
  PosType length(){
    return sqrt(x[0]*x[0] + x[1]*x[1]);
  }
  
  void rotate(PosType theta){
    PosType c = cos(theta),s = sin(theta);
    PosType tmp = x[0];
    x[0] = c*tmp + s*x[1];
    x[1] = c*x[1] + s*tmp;
  }
  
  PosType x[2];
  PosType & operator[](size_t i){return x[i];}
};

std::ostream &operator<<(std::ostream &os, Point_2d const &p);


//struct branchstruct;
struct Branch;


/** \brief A point on the source or image plane that contains a position and the lensing quantities */

struct Point: public Point_2d{
    
  Point();
  Point *next;    // pointer to next point in linked list
  Point *prev;
  Point *image;  // pointer to point on image or source plane
  unsigned long id;
  //double x[2];         // the position of the point
  unsigned long head;         // marks beginning of allocated array of points for easy deallocation
  Boo in_image; // marks if point is in image

  //PosType operator[](int i){return x[i];}
  
  // redundant information in image and source points
  KappaType kappa;           // surface density
  KappaType gamma[3];        // shear, third component is the rotation quantity that is only non-zero for multi-plane lensing
  double dt;                 // time delay : double implies permanent precision independently from DOUBLE_PRECISION
  KappaType invmag = 1;          // inverse of magnification
    
  double gridsize;           // the size of the most refined grid the point is in
  float surface_brightness;  // the surface brightness at this points

  Branch *leaf;
  bool flag;

  void Print();
  
  static bool orderX(Point *p1,Point *p2){
    return (p1->x[0] < p2->x[0]);
  }
  static bool orderXrev(Point *p1,Point *p2){
    return (p1->x[0] > p2->x[0]);
  }
  static bool orderY(Point *p1,Point *p2){
    return (p1->x[1] < p2->x[1]);
  }
  static bool orderYrev(Point *p1,Point *p2){
    return (p1->x[1] > p2->x[1]);
  }

  /// returns true if the image is double inverted
  bool inverted(){ return 0 > (1 - kappa + sqrt( fabs((1-kappa)*(1-kappa) - invmag) ) ); }
  
private:
  // make a point uncopyable
  //Point(const Point &p);
  //Point &operator=(Point &p);
  //Point &operator=(const Point &p);
};

std::ostream &operator<<(std::ostream &os, Point const &p);

/// The box representing a branch of a binary tree structure.  Used specifically in TreeStruct for organizing points in the grid.
struct Branch{
	Branch(Point *my_points,unsigned long my_npoints
			  ,double my_boundary_p1[2],double my_boundary_p2[2]
			  ,double my_center[2],int my_level);
	~Branch();

  struct Point *points;        /// pointer to first points in Branch
  //Kist<Point>::iterator pointit;       /// Kist iterator pointing to first point in branch
  
  unsigned long npoints;
  double center[2];
  int level;
  unsigned long number;
  double boundary_p1[2];
  double boundary_p2[2];
  Branch *child1;
  Branch *child2;
  Branch *brother;
  Branch *prev;
  /// Marks point as start of a level of refinement
  bool refined;

  void print();

  PosType area(){return (boundary_p2[0]-boundary_p1[0])*(boundary_p2[1]-boundary_p1[1]);}
    
  std::list<Branch *> neighbors;
private:
  static unsigned long countID;

  // make a Branch uncopyable
  Branch(const Branch &p);
  Branch &operator=(Branch &p);
  Branch &operator=(const Branch &p);

} ;

//typedef struct branchstruct Branch;

/** \brief link list for points, uses the linking pointers within the Point type unlike  Kist */
struct PointList{
  PointList(){
    top=NULL;
    Npoints=0;
    bottom = top;
  }
  ~PointList(){EmptyList();}
  
  class iterator{
  private:
    Point *current;
  public:
    
    iterator(){
      current = NULL;
    }
    iterator(PointList &list){
      current = list.top;
    }
    iterator(iterator &it){
      current = *it;
    }
    iterator(Point *p){
      current = p;
    }
    
    Point *operator*(){return current;}
    PointList::iterator &operator=(PointList::iterator &p){
      if(&p == this) return *this;
      current = p.current;
      return *this;
    }
    
    PointList::iterator &operator=(Point *point){
      current = point;
      return *this;
    }
    
    bool operator++(){
      assert(current);
      if(current->prev == NULL) return false;
      current=current->prev;
      return true;
    }
    
    /// Same as Up()
    bool operator++(int){
      assert(current);
      if(current->prev == NULL) return false;
      current=current->prev;
      return true;
    }
    
    /// Same as Down()
    bool operator--(){
      assert(current);
      if(current->next == NULL) return false;
      current=current->next;
      return true;
     }
    
    /// Same as Down()
    bool operator--(int){
      assert(current);
      if(current->next == NULL) return false;
      current=current->next;
      return true;
    }
    
    void JumpDownList(int jump){
      int i;
      
      if(jump > 0) for(i=0;i<jump;++i) --(*this);
      if(jump < 0) for(i=0;i<abs(jump);++i) ++(*this);
    }

    
    bool operator==(const iterator my_it){
      return (current == my_it.current);
    }
    
    bool operator!=(const iterator my_it){
      return (current != my_it.current);
    }
    
  };
  
  unsigned long size() const {return Npoints;}

  inline bool IsTop(PointList::iterator &it) const{
    return *it == top;
  };
  inline bool IsBottom(PointList::iterator &it) const{
    return *it == bottom;
  };
  
  Point *Top() const {return top;}
  Point *Bottom() const {return bottom;}
  
  void EmptyList();
  void InsertAfterCurrent(iterator &current,double *x,unsigned long id,Point *image);
  void InsertBeforeCurrent(iterator &current,double *x,unsigned long id,Point *image);
  void InsertPointAfterCurrent(iterator &current,Point *);
  void InsertPointBeforeCurrent(iterator &current,Point *);
  
  void MoveCurrentToBottom(iterator &current);
  Point *TakeOutCurrent(iterator &current);

  void InsertListAfterCurrent(iterator &current,PointList *list2);
  void InsertListBeforeCurrent(iterator &current,PointList *list2);
  void MergeLists(PointList* list2);
  void ShiftList(iterator &current);

  void FillList(double **x,unsigned long N
                ,unsigned long idmin);
  void PrintList();

  void setN(unsigned long N){Npoints = N;}
  void setTop(Point *p){top = p;}
  void setBottom(Point *p){bottom = p;}

private:
  Point *top;
  Point *bottom;
  //Point *current;
  unsigned long Npoints;
  
private:
  // make a point uncopyable
  //PointList(const PointList &p);
  PointList &operator=(const PointList &p);

};

typedef struct PointList *ListHndl;
//bool AtTopList(ListHndl list);
//bool AtBottomList(ListHndl list);


//inline void MoveToTopList(ListHndl list);
//inline void MoveToBottomList(ListHndl list);
/*inline void MoveToTopList(ListHndl list){
  list->current=list->top;
};
inline void MoveToBottomList(ListHndl list){
  list->current=list->bottom;
};*/

// ***********************************************************
//   routines for linked list of points
// ************************************************************

Point *NewPoint(double *x,unsigned long id);

//void JumpDownList(ListHndl list,int jump);
//bool MoveDownList(ListHndl list);
//bool MoveUpList(ListHndl list);

void SwapPointsInList(ListHndl list,Point *p1,Point *p2);
Point *sortList(long n, double arr[],ListHndl list,Point *firstpoint);

//void UnionList(ListHndl list1,ListHndl list2);
//bool ArePointsUniqueList(ListHndl list);
//bool IntersectionList(ListHndl list1,ListHndl list2,ListHndl intersection);

/**  \brief Class for representing points or vectors in 3 dimensions.  Not that the dereferencing operator is overridden.
 
 */
struct Point_3d{
  Point_3d(){
    x[0]=x[1]=x[2]=0.0;
  }
  ~Point_3d(){};
  
  Point_3d(const Point_3d &p){
    x[0]=p.x[0];
    x[1]=p.x[1];
    x[2]=p.x[2];
  }
  Point_3d & operator=(const Point_3d &p){
    if(this == &p) return *this;
    x[0]=p.x[0];
    x[1]=p.x[1];
    x[2]=p.x[2];
    return *this;
  }
  Point_3d  operator+(const Point_3d &p) const{
    Point_3d tmp;
    tmp.x[0] = x[0] + p.x[0];
    tmp.x[1] = x[1] + p.x[1];
    tmp.x[2] = x[2] + p.x[2];
    return tmp;
  }
  Point_3d  operator-(const Point_3d &p) const{
    Point_3d tmp;
    tmp.x[0] = x[0] - p.x[0];
    tmp.x[1] = x[1] - p.x[1];
    tmp.x[2] = x[2] - p.x[2];
    return tmp;
  }
  Point_3d & operator+=(const Point_3d &p){
    x[0]+=p.x[0];
    x[1]+=p.x[1];
    x[2]+=p.x[2];
    return *this;
  }
  Point_3d & operator/=(PosType value){
    x[0]/=value;
    x[1]/=value;
    x[2]/=value;
    return *this;
  }
  Point_3d & operator/(PosType value){
    Point_3d tmp;
    
    tmp[0] = x[0]/value;
    tmp[1] = x[1]/value;
    tmp[2] = x[2]/value;
    
    return *this;
  }
  Point_3d & operator*=(PosType value){
    x[0] *=value;
    x[1] *=value;
    x[2] *=value;
    return *this;
  }
  /// scalar product
  PosType operator*(const Point_3d &p){
    return x[0]*p.x[0] + x[1]*p.x[1] + x[2]*p.x[2];
  }
  
  /// length
  PosType length(){
    return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  }
  
  void rotate(PosType theta){
    PosType c = cos(theta),s = sin(theta);
    PosType tmp = x[0];
    x[0] = c*tmp + s*x[1];
    x[1] = c*x[1] + s*tmp;
  }
  
  PosType x[3];
  PosType & operator[](size_t i){return x[i];}
};


std::ostream &operator<<(std::ostream &os, Point_3d const &p);

inline double pointx(Point &p){return p.x[0];}
inline double pointy(Point &p){return p.x[1];}

#endif
