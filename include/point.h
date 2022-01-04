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

#ifndef PI
#define PI  3.141593
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
  Point_2d(const PosType *p){
    x[0]=p[0];
    x[1]=p[1];
  }

  
  bool operator==(const Point_2d &p) const{
    return (x[0] == p.x[0])*(x[1] == p.x[1]);
  }
  bool operator!=(const Point_2d &p) const{
    return (x[0] != p.x[0]) || (x[1] != p.x[1]);
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
  Point_2d & operator-=(const Point_2d &p){
    x[0]-=p.x[0];
    x[1]-=p.x[1];
    return *this;
  }
  Point_2d & operator/=(PosType value){
    x[0]/=value;
    x[1]/=value;
    return *this;
  }
  Point_2d operator/(PosType value) const{
    Point_2d tmp;
    tmp[0] = x[0]/value;
    tmp[1] = x[1]/value;
    return tmp;
  }
  Point_2d & operator*=(PosType value){
    x[0]*=value;
    x[1]*=value;
    return *this;
  }
  Point_2d operator*(PosType value) const{
    return Point_2d(x[0]*value,x[1]*value);
  }

  /// scalar product
  PosType operator*(const Point_2d &p) const{
    return x[0]*p.x[0] + x[1]*p.x[1];
  }
  /// outer product
  PosType operator^(const Point_2d &p) const{
    return x[0]*p.x[1] - x[1]*p.x[0];
  }
  
  /// length
  PosType length() const{
    return sqrt(x[0]*x[0] + x[1]*x[1]);
  }

  /// length^2
  PosType length_sqr() const{
    return x[0]*x[0] + x[1]*x[1];
  }
  
  // rotates the point
  void rotate(PosType theta){
    PosType c = cos(theta),s = sin(theta);
    PosType tmp = x[0];
    x[0] = c*tmp - s*x[1];
    x[1] = c*x[1] + s*tmp;
  }
  
  /// returns a copy of the point that it rotated
  Point_2d rotated(PosType theta) const{
    Point_2d p;
    PosType c = cos(theta),s = sin(theta);
    p[0] = c*x[0] - s*x[1];
    p[1] = c*x[1] + s*x[0];
    
    return p;
  }
  
  /// rescale to make a unit length vector
  void unitize(){
    PosType s = length();
    x[0] /= s;
    x[1] /= s;
  }
  
  PosType* data(){return x;}
  
  PosType x[2];
  
  PosType & operator[](size_t i) {return x[i];}
  const PosType & operator[](size_t i) const {return x[i];}
};

std::ostream &operator<<(std::ostream &os, Point_2d const &p);

template <typename T>
struct Matrix2x2{
  
  Matrix2x2(){
  }
  
  Matrix2x2(const Matrix2x2<T> &F){
    a[0] = F.a[0];
    a[1] = F.a[1];
    a[2] = F.a[2];
    a[3] = F.a[3];
  }
  
  template <typename B>
  Matrix2x2<T> operator=(const Matrix2x2<B> &F){
    a[0] = F.a[0];
    a[1] = F.a[1];
    a[2] = F.a[2];
    a[3] = F.a[3];
    
    return *this;
  }
  
 
  Matrix2x2<T> operator*=(T f){
    a[0] *= f;
    a[1] *= f;
    a[2] *= f;
    a[3] *= f;
    
    return *this;
  }
 
  Matrix2x2<T> operator/=(T f){
    a[0] /= f;
    a[1] /= f;
    a[2] /= f;
    a[3] /= f;
    
    return *this;
  }

  Matrix2x2<T> operator*(T f) const{
    Matrix2x2<T> m = *this;
    m *= f;
    return m;
  }

  Point_2d operator*(const Point_2d &v) const{
    Point_2d v2;
    v2[0] = a[0]*v[0] + a[1]*v[1];
    v2[1] = a[2]*v[0] + a[3]*v[1];
    return v2;
  }


  Matrix2x2<T> operator/(T f) const{
    Matrix2x2<T> m = *this;
    m /= f;
    
    return m;
  }
  
  template <typename B>
  Matrix2x2<T> operator*(const Matrix2x2<B> &F) const{
    Matrix2x2<T> m;
    
    m.a[0] = a[0] * F.a[0] + a[1] * F.a[2];
    m.a[1] = a[0] * F.a[1] + a[1] * F.a[3];
    m.a[2] = a[2] * F.a[0] + a[3] * F.a[2];
    m.a[3] = a[2] * F.a[1] + a[3] * F.a[3];
    
    return m;
  }
 
  template <typename B>
  Matrix2x2<T> operator+=(const Matrix2x2<B> &F){
    a[0] += F.a[0];
    a[1] += F.a[1];
    a[2] += F.a[2];
    a[3] += F.a[3];
    
    return *this;
  }
  
  template <typename B>
  Matrix2x2<T> operator+(const Matrix2x2<B> &F) const{
    Matrix2x2<T> m = *this;
    m += F;

    return m;
  }
  template <typename B>
  Matrix2x2<T> operator-=(const Matrix2x2<B> &F){
    a[0] -= F.a[0];
    a[1] -= F.a[1];
    a[2] -= F.a[2];
    a[3] -= F.a[3];
    
    return *this;
  }
  
  template <typename B>
  Matrix2x2<T> operator-(const Matrix2x2<B> &F)const{
    Matrix2x2<T> m = *this;
    m -= F;

    return m;
  }

  ///  column , row
  T & operator()(int i,int j){
    return a[ i + 2*j ];
  }
  
  T operator()(int i,int j) const{
    return a[ i + 2*j ];
  }
 
  
  T & operator[](int i){
    return a[i];
  }
 
  T det() const{
    return a[0]*a[3] - a[1]*a[2];
  }

  Matrix2x2<T> inverse() const{
    Matrix2x2<T> m;
    
    m.a[0] = a[3];
    m.a[3] = a[0];
    m.a[1] = -a[1];
    m.a[2] = -a[2];
 
    m /= det();
    return m;
  }
  
  void invert(){
    
    std::swap(a[0],a[3]);
    a[1] *= -1;
    a[2] *= -1;
    
    *this /= det();
  }

  T a[4];
  
  // make the matrix the identity
  void setToI(){
    a[0] = a[3] = 1;
    a[1] = a[2] = 0;
  }
  
  static Matrix2x2<T> I(){
    Matrix2x2<T> m;
    m.a[0] = m.a[3] = 1;
    m.a[1] = m.a[2] = 0;
    return m;
  }
  
  static Matrix2x2<T> sig1(){
    Matrix2x2<T> m;
    m.a[0] = -1;
    m.a[3] = 1;
    m.a[1] = m.a[2] = 0;
    return m;
  }
  static Matrix2x2<T> sig2(){
    Matrix2x2<T> m;
    m.a[0] = m.a[3] = 0;
    m.a[1] = m.a[2] = -1;
    return m;
  }

  static Matrix2x2<T> sig3(){
    Matrix2x2<T> m;
    m.a[0] = m.a[3] = 0;
    m.a[1] = 1;
    m.a[2] = -1;
    return m;
  }

  T kappa() const {return 1-(a[0] + a[3])/2; }
  // defined according to GLAMER II convention
  T gamma1() const{ return (a[0] - a[3])/2;}
  T gamma2() const{ return (a[1] + a[2])/2;}
  T gamma3() const{ return (a[2] - a[1])/2;}
  
  void set(T kappa,T gamma[3]){
    a[0] = 1 - kappa + gamma[0];
    a[3] = 1 - kappa - gamma[0];
    
    a[1] = gamma[1] - gamma[3];
    a[2] = gamma[1] + gamma[3];
  }

  void gamma(T *g) const{
    g[0] = gamma1();
    g[1] = gamma2();
    g[2] = gamma3();
  }
  

  void print() const {
    std::cout << std::endl;
    std::cout << a[0] << "  " << a[1] << std::endl;
    std::cout << a[2] << "  " << a[3] << std::endl;
  }
};


//struct branchstruct;
struct Branch;


/** \brief A point on the source or image plane that contains a position and the lensing quantities */

struct Point: public Point_2d{
    
  Point();
  Point(const Point_2d &p);
  Point(PosType x,PosType y);
  Point *next;    // pointer to next point in linked list
  Point *prev;
  Point *image;  // pointer to point on image or source plane
  unsigned long id;
  unsigned long head;         // marks beginning of allocated array of points for easy deallocation
  Boo in_image; // marks if point is in image

  Point operator=(const Point_2d &p){
    Point_2d::operator=(p);
    
    return *this;
  }
 
  double dt;                 // time delay : double implies permanent precision independently from DOUBLE_PRECISION
 
  Matrix2x2<KappaType> A;
  
  KappaType invmag() const{
    return A.det();
  }
  KappaType gamma1() const{
    return A.gamma1();
  }
  KappaType gamma2() const{
    return A.gamma2();
  }
  KappaType gamma3() const{
    return A.gamma3();
  }
  KappaType kappa() const{
    return A.kappa();
  }

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

  /// returns true if the image is double inverted,  At very low magnification this can fail.
  bool inverted(){ return 0 > (A.a[0] + sqrt( fabs(A.a[0]*A.a[0] - invmag()) ) ); }
};



/** \brief A point that automatically has an image point.
 
 This does not produce a stack of source plan points that are contigious in memory.
 */
struct LinkedPoint : public Point
{
  LinkedPoint(){
    image = &im;
    im.image = this;
  }
private:
  Point im;
};

/** \brief Simple representaion of a light path giving position on the image and source planes and lensing quantities.
*/
struct RAY{
  RAY(){
    dt = 0.0;
    A = Matrix2x2<KappaType>::I();
  };
  
  RAY(const Point &p){
    x[0] = p.x[0];
    x[1] = p.x[1];
    y[0] = p.image->x[0];
    y[1] = p.image->x[1];

    dt = p.dt;
        
    A = p.A;
  };
  RAY(const LinkedPoint &p){
    x[0] = p.x[0];
    x[1] = p.x[1];
    y[0] = p.image->x[0];
    y[1] = p.image->x[1];

    dt = p.dt;
        
    A = p.A;
  };
  RAY(const RAY &p){
    x = p.x;
    y = p.y;

    dt = p.dt;
    
    A = p.A;
  };

  RAY & operator=(const Point &p){
    x[0] = p.x[0];
    x[1] = p.x[1];
    y[0] = p.image->x[0];
    y[1] = p.image->x[1];

    dt = p.dt;
        
    A = p.A;
    return *this;
  };
  
  RAY & operator=(const RAY &p){
    x = p.x;
    y = p.y;
    
    dt = p.dt;
        
    A = p.A;
    return *this;
  };
  
  ~RAY(){};
  
  /// image position
  Point_2d x;
  /// source position
  Point_2d y;
  
  Matrix2x2<KappaType> A;
  
  KappaType invmag() const{
    return A.det();
  }
  KappaType gamma1() const{
    return A.gamma1();
  }
  KappaType gamma2() const{
    return A.gamma2();
  }
  KappaType gamma3() const{
    return A.gamma3();
  }
  KappaType kappa() const{
    return A.kappa();
  }
  
  /// time-delay
  KappaType dt;
  
  /// inverse of the magnification
  KappaType invmag(){return A.det();}
  
  /// deflection angle
  Point_2d alpha(){return x - y;}
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
  
  struct iterator{

    Point *current;
    
    iterator():current(NULL){ }
    
/*    iterator(PointList &list){
      current = list.top;
    }
    iterator(iterator &it){
      current = *it;
    }
 */
    iterator(Point *p){
      current = p;
    }
    
 /*
    PointList::iterator &operator=(Point *point){
      current = point;
      return *this;
    }
*/
    Point *operator*(){return current;}
    
    /*iterator &operator=(iterator &p){
      if(&p == this) return *this;
      current = p.current;
      return *this;
    }*/

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
  
  /* Many changes need to be made to implement this correctly
  PointList::iterator begin() const{
    PointList::iterator it;
    it.current = bottom;
    return it;
  }

  PointList::iterator end() const{
    PointList::iterator it;
    it.current = top->prev;
    return it;
  }
  */
  
  
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
  unsigned long Npoints;
  
  // make a point uncopyable
  //PointList(const PointList &p);
  PointList &operator=(const PointList &p);

};

typedef struct PointList *ListHndl;

// ***********************************************************
//   routines for linked list of points
// ************************************************************

Point *NewPoint(double *x,unsigned long id);
void SwapPointsInList(ListHndl list,Point *p1,Point *p2);
Point *sortList(long n, double arr[],ListHndl list,Point *firstpoint);

/**  \brief Class for representing points or vectors in 3 dimensions.  Not that the dereferencing operator is overridden.
 
 */
template <typename T = PosType>
struct Point_3d{
  Point_3d(){
    x[0]=x[1]=x[2]=0.0;
  }
  Point_3d(T xx,T yy,T zz){
    x[0]=xx;
    x[1]=yy;
    x[2]=zz;
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
  Point_3d & operator-=(const Point_3d &p){
    x[0]-=p.x[0];
    x[1]-=p.x[1];
    x[2]-=p.x[2];
    return *this;
  }
  Point_3d & operator/=(T value){
    x[0]/=value;
    x[1]/=value;
    x[2]/=value;
    return *this;
  }
  Point_3d operator/(T value) const{
    Point_3d tmp;
    tmp[0] = x[0]/value;
    tmp[1] = x[1]/value;
    tmp[2] = x[2]/value;
    
    return tmp;
  }
  Point_3d & operator*=(T value){
    x[0] *=value;
    x[1] *=value;
    x[2] *=value;
    return *this;
  }
  Point_3d operator*(PosType value) const{
    Point_3d tmp;
    tmp[0] = x[0]*value;
    tmp[1] = x[1]*value;
    tmp[2] = x[2]*value;
    
    return tmp;
  }
  
  /// scalar product
  T operator*(const Point_3d &p) const {
    return x[0]*p.x[0] + x[1]*p.x[1] + x[2]*p.x[2];
  }

  /// outer product
  Point_3d<T> operator^(const Point_3d<T> &p) const{
    Point_3d<T> v;
    v[0] = x[1]*p[2] - x[2]*p[1];
    v[1] = x[2]*p[0] - x[0]*p[2];
    v[2] = x[0]*p[1] - x[1]*p[0];
    return v;
  }

  /// length
  T length() const{
    return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  }

  T length_sqr() const{
    return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
  }
  
  /// a rotation theta around the z-axis followed by a rotation phi around the y-axis
  void rotate(T theta,T phi){
    T c = cos(theta),s = sin(theta);
    T tmp = c*x[0] - s*x[1];
    x[1] = c*x[1] + s*x[0];
    
    c = cos(phi);
    s = sin(phi);;
    x[0]  = c*tmp - s*x[2];
    x[2] = c*x[2] + s*tmp;
  }
  
  /// rescale to make a unit length vector
  void unitize(){
    T s = length();
    x[0] /= s;
    x[1] /= s;
    x[2] /= s;
  }
  
  /// returns the unit vector in the direction of the right handed spherical coordinate Phi
  Point_3d<T> unitPhi(){
    double s = sqrt(x[0]*x[0] + x[1]*x[1]);
    return Point_3d<T>(-x[1],x[0],0) / s;
  }

  /// returns the unit vector in the direction of the right handed spherical coordinate Theta, Theta= pi/2 is the north
  /// pole and Theta=0 is the equater (z=0).
  Point_3d<T> unitTheta(){
    Point_3d<T> v(-x[0]*x[2],-x[1]*x[2] ,x[0]*x[0] + x[1]*x[1]);
    v.unitize();
    return v;
  }

  T* data(){return x;}
  
  T x[3];
  T & operator[](size_t i){return x[i];}
  const T &  operator[](size_t i) const {return x[i];}
};

template <typename T>
std::ostream &operator<<(std::ostream &os, Point_3d<T> const &p) {
  return os << p.x[0] << " " << p.x[1] << " " << p.x[2];
}

inline double pointx(Point &p){return p.x[0];}
inline double pointy(Point &p){return p.x[1];}


#endif
