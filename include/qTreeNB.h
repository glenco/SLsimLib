//
//  qTreeNB.h
//  GLAMER
//
//  Created by Ben Metcalf on 11/11/2018.
//

#ifndef qTreeNB_h
#define qTreeNB_h

#ifndef PI
#define PI  3.141593
#endif

#ifndef treeNBdim
#define treeNBdim 2  // dimension of boxes in tree
#endif

/** type for particle positions and boundaries etc **/

#ifndef IndexType_declare
#define IndexType_declare
typedef unsigned long IndexType;
#endif

/** \brief Box representing a branch in a tree.  It has four children.  Used in QTreeNB which is used in TreeQuad.
 */
struct QBranchNB{
  QBranchNB(){
    static unsigned long n=0;
    child0 = NULL;
    child1 = NULL;
    child2 = NULL;
    child3 = NULL;
    prev = NULL;
    brother = NULL;
    particles = NULL;
    big_particles = NULL;
    nparticles = 0;
    number = n;
    ++n;
  };
  
  ~QBranchNB(){};
  
  /// array of particles in QBranchNB
  IndexType *particles;
  IndexType *big_particles;
  IndexType nparticles;
  /// the number of particles that aren't in children
  IndexType Nbig_particles;
  /// Size of largest particle in branch
  PosType maxrsph;
  /// center of mass
  PosType center[2];
  PosType mass;
  /// level in tree
  int level;
  unsigned long number;
  /// bottom, left, back corner of box
  PosType boundary_p1[2];
  /// top, right, front corner of box
  PosType boundary_p2[2];
  QBranchNB *child0;
  QBranchNB *child1;
  QBranchNB *child2;
  QBranchNB *child3;
  
  /// father of branch
  QBranchNB *prev;
  /// Either child2 of father is branch is child1 and child2 exists or the brother of the father.
  /// Used for iterative tree walk.
  QBranchNB *brother;
  
  /* projected quantities */
  /// quadropole moment of branch
  PosType quad[3];
  /// largest dimension of box
  PosType rmax;
  /// the critical distance below which a branch is opened in the
  PosType rcrit_angle;
  /* force calculation */
  PosType rcrit_part;
  //PosType cm[2]; /* projected center of mass */
  
};

/** \brief
 * QTreeNB: Tree structure used for force calculation with particles (i.e. stars, Nbody or halos).
 *
 * The tree also contains pointers to the list of positions, sizes and masses of the particles.
 * Also flags for the number of dimensions the tree is defined in (2 or 3), and if multiple
 * masses and sizes should be used.
 */
template<typename PType = double*>
struct QTreeNB{
  QTreeNB(PType *xp,IndexType *particles,IndexType nparticles
          ,PosType boundary_p1[],PosType boundary_p2[]);
  ~QTreeNB();
  
  const bool isEmpty();
  const bool atTop();
  const bool noChild();
  const bool offEnd();
  const unsigned long getNbranches();
  const bool atLeaf(QBranchNB *branch){
    return (branch->child0 == NULL)*(branch->child1 == NULL)
    *(branch->child2 == NULL)*(branch->child3 == NULL);
  }
  
  
  void getCurrent(IndexType *particles,IndexType *nparticles);
  void moveTop();
  void moveUp();
  void moveToChild(int child);
  bool WalkStep(bool allowDescent);
  
  short empty();
  void attachChildrenToCurrent(QBranchNB *branch0,QBranchNB *branch1
                               ,QBranchNB *branch2,QBranchNB *branch3);
  
  QBranchNB *top;
  QBranchNB *current;
  /// Array of particle positions
  PType *xxp;
  
private:
  /// number of branches in tree
  unsigned long Nbranches;
  void _freeQTree(short child);
};

/************************************************************************
 * This constructor just creates a root branch.  To construct a full tree
 * an external function must be used.  For example, TreeQuad uses QTreeNB
 ************************************************************************/
template<typename PType>
QTreeNB<PType>::QTreeNB(PType *xp,IndexType *particles,IndexType nparticles
                          ,PosType boundary_p1[],PosType boundary_p2[]): xxp(xp) {
  
  top = new QBranchNB();
  
  top->boundary_p1[0] = boundary_p1[0];
  top->boundary_p1[1] = boundary_p1[1];
  top->boundary_p2[0] = boundary_p2[0];
  top->boundary_p2[1] = boundary_p2[1];
  
  top->nparticles = nparticles;
  top->level = 0;
  top->particles = particles;
  
  top->center[0] = (boundary_p1[0] + boundary_p2[0])/2;
  top->center[1] = (boundary_p1[1] + boundary_p2[1])/2;
  
  Nbranches = 1;
  current = top;
  
}


/**
 *  \brief A iterator class fore TreeStruct that allows for movement through the tree without changing
 *      anything in the tree itself.
 *
 *   This class should be able to preform all of the constant movements within the tree without causing
 *   any change to the tree.
 */
template<typename PType>
class QBiterator{
  
private:
  QBranchNB *current;
  QBranchNB *top;
  
public:
  /// Sets the top or root to the top of "tree".
  QBiterator(QTreeNB<PType> * tree){current = top = tree->top;}
  /// Sets the root to the input branch so that this will be a subtree in branch is not the real root.
  QBiterator(QBranchNB *branch){current = top = branch;}
  
  /// Returns a pointer to the current Branch.
  QBranchNB *operator*(){return current;}
  
  void movetop(){current = top;}
  
  /// Same as up()
  bool operator++(){ return up();}
  
  /// Same as up()
  bool operator++(int){ return up();}
  
  bool up();
  /// Move to child
  bool down(short child);
  const bool atLeaf(){
    return (current->child0 == NULL)*(current->child1 == NULL)
    *(current->child2 == NULL)*(current->child3 == NULL);
  }
  bool TreeWalkStep(bool allowDescent);
};



/// Free treeNB. Does not free the particle positions, masses or sizes
template<typename PType>
QTreeNB<PType>::~QTreeNB(){
  //  void TreeQuad::freeQTreeNB(QTreeNBHndl tree){
  
  empty();
  delete top;
  
  return;
}

template<typename PType>
short QTreeNB<PType>::empty(){
  
  if(Nbranches <= 1) return 1;
  moveTop();
  _freeQTree(0);
  
  assert(Nbranches == 1);
  
  return 1;
}

template<typename PType>
void QTreeNB<PType>::_freeQTree(short child){
  QBranchNB *branch;
  
  assert(current);
  if(current->particles != current->big_particles
     && current->big_particles != NULL) delete[] current->big_particles;
  
  if(current->child0 != NULL){
    moveToChild(0);
    _freeQTree(0);
  }
  
  if(current->child1 != NULL){
    moveToChild(1);
    _freeQTree(1);
  }
  
  if(current->child2 != NULL){
    moveToChild(2);
    _freeQTree(2);
  }
  
  if(current->child3 != NULL){
    moveToChild(3);
    _freeQTree(3);
  }
  
  if( atLeaf(current) ){
    
    if(atTop()) return;
    
    branch = current;
    moveUp();
    delete branch;
    
    /*printf("*** removing branch %i number of branches %i\n",branch->number
     ,Nbranches-1);*/
    
    if(child==0) current->child0 = NULL;
    if(child==1) current->child1 = NULL;
    if(child==2) current->child2 = NULL;
    if(child==3) current->child3 = NULL;
    
    --Nbranches;
    
    return;
  }
  
  return;
}

/************************************************************************
 * isEmpty
 * Returns "true" if the QTreeNB is empty and "false" otherwise.  Exported.
 ************************************************************************/
template<typename PType>
const bool QTreeNB<PType>::isEmpty(){
  
  return(Nbranches == 0);
}

/************************************************************************
 * atTop
 * Returns "true" if current is the same as top and "false" otherwise.
 * Exported.
 * Pre: !isEmptyNB(tree)
 ************************************************************************/
template<typename PType>
const bool QTreeNB<PType>::atTop(){
  
  if( isEmpty() ){
    ERROR_MESSAGE();
    std::cout << "QTreeNB Error: calling atTop() on empty tree" << std::endl;
    exit(1);
  }
  return(current == top);
}

/************************************************************************
 * offEndNB
 * Returns "true" if current is off end and "false" otherwise.  Exported.
 ************************************************************************/
template<typename PType>
const bool QTreeNB<PType>::offEnd(){
  
  return(current == NULL);
}

/************************************************************************
 * getCurrentNB
 * Returns the particles of current.  Exported.
 * Pre: !offEnd()
 ************************************************************************/
template<typename PType>
void QTreeNB<PType>::getCurrent(IndexType *particles,IndexType *nparticles){
  
  if( offEnd() ){
    ERROR_MESSAGE();
    std::cout << "QTreeNB Error: calling getCurrent() when current is off end" << std::endl;
    exit(1);
  }
  
  *nparticles = current->nparticles;
  particles = current->particles;
  
  return;
}

/************************************************************************
 * getNbranches
 * Returns the NbranchNBes of tree.  Exported.
 ************************************************************************/
template<typename PType>
const unsigned long QTreeNB<PType>::getNbranches(){
  
  return(Nbranches);
}

/***** Manipulation procedures *****/

/************************************************************************
 * moveTop
 * Moves current to the front of tree.  Exported.
 * Pre: !isEmpty(tree)
 ************************************************************************/
template<typename PType>
void QTreeNB<PType>::moveTop(){
  //std::cout << tree << std::endl;
  //std::cout << tree->current << std::endl;
  //std::cout << tree->top << std::endl;
  
  if( isEmpty() ){
    ERROR_MESSAGE();
    std::cout << "QTreeNB Error: calling moveTop() on empty tree" << std::endl;
    exit(1);
  }
  
  current = top;
  
  return;
}

/************************************************************************
 * movePrev
 * Moves current to the branchNB before it in tree.  This can move current
 * off end.  Exported.
 * Pre: !offEndNB(tree)
 ************************************************************************/
template<typename PType>
void QTreeNB<PType>::moveUp(){
  
  if( offEnd() ){
    ERROR_MESSAGE();
    std::cout << "QTreeNB Error: call to moveUp() when current is off end" << std::endl;
    exit(1);
  }
  if( current == top ){
    ERROR_MESSAGE();
    std::cout << "QTreeNB Error: call to moveUp() tried to move off the top" << std::endl;
    exit(1);
  }
  
  current = current->prev;  /* can move off end */
  return;
}

/************************************************************************
 * moveToChild
 * Moves current to child branchNB after it in tree.  This can move current off
 * end.  Exported.
 * Pre: !offEnd(tree)
 ************************************************************************/
template<typename PType>
void QTreeNB<PType>::moveToChild(int child){
  
  if( offEnd() ){
    ERROR_MESSAGE(); std::cout << "QTreeNB Error: calling moveChildren() when current is off end" << std::endl;
    exit(1);
  }
  if(child==0){
    if( current->child0 == NULL ){
      ERROR_MESSAGE(); std::cout << "QTreeNB Error: moveToChild() typing to move to child1 when it doesn't exist" << std::endl;
      exit(1);
    }
    current = current->child0;
  }
  if(child==1){
    if( current->child1 == NULL ){
      ERROR_MESSAGE(); std::cout << "QTreeNB Error: moveToChild() typing to move to child1 when it doesn't exist" << std::endl;
      exit(1);
    }
    current = current->child1;
  }
  if(child==2){
    if( current->child2 == NULL ){
      ERROR_MESSAGE(); std::cout << "QTreeNB Error: moveToChild() typing to move to child2 when it doesn't exist" << std::endl;
      exit(1);
    }
    current = current->child2;
  }
  if(child==3){
    if( current->child3 == NULL ){
      ERROR_MESSAGE(); std::cout << "QTreeNB Error: moveToChild() typing to move to child2 when it doesn't exist" << std::endl;
      exit(1);
    }
    current = current->child3;
  }
  return;
}

template<typename PType>
void QTreeNB<PType>::attachChildrenToCurrent(QBranchNB *branch0,QBranchNB *branch1
                                      ,QBranchNB *branch2,QBranchNB *branch3){
  
  Nbranches += 4;
  current->child0 = branch0;
  current->child1 = branch1;
  current->child2 = branch2;
  current->child3 = branch3;
  
  int level = current->level+1;
  branch0->level = level;
  branch1->level = level;
  branch2->level = level;
  branch3->level = level;
  
  // work out brothers for children
  branch0->brother = branch1;
  branch1->brother = branch2;
  branch2->brother = branch3;
  branch3->brother = current->brother;
  
  branch0->prev = current;
  branch1->prev = current;
  branch2->prev = current;
  branch3->prev = current;
  
  return;
}

// step for walking tree by iteration instead of recursion
template<typename PType>
bool QTreeNB<PType>::WalkStep(bool allowDescent){
  if(allowDescent && current->child0 != NULL){
    moveToChild(0);
    return true;
  }
  
  if(current->brother != NULL){
    current=current->brother;
    return true;
  }
  return false;
}

template<typename PType>
bool QBiterator<PType>::up(){
  if( current == top ){
    return false;
  }
  current = current->prev;  /* can move off end */
  return true;
}

/// Move to child
template<typename PType>
bool QBiterator<PType>::down(short child){
  
  if(child==0){
    if( current->child0 == NULL ){
      ERROR_MESSAGE(); std::cout << "QTreeNB Error: moveToChild() typing to move to child1 when it doesn't exist" << std::endl;
      return false;
    }
    current = current->child0;
    return true;
  }
  if(child==1){
    if( current->child1 == NULL ){
      ERROR_MESSAGE(); std::cout << "QTreeNB Error: moveToChild() typing to move to child1 when it doesn't exist" << std::endl;
      return false;
    }
    current = current->child1;
    return true;
  }
  if(child==2){
    if( current->child2 == NULL ){
      ERROR_MESSAGE(); std::cout << "QTreeNB Error: moveToChild() typing to move to child2 when it doesn't exist" << std::endl;
      return false;
    }
    current = current->child2;
    return true;
  }
  if(child==3){
    if( current->child3 == NULL ){
      ERROR_MESSAGE(); std::cout << "QTreeNB Error: moveToChild() typing to move to child2 when it doesn't exist" << std::endl;
      return false;
    }
    current = current->child3;
    return true;
  }
  return true;
}

template<typename PType>
bool QBiterator<PType>::TreeWalkStep(bool allowDescent){
  if(allowDescent && current->child0 != NULL){
    down(0);
    return true;
  }
  
  if(current->brother != NULL){
    current=current->brother;
    return true;
  }
  return false;
}

#endif /* qTreeNB_h */
