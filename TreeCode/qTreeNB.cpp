/*
 * qTreeNB.cpp
 *
 *  Created on: Jun 6, 2012
 *      Author: bmetcalf
 */

#include "slsimlib.h"

QBranchNB::QBranchNB(){
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
}
QBranchNB::~QBranchNB(){};


/************************************************************************
 * This constructor just creates a root branch.  To construct a full tree
 * an external function must be used.  For example, TreeQuad uses QTreeNB
 ************************************************************************/
QTreeNB::QTreeNB(PosType **xp,IndexType *particles,IndexType nparticles
		 ,PosType boundary_p1[],PosType boundary_p2[]): xp(xp) {

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

/// Free treeNB. Does not free the particle positions, masses or sizes
QTreeNB::~QTreeNB(){
//	void TreeQuad::freeQTreeNB(QTreeNBHndl tree){

	empty();
  	delete top;

	return;
}

short QTreeNB::empty(){

	moveTop();
	_freeQTree(0);

	assert(Nbranches == 1);

	return 1;
}

void QTreeNB::_freeQTree(short child){
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

    if( atLeaf() ){

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
const bool QTreeNB::isEmpty(){

    return(Nbranches == 0);
}

/************************************************************************
 * atTop
 * Returns "true" if current is the same as top and "false" otherwise.
 * Exported.
 * Pre: !isEmptyNB(tree)
 ************************************************************************/
const bool QTreeNB::atTop(){

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
const bool QTreeNB::offEnd(){

    return(current == NULL);
}

/************************************************************************
 * getCurrentNB
 * Returns the particles of current.  Exported.
 * Pre: !offEnd()
 ************************************************************************/
void QTreeNB::getCurrent(IndexType *particles,IndexType *nparticles){

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
const unsigned long QTreeNB::getNbranches(){

    return(Nbranches);
}

/***** Manipulation procedures *****/

/************************************************************************
 * moveTop
 * Moves current to the front of tree.  Exported.
 * Pre: !isEmpty(tree)
 ************************************************************************/
void QTreeNB::moveTop(){
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
void QTreeNB::moveUp(){

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
void QTreeNB::moveToChild(int child){

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


void QTreeNB::attachChildrenToCurrent(QBranchNB *branch0,QBranchNB *branch1
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
bool QTreeNB::WalkStep(bool allowDescent){
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
