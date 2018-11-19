/**
 * simpleTree.h
 *
 *  Created on: Oct 14, 2011
 *      Author: bmetcalf
 */

#ifndef SIMP_TREE_H_
#define SIMP_TREE_H_

#include "standard.h"
#include "Tree.h"
//#include "lens_halos.h"

/** \brief Box representing a branch in a tree.  It has four children.
 Used in TreeNBStruct which is used in TreeForce.
 */
struct BranchNB{
  
  BranchNB(int Ndim){
    center = new PosType[Ndim*3];
    boundary_p1 = &center[Ndim];
    boundary_p2 = &center[2*Ndim];
  }
  ~BranchNB(){
    delete center;
  }

	/// array of particles in BranchNB
  IndexType *particles;
  IndexType nparticles;
  /// the number of particles that aren't in children
  IndexType big_particle;
  /// Size of largest particle in branch
  PosType maxrsph;
  /// center of mass
  PosType *center;
  PosType mass;
  /// level in tree
  int level;
  unsigned long number;
  /// bottom, left, back corner of box
  PosType *boundary_p1;
  /// top, right, front corner of box
  PosType *boundary_p2;
  BranchNB *child1;
  BranchNB *child2;
  /// father of branch
  BranchNB *prev;
  /// Either child2 of father is branch is child1 and child2 exists or the brother of the father.
  /// Used for iterative tree walk.
  BranchNB *brother;

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
 * TreeNBStruct: Tree structure used for force calculation with particles (i.e. stars, Nbody or halos).
 *
 * The tree also contains pointers to the list of positions, sizes and masses of the particles.
 * Also flags for the number of dimensions the tree is defined in (2 or 3), and if multiple
 * masses and sizes should be used.
 */
template<typename PType>
struct TreeNBStruct{
  BranchNB *top;
  BranchNB *current;
  /// number of branches in tree
  unsigned long Nbranches;
  /// Dimension of tree, 2 or 3.  This will dictate how the force is calculated.
  short Ndimensions;
  /// Array of particle positions
  PType *pp;
};

//typedef struct TreeNBStruct * TreeNBHndl;

/**
 * \brief A C++ class wrapper for the bianary treeNB used in the Nobody
 *  force calculation, but also useful for general purpose searches.
 *
 * Most of the code in TreeNB.c and TreeDriverNB.c is duplicated here as private
 * methods and a few public ones.
 */
template<typename PType>
class TreeSimple {
public:
	TreeSimple(PType *xp,IndexType Npoints,int bucket = 5,int dimensions = 2,bool median = true);
	virtual ~TreeSimple();

	/// Finds the points within a circle around center and puts their index numbers in a list
	template<typename T>
  void PointsWithinCircle(T center[2],float radius,std::list<unsigned long> &neighborkist);
  
	/// Finds the points within an ellipse around center and puts their index numbers in a list
  template<typename T>
	void PointsWithinEllipse(T center[2],float a_max,float a_min,float posangle,std::list<unsigned long> &neighborkist);
  
	/// Finds the nearest N neighbors and puts their index numbers in an array, also returns the distance to the Nth neighbor for calculating smoothing
  //template<typename T>
  //void NearestNeighbors(T *ray,int Nneighbors,float *rsph,IndexType *neighbors) const;
  
  template<typename T>
  T NNDistance(T *ray ,int Nneighbors) const;
  
  class iterator{
  public:
    
    iterator(BranchNB *branch){
      current = top = branch;
    }
    iterator(const iterator &it){
      current = it.current;
      top = it.top;
    }
    
    iterator & operator=(const iterator &it){
      if(&it == this) return *this;
      current = it.current;
      top = it.top;
      
      return *this;
    }

    bool operator==(const iterator &it){
      return (current == it.current)*(top == it.top);
    }

    BranchNB *operator*(){return current;}
    
    /// walk tree below branch last assigned to iterator.  Returnes false if it has completed walking subtree.
    bool walk(bool decend){

      if(decend && current->child1 != NULL){
        current = current->child1;
        return true;
      }
      if(decend && current->child2 != NULL){
        current = current->child2;
        return true;
      }
      
      if(current->brother == top->brother) return false;

      current = current->brother;
      return true;
    }
    
    bool up(){
      if(current == top) return false;
      current = current->prev;
      return true;
    }
    
    bool down(int child){
      if(child == 1){
        if(current->child1 == NULL) return false;
        current = current->child1;
        return true;
      }
      if(child == 2){
        if(current->child2 == NULL) return false;
        current = current->child2;
        return true;
      }
      return false;
    }
    
    bool atleaf(){
      return (current->child1 == NULL)*(current->child2 == NULL);
    }
    bool atTop(){
      return (current == top);
    }

  private:
         BranchNB *current;
         BranchNB *top;
  };
protected:

	int Ndim;//,incell2;
  int incell;
  TreeNBStruct<PType>* tree;
	IndexType *index;
	IndexType Nparticles;
	bool median_cut;
	int Nbucket;
	PosType realray[3];
	PType *xxp;

	TreeNBStruct<PType>* BuildTreeNB(PType *xxp,IndexType Nparticles,IndexType *particles,int Ndimensions,PosType theta);
	void _BuildTreeNB(TreeNBStruct<PType>* tree,IndexType nparticles,IndexType *particles);
  
  template <typename T>
  void _findleaf(T *ray,TreeSimple::iterator &it) const;
  

  template<typename T>
	void _PointsWithin(T *ray,float *rmax,std::list<unsigned long> &neighborkist);
  
  template<typename T>
	void _NearestNeighbors(T *ray,int Nneighbors,unsigned long *neighbors,PosType *rneighbors);

  //template<typename T>
  //void _NearestNeighbors(T *ray,PosType *realray,int Nneighbors,unsigned long *neighbors,PosType *rneighbors
  //                                   ,TreeSimple<PType>::iterator &it,int &bool) const ;

	BranchNB *NewBranchNB(IndexType *particles,IndexType nparticles
			  ,PosType boundary_p1[],PosType boundary_p2[]
			  ,PosType center[],int level,unsigned long branchNBnumber);
	void FreeBranchNB(BranchNB *branchNB);
	TreeNBStruct<PType> * NewTreeNB(IndexType *particles,IndexType nparticles
			 ,PosType boundary_p1[],PosType boundary_p2[],
			     PosType center[],short Ndimensions);
	void freeTreeNB(TreeNBStruct<PType> * tree);
	short emptyTreeNB(TreeNBStruct<PType> * tree);
	void _freeTreeNB(TreeNBStruct<PType> * tree,short child);
	bool isEmptyNB(TreeNBStruct<PType> * tree);
	bool atTopNB(TreeNBStruct<PType> * tree);
	bool noChildNB(TreeNBStruct<PType> * tree);
	bool offEndNB(TreeNBStruct<PType> * tree);
	void getCurrentNB(TreeNBStruct<PType> * tree,IndexType *particles,IndexType *nparticles);
	unsigned long getNbranchesNB(TreeNBStruct<PType> * tree);
	void moveTopNB(TreeNBStruct<PType> * tree);
	void moveUpNB(TreeNBStruct<PType> * tree);
	void moveToChildNB(TreeNBStruct<PType> * tree,int child);
	void insertChildToCurrentNB(TreeNBStruct<PType> * tree, IndexType *particles,IndexType nparticles
				  ,PosType boundary_p1[],PosType boundary_p2[]
				  ,PosType center[],int child);
	void attachChildToCurrentNB(TreeNBStruct<PType> * tree,BranchNB &data,int child);
	bool TreeNBWalkStep(TreeNBStruct<PType> * tree,bool allowDescent);

	inline bool atLeaf(){
		return (tree->current->child1 == NULL)*(tree->current->child2 == NULL);
	}
  template<typename T>
  inline bool inbox(const T* center,PosType *p1,PosType *p2) const{
    return (center[0]>=p1[0])*(center[0]<=p2[0])*(center[1]>=p1[1])*(center[1]<=p2[1]);
  }

	/*TreeNBStruct<PType> * rotate_simulation(PosType **xp,IndexType Nparticles,IndexType *particles
			,PosType **coord,PosType theta,float *rsph,float *mass
			,bool MultiRadius,bool MultiMass);
	TreeNBStruct<PType> * rotate_project(PosType **xp,IndexType Nparticles,IndexType *particles
			,PosType **coord,PosType theta,float *rsph,float *mass
			,bool MultiRadius,bool MultiMass);
	TreeNBStruct<PType> * spread_particles(PosType **xp,IndexType Nparticles,IndexType *particles
			,PosType theta,float *rsph,float *mass,bool MultiRadius,bool MultiMass);
	 void cuttoffscale(TreeNBStruct<PType> * tree,PosType *theta);*/
};


template<typename PType>
TreeSimple<PType>::TreeSimple(PType *xpt,IndexType Npoints,int bucket,int Ndimensions,bool median){
  index = new IndexType[Npoints];
  IndexType ii;
  
  Nbucket = bucket;
  Ndim = Ndimensions;
  median_cut = median;
  Nparticles = Npoints;
  
  xxp = xpt;
  
  for(ii=0;ii<Npoints;++ii) index[ii] = ii;
  
  tree = TreeSimple::BuildTreeNB(xxp,Npoints,index,Ndimensions,0);
  
  return;
}

template<typename PType>
TreeSimple<PType>::~TreeSimple()
{
  freeTreeNB(tree);
  delete[] index;
  return;
}

template<typename PType>
template<typename T>
void TreeSimple<PType>::PointsWithinEllipse(
                                            T *ray     /// center of ellipse
                                            ,float rmax      /// major axis
                                            ,float rmin     /// minor axis
                                            ,float posangle  /// position angle of major axis, smallest angle between the x-axis and the long axis
                                            ,std::list <unsigned long> &neighborlist  /// output neighbor list, will be emptied if it contains anything on entry
){
  PosType x,y,cs,sn;
  
  
  if(rmax < rmin){float tmp = rmax; rmax=rmin; rmin=tmp;}
  
  if(rmax <=0.0 || rmin <= 0.0){ neighborlist.clear(); return; }
  
  // find point within a circle circumscribes the ellipse
  PointsWithinCircle(ray,rmax,neighborlist);
  
  cs = cos(posangle);
  sn = sin(posangle);
  // go through points within the circle and reject points outside the ellipse
  for(  std::list<unsigned long>::iterator it = neighborlist.begin();it != neighborlist.end();){
    x = xxp[*it][0]*cs - xxp[*it][1]*sn;
    y = xxp[*it][0]*sn + xxp[*it][1]*cs;
    if( ( pow(x/rmax,2) + pow(y/rmin,2) ) > 1) it = neighborlist.erase(it);
    else ++it;
  }
  return;
}

template<typename PType>
template<typename T>
void TreeSimple<PType>::PointsWithinCircle(
                                           T *ray     /// center of circle
                                           ,float rmax      /// radius of circle
                                           ,std::list <unsigned long> &neighborlist  /// output neighbor list, will be emptied if it contains anything on entry
){
  
  neighborlist.clear();
  
  realray[0]=ray[0];
  realray[1]=ray[1];
  
  //cout << "ray = " << ray[0] << " " << ray[1] << " rmax = " << rmax << endl;
  moveTopNB(tree);
  if( inbox(ray,tree->current->boundary_p1,tree->current->boundary_p2) == 0 ){
    
    ray[0] = (ray[0] > tree->current->boundary_p1[0]) ? ray[0] : tree->current->boundary_p1[0];
    ray[0] = (ray[0] < tree->current->boundary_p2[0]) ? ray[0] : tree->current->boundary_p2[0];
    
    ray[1] = (ray[1] > tree->current->boundary_p1[1]) ? ray[1] : tree->current->boundary_p1[1];
    ray[1] = (ray[1] < tree->current->boundary_p2[1]) ? ray[1] : tree->current->boundary_p2[1];
    
    //cout << "ray = " << ray[0] << " " << ray[1] << endl;
    //ray[0]=DMAX(ray[0],tree->current->boundary_p1[0]);
    //ray[0]=DMIN(ray[0],tree->current->boundary_p2[0]);
    
    //ray[1]=DMAX(ray[1],tree->current->boundary_p1[1]);
    //ray[1]=DMIN(ray[1],tree->current->boundary_p2[1]);
  }
  incell=1;
  
  _PointsWithin<T>(ray,&rmax,neighborlist);
  
  return;
}
/**
 * Used in PointsWithinKist() to walk tree.*/
template<typename PType>
template<typename T>
void TreeSimple<PType>::_PointsWithin(T *ray,float *rmax,std::list <unsigned long> &neighborlist){
  
  int j,incell2=1;
  unsigned long i;
  PosType radius;
  short pass;
  
  if(incell){  // not found cell yet
    
    if( inbox(ray,tree->current->boundary_p1,tree->current->boundary_p2) ){
      //cout << "   In box" << endl;
      
      // found the box small enough
      if( Utilities::cutbox(ray,tree->current->boundary_p1,tree->current->boundary_p2,*rmax)==1
         || atLeaf() ){
        //cout << "   Found cell" << endl;
        
        // whole box in circle or a leaf with ray in it
        
        incell=0;
        
        ray[0]=realray[0];
        ray[1]=realray[1];
        
        if( atLeaf() ){
          // if leaf calculate the distance to all the points in cell
          for(i=0;i<tree->current->nparticles;++i){
            for(j=0,radius=0.0;j<2;++j) radius+=pow(xxp[tree->current->particles[i]][j]-ray[j],2);
            if( radius < *rmax**rmax ){
              neighborlist.push_back(tree->current->particles[i]);
              //cout << "add point to list" << neighborlist.size() << endl;
              //InsertAfterCurrentKist(neighborlist,tree->current->particles[i]);
            }
          }
        }else{ // put all of points in box into getCurrentKist(imagelist)
          for(i=0;i<tree->current->nparticles;++i){
            neighborlist.push_back(tree->current->particles[i]);
            //InsertAfterCurrentKist(neighborlist,tree->current->particles[i]);
            //cout << "add point to list" << neighborlist.size() << endl;
            
          }
        }
        
      }else{ // keep going down the tree
        
        if(tree->current->child1 !=NULL){
          moveToChildNB(tree,1);
          _PointsWithin(ray,rmax,neighborlist);
          moveUpNB(tree);
          
          incell2=incell;
        }
        
        if(tree->current->child2 !=NULL){
          moveToChildNB(tree,2);
          _PointsWithin(ray,rmax,neighborlist);
          moveUpNB(tree);
        }
        
        // if ray found in second child go back to first to search for neighbors
        if( (incell2==1) && (incell==0) ){
          if(tree->current->child1 !=NULL){
            moveToChildNB(tree,1);
            _PointsWithin(ray,rmax,neighborlist);
            moveUpNB(tree);
          }
        }
      }
    }  // not in the box
    
  }else{    // found cell
    
    pass=Utilities::cutbox(ray,tree->current->boundary_p1,tree->current->boundary_p2,*rmax);
    // does radius cut into the box
    if( pass ){
      //cout << "   Cell found searching other cells" << endl;
      if( atLeaf()  ){  /* leaf case */
        
        for(i=0;i<tree->current->nparticles;++i){
          for(j=0,radius=0.0;j<2;++j) radius+=pow(xxp[tree->current->particles[i]][j]-ray[j],2);
          if( radius < *rmax**rmax ){
            neighborlist.push_back(tree->current->particles[i]);
            //InsertAfterCurrentKist(neighborlist,tree->current->particles[i]);
            //cout << "add point to list" << neighborlist.size() << endl;
            
          }
        }
      }else if(pass==1){ // whole box is inside radius
        
        for(i=0;i<tree->current->nparticles;++i){
          neighborlist.push_back(tree->current->particles[i]);
          //InsertAfterCurrentKist(neighborlist,tree->current->particles[i]);
          //cout << "add point to list" << neighborlist.size() << endl;
          
        }
      }else{
        if(tree->current->child1 !=NULL){
          moveToChildNB(tree,1);
          _PointsWithin(ray,rmax,neighborlist);
          moveUpNB(tree);
        }
        
        if(tree->current->child2 !=NULL){
          moveToChildNB(tree,2);
          _PointsWithin(ray,rmax,neighborlist);
          moveUpNB(tree);
        }
      }
      
    }
  }
  
  return;
}

/**
 *  \brief Nearest neighbor distance
 *
 * Finds the distance to the Nth nearest neighbors in whatever dimensions tree is defined in.
 *  */
template<typename PType>
template<typename T>
T TreeSimple<PType>::NNDistance(
                                T * ray       /// position
                                ,int Nneighbors    /// number of neighbors to be found
) const{
  
  if(tree->top->nparticles <= Nneighbors){
    ERROR_MESSAGE();
    printf("ERROR: in TreeSimple::NearestNeighbors, number of neighbors > total number of particles\n");
    exit(1);
  }
  
  std::vector<PosType> tmp(Ndim);
  
  TreeSimple::iterator iter(tree->top);
  
  // if the real ray is outside of root box find the closest boundary point that is inside
    for(int j=0;j<tree->Ndimensions;++j){
      tmp[j] = (ray[j] > tree->current->boundary_p1[j]) ? ray[j] : tree->current->boundary_p1[j];
      tmp[j] = (ray[j] < tree->current->boundary_p2[j]) ? ray[j] : tree->current->boundary_p2[j];
    }
    
  //std::cout << tmp[0] << " " << tmp[1] << " " << tmp[2] << std::endl;
  
  // find leaf
  _findleaf(tmp.data(),iter);

  //PosType tmp = (10*(tree->top->boundary_p2[0]-tree->top->boundary_p1[0]));
  
  while((*iter)->nparticles < Nneighbors) iter.up();
  auto pass = iter;
  
  std::vector<PosType> r2neighbors((*iter)->nparticles);
  for(int i=0 ; i < (*iter)->nparticles ;++i){
    for(int j=0;j<tree->Ndimensions;++j){
      r2neighbors[i] += pow(tree->pp[(*iter)->particles[i]][j]-ray[j],2);
    }
  }
  
  std::sort(r2neighbors.begin(),r2neighbors.end());
  r2neighbors.resize(Nneighbors);
  
  PosType bestr = sqrt( r2neighbors.back() );
  
  int cutbox = 1;
  for(int dim = 0 ; dim < Ndim ; ++dim){
    cutbox *= ( (ray[dim] - bestr) > (*iter)->boundary_p1[dim] )
    *( (ray[dim] + bestr) < (*iter)->boundary_p2[dim] );
  }
  
  while(!cutbox && !iter.atTop()){
    iter.up();
    cutbox = 1;
    for(int dim = 0 ; dim < Ndim ; ++dim){
      cutbox *= ( (ray[dim] - bestr) > (*iter)->boundary_p1[dim] )
      *( (ray[dim] + bestr) < (*iter)->boundary_p2[dim] );
    }
  }
  
  auto end = (*iter)->brother;
    
  /// walk tree from the top
  bool decend = true;
  do{
    
    if(iter == pass){
      decend =false;
    }else{
      int cutbox = 1;
      for(int dim = 0 ; dim < Ndim ; ++dim){
        cutbox *= ( (ray[dim] + bestr) > (*iter)->boundary_p1[dim] )
        *( (ray[dim] - bestr) < (*iter)->boundary_p2[dim] );
      }

      if(cutbox == 1){
        decend = true;
      
        if(iter.atleaf()){
        
          for(int i=0 ; i<(*iter)->nparticles ;++i){
            PosType r2 = 0;
            for(int j=0;j<tree->Ndimensions;++j){
              r2 += pow(tree->pp[(*iter)->particles[i]][j]-ray[j],2);
            }
          
            if(r2 < bestr*bestr){
              r2neighbors.back() = r2;
              int ii = r2neighbors.size() - 1;
              while(r2neighbors[ii] < r2neighbors[ii-1] && ii > 0){
                std::swap(r2neighbors[ii],r2neighbors[ii-1]);
                --ii;
              }
              bestr = sqrt( r2neighbors.back() );
            }
          }
          decend = false;
        }
      }else{
        decend = false;
      }
    }
  }while(iter.walk(decend) && !(*iter == end) );
  
  return bestr;
}

template<typename PType>
template<typename T>
void TreeSimple<PType>::_findleaf(T *ray,TreeSimple::iterator &it) const {
  
  if(it.atleaf() ) return;
  
  if( inbox(ray,(*it)->child1->boundary_p1,(*it)->child1->boundary_p2) ){
    it.down(1);
    _findleaf(ray,it);
  }else{
    it.down(2);
    _findleaf(ray,it);
  }

  return;
}

/*
template<typename PType>
template<typename T>
void TreeSimple<PType>::_NearestNeighbors(T *ray,PosType *rlray,int Nneighbors
                                          ,unsigned long *neighbors,PosType *rneighbors
                                          ,TreeSimple::iterator &it,bool &notfound) const {
  
  int incellNB2=1;
  IndexType i;
  short j;
  
  if(notfound){  // not found cell yet
    
    if( inbox(ray,(*it)->boundary_p1,(*it)->boundary_p2) ){
      
      // found the box small enough
      if( (*it)->nparticles <= Nneighbors+Nbucket ){
        notfound=0;
        for(j=0;j<tree->Ndimensions;++j) ray[j]=rlray[j];
        
        // calculate the distance to all the particles in cell
        for(i=0;i<(*it)->nparticles;++i){
          for(j=0,rneighbors[i]=0.0;j<tree->Ndimensions;++j){
            rneighbors[i] += pow(tree->pp[(*it)->particles[i]][j]-ray[j],2);
          }
          rneighbors[i]=sqrt( rneighbors[i] );
          assert(rneighbors[i] < 10);
          neighbors[i]=(*it)->particles[i];
        }
        
        Utilities::quicksort(neighbors,rneighbors,(*it)->nparticles);
        
        
      }else{ // keep going down the tree
        
        if(it.down(1)){
          //moveToChildNB(tree,1);
          _NearestNeighbors(ray,rlray,Nneighbors,neighbors,rneighbors,it,notfound);
          //printf("moving up from level %i\n",(*it)->level);
          //moveUpNB(tree);
          it.up();
          incellNB2=notfound;
        }
        
        if(it.down(2)){
          //printf("moving to child2 from level %i\n",(*it)->level);
          //moveToChildNB(tree,2);
          _NearestNeighbors(ray,rlray,Nneighbors,neighbors,rneighbors,it,notfound);
          //printf("moving up from level %i\n",(*it)->level);
          //moveUpNB(tree);
          it.up();
        }
        
        // if ray found in second child go back to first to search for neighbors
        if( (incellNB2==1) && (notfound==0) ){
          if(it.down(1)){
            //moveToChildNB(tree,1);
            _NearestNeighbors(ray,rlray,Nneighbors,neighbors,rneighbors,it,notfound);
            //moveUpNB(tree);
            it.up();
          }
        }
      }
    }
  }else{ // found cell
    // does radius cut into the box
    if( Utilities::cutbox(ray,(*it)->boundary_p1,(*it)->boundary_p2,rneighbors[Nneighbors-1]) ){
      
      if( ((*it)->child1 == NULL)*((*it)->child2 == NULL)){  // leaf case
        
        // combine found neighbors with particles in box and resort
        for(i=Nneighbors;i<((*it)->nparticles+Nneighbors);++i){
          for(j=0,rneighbors[i]=0.0;j<tree->Ndimensions;++j){
            rneighbors[i]+=pow(tree->pp[(*it)->particles[i-Nneighbors]][j]-ray[j],2);
          }
          rneighbors[i]=sqrt( rneighbors[i] );
          assert(rneighbors[i] < 10);
          neighbors[i]=(*it)->particles[i-Nneighbors];
        }
        
        Utilities::quicksort(neighbors,rneighbors,Nneighbors+Nbucket);
        
      }else{
        
        if(it.down(1)){
          //moveToChildNB(tree,1);
          _NearestNeighbors(ray,rlray,Nneighbors,neighbors,rneighbors,it,notfound);
          //moveUpNB(tree);
          it.up();
        }
        
        if(it.down(2)){
          //moveToChildNB(tree,2);
          _NearestNeighbors(ray,rlray,Nneighbors,neighbors,rneighbors,it,notfound);
          //moveUpNB(tree);
          it.up();
        }
      }
    }
  }
  return;
}
*/

template<typename PType>
template<typename T>
void TreeSimple<PType>::_NearestNeighbors(T *ray,int Nneighbors,unsigned long *neighbors,PosType *rneighbors){
  
  int incellNB2=1;
  IndexType i;
  short j;
  
  if(incell){  /* not found cell yet */
    
    if( inbox(ray,tree->current->boundary_p1,tree->current->boundary_p2) ){
      
      /* found the box small enough */
      if( tree->current->nparticles <= Nneighbors+Nbucket ){
        incell=0;
        for(j=0;j<tree->Ndimensions;++j) ray[j]=realray[j];
        
        /* calculate the distance to all the particles in cell */
        for(i=0;i<tree->current->nparticles;++i){
          for(j=0,rneighbors[i]=0.0;j<tree->Ndimensions;++j){
            rneighbors[i] += pow(tree->xp[tree->current->particles[i]][j]-ray[j],2);
          }
          rneighbors[i]=sqrt( rneighbors[i] );
          assert(rneighbors[i] < 10);
          neighbors[i]=tree->current->particles[i];
        }
        
        Utilities::quicksort(neighbors,rneighbors,tree->current->nparticles);
        
      }else{ /* keep going down the tree */
        
        if(tree->current->child1 !=NULL){
          moveToChildNB(tree,1);
          _NearestNeighbors(ray,Nneighbors,neighbors,rneighbors);
          /*printf("moving up from level %i\n",tree->current->level);*/
          moveUpNB(tree);
          
          incellNB2=incell;
        }
        
        if(tree->current->child2 !=NULL){
          /*printf("moving to child2 from level %i\n",tree->current->level);*/
          moveToChildNB(tree,2);
          _NearestNeighbors(ray,Nneighbors,neighbors,rneighbors);
          /*printf("moving up from level %i\n",tree->current->level);*/
          moveUpNB(tree);
        }
        
        /** if ray found in second child go back to first to search for neighbors **/
        if( (incellNB2==1) && (incell==0) ){
          if(tree->current->child1 !=NULL){
            moveToChildNB(tree,1);
            _NearestNeighbors(ray,Nneighbors,neighbors,rneighbors);
            moveUpNB(tree);
          }
        }
      }
    }
  }else{ // found cell
    /* does radius cut into the box */
    if( Utilities::cutbox(ray,tree->current->boundary_p1,tree->current->boundary_p2,rneighbors[Nneighbors-1]) ){
      
      if( (tree->current->child1 == NULL)*(tree->current->child2 == NULL)){  /* leaf case */
        
        /* combine found neighbors with particles in box and resort */
        for(i=Nneighbors;i<(tree->current->nparticles+Nneighbors);++i){
          for(j=0,rneighbors[i]=0.0;j<tree->Ndimensions;++j){
            rneighbors[i]+=pow(tree->xp[tree->current->particles[i-Nneighbors]][j]-ray[j],2);
          }
          rneighbors[i]=sqrt( rneighbors[i] );
          assert(rneighbors[i] < 10);
          neighbors[i]=tree->current->particles[i-Nneighbors];
        }
        
        Utilities::quicksort(neighbors,rneighbors,Nneighbors+Nbucket);
        
      }else{
        
        if(tree->current->child1 !=NULL){
          moveToChildNB(tree,1);
          _NearestNeighbors(ray,Nneighbors,neighbors,rneighbors);
          moveUpNB(tree);
        }
        
        if(tree->current->child2 !=NULL){
          moveToChildNB(tree,2);
          _NearestNeighbors(ray,Nneighbors,neighbors,rneighbors);
          moveUpNB(tree);
        }
      }
    }
  }
  return;
}

template<typename PType>
TreeNBStruct<PType> * TreeSimple<PType>::BuildTreeNB(PType *xp,IndexType Nparticles,IndexType *particles,int Ndims
                                          ,PosType theta){
  TreeNBStruct<PType> * tree;
  IndexType i;
  short j;
  PosType p1[3],p2[3],center[3];
  
  assert(Ndim == Ndims);
  
  for(j=0;j<Ndim;++j){
    p1[j]=xp[0][j];
    p2[j]=xp[0][j];
  }
  
  for(i=0;i<Nparticles;++i){
    for(j=0;j<Ndim;++j){
      if(xp[i][j] < p1[j] ) p1[j]=xp[i][j];
      if(xp[i][j] > p2[j] ) p2[j]=xp[i][j];
    }
  }
  
  for(j=0;j<Ndim;++j) center[j]=(p1[j]+p2[j])/2;
  
  /* Initialize tree root */
  tree=NewTreeNB(particles,Nparticles,p1,p2,center,Ndim);
  
  tree->pp = xp;
  
  /* build the tree */
  _BuildTreeNB(tree,Nparticles,particles);
  
  /* visit every branch to find center of mass and cutoff scale */
  moveTopNB(tree);
  
  return tree;
}


// tree must be created and first branch must be set before start
template<typename PType>
void TreeSimple<PType>::_BuildTreeNB(TreeNBStruct<PType> * tree,IndexType nparticles,IndexType *particles){
  IndexType i,cut,dimension;
  short j;
  BranchNB *cbranch,branch1(Ndim),branch2(Ndim);
  PosType xcut;
  PosType *x;
  
  cbranch=tree->current; /* pointer to current branch */
  
  cbranch->center[0] = (cbranch->boundary_p1[0] + cbranch->boundary_p2[0])/2;
  cbranch->center[1] = (cbranch->boundary_p1[1] + cbranch->boundary_p2[1])/2;
  cbranch->quad[0]=cbranch->quad[1]=cbranch->quad[2]=0;
  
  /* leaf case */
  if(cbranch->nparticles <= Nbucket){
    cbranch->big_particle=0;
    return;
  }
  
  /* initialize boundaries to old boundaries */
  for(j=0;j<Ndim;++j){
    branch1.boundary_p1[j]=cbranch->boundary_p1[j];
    branch1.boundary_p2[j]=cbranch->boundary_p2[j];
    
    branch2.boundary_p1[j]=cbranch->boundary_p1[j];
    branch2.boundary_p2[j]=cbranch->boundary_p2[j];
  }
  cbranch->big_particle=0;
  
  // **** makes sure force does not require nbucket at leaf
  
  /* set dimension to cut box */
  dimension=(cbranch->level % Ndim);
  
  x=(PosType *)malloc((cbranch->nparticles-cbranch->big_particle)*sizeof(PosType));
  for(i=cbranch->big_particle;i<cbranch->nparticles;++i) x[i]=tree->pp[particles[i]][dimension];
  
  if(median_cut){
    Utilities::quicksort(particles,x,cbranch->nparticles-cbranch->big_particle);
    
    cut=(cbranch->nparticles-cbranch->big_particle)/2;
    branch1.boundary_p2[dimension]=x[cut];
    branch2.boundary_p1[dimension]=x[cut];
  }else{
    xcut=(cbranch->boundary_p1[dimension]+cbranch->boundary_p2[dimension])/2;
    branch1.boundary_p2[dimension]=xcut;
    branch2.boundary_p1[dimension]=xcut;
    
    Utilities::quickPartition(xcut,&cut,particles
                              ,x,cbranch->nparticles-cbranch->big_particle);
  }
  
  /* set particle numbers and pointers to particles */
  branch1.prev=cbranch;
  branch1.nparticles=cut;
  branch1.particles=particles+cbranch->big_particle;
  
  branch2.prev=cbranch;
  branch2.nparticles=cbranch->nparticles-cbranch->big_particle - cut;
  if(cut < (cbranch->nparticles-cbranch->big_particle) )
    branch2.particles=&particles[cut+cbranch->big_particle];
  else branch2.particles=NULL;
  
  free(x);
  
  if(branch1.nparticles > 0) attachChildToCurrentNB(tree,branch1,1);
  if(branch2.nparticles > 0) attachChildToCurrentNB(tree,branch2,2);
  
  // work out brothers for children
  if( (cbranch->child1 != NULL) && (cbranch->child2 != NULL) ){
    cbranch->child1->brother = cbranch->child2;
    cbranch->child2->brother = cbranch->brother;
  }
  if( (cbranch->child1 == NULL) && (cbranch->child2 != NULL) )
    cbranch->child2->brother = cbranch->brother;
  if( (cbranch->child1 != NULL) && (cbranch->child2 == NULL) )
    cbranch->child1->brother = cbranch->brother;
  
  
  if( branch1.nparticles > 0 ){
    moveToChildNB(tree,1);
    _BuildTreeNB(tree,branch1.nparticles,branch1.particles);
    moveUpNB(tree);
  }
  
  if(branch2.nparticles > 0 ){
    moveToChildNB(tree,2);
    _BuildTreeNB(tree,branch2.nparticles,branch2.particles);
    moveUpNB(tree);
  }
  
  return;
}

/**
 * visits every branch in tree to calculate
 * two critical lengths rcrit_angle and rcrit_part
 *   and mark largest particle in each node subject
 *   to some conditions
 *
 *   if the sph[] are negative rcrit_part = 0
 */

/***** Constructors/Destructors *****/

/************************************************************************
 * NewBranchNB
 * Returns pointer to new BranchNB struct.  Initializes children pointers to NULL,
 * and sets data field to input.  Private.
 ************************************************************************/
template<typename PType>
BranchNB *TreeSimple<PType>::NewBranchNB(IndexType *particles,IndexType nparticles
                                         ,PosType boundary_p1[],PosType boundary_p2[]
                                         ,PosType center[],int level,unsigned long branchNBnumber){
  
  BranchNB *branchNB = new BranchNB(Ndim);
  int i;
  
  if (!branchNB){
    ERROR_MESSAGE(); fprintf(stderr,"allocation failure in NewBranchNB()\n");
    exit(1);
  }
  branchNB->particles = particles;
  branchNB->nparticles = nparticles;
  
  for(i=0;i<Ndim;++i) branchNB->center[i]=center[i];
  branchNB->level=level;
  
  for(i=0;i<Ndim;++i){
    branchNB->boundary_p1[i]= boundary_p1[i];
    branchNB->boundary_p2[i]= boundary_p2[i];
  }
  
  branchNB->number=branchNBnumber;
  
  branchNB->child1 = NULL;
  branchNB->child2 = NULL;
  branchNB->prev = NULL;
  branchNB->brother = NULL;
  
  return(branchNB);
}

/************************************************************************
 * FreeBranchNB
 * Frees memory pointed to by branchNB.  Private.
 ************************************************************************/
template<typename PType>
void TreeSimple<PType>::FreeBranchNB(BranchNB *branchNB){
  
  assert( branchNB != NULL);
  delete branchNB;
  
  return;
}

/************************************************************************
 * NewTreeNB
 * Returns pointer to new TreeNB struct.  Initializes top, last, and
 * current pointers to NULL.  Sets NbranchNBes field to 0.  Exported.
 ************************************************************************/
template<typename PType>
TreeNBStruct<PType> * TreeSimple<PType>::NewTreeNB(IndexType *particles,IndexType nparticles
                                        ,PosType boundary_p1[],PosType boundary_p2[],
                                        PosType center[],short Ndimensions){
  
  TreeNBStruct<PType> * tree;
  
  tree = new TreeNBStruct<PType>;// *)malloc(sizeof(TreeNBStruct));
  if (!tree){
    ERROR_MESSAGE(); fprintf(stderr,"allocation failure in NewTreeNB()\n");
    exit(1);
  }
  
  tree->top= NewBranchNB(particles,nparticles,boundary_p1,boundary_p2,center,0,0);
  if (!(tree->top)){
    ERROR_MESSAGE(); fprintf(stderr,"allocation failure in NewTreeNB()\n");
    exit(1);
  }
  
  tree->Nbranches = 1;
  tree->current = tree->top;
  tree->Ndimensions=Ndimensions;
  
  return(tree);
}

template<typename PType>
void TreeSimple<PType>::freeTreeNB(TreeNBStruct<PType> * tree){
  /* free treeNB
   *  does not free the particle positions, masses or rsph
   */
  
  if(tree == NULL) return;
  
  emptyTreeNB(tree);
  FreeBranchNB(tree->top);
  delete tree;
  
  return;
}

template<typename PType>
short TreeSimple<PType>::emptyTreeNB(TreeNBStruct<PType> * tree){
  
  moveTopNB(tree);
  _freeTreeNB(tree,0);
  
  assert(tree->Nbranches == 1);
  
  return 1;
}

template<typename PType>
void TreeSimple<PType>::_freeTreeNB(TreeNBStruct<PType> * tree,short child){
  BranchNB *branch;
  
  assert( tree );
  assert( tree->current);
  
  if(tree->current->child1 != NULL){
    moveToChildNB(tree,1);
    _freeTreeNB(tree,1);
  }
  
  if(tree->current->child2 != NULL){
    moveToChildNB(tree,2);
    _freeTreeNB(tree,2);
  }
  
  if( (tree->current->child1 == NULL)*(tree->current->child2 == NULL) ){
    
    if(atTopNB(tree)) return;
    
    branch = tree->current;
    moveUpNB(tree);
    FreeBranchNB(branch);
    
    /*printf("*** removing branch %i number of branches %i\n",branch->number
     ,tree->Nbranches-1);*/
    
    if(child==1) tree->current->child1 = NULL;
    if(child==2) tree->current->child2 = NULL;
    
    --tree->Nbranches;
    
    return;
  }
  
  return;
}

/************************************************************************
 * isEmptyNB
 * Returns "true" if the TreeNB is empty and "false" otherwise.  Exported.
 ************************************************************************/
template<typename PType>
bool TreeSimple<PType>::isEmptyNB(TreeNBStruct<PType> * tree){
  
  assert(tree != NULL);
  return(tree->Nbranches == 0);
}

/************************************************************************
 * atTop
 * Returns "true" if current is the same as top and "false" otherwise.
 * Exported.
 * Pre: !isEmptyNB(tree)
 ************************************************************************/
template<typename PType>
bool TreeSimple<PType>::atTopNB(TreeNBStruct<PType> * tree){
  
  assert(tree != NULL);
  if( isEmptyNB(tree) ){
    ERROR_MESSAGE();
    fprintf(stderr, "TreeNB Error: calling atTop() on empty tree\n");
    exit(1);
  }
  return(tree->current == tree->top);
}

/************************************************************************
 * noChild
 * Returns "true" if the child of the current branchNB does not exist and "false" otherwise.
 * Exported.
 * Pre: !isEmptyNB(tree)
 ************************************************************************/
template<typename PType>
bool TreeSimple<PType>::noChildNB(TreeNBStruct<PType> * tree){
  
  assert(tree != NULL);
  if( isEmptyNB(tree) ){
    
    ERROR_MESSAGE(); fprintf(stderr, "TreeNB Error: calling atTop() on empty tree\n");
    exit(1);
  }
  
  if( (tree->current->child1 == NULL) || (tree->current->child2 == NULL) ) return true;
  return false;
}

/************************************************************************
 * offEndNB
 * Returns "true" if current is off end and "false" otherwise.  Exported.
 ************************************************************************/
template<typename PType>
bool TreeSimple<PType>::offEndNB(TreeNBStruct<PType> * tree){
  
  assert(tree != NULL);
  return(tree->current == NULL);
}

/************************************************************************
 * getCurrentNB
 * Returns the particuls of current.  Exported.
 * Pre: !offEndNB(tree)
 ************************************************************************/
template<typename PType>
void TreeSimple<PType>::getCurrentNB(TreeNBStruct<PType> * tree,IndexType *particles,IndexType *nparticles){
  
  assert(tree != NULL);
  if( offEndNB(tree) ){
    
    ERROR_MESSAGE(); fprintf(stderr, "TreeNB Error: calling getCurrent() when current is off end\n");
    exit(1);
  }
  
  *nparticles=tree->current->nparticles;
  particles=tree->current->particles;
  
  return;
}

/************************************************************************
 * getNbranchesNB
 * Returns the NbranchNBes of tree.  Exported.
 ************************************************************************/
template<typename PType>
unsigned long TreeSimple<PType>::getNbranchesNB(TreeNBStruct<PType> * tree){
  
  assert(tree != NULL);
  return(tree->Nbranches);
}

/***** Manipulation procedures *****/

/************************************************************************
 * moveTopNB
 * Moves current to the front of tree.  Exported.
 * Pre: !isEmptyNB(tree)
 ************************************************************************/
template<typename PType>
void TreeSimple<PType>::moveTopNB(TreeNBStruct<PType> * tree){
  //std::cout << tree << std::endl;
  //std::cout << tree->current << std::endl;
  //std::cout << tree->top << std::endl;
  
  assert(tree != NULL);
  if( isEmptyNB(tree) ){
    
    ERROR_MESSAGE(); fprintf(stderr, "TreeNB Error: calling moveTopNB() on empty tree\n");
    exit(1);
  }
  
  
  tree->current = tree->top;
  
  return;
}

/************************************************************************
 * movePrev
 * Moves current to the branchNB before it in tree.  This can move current
 * off end.  Exported.
 * Pre: !offEndNB(tree)
 ************************************************************************/
template<typename PType>
void TreeSimple<PType>::moveUpNB(TreeNBStruct<PType> * tree){
  
  assert(tree != NULL);
  if( offEndNB(tree) ){
    ERROR_MESSAGE(); fprintf(stderr, "TreeNB Error: call to moveUpNB() when current is off end\n");
    exit(1);
  }
  if( tree->current == tree->top ){
    ERROR_MESSAGE(); fprintf(stderr, "TreeNB Error: call to moveUpNB() tried to move off the top\n");
    exit(1);
  }
  
  tree->current = tree->current->prev;  /* can move off end */
  return;
}

/************************************************************************
 * moveToChildNB
 * Moves current to child branchNB after it in tree.  This can move current off
 * end.  Exported.
 * Pre: !offEndNB(tree)
 ************************************************************************/
template<typename PType>
void TreeSimple<PType>::moveToChildNB(TreeNBStruct<PType> * tree,int child){
  
  assert(tree != NULL);
  if( offEndNB(tree) ){
    
    ERROR_MESSAGE(); fprintf(stderr, "TreeNB Error: calling moveChildren() when current is off end\n");
    exit(1);
  }
  if(child==1){
    if( tree->current->child1 == NULL ){
      ERROR_MESSAGE(); fprintf(stderr, "TreeNB Error: moveToChildNB() typing to move to child1 when it doesn't exist\n");
      exit(1);
    }
    tree->current = tree->current->child1;
  }
  if(child==2){
    if( tree->current->child2 == NULL ){
      ERROR_MESSAGE(); fprintf(stderr, "TreeNB Error: moveToChildNB() typing to move to child2 when it doesn't exist\n");
      exit(1);
    }
    tree->current = tree->current->child2;
  }
  return;
}

/************************************************************************
 * insertAfterCurrent
 * Inserts a new BranchNB after the current branchNB in the tree and sets the
 * data field of the new BranchNB to input.  Exported.
 * Pre: !offEndNB(tree)
 ************************************************************************/
template<typename PType>
void TreeSimple<PType>::insertChildToCurrentNB(TreeNBStruct<PType> * tree, IndexType *particles,IndexType nparticles
                                               ,PosType boundary_p1[],PosType boundary_p2[]
                                               ,PosType center[],int child){
  
  BranchNB *branchNB;
  
  /*printf("attaching child%i  current paricle number %i\n",child,tree->current->nparticles);*/
  
  branchNB = NewBranchNB(particles,nparticles,boundary_p1,boundary_p2,center
                         ,tree->current->level+1,tree->Nbranches);
  
  assert(tree != NULL);
  
  if( offEndNB(tree) ){
    
    ERROR_MESSAGE(); fprintf(stderr, "TreeNB Error: calling insertChildToCurrentNB() when current is off end\n");
    exit(1);
  }
  
  branchNB->prev = tree->current;
  
  if(child==1){
    if(tree->current->child1 != NULL){
      ERROR_MESSAGE(); fprintf(stderr, "TreeNB Error: calling insertChildToCurrentNB() when child1 alread exists\n");
      exit(1);
    }
    tree->current->child1 = branchNB;
  }
  if(child==2){
    if(tree->current->child2 != NULL){
      ERROR_MESSAGE();
      fprintf(stderr, "TreeNB Error: calling insertChildToCurrentNB() when child2 alread exists\n  current level=%i Nbranches=%li\n"
              ,tree->current->level,tree->Nbranches);
      exit(1);
    }
    tree->current->child2 = branchNB;
  }
  
  tree->Nbranches++;
  
  return;
}

/* same as above but takes a branchNB structure */
template<typename PType>
void TreeSimple<PType>::attachChildToCurrentNB(TreeNBStruct<PType> * tree,BranchNB &data,int child){
  
  insertChildToCurrentNB(tree,data.particles,data.nparticles
                         ,data.boundary_p1,data.boundary_p2,data.center,child);
  return;
}

// step for walking tree by iteration instead of recursion
template<typename PType>
bool TreeSimple<PType>::TreeNBWalkStep(TreeNBStruct<PType> * tree,bool allowDescent){
  if(allowDescent && tree->current->child1 != NULL){
    moveToChildNB(tree,1);
    return true;
  }
  if(allowDescent && tree->current->child2 != NULL){
    moveToChildNB(tree,2);
    return true;
  }
  
  if(tree->current->brother != NULL){
    tree->current=tree->current->brother;
    return true;
  }
  return false;
}

#endif /* SIMP_TREE_H_ */
