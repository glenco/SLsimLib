/**
 * simpleTreeVec.h
 *
 *  Created on: Oct 14, 2011
 *      Author: bmetcalf
 */

#ifndef SIMP_TREE_V_
#define SIMP_TREE_V_

#include <memory>
#include "standard.h"
#include "Tree.h"

    /** \brief
     * A tree for doing quick searches in multidimensional space.  A pointer to an array of
     *  objects type T is provided.  If they are in a vector, vector.data() can be put in the constructor.
     *  The order of objects in the array is not changed.  If Mypos is not specified the method D * T::x() must exist.
     *  It is this coordinate that is used for the object's position.  Another definition for the position of and object can made by specifying Mypos.
     */
template <typename T,typename D = PosType>
class TreeSimpleVec {
public:
  
  TreeSimpleVec(
                T *xpt                 /// array of object, xpt[i].x[0...dimensions-1] must exist or Mypos must be specified
                ,IndexType Npoints      /// number of object in array
                ,int bucket = 5         /// number of points in leaves of tree
                ,int dimensions = 2     /// dimension of space
                ,bool median = true     /// whether to use a median cut or a space cut in splitting branches
                //  ,PosType *(*Mypos)(T&) = defaultposition  /// function that takes
                ,D *(*Mypos)(T&) = [](T& in){return in.x;}  /// function that takes the object T and returns a pointer to its position, default is t.x[]
  ):realray(dimensions),Nbranches(0),position(Mypos),Nbucket(bucket),Ndimensions(dimensions),
  median_cut(median),Nparticles(Npoints),points(xpt)
  {
    index.resize(Npoints);
    for(IndexType ii=0;ii<Npoints;++ii) index[ii] = ii;
    
    BuildTree();
  }
  
  virtual ~TreeSimpleVec()
  {
    //delete top_ptr;
    //freeTree();
    //delete[] index;
    //assert(Nbranches == 0);
    return;
  };
  
  /// \brief Finds the points within a circle around center and puts their index numbers in a list
  void PointsWithinCircle(D center[2],float radius,std::list<unsigned long> &neighborkist);
  /// \brief Finds the points within an ellipse around center and puts their index numbers in a list
  void PointsWithinEllipse(D center[2],float a_max,float a_min,float posangle,std::list<unsigned long> &neighborkist);
  /// \brief Finds the nearest N neighbors and puts their index numbers in an array, also returns the distance to the Nth neighbor for calculating smoothing
  void NearestNeighbors(D *ray,int Nneighbors,std::vector<D> &radii
                        ,std::vector<IndexType> &neighbors);
  void NearestNeighbor(D *ray,D &radius,IndexType &neighbor);
  
  /// pop this point out of the tree
  void pop(IndexType index_of_point){
    
    //print();
    
    size_t n = 0;
    while(index[n] != index_of_point && n < Nparticles) ++n;
    if(n == Nparticles) return; // particle not in tree
    
    std::shared_ptr<BranchV> leaf;
    
    //bool leaf_found = false;
    moveTop();
    do{
      long p = current->branch_index - index.data();
      if(atLeaf() && n >= p && n < p + current->nparticles ){
        leaf=current;
      }
      if(p > n){
        current->branch_index -= 1;
      }
    }while(TreeWalkStep(true));
    
    assert(leaf->child1==nullptr);
    assert(leaf->child2==nullptr);
    
    //BranchV *tmp = leaf;
    while(leaf != top_ptr){
      if(leaf->nparticles > 0) leaf->nparticles -= 1;
      leaf = leaf->prev_ptr;
    }
    leaf->nparticles -= 1;
    
    // remove from index
    for(size_t i = n ; i < Nparticles - 1 ; ++i){
      index[i] = index[i+1];
    }
    --Nparticles;
    index.pop_back();
    
    //print();
    
  }
  
  /// returns the index numbers of the remaining points which are less than the original if pop() was used.
  std::vector<IndexType> get_index() const{
    return index;
  }
  
  void print(){
    moveTop();
    do{
      for(int i=0 ; i<current->level ; ++i ) std::cout << "    ";
      std::cout << current->level << " " << current->nparticles << " : ";
      for(int i=0 ; i<current->nparticles ; ++i){
        std::cout << *(current->branch_index+i) << " " ;
      }
      std::cout << std::endl;
    }while(TreeWalkStep(true));
    std::cout << std::endl;
  }
  
  /** \brief Box representing a branch in a tree.  It has four children.  Used in TreeNBStruct which is used in TreeForce.
   */
  
  struct BranchV{
    
    int Ndim;
    /// array of branch_index in BranchV
    IndexType *branch_index;
    IndexType nparticles;
    /// level in tree
    int level;
    unsigned long number;
    /// bottom, left, back corner of box
    D *boundary_p1 = nullptr;
    /// top, right, front corner of box
    D *boundary_p2 = nullptr;
    std::vector<D> bank;
    
    std::shared_ptr<BranchV> child1_ptr;
    std::shared_ptr<BranchV> child2_ptr;
    /// father of branch
    std::shared_ptr<BranchV> prev_ptr;
    /// Either child2 of father is branch is child1 and child2 exists or the brother of the father.
    /// Used for iterative tree walk.
    std::shared_ptr<BranchV> brother_ptr;
    
    
    BranchV(int my_Ndim,IndexType *my_branch_index,IndexType my_nparticles
            ,D my_boundary_p1[],D my_boundary_p2[]
            ,int my_level,unsigned long my_BranchVnumber)
    :Ndim(my_Ndim),branch_index(my_branch_index),nparticles(my_nparticles)
    ,level(my_level),number(my_BranchVnumber)
    {
      
      bank.resize(Ndim*2);
      boundary_p1 = bank.data();
      boundary_p2 = bank.data() + Ndim;
      
      for(size_t i=0;i<Ndim;++i){
        boundary_p1[i]= my_boundary_p1[i];
        boundary_p2[i]= my_boundary_p2[i];
      }
      
      child1_ptr = nullptr;
      child2_ptr = nullptr;
      prev_ptr = nullptr;
      brother_ptr = nullptr;
    };
    BranchV(BranchV &branch){
      
      Ndim = branch.Ndim;
      branch_index = branch.branch_index;
      nparticles = branch.nparticles;
      level = branch.level;
      number = branch.number;
      
      bank.resize(Ndim*2);
      boundary_p1 = bank.data();
      boundary_p2 = bank.data() + Ndim;
      
      for(size_t i=0;i<Ndim;++i){
        boundary_p1[i]= branch.boundary_p1[i];
        boundary_p2[i]= branch.boundary_p2[i];
      }
      
      child1_ptr = nullptr;
      child2_ptr = nullptr;
      prev_ptr = nullptr;
      brother_ptr = nullptr;
    };
    
    ~BranchV(){
      //delete child1_ptr;
      //delete child2_ptr;
    };
  };
  
protected:
  
  int incell,incell2;
  std::vector<IndexType> index;
  IndexType Nparticles;
  bool median_cut;
  int Nbucket;
  Point_nd<D> realray;
  T *points;
  D *(*position)(T&);
  
  std::vector<D> workspace;
  
  std::shared_ptr<BranchV> top_ptr;
  std::shared_ptr<BranchV> current;
  /// number of branches in tree
  unsigned long Nbranches;
  /// Dimension of tree, 2 or 3.  This will dictate how the force is calculated.
  short Ndimensions;
  
  void BuildTree();
  
  void _BuildTree(IndexType nparticles,IndexType *tmp_index);
  
  void _PointsWithin(D *ray,float *rmax,std::list<unsigned long> &neighborkist);
  void _NearestNeighbors(D *ray,int Nneighbors,unsigned long *neighbors,D *rneighbors);
  
  void freeTree();
  
  //short clearTree();
  //void _freeTree(short child);
  bool isEmpty();
  bool atTop();
  bool noChild();
  bool offEnd();
  void getCurrent(IndexType *branch_index,IndexType *nparticles);
  unsigned long getNbranches();
  void moveTop();
  void moveUp();
  void moveToChild(int child);
  void attachChildToCurrent(IndexType *branch_index,IndexType nparticles
                            ,D boundary_p1[],D boundary_p2[]
                            ,int child);
  void attachChildToCurrent(BranchV &data,int child);
  bool TreeWalkStep(bool allowDescent);
  
  inline bool atLeaf(){
    return (current->child1_ptr == nullptr)*(current->child2_ptr == nullptr);
  }
  inline bool inbox(const D* center,D *p1,D *p2){
    int tt=1;
    for(int i=0;i<Ndimensions;++i) tt *= (center[i]>=p1[i]);
    
    return tt;
    //        return (center[0]>=p1[0])*(center[0]<=p2[0])*(center[1]>=p1[1])*(center[1]<=p2[1]);
  }
  
  
  int cutbox(const D *center,D *p1,D *p2,float rmax);
};
/*
template <class T>
TreeSimpleVec<T>::TreeSimpleVec<T>(T *xpt,IndexType Npoints,int bucket,int dimensions,bool median){
        index = new IndexType[Npoints];
        IndexType ii;
        
        Nbucket = bucket;
        Ndimensions = dimensions;
        median_cut = median;
        Nparticles = Npoints;
        
        points = xpt;
        
        for(ii=0;ii<Npoints;++ii) index[ii] = ii;
        
        BuildTree();
};*/

/************************************************************************
 * NewBranchV
 * Returns pointer to new BranchV struct.  Initializes children pointers to nullptr,
 * and sets data field to input.  Private.
 ************************************************************************/



template <typename T,typename D>
void TreeSimpleVec<T,D>::PointsWithinEllipse(
                                        D *ray     /// center of ellipse
                                        ,float rmax      /// major axis
                                        ,float rmin     /// minor axis
                                        ,float posangle  /// position angle of major axis, smallest angle between the x-axis and the long axis
                                        ,std::list <unsigned long> &neighborlist  /// output neighbor list, will be emptied if it contains anything on entry
                                        ){
	D x,y,cs,sn;
    
    
	if(rmax < rmin){float tmp = rmax; rmax=rmin; rmin=tmp;}
    
	if(rmax <=0.0 || rmin <= 0.0){ neighborlist.clear(); return; }
    
	// find point within a circle circumscribes the ellipse
	PointsWithinCircle(ray,rmax,neighborlist);
    
	cs = cos(posangle);
	sn = sin(posangle);
	// go through points within the circle and reject points outside the ellipse
	for(  std::list<unsigned long>::iterator it = neighborlist.begin();it != neighborlist.end();){
		x = points[*it][0]*cs - points[*it][1]*sn;
		y = points[*it][0]*sn + points[*it][1]*cs;
		if( ( pow(x/rmax,2) + pow(y/rmin,2) ) > 1) it = neighborlist.erase(it);
		else ++it;
	}
	return;
}

template <typename T,typename D>
void TreeSimpleVec<T,D>::PointsWithinCircle(
                                       D *ray     /// center of circle
                                       ,float rmax      /// radius of circle
                                       ,std::list <unsigned long> &neighborlist  /// output neighbor list, will be emptied if it contains anything on entry
                                       ){
    
    neighborlist.clear();
    
    //realray[0]=ray[0];
    //realray[1]=ray[1];
    for(int j=0 ; j<Ndimensions ; ++j) realray[j]=ray[j];
    
    //cout << "ray = " << ray[0] << " " << ray[1] << " rmax = " << rmax << endl;
    moveTop();
    if( inbox(ray,current->boundary_p1,current->boundary_p2) == 0 ){
        
      for(int i=0 ; i<Ndimensions ; ++i){
        ray[i] = (ray[i] > current->boundary_p1[i]) ? ray[i] : current->boundary_p1[i];
        ray[i] = (ray[i] < current->boundary_p2[i]) ? ray[i] : current->boundary_p2[i];
      }
      
    }
    incell=1;
    
    _PointsWithin(ray,&rmax,neighborlist);
    
    //for(auto i : neighborlist){
    //  assert( rmax*rmax > (position(points[i])[0]-ray[0])*(position(points[i])[0]-ray[0]) + (position(points[i])[1]-ray[1])*(position(points[i])[1]-ray[1]) );
    //}
    return;
}
/**
 * Used in PointsWithinKist() to walk tree.*/
template <typename T,typename D>
void TreeSimpleVec<T,D>::_PointsWithin(D *ray,float *rmax,std::list <unsigned long> &neighborlist){
    
    int j,incell2=1;
    unsigned long i;
    D radius;
    short pass;
    
    if(incell){  // not found cell yet
        
        if( inbox(ray,current->boundary_p1,current->boundary_p2) ){
            //cout << "   In box" << endl;
            
            // found the box small enough
            //if( cutbox(ray,current->boundary_p1,current->boundary_p2,*rmax)==1
            //   || atLeaf() ){
            if( atLeaf() ){
               //cout << "   Found cell" << endl;
                
                // whole box in circle or a leaf with ray in it
                
                incell=0;
                
                //ray[0]=realray[0];
                //ray[1]=realray[1];
                for(int j=0;j<Ndimensions;++j) ray[j]=realray[j];
                if( atLeaf() ){
                    // if leaf calculate the distance to all the points in cell
                    for(i=0;i<current->nparticles;++i){
                        for(j=0,radius=0.0;j<2;++j) radius+=pow(position(points[current->branch_index[i]])[j]-ray[j],2);
                        if( radius < *rmax**rmax ){
                            neighborlist.push_back(current->branch_index[i]);
                            //cout << "add point to list" << neighborlist.size() << endl;
                            //InsertAfterCurrentKist(neighborlist,current->branch_index[i]);
                        }
                    }
                }else{ // put all of points in box into getCurrentKist(imagelist)
                    for(i=0;i<current->nparticles;++i){
                        neighborlist.push_back(current->branch_index[i]);
                        //InsertAfterCurrentKist(neighborlist,current->branch_index[i]);
                        //cout << "add point to list" << neighborlist.size() << endl;
                        
                    }
                }
                
            }else{ // keep going down the tree
                
                if(current->child1_ptr != nullptr){
                    moveToChild(1);
                    _PointsWithin(ray,rmax,neighborlist);
                    moveUp();
                    
                    incell2=incell;
                }
                
                if(current->child2_ptr != nullptr){
                    moveToChild(2);
                    _PointsWithin(ray,rmax,neighborlist);
                    moveUp();
                }
                
                // if ray found in second child go back to first to search for neighbors
                if( (incell2==1) && (incell==0) ){
                    if(current->child1_ptr != nullptr){
                        moveToChild(1);
                        _PointsWithin(ray,rmax,neighborlist);
                        moveUp();
                    }
                }
            }
        }  // not in the box
        
    }else{    // found cell
        
        pass = cutbox(ray,current->boundary_p1,current->boundary_p2,*rmax);
        // does radius cut into the box
        if( pass ){
            //cout << "   Cell found searching other cells" << endl;
            if( atLeaf()  ){  /* leaf case */
                
                for(i=0;i<current->nparticles;++i){
                    for(j=0,radius=0.0;j<2;++j) radius+=pow(position(points[current->branch_index[i]])[j]-ray[j],2);
                    if( radius < *rmax**rmax ){
                        neighborlist.push_back(current->branch_index[i]);
                    }
                }
            //}else if(pass==1){ // whole box is inside radius
                
            //    for(i=0;i<current->nparticles;++i){
            //        neighborlist.push_back(current->branch_index[i]); ?????
            //    }
            }else{
                if(current->child1_ptr != nullptr){
                    moveToChild(1);
                    _PointsWithin(ray,rmax,neighborlist);
                    moveUp();
                }
                
                if(current->child2_ptr != nullptr){
                    moveToChild(2);
                    _PointsWithin(ray,rmax,neighborlist);
                    moveUp();
                }
            }
            
        }
    }
    
    return;
}

/**
 *  \brief finds the nearest neighbors in whatever dimensions tree is defined in
 *  */
template <typename T,typename D>
void TreeSimpleVec<T,D>::NearestNeighbors(
                                     D *ray       /// position
                                     ,int Nneighbors    /// number of neighbors to be found
                                     ,std::vector<D> &radii     /// distance to neighbors
                                     ,std::vector<IndexType> &neighbors  /// list of the indexes of the neighbors
                                     ){
    IndexType i;
    short j;

    
    if(top_ptr->nparticles <= Nneighbors){
  
      std::vector<D> tmp(Nparticles);
      std::vector<size_t> sort_index(Nparticles);
      for(int i = 0 ; i< Nparticles ; ++i){
        tmp[i] = pow(position(points[i])[0]-ray[0],2) + pow(position(points[i])[1]-ray[1],2);
        sort_index[i]=i;
      }
      
      std::sort(sort_index.begin(),sort_index.end(),[&tmp](size_t i,size_t j){return tmp[i] < tmp[j];});
      
      radii.resize(Nparticles);
      neighbors.resize(Nparticles);

      for(int i = 0 ; i< Nparticles ; ++i){
        neighbors[i] = index[sort_index[i]];
        radii[i] = tmp[sort_index[i]];
      }
      
      return;
    }
                                       
    radii.resize(Nneighbors+Nbucket);
    neighbors.resize(Nneighbors+Nbucket);

    
    /* initalize distance to neighbors to a large number */
    for(i=0;i<Nbucket+Nneighbors;++i){
        radii[i] = (10*(top_ptr->boundary_p2[0]-top_ptr->boundary_p1[0]));
        neighbors[i] = 0;
    }
    
    for(int j=0;j<Ndimensions;++j) realray[j]=ray[j];
    
    moveTop();
    if( inbox(ray,current->boundary_p1,current->boundary_p2) == 0 ){
        //ERROR_MESSAGE();
        
        for(j=0;j<Ndimensions;++j){
          ray[j] = (ray[j] > current->boundary_p1[j]) ? ray[j] : current->boundary_p1[j];
          ray[j] = (ray[j] < current->boundary_p2[j]) ? ray[j] : current->boundary_p2[j];
        }
    }
    incell = 1;
    
    _NearestNeighbors(ray,Nneighbors,neighbors.data(),radii.data());
    
    neighbors.resize(Nneighbors);
    radii.resize(Nneighbors);
                                       
    return;
}

/**
 *  \brief finds the nearest neighbors in whatever dimensions tree is defined in
 *  */
template <typename T,typename D>
void TreeSimpleVec<T,D>::NearestNeighbor(
                                     D *ray       /// position
                                     ,D &radius     /// distance of furthest neighbor found from ray[]
                                     ,IndexType &neighbor  /// list of the indexes of the neighbors
                                     ){
    IndexType i;
    //static int count=0,oldNneighbors=-1;
    //short j;
    
    D rneighbors[1+Nbucket];
    IndexType neighbors[1+Nbucket];
    
    if(top_ptr->nparticles < 1){
        ERROR_MESSAGE();
        printf("ERROR: in NearestNeighbors, number of neighbors > total number of particles\n");
      throw std::runtime_error("Asked for too many neighbors");
    }
    
    /* initalize distance to neighbors to a large number */
    for(i=0;i<Nbucket+1;++i){
        rneighbors[i] = (10*(top_ptr->boundary_p2[0]-top_ptr->boundary_p1[0]));
        neighbors[i] = 0;
    }
    
    for(int j=0;j<Ndimensions;++j) realray.x[j]=ray[j];
    
    moveTop();
    if( inbox(ray,current->boundary_p1,current->boundary_p2) == 0 ){
        //ERROR_MESSAGE();
        
        for(int j=0;j<Ndimensions;++j){
          ray[j] = (ray[j] > current->boundary_p1[j]) ? ray[j] : current->boundary_p1[j];
          ray[j] = (ray[j] < current->boundary_p2[j]) ? ray[j] : current->boundary_p2[j];
        }
    }
    incell = 1;
    
    _NearestNeighbors(ray,1,neighbors,rneighbors);
    
    neighbor = neighbors[0];
    radius = rneighbors[0];
    
    return;
}

template <typename T,typename D>
void TreeSimpleVec<T,D>::_NearestNeighbors(D *ray,int Nneighbors,unsigned long *neighbors,D *rneighbors){
    
    int incellNB2=1;
    IndexType i;
    short j;
    
    if(incell){  /* not found cell yet */
        
        if( inbox(ray,current->boundary_p1,current->boundary_p2) ){
            
            /* found the box small enough */
            if( current->nparticles <= Nneighbors+Nbucket ){
                incell=0;
                for(j=0;j<Ndimensions;++j) ray[j]=realray[j];
                
                /* calculate the distance to all the particles in cell */
                for(i=0;i<current->nparticles;++i){
                    for(j=0,rneighbors[i]=0.0;j<Ndimensions;++j){
                        rneighbors[i] += pow(position(points[current->branch_index[i]])[j]-ray[j],2);
                    }
                    rneighbors[i]=sqrt( rneighbors[i] );
                    //assert(rneighbors[i] < 10);
                    neighbors[i]=current->branch_index[i];
                }
                
                Utilities::quicksort(neighbors,rneighbors,current->nparticles);
                
            }else{ /* keep going down the tree */
                
                if(current->child1_ptr !=nullptr){
                    moveToChild(1);
                    _NearestNeighbors(ray,Nneighbors,neighbors,rneighbors);
                    /*printf("moving up from level %i\n",current->level);*/
                    moveUp();
                    
                    incellNB2=incell;
                }
                
                if(current->child2_ptr !=nullptr){
                    /*printf("moving to child2 from level %i\n",current->level);*/
                    moveToChild(2);
                    _NearestNeighbors(ray,Nneighbors,neighbors,rneighbors);
                    /*printf("moving up from level %i\n",current->level);*/
                    moveUp();
                }
                
                /** if ray found in second child go back to first to search for neighbors **/
                if( (incellNB2==1) && (incell==0) ){
                    if(current->child1_ptr !=nullptr){
                        moveToChild(1);
                        _NearestNeighbors(ray,Nneighbors,neighbors,rneighbors);
                        moveUp();
                    }
                }
            }
        }
    }else{ // found cell
		// does radius cut into the box
        if( cutbox(ray,current->boundary_p1,current->boundary_p2,rneighbors[Nneighbors-1]) ){
            
            if( (current->child1_ptr == nullptr)*(current->child2_ptr == nullptr)){  /* leaf case */
                
                /* combine found neighbors with particles in box and resort */
                for(i=Nneighbors;i<(current->nparticles+Nneighbors);++i){
                    for(j=0,rneighbors[i]=0.0;j<Ndimensions;++j){
                        rneighbors[i]+=pow(position(points[current->branch_index[i-Nneighbors]])[j]-ray[j],2);
                    }
                    rneighbors[i]=sqrt( rneighbors[i] );
                    //assert(rneighbors[i] < 10);
                    neighbors[i]=current->branch_index[i-Nneighbors];
                }
                
                Utilities::quicksort(neighbors,rneighbors,Nneighbors+Nbucket);
                
            }else{
                
                if(current->child1_ptr !=nullptr){
                    moveToChild(1);
                    _NearestNeighbors(ray,Nneighbors,neighbors,rneighbors);
                    moveUp();
                }
                
                if(current->child2_ptr !=nullptr){
                    moveToChild(2);
                    _NearestNeighbors(ray,Nneighbors,neighbors,rneighbors);
                    moveUp();
                }
            }
        }
    }
    return;
}

template <typename T,typename D>
void TreeSimpleVec<T,D>::BuildTree(){
    IndexType i;
    short j;
  std::vector<D> p1(Ndimensions),p2(Ndimensions);
  
  if(Nparticles > 0){
    for(j=0;j<Ndimensions;++j){
        p1[j] = position(points[0])[j];
        p2[j] = position(points[0])[j];
    }
    
    for(i=0;i<Nparticles;++i){
        for(j=0;j<Ndimensions;++j){
            if(position(points[i])[j] < p1[j] ) p1[j] = position(points[i])[j];
            if(position(points[i])[j] > p2[j] ) p2[j] = position(points[i])[j];
        }
    }
    
    /* Initialize tree root */
    top_ptr.reset(new BranchV(Ndimensions,index.data(),Nparticles,p1.data(),p2.data(),0,0));
    current = top_ptr;
    ++Nbranches;
    if (!(top_ptr)){
      ERROR_MESSAGE(); fprintf(stderr,"allocation failure in NewTree()\n");
      exit(1);
    }
  
  
    /* build the tree */
    workspace.resize(Nparticles);
    _BuildTree(Nparticles,index.data());
    workspace.clear();
    workspace.shrink_to_fit();
  }else{
    
    for(j=0;j<Ndimensions;++j){
      p1[j] *= 0;
      p2[j] *= 0;
    }

    top_ptr.reset(new BranchV(Ndimensions,index.data(),0,p1.data(),p2.data(),0,0));
    current = top_ptr;
    ++Nbranches;
    if (!(top_ptr)){
      ERROR_MESSAGE(); fprintf(stderr,"allocation failure in NewTree()\n");
      exit(1);
    }
  }
    /* visit every branch to find center of mass and cutoff scale */
    moveTop();
}


// tree must be created and first branch must be set before start
template <typename T,typename D>
void TreeSimpleVec<T,D>::_BuildTree(IndexType nparticles,IndexType *tmp_index){
    IndexType i,cut,dimension;
    BranchV branch1(*current),branch2(*current);
    D xcut;
  std::shared_ptr<BranchV> cbranch;
    
    cbranch=current; // pointer to current branch
  
    // leaf case
    if(cbranch->nparticles <= Nbucket){
        return;
    }
       //cbranch->big_particle=0;
    
    // **** makes sure force does not require nbucket at leaf
    
    // set dimension to cut box
    dimension=(cbranch->level % Ndimensions);
  
   // D *x = new D[cbranch->nparticles];
    D *x = workspace.data();
    for(i=0;i<cbranch->nparticles;++i) x[i] = position(points[tmp_index[i]])[dimension];
    
    if(median_cut){
        Utilities::quicksort(tmp_index,x,cbranch->nparticles);
        
        cut=(cbranch->nparticles)/2;
        branch1.boundary_p2[dimension]=x[cut];
        branch2.boundary_p1[dimension]=x[cut];
    }else{
        xcut=(cbranch->boundary_p1[dimension]+cbranch->boundary_p2[dimension])/2;
        branch1.boundary_p2[dimension]=xcut;
        branch2.boundary_p1[dimension]=xcut;
        
        Utilities::quickPartition(xcut,&cut,tmp_index
                                  ,x,cbranch->nparticles);
    }
    
    // set particle numbers and pointers to my_index
    branch1.prev_ptr = cbranch;
    branch1.nparticles=cut;
  
    branch2.prev_ptr = cbranch;
    branch2.nparticles=cbranch->nparticles - cut;
    if(cut < (cbranch->nparticles) )
        branch2.branch_index = &tmp_index[cut];
    else branch2.branch_index=nullptr;
  
    if(branch1.nparticles > 0) attachChildToCurrent(branch1,1);
    if(branch2.nparticles > 0) attachChildToCurrent(branch2,2);
    
    // work out brothers for children
    if( (cbranch->child1_ptr != nullptr) && (cbranch->child2_ptr != nullptr) ){
        cbranch->child1_ptr->brother_ptr = cbranch->child2_ptr;
        cbranch->child2_ptr->brother_ptr = cbranch->brother_ptr;
    }
    if( (cbranch->child1_ptr == nullptr) && (cbranch->child2_ptr != nullptr) )
        cbranch->child2_ptr->brother_ptr = cbranch->brother_ptr;
    if( (cbranch->child1_ptr != nullptr) && (cbranch->child2_ptr == nullptr) )
        cbranch->child1_ptr->brother_ptr = cbranch->brother_ptr;
    
    
    if( branch1.nparticles > 0 ){
        moveToChild(1);
        _BuildTree(branch1.nparticles,branch1.branch_index);
        moveUp();
    }
    
    if(branch2.nparticles > 0 ){
        moveToChild(2);
        _BuildTree(branch2.nparticles,branch2.branch_index);
        moveUp();
    }
    
    return;
}

/// Free all the tree branches
template <typename T,typename D>
void TreeSimpleVec<T,D>::freeTree(){
    
	//clearTree();
  delete top_ptr;
  Nbranches=0;

	return;
}

//template <typename T,typename D>
//short TreeSimpleVec<T,D>::clearTree(){
//
//	moveTop();
//	_freeTree(0);
//
//	assert(Nbranches == 1);
//
//	return 1;
//}

//template <typename T,typename D>
//void TreeSimpleVec<T,D>::_freeTree(short child){
//	BranchV *branch;
//
//	assert( current);
//
//	if(current->child1 != nullptr){
//		moveToChild(1);
//		_freeTree(1);
//	}
//
//    if(current->child2 != nullptr){
//        moveToChild(2);
//        _freeTree(2);
//    }
//
//    if( (current->child1 == nullptr)*(current->child2 == nullptr) ){
//
//    	if(atTop()) return;
//
//    	branch = current;
//    	moveUp();
//      delete branch;
//      --Nbranches;
//
//    	/*printf("*** removing branch %i number of branches %i\n",branch->number
//         ,Nbranches-1);*/
//
//      if(child==1) current->child1 = nullptr;
//    	if(child==2) current->child2 = nullptr;
//
//    	return;
//    }
//
//    return;
//}

/************************************************************************
 * isEmpty
 * Returns "true" if the Tree is empty and "false" otherwise.  Exported.
 ************************************************************************/
template <typename T,typename D>
bool TreeSimpleVec<T,D>::isEmpty(){
    return(Nbranches == 0);
}

/************************************************************************
 * atTop
 * Returns "true" if current is the same as top and "false" otherwise.
 * Exported.
 * Pre: !isEmpty()
 ************************************************************************/
template <typename T,typename D>
bool TreeSimpleVec<T,D>::atTop(){
    
    if( isEmpty() ){
    	ERROR_MESSAGE();
    	fprintf(stderr, "Tree Error: calling atTop() on empty tree\n");
    	exit(1);
    }
    return(current == top_ptr);
}

/************************************************************************
 * noChild
 * Returns "true" if the child of the current BranchV does not exist and "false" otherwise.
 * Exported.
 * Pre: !isEmpty()
 ************************************************************************/
template <typename T,typename D>
bool TreeSimpleVec<T,D>::noChild(){
    
    if( isEmpty() ){
        
        ERROR_MESSAGE(); fprintf(stderr, "Tree Error: calling atTop() on empty tree\n");
        exit(1);
    }
    
    if( (current->child1 == nullptr) || (current->child2 == nullptr) ) return true;
    return false;
}

/************************************************************************
 * offEnd
 * Returns "true" if current is off end and "false" otherwise.  Exported.
 ************************************************************************/
template <typename T,typename D>
bool TreeSimpleVec<T,D>::offEnd(){
    
    return(current == nullptr);
}

/************************************************************************
 * getCurrent
 * Returns the particuls of current.  Exported.
 * Pre: !offEnd()
 ************************************************************************/
template <typename T,typename D>
void TreeSimpleVec<T,D>::getCurrent(IndexType *branch_index,IndexType *nparticles){
    
    if( offEnd() ){
        
        ERROR_MESSAGE(); fprintf(stderr, "Tree Error: calling getCurrent() when current is off end\n");
        exit(1);
    }
    
    *nparticles = current->nparticles;
    branch_index = current->branch_index;
    
    return;
}

/************************************************************************
 * getNbranches
 * Returns the NBranchVes of tree.  Exported.
 ************************************************************************/
template <typename T,typename D>
size_t TreeSimpleVec<T,D>::getNbranches(){
    
    return(Nbranches);
}

/***** Manipulation procedures *****/

/************************************************************************
 * moveTop
 * Moves current to the front of tree.  Exported.
 * Pre: !isEmpty()
 ************************************************************************/
template <typename T,typename D>
void TreeSimpleVec<T,D>::moveTop(){
   
    if( isEmpty() ){
        
        ERROR_MESSAGE(); fprintf(stderr, "Tree Error: calling moveTop() on empty tree\n");
        exit(1);
    }
    
    current = top_ptr;
    
	return;
}

/************************************************************************
 * movePrev
 * Moves current to the BranchV before it in tree.  This can move current
 * off end.  Exported.
 * Pre: !offEnd()
 ************************************************************************/
template <typename T,typename D>
void TreeSimpleVec<T,D>::moveUp(){
    
    if( offEnd() ){
        ERROR_MESSAGE(); fprintf(stderr, "Tree Error: call to moveUp() when current is off end\n");
        exit(1);
    }
    if( current == top_ptr ){
        ERROR_MESSAGE(); fprintf(stderr, "Tree Error: call to moveUp() tried to move off the top\n");
        exit(1);
    }
    
    current = current->prev_ptr;  /* can move off end */
    return;
}

/************************************************************************
 * moveToChild
 * Moves current to child BranchV after it in tree.  This can move current off
 * end.  Exported.
 * Pre: !offEnd()
 ************************************************************************/
template <typename T,typename D>
void TreeSimpleVec<T,D>::moveToChild(int child){
    
    if( offEnd() ){
        
        ERROR_MESSAGE(); fprintf(stderr, "Tree Error: calling moveChildren() when current is off end\n");
        exit(1);
    }
    if(child==1){
        if( current->child1_ptr == nullptr ){
            ERROR_MESSAGE(); fprintf(stderr, "Tree Error: moveToChild() typing to move to child1 when it doesn't exist\n");
            exit(1);
        }
        current = current->child1_ptr;
    }
    if(child==2){
        if( current->child2_ptr == nullptr ){
            ERROR_MESSAGE(); fprintf(stderr, "Tree Error: moveToChild() typing to move to child2 when it doesn't exist\n");
            exit(1);
        }
        current = current->child2_ptr;
    }
    return;
}

/************************************************************************
 * insertAfterCurrent
 * Inserts a new BranchV after the current BranchV in the tree and sets the
 * data field of the new BranchV to input.  Exported.
 * Pre: !offEnd()
 ************************************************************************/
template <typename T,typename D>
void TreeSimpleVec<T,D>::attachChildToCurrent(IndexType *branch_index,IndexType nparticles
                                           ,D boundary_p1[],D boundary_p2[]
                                           ,int child){
    
    /*printf("attaching child%i  current paricle number %i\n",child,current->nparticles);*/
  std::shared_ptr<BranchV> branchV(new BranchV(Ndimensions,branch_index,nparticles,boundary_p1,boundary_p2
                           ,current->level+1,Nbranches));
  
  ++Nbranches;
  
    
    if( offEnd() ){
        
        ERROR_MESSAGE(); fprintf(stderr, "Tree Error: calling attachChildToCurrent() when current is off end\n");
        exit(1);
    }
    
    branchV->prev_ptr = current;
    
    if(child==1){
        if(current->child1_ptr != nullptr){
            ERROR_MESSAGE(); fprintf(stderr, "Tree Error: calling attachChildToCurrent() when child1 alread exists\n");
            exit(1);
        }
        current->child1_ptr = branchV;
    }
    if(child==2){
        if(current->child2_ptr != nullptr){
            ERROR_MESSAGE();
            fprintf(stderr, "Tree Error: calling attachChildToCurrent() when child2 alread exists\n  current level=%i Nbranches=%li\n"
                    ,current->level,Nbranches);
            exit(1);
        }
        current->child2_ptr = branchV;
    }
  
    return;
}

/* same as above but takes a BranchV structure */
template <typename T,typename D>
void TreeSimpleVec<T,D>::attachChildToCurrent(BranchV &data,int child){
    
    attachChildToCurrent(data.branch_index,data.nparticles
                           ,data.boundary_p1,data.boundary_p2,child);
    return;
}

// step for walking tree by iteration instead of recursion
template <typename T,typename D>
bool TreeSimpleVec<T,D>::TreeWalkStep(bool allowDescent){
	if(allowDescent && current->child1_ptr != nullptr){
		moveToChild(1);
		return true;
	}
	if(allowDescent && current->child2_ptr != nullptr){
		moveToChild(2);
		return true;
	}
    
	if(current->brother_ptr != nullptr){
		current = current->brother_ptr;
		return true;
	}
	return false;
}

// = 0 if not possibly intersection
// > 0 sphere might intersect box
// = 2 sphere completely inside box
template <typename T,typename D>
int TreeSimpleVec<T,D>::cutbox(const D *center,D *p1,D *p2,float rmax){
  int intersection = 1,inside=1;
  for(int i=0 ; i<Ndimensions ; ++i){
    intersection *= ( (center[i] + rmax) > p1[i] )*( (center[i] - rmax) < p2[i] );
    inside *= ( (center[i] - rmax) > p1[i] )*( (center[i] + rmax) < p2[i] );
  }
  
  return intersection + inside;
}
#endif /* SIMP_TREE_V_ */
