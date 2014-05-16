
/*
 * Code Name:     tree.c                                       
 * Programmer:    R Ben Metcalf
 * Last Revised:  Nov, 2005                                   
 * Discription:  2 way tree data structure with 2 branchNBes
 * Comments:                           
 */

#include "slsimlib.h"


//PosType dummy;

BranchNB::BranchNB(int Ndim){
	center = new PosType[Ndim*3];
	boundary_p1 = &center[Ndim];
	boundary_p2 = &center[2*Ndim];
}
BranchNB::~BranchNB(){
	delete center;
}

TreeSimple::TreeSimple(PosType **xpt,IndexType Npoints,int bucket,int Ndimensions,bool median){
	index = new IndexType[Npoints];
	IndexType ii;

	Nbucket = bucket;
	Ndim = Ndimensions;
	median_cut = median;
	Nparticles = Npoints;

	xp = xpt;

	for(ii=0;ii<Npoints;++ii) index[ii] = ii;

	tree = TreeSimple::BuildTreeNB(xp,Npoints,index,Ndimensions,0);

	return;
}

TreeSimple::~TreeSimple()
{
	freeTreeNB(tree);
	delete[] index;
	return;
}

void TreeSimple::PointsWithinEllipse(
	PosType *ray     /// center of ellipse
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
		x = xp[*it][0]*cs - xp[*it][1]*sn;
		y = xp[*it][0]*sn + xp[*it][1]*cs;
		if( ( pow(x/rmax,2) + pow(y/rmin,2) ) > 1) it = neighborlist.erase(it);
		else ++it;
	}
	return;
}
void TreeSimple::PointsWithinCircle(
		PosType *ray     /// center of circle
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

  _PointsWithin(ray,&rmax,neighborlist);

  return;
}
/**
 * Used in PointsWithinKist() to walk tree.*/
void TreeSimple::_PointsWithin(PosType *ray,float *rmax,std::list <unsigned long> &neighborlist){

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
    			  for(j=0,radius=0.0;j<2;++j) radius+=pow(xp[tree->current->particles[i]][j]-ray[j],2);
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
				  for(j=0,radius=0.0;j<2;++j) radius+=pow(xp[tree->current->particles[i]][j]-ray[j],2);
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
 *  \brief finds the nearest neighbors in whatever dimensions tree is defined in
 *  */
void TreeSimple::NearestNeighbors(
            PosType *ray       /// position
            ,int Nneighbors    /// number of neighbors to be found
            ,float *radius     /// distance of furthest neighbor found from ray[]
            ,IndexType *neighborsout  /// list of the indexes ofx the neighbors
                                  ){
  IndexType i;
  //static int count=0,oldNneighbors=-1;
   short j;

  Ndim = tree->Ndimensions;

  PosType rneighbors[Nneighbors+Nbucket];
  IndexType neighbors[Nneighbors+Nbucket];

  if(tree->top->nparticles <= Nneighbors){
	ERROR_MESSAGE();
    printf("ERROR: in TreeSimple::NearestNeighbors, number of neighbors > total number of particles\n");
    exit(1);
  }

  /* initalize distance to neighbors to a large number */
  for(i=0;i<Nbucket+Nneighbors;++i){
    rneighbors[i] = (10*(tree->top->boundary_p2[0]-tree->top->boundary_p1[0]));
    neighbors[i] = 0;
  }

  for(j=0;j<Ndim;++j) realray[j]=ray[j];

  moveTopNB(tree);
  if( inbox(ray,tree->current->boundary_p1,tree->current->boundary_p2) == 0 ){
	  ERROR_MESSAGE();

    for(j=0;j<Ndim;++j){
        ray[j] = (ray[j] > tree->current->boundary_p1[j]) ? ray[j] : tree->current->boundary_p1[j];
        ray[j] = (ray[j] < tree->current->boundary_p2[j]) ? ray[j] : tree->current->boundary_p2[j];

        //ray[j]=DMAX(ray[j],tree->current->boundary_p1[j]);
        //ray[j]=DMIN(ray[j],tree->current->boundary_p2[j]);
    }
  }
  incell = 1;

  _NearestNeighbors(ray,Nneighbors,neighbors,rneighbors);

  for(i=0;i<Nneighbors;++i) neighborsout[i] = neighbors[i];
  *radius = rneighbors[Nneighbors-1];

  return;
}


void TreeSimple::_NearestNeighbors(PosType *ray,int Nneighbors,unsigned long *neighbors,PosType *rneighbors){

  int incellNB2=1;
  IndexType i;
  short j;

  if(incell){  /* not found cell yet */

    if( inbox(ray,tree->current->boundary_p1,tree->current->boundary_p2) ){

      /* found the box small enough */
    	if( tree->current->nparticles <= Nneighbors+Nbucket ){
    		incell=0;
    		for(j=0;j<Ndim;++j) ray[j]=realray[j];

    		/* calculate the distance to all the particles in cell */
    		for(i=0;i<tree->current->nparticles;++i){
    			for(j=0,rneighbors[i]=0.0;j<Ndim;++j){
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
    			for(j=0,rneighbors[i]=0.0;j<Ndim;++j){
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


TreeNBHndl TreeSimple::BuildTreeNB(PosType **xp,IndexType Nparticles,IndexType *particles,int Ndims
	   ,PosType theta){
  TreeNBHndl tree;
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

  tree->xp = xp;

  /* build the tree */
  _BuildTreeNB(tree,Nparticles,particles);

  /* visit every branch to find center of mass and cutoff scale */
   moveTopNB(tree);

  return tree;
}


// tree must be created and first branch must be set before start
void TreeSimple::_BuildTreeNB(TreeNBHndl tree,IndexType nparticles,IndexType *particles){
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
  for(i=cbranch->big_particle;i<cbranch->nparticles;++i) x[i]=tree->xp[particles[i]][dimension];

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
BranchNB *TreeSimple::NewBranchNB(IndexType *particles,IndexType nparticles
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
void TreeSimple::FreeBranchNB(BranchNB *branchNB){

    assert( branchNB != NULL);
    delete branchNB;

    return;
}

/************************************************************************
 * NewTreeNB
 * Returns pointer to new TreeNB struct.  Initializes top, last, and
 * current pointers to NULL.  Sets NbranchNBes field to 0.  Exported.
 ************************************************************************/
TreeNBHndl TreeSimple::NewTreeNB(IndexType *particles,IndexType nparticles
		 ,PosType boundary_p1[],PosType boundary_p2[],
		     PosType center[],short Ndimensions){

    TreeNBHndl tree;
    
    tree = (TreeNBStruct *)malloc(sizeof(TreeNBStruct));
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

void TreeSimple::freeTreeNB(TreeNBHndl tree){
	/* free treeNB
	 *  does not free the particle positions, masses or rsph
	 */

	if(tree == NULL) return;

	emptyTreeNB(tree);
  	FreeBranchNB(tree->top);
	free(tree);

	return;
}

short TreeSimple::emptyTreeNB(TreeNBHndl tree){

	moveTopNB(tree);
	_freeTreeNB(tree,0);

	assert(tree->Nbranches == 1);

	return 1;
}

void TreeSimple::_freeTreeNB(TreeNBHndl tree,short child){
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
bool TreeSimple::isEmptyNB(TreeNBHndl tree){

    assert(tree != NULL);
    return(tree->Nbranches == 0);
}

/************************************************************************
 * atTop
 * Returns "true" if current is the same as top and "false" otherwise.
 * Exported.
 * Pre: !isEmptyNB(tree)
 ************************************************************************/
bool TreeSimple::atTopNB(TreeNBHndl tree){

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
bool TreeSimple::noChildNB(TreeNBHndl tree){

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
bool TreeSimple::offEndNB(TreeNBHndl tree){

    assert(tree != NULL);
    return(tree->current == NULL);
}

/************************************************************************
 * getCurrentNB
 * Returns the particuls of current.  Exported.
 * Pre: !offEndNB(tree)
 ************************************************************************/
void TreeSimple::getCurrentNB(TreeNBHndl tree,IndexType *particles,IndexType *nparticles){

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
unsigned long TreeSimple::getNbranchesNB(TreeNBHndl tree){

    assert(tree != NULL);
    return(tree->Nbranches);
}

/***** Manipulation procedures *****/

/************************************************************************
 * moveTopNB
 * Moves current to the front of tree.  Exported.
 * Pre: !isEmptyNB(tree)
 ************************************************************************/
void TreeSimple::moveTopNB(TreeNBHndl tree){
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
void TreeSimple::moveUpNB(TreeNBHndl tree){
    
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
void TreeSimple::moveToChildNB(TreeNBHndl tree,int child){

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
void TreeSimple::insertChildToCurrentNB(TreeNBHndl tree, IndexType *particles,IndexType nparticles
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

void TreeSimple::attachChildToCurrentNB(TreeNBHndl tree,BranchNB &data,int child){

  insertChildToCurrentNB(tree,data.particles,data.nparticles
		  ,data.boundary_p1,data.boundary_p2,data.center,child);
  return;
}

// step for walking tree by iteration instead of recursion
bool TreeSimple::TreeNBWalkStep(TreeNBHndl tree,bool allowDescent){
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

