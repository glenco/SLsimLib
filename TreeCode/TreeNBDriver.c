#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
/*#include "double_sort.c"*/
#include <nrutil.h>
#include <nr.h>
#include <nrD.h>
#include "TreeNB.h"

#define Nbucket 5


/*
 *  median_cutNB determins how the cells are subdivided */
/*    if ==0  equal volume cuts */
/*    if ==1  median particle cuts */
static int incellNB,median_cutNB=0,Ndim;
static PosType realrayNB[treeNBdim];

/** \ingroup DeflectionL2
 *
 */
TreeNBHndl BuildTreeNB(PosType **xp,float *rsph,float *mass,bool MultiRadius
	   ,bool MultiMass,IndexType Nparticles,IndexType *particles,int Ndimensions
	   ,double theta){
  TreeNBHndl tree;
  IndexType i,j;
  PosType p1[3],p2[3],center[3];
  void _BuildTreeNB(TreeNBHndl tree,IndexType nparticls,IndexType *particles);

  Ndim=Ndimensions;
  if(Ndim > treeNBdim){
	  ERROR_MESSAGE();
	  printf("ERROR: BuildTree, need to change treeNBdim in TreeNB.h to make 3-d tree\n");
	  exit(0);
  }

  for(j=0;j<Ndim;++j){
    p1[j]=xp[0][j];
    p2[j]=xp[0][j];
  }

  for(i=0;i<Nparticles;++i){
    particles[i]=i;
    
    for(j=0;j<Ndim;++j){
      if(xp[i][j] < p1[j] ) p1[j]=xp[i][j];
      if(xp[i][j] > p2[j] ) p2[j]=xp[i][j];
    }
  }

  for(j=0;j<Ndim;++j) center[j]=(p1[j]+p2[j])/2;

  /* Initialize tree root */
  tree=NewTreeNB(particles,Nparticles,p1,p2,center,Ndim);

  tree->MultiMass = MultiMass;
  tree->MultiRadius = MultiRadius;
  tree->xp = xp;
  tree->rsph = rsph;
  tree->mass = mass;

  /* build the tree */
  _BuildTreeNB(tree,Nparticles,particles);

  /* visit every branch to find center of mass and cutoff scale */
   moveTopNB(tree);
   cuttoffscale(tree,&theta);

  return tree;
}


// tree must be created and first branch must be set before start
void _BuildTreeNB(TreeNBHndl tree,IndexType nparticles,IndexType *particles){
  IndexType i,j,cut,dimension;
  BranchNB *cbranch,branch1,branch2;
  PosType xcut,xcm[2];
  double *x;

  cbranch=tree->current; /* pointer to current branch */

  cbranch->rmax = sqrt( pow(cbranch->boundary_p2[0]-cbranch->boundary_p1[0],2)
  		     + pow(cbranch->boundary_p2[1]-cbranch->boundary_p1[1],2) );

  for(i=0,cbranch->mass=0;i<cbranch->nparticles;++i)
 	  cbranch->mass+=tree->mass[cbranch->particles[i]*tree->MultiMass];
  cbranch->center[0]=cbranch->center[1]=0;
  for(i=0;i<cbranch->nparticles;++i){
	  cbranch->center[0] += tree->mass[cbranch->particles[i]*tree->MultiMass]
	                       *tree->xp[cbranch->particles[i]][0]/cbranch->mass;
	  cbranch->center[1] += tree->mass[cbranch->particles[i]*tree->MultiMass]
	                       *tree->xp[cbranch->particles[i]][1]/cbranch->mass;
  }
  //////////////////////////////////////////////
  // calculate quadropole moment of branch
  //////////////////////////////////////////////
  cbranch->quad[0]=cbranch->quad[1]=cbranch->quad[2]=0;
  for(i=0;i<cbranch->nparticles;++i){
	  xcm[0]=tree->xp[cbranch->particles[i]][0]-cbranch->center[0];
	  xcm[1]=tree->xp[cbranch->particles[i]][1]-cbranch->center[1];
	  xcut=pow(xcm[0],2) + pow(xcm[1],2);
 	  cbranch->quad[0] += (xcut-2*xcm[0]*xcm[0])*tree->mass[cbranch->particles[i]];
 	  cbranch->quad[1] += (xcut-2*xcm[1]*xcm[1])*tree->mass[cbranch->particles[i]];
 	  cbranch->quad[2] += -2*xcm[0]*xcm[1]*tree->mass[cbranch->particles[i]];
  }


  /* leaf case */
  if(cbranch->nparticles <= Nbucket){
	  cbranch->big_particle=0;
	  return;
  }
 
  /* initialize boundaries to old boundaries */
  for(i=0;i<Ndim;++i){
      branch1.boundary_p1[i]=cbranch->boundary_p1[i];
      branch1.boundary_p2[i]=cbranch->boundary_p2[i];

      branch2.boundary_p1[i]=cbranch->boundary_p1[i];
      branch2.boundary_p2[i]=cbranch->boundary_p2[i];
  }
/*
  for(i=cbranch->nparticles-1,cbranch->big_particle=0
		 ; i > cbranch->big_particle-1 ; --i){
	  //printf("i=%li  %li\n",i,cbranch->big_particle);
	  if(2*tree->rsph[tree->MultiRadius*particles[i]] > cbranch->rmax){

		  cut=particles[cbranch->big_particle];
		  particles[cbranch->big_particle]=particles[i];
		  particles[i]=cut;

		  ++i;
		  ++(cbranch->big_particle);
	  }
  }

  if(cbranch->big_particle==cbranch->nparticles) return;
*/
  cbranch->big_particle=0;

  //**** makes sure force does not require nbucket at leaf

  /* set dimension to cut box */
  dimension=(cbranch->level % Ndim);

  x=(double *)malloc((cbranch->nparticles-cbranch->big_particle)*sizeof(double));
  for(i=cbranch->big_particle;i<cbranch->nparticles;++i) x[i]=tree->xp[particles[i]][dimension];

 //double_sort(cbranch->nparticles-cbranch->big_particle,x-1,particles-1+cbranch->big_particle);

  if(median_cutNB){
	  quicksort(particles,x,cbranch->nparticles-cbranch->big_particle);

	  cut=(cbranch->nparticles-cbranch->big_particle)/2;
      branch1.boundary_p2[dimension]=x[cut];
      branch2.boundary_p1[dimension]=x[cut];
  }else{
      xcut=(cbranch->boundary_p1[dimension]+cbranch->boundary_p2[dimension])/2;
      branch1.boundary_p2[dimension]=xcut;
      branch2.boundary_p1[dimension]=xcut;

      quickPartition(xcut,&cut,particles
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

  // centers of mass
/*  for(i=0;i<Ndim;++i) branch1.center[i]=0;
  for(i=0;i<cut; ++i) for(j=0;j<Ndim;++j)
 		       branch1.center[j]+=xp[particles[i]][j]/branch1.nparticles;
  for(i=0;i<Ndim;++i) branch2.center[i]=0;
  for(i=cut;i<cbranch->nparticles; ++i) for(j=0;j<Ndim;++j)
 		       branch2.center[j]+=xp[particles[i]][j]/branch2.nparticles;
*/

  // centers of mass

  //cbranch->child1->mass+=masses[branch1.particles[i]];

  for(i=0;i<Ndim;++i) branch1.center[i]=0;
  for(i=0;i<branch1.nparticles; ++i) for(j=0;j<Ndim;++j)
 		       branch1.center[j]+=tree->xp[branch1.particles[i]][j]/branch1.nparticles;

  for(i=0;i<Ndim;++i) branch2.center[i]=0;
  for(i=0;i<branch2.nparticles; ++i) for(j=0;j<Ndim;++j)
 		       branch2.center[j]+=tree->xp[branch2.particles[i]][j]/branch2.nparticles;


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
      //attachChildToCurrentNB(tree,branch1,1);
      moveToChildNB(tree,1);
     _BuildTreeNB(tree,branch1.nparticles,branch1.particles);
     moveUpNB(tree);
 }

 if(branch2.nparticles > 0 ){ 
     //attachChildToCurrentNB(tree,branch2,2);
     moveToChildNB(tree,2);
     _BuildTreeNB(tree,branch2.nparticles,branch2.particles);
     moveUpNB(tree);
 }


 /*printf("reached end of _BuildTreeNB level=%i\n",tree->current->level);*/
 return;
}

/** \ingroup DeflectionL2
 *  finds the nearest neighbors in whatever dimensions tree is defined in
 *  */
IndexType *NearestNeighborNB(TreeNBHndl tree,double *ray,int Nneighbors
				 ,float *rsph){
  IndexType i;
  void _NearestNeighborNB(TreeNBHndl tree,double *ray,int Nneighbors,IndexType *neighbors,double *rneighbor);
  static int count=0,oldNneighbors=-1;
  static double *rneighbors;
  static IndexType *neighbors;

  Ndim=tree->Ndimensions;

  if(tree->top->nparticles <= Nneighbors){
	ERROR_MESSAGE();
    printf("ERROR: in NearestNeighborNB, number of neighbors > total number of particles\n");
    exit(1);
  }

  if(count==0){
    /*printf("allocating memory\n");*/
    neighbors=(IndexType *)malloc((Nneighbors+Nbucket)*sizeof(IndexType));
    rneighbors=(double *)malloc((Nneighbors+Nbucket)*sizeof(double));
    ++count;
    oldNneighbors=Nneighbors;
  }else if(oldNneighbors < Nneighbors){ /* if the number of nearest neighbors goes up get more mem */
    /*printf("re-allocating memory\n");*/
    neighbors=(IndexType *)realloc(neighbors,(Nneighbors+Nbucket)*sizeof(IndexType));
    rneighbors=(double *)realloc(rneighbors,(Nneighbors+Nbucket)*sizeof(double));
    oldNneighbors=Nneighbors;
  }

  /*   printf("Nneighbors=%i\n",Nneighbors); */
  /*   printf("array sizes=%i\n",Nneighbors+Nbucket); */

  /* initalize distance to neighbors to a large number */
  for(i=0;i<Nbucket+Nneighbors;++i){
    rneighbors[i]=10*(tree->top->boundary_p2[0]-tree->top->boundary_p1[0]);
    neighbors[i]=0;
  }

  for(i=0;i<Ndim;++i) realrayNB[i]=ray[i];

  moveTopNB(tree);
  if( inboxNB(ray,tree->current->boundary_p1,tree->current->boundary_p2) == 0 ){
	  ERROR_MESSAGE();
    /*printf("ERROR: in NearestNeighborNB, ray is not inside the simulation box\n      ray=%e %e  boundary= %e %e   %e %e\n",ray[0],ray[1],tree->current->boundary_p1[0],tree->current->boundary_p1[1]
      ,tree->current->boundary_p1[0],tree->current->boundary_p2[1]);*/

    for(i=0;i<Ndim;++i){
      ray[i]=DMAX(ray[i],tree->current->boundary_p1[i]);
      ray[i]=DMIN(ray[i],tree->current->boundary_p2[i]);
    }
  }
  incellNB=1;

  _NearestNeighborNB(tree,ray,Nneighbors,neighbors,rneighbors);

  /*for(i=0;i<Nneighbors;++i) printf("%e %e %e\n",xp[neighbors[i]][0],xp[neighbors[i]][1],xp[neighbors[i]][2]);*/

  *rsph=(float)(rneighbors[Nneighbors-1]);
  /*printf("hi\n");*/
  return neighbors;
}


void _NearestNeighborNB(TreeNBHndl tree,double *ray,int Nneighbors
			,IndexType *neighbors,double *rneighbors){

  int i,j,incellNB2=1;

  if(incellNB){  /* not found cell yet */

    if( inboxNB(ray,tree->current->boundary_p1,tree->current->boundary_p2) ){

      /* found the box small enough */
    	if( tree->current->nparticles <= Nneighbors+Nbucket ){
    		incellNB=0;
    		for(i=0;i<Ndim;++i) ray[i]=realrayNB[i];

    		/* calculate the distance to all the particles in cell */
    		for(i=0;i<tree->current->nparticles;++i){
    			for(j=0,rneighbors[i]=0.0;j<Ndim;++j){
    				rneighbors[i]+=pow(tree->xp[tree->current->particles[i]][j]-ray[j],2);
    			}
    			rneighbors[i]=sqrt( rneighbors[i] );
    			neighbors[i]=tree->current->particles[i];
    		}

    		/*printf("first sort at level =%i\n",tree->current->level);*/
    		/*printf("N=%i\n",tree->current->nparticles);*/
    		/*for(i=0;i<tree->current->nparticles;++i) printf("  %i  %e\n",neighbors[i],rneighbors[i]);*/
    		double_sort(tree->current->nparticles,rneighbors-1,neighbors-1);
    		/*printf("end sort\n");*/
       
    		/*for(i=0;i<tree->current->nparticles;++i) printf("  %i  %e\n",neighbors[i],rneighbors[i]);*/

      }else{ /* keep going down the tree */

    	  /*printf("moving to child1 from level %i\n",tree->current->level);*/
    	  if(tree->current->child1 !=NULL){
    		  moveToChildNB(tree,1);
    		  _NearestNeighborNB(tree,ray,Nneighbors,neighbors,rneighbors);
    		  /*printf("moving up from level %i\n",tree->current->level);*/
    		  moveUpNB(tree);

    		  incellNB2=incellNB;
    	  }

    	  if(tree->current->child2 !=NULL){
    		  /*printf("moving to child2 from level %i\n",tree->current->level);*/
    		  moveToChildNB(tree,2);
    		  _NearestNeighborNB(tree,ray,Nneighbors,neighbors,rneighbors);
    		  /*printf("moving up from level %i\n",tree->current->level);*/
    		  moveUpNB(tree);
    	  }

    	  /** if ray found in second child go back to first to search for neighbors **/
    	  if( (incellNB2==1) && (incellNB==0) ){
    		  if(tree->current->child1 !=NULL){
    			  /*printf("moving to child1 again from level %i\n",tree->current->level);*/
    			  moveToChildNB(tree,1);
    			  _NearestNeighborNB(tree,ray,Nneighbors,neighbors,rneighbors);
    			  /*printf("moving up from level %i\n",tree->current->level);*/
    			  moveUpNB(tree);
    		  }
    	  }
      }
	}/*else{ // not in the box
		//printf("not in box \n");}
	  }else{   found cell */

    /*printf("finding neighboring boxes at level = %i\n",tree->current->level);*/

		/* does radius cut into the box */
    if( cutboxNB(ray,tree->current->boundary_p1,tree->current->boundary_p2,rneighbors[Nneighbors-1]) ){

    	if( (tree->current->child1 == NULL)*(tree->current->child2 == NULL)){  /* leaf case */

    		/* combine found neighbors with particles in box and resort */
    		for(i=Nneighbors;i<(tree->current->nparticles+Nneighbors);++i){
    			for(j=0,rneighbors[i]=0.0;j<Ndim;++j){
    				rneighbors[i]+=pow(tree->xp[tree->current->particles[i-Nneighbors]][j]-ray[j],2);
	  			}
    			rneighbors[i]=sqrt( rneighbors[i] );
    			neighbors[i]=tree->current->particles[i-Nneighbors];
    		}

    		double_sort(Nneighbors+Nbucket,rneighbors-1,neighbors-1);

    	}else{

    		if(tree->current->child1 !=NULL){
    			moveToChildNB(tree,1);
    			_NearestNeighborNB(tree,ray,Nneighbors,neighbors,rneighbors);
    			moveUpNB(tree);
    		}

    		if(tree->current->child2 !=NULL){
    			moveToChildNB(tree,2);
    			_NearestNeighborNB(tree,ray,Nneighbors,neighbors,rneighbors);
    			moveUpNB(tree);
    		}
    	}
	}
  }
  return;
}

/* return 1 (0) if ray is (not) in the cube */
int inboxNB(double *ray,PosType *p1,PosType *p2){
  short i,ans;

  for(i=1,ans=(ray[0]>=p1[0])*(ray[0]<=p2[0]);i<Ndim;++i) 
    ans*=(ray[i]>=p1[i])*(ray[i]<=p2[i]);

  return ans;
}

/* return 1 (0) if box is (not) within rmax of ray */
int cutboxNB(double *ray,PosType *p1,PosType *p2,PosType rmax){
  int i;
  PosType close[3],rtmp;
  
  for(i=0;i<Ndim;++i){

    if( ray[i] < p1[i] ){
      close[i]=p1[i];
    }else if(ray[i] > p2[i]){
      close[i]=p2[i];
    }else{
      close[i]=ray[i];
    }
  }
  
  for(i=0,rtmp=0;i<Ndim;++i) rtmp+= pow(ray[i] - close[i],2);
 
  if(rtmp<rmax*rmax) return 1;
  return 0;
}

/** \ingroup DeflectionL2
 *  rotate particle positions and build 3d tree
 *
 */
TreeNBHndl rotate_simulation(PosType **xp,IndexType Nparticles,IndexType *particles
		,double **coord,double theta,float *rsph,float *mass
		,bool MultiRadius,bool MultiMass){
 	TreeNBHndl tree;
	IndexType j;
	int i;
	PosType tmp[3];
  
  /* rotate particle positions */
  for(i=0;i<Nparticles;++i){
    for(j=0;j<3;++j) tmp[j]=0.0;
    for(j=0;j<3;++j){
      tmp[0]+=coord[0][j]*xp[i][j];
      tmp[1]+=coord[1][j]*xp[i][j];
      tmp[2]+=coord[2][j]*xp[i][j];
    }
    for(j=0;j<3;++j) xp[i][j]=tmp[j];
  }

  // need to erase tree !!!!
  tree=BuildTreeNB(xp,rsph,mass,MultiRadius,MultiMass,Nparticles,particles,3,theta);

  return tree;
}

/** \ingroup DeflectionL2
 * rotate particle positions and build 2d tree in x-y plane
 */
TreeNBHndl rotate_project(PosType **xp,IndexType Nparticles,IndexType *particles
		,double **coord,double theta,float *rsph,float *mass
		,bool MultiRadius,bool MultiMass){
   IndexType j;
  int i;
  PosType tmp[3];
  TreeNBHndl tree;

  /* rotate particle positions */
   for(i=0;i<Nparticles;++i){
     for(j=0;j<3;++j) tmp[j]=0.0;
     for(j=0;j<3;++j){
       tmp[0]+=coord[0][j]*xp[i][j];
       tmp[1]+=coord[1][j]*xp[i][j];
       tmp[2]+=coord[2][j]*xp[i][j];
     }
     for(j=0;j<3;++j) xp[i][j]=tmp[j];
   }

  tree=BuildTreeNB(xp,rsph,mass,MultiRadius,MultiMass,Nparticles,particles,2,theta);
  printf("projected with 2D tree\n");

  return tree;
}

/** \ingroup DeflectionL2
\brief !!!! Not yet come up with a good way of doing this !!!!

	 * gives the particles a third dimension depending on their size
	 *  This is used to spread subhalos out into 3d so that their density
	 *  is about one per smoothing length.  This is to avoid the big particle
	 *  problem that slows the tree-force calculation.
	 *
	 *  Warning: This will erase the third coordinate of the particles.
	 */
TreeNBHndl spread_particles(PosType **xp,IndexType Nparticles,IndexType *particles
		,double theta,float *rsph,float *mass
		,bool MultiRadius,bool MultiMass){


	IndexType i,*dummy;
	TreeNBHndl tree;
	float tmp=0;

	//Build 2d tree for
	tree = BuildTreeNB(xp,&tmp,mass,false,false,Nparticles,particles,2,theta);

	for(i=0;i<Nparticles;++i){
		dummy = NearestNeighborNB(tree,xp[i],1,&tmp);
		xp[i][2] = 4*pi*pow(rsph[i],3)/pow(tmp,2)/3;
	}

	// free tree
	freeTreeNB(tree);

	// Build 3d tree for force calculations
	tree = BuildTreeNB(xp,rsph,mass,MultiRadius,MultiMass,Nparticles,particles,3,theta);

	return tree;
}

/** \ingroup DeflectionL2
 * visits every branch in tree to calculate
 * two critical lengths rcrit_angle and rcrit_part
 *   and mark largest particle in each node subject
 *   to some conditions
 *
 *   if the sph[] are negative rcrit_part = 0
 */
void cuttoffscale(TreeNBHndl tree,double *theta){

	IndexType i,prev_big;
	PosType rcom,maxrsph,prev_size;
	BranchNB *cbranch;

	moveTopNB(tree);
	do{
		cbranch=tree->current;

		// largest distance from center of mass of cell
		for(i=0,rcom=0.0;i<2;++i) rcom+=DMAX( pow(cbranch->center[i]-cbranch->boundary_p1[i],2)
 				       ,pow(cbranch->center[i]-cbranch->boundary_p2[i],2) );
		rcom=sqrt(rcom);

		if(*theta > 0.0) cbranch->rcrit_angle=1.15470*rcom/(*theta);
		else  cbranch->rcrit_angle=1.0e100;

		//////////////////////////////////////////////
		/* find the biggest particle that is not the
		 * biggest particle of any of its parents
		///////////////////////////////////////////////*/

		if(tree->MultiRadius){

			if(atTopNB(tree) || cbranch->prev->rcrit_part==0.0){
				prev_big = -1;
				prev_size = 1.0e100;
			}else{
				prev_big = cbranch->prev->big_particle;
				prev_size = cbranch->prev->maxrsph;
			}

			//assert(tree->MultiRadius == 0);

			// find largest rsph smoothing in cell  that is not the largest in a previous cell
			cbranch->big_particle = cbranch->particles[0];
			for(i=0,maxrsph=0.0,cbranch->big_particle=0;i<cbranch->nparticles;++i){

				assert(i < tree->top->nparticles);

				if(maxrsph <= tree->rsph[tree->MultiRadius*cbranch->particles[i]] ){
					maxrsph = tree->rsph[tree->MultiRadius*cbranch->particles[i]];
					cbranch->big_particle = cbranch->particles[i];
				}
			}
			// if it is possible for the largest particle size to be big enough
			//    cause the opening of a box then tag it
			cbranch->maxrsph = maxrsph;
			//printf("set maxrsp=%e\n",cbranch->maxrsph);
			//if(maxrsph > rcom*(1.0/(*theta)-1)/2){

			if(maxrsph > cbranch->rcrit_angle/2){
				tree->rsph[tree->MultiRadius*cbranch->big_particle] = 0.0; // zero out so it wont show up in leaf
				cbranch->rcrit_part = rcom + 2*maxrsph;
			}else{
				cbranch->rcrit_part = 0;
			}

		}else{
			// point mass case
			cbranch->maxrsph = tree->rsph[0];
			cbranch->rcrit_part = rcom + 2*cbranch->maxrsph;
			cbranch->big_particle = cbranch->particles[0];
		}
		cbranch->rcrit_angle += cbranch->rcrit_part;

	}while(TreeNBWalkStep(tree,true));

	return;
}
