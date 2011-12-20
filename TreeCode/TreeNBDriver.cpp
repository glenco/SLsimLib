/*#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <nrutil.h>
#include <nr.h>
#include <nrD.h>
#include "TreeNB.h"*/

#include <slsimlib.h>

#define Nbucket 5


/*
 *  median_cutNB determins how the cells are subdivided */
/*    if ==0  equal volume cuts */
/*    if ==1  median particle cuts */
static int incellNB,median_cutNB=0,Ndim;
static PosType realrayNB[treeNBdim];



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
  std::printf("projected with 2D tree\n");

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
