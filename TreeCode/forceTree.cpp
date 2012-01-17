/*
 * forceTree.cpp
 *
 *  Created on: Oct 18, 2011
 *      Author: bmetcalf
 */
#include <assert.h>
#include <math.h>
#include "forceTree.h"

const float pi = 3.141593;

/// Gaussian particles
double alpha_o(double r2,float sigma){
  if(r2==0) return 0.0;
  if(sigma == 0.0 ) return -1.0/r2/pi;
  if(sigma < 0.0)  return -exp(-r2/sigma/sigma/2)/r2/pi;  // screened potential
  return ( exp(-r2/sigma/sigma/2) - 1.0 )/r2/pi;
}
double kappa_o(double r2,float sigma){
  if(sigma == 0.0) return 0.0;
  if(sigma < 0.0) return -exp(-r2/sigma/sigma/2)/2/pi/sigma/sigma;
  return exp(-r2/sigma/sigma/2)/2/pi/sigma/sigma;
}
double gamma_o(double r2,float sigma){
  if(r2==0) return 0.0;
  if(sigma == 0.0) return -2.0/pi/pow(r2,2);
  if(sigma < 0.0) return -(-2.0 + (2.0+(r2/sigma/sigma))*exp(-r2/sigma/sigma/2) )/pi/pow(r2,2);
  return (-2.0 + (2.0 + r2/sigma/sigma)*exp(-r2/sigma/sigma/2) )/pi/pow(r2,2);
}

ForceTree::ForceTree(PosType **xp,IndexType Npoints,float *masses,float *rsph,bool Multimass,bool Multisize
		,int bucket,int dimension,bool median,double theta) :
	SimpleTree(xp,Npoints,bucket,dimension,median)
{
	tree->MultiMass = Multimass;
	tree->MultiRadius = Multisize;
	tree->masses = masses;
	tree->rsph = rsph;
	force_theta = theta;

	// This should be changed so other particle profiles can be used.
	alpha_internal = alpha_o;
	kappa_internal = kappa_o;
	gamma_internal = gamma_o;

	CalcMoments();
}

ForceTree::~ForceTree(){
}

// calculates moments of the mass in the squares and the cutoff scale for each box
void ForceTree::CalcMoments(){

	IndexType i,prev_big;
	PosType rcom,maxrsph,prev_size,xcm[2],xcut;
	BranchNB *cbranch;

	moveTopNB(tree);
	do{
		cbranch=tree->current; /* pointer to current branch */

		cbranch->rmax = sqrt( pow(cbranch->boundary_p2[0]-cbranch->boundary_p1[0],2)
				+ pow(cbranch->boundary_p2[1]-cbranch->boundary_p1[1],2) );

		// calculate mass
		for(i=0,cbranch->mass=0;i<cbranch->nparticles;++i)
			cbranch->mass += tree->masses[cbranch->particles[i]*tree->MultiMass];
		  // calculate center of mass
		cbranch->center[0]=cbranch->center[1]=0;
		for(i=0;i<cbranch->nparticles;++i){
			cbranch->center[0] += tree->masses[cbranch->particles[i]*tree->MultiMass]
			                       *tree->xp[cbranch->particles[i]][0]/cbranch->mass;
			cbranch->center[1] += tree->masses[cbranch->particles[i]*tree->MultiMass]
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
			cbranch->quad[0] += (xcut-2*xcm[0]*xcm[0])*tree->masses[cbranch->particles[i]*tree->MultiMass];
			cbranch->quad[1] += (xcut-2*xcm[1]*xcm[1])*tree->masses[cbranch->particles[i]*tree->MultiMass];
			cbranch->quad[2] += -2*xcm[0]*xcm[1]*tree->masses[cbranch->particles[i]*tree->MultiMass];
		}

		// largest distance from center of mass of cell
		for(i=0,rcom=0.0;i<2;++i) rcom += ( pow(cbranch->center[i]-cbranch->boundary_p1[i],2) > pow(cbranch->center[i]-cbranch->boundary_p2[i],2) ) ?
		pow(cbranch->center[i]-cbranch->boundary_p1[i],2) : pow(cbranch->center[i]-cbranch->boundary_p2[i],2);

		rcom=sqrt(rcom);

		if(force_theta > 0.0) cbranch->rcrit_angle = 1.15470*rcom/(force_theta);
		else  cbranch->rcrit_angle=1.0e100;


		if(tree->MultiRadius){

			for(i=0,cbranch->maxrsph=0.0;i<cbranch->nparticles;++i){
				if(cbranch->maxrsph <= tree->rsph[cbranch->particles[i]*tree->MultiRadius] ){
					cbranch->maxrsph = tree->rsph[cbranch->particles[i]*tree->MultiRadius];
				}
			}

			cbranch->rcrit_part = rcom + 2*cbranch->maxrsph;

			//////////////////////////////////////////////
			/* find the biggest particle that is not the
			 * biggest particle of any of its parents
			 *
			 * This was to implement the big particle subtraction scheme that
			 * has been semi-abandoned.
			///////////////////////////////////////////////

			if(atTopNB(tree) || cbranch->prev->rcrit_part==0.0){
				prev_big = -1;
				prev_size = 1.0e100;
			}else{
				prev_big = cbranch->prev->big_particle;
				prev_size = cbranch->prev->maxrsph;
			}

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

			if(maxrsph > cbranch->rcrit_angle/2){
				tree->rsph[tree->MultiRadius*cbranch->big_particle] = 0.0; // zero out so it wont show up in leaf
				cbranch->rcrit_part = rcom + 2*maxrsph;
			}else{
				cbranch->rcrit_part = 0;
			}
			 */

		}else{
			// single size case
			cbranch->maxrsph = tree->rsph[0];
			cbranch->rcrit_part = rcom + 2*cbranch->maxrsph;
			cbranch->big_particle = cbranch->particles[0];
		}
		cbranch->rcrit_angle += cbranch->rcrit_part;

	}while(TreeNBWalkStep(tree,true));

	return;
}

/// simple rotates the coordinates in the xp array
void ForceTree::rotate_coordinates(double **coord){
	IndexType i;
	short j;
	PosType tmp[3];

  /* rotate particle positions */
  for(i=0;i<tree->top->nparticles;++i){
    for(j=0;j<3;++j) tmp[j]=0.0;
    for(j=0;j<3;++j){
      tmp[0]+=coord[0][j]*xp[i][j];
      tmp[1]+=coord[1][j]*xp[i][j];
      tmp[2]+=coord[2][j]*xp[i][j];
    }
    for(j=0;j<3;++j) xp[i][j]=tmp[j];
  }

  return;
}

float * ForceTree::CalculateSPHsmoothing(int N){
	IndexType i;//,neighbors[N];
	IndexType *neighbors = new IndexType[N];

	for(i=0;i<tree->top->nparticles;++i){
		NearestNeighbors(xp[i],N,&(tree->rsph[i]),neighbors);
	}

	// recalculate cutoff scales for each branch to agree with new sizes
	CalcMoments();
	return tree->rsph;
}

/** TreeNBForce2D calculates the defection, convergence and shear using
 *   the plane-lens approximation with 3D SPH smoothing of the density
 *   rsph must be calculated before doing this with FindRSPH
 *   tangent[3] - the direction of light rays or orientation of the simulation
 *   tree can be either a 3d or 2d tree although 2d is more efficient
 *       need to change the projected cm in _TreeNBForce to us 3d tree
 *
 *   negative rsph give a screened potential
 * */

void ForceTree::force2D(double *ray,double *alpha,double *kappa,double *gamma,bool no_kappa){

  PosType xcm,ycm,rcm2,tmp;
  int OpenBox(TreeNBHndl tree,PosType r);
  IndexType i;
  bool allowDescent=true;
  unsigned long count=0;

  assert(tree);
  moveTopNB(tree);

  alpha[0]=alpha[1]=gamma[0]=gamma[1]=0.0;
  *kappa=0.0;

  do{

	  ++count;
	  xcm=tree->current->center[0]-ray[0];
	  ycm=tree->current->center[1]-ray[1];

	  rcm2 = xcm*xcm + ycm*ycm;

	  /* if the box is close enough that the smoothing scale could be important
	   * add those particles individually and subtract their point particle contribution
	   * that will be added later
	   */

	  if( rcm2 < pow(tree->current->rcrit_angle,2) ){
		  // includes rcrit_particle constraint
		  allowDescent=true;

		  //printf("right place\n");
		  if( atLeaf() ){
			  // leaf case
			  // particle treatment

			  for(i=0;i<tree->current->nparticles;++i){

				  xcm = tree->xp[tree->current->particles[i]][0] - ray[0];
				  ycm = tree->xp[tree->current->particles[i]][1] - ray[1];

				  rcm2 = xcm*xcm + ycm*ycm;

				  tmp =  alpha_internal(rcm2,tree->rsph[tree->MultiRadius*tree->current->particles[i]])
						  *tree->masses[tree->MultiMass*tree->current->particles[i]];
				  alpha[0] += tmp*xcm;
				  alpha[1] += tmp*ycm;

				  // can turn off kappa and gamma calculations to save times
				  if(!no_kappa){
					  *kappa+=kappa_internal(rcm2,tree->rsph[tree->MultiRadius*tree->current->particles[i]])
		                  *tree->masses[tree->MultiMass*tree->current->particles[i]];
					  tmp= gamma_internal(rcm2,tree->rsph[tree->MultiRadius*tree->current->particles[i]])
						  *tree->masses[tree->MultiMass*tree->current->particles[i]];
					  gamma[0] += 0.5*(xcm*xcm-ycm*ycm)*tmp;
					  gamma[1] += xcm*ycm*tmp;
				  }
			  }
		  }
	  }else{ // use whole cell
		  allowDescent=false;

		  tmp = -1.0*tree->current->mass/rcm2/pi;

		  alpha[0] += tmp*xcm;
		  alpha[1] += tmp*ycm;

		  if(!no_kappa){      //  taken out to speed up
			  tmp=-2.0*tree->current->mass/pi/rcm2/rcm2;
			  gamma[0] += 0.5*(xcm*xcm-ycm*ycm)*tmp;
			  gamma[1] += xcm*ycm*tmp;
		  }

		  // quadrapole contribution
		  //   the kappa and gamma are not calculated to this order
		  alpha[0] -= (tree->current->quad[0]*xcm + tree->current->quad[2]*ycm)
    				  /pow(rcm2,2)/pi;
		  alpha[1] -= (tree->current->quad[1]*ycm + tree->current->quad[2]*xcm)
    				  /pow(rcm2,2)/pi;

		  tmp = 4*(tree->current->quad[0]*xcm*xcm + tree->current->quad[1]*ycm*ycm
				  + 2*tree->current->quad[2]*xcm*ycm)/pow(rcm2,3)/pi;

		  alpha[0] += tmp*xcm;
		  alpha[1] += tmp*ycm;
	  }

  }while(TreeNBWalkStep(tree,allowDescent));

  return;
}

void ForceTree::ChangeParticleProfile(PartProf partprof){
	switch(partprof){
	case gaussian:
		alpha_internal = alpha_o;
		kappa_internal = kappa_o;
		gamma_internal = gamma_o;
		break;
	default:
		alpha_internal = alpha_o;
		kappa_internal = kappa_o;
		gamma_internal = gamma_o;
		break;
	}
}
