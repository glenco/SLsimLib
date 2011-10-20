#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "TreeNB.h"

/* find the SPH smoothing length for each particle */
float *FindRSPH(TreeNBHndl tree,int Nsph){
  float *rsph;
  IndexType *neighbors,i;

  if(sizeof(tree->xp[0]) != sizeof(double)){
	  ERROR_MESSAGE();
	  printf("ERROR: FindRSPH, xp is not a double, need to change PosType\n");
	  exit(0);
  }
  rsph=(float *)malloc((tree->top->nparticles)*sizeof(float));

  for(i=0;i<tree->top->nparticles;++i){
	if(i<0){ERROR_MESSAGE(); printf("ERROR: findRSPH, not enough bits to address particle number\n"); exit(1);}
    if(Nsph > 0) neighbors=NearestNeighborNB(tree,tree->xp[i],Nsph,&(tree->rsph[i]));
    else tree->rsph[i]=0.0;
    if(i % 50000 == 0) printf("         %.2f percent done\n",100.*i/tree->top->nparticles);
  }

  return rsph;
}

/* TreeNBForce2D calculates the defection, convergence and shear using
 *   the plane-lens approximation with 3D SPH smoothing of the density
 *   rsph must be calculated before doing this with FindRSPH
 *   tangent[3] - the direction of light rays or orientation of the simulation
 *   tree can be either a 3d or 2d tree although 2d is more efficient
 *       need to change the projected cm in _TreeNBForce to us 3d tree
 *
 *   negative rsph give a screened potential
 * */

void TreeNBForce2D(TreeNBHndl tree,double *ray
		,double *alpha,double *kappa,double *gamma,bool no_kappa
		   ){

	assert(tree);

  alpha[0]=0.0; alpha[1]=0.0;
  *kappa=0.0; gamma[0]=0.0; gamma[1]=0.0;

  TreeNBParticleForce2Diter(tree,ray,alpha,kappa,gamma,no_kappa,
		  alpha_o,kappa_o,gamma_o);
  //moveTopNB(tree);
//   _TreeNBParticleForce2D(tree,xp_2d,rsph,ray,alpha,kappa,gamma,no_kappa);
  //_TreeNBForce2D(tree,xp_2d,rsph,ray,alpha,kappa,gamma,no_kappa);

  return ;
}

void _TreeNBForce2D(TreeNBHndl tree,double *ray
		  ,double *alpha,double *kappa,double *gamma,bool no_kappa){
  PosType xcm,ycm,rcm2,tmp;
  int OpenBox(TreeNBHndl tree,PosType r);
  int i;

/*   xcm=tree->current->cm[0]-ray[0]; */
/*   ycm=tree->current->cm[1]-ray[1]; */

  xcm=tree->current->center[0]-ray[0];
  ycm=tree->current->center[1]-ray[1];

  rcm2 = xcm*xcm + ycm*ycm;
  //rcrit2=pow(tree->current->rcrit_angle+tree->current->rcrit_part,2);
  //printf("rsm=%e rcrit=%e\n",sqrt(rcm2),tree->current->rcrit_angle);

  if( rcm2 < pow(tree->current->rcrit_angle,2) ){
	  if( (tree->current->child1==NULL)*(tree->current->child2==NULL) ){ /* leaf case */
	  // particle treatment

		  for(i=0;i<tree->current->nparticles;++i){

			  xcm=tree->xp[tree->current->particles[i]][0]-ray[0];
			  ycm=tree->xp[tree->current->particles[i]][1]-ray[1];

			  rcm2=xcm*xcm + ycm*ycm;

			  tmp=alpha_o(rcm2,tree->rsph[tree->MultiRadius*tree->current->particles[i]]);
			  alpha[0]+=tmp*xcm;
			  alpha[1]+=tmp*ycm;

			  // can turn off kappa and gamma calculations to save times
			  if(!no_kappa){
				  *kappa+=kappa_o(rcm2,tree->rsph[tree->MultiRadius*tree->current->particles[i]]);
				  tmp=gamma_o(rcm2,tree->rsph[tree->MultiRadius*tree->current->particles[i]]);
				  gamma[0]+=0.5*(xcm*xcm-ycm*ycm)*tmp;
				  gamma[1]+=xcm*ycm*tmp;
			  }
		  }

	  }else{

		  if(tree->current->child1 != NULL){
			  moveToChildNB(tree,1);
			  _TreeNBForce2D(tree,ray,alpha,kappa,gamma,no_kappa);
			  moveUpNB(tree);
		  }
		  if(tree->current->child2 != NULL){
			  moveToChildNB(tree,2);
			  _TreeNBForce2D(tree,ray,alpha,kappa,gamma,no_kappa);
			  moveUpNB(tree);
		  }

	  }

  }else{ /* use whole cell  */

      tmp=-1.0*tree->current->nparticles/rcm2/pi;
      alpha[0]+=tmp*xcm;
      alpha[1]+=tmp*ycm;

      if(!no_kappa){      //  taken out to speed up
    	  //*kappa+=0;
    	  tmp=-2.0*tree->current->nparticles/pi/rcm2/rcm2;
    	  gamma[0]+=0.5*(xcm*xcm-ycm*ycm)*tmp;
    	  gamma[1]+=xcm*ycm*tmp;
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

      /*
      if(isnan(*kappa)){
    	  printf(" kappa=nan from cell tmp=%e xcm=%e ycm=%e \nrsph=%e  rcrit=%e  level=%i\n",tmp,xcm,ycm
    			  ,tree->rsph[tree->MultiRadius*tree->current->particles[0]],tree->current->rcrit,tree->current->level);
    	  printBranchNB(tree->current,xp_2d);
    	  exit(1);
      }
      */
  }
}

void _TreeNBParticleForce2D(TreeNBHndl tree,double *ray
		  ,double *alpha,double *kappa,double *gamma,bool no_kappa){
  PosType xcm,ycm,rcm2,tmp,rsphmax;
  int OpenBox(TreeNBHndl tree,PosType r);
  int i;

/*   xcm=tree->current->cm[0]-ray[0]; */
/*   ycm=tree->current->cm[1]-ray[1]; */

  xcm=tree->current->center[0]-ray[0];
  ycm=tree->current->center[1]-ray[1];

  rcm2 = xcm*xcm + ycm*ycm;
  rsphmax=tree->rsph[tree->MultiRadius*tree->current->big_particle];

  /* if the box is close enough that the smoothing scale could be important
   * add those particles individually and subtract their point particle contribution
   * that will be added later
   */

  if( rcm2 < pow(tree->current->rcrit_angle,2) ){
	  // includes rcrit_particle constraint

	  if( (tree->current->child1==NULL)*(tree->current->child2==NULL) ){
	   // leaf case
	  // particle treatment

		  for(i=0;i<tree->current->nparticles;++i){

			  xcm=tree->xp[tree->current->particles[i]][0]-ray[0];
			  ycm=tree->xp[tree->current->particles[i]][1]-ray[1];

			  rcm2=xcm*xcm + ycm*ycm;

			  tmp=alpha_o(rcm2,tree->rsph[tree->MultiRadius*tree->current->particles[i]]);
			  alpha[0]+=tmp*xcm;
			  alpha[1]+=tmp*ycm;

			  // can turn off kappa and gamma calculations to save times
			  if(!no_kappa){
				  *kappa+=kappa_o(rcm2,tree->rsph[tree->MultiRadius*tree->current->particles[i]]);
				  tmp=gamma_o(rcm2,tree->rsph[tree->MultiRadius*tree->current->particles[i]]);
				  gamma[0]+=0.5*(xcm*xcm-ycm*ycm)*tmp;
				  gamma[1]+=xcm*ycm*tmp;
			  }
		  }

	  }else{

		  // case where largest particle needs to be added singly
		  //   and a point mass needs to be subtracted
		  if(rcm2 < tree->current->rcrit_part*tree->current->rcrit_part){
			  ERROR_MESSAGE();
			  assert(tree->current->big_particle != -1);
			  xcm=tree->xp[tree->current->big_particle][0]-ray[0];
			  ycm=tree->xp[tree->current->big_particle][1]-ray[1];

			  rcm2=xcm*xcm + ycm*ycm;
			  tmp=alpha_o(rcm2,tree->current->maxrsph) - alpha_o(rcm2,0);
			  alpha[0]+=xcm*tmp;
			  alpha[1]+=ycm*tmp;

			  // can turn off kappa and gamma calculations to save times
			  if(!no_kappa){
				  *kappa+=kappa_o(rcm2,tree->current->maxrsph)
						  - kappa_o(rcm2,0);
				  tmp=gamma_o(rcm2,tree->current->maxrsph)
						  - gamma_o(rcm2,0);
				  gamma[0]+=0.5*(xcm*xcm-ycm*ycm)*tmp;
				  gamma[1]+=xcm*ycm*tmp;
			  }
		  }

		  if(tree->current->child1 != NULL){
			  moveToChildNB(tree,1);
			  _TreeNBForce2D(tree,ray,alpha,kappa,gamma,no_kappa);
			  moveUpNB(tree);
		  }
		  if(tree->current->child2 != NULL){
			  moveToChildNB(tree,2);
			  _TreeNBForce2D(tree,ray,alpha,kappa,gamma,no_kappa);
			  moveUpNB(tree);
		  }
	  }

  }else{ /* use whole cell  */

	  if(rsphmax<0) tmp=tree->current->nparticles*alpha_o(rcm2,rsphmax);  //case of screened potential
	  else tmp=-1.0*tree->current->nparticles/rcm2/pi;

      alpha[0]+=tmp*xcm;
      alpha[1]+=tmp*ycm;

      if(!no_kappa){      //  taken out to speed up
    	  if(rsphmax<0){
    		  //case of screened potential
    		  *kappa+=tree->current->nparticles*kappa_o(rcm2,rsphmax);
    		  tmp=tree->current->nparticles*gamma_o(rcm2,rsphmax);
    	  }else tmp=-2.0*tree->current->nparticles/pi/rcm2/rcm2;
    	  gamma[0]+=0.5*(xcm*xcm-ycm*ycm)*tmp;
    	  gamma[1]+=xcm*ycm*tmp;
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
}

void TreeNBParticleForce2Diter(TreeNBHndl tree
		  ,double *ray
		  ,double *alpha,double *kappa,double *gamma,bool no_kappa
		  ,double (*alpha_internal)(double r,float rmax)
		  ,double (*kappa_internal)(double r,float rmax)
		  ,double (*gamma_internal)(double r,float rmax)
){

  PosType xcm,ycm,rcm2,tmp,rsphmax;
  int OpenBox(TreeNBHndl tree,PosType r);
  int i;
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
	  rsphmax=tree->current->maxrsph;

	  /* if the box is close enough that the smoothing scale could be important
	   * add those particles individually and subtract their point particle contribution
	   * that will be added later
	   */

	  //printf("rsm %e rcrit=%e\n",sqrt(rcm2),tree->current->rcrit_angle);
	  //printf("ray = %e %e\n",ray[0],ray[1]);
	  //printBranchNB(tree->current,xp_2d,2);

	  if( rcm2 < pow(tree->current->rcrit_angle,2) ){
		  // includes rcrit_particle constraint
		  allowDescent=true;

		  //printf("right place\n");
		  if( (tree->current->child1==NULL)*(tree->current->child2==NULL) ){
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
		  }/*else{

			  // case where largest particle needs to be added singly
			  //   and a point mass needs to be subtracted
			  if(rcm2 < tree->current->rcrit_part*tree->current->rcrit_part){
				  //ERROR_MESSAGE();
				  assert(tree->current->big_particle != -1);

				  xcm = tree->xp[tree->current->big_particle][0] - ray[0];
				  ycm = tree->xp[tree->current->big_particle][1] - ray[1];

				  rcm2 = xcm*xcm + ycm*ycm;
				  tmp =  (alpha_internal(rcm2,tree->current->maxrsph)
						  + tree->current->mass/rcm2/pi)
						  *tree->mass[tree->MultiMass*tree->current->big_particle];
				  alpha[0] += xcm*tmp;
				  alpha[1] += ycm*tmp;

				  // can turn off kappa and gamma calculations to save times
				  if(!no_kappa){
					  *kappa+=kappa_internal(rcm2,tree->current->maxrsph)
						    *tree->mass[tree->MultiMass*tree->current->big_particle];
					  tmp= (gamma_internal(rcm2,tree->current->maxrsph)
										  + 2.0*tree->current->mass/pi/rcm2/rcm2 )
							*tree->mass[tree->MultiMass*tree->current->big_particle];
					  gamma[0] += 0.5*(xcm*xcm-ycm*ycm)*tmp;
					  gamma[1] += xcm*ycm*tmp;
				  }
			  }
		  }*/
	  }else{ // use whole cell
		  allowDescent=false;

		  tmp = -1.0*tree->current->mass/rcm2/pi;
		  //if(rsphmax <  0) tmp*=exp(-rcm2/rsphmax/rsphmax/2);
		  alpha[0] += tmp*xcm;
		  alpha[1] += tmp*ycm;

		  if(!no_kappa){      //  taken out to speed up
			  /*if(rsphmax<0){
				  //case of screened potential
				  *kappa+=tree->current->mass*kappa_internal(rcm2,rsphmax);
				  tmp=tree->current->mass* gamma_internal(rcm2,rsphmax);
			  }else */
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
