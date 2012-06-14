/*
 * forceTree.cpp
 *
 *  Created on: Oct 18, 2011
 *      Author: bmetcalf
 */
#include <assert.h>
#include <math.h>
#include "forceTree.h"

//const float pi = 3.141593;

ForceTree::ForceTree(
		PosType **xp
		,IndexType Npoints
		,float *Masses
		,float *Rsphs
		,bool Multimass
		,bool Multisize
		,double my_kappa_background /// background kappa that is subtracted
		,int bucket
		,int dimension
		,bool median
		,double theta
		) :
	SimpleTree(xp,Npoints,bucket,dimension,median)
	, MultiMass(Multimass), MultiRadius(Multisize), masses(Masses), rsph(Rsphs), force_theta(theta), kappa_background(my_kappa_background)
{
/*	tree->MultiMass = Multimass;
	tree->MultiRadius = Multisize;
	tree->masses = masses;
	tree->rsph= rsph;
	force_theta = theta;*/

	haloON = false; // don't use internal halo parameters
	halo_params = NULL;
	init = false;


	// This should be changed so other particle profiles can be used.
	//alpha_particle = alpha_o;
	//kappa_particle = kappa_o;
	//gamma_particle = gamma_o;
}

ForceTree::~ForceTree(){
}


// calculates moments of the mass and the cutoff scale for each box
void ForceTree::CalcMoments(){

	//*** make compatable
	IndexType i;
	PosType rcom,xcm[2],xcut;
	BranchNB *cbranch;
	double tmp;

	moveTopNB(tree);
	do{
		cbranch=tree->current; /* pointer to current branch */

		cbranch->rmax = sqrt( pow(cbranch->boundary_p2[0]-cbranch->boundary_p1[0],2)
				+ pow(cbranch->boundary_p2[1]-cbranch->boundary_p1[1],2) );

		// calculate mass
		for(i=0,cbranch->mass=0;i<cbranch->nparticles;++i)
			cbranch->mass +=  haloON ? halo_params[cbranch->particles[i]*MultiRadius].mass : masses[cbranch->particles[i]*MultiMass];

		// calculate center of mass
		cbranch->center[0]=cbranch->center[1]=0;
		for(i=0;i<cbranch->nparticles;++i){
			tmp = haloON ? halo_params[cbranch->particles[i]*MultiRadius].mass : masses[cbranch->particles[i]*MultiMass];
			cbranch->center[0] += tmp*tree->xp[cbranch->particles[i]][0]/cbranch->mass;
			cbranch->center[1] += tmp*tree->xp[cbranch->particles[i]][1]/cbranch->mass;
		}
		//////////////////////////////////////////////
		// calculate quadropole moment of branch
		//////////////////////////////////////////////
		cbranch->quad[0]=cbranch->quad[1]=cbranch->quad[2]=0;
		for(i=0;i<cbranch->nparticles;++i){
			xcm[0]=tree->xp[cbranch->particles[i]][0]-cbranch->center[0];
			xcm[1]=tree->xp[cbranch->particles[i]][1]-cbranch->center[1];
			xcut=pow(xcm[0],2) + pow(xcm[1],2);
			tmp = haloON ? halo_params[cbranch->particles[i]*MultiRadius].mass : masses[cbranch->particles[i]*MultiMass];

			cbranch->quad[0] += (xcut-2*xcm[0]*xcm[0])*tmp;
			cbranch->quad[1] += (xcut-2*xcm[1]*xcm[1])*tmp;
			cbranch->quad[2] += -2*xcm[0]*xcm[1]*tmp;
		}

		// largest distance from center of mass of cell
		for(i=0,rcom=0.0;i<2;++i) rcom += ( pow(cbranch->center[i]-cbranch->boundary_p1[i],2) > pow(cbranch->center[i]-cbranch->boundary_p2[i],2) ) ?
		pow(cbranch->center[i]-cbranch->boundary_p1[i],2) : pow(cbranch->center[i]-cbranch->boundary_p2[i],2);

		rcom=sqrt(rcom);

		if(force_theta > 0.0) cbranch->rcrit_angle = 1.15470*rcom/(force_theta);
		else  cbranch->rcrit_angle=1.0e100;


		if(MultiRadius){

			for(i=0,cbranch->maxrsph=0.0;i<cbranch->nparticles;++i){
				tmp = haloON ? halo_params[cbranch->particles[i]*MultiRadius].Rmax : rsph[cbranch->particles[i]*MultiRadius];
				if(cbranch->maxrsph <= tmp ){
					cbranch->maxrsph = tmp;
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

				if(maxrsph <= rsph[MultiRadius*cbranch->particles[i]] ){
					maxrsph = rsph[MultiRadius*cbranch->particles[i]];
					cbranch->big_particle = cbranch->particles[i];
				}
			}
			// if it is possible for the largest particle size to be big enough
			//    cause the opening of a box then tag it
			cbranch->maxrsph = maxrsph;

			if(maxrsph > cbranch->rcrit_angle/2){
				rsph[MultiRadius*cbranch->big_particle] = 0.0; // zero out so it wont show up in leaf
				cbranch->rcrit_part = rcom + 2*maxrsph;
			}else{
				cbranch->rcrit_part = 0;
			}
			 */

		}else{
			// single size case
			cbranch->maxrsph = haloON ? halo_params[0].Rmax : rsph[0];
			cbranch->rcrit_part = rcom + 2*cbranch->maxrsph;
			cbranch->big_particle = cbranch->particles[0];
		}
		cbranch->rcrit_angle += cbranch->rcrit_part;

	}while(TreeNBWalkStep(tree,true));

	init = true;
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
		NearestNeighbors(xp[i],N,&(rsph[i]),neighbors);
	}

	// recalculate cutoff scales for each branch to agree with new sizes
	CalcMoments();
	return rsph;
}

/** TreeNBForce2D calculates the defection, convergence and shear using
 *   the plane-lens approximation with 3D SPH smoothing of the density
 *   rsph must be calculated before doing this with FindRSPH or by other means.
 *   tangent[3] - the direction of light rays or orientation of the simulation
 *   tree can be either a 3d or 2d tree although 2d is more efficient
 *       need to change the projected cm in _TreeNBForce to us 3d tree
 *
 *       The output alpha[] is in units of mass_scale/Mpc, ie it needs to be
 *       divided by Sigma_crit and multiplied by mass_scale to be the defelction
 *       in the lens equation expressed on the lens plane or multiplied by
 *       4*pi*G*mass_scale to get the deflection angle caused by the plane lens.
 * */

void ForceTree::force2D(double *ray,double *alpha,double *kappa,double *gamma,bool no_kappa){

  PosType xcm,ycm,rcm2,tmp;
  int OpenBox(TreeNBHndl tree,PosType r);
  IndexType i;
  bool allowDescent=true;
  unsigned long count=0,index;

  assert(tree);
  moveTopNB(tree);

  // This ensures the moments are calculated even if they are not calculated on construction
  if(!init) CalcMoments();

  alpha[0]=alpha[1]=gamma[0]=gamma[1]=gamma[2]=0.0;
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

				  index = MultiRadius*tree->current->particles[i];
				  if(haloON) tmp = alpha_h(rcm2,halo_params[index]);
				  else tmp =  alpha_o(rcm2,rsph[index])*masses[MultiMass*tree->current->particles[i]];
				  alpha[0] += tmp*xcm;
				  alpha[1] += tmp*ycm;

				  // can turn off kappa and gamma calculations to save times
				  if(!no_kappa){
					  if(haloON) *kappa += kappa_h(rcm2,halo_params[index]);
					  else *kappa += kappa_o(rcm2,rsph[index])*masses[MultiMass*tree->current->particles[i]];

					  if(haloON) tmp = gamma_h(rcm2,halo_params[index]);
					  else tmp = gamma_o(rcm2,rsph[index])*masses[MultiMass*tree->current->particles[i]];
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

  // Subtract off uniform mass sheet to compensate for the extra mass
  //  added to the universe in the halos.
  alpha[0] += ray[0]*kappa_background;
  alpha[1] += ray[1]*kappa_background;
  if(!no_kappa){      //  taken out to speed up
	  *kappa -= kappa_background;
  }

  return;
}
/*
void ForceTree::ChangeParticleProfile(PartProf partprof){
	switch(partprof){
	case gaussian:
		alpha_particle = alpha_o;
		kappa_particle = kappa_o;
		gamma_particle = gamma_o;
		break;
	default:
		alpha_particle = alpha_o;
		kappa_particle = kappa_o;
		gamma_particle = gamma_o;
		break;
	}
}
*/
/** \ingroup DeflectionL2
\brief !!!! Not yet come up with a good way of doing this !!!!

	 * gives the particles a third dimension depending on their size
	 *  This is used to spread subhalos out into 3d so that their density
	 *  is about one per smoothing length.  This is to avoid the big particle
	 *  problem that slows the tree-force calculation.
	 *
	 *  Warning: This will erase the third coordinate of the particles.
	 */
void ForceTree::spread_particles(){

	IndexType i;
	TreeNBHndl tree;
	float tmp=0;

	for(i=0;i<Nparticles;++i){
		//TODO: this needs to be fixed!!!
		//dummy = NearestNeighborNB(tree,xp[i],1,&tmp);
		tree->xp[i][2] = 4*pi*pow(rsph[i*MultiRadius],3)/pow(tmp,2)/3;
	}

	return;
}

ForceTreePowerLaw::ForceTreePowerLaw(
		float beta                 /// slop of mass profile \f$ \Sigma \propto r^\beta \f$
		,PosType **xp              /// positions of the halos xp[0..Npoints-1][0..1 or 2]
		,IndexType Npoints         /// number of halos
		,HaloStructure *h_params   /// array with internal properties of halos
		,bool Multisize            /// flag false if only one halo size and structure should be used, default is true
		,double my_kappa_bk        /// Background convergence to be subtracted
		,int bucket                /// maximum number of halos in a leaf of the tree
		,int dimension             /// 2 or 3, dimension of tree, default 2
		,bool median               /// If true will divide branches at the median position of the particles, if false an equal area cut is used, default false
		,PosType theta             /// Opening angle used in tree force calculation
		) :
		ForceTree(xp,Npoints,NULL,NULL,false,Multisize,my_kappa_bk,bucket,dimension,median,theta), beta(beta)
{

	haloON = true;
	halo_params = h_params;

	CalcMoments();
}

ForceTreePowerLaw::~ForceTreePowerLaw(){
}

ForceTreeNFW::ForceTreeNFW(
		PosType **xp               /// positions of the halos xp[0..Npoints-1][0..1 or 2]
		,IndexType Npoints         /// number of halos
		,HaloStructure *h_params   /// array with internal properties of halos
		,bool Multisize            /// flag false if only one halo size and structure should be used, default is true
		,double my_kappa_bk        /// Background convergence to be subtracted
		,int bucket                /// maximum number of halos in a leaf of the tree
		,int dimension             /// 2 or 3, dimension of tree, default 2
		,bool median               /// If true will divide branches at the median position of the particles, if false an equal area cut is used, default false
		,PosType theta             /// Opening angle used in tree force calculation
		) :
		ForceTree(xp,Npoints,NULL,NULL,false,Multisize,my_kappa_bk,bucket,dimension,median,theta)
{
	for(unsigned long i=0;i<Npoints;++i){
		if(h_params[i].Rmax <= 0.0 || h_params[i].rscale <= 0.0){
			ERROR_MESSAGE();
			cout << "Illegal values for halo internal valuables." << endl;
			exit(1);
		}
	}

	haloON = true;
	halo_params = h_params;

	CalcMoments();
}

ForceTreeNFW::~ForceTreeNFW(){
}

ForceTreePseudoNFW::ForceTreePseudoNFW(
		int my_beta                 /// outer slope of profile is \f$ \Sigma \propto r^{-\beta} \f$
		,PosType **xp              /// positions of the halos xp[0..Npoints-1][0..1 or 2]
		,IndexType Npoints         /// number of halos
		,HaloStructure *h_params   /// array with internal properties of halos
		,bool Multisize            /// flag false if only one halo size and structure should be used, default is true
		,double my_kappa_bk       /// Background convergence to be subtracted
		,int bucket                /// maximum number of halos in a leaf of the tree
		,int dimension             /// 2 or 3, dimension of tree, default 2
		,bool median               /// If true will divide branches at the median position of the particles, if false an equal area cut is used, default false
		,PosType theta             /// Opening angle used in tree force calculation
		) :
		ForceTree(xp,Npoints,NULL,NULL,false,Multisize,my_kappa_bk,bucket,dimension,median,theta), beta(my_beta)
{

	if(beta == 0){
		cout << "The slope can not be zero!" << endl;
		exit(1);
	}

	// Check for values that would make the rayshooter return nan.
	for(unsigned long i=0;i<Npoints;++i){
		if(h_params[i].Rmax <= 0.0 || h_params[i].rscale <= 0.0){
			ERROR_MESSAGE();
			cout << "Illegal values for halo internal valuables." << endl;
			exit(1);
		}
	}
	haloON = true;
	halo_params = h_params;

	CalcMoments();
}

ForceTreePseudoNFW::~ForceTreePseudoNFW(){
}
