/*
 * forceTree.cpp
 *
 *  Created on: Oct 18, 2011
 *      Author: bmetcalf
 */

#include "slsimlib.h"

TreeForce::TreeForce(
		PosType **xp
		,IndexType Npoints
		,float *Masses
		,float *Rsphs
		,bool Multimass
		,bool Multisize
		,PosType my_kappa_background /// background kappa that is subtracted
		,int bucket
		,int dimension
		,bool median
		,PosType theta
		) :
	TreeSimple(xp,Npoints,bucket,dimension,median)
	, MultiMass(Multimass), MultiRadius(Multisize), masses(Masses), rsph(Rsphs), kappa_background(my_kappa_background), force_theta(theta)
{

	haloON = false; // don't use internal halo parameters
	halos = NULL;
	init = false;

	CalcMoments();
}

TreeForce::TreeForce(
		PosType **xp              /// positions of the halos xp[0..Npoints-1][0..1 or 2]
		,IndexType Npoints         /// number of halos
		,LensHalo *my_halos   /// array with internal properties of halos
		,bool Multisize            /// flag false if only one halo size and structure should be used, default is true
		,PosType my_kappa_bk        /// Background convergence to be subtracted
		,int bucket                /// maximum number of halos in a leaf of the tree
		,int dimension             /// 2 or 3, dimension of tree, default 2
		,bool median               /// If true will divide branches at the median position of the particles, if false an equal area cut is used, default false
		,PosType theta             /// Opening angle used in tree force calculation
		) :
		TreeSimple(xp,Npoints,bucket,dimension,median)
, MultiMass(true), MultiRadius(true), masses(NULL), rsph(NULL), kappa_background(my_kappa_bk), force_theta(theta)
{

	haloON = true;
	halos = my_halos;

	CalcMoments();
}

TreeForce::~TreeForce(){
}


// calculates moments of the mass and the cutoff scale for each box
void TreeForce::CalcMoments(){

	//*** make compatable
	IndexType i;
	PosType rcom,xcm[2],xcut;
	BranchNB *cbranch;
	PosType tmp;

	moveTopNB(tree);
	do{
		cbranch=tree->current; /* pointer to current branch */

		cbranch->rmax = sqrt( pow(cbranch->boundary_p2[0]-cbranch->boundary_p1[0],2)
				+ pow(cbranch->boundary_p2[1]-cbranch->boundary_p1[1],2) );

		// calculate mass
		for(i=0,cbranch->mass=0;i<cbranch->nparticles;++i)
			cbranch->mass +=  haloON ? halos[cbranch->particles[i]*MultiRadius].get_mass() : masses[cbranch->particles[i]*MultiMass];

		// calculate center of mass
		cbranch->center[0]=cbranch->center[1]=0;
		for(i=0;i<cbranch->nparticles;++i){
			tmp = haloON ? halos[cbranch->particles[i]*MultiRadius].get_mass() : masses[cbranch->particles[i]*MultiMass];
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
			tmp = haloON ? halos[cbranch->particles[i]*MultiRadius].get_mass() : masses[cbranch->particles[i]*MultiMass];

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
				tmp = haloON ? halos[cbranch->particles[i]*MultiRadius].get_Rmax() : rsph[cbranch->particles[i]*MultiRadius];
				if(cbranch->maxrsph <= tmp ){
					cbranch->maxrsph = tmp;
				}
			}

			cbranch->rcrit_part = rcom + 2*cbranch->maxrsph;

		}else{
			// single size case
			cbranch->maxrsph = haloON ? halos[0].get_Rmax() : rsph[0];
			cbranch->rcrit_part = rcom + 2*cbranch->maxrsph;
			cbranch->big_particle = cbranch->particles[0];
		}
		cbranch->rcrit_angle += cbranch->rcrit_part;

	}while(TreeNBWalkStep(tree,true));

	init = true;
	return;
}

/// simple rotates the coordinates in the xp array
void TreeForce::rotate_coordinates(PosType **coord){
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

float * TreeForce::CalculateSPHsmoothing(int N){
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

void TreeForce::force2D(PosType const *ray
                        ,PosType *alpha
                        ,KappaType *kappa
                        ,KappaType *gamma
                        ,KappaType *phi
                        ,bool no_kappa)
{

  PosType xcm[2],rcm2,tmp;
  int OpenBox(TreeNBHndl tree,PosType r);
  IndexType i;
  bool allowDescent=true;
  unsigned long count=0,index;
  PosType rcm, arg1, arg2, prefac;

  assert(tree);
  moveTopNB(tree);

  // This ensures the moments are calculated even if they are not calculated on construction
  if(!init) CalcMoments();

  alpha[0]=alpha[1]=gamma[0]=gamma[1]=gamma[2]=0.0;
  *kappa=0.0;

  do{

	  ++count;
	  xcm[0]=tree->current->center[0]-ray[0];
	  xcm[1]=tree->current->center[1]-ray[1];

	  rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];

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

				  xcm[0] = tree->xp[tree->current->particles[i]][0] - ray[0];
				  xcm[1] = tree->xp[tree->current->particles[i]][1] - ray[1];

				  index = MultiRadius*tree->current->particles[i];

				  rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
				  if(rcm2 < 1e-20) rcm2 = 1e-20;
				  rcm = sqrt(rcm2);

				  prefac = haloON ? halos[index].get_mass()/rcm2/pi : masses[MultiMass*index]/rcm2/pi;

				  tmp = -1.0*prefac;

				  alpha[0] += tmp*xcm[0];
				  alpha[1] += tmp*xcm[1];

				  // can turn off kappa and gamma calculations to save times
				  if(!no_kappa){
					  tmp = -2.0*prefac/rcm2;

					  gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
					  gamma[1] += xcm[0]*xcm[1]*tmp;
				  }

				  if(haloON){
					  halos[index].force_halo(alpha,kappa,gamma,phi,xcm,no_kappa,true); // PHI BY Fabien
				  }else{  // case of no halos just particles and no class derived from TreeQuad

					  arg1 = rcm2/(rsph[index*MultiRadius]*rsph[index*MultiRadius]);
					  arg2 = rsph[index*MultiRadius];
					  tmp = rsph[index*MultiRadius];

					  /// intersecting, subtract the point particle
					  if(rcm2 < tmp*tmp){
						  tmp = (alpha_h(arg1,arg2)+1.0)*prefac;
						  alpha[0] += tmp*xcm[0];
						  alpha[1] += tmp*xcm[1];

						  // can turn off kappa and gamma calculations to save times
						  if(!no_kappa){
							  *kappa += kappa_h(arg1,arg2)*prefac;

							  tmp = (gamma_h(arg1,arg2)+2.0)*prefac/rcm2;

							  gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
							  gamma[1] += xcm[0]*xcm[1]*tmp;
						  }
					  }
				  }
			  }
		  }
	  }else{ // use whole cell
		  allowDescent=false;

		  tmp = -1.0*tree->current->mass/rcm2/pi;

		  alpha[0] += tmp*xcm[0];
		  alpha[1] += tmp*xcm[1];

		  if(!no_kappa){      //  taken out to speed up
			  tmp=-2.0*tree->current->mass/pi/rcm2/rcm2;
			  gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
			  gamma[1] += xcm[0]*xcm[1]*tmp;
		  }

		  // quadrapole contribution
		  //   the kappa and gamma are not calculated to this order
		  alpha[0] -= (tree->current->quad[0]*xcm[0] + tree->current->quad[2]*xcm[1])
    				  /pow(rcm2,2)/pi;
		  alpha[1] -= (tree->current->quad[1]*xcm[1] + tree->current->quad[2]*xcm[0])
    				  /pow(rcm2,2)/pi;

		  tmp = 2*(tree->current->quad[0]*xcm[0]*xcm[0] + tree->current->quad[1]*xcm[1]*xcm[1]
				  + 2*tree->current->quad[2]*xcm[0]*xcm[1])/pow(rcm2,3)/pi;

		  alpha[0] += tmp*xcm[0];
		  alpha[1] += tmp*xcm[1];
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

