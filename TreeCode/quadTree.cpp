/*
 * quadTree.cpp
 *
 *  Created on: Sep 4, 2012
 *      Author: mpetkova
 */

/*
 * Programmer:    R Ben Metcalf
 */

#include "slsimlib.h"

/** \brief Constructor meant for point particles, simulation particles
 */
QuadTree::QuadTree(
		PosType **xpt
		,float *my_masses
		,float *my_sizes
		,IndexType Npoints
		,bool Multimass
		,bool Multisize
		,double my_sigma_background /// background kappa that is subtracted
		,int bucket
		,double theta_force
		):
	xp(xpt),MultiMass(Multimass), MultiRadius(Multisize), masses(my_masses),sizes(my_sizes),Nparticles(Npoints),
	sigma_background(my_sigma_background),Nbucket(bucket),force_theta(theta_force)
{
	index = new IndexType[Npoints];
	IndexType ii;

	for(ii=0;ii<Npoints;++ii) index[ii] = ii;

	haloON = false; // don't use internal halo parameters
	halos = NULL;

	tree = BuildQTreeNB(xp,Npoints,index);

	CalcMoments();

	return;
}
/** \brief Constructor meant for halos with internal structure parameters.  This is a protected constructor because
 * it should only be invoked from the derived classes that have specific defined halo models.
 */
QuadTree::QuadTree(
		PosType **xpt
		,LensHaloHndl *my_halos
		,IndexType Npoints
		,double my_sigma_background /// background kappa that is subtracted
		,int bucket
		,double theta_force
		):
	xp(xpt),MultiMass(true),MultiRadius(true),masses(NULL),sizes(NULL),Nparticles(Npoints),
	sigma_background(my_sigma_background),Nbucket(bucket),force_theta(theta_force),halos(my_halos)
{
	index = new IndexType[Npoints];
	IndexType ii;

	for(ii=0;ii<Npoints;++ii) index[ii] = ii;

	haloON = true; //use internal halo parameters

	tree = BuildQTreeNB(xp,Npoints,index);

	CalcMoments();

	return;
}


QuadTree::~QuadTree()
{
	delete tree;
	delete[] index;
	return;
}

QTreeNBHndl QuadTree::BuildQTreeNB(PosType **xp,IndexType Nparticles,IndexType *particles){
  IndexType i;
  short j;
  PosType p1[2],p2[2];

  for(j=0;j<2;++j){
    p1[j]=xp[0][j];
    p2[j]=xp[0][j];
  }

  // Find borders that enclose all particles
  for(i=0;i<Nparticles;++i){
    for(j=0;j<2;++j){
      if(xp[i][j] < p1[j] ) p1[j]=xp[i][j];
      if(xp[i][j] > p2[j] ) p2[j]=xp[i][j];
    }
  }

  if(Nparticles == 1){
	  p1[0]=-0.25;
	  p1[1]=-0.25;
	  p2[0]=0.25;
	  p2[1]=0.25;
  }

  /* Initialize tree root */
  tree = new QTreeNB(xp,particles,Nparticles,p1,p2);

  /* build the tree */
  _BuildQTreeNB(Nparticles,particles);

  /* visit every branch to find center of mass and cutoff scale */
  tree->moveTop();

  return tree;
}

/// returns an index for which of the four quadrangles of the branch the point x[] is in
inline short QuadTree::WhichQuad(double *x,QBranchNB &branch){
	return (x[0] < branch.center[0]) + 2*(x[1] < branch.center[1]);
}

/// tree must be created and first branch must be set before start
void QuadTree::_BuildQTreeNB(IndexType nparticles,IndexType *particles){

	QBranchNB *cbranch = tree->current; /* pointer to current branch */
	IndexType i,j,cut,cut2,jt;

	cbranch->center[0] = (cbranch->boundary_p1[0] + cbranch->boundary_p2[0])/2;
	cbranch->center[1] = (cbranch->boundary_p1[1] + cbranch->boundary_p2[1])/2;
	cbranch->quad[0]=cbranch->quad[1]=cbranch->quad[2]=0;
	cbranch->mass = 0.0;


	// leaf case
	if(cbranch->nparticles <= Nbucket){
		PosType r;
		cbranch->Nbig_particles = 0;
		for(i=0;i<cbranch->nparticles;++i){
			jt = particles[i]*MultiRadius;
			if(haloON){
				r = halos[jt]->get_Rmax();
			}else{
				r = sizes[jt];
			}
			//r = haloON ? halos[particles[i]*MultiRadius].Rmax	: sizes[particles[i]*MultiRadius];
			if(r < (cbranch->boundary_p2[0]-cbranch->boundary_p1[0])) ++cbranch->Nbig_particles;
		}
		if(cbranch->Nbig_particles){
			cbranch->big_particles = new IndexType[cbranch->Nbig_particles];
			for(i=0,j=0;i<cbranch->nparticles;++i){
				jt = particles[i]*MultiRadius;
				if(haloON){
					r = halos[jt]->get_Rmax();
				}else{
					r = sizes[jt];
				}
				//r = haloON ? halos[particles[i]*MultiRadius].Rmax	: sizes[particles[i]*MultiRadius];
				if(r < (cbranch->boundary_p2[0]-cbranch->boundary_p1[0])) cbranch->big_particles[j++] = particles[i];
			}
		}
		else{
			cbranch->big_particles = NULL;
		}

		return;
	}

	// find particles too big to be in children

	double *x = new double[cbranch->nparticles];

	cbranch->Nbig_particles=0;

	if(MultiRadius){
		// store the particles that are too large to be in a child at the end of the list of particles in cbranch
		for(i=0;i<cbranch->nparticles;++i){
			if(haloON){
				x[i] = halos[particles[i]]->get_Rmax();
			}else{
				x[i] = sizes[particles[i]];
			}
			//x[i] =  haloON ? halos[particles[i]].Rmax : sizes[particles[i]];
		}
		double maxsize = (cbranch->boundary_p2[0]-cbranch->boundary_p1[0])/2;

		// sort particles in size
		//quicksort(particles,x,cbranch->nparticles);
		Utilities::quickPartition(maxsize,&cut,particles,x,cbranch->nparticles);

		if(cut < cbranch->nparticles){
			if(tree->atTop()){
				cbranch->Nbig_particles = cut2 = cbranch->nparticles-cut;
			}else{
				Utilities::quickPartition(2*maxsize,&cut2,&particles[cut],&x[cut],cbranch->nparticles-cut);
				cbranch->Nbig_particles = cut2;
			}

			cbranch->big_particles = new IndexType[cbranch->Nbig_particles];
			for(i=cut;i<(cut+cut2);++i) cbranch->big_particles[i-cut] = particles[i];
		}

	}else{
		if(haloON){
			x[0] = halos[0]->get_Rmax();
		}else{
			x[0] = sizes[0];
		}
		//x[0] =  haloON ? halos[0].Rmax : sizes[0];
		if(x[0] > (cbranch->boundary_p2[0]-cbranch->boundary_p1[0])/2
				&& x[0] < (cbranch->boundary_p2[0]-cbranch->boundary_p1[0])){
			cbranch->Nbig_particles = cbranch->nparticles;
			cbranch->big_particles = cbranch->particles;
		}
	}

	assert( cbranch->nparticles >= cbranch->Nbig_particles );
	IndexType NpInChildren = cbranch->nparticles;// - cbranch->Nbig_particles;
	assert(NpInChildren >= 0);

	if(NpInChildren == 0){
		delete[] x;
		return;
	}

	IndexType cutx,cuty;
	PosType xcut,ycut;

	QBranchNB *child0 = new QBranchNB();
	QBranchNB *child1 = new QBranchNB();
	QBranchNB *child2 = new QBranchNB();
	QBranchNB *child3 = new QBranchNB();

	tree->attachChildrenToCurrent(child0,child1,child2,child3);

	// initialize boundaries of children

	child0->boundary_p1[0]=cbranch->center[0];
	child0->boundary_p1[1]=cbranch->center[1];
	child0->boundary_p2[0]=cbranch->boundary_p2[0];
	child0->boundary_p2[1]=cbranch->boundary_p2[1];

	child1->boundary_p1[0]=cbranch->boundary_p1[0];
	child1->boundary_p1[1]=cbranch->center[1];
	child1->boundary_p2[0]=cbranch->center[0];
	child1->boundary_p2[1]=cbranch->boundary_p2[1];

	child2->boundary_p1[0]=cbranch->center[0];
	child2->boundary_p1[1]=cbranch->boundary_p1[1];
	child2->boundary_p2[0]=cbranch->boundary_p2[0];
	child2->boundary_p2[1]=cbranch->center[1];

	child3->boundary_p1[0]=cbranch->boundary_p1[0];
	child3->boundary_p1[1]=cbranch->boundary_p1[1];
	child3->boundary_p2[0]=cbranch->center[0];
	child3->boundary_p2[1]=cbranch->center[1];


	// **** makes sure force does not require nbucket at leaf

	xcut = cbranch->center[0];
	ycut = cbranch->center[1];

	// divide along y direction
	for(i=0;i<NpInChildren;++i) x[i] = tree->xp[particles[i]][1];
	Utilities::quickPartition(ycut,&cuty,particles,x,NpInChildren);

	if(cuty > 0){
	  // divide first group in the x direction
      for(i=0;i<cuty;++i) x[i] = tree->xp[particles[i]][0];
      Utilities::quickPartition(xcut,&cutx,particles,x,cuty);

      child3->nparticles=cutx;
      assert(child3->nparticles >= 0);
      if(child3->nparticles > 0)
    	  child3->particles = particles;
      else child3->particles = NULL;

      child2->nparticles = cuty - cutx;
      assert(child2->nparticles >= 0);
      if(child2->nparticles > 0)
    	  child2->particles = &particles[cutx];
      else child2->particles = NULL;

	}else{
	  child3->nparticles = 0;
	  child3->particles = NULL;

	  child2->nparticles = 0;
	  child2->particles = NULL;
	}

	if(cuty < NpInChildren){
	  // divide second group in the x direction
	  for(i=cuty;i<NpInChildren;++i) x[i-cuty] = tree->xp[particles[i]][0];
	  Utilities::quickPartition(xcut,&cutx,&particles[cuty],x,NpInChildren-cuty);

	  child1->nparticles=cutx;
	  assert(child1->nparticles >= 0);
	  if(child1->nparticles > 0)
		  child1->particles = &particles[cuty];
	  else child1->particles = NULL;

	  child0->nparticles=NpInChildren - cuty - cutx;
	  assert(child0->nparticles >= 0);
	  if(child0->nparticles > 0)
		  child0->particles = &particles[cuty + cutx];
	  else child0->particles = NULL;

	}else{
	  child1->nparticles = 0;
	  child1->particles = NULL;

	  child0->nparticles = 0;
	  child0->particles = NULL;
	}

	delete[] x;

	tree->moveToChild(0);
	_BuildQTreeNB(child0->nparticles,child0->particles);
	tree->moveUp();

	tree->moveToChild(1);
	_BuildQTreeNB(child1->nparticles,child1->particles);
	tree->moveUp();

	tree->moveToChild(2);
	_BuildQTreeNB(child2->nparticles,child2->particles);
	tree->moveUp();

	tree->moveToChild(3);
	_BuildQTreeNB(child3->nparticles,child3->particles);
	tree->moveUp();

	return;
}

// calculates moments of the mass and the cutoff scale for each box
void QuadTree::CalcMoments(){

	//*** make compatable
	IndexType i;
	PosType rcom,xcm[2],xcut;
	QBranchNB *cbranch;
	double tmp;

	tree->moveTop();
	do{
		cbranch=tree->current; /* pointer to current branch */

		cbranch->rmax = sqrt( pow(cbranch->boundary_p2[0]-cbranch->boundary_p1[0],2)
				+ pow(cbranch->boundary_p2[1]-cbranch->boundary_p1[1],2) );

		// calculate mass
		for(i=0,cbranch->mass=0;i<cbranch->nparticles;++i)
			if(haloON ){
				cbranch->mass +=  halos[cbranch->particles[i]]->get_mass();
			}else{
				cbranch->mass += masses[cbranch->particles[i]*MultiMass];
			}
			//cbranch->mass +=  haloON ? halos[cbranch->particles[i]*MultiMass].mass : masses[cbranch->particles[i]*MultiMass];

		// calculate center of mass
		cbranch->center[0]=cbranch->center[1]=0;
		for(i=0;i<cbranch->nparticles;++i){
			if(haloON ){
				tmp = halos[cbranch->particles[i]]->get_mass();
			}else{
				tmp = masses[cbranch->particles[i]*MultiMass];
			}
			//tmp = haloON ? halos[cbranch->particles[i]*MultiMass].mass : masses[cbranch->particles[i]*MultiMass];
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
			xcut = pow(xcm[0],2) + pow(xcm[1],2);
			if(haloON ){
				tmp = halos[cbranch->particles[i]]->get_mass();
			}else{
				tmp = masses[cbranch->particles[i]*MultiMass];
			}
			//tmp = haloON ? halos[cbranch->particles[i]*MultiMass].mass : masses[cbranch->particles[i]*MultiMass];

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

		/*if(MultiRadius){

			for(i=0,cbranch->maxrsph=0.0;i<cbranch->nparticles;++i){
				tmp = haloON ? halos[cbranch->particles[i]].Rmax : sizes[cbranch->particles[i]];
				if(cbranch->maxrsph <= tmp ){
					cbranch->maxrsph = tmp;
				}
			}

			cbranch->rcrit_part = rcom + 2*cbranch->maxrsph;

		}else{
			// single size case
			cbranch->maxrsph = haloON ? halos[0].Rmax : sizes[0];
			cbranch->rcrit_part = rcom + 2*cbranch->maxrsph;
		}*/
		//cbranch->rcrit_angle += cbranch->rcrit_part;

	}while(tree->WalkStep(true));

	return;
}

/// simple rotates the coordinates in the xp array
void QuadTree::rotate_coordinates(double **coord){
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

/** \brief Force2D calculates the defection, convergence and shear using
 *   the plane-lens approximation.
 *
 *       The output alpha[] is in units of mass_scale/Mpc, ie it needs to be
 *       divided by Sigma_crit and multiplied by mass_scale to be the deflection
 *       in the lens equation expressed on the lens plane or multiplied by
 *       4*pi*G*mass_scale to get the deflection angle caused by the plane lens.
 *
 *       kappa and gamma need to by multiplied by mass_scale/Sigma_crit to get
 *       the traditional units for these where Sigma_crit are in the mass/units(ray)^2
 *       NB : the units of sigma_backgound need to be mass/units(ray)^2
 * */

void QuadTree::force2D(double *ray,double *alpha,KappaType *kappa,KappaType *gamma,bool no_kappa){

  PosType xcm[2],rcm2cell,rcm2,tmp,boxsize2;
  IndexType i;
  bool allowDescent=true;
  unsigned long count=0,index;
  double rcm, arg1, arg2, prefac;

  assert(tree);
  tree->moveTop();

  alpha[0]=alpha[1]=gamma[0]=gamma[1]=gamma[2]=0.0;

  *kappa=0.0;

  do{
	  ++count;
	  allowDescent=false;
	  if(tree->current->nparticles > 0){

		  xcm[0]=tree->current->center[0]-ray[0];
		  xcm[1]=tree->current->center[1]-ray[1];

		  rcm2cell = xcm[0]*xcm[0] + xcm[1]*xcm[1];

		  boxsize2 = (tree->current->boundary_p2[0]-tree->current->boundary_p1[0])*(tree->current->boundary_p2[0]-tree->current->boundary_p1[0]);

		  if( rcm2cell < (tree->current->rcrit_angle)*(tree->current->rcrit_angle) || rcm2cell < 5.83*boxsize2){

			  // includes rcrit_particle constraint
			  allowDescent=true;


			  // Treat all particles in a leaf as a point particle
			  if(tree->atLeaf()){

				  for(i = 0 ; i < tree->current->nparticles ; ++i){

					  xcm[0] = tree->xp[tree->current->particles[i]][0] - ray[0];
					  xcm[1] = tree->xp[tree->current->particles[i]][1] - ray[1];

					  rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
					  if(rcm2 < 1e-20) rcm2 = 1e-20;

					  index = MultiMass*tree->current->particles[i];

					  if(haloON ){
						  prefac = halos[index]->get_mass();
					  }else{
						  prefac = masses[index];
					  }
					  prefac /= rcm2*pi;

					  tmp = -1.0*prefac;

					  alpha[0] += tmp*xcm[0];
					  alpha[1] += tmp*xcm[1];

					  // can turn off kappa and gamma calculations to save times
					  if(!no_kappa){
						  tmp = -2.0*prefac/rcm2;

						  gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
						  gamma[1] += xcm[0]*xcm[1]*tmp;
					  }
				  }
			  }

			  // Fined the particles that are intersect with ray and add them individually.
			  if(rcm2cell < 5.83*boxsize2){
				  for(i = 0 ; i < tree->current->Nbig_particles ; ++i){

					  index = tree->current->big_particles[i];

					  xcm[0] = tree->xp[index][0] - ray[0];
					  xcm[1] = tree->xp[index][1] - ray[1];

					  if(haloON){
						  halos[index]->force_halo(alpha,kappa,gamma,xcm,no_kappa);
					  }else{  // case of no halos just particles and no class derived from QuadTree

						  rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
						  if(rcm2 < 1e-20) rcm2 = 1e-20;
						  rcm = sqrt(rcm2);

						  prefac = masses[MultiMass*index]/rcm2/pi;
						  arg1 = rcm2/(sizes[index*MultiRadius]*sizes[index*MultiRadius]);
						  arg2 = sizes[index*MultiRadius];
						  tmp = sizes[index*MultiRadius];

						  /// intersecting, subtract the point particle
						  if(rcm2 < tmp*tmp){
							  tmp = (alpha_h(arg1,arg2) + 1.0)*prefac;
							  alpha[0] += tmp*xcm[0];
							  alpha[1] += tmp*xcm[1];

							  // can turn off kappa and gamma calculations to save times
							  if(!no_kappa){
								  *kappa += kappa_h(arg1,arg2)*prefac;

								  tmp = (gamma_h(arg1,arg2) + 2.0)*prefac/rcm2;

								  gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
								  gamma[1] += xcm[0]*xcm[1]*tmp;
							  }
						  }
					  }
				  }
			  }

		  }else{ // use whole cell
			  allowDescent=false;

			  tmp = -1.0*tree->current->mass/rcm2cell/pi;

			  alpha[0] += tmp*xcm[0];
			  alpha[1] += tmp*xcm[1];

			  if(!no_kappa){      //  taken out to speed up
				  tmp=-2.0*tree->current->mass/pi/rcm2cell/rcm2cell;
				  gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
				  gamma[1] += xcm[0]*xcm[1]*tmp;
			  }

			  // quadrapole contribution
			  //   the kappa and gamma are not calculated to this order
			  alpha[0] -= (tree->current->quad[0]*xcm[0] + tree->current->quad[2]*xcm[1])
    				  /(rcm2cell*rcm2cell)/pi;
			  alpha[1] -= (tree->current->quad[1]*xcm[1] + tree->current->quad[2]*xcm[0])
    				  /(rcm2cell*rcm2cell)/pi;

			  tmp = 4*(tree->current->quad[0]*xcm[0]*xcm[0] + tree->current->quad[1]*xcm[1]*xcm[1]
				  + 2*tree->current->quad[2]*xcm[0]*xcm[1])/(rcm2cell*rcm2cell*rcm2cell)/pi;

			  alpha[0] += tmp*xcm[0];
			  alpha[1] += tmp*xcm[1];
		  }
	  }
  }while(tree->WalkStep(allowDescent));


  // Subtract off uniform mass sheet to compensate for the extra mass
  //  added to the universe in the halos.
  alpha[0] -= ray[0]*sigma_background;
  alpha[1] -= ray[1]*sigma_background;
  if(!no_kappa){      //  taken out to speed up
	  *kappa -= sigma_background;
  }

  return;
}

/** \brief Force2D_recur calculates the defection, convergence and shear using
 *   the plane-lens approximation.
 *
 *  This function should do the same work as QuadTree::force2D() accept it is
 *  done recursively instead of iteratively.  This is done to enable multi-threading
 *  of the force calculation.
 *
 *       The output alpha[] is in units of mass_scale/Mpc, ie it needs to be
 *       divided by Sigma_crit and multiplied by mass_scale to be the deflection
 *       in the lens equation expressed on the lens plane or multiplied by
 *       4*pi*G*mass_scale to get the deflection angle caused by the plane lens.
 *
 *       kappa and gamma need to by multiplied by mass_scale/Sigma_crit to get
 *       the traditional units for these where Sigma_crit are in the mass/units(ray)^2
 *       NB : the units of sigma_backgound need to be mass/units(ray)^2
 * */

void QuadTree::force2D_recur(double *ray,double *alpha,KappaType *kappa,KappaType *gamma,bool no_kappa){

  assert(tree);

  alpha[0]=alpha[1]=gamma[0]=gamma[1]=gamma[2]=0.0;
  *kappa=0.0;

  walkTree_recur(tree->top,&ray[0],&alpha[0],kappa,&gamma[0],no_kappa);

  // Subtract off uniform mass sheet to compensate for the extra mass
  //  added to the universe in the halos.
  alpha[0] -= ray[0]*sigma_background;
  alpha[1] -= ray[1]*sigma_background;
  if(!no_kappa){      //  taken out to speed up
	  *kappa -= sigma_background;
  }

  return;
}

void QuadTree::walkTree_recur(QBranchNB *branch,double *ray,double *alpha,KappaType *kappa,KappaType *gamma,bool no_kappa){

	PosType xcm[2],rcm2cell,rcm2,tmp,boxsize2;
	IndexType i;
	unsigned long count=0,index;
	double rcm, arg1, arg2, prefac;

	if(branch->nparticles > 0){
		xcm[0]=branch->center[0]-ray[0];
		xcm[1]=branch->center[1]-ray[1];

		rcm2cell = xcm[0]*xcm[0] + xcm[1]*xcm[1];

		boxsize2 = (branch->boundary_p2[0]-branch->boundary_p1[0])*(branch->boundary_p2[0]-branch->boundary_p1[0]);

		if( rcm2cell < (branch->rcrit_angle)*(branch->rcrit_angle) || rcm2cell < 5.83*boxsize2){

			// Treat all particles in a leaf as a point particle
			if(tree->atLeaf(branch)){

				for(i = 0 ; i < branch->nparticles ; ++i){

					xcm[0] = tree->xp[branch->particles[i]][0] - ray[0];
					xcm[1] = tree->xp[branch->particles[i]][1] - ray[1];

					rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
					if(rcm2 < 1e-20) rcm2 = 1e-20;

					index = MultiMass*branch->particles[i];


					  if(haloON ){
						  prefac = halos[index]->get_mass();
					  }else{
						  prefac = masses[index];
					  }
					  prefac /= rcm2*pi;

					alpha[0] += -1.0*prefac*xcm[0];
					alpha[1] += -1.0*prefac*xcm[1];

					// can turn off kappa and gamma calculations to save times
					if(!no_kappa){
						tmp = -2.0*prefac/rcm2;

						gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
						gamma[1] += xcm[0]*xcm[1]*tmp;
					}
				}
			}

			// Fined the particles that intersect with ray and add them individually.
			if(rcm2cell < 5.83*boxsize2){
				for(i = 0 ; i < branch->Nbig_particles ; ++i){

					index = branch->big_particles[i];

					xcm[0] = tree->xp[index][0] - ray[0];
					xcm[1] = tree->xp[index][1] - ray[1];

					/////////////////////////////////////////
					if(haloON){
						halos[index]->force_halo(alpha,kappa,gamma,xcm,no_kappa);
					}else{  // case of no halos just particles and no class derived from QuadTree

						rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
						if(rcm2 < 1e-20) rcm2 = 1e-20;
						rcm = sqrt(rcm2);

						prefac = masses[MultiMass*branch->particles[i]]/rcm2/pi;
						arg1 = rcm2/(sizes[index*MultiRadius]*sizes[index*MultiRadius]);
						arg2 = sizes[index*MultiRadius];
						tmp = sizes[index*MultiRadius];

						/// intersecting, subtract the point particle
						if(rcm2 < tmp*tmp){
							tmp = (alpha_h(arg1,arg2) + 1.0)*prefac;
							alpha[0] += tmp*xcm[0];
							alpha[1] += tmp*xcm[1];

							// can turn off kappa and gamma calculations to save times
							if(!no_kappa){
								*kappa += kappa_h(arg1,arg2)*prefac;
								tmp = (gamma_h(arg1,arg2) + 2.0)*prefac/rcm2;

								gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
								gamma[1] += xcm[0]*xcm[1]*tmp;
							}
						}
					}
				}
			}

			if(branch->child0 != NULL)
				walkTree_recur(branch->child0,&ray[0],&alpha[0],kappa,&gamma[0],no_kappa);
			if(branch->child1 != NULL)
				walkTree_recur(branch->child1,&ray[0],&alpha[0],kappa,&gamma[0],no_kappa);
			if(branch->child2 != NULL)
				walkTree_recur(branch->child2,&ray[0],&alpha[0],kappa,&gamma[0],no_kappa);
			if(branch->child3 != NULL)
				walkTree_recur(branch->child3,&ray[0],&alpha[0],kappa,&gamma[0],no_kappa);

		}else{ // use whole cell
			tmp = -1.0*branch->mass/rcm2cell/pi;

			alpha[0] += tmp*xcm[0];
			alpha[1] += tmp*xcm[1];

			if(!no_kappa){      //  taken out to speed up
				tmp=-2.0*branch->mass/pi/rcm2cell/rcm2cell;
				gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
				gamma[1] += xcm[0]*xcm[1]*tmp;
			}

			// quadrapole contribution
			//   the kappa and gamma are not calculated to this order
			alpha[0] -= (branch->quad[0]*xcm[0] + branch->quad[2]*xcm[1])
	    								  /(rcm2cell*rcm2cell)/pi;
			alpha[1] -= (branch->quad[1]*xcm[1] + branch->quad[2]*xcm[0])
	    								  /(rcm2cell*rcm2cell)/pi;

			tmp = 4*(branch->quad[0]*xcm[0]*xcm[0] + branch->quad[1]*xcm[1]*xcm[1]
					+ 2*branch->quad[2]*xcm[0]*xcm[1])/(rcm2cell*rcm2cell*rcm2cell)/pi;

			alpha[0] += tmp*xcm[0];
			alpha[1] += tmp*xcm[1];

			return;
		}
	}
	return;
}

/** This is a diagnostic routine that prints the position of every point in a
 * given branch of the tree.
 */
void QuadTree::printParticlesInBranch(unsigned long number){
	unsigned long i;

	tree->moveTop();
	do{
		if(tree->current->number == number){
			std::cout << tree->current->nparticles << std::endl;
			for(i=0;i<tree->current->nparticles;++i){
				std::cout << xp[tree->current->particles[i]][0] << "  " << xp[tree->current->particles[i]][1] << std::endl;
			}
			return;
		}
	}while(tree->WalkStep(true));

	return;
}
/**
 * Prints to stdout the borders of each branch in the tree below level.
 * If level < 0 or not specified the whole tree will be printed.
 */
void QuadTree::printBranchs(int level){

	bool decend = true;
	tree->moveTop();
	do{
		std::cout << tree->current->boundary_p1[0] << "  " << tree->current->boundary_p1[1] << "   "
		     << tree->current->boundary_p2[0] << "  " << tree->current->boundary_p2[1] << std::endl;
		if(tree->current->level == level) decend = false;
		else decend = true;
	}while(tree->WalkStep(decend));

	return;
}
