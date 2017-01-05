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
#include "Tree.h"


/** \brief Constructor meant for point particles, simulation particles
 */
TreeQuad::TreeQuad(
                   PosType **xpt
                   ,float *my_masses
                   ,float *my_sizes
                   ,IndexType Npoints
                   ,bool Multimass
                   ,bool Multisize
                   ,PosType my_sigma_background /// background kappa that is subtracted
                   ,int bucket
                   ,PosType theta_force
                   ,bool my_periodic_buffer  /// if true a periodic buffer will be imposed in the force calulation.  See documentation on TreeQuad::force2D() for details.  See note for TreeQuad::force2D_recur().
                   ,PosType my_inv_screening_scale   /// the inverse of the square of the sreening length. See note for TreeQuad::force2D_recur().
                   ):
	xp(xpt),MultiMass(Multimass), MultiRadius(Multisize), masses(my_masses)
  ,sizes(my_sizes),Nparticles(Npoints),sigma_background(my_sigma_background)
  ,Nbucket(bucket),force_theta(theta_force),periodic_buffer(my_periodic_buffer)
  ,inv_screening_scale2(my_inv_screening_scale*my_inv_screening_scale)
{
  
  
	index = new IndexType[Npoints];
	IndexType ii;

	for(ii=0;ii<Npoints;++ii) index[ii] = ii;

	haloON = false; // don't use internal halo parameters
	halos = NULL;

  tree = BuildQTreeNB_iter(xp,Npoints,index);
  //tree = BuildQTreeNB(xp,Npoints,index);  // ???
  //std::cout << "Nbranches " << tree->getNbranches() << std::endl;

	CalcMoments();

  // test lines ????
  /*{
    size_t i=0;
    tree->moveTop();
    QBranchNB *tmp;
    do{
      ++i;
      tmp = tree->current;
      assert(tree->current->constructed);
      
      std::cout << tree->current->level << "," << tree->current->nparticles << ","
      << tree->current->Nbig_particles << "," << tree->current->number << ","
      << tree->current->mass << "|"
      << tree->current->boundary_p1[0] << "," << tree->current->boundary_p1[1] << ","
      << tree->current->boundary_p2[0] << "," << tree->current->boundary_p2[1]
      << std::endl;
      
    }while( tree->WalkStep(true) );
    assert(i == tree->getNbranches());
    exit(0);
  }*/

  phiintconst = (120*log(2.) - 631.)/840 + 19./70;
	return;
}
/** \brief Constructor meant for halos with internal structure parameters.  This is a protected constructor because
 * it should only be invoked from the derived classes that have specific defined halo models.
 */
TreeQuad::TreeQuad(
                   LensHaloHndl *my_halos     /// list of halos to be put in tree
                   ,IndexType Npoints          /// number of halos
                   ,PosType my_sigma_background /// background kappa that is subtracted
                   ,int bucket                  /// maximum number of halos in each leaf of the tree
                   ,PosType theta_force         /// the openning angle used in the criterion for decending into a subcell
                   ,bool periodic_buffer        /// if true a periodic buffer will be imposed in the force calulation.  See documentation on TreeQuad::force2D() for details. See note for TreeQuad::force2D_recur().                   
                   ,PosType my_inv_screening_scale   /// the inverse of the square of the sreening length. See note for TreeQuad::force2D_recur().
                   ):
MultiMass(true),MultiRadius(true),masses(NULL),sizes(NULL)
,Nparticles(Npoints),sigma_background(my_sigma_background),Nbucket(bucket)
,force_theta(theta_force),halos(my_halos),periodic_buffer(periodic_buffer)
,inv_screening_scale2(my_inv_screening_scale*my_inv_screening_scale)
{
	index = new IndexType[Npoints];
	IndexType ii;

	for(ii=0;ii<Npoints;++ii) index[ii] = ii;
  
  // copy halo positions into xp array to make compatable with QTreeNB
  xp = Utilities::PosTypeMatrix(Npoints,2);
  for(ii=0;ii<Npoints;++ii) halos[ii]->getX(xp[ii]);

	haloON = true;  //use internal halo parameters

  tree = BuildQTreeNB_iter(xp,Npoints,index);
  //tree = BuildQTreeNB(xp,Npoints,index);
  
	CalcMoments();

  phiintconst = (120*log(2.) - 631.)/840 + 19./70;
	return;
}

/// Particle positions and other data are not destroyed.
TreeQuad::~TreeQuad()
{
	delete tree;
	delete[] index;
  if(haloON) Utilities::free_PosTypeMatrix(xp,Nparticles,2);
	return;
}

QTreeNBHndl TreeQuad::BuildQTreeNB(PosType **xp,IndexType Nparticles,IndexType *particles){
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
  
  if(Nparticles <= 1){
    p1[0]=-0.25;
    p1[1]=-0.25;
    p2[0]=0.25;
    p2[1]=0.25;
  }
  
  // store true dimensions of simulation
  PosType lengths[2] = {p2[0]-p1[0],p2[1]-p1[1]};
  original_xl = lengths[0];
  original_yl = lengths[1];
  
  // If region is not square, make it square.
  j = lengths[0] > lengths[1] ? 1 : 0;
  p2[j] = p1[j] + lengths[!j];
  
  /* Initialize tree root */
  tree = new QTreeNB(xp,particles,Nparticles,p1,p2);
  
  /* build the tree */
  _BuildQTreeNB();
  
  /* visit every branch to find center of mass and cutoff scale */
  tree->moveTop();
  
  return tree;
}

// This is a iterative replacement for the recursive combination of BuildQTreeNB() and _BuildQTreeNB()
QTreeNBHndl TreeQuad::BuildQTreeNB_iter(PosType **xp,IndexType Nparticles,IndexType *particles){
  
  //std::cout << "In TreeQuad::BuildQTreeNB_iter()" << std::endl;
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
  
  if(Nparticles <= 1){
    p1[0]=-0.25;
    p1[1]=-0.25;
    p2[0]=0.25;
    p2[1]=0.25;
  }
  
  // store true dimensions of simulation
  PosType lengths[2] = {p2[0]-p1[0],p2[1]-p1[1]};
  original_xl = lengths[0];
  original_yl = lengths[1];
  
  // If region is not square, make it square.
  j = lengths[0] > lengths[1] ? 1 : 0;
  p2[j] = p1[j] + lengths[!j];
  
  // Initialize tree root
  tree = new QTreeNB(xp,particles,Nparticles,p1,p2);
  QBranchNB *branch = tree->top;
  for(;;){
    
    if(branch->constructed){
      if(branch->child0->constructed == false){
        branch = branch->child0;
      }else if(branch->child1->constructed == false){
        branch = branch->child1;
      }else if(branch->child2->constructed == false){
        branch = branch->child2;
      }else if(branch->child3->constructed == false){
        branch = branch->child3;
      }else{
        if(branch == tree->top) break;
        branch = branch->prev;
      }
    }else{
      // make the children of branch
      if(!MakeChildren(branch)){
        if(branch == tree->top) break;
        branch = branch->prev;
      }
    }
  }
  
  //std::cout << "Out TreeQuad::BuildQTreeNB_iter(), " << tree->getNbranches() << " Branches" << std::endl;


  return tree;
}

/// returns an index for which of the four quadrangles of the branch the point x[] is in
inline short TreeQuad::WhichQuad(PosType *x,QBranchNB &branch){
	return (x[0] < branch.center[0]) + 2*(x[1] < branch.center[1]);
}

/// tree must be created and first branch must be set before start
//void TreeQuad::_BuildQTreeNB(IndexType nparticles,IndexType *particles){
void TreeQuad::_BuildQTreeNB(){
  
  if(MakeChildren(tree->current)){
  
    tree->moveToChild(0);
    //_BuildQTreeNB(child0->nparticles,child0->particles);
    _BuildQTreeNB();
    tree->moveUp();
  
    tree->moveToChild(1);
    //_BuildQTreeNB(child1->nparticles,child1->particles);
    _BuildQTreeNB();
    tree->moveUp();
  
    tree->moveToChild(2);
    //_BuildQTreeNB(child2->nparticles,child2->particles);
    _BuildQTreeNB();
    tree->moveUp();
  
    tree->moveToChild(3);
    //_BuildQTreeNB(child3->nparticles,child3->particles);
    _BuildQTreeNB();
    tree->moveUp();
    
  }
}

// return true if branch cbranch is divided and false otherwise
bool TreeQuad::MakeChildren(QBranchNB *cbranch){

  if(cbranch->constructed) return false;
  cbranch->constructed = true;
	//QBranchNB *cbranch = tree->current; // pointer to current branch
  IndexType *particles = cbranch->particles;

  IndexType i,j,cut,cut2,jt;

	cbranch->center[0] = (cbranch->boundary_p1[0] + cbranch->boundary_p2[0])/2;
	cbranch->center[1] = (cbranch->boundary_p1[1] + cbranch->boundary_p2[1])/2;
	cbranch->quad[0] = cbranch->quad[1] = cbranch->quad[2] = 0;
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
				if(r < (cbranch->boundary_p2[0]-cbranch->boundary_p1[0])) cbranch->big_particles[j++] = particles[i];
			}
		}
		else{
			cbranch->big_particles = NULL;
		}

		return false;
	}

	// find particles too big to be in children

	//PosType *x = new PosType[cbranch->nparticles];

	cbranch->Nbig_particles = 0;

	if(MultiRadius){
		// store the particles that are too large to be in a child at the end of the list of particles in cbranch
		/*for(i=0;i<cbranch->nparticles;++i){
			if(haloON){
				x[i] = halos[particles[i]]->get_Rmax();
			}else{
				x[i] = sizes[particles[i]];
			}
		}*/
		PosType maxsize = (cbranch->boundary_p2[0]-cbranch->boundary_p1[0])/2;

		// sort particles in size
		//Utilities::quickPartition(maxsize,&cut,particles,x,cbranch->nparticles);
    
    if(haloON){
      LensHalo **tmp = halos;
      cut = Utilities::Partition(maxsize,particles,cbranch->nparticles
                                  ,[tmp](IndexType i){return tmp[i]->get_Rmax();}  );
    }else{
      float *tmp = sizes;
      cut = Utilities::Partition(maxsize,particles,cbranch->nparticles
                                 ,[tmp](IndexType i){return tmp[i];}  );
    }


		if(cut < cbranch->nparticles){
			if(cbranch == tree->top){
				cbranch->Nbig_particles = cut2 = cbranch->nparticles-cut;
			}else{
				//Utilities::quickPartition(2*maxsize,&cut2,&particles[cut],&x[cut],cbranch->nparticles-cut);
        
        if(haloON){
          LensHalo **tmp = halos;
          cut2 = Utilities::Partition(2*maxsize,&particles[cut],cbranch->nparticles - cut
                                     ,[tmp](IndexType i){return tmp[i]->get_Rmax();}  );
        }else{
          float *tmp = sizes;
          cut2 = Utilities::Partition(2*maxsize,&particles[cut],cbranch->nparticles - cut
                                     ,[tmp](IndexType i){return tmp[i];}  );
        }
        
				cbranch->Nbig_particles = cut2;
			}

			cbranch->big_particles = new IndexType[cbranch->Nbig_particles];
			for(i=cut;i<(cut+cut2);++i) cbranch->big_particles[i-cut] = particles[i];
		}

	}else{
    float tmp;
		if(haloON){
			tmp = halos[0]->get_Rmax();
		}else{
			tmp = sizes[0];
		}
		if(tmp > (cbranch->boundary_p2[0]-cbranch->boundary_p1[0])/2
				&& tmp < (cbranch->boundary_p2[0]-cbranch->boundary_p1[0])){
			cbranch->Nbig_particles = cbranch->nparticles;
			cbranch->big_particles = cbranch->particles;
		}
	}

	assert( cbranch->nparticles >= cbranch->Nbig_particles );
	IndexType NpInChildren = cbranch->nparticles;// - cbranch->Nbig_particles;
	assert(NpInChildren >= 0);

	if(NpInChildren == 0){
		//delete[] x;
		return false;
	}

	IndexType cutx,cuty;
	PosType xcut,ycut;

	QBranchNB *child0;
	QBranchNB *child1;
	QBranchNB *child2;
	QBranchNB *child3;

	//tree->attachChildrenToCurrent(child0,child1,child2,child3);
  {  // construct children
    tree->Nbranches += 4;
    cbranch->child0 = child0 = new QBranchNB();
    assert(child0->constructed == false);
    cbranch->child1 = child1 = new QBranchNB();
        assert(child1->constructed == false);
    cbranch->child2 = child2 = new QBranchNB();
        assert(child2->constructed == false);
    cbranch->child3 = child3 = new QBranchNB();
        assert(child3->constructed == false);
    
    int level = cbranch->level + 1;
    child0->level = level;
    child1->level = level;
    child2->level = level;
    child3->level = level;
  
    // work out brothers for children
    child0->brother = child1;
    child1->brother = child2;
    child2->brother = child3;
    child3->brother = cbranch->brother;
  
    child0->prev = cbranch;
    child1->prev = cbranch;
    child2->prev = cbranch;
    child3->prev = cbranch;
  }

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
	//for(i=0;i<NpInChildren;++i) x[i] = tree->xp[particles[i]][1];
	//Utilities::quickPartition(ycut,&cuty,particles,x,NpInChildren);
  
  {
    PosType **tmp = tree->xp;
    cuty = Utilities::Partition(ycut,particles,NpInChildren
                       ,[tmp](IndexType i){return tmp[i][1];}  );
  }

	if(cuty > 0){
	  // divide first group in the x direction
      //for(i=0;i<cuty;++i) x[i] = tree->xp[particles[i]][0];
      //Utilities::quickPartition(xcut,&cutx,particles,x,cuty);

    {
      PosType **tmp = tree->xp;
      cutx = Utilities::Partition(xcut,particles,cuty
                                  ,[tmp](IndexType i){return tmp[i][0];}  );
    }

    
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
	  //for(i=cuty;i<NpInChildren;++i) x[i-cuty] = tree->xp[particles[i]][0];
	  //Utilities::quickPartition(xcut,&cutx,&particles[cuty],x,NpInChildren-cuty);
    {
      PosType **tmp = tree->xp;
      cutx = Utilities::Partition(xcut,&particles[cuty],NpInChildren-cuty
                                  ,[tmp](IndexType i){return tmp[i][0];}  );
    }

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

	//delete[] x;

	return true;
}

// calculates moments of the mass and the cutoff scale for each box
void TreeQuad::CalcMoments(){

	//*** make compatable
	IndexType i;
	PosType rcom,xcm[2],xcut;
	QBranchNB *cbranch;
	PosType tmp;
  size_t Ncheck = 0;

  
	tree->moveTop();
	do{
		cbranch=tree->current; /* pointer to current branch */

		cbranch->rmax = sqrt( pow(cbranch->boundary_p2[0]-cbranch->boundary_p1[0],2)
				+ pow(cbranch->boundary_p2[1]-cbranch->boundary_p1[1],2) );

		// calculate mass
		for(i=0,cbranch->mass=0;i<cbranch->nparticles;++i)
			if(haloON ){
				cbranch->mass += halos[cbranch->particles[i]]->get_mass();
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

    ++Ncheck;
	}while(tree->WalkStep(true));

  assert(Ncheck == tree->getNbranches());
	return;
}

/// simple rotates the coordinates in the xp array
void TreeQuad::rotate_coordinates(PosType **coord){
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
 *   the plane-lens approximation.  The tree is walked iteratively.
 *
 *       The output alpha[] is in units of mass_scale/Mpc, ie it needs to be
 *       divided by Sigma_crit and multiplied by mass_scale to be the deflection
 *       in the lens equation expressed on the lens plane or multiplied by
 *       4*pi*G*mass_scale to get the deflection angle caused by the plane lens.
 *
 *       kappa and gamma need to by multiplied by mass_scale/Sigma_crit to get
 *       the traditional units for these where Sigma_crit are in the mass/units(ray)^2
 *       NB : the units of sigma_backgound need to be mass/units(ray)^2
 *
 *   If periodic_buffer == true a periodic buffer is included.  The ray is calculated as if there
 *   where identical copies of it on the borders of the tree's root branch.  This is not the same thing as periodic boundary conditions, because
 *   there are not an infinite number of copies in all direction as for the DFT force solver.
 *
 *   If inv_screening_scale2 != 0 the mass of cells are reduced by a factor of exp(-|ray - center of mass|^2*inv_screening_scale2) which
 *   screens the large scale geometry of the simulation on the sky.  This is useful when the region is rectangular instead of circular.
 * */

void TreeQuad::force2D(const PosType *ray,PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi) const{
  
  assert(tree);
  QTreeNB::iterator it(tree);
  
  alpha[0]=alpha[1]=gamma[0]=gamma[1]=gamma[2]=0.0;
  *kappa=*phi=0.0;
    
  if(periodic_buffer){
    PosType tmp_ray[2];
        
    for(int i = -1;i<2;++i){
      for(int j = -1;j<2;++j){
        tmp_ray[0] = ray[0] + i*original_xl;
        tmp_ray[1] = ray[1] + j*original_yl;
        walkTree_iter(it,tmp_ray,alpha,kappa,gamma,phi);
      }
    }
  }else{
    walkTree_iter(it,ray,alpha,kappa,gamma,phi);
  }
  
  // Subtract off uniform mass sheet to compensate for the extra mass
  //  added to the universe in the halos.
  alpha[0] -= ray[0]*sigma_background*(inv_screening_scale2 == 0);
  alpha[1] -= ray[1]*sigma_background*(inv_screening_scale2 == 0);
  
  *kappa -= sigma_background;
  
  return;
}

/** \brief Returns the halos that are within rmax of ray[]
 */
void TreeQuad::neighbors(PosType ray[],PosType rmax,std::list<IndexType> &neighbors) const{
  QTreeNB::iterator it(tree);
  neighbors.clear();
  
  PosType r2 = rmax*rmax;
  bool decend=true;
  it.movetop();
  do{
    int cut = Utilities::cutbox(ray,(*it)->boundary_p1,(*it)->boundary_p2,rmax);
    decend = true;
    
    if(cut == 0){      // box outside circle
      decend = false;
    }else if(cut == 1){  // whole box inside circle
      for(int i=0;i<(*it)->nparticles;++i) neighbors.push_back((*it)->particles[i]);
      decend = false;
    }else if((*it)->child0 == NULL){  // at leaf
      for(int i=0;i<(*it)->nparticles;++i){
        if( (tree->xp[(*it)->particles[i]][0] - ray[0])*(tree->xp[(*it)->particles[i]][0] - ray[0])
           + (tree->xp[(*it)->particles[i]][1] - ray[1])*(tree->xp[(*it)->particles[i]][1] - ray[1]) < r2) neighbors.push_back((*it)->particles[i]);
      }
      decend = false;
    }
  }while(it.TreeWalkStep(decend));
  
}
/** \brief Returns the halos that are within rmax of ray[]
 */
void TreeQuad::neighbors(PosType ray[],PosType rmax,std::vector<LensHalo *> &neighbors) const{
  neighbors.clear();
  
  if(halos == NULL){
    std::cerr << "TreeQuad::neighbors - The are no halos in this tree use other version of this function" << std::endl;
    throw std::runtime_error("no halos");
    return;
  }
  
  std::list<IndexType> neighbor_indexes;
  TreeQuad::neighbors(ray,rmax,neighbor_indexes);
  
  neighbors.resize(neighbor_indexes.size());
  std::list<IndexType>::iterator it = neighbor_indexes.begin();
  for(size_t i=0;i<neighbors.size();++i,++it){
      neighbors[i] = halos[*it];
  }
  
  return;
}

void TreeQuad::walkTree_iter(
                             QTreeNB::iterator &treeit,
                             const PosType *ray
                             ,PosType *alpha
                             ,KappaType *kappa
                             ,KappaType *gamma
                             ,KappaType *phi
                             ) const {
  
  PosType xcm[2],rcm2cell,rcm2,tmp,boxsize2;
  IndexType i;
  bool allowDescent=true;
  unsigned long count=0,tmp_index;
  PosType prefac;
  //bool notscreen = true;
  
  assert(tree);
  
  treeit.movetop();
  //tree->moveTop();
  do{
	  ++count;
	  allowDescent=false;
    
    /*
    if(inv_screening_scale2 != 0) notscreen = BoxIntersectCircle(ray,3*sqrt(1.0/inv_screening_scale2)
                                                                 ,(*treeit)->boundary_p1
                                                                 ,(*treeit)->boundary_p2);
    else notscreen = true;
*/
	  //if(tree->current->nparticles > 0)
    if((*treeit)->nparticles > 0)
    {
      
		  xcm[0]=(*treeit)->center[0]-ray[0];
		  xcm[1]=(*treeit)->center[1]-ray[1];
      
		  rcm2cell = xcm[0]*xcm[0] + xcm[1]*xcm[1];
      
		  boxsize2 = ((*treeit)->boundary_p2[0]-(*treeit)->boundary_p1[0])*((*treeit)->boundary_p2[0]-(*treeit)->boundary_p1[0]);
      
		  if( rcm2cell < ((*treeit)->rcrit_angle)*((*treeit)->rcrit_angle) || rcm2cell < 5.83*boxsize2)
      {
        
			  // includes rcrit_particle constraint
			  allowDescent=true;
        
        
			  // Treat all particles in a leaf as a point particle
			  if(treeit.atLeaf())
        {
          
				  for(i = 0 ; i < (*treeit)->nparticles ; ++i)
          {
            
					  xcm[0] = tree->xp[(*treeit)->particles[i]][0] - ray[0];
					  xcm[1] = tree->xp[(*treeit)->particles[i]][1] - ray[1];
            

					  rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
            double screening = exp(-rcm2*inv_screening_scale2);
					  if(rcm2 < 1e-20) rcm2 = 1e-20;
            
					  tmp_index = MultiMass*(*treeit)->particles[i];
            
					  if(haloON){ prefac = halos[tmp_index]->get_mass(); }
            else{ prefac = masses[tmp_index]; }
					  prefac /= rcm2*pi/screening;
            
					  tmp = -1.0*prefac;
            
					  alpha[0] += tmp*xcm[0];
					  alpha[1] += tmp*xcm[1];
            
					  
            tmp = -2.0*prefac/rcm2;
              
            gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
            gamma[1] += xcm[0]*xcm[1]*tmp;
              
            *phi += prefac*rcm2*0.5*log(rcm2);
          
				  } // end of for
			  } // end of if(tree->atLeaf())
        
        
			  // Find the particles that intersect with ray and add them individually.
			  if(rcm2cell < 5.83*boxsize2)
        {
          
				  for(i = 0 ; i < (*treeit)->Nbig_particles ; ++i)
          {
            
					  tmp_index = (*treeit)->big_particles[i];
            
					  xcm[0] = tree->xp[tmp_index][0] - ray[0];
					  xcm[1] = tree->xp[tmp_index][1] - ray[1];
            
            rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
					  if(haloON){
              double screening = exp(-rcm2*inv_screening_scale2);

						  halos[tmp_index]->force_halo(alpha,kappa,gamma,phi,xcm,true,screening);
					  }else{  // case of no halos just particles and no class derived from TreeQuad
              
						  if(rcm2 < 1e-20) rcm2 = 1e-20;
						  //rcm = sqrt(rcm2);
              
						  //prefac = masses[MultiMass*tmp_index]/rcm2/pi;
						  //arg1 = rcm2/(sizes[tmp_index*MultiRadius]*sizes[tmp_index*MultiRadius]);
						  //arg2 = sizes[tmp_index*MultiRadius];
						  //tmp = sizes[tmp_index*MultiRadius];
              
              PosType size = sizes[tmp_index*MultiRadius];
              
						  // intersecting, subtract the point particle
						  if(rcm2 < 4*size*size)
              {
                
                b_spline_profile(xcm,sqrt(rcm2),masses[MultiMass*tmp_index],size,alpha,kappa,gamma,phi);
						  }
					  }
				  }
			  }
        
		  }
      else
      { // use whole cell
        
			  allowDescent=false;
                
        double screening = exp(-rcm2cell*inv_screening_scale2);
        double tmp = -1.0*(*treeit)->mass/rcm2cell/pi*screening;
        
        alpha[0] += tmp*xcm[0];
        alpha[1] += tmp*xcm[1];
        
        {      //  taken out to speed up
          tmp=-2.0*(*treeit)->mass/pi/rcm2cell/rcm2cell*screening;
          gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
          gamma[1] += xcm[0]*xcm[1]*tmp;
          
          *phi += 0.5*(*treeit)->mass*log( rcm2cell )/pi*screening;
          *phi -= 0.5*( (*treeit)->quad[0]*xcm[0]*xcm[0] + (*treeit)->quad[1]*xcm[1]*xcm[1] + 2*(*treeit)->quad[2]*xcm[0]*xcm[1] )/(pi*rcm2cell*rcm2cell)*screening;
        }
        
        // quadrapole contribution
        //   the kappa and gamma are not calculated to this order
        alpha[0] -= ((*treeit)->quad[0]*xcm[0] + (*treeit)->quad[2]*xcm[1])
        /(rcm2cell*rcm2cell)/pi*screening;
        alpha[1] -= ((*treeit)->quad[1]*xcm[1] + (*treeit)->quad[2]*xcm[0])
        /(rcm2cell*rcm2cell)/pi*screening;
        
        tmp = 4*((*treeit)->quad[0]*xcm[0]*xcm[0] + (*treeit)->quad[1]*xcm[1]*xcm[1]
                 + 2*(*treeit)->quad[2]*xcm[0]*xcm[1])/(rcm2cell*rcm2cell*rcm2cell)/pi*screening;
        
        alpha[0] += tmp*xcm[0];
        alpha[1] += tmp*xcm[1];

		  }
	  }
  }while(treeit.TreeWalkStep(allowDescent));
  
  
  // Subtract off uniform mass sheet to compensate for the extra mass
  //  added to the universe in the halos.
  alpha[0] -= ray[0]*sigma_background*(inv_screening_scale2 == 0.0);
  alpha[1] -= ray[1]*sigma_background*(inv_screening_scale2 == 0.0);
  {      //  taken out to speed up
	  *kappa -= sigma_background;
    // *phi -= sigma_background * sigma_background ;
  }
  
  return;
}

/** \brief Force2D_recur calculates the defection, convergence and shear using
 *   the plane-lens approximation.
 *
 *  This function should do the same work as TreeQuad::force2D() except it is
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
 *
 *   If periodic_buffer == true a periodic buffer is included.  The ray is calculated as if there
 *   where identical copies of it on the borders of the tree's root branch.  This is not the same thing as periodic boundary conditions, because
 *   there are not an infinite number of copies in all direction as for the DFT force solver.
 *
 *   If inv_screening_scale2 != 0 the mass of cells are reduced by a factor of exp(-|ray - center of mass|^2*inv_screening_scale2) which
 *   screens the large scale geometry of the simulation on the sky.  This is useful when the region is rectangular instead of circular.
 * */

void TreeQuad::force2D_recur(const PosType *ray,PosType *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi){
  
  assert(tree);
  
  alpha[0]=alpha[1]=gamma[0]=gamma[1]=gamma[2]=0.0;
  *kappa=*phi=0.0;
  
  //walkTree_recur(tree->top,ray,&alpha[0],kappa,&gamma[0],phi);
  
  if(periodic_buffer){
    PosType tmp_ray[2],tmp_c[2];
    
    // for points outside of primary zone shift to primary zone
    //if(inbox(ray,tree->top->boundary_p1,tree->top->boundary_p2)){
      tmp_c[0] = ray[0];
      tmp_c[1] = ray[1];
    /*}else{
      int di[2];
      PosType center[2] = { (tree->top->boundary_p1[0] + tree->top->boundary_p2[0])/2 ,
        (tree->top->boundary_p1[1] + tree->top->boundary_p2[1])/2 };
      
      di[0] = (int)( 2*(ray[0] - center[0])/lx );
      di[1] = (int)( 2*(ray[1] - center[1])/ly );
      
      di[0] = (di[0] > 0) ? (di[0] + 1)/2 : (di[0] - 1)/2;
      di[1] = (di[1] > 0) ? (di[1] + 1)/2 : (di[1] - 1)/2;
      
      tmp_c[0] = ray[0] - lx * di[0];
      tmp_c[1] = ray[1] - ly * di[1];
      
      assert(inbox(tmp_c,tree->top->boundary_p1,tree->top->boundary_p2));
    }*/
    
    for(int i = -1;i<2;++i){
      for(int j = -1;j<2;++j){
        tmp_ray[0] = tmp_c[0] + i*original_xl;
        tmp_ray[1] = tmp_c[1] + j*original_yl;
        walkTree_recur(tree->top,tmp_ray,alpha,kappa,gamma,phi);
      }
    }
  }else{
    walkTree_recur(tree->top,ray,alpha,kappa,gamma,phi);
  }
  
  // Subtract off uniform mass sheet to compensate for the extra mass
  //  added to the universe in the halos.
  alpha[0] -= ray[0]*sigma_background*(inv_screening_scale2 == 0);
  alpha[1] -= ray[1]*sigma_background*(inv_screening_scale2 == 0);
  
  *kappa -= sigma_background;
    
  return;
}



void TreeQuad::walkTree_recur(QBranchNB *branch,PosType const *ray,PosType *alpha,KappaType *kappa,KappaType *gamma, KappaType *phi){
  
	PosType xcm[2],rcm2cell,rcm2,tmp,boxsize2;
	IndexType i;
	std::size_t tmp_index;
	PosType prefac;
  /*bool notscreen = true;
  
  if(inv_screening_scale2 != 0) notscreen = BoxIntersectCircle(ray,3*sqrt(1.0/inv_screening_scale2), branch->boundary_p1, branch->boundary_p2);
    */
	if(branch->nparticles > 0)
  {
		xcm[0]=branch->center[0]-ray[0];
		xcm[1]=branch->center[1]-ray[1];
    
		rcm2cell = xcm[0]*xcm[0] + xcm[1]*xcm[1];
    
    boxsize2 = (branch->boundary_p2[0]-branch->boundary_p1[0])*(branch->boundary_p2[0]-branch->boundary_p1[0]);
    
    if( rcm2cell < (branch->rcrit_angle)*(branch->rcrit_angle) || rcm2cell < 5.83*boxsize2)
    {
      
      // Treat all particles in a leaf as a point particle
      if(tree->atLeaf(branch))
      {
        
        for(i = 0 ; i < branch->nparticles ; ++i){
          
          xcm[0] = tree->xp[branch->particles[i]][0] - ray[0];
          xcm[1] = tree->xp[branch->particles[i]][1] - ray[1];
          
          rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
          double screening = exp(-rcm2*inv_screening_scale2);
          if(rcm2 < 1e-20) rcm2 = 1e-20;
          
          tmp_index = MultiMass*branch->particles[i];
          
          
          if(haloON ) { prefac = halos[tmp_index]->get_mass(); }
          else{ prefac = masses[tmp_index]; }
          prefac /= rcm2*pi/screening;
          
          
          alpha[0] += -1.0*prefac*xcm[0];
          alpha[1] += -1.0*prefac*xcm[1];
          
          {
            tmp = -2.0*prefac/rcm2;
            
            
            gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
            gamma[1] += xcm[0]*xcm[1]*tmp;
            
            *phi += prefac*rcm2*0.5*log(rcm2);
          }
        }
      }
      
			// Find the particles that intersect with ray and add them individually.
			if(rcm2cell < 5.83*boxsize2)
      {
        
				for(i = 0 ; i < branch->Nbig_particles ; ++i)
        {
          
					tmp_index = branch->big_particles[i];
          
					xcm[0] = tree->xp[tmp_index][0] - ray[0];
					xcm[1] = tree->xp[tmp_index][1] - ray[1];
          rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];

					/////////////////////////////////////////
					if(haloON){
            double screening = exp(-rcm2*inv_screening_scale2);
						halos[tmp_index]->force_halo(alpha,kappa,gamma,phi,xcm,true,screening);
 					}else{  // case of no halos just particles and no class derived from TreeQuad
            
						if(rcm2 < 1e-20) rcm2 = 1e-20;
            
            PosType size = sizes[tmp_index*MultiRadius];
            
						// intersecting, subtract the point particle
						if(rcm2 < 4*size*size)
            {
              
              b_spline_profile(xcm,sqrt(rcm2),masses[MultiMass*tmp_index],size,alpha,kappa,gamma,phi);
              
						}
					}
				}
			}
      
			if(branch->child0 != NULL)
				walkTree_recur(branch->child0,&ray[0],&alpha[0],kappa,&gamma[0],&phi[0]);
			if(branch->child1 != NULL)
				walkTree_recur(branch->child1,&ray[0],&alpha[0],kappa,&gamma[0],&phi[0]);
			if(branch->child2 != NULL)
				walkTree_recur(branch->child2,&ray[0],&alpha[0],kappa,&gamma[0],&phi[0]);
			if(branch->child3 != NULL)
				walkTree_recur(branch->child3,&ray[0],&alpha[0],kappa,&gamma[0],&phi[0]);
      
		}
    else
    { // use whole cell
      
      double screening = exp(-rcm2cell*inv_screening_scale2);
			double tmp = -1.0*branch->mass/rcm2cell/pi*screening;
      
			alpha[0] += tmp*xcm[0];
			alpha[1] += tmp*xcm[1];
      
			{      //  taken out to speed up
				tmp=-2.0*branch->mass/pi/rcm2cell/rcm2cell*screening;
				gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
				gamma[1] += xcm[0]*xcm[1]*tmp;
        
        *phi += 0.5*tree->current->mass*log( rcm2cell )/pi*screening;
        *phi -= 0.5*( tree->current->quad[0]*xcm[0]*xcm[0] + tree->current->quad[1]*xcm[1]*xcm[1] + 2*tree->current->quad[2]*xcm[0]*xcm[1] )/(pi*rcm2cell*rcm2cell)*screening;
			}
      
			// quadrapole contribution
			//   the kappa and gamma are not calculated to this order
			alpha[0] -= (branch->quad[0]*xcm[0] + branch->quad[2]*xcm[1])
      /(rcm2cell*rcm2cell)/pi*screening;
			alpha[1] -= (branch->quad[1]*xcm[1] + branch->quad[2]*xcm[0])
      /(rcm2cell*rcm2cell)/pi*screening;
      
			tmp = 4*(branch->quad[0]*xcm[0]*xcm[0] + branch->quad[1]*xcm[1]*xcm[1]
               + 2*branch->quad[2]*xcm[0]*xcm[1])/(rcm2cell*rcm2cell*rcm2cell)/pi*screening;
      
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
void TreeQuad::printParticlesInBranch(unsigned long number){
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
void TreeQuad::printBranchs(int level){

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
