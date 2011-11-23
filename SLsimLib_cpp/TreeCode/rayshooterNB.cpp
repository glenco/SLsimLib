/*#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <nrutil.h>
#include <nr.h>
#include <nrD.h>
#include <cosmo.h>
#include "TreeNB.h"*/

#include <slsimlib.h>

COSMOLOGY cosmo;
int kmax,kount; 
PosType *xp,**yp,dxsav;

/** \ingroup DeflectionL2
 *  \brief Wrapper function for initializing and calculating the deflection from particles.
 *  */

void rayshooterNB(unsigned long Nrays,Point *points,TreeHndl i_tree,char *paramfile,bool no_kappa){
 
  int i,j,k;
  float *mass;
  //double r[9],rmax,dkappa=0;
  static SimLens *lens;
  static unsigned long counter=0;
  static double convert_factor=0,tmp=0;
  static ListHndl neighbors;
  static short initialize=1;
  TreeNBStruct tree;
  time_t to,t1;
  //Point *pointo;

  ++counter;

  /* read in simulation data on first use */
  if(initialize){
	initialize=0;
    if(paramfile == NULL){
      ERROR_MESSAGE();
      printf("ERROR: in rayshooter, no paramfile\n");
      exit(0);
    }
    /* read parameter file and simulation */

    printf("READING PARAMETERS FROM %s\n\n",paramfile);
    lens=readparams(paramfile,&cosmo);

    /*printf("READING PARAMETERS FROM %s\n\n",teststring);
      readparams(teststring,lens,&cosmo);*/

    printf("\nreading particle positions from %s\n\n",lens->simfilename);
    for(i=0;i<lens->Nspecies;++i){

    	readpositions(&lens[i]);

    	printf("*********    loaded    *********************\n\n");
    	printf("number of particles=%li\n",lens[i].Nparticles);

    	/*for(i=0;i<Nparticles;++i) printf("%e %e %e\n",xp[i][0],xp[i][1],xp[i][2]);*/
    	/*     printf("\n***********************************************\n"); */
    	/*     printf("********* tree from %s *********************\n\n",lens[i].treefilenames); */
    	/*lens[i].tree=readTreeNB(lens[i].particles,lens[i].rsph,lens[i].Nparticles,lens[i].treefilenames);*/

    	printf("\n***********************************************\n");
    	printf("********* reading rsph from %s *********************\n\n",lens[i].treefilenames);

    	readSmoothingNB(lens[i].rsph,lens[i].treefilenames);

    	// rescaling for Steffan's
    	float factor=4.*sqrt(103./1120.)/(1+lens[i].zlens)/cosmo.h;
    	for(j=0;j<lens[i].Nparticles;++j) lens[i].rsph[j]*=factor;
    	factor=1./(1+lens[i].zlens)/cosmo.h;

    	PosType tmp[3];
    	for(j=0;j<lens[i].Nparticles;++j){
    		for(k=0;k<3;++k) lens[i].xp[j][k] *= factor;
    		tmp[0]=lens[i].xp[j][0];
    		tmp[1]=lens[i].xp[j][1];
    		tmp[2]=lens[i].xp[j][2];

    		lens[i].xp[j][0]=tmp[1];
       		lens[i].xp[j][1]=tmp[2];
       		lens[i].xp[j][2]=tmp[0];
    	}


    	printf("\n***********************************************\n");
    	printf("*********    projecting lens    *********************\n\n");

    	mass=(float*)calloc(1,sizeof(float));
    	*mass=1.0;
    	lens[i].tree=rotate_project(lens[i].xp,lens[i].Nparticles,lens[i].particles
    		      ,lens[i].coord,lens[i].theta_force,lens[i].rsph,mass,true,false);
    	lens[i].xp_2d=lens[i].xp;

    	//lens[i].tree=rotate_simulation(lens[i].xp,lens[i].Nparticles,lens[i].particles
    	//		,lens[i].coord,lens[i].theta_force,lens[i].rsph);
    	// 2D smoothing
    	//lens[i].rsph=FindRSPH(lens[i].tree,lens[i].xp,lens[i].Nsph);

    	printf("*********    projected    ***************************\n\n");
    	//exit(0);
    }

    convert_factor=lens->mass_units/lens->Sigma_crit;
    neighbors=NewList();
 }

  // special test case of massless particles
  /*
  if(lens[0].mass_units == 0.0){
	  for(i=0;i<Nrays;++i){
		  points[i].image->x[0]=points[i].x[0];
		  points[i].image->x[1]=points[i].x[1];

		  points[i].kappa*=0;
		  points[i].gamma[0]*=0;
		  points[i].gamma[1]*=0;
	  }
	  return ;
  }
  */

//  PrintSimLens(lens);
  /*** loop through rays ***/

  time(&to);
  for(k=0;k<lens->Nspecies;++k){
	  tree=*(lens[k].tree);
	  //#pragma omp parallel for private(neighbors,tmp,rmax,pointo) firstprivate(tree) firstprivate(i_tree)
	  for(i=0;i<Nrays;++i){

		  /*     printf(" i=%i ray = %e %e  ray image = %e %e\n",i */
		  /* 	   ,points[i].x[0],points[i].x[1] */
		  /* 	   ,points[i].image->x[0],points[i].image->x[1]); */

		  /*     if( i == 0 ) {
		   nthreads = omp_get_num_threads();
			   printf("There are %d threads Nrays=%i\n",nthreads,Nrays);
       }
		   */
		  // determine if interpolation should be used or direct tree summation

		  // grid size is below interpolation scale find neighbors and
		  // see if point is on a grid boundary
/*		  if(points[i].gridsize < lens[k].interpolation_scale){
			  pointo=NearestNeighbor(i_tree,points[i].x,9,neighbors,0);
			  // convert to nearest neighbor centered coordinates
			  MoveToTopList(neighbors);
			  for(j=0,rmax=0,dkappa=0;j<neighbors->Npoints;++j){
				  r[j]=sqrt(pow(neighbors->current->x[0]-pointo->x[0],2)
						  + pow(neighbors->current->x[1]-pointo->x[1],2))/pointo->gridsize;
				  if(rmax < r[j]) rmax=r[j];
				  if(dkappa < fabs( (neighbors->current->kappa-pointo->kappa)/pointo->kappa) )
					  dkappa = fabs( (neighbors->current->kappa-pointo->kappa)/pointo->kappa);
				  MoveDownList(neighbors);
			  }

		  }else rmax=0.0;
*/
		  tmp=points[i].gridsize;

		  //if( (points[i].gridsize > lens[k].interpolation_scale)
				  //|| (fabs(rmax - 4.2426) > 0.01) || (dkappa>1.0e-2)
		//		  ){
			  // do full tree calculation

			  // 2d tree option
			  //      TreeNBForce2D(&tree,lens[k].xp_2d,lens[k].rsph,points[i].x,lens[k].theta
			  //		    ,points[i].image->x,&(points[i].kappa),points[i].gamma);

			  //3d tree option
			  TreeNBForce2D(&tree,points[i].x,points[i].image->x
					  ,&(points[i].kappa),points[i].gamma,no_kappa);

			  points[i].image->x[0]=points[i].x[0] - convert_factor*points[i].image->x[0];
			  points[i].image->x[1]=points[i].x[1] - convert_factor*points[i].image->x[1];
    
			  points[i].kappa*=convert_factor;
			  points[i].gamma[0]*=convert_factor;
			  points[i].gamma[1]*=convert_factor;

		 /* }else{ // interpolate from closest points

			  dx[0]=points[i].x[0] - pointo->x[0];
			  dx[1]=points[i].x[1] - pointo->x[1];

			  points[i].kappa = pointo->kappa;
			  points[i].gamma[0] = pointo->gamma[0];
			  points[i].gamma[1] = pointo->gamma[1];
			  MoveToTopList(neighbors);
			  for(j=0;j<9;++j){
				  if(fabs(r[j]-3) < 0.01){
					  tmp=( (neighbors->current->x[0]-pointo->x[0])*dx[0]
                        + (neighbors->current->x[1]-pointo->x[1])*dx[1] )
                        /pow(3*pointo->gridsize,2)/2;
					  //printf("%e %e\n",(neighbors->current->x[0]-pointo->x[0])/3/pointo->gridsize
					  //   ,(neighbors->current->x[1]-pointo->x[1])/3/pointo->gridsize);
					  points[i].kappa += neighbors->current->kappa*tmp;
					  points[i].gamma[0] += neighbors->current->gamma[0]*tmp;
					  points[i].gamma[1] += neighbors->current->gamma[1]*tmp;
				  }
				  MoveDownList(neighbors);
			  }
			  points[i].image->x[0]= pointo->image->x[0]
			             + (1-points[i].kappa-points[i].gamma[0])*dx[0]
			             - points[i].gamma[1]*dx[1];
			  points[i].image->x[1]= pointo->image->x[1]
			             + (1-points[i].kappa+points[i].gamma[0])*dx[1]
        	             - points[i].gamma[1]*dx[0];

			  if(fabs(points[i].kappa-pointo->kappa) > 1){
			      ERROR_MESSAGE();
				  printf("error\n");
				  printf("D kappa=%e   dx= %e  %e  gridsize=%e\n",points[i].kappa-pointo->kappa
						  ,dx[0],dx[1],pointo->gridsize);
				  PrintPoint(pointo);
				  PrintPoint(&points[i]);

			  }
		  }*/

		  points[i].invmag=pow(1-points[i].kappa,2)
    		  - pow(points[i].gamma[0],2) - pow(points[i].gamma[1],2);
     
		  points[i].image->invmag=points[i].invmag;
		  points[i].image->kappa=points[i].kappa;
		  points[i].image->gamma[0]=points[i].gamma[0];
		  points[i].image->gamma[1]=points[i].gamma[1];
 
		  /*printf("   alpha = %e %e\n   gamma=%e %e\n   kappa = %e\n   invmag = %e\n"
	   ,alpha[2*i],alpha[2*i+1],gamma[2*i],gamma[2*i+1],kappa[i],invmag[i]);*/
	  }
  }
  time(&t1);
	  //printf("rayshooter: time for calculating %i points %f sec\n            %f sec per point\n",Nrays
//		  ,difftime(t1,to),difftime(t1,to)/Nrays);

}
