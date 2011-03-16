
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <nr.h>
#include "../../Library/Recipes/nrutil.h"
#include "../../Library/Recipes/nrutil.c"
#include "../../Library/Recipes/ran2.c"
#include <nrD.h>
#include "../../Library/RecipesD/locateD.c"
#include "../../Library/Recipes/ran2.c"
#include "../../Library/RecipesD/powellD.c"
#include "../../Library/RecipesD/odeintD.c"
#include "../../Library/RecipesD/bsstepD.c"
#include "../../Library/RecipesD/mmidD.c"
#include "../../Library/RecipesD/pzextrD.c"
#include "../../Library/RecipesD/polintD.c"
#include "../../Library/RecipesD/dfridrD.c"

#include "../../Library/cosmo.h"
#include "../../Library/cosmo.c"
#include "../../Library/powerCDM.c"

#include "../TreeCode_link/Tree.h"
//#include "../TreeCode_link/Tree.c"
//#include "../TreeCode_link/double_sort.c"
//#include "../TreeCode_link/TreeDriver.c"
//#include "../TreeCode_link/image_finder.c"

#include "../TreeCode/TreeNB.h"/**/
//#include "../TreeCode/rayshooter.c"/**/

char *paramfile;

int main(int arg,char **argv){
  TreeHndl i_tree,s_tree;
  Point *i_points,*s_points,*source_centers,*image_centers;
  double range,r_source,center[2],linkinglength,rmin,rmax,*times,x[2];
  unsigned long i,j,k,Ngrid,Nimages,NImagePoints,Ngrid_sources;
  ImageInfo *imageinfo;
  time_t to,t1,t2,now,t3,t4;
  Boolean verbose,nocrit;

  verbose=False;
  //verbose=True;
  nocrit=True;

  if(arg==2) paramfile=argv[1];
  else paramfile=NULL;

  time(&to);
  range=0.3;
  Ngrid=21;
  Ngrid_sources=10;

  /* select source */
  r_source=1.0e-5;

  // print out grid of for testing rayshooting
/*   double alpha[2],gamma[2],kappa,invmag,r; */
/*   for(r=0.0001;r<0.1;r*=pow(10,1.0/30)){ */
/*     for(tmp=0;tmp<2*pi;tmp+=2*pi/10){ */
/*       x[0]=r*cos(tmp); */
/*       x[1]=r*sin(tmp); */
/*       rayshooterInternal(x,alpha,gamma,&kappa,&invmag); */
/*       printf("%e  %e %e  %e %e  %e %e  %e  %e\n",r,x[0],x[1],alpha[0],alpha[1] */
/* 	     ,gamma[0],gamma[1],kappa,invmag); */
/*     } */
/*   } */
/*   exit(0); */

  center[0]=0.0; center[1]=0.0;

  // define source positions on polar grid
  //   make grid

  rmax=3.0e-2;
  rmin=1.0e-4;

  rmax=1.0e-3;
  rmin=3.0e-2;

  image_centers=NewPointArray(Ngrid_sources*Ngrid_sources,True);
  log_polar_grid(image_centers,rmax,rmin,center,Ngrid_sources);
  source_centers=LinkToSourcePoints(image_centers,Ngrid_sources*Ngrid_sources);
  time(&to);
  rayshooterInternal(Ngrid_sources*Ngrid_sources,image_centers,NULL,False);
  time(&t1);
  printf("time for calculating %i points %f sec\n       %f sec per point\n"
		  ,Ngrid_sources*Ngrid_sources
		  ,difftime(t1,to),difftime(t1,to)/Ngrid_sources/Ngrid_sources);

  //   make initial grid
  time(&to);
  //initialize_grid(center,range,Ngrid,s_tree,i_tree);

  i_points=NewPointArray(Ngrid*Ngrid,True);
  xygridpoints(i_points,range,center,Ngrid,0);
  s_points=LinkToSourcePoints(i_points,Ngrid*Ngrid);
  rayshooterInternal(Ngrid*Ngrid,i_points,i_tree,False);
  // build trees
  i_tree=BuildTree(i_points,Ngrid*Ngrid);
  s_tree=BuildTree(s_points,Ngrid*Ngrid);

//  PrintList(i_tree->pointlist);
//  PrintList(s_tree->pointlist);

  time(&now);
  if(verbose) printf("set up initial grid in %f min\n     %f sec per point\n\n"
	 ,difftime(now,to)/60.,difftime(now,to)/Ngrid/Ngrid);

  /* find critical curves */
  time(&t1);

  imageinfo=NewImageInfo(200);

  time(&t2);

  times=(double *) malloc(Ngrid_sources*Ngrid_sources*sizeof(double));
  printf("\n%li\n",Ngrid_sources*Ngrid_sources);
  for(i=0;i<Ngrid_sources*Ngrid_sources;++i){

      printf("\n%li  %e %e\n",i
    		  ,source_centers[i].x[0],source_centers[i].x[1]);

      printf("%e %e\n",image_centers[i].x[0],image_centers[i].x[1]);
      // find images
      if(verbose) printf("finding images\n");

      // refine grid
      linkinglength=0.001;

      time(&t3);
      find_images(source_centers[i].x,r_source,s_tree,i_tree,&Nimages
    		  ,imageinfo,&NImagePoints,0,True,2,verbose,False);

      time(&t4);
      times[i]=difftime(t4,t3);

   /*****************************/
    /*printf("end of refinements\n %e %e  \n",source_centers[i].x[0],source_centers[i].x[1]);*/
    /* print centers of images and magnifications */
    /*printf(" image centers and magnifications\n");*/

    // find the point source magnification for each image
    if(verbose) printf("%li images\n",Nimages);
    else printf("%li\n",Nimages);

    for(k=0,rmin=1.0e99;k<Nimages;++k){
			//printf("%i\n",imageinfo[k].Npoints);
    	x[0]=x[1]=0;
    	for(j=0;j<imageinfo[k].Npoints;++j){
    		x[0]+=imageinfo[k].points[j].x[0]/imageinfo[k].Npoints;
    		x[1]+=imageinfo[k].points[j].x[1]/imageinfo[k].Npoints;
    	}
    	printf("%li  %e %e  %e\n",k,x[0],x[1],imageinfo[k].area/(pi*r_source*r_source));
    }

    /* print all points in images */
    /*
	if(!verbose){
		for(k=0;k<Nimages;++k){
			printf("%i\n",imageinfo[k].Npoints);
			for(j=0;j<imageinfo[k].Npoints;++j){
				printf(" %e %e %e\n",imageinfo[k].points[j].x[0]
				       ,imageinfo[k].points[j].x[1],imageinfo[k].points[j].gridsize);
			}
			PrintList(imageinfo[k].outerborder);
		}
	}
	*/
  }

  printf("\n%li\n",NImagePoints);

//  for(i=0;i<Nimages;++i){
//	printf("%i\n",imageinfo[i].Npoints);
//    for(j=0;j<imageinfo[i].Npoints;++j){
//      printf("%e  %e  %e\n",imageinfo[i].points[j].x[0],imageinfo[i].points[j].x[1],imageinfo[i].points[j].invmag);
//    }
//  }

  printf("\nNImagePoints = %li Nimages = %li\n",NImagePoints,Nimages);

  time(&now);
  printf("found images in %f min\n     %f sec per source\n    %f sec per point\n"
	 ,difftime(now,t2)/60.
	 ,difftime(now,t2)/Ngrid_sources/Ngrid_sources
	 ,difftime(now,to)/i_tree->pointlist->Npoints
	 );
  printf("times for each source\n");
  for(i=0;i<Ngrid_sources*Ngrid_sources;++i) printf("    %f sec\n",times[i]);

  printf("\nNumber of rays shot: %li\n",i_tree->pointlist->Npoints);

  return 1;
}
