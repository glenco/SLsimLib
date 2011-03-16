#include <math.h>
#include <stdio.h>
#include <time.h>
#include <nr.h>
#include "../../Library/Recipes/nrutil.h"
#include "../../Library/Recipes/nrutil.c"
#include "../../Library/Recipes/ran2.c"
#include <nrD.h>
#include "../../Library/RecipesD/locateD.c"
#include "../../Library/RecipesD/ran2D.c"
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

int main(int arg,char **argv){
  TreeHndl i_tree,s_tree;
  Point *i_points,*s_points,*point;
  double range,source[2],r_source,center[2],linkinglength,tmp,x[2],xrange[2],yrange[2],rmin,r,critrange[2],*times;
  unsigned long i,j,k,Ngrid,Nimages,NImagePoints;
  long seed=172930;
  ImageInfo *imageinfo,*critical;
  FILE *file;
  int Ncrit,tange_caustic,Nsources;
  time_t to,t1,t2,now,t3,t4;
  Boolean verbose;

  verbose=False;
  //verbose=True;

  time(&to);
  range=0.3;
  Ngrid=21;
  Nsources=1;

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

  times=(double *) malloc(Nsources*sizeof(double));

  //   make initial grid
  center[0]=0.0; center[1]=0.0;

  time(&to);
  //initialize_grid(center,range,Ngrid,s_tree,i_tree);

  i_points=NewPointArray(Ngrid*Ngrid,True);
  xygridpoints(i_points,range,center,Ngrid,0);
  s_points=LinkToSourcePoints(i_points,Ngrid*Ngrid);
  rayshooterInternal(Ngrid*Ngrid,i_points,i_tree,False);
  // build trees
  i_tree=BuildTree(i_points,Ngrid*Ngrid);
  s_tree=BuildTree(s_points,Ngrid*Ngrid);

  time(&now);
  if(verbose) printf("set up initial grid in %f min\n     %f sec per point\n\n"
	 ,difftime(now,to)/60.,difftime(now,to)/Ngrid/Ngrid);

  /* find critical curves */
  time(&t1);

   // critical=find_crit(s_tree,i_tree,&Ncrit,5.0e-5);
  critical=find_crit(s_tree,i_tree,&Ncrit,9.0e-5);

  time(&now);
  if(verbose) printf("found caustics in %f min\n",difftime(now,t1)/60.);
  if(verbose) printf(" %f sec per point \n",difftime(now,t1)/(i_tree->pointlist->Npoints-Ngrid*Ngrid));

  printf(" output\n");

  // output critical curves
  if(!verbose){
	  printf("%i\n",Ncrit);
	  for(i=0;i<Ncrit;++i){
		  printf("%i\n",critical[i].Npoints);
		  for(j=0;j<critical[i].Npoints;++j){
			  printf("%e %e\n",critical[i].points[j].x[0]
			         ,critical[i].points[j].x[1]);
		  }
	  }
  }
  // output caustic curves
  printf("%i\n",Ncrit);
  for(i=0;i<Ncrit;++i){
    printf("%i\n",critical[i].Npoints);
    for(j=0;j<critical[i].Npoints;++j){
      printf("%e %e\n",critical[i].points[j].image->x[0]
	     ,critical[i].points[j].image->x[1]);
    }
  }

  /** find range in tangential caustics, uses caustic of largest area critical curve  **/

  for(i=0,tmp=0;i<Ncrit;++i){
    if(tmp< critical[i].area){
      tmp = critical[i].area;
      tange_caustic=i;
    }
  }

  if(verbose) for(i=0;i<Ncrit;++i)  printf("caustic %i area %e\n",i,critical[i].area);

  xrange[0]=xrange[1]=critical[tange_caustic].points[0].image->x[0];
  yrange[0]=yrange[1]=critical[tange_caustic].points[0].image->x[1];

  critrange[0]=critrange[1]=critical[tange_caustic].points[0].x[0];

  for(j=0;j<critical[tange_caustic].Npoints;++j){

	if(critical[tange_caustic].points[j].image->x[0] < xrange[0])
      xrange[0]=critical[tange_caustic].points[j].image->x[0];
    if(critical[tange_caustic].points[j].image->x[0] > xrange[1]) 
      xrange[1]=critical[tange_caustic].points[j].image->x[0];
    
    if(critical[tange_caustic].points[j].image->x[1] < yrange[0]) 
      yrange[0]=critical[tange_caustic].points[j].image->x[1];
    if(critical[tange_caustic].points[j].image->x[1] > yrange[1]) 
      yrange[1]=critical[tange_caustic].points[j].image->x[1];

    if(critical[tange_caustic].points[j].x[0] < critrange[0]) 
      critrange[0]=critical[tange_caustic].points[j].x[0];
    if(critical[tange_caustic].points[j].x[0] > critrange[1]) 
      critrange[1]=critical[tange_caustic].points[j].x[0];
  }

  if(tmp==0) printf("error: no tangential caustic found\n");
  if(verbose) printf(" xrange = %e %e\n  yrange = %e %e\n",xrange[0],xrange[1],yrange[0],yrange[1]);

  imageinfo=NewImageInfo(200);

  /* select source */
  r_source=1.0e-5;

  time(&t2);

  printf("\n%i\n",Nsources);
  for(i=0;i<Nsources;++i){
    do{
      source[0]=xrange[0] + ran2D(&seed)*(xrange[1]-xrange[0]);
      source[1]=yrange[0] + ran2D(&seed)*(yrange[1]-yrange[0]);
      /*printf("source = %e %e\n",source[0],source[1]);*/
    }while(abs(windings(source,critical[tange_caustic].points,critical[tange_caustic].Npoints,&tmp,1)) < 1);

    printf("\n%i  %e %e\n",i,source[0],source[1]);

    // find images
    if(verbose) printf("finding images\n");

    // refine grid
    linkinglength=fabs(critrange[1]-critrange[0])/10;

    time(&t3);
    find_images(source,r_source,linkinglength,s_tree,i_tree,&Nimages
    		   ,imageinfo,&NImagePoints,fabs(critrange[1]-critrange[0])/3,verbose);
    time(&t4);
    times[i]=difftime(t4,t3);

    MoveToTopList(i_tree->pointlist);
    printf("%i\n",i_tree->pointlist->Npoints);
    for(k=0;k<i_tree->pointlist->Npoints;++k){
    	printf("%e  %e   %e  %e\n",i_tree->pointlist->current->x[0],i_tree->pointlist->current->x[1]
    	                   ,i_tree->pointlist->current->kappa,i_tree->pointlist->current->invmag);
    	MoveDownList(i_tree->pointlist);
    }
//    PrintList(i_tree->pointlist);
    PrintList(s_tree->pointlist);

    /*****************************/
    /*printf("end of refinements\n %e %e  \n",source[0],source[1]);*/
    /* print centers of images and magnifications */
    /*printf(" image centers and magnifications\n");*/

    // find the point source magnification for each image
    printf("%i images\n",Nimages);
    for(k=0;k<Nimages;++k){
    	//printf("%i\n",imageinfo[k].Npoints);
      x[0]=0.0; x[1]=0.0;
      for(j=0,rmin=1.0e99;j<imageinfo[k].Npoints;++j){
    	  r=sqrt( pow(imageinfo[k].points[j].image->x[0] - source[0],2)
    			  + pow(imageinfo[k].points[j].image->x[1] - source[1],2) );
    	  if(rmin>r){
				 rmin=r;
				 point=&(imageinfo[k].points[j]);
    	  }
    	  x[0]+=imageinfo[k].points[j].x[0]/imageinfo[k].Npoints;
    	  x[1]+=imageinfo[k].points[j].x[1]/imageinfo[k].Npoints;
      }
      printf("%e %e  %e %e\n",x[0],x[1],imageinfo[k].area/(pi*r_source*r_source)
	     ,1.0/(pow(1-point->kappa,2) - pow(point->image->gamma[0],2) - pow(point->image->gamma[1],2)) );
    }

    /* print all points in images */
	if(!verbose){
		for(k=0;k<Nimages;++k){
			printf("%i\n",imageinfo[k].Npoints);
			for(j=0;j<imageinfo[k].Npoints;++j){
				printf(" %e %e\n",imageinfo[k].points[j].x[0]
				                 ,imageinfo[k].points[j].x[1]);
			}
			PrintList(imageinfo[k].outerborder);
		}
	}
  }

  // find caustic again
/*  critical=find_crit(s_tree,i_tree,&Ncrit,1.0e-4);

  // output critical curves
  printf("%i\n",Ncrit);
  for(i=0;i<Ncrit;++i){
    printf("%i\n",critical[i].Npoints);
    for(j=0;j<critical[i].Npoints;++j){
      printf("%e %e\n",critical[i].points[j].x[0]
	     ,critical[i].points[j].x[1]);
    }
  }

  // output cautic curves
  printf("%i\n",Ncrit);
  for(i=0;i<Ncrit;++i){
    printf("%i\n",critical[i].Npoints);
    for(j=0;j<critical[i].Npoints;++j){
      printf("%e %e\n",critical[i].points[j].image->x[0]
	     ,critical[i].points[j].image->x[1]);
    }
  }
*/

  printf("\n%i\n",NImagePoints);

//  for(i=0;i<Nimages;++i){
//	printf("%i\n",imageinfo[i].Npoints);
//    for(j=0;j<imageinfo[i].Npoints;++j){
//      printf("%e  %e  %e\n",imageinfo[i].points[j].x[0],imageinfo[i].points[j].x[1],imageinfo[i].points[j].invmag);
//    }
//  }

  printf("\nNImagePoints = %i Nimages = %i\n",NImagePoints,Nimages);

  time(&now);
  printf("found images in %f min\n     %f sec per source\n    %f sec per point\n"
	 ,difftime(now,t2)/60.
	 ,difftime(now,t2)/Nsources
	 ,difftime(now,to)/i_tree->pointlist->Npoints
	 );
  printf("times for each source\n");
  for(i=0;i<Nsources;++i) printf("    %f sec\n",times[i]);

  printf("\nNumber of rays shot: %i\n",i_tree->pointlist->Npoints);


  exit(0);
}

void rayshooterInternal(unsigned long Npoints,Point *i_points,TreeHndl i_tree){
  /* i_points need to be already linked to s_points */
  static char paramfile[50];
  static short i=0;

  if(i==0){
	  sprintf(paramfile,"ParameterFiles/paramfileSIE1M_vol_512");
	  ++i;
  }
  /*rayshooter(x,1,alpha,gamma,kappa,invmag,paramfile);*/
  rayshooter(Npoints,i_points,i_tree,paramfile);
}

void in_source(double *y_source,ListHndl sourcelist){
  return;
}

