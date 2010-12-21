
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include "../../Library/Recipes/nr.h"
#include "../../Library/Recipes/nrutil.h"
#include "../../Library/Recipes/nrutil.c"
#include "../../Library/Recipes/ran2.c"
#include "../../Library/RecipesD/nrD.h"
#include "../../Library/RecipesD/locateD.c"
#include "../../Library/RecipesD/ran2D.c"
#include "../../Library/RecipesD/powellD.c"
#include "../../Library/RecipesD/odeintD.c"
#include "../../Library/RecipesD/bsstepD.c"
#include "../../Library/RecipesD/mmidD.c"
#include "../../Library/RecipesD/pzextrD.c"
#include "../../Library/RecipesD/polintD.c"
#include "../../Library/RecipesD/dfridrD.c"
#include "../../Library/Recipes/gasdev.c"

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

#include "../ImageProcessing/image_processing.h"
#include "../Galaxies/galaxies.h"

#define MaxNimages 10

char *paramfile,*outputfile;;

int main(int arg,char **argv){
  TreeHndl i_tree,s_tree;
  ListHndl plist;
  Point *i_points,*s_points,*caustic_points,*tmp_point;
  double range,source[2],r_source,center[2],tmp,xrange[2],yrange[2],*image;
  double **x_centers,pixelsize,pixelrange,position_angle;
  unsigned long i,j,Ngrid,Nimages,Nimages_tmp,NImagePoints,NImagePoints_tmp
	  ,Ncritpoints;
  unsigned long Npixels=0;
  long seed=182993;
  ImageInfo *imageinfo,*critical;
  FILE *file;
  int Ncrit,tange_caustic,Nsources;
  Boolean verbose,nocrit,just_mags;

  verbose=False;
  //verbose=True;
  nocrit=False;

  just_mags=True;

  if(arg==2) paramfile=argv[1];
  else paramfile=NULL;


  range=0.3;
  Ngrid=64;

  tmp_point=NewPointArray(1,True);
  tmp_point->image=NewPointArray(1,True);

  //   make initial grid
  center[0]=0.0; center[1]=0.0;

  i_points=NewPointArray(Ngrid*Ngrid,True);
  xygridpoints(i_points,range,center,Ngrid,0);
  printf("Ngrid=%i\n",Ngrid);
  s_points=LinkToSourcePoints(i_points,Ngrid*Ngrid);
  rayshooterInternal(Ngrid*Ngrid,i_points,i_tree,False);
  // build trees
  i_tree=BuildTree(i_points,Ngrid*Ngrid);
  s_tree=BuildTree(s_points,Ngrid*Ngrid);

  /*
  PrintPoint(s_points);
   PrintPoint(&s_points[1]);
   exit(0);
   */

  /* find critical curves */

  if(!nocrit){
    critical=find_crit(s_tree,i_tree,&Ncrit,3.0e-5,verbose);
    printf("\n\n");

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

    /* select critical curve with largest area
     *     this will generally be the tangential critical curve if
     */
    for(i=0,tmp=0,tange_caustic=0;i<Ncrit;++i){
      if(tmp< critical[i].area){
    	  tmp = critical[i].area;
    	  tange_caustic=i;
      }
    }

   //printf("%e\n",critical[tange_caustic].area);
   if(verbose) for(i=0;i<Ncrit;++i)  printf("caustic %i area %e\n",i,critical[i].area);


   // re-assign caustic curve
   Ncritpoints=critical[tange_caustic].Npoints;
   caustic_points=NewPointArray(Ncritpoints,True);
   for(j=0;j<Ncritpoints;++j){
	   caustic_points[j].x[0]=critical[tange_caustic].points[j].image->x[0];
	   caustic_points[j].x[1]=critical[tange_caustic].points[j].image->x[1];
   }

   xrange[0]=xrange[1]=caustic_points->x[0];
   yrange[0]=yrange[1]=caustic_points->x[1];

   for(j=0;j<critical[tange_caustic].Npoints;++j){
	   if(caustic_points[j].x[0] < xrange[0])
		   xrange[0]=caustic_points[j].x[0];
	   if(caustic_points[j].x[0] > xrange[1])
		   xrange[1]=caustic_points[j].x[0];
    
	   if(caustic_points[j].x[1] < yrange[0])
		   yrange[0]=caustic_points[j].x[1];
	   if(caustic_points[j].x[1] > yrange[1])
		   yrange[1]=caustic_points[j].x[1];
   }

   if(tmp==0) printf("error: no tangential caustic found\n");
   if(verbose) printf(" xrange = %e %e\n  yrange = %e %e\n",xrange[0],xrange[1],yrange[0],yrange[1]);
  }else printf("0\n");

  do{
	  source[0]=xrange[0] + ran2D(&seed)*(xrange[1]-xrange[0]);
	  source[1]=yrange[0] + ran2D(&seed)*(yrange[1]-yrange[0]);
	  //printf("source = %e %e\n",source[0],source[1]);
  }while(abs(windings(source,caustic_points,Ncritpoints,&tmp,0)) < 1 );
		  //&& sqrt(pow(source[0],2)+pow(source[1],2)) < 3.0e-4);

  source[0]*=3;
  source[0]=1.672116e-04; source[1]=-5.051015e-04;
  // create galaxy source
  Nsources=600;
  imageinfo=NewImageInfo(MaxNimages);
  x_centers=dmatrix(0,Nsources-1,0,1);
  position_angle=30*pi/180;
  //printf("galaxy center: %e %e\n",source[0],source[1]); exit(1);
  create_sersic(1,1.0e-3,0.5,source,position_angle,x_centers,Nsources);
  source[0]/=3;
  source[0]+=1.0e-3*cos(position_angle+50*pi/180);
  source[1]+=1.0e-3*sin(position_angle+50*pi/180);
  create_sersic(1,7.0e-5,1.0,source,position_angle,&x_centers[Nsources - 1 - Nsources/5],Nsources/5);

  r_source=1.0e-4;

  file=fopen(outputfile,"w");

  // setup pixel map
  //printf("pixelizing\n");
  pixelsize=1.0e-4;
  pixelrange=1.5e-2;
  Npixels=pow(2, (int)(log(pixelrange/pixelsize)/log(2)) + 1.0);
  pixelrange=(Npixels-1)*pixelsize;

  image=(double *)calloc(Npixels*Npixels,sizeof(double));

  // calculate the area within the caustic
  windings(caustic_points->x,caustic_points,Ncritpoints,&tmp,0);
  //fprintf(file,"%e\n",tmp);
  printf("%e\n",tmp);

  //fprintf(file,"\n%i\n",Nsources);
  printf("\n%i\n",Nsources);
  for(i=0;i<Nsources;++i){
	  //fprintf(file,"%i  %e %e  %e\n",i,x_centers[i][0],x_centers[i][1],r_source);
	  printf("%i  %e %e  %e\n",i,x_centers[i][0],x_centers[i][1],r_source);
  }

  for(i=0,Nimages=0,NImagePoints=0;i<Nsources;++i){
	  // find images
	  find_images(x_centers[i],r_source,s_tree,i_tree,&Nimages_tmp
			  ,imageinfo,&NImagePoints_tmp,0,True,1,verbose);

	  pixelize(image,Npixels,pixelrange,center,imageinfo,Nimages_tmp,True,False);

	  Nimages += Nimages_tmp;
	  NImagePoints += NImagePoints_tmp;
	  //printf("i=%i\n",i);
  }
/*
  printf("%i\n",NImagePoints);
  for(i=0;i<Nimages;++i){
	  //printf("%i\n",imageinfo[i].Npoints);
	  for(j=0;j<imageinfo[i].Npoints;++j){
		  printf("%e  %e  %e\n",imageinfo[i].points[j].x[0],imageinfo[i].points[j].x[1]
		                                     ,imageinfo[i].points[j].gridsize);
	  }
  }
*/

  printf("%e\n%i  %i\n",pixelsize,Npixels,Npixels);
  for(i=0;i<Npixels*Npixels;++i) printf("%e\n",image[i]);

  //PrintList(plist);

  printf("\nNImagePoints = %i Nimages = %i\n",NImagePoints,Nimages);
  printf("\nNumber of rays shot: %i\n",i_tree->pointlist->Npoints);
  printf("param file : %s\n",paramfile);

  free(imageinfo->points);
  freeImageInfo(imageinfo,MaxNimages);
  FreePointArray(caustic_points);
  FreePointArray(tmp_point->image);
  FreePointArray(tmp_point);
  freeTree(i_tree);
  freeTree(s_tree);
  if(!nocrit){
	  free(critical->points);
	  freeImageInfo(critical,Ncrit);
  }
  free(image);
  free_dmatrix(x_centers,0,Nsources-1,0,1);

  fclose(file);
}
