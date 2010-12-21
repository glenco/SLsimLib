#include <math.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include "../../Library/Recipes/nrutil.h"
#include "../../Library/Recipes/nrutil.c"
/*#include "../Include/nr.h"*/
#include "../../Library/Recipes/ran2.c"
#include "../../Library/RecipesD/locateD.c"
#include "../../Library/RecipesD/bsstepD.c"
#include "../../Library/RecipesD/polintD.c"
#include "../../Library/RecipesD/dfridrD.c"
#include "../../Library/RecipesD/odeintD.c"
#include "../../Library/RecipesD/mmidD.c"
#include "../../Library/RecipesD/pzextrD.c"


#include "../../Library/cosmo.h"
#include "../../Library/cosmo.c"
#include "../../Library/powerCDM.c"

#include "../TreeCode_link/Tree.h"
//#include "../TreeCode_link/double_sort.c"
//#include "../TreeCode_link/Tree.c"
//#include "../TreeCode_link/TreeDriver.c"
//#include "../TreeCode_link/image_finder.c"

#include "../TreeCode/TreeNB.h"
//#include "rayshooter_stars.c"
//#include "../TreeCode/rayshooter.c"

int main(int arg,char **argv){
  TreeHndl i_tree,s_tree;
  Point *i_points,*s_points;
  double range,source[2],r_source,center[2],linkinglength,area,step;
  unsigned long i,j,k,m,Ngrid,Nimages,NImagePoints;
  ImageInfo *imageinfo;
  FILE *file;
  long Nframes;

  // set number of threads to be used in parallel sections
  omp_set_num_threads(4);

  Nframes=(long)(atol(argv[1]));
  file=fopen(argv[2],"w");
  fwrite(&Nframes,sizeof(long),1,file);
  fclose(file);

  step=-1.0e-8/2/2;
  source[0]=fabs(step)*Nframes/2;
  source[1]=0.0;

  imageinfo=NewImageInfo(300);

  //*** make initial grid
  Ngrid=50;
  range=5.0e-5;
  printf("grid range=%e Mpc\n",range);
  center[0]=0.0; center[1]=0.0;

  printf("finding images\n");
  for(i=0;i<Nframes;++i){
	printf("frame %i\n",i);

    // make new course grid
    i_points=NewPointArray(Ngrid*Ngrid,True);
    xygridpoints(i_points,range,center,Ngrid,0);
    s_points=LinkToSourcePoints(i_points,Ngrid*Ngrid);

    rayshooterInternal(Ngrid*Ngrid,i_points,NULL,False);

    // build trees
    i_tree=BuildTree(i_points,Ngrid*Ngrid);
    s_tree=BuildTree(s_points,Ngrid*Ngrid);

    // refine grid to smallest source size to find holes in image
//	r_source=2*5.0e-8/4;
//	linkinglength=1.1*sqrt(2.)*range/(Ngrid-1);
//	find_images(source,r_source,linkinglength,s_tree,i_tree,&Nimages
//			,imageinfo,&NImagePoints,0,False,False);

	r_source=1.0e-7;
    for(m=0;m<3;++m){ /* source sizes */

		linkinglength=1.1*sqrt(2.)*range/(Ngrid-1);
    	find_images(source,r_source,linkinglength,s_tree,i_tree,&Nimages
    			,imageinfo,&NImagePoints,0,False,False);

    	/********************/
    	/* output to files */
    	/*******************/
		file=fopen(argv[2],"a");

      /* print image plane grid */
		/*       fwrite(&(i_tree->pointlist->Npoints),sizeof(unsigned long),1,file); */
		/*       MoveToTopList(i_tree->pointlist); */
		/*       for(j=0;j<i_tree->pointlist->Npoints;++j){ */
		/* 	fwrite(i_tree->pointlist->current->x,sizeof(double),2,file); */
		/* 	MoveDownList(i_tree->pointlist); */
		/*       } */
		/* print source plane grid */
		/*       MoveToTopList(s_tree->pointlist); */
		/*       for(j=0;j<s_tree->pointlist->Npoints;++j){ */
		/* 	fwrite(s_tree->pointlist->current->x,sizeof(double),2,file); */
		/* 	MoveDownList(s_tree->pointlist); */
		/*       } */
		/* print magnification on source plane */
		/*       MoveToTopList(s_tree->pointlist); */
		/*       for(j=0;j<s_tree->pointlist->Npoints;++j){ */
		/* 	fwrite(&(s_tree->pointlist->current->invmag),sizeof(double),1,file); */
		/* 	MoveDownList(s_tree->pointlist); */
		/*       } */

      /* print image points */
    	fwrite(&NImagePoints,sizeof(unsigned long),1,file);
    	for(k=0;k<Nimages;++k){
    		for(j=0;j<imageinfo[k].Npoints;++j){
    			fwrite(imageinfo[k].points[j].x,sizeof(double),2,file);
    		}
    	}
		for(k=0;k<Nimages;++k){
			for(j=0;j<imageinfo[k].Npoints;++j){
				fwrite(&(imageinfo[k].points[j].gridsize),sizeof(double),1,file);
			}
		}

      /* print image magnifications */
      fwrite(&Nimages,sizeof(unsigned long),1,file);
      for(k=0;k<Nimages;++k){
    	  area=imageinfo[k].area;
    	  area/=pi*r_source*r_source;
    	  fwrite(&area,sizeof(double),1,file);
      }
    
      fclose(file);

      /*****************************/

      r_source/=2.0;
    }
    source[0]+=step;
    source[1]=0.0;

    freeTree(i_tree);
    freeTree(s_tree);
  }

/*   PrintList(i_tree->pointlist); */
/*   PrintList(s_tree->pointlist); */
 /*   PrintKist(imageinfo[5].innerborderkist); */
/*    PrintList(imageinfo[5].outerborder); */

/*   printf("\n%i\n",NImagePoints);   */
/*   for(i=0;i<Nimages;++i){ */
/*     for(j=0;j<imageinfo[i].Npoints;++j){ */
/*       printf("%e  %e  %e\n",imageinfo[i].points[j].x[0] */
/* 	     ,imageinfo[i].points[j].x[1],imageinfo[i].points[j].invmag); */
/*     } */
/*   } */

  printf("\nNImagePoints = %i Nimages = %i\n",NImagePoints,Nimages);

  return 1;
}

