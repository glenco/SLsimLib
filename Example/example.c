
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
#include "../../Library/RecipesD/powellD.c"
#include "../../Library/RecipesD/odeintD.c"
#include "../../Library/RecipesD/bsstepD.c"
#include "../../Library/RecipesD/mmidD.c"
#include "../../Library/RecipesD/pzextrD.c"
#include "../../Library/RecipesD/polintD.c"
#include "../../Library/Recipes/poidev.c"
#include "../../Library/Recipes/gammln.c"
#include "../../Library/RecipesD/dfridrD.c"
#include "../../Library/Recipes/gasdev.c"

#include "../../Library/cosmo.h"
#include "../../Library/powerCDM.c"
#include "../../Library/cosmo.c"

#include "../TreeCode_link/Tree.h"
#include "../TreeCode/TreeNB.h"
#include "../AnalyticNSIE/analytic_lens.h"

char *paramfile,*outputfile;

int main(int arg,char **argv){
  TreeHndl i_tree,s_tree;
  Point *i_points,*s_points;
  double range,source[2],center[2],z;
  unsigned long i,j,Ngrid,NImagePoints;
  int Nimages;
  ImageInfo *imageinfo,*critical;
  //FILE *file;
  int Ncrit;
  time_t to;
  Boolean verbose,success;

  verbose=False;
  //verbose=True;


  // if parameter file if not given on execution a default parameter file is used
  if(arg==2) paramfile=argv[1];
  else paramfile=NULL;

  time(&to);
  range=4.*0.25;  // size of initial grid
  Ngrid=64;   // 1-d number of grid points

  /* this data structure holds all the information
   *   about the images
   */
  imageinfo=NewImageInfo(200);


  //   make initial grid
  /****************************/

  center[0]=0.0; center[1]=0.0;                          // center of grid
  i_points=NewPointArray(Ngrid*Ngrid,True);                   // make a new array for image points
  xygridpoints(i_points,range,center,Ngrid,0);           // lay down grid points
  s_points=LinkToSourcePoints(i_points,Ngrid*Ngrid);     // make and link source points
  i_tree=NULL;                                           // to make sure interpolation isn't used
  rayshooterInternal(Ngrid*Ngrid,i_points,i_tree,False); // lens calculate the deflection angles & shear etc

  /*
   * the parameters for the lens are read on the first use of rayshooterInternal
   *   they can be changed by hand in the structure lens
   */

  // build trees
  i_tree=BuildTree(i_points,Ngrid*Ngrid);              //  build tree with image points
  s_tree=BuildTree(s_points,Ngrid*Ngrid);              //  build tree with source points

  // print some information about the lens
  //printf("%.2f  %.3f  %.3f  %.3f\n",lens->sigma,lens->axis_ratio,lens->zlens,lens->zsource);

  /*
   *   caustics and critical curves
   */

  critical=find_crit(s_tree,i_tree,&Ncrit,5.0e-3,&success,False,verbose);  // find critical curves and caustics

  // print critical curves
   printf("%i\n",Ncrit);
   for(j=0;j<Ncrit;++j){
 	  printf("%li\n",critical[j].Npoints);
 	  for(i=0;i<critical[j].Npoints;++i){
 		  printf("%e %e\n",critical[j].points[i].x[0],critical[j].points[i].x[1]);
 	  }
   }
   // print caustic curves
    printf("%i\n",Ncrit);
    for(j=0;j<Ncrit;++j){
  	  printf("%li\n",critical[j].Npoints);
  	  for(i=0;i<critical[j].Npoints;++i){
  		  printf("%e %e\n",critical[j].points[i].image->x[0],critical[j].points[i].image->x[1]);
  	  }
    }

    /*
     * set a source
     */

    source[0] = 1.0e-3;    // set source position, same units
    source[1] = 1.3e-3;

    source[0] = 0.0908*1.05374078;    // set source position, same units
    source[1] =-0.0565*1.05374078;

    /*
     * find images
     */

   find_images(source,1.05*4.0e-3*1.05374078,s_tree,i_tree,&Nimages
 			  ,imageinfo,&NImagePoints,0,False,1,verbose,True);

    // print out image points

    printf("%li\n",Nimages);
    for(i=0;i<Nimages;++i){
    	printf("%li\n",imageinfo[i].Npoints);
    	for(j=0;j<imageinfo[i].Npoints;++j){
    		// print the position and cell size
    		printf("%e %e  %e\n",imageinfo[i].points[j].x[0],imageinfo[i].points[j].x[1],imageinfo[i].points[j].gridsize);
    	}
    }
    
    printf("%li\n",Nimages);
    for(i=0;i<Nimages;++i){
    	printf("%li\n",imageinfo[i].Npoints);
    	for(j=0;j<imageinfo[i].Npoints;++j){
    		// print the position and cell size
    		printf("%e %e  %e\n",imageinfo[i].points[j].image->x[0],imageinfo[i].points[j].image->x[1]
    		                                                    ,imageinfo[i].points[j].image->gridsize);
    	}
    }

    PrintList(i_tree->pointlist);
    //PrintAnaLens(lens,True,False);  // print out lens parameters


    // free memory
    free(critical->points);
    freeImageInfo(critical,Ncrit);

    free(imageinfo->points);
    freeImageInfo(imageinfo,200);
    freeTree(i_tree);
    freeTree(s_tree);

    return 1;
}
