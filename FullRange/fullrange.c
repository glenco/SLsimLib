
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <string.h>
#include "../../Library/Recipes/nr.h"
#include "../../Library/Recipes/nrutil.h"
#include "../../Library/Recipes/nrutil.c"
#include "../../Library/Recipes/ran2.c"
#include "../../Library/RecipesD/nrD.h"
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

#include "../AnalyticNSIE/analytic_lens.h"
#include "../TreeCode_link/Tree.h"
#include "../TreeCode/TreeNB.h"

//char *paramfile,*outputfile;
AnaLens *lens=0;

/*
 *   fullrange.c calculates the magnifications of circular images as a function of
 *   source size including substructures and microlensing by stars
 */

int main(int arg,char **argv){
  TreeHndl i_tree,s_tree;
  Point *i_points,*s_points,*caustic_points,*tmp_point,*point=0,*centers;
  double range,r_source_phys,center[2]
                          ,xrange[2],yrange[2],rmin,r,*times;
  unsigned long i,j,k,m,Ngrid,NImagePoints,Ncritpoints;
  int Nimages,NewNimages;
  long seed=282923;
  //long seed=28374;
  //time(&seed);
  const int NimageMax = 750;

  printf("seed = %li\n",seed);

  double min_image_seporation;
  ImageInfo *imageinfo,*critical;
  //FILE *file;
  int Ncrit,Nsources,Ntotal_lenses;
  time_t to,now,t4;
  Boolean verbose,nocrit,just_mags,success;
  char *paramfile;
  FILE *file;


  verbose=False;
  //verbose=True;
  nocrit=False;

  just_mags=True;

  if(arg==2) paramfile=argv[1];
  else paramfile=NULL;

  time(&to);
  range=0.5;
  Ngrid=64;
  Ntotal_lenses=1;
  //Ntotal_lenses=500;
  //Nsources=30;
  Nsources=1;
  // select source
  r_source_phys=1.0e-5;
  //r_source_phys=2.0e-6;
  min_image_seporation=0.05;  // in arcsec

  times=(double *) malloc(Nsources*sizeof(double));
  tmp_point=NewPointArray(1,True);
  tmp_point->image=NewPointArray(1,True);

  imageinfo=NewImageInfo(NimageMax);

  //   make initial grid
  center[0]=0.0; center[1]=0.0;

  time(&to);

   // read parameter file in
  if(paramfile==NULL){
	  paramfile=(char *)malloc(60*sizeof(char));
	  sprintf(paramfile,"FullRange/paramfile");
  }
  lens=(AnaLens *)calloc(1,sizeof(AnaLens));
  readparams_ana(paramfile,&cosmo,lens);

  //RandomizeSubstructure2(lens,2,&seed);

  // initialize tree and lens model
  i_points=NewPointArray(Ngrid*Ngrid,True);
  xygridpoints(i_points,range,center,Ngrid,0);
  s_points=LinkToSourcePoints(i_points,Ngrid*Ngrid);
  i_tree=NULL;
  rayshooterInternal(Ngrid*Ngrid,i_points,i_tree,False);

  //exit(0);
  // build trees
  i_tree=BuildTree(i_points,Ngrid*Ngrid);
  s_tree=BuildTree(s_points,Ngrid*Ngrid);

   // find points on critical curves

  critical=find_crit(s_tree,i_tree,&Ncrit,5.0e-5,&success,False,verbose);

  if(Ncrit == 0 || Ncrit > 1){
	  PrintAnaLens(lens,False,False);
	  ERROR_MESSAGE();
	  printf("ERROR: Ncrit=%i \n",Ncrit);
	  exit(0);
  }
  printf("\n\n");

   // re-assign caustic curve
  for(i=0,Ncritpoints=0;i<Ncrit;++i) Ncritpoints+=critical[i].Npoints;
  caustic_points=NewPointArray(Ncritpoints,True);
  for(i=0,j=0;i<Ncrit;++i){
	  for(k=0;k<critical[i].Npoints;++k){
		  caustic_points[j].x[0] = critical[i].points[k].image->x[0];
		  caustic_points[j].x[1] = critical[i].points[k].image->x[1];
		  ++j;
	  }
  }

  // find extreme values of caustic region
  xrange[0]=xrange[1]=caustic_points->x[0];
  yrange[0]=yrange[1]=caustic_points->x[1];
  for(j=0;j<Ncritpoints;++j){
	  if(caustic_points[j].x[0] < xrange[0])
		  xrange[0]=caustic_points[j].x[0];
	  if(caustic_points[j].x[0] > xrange[1])
		  xrange[1]=caustic_points[j].x[0];
    
	  if(caustic_points[j].x[1] < yrange[0])
		  yrange[0]=caustic_points[j].x[1];
	  if(caustic_points[j].x[1] > yrange[1])
		  yrange[1]=caustic_points[j].x[1];
  }
  xrange[0]*=1.01;
  xrange[1]*=1.01;
  yrange[0]*=1.01;
  yrange[1]*=1.01;

/*
  MoveToTopList(i_tree->pointlist);
  xrange[0]=yrange[0]=1.0e10;
  xrange[1]=yrange[1]=-1.0e10;
  do{
	  if(i_tree->pointlist->current->invmag < 0){
		  source[0]=i_tree->pointlist->current->image->x[0];
		  source[1]=i_tree->pointlist->current->image->x[1];

		  if(source[0]<xrange[0]) xrange[0]=source[0];
		  if(source[0]>xrange[1]) xrange[1]=source[0];
		  if(source[1]<yrange[0]) yrange[0]=source[1];
		  if(source[1]>yrange[1]) yrange[1]=source[1];

		  //printf("%e %e\n",source[0],source[1]);
		  break;
	  }
  }while(MoveDownList(i_tree->pointlist));
  //source[0]=xrange[0] + ran2(&seed)*(xrange[1]-xrange[0]);
  //source[1]=yrange[0] + ran2(&seed)*(yrange[1]-yrange[0]);
*/

  printf("%.2f  %.3f  %.3f  %.3f   %i  %e\n",lens->host_sigma,lens->host_axis_ratio,lens->zlens,lens->zsource
		  ,lens->sub_N,lens->source_r);

  // randomize lens
    //RandomizeHost(lens,r_source_phys,&seed,True);
  //RandomizeSubstructure2(lens,2,&seed);
   //lens->NSubstruct=10;
   //PrintAnaLens(lens,False,False);

  // free old tree to speed up image finding
  emptyTree(i_tree);
  emptyTree(s_tree);
  // build new inital grid
  i_points=NewPointArray(Ngrid*Ngrid,True);
  xygridpoints(i_points,range,center,Ngrid,0);
  s_points=LinkToSourcePoints(i_points,Ngrid*Ngrid);
  rayshooterInternal(Ngrid*Ngrid,i_points,i_tree,True);
  // fill trees
  FillTree(i_tree,i_points,Ngrid*Ngrid);
  FillTree(s_tree,s_points,Ngrid*Ngrid);

  //strcpy(filename,lens->outputfile);
  //sprintf(string,"%ld",k);
  //strcat(filename,string);
  file=fopen(lens->outputfile,"w");


  printf("\n %i\n",Nsources);
  for(i=0; i < Nsources ;++i){
		  //printf("source %i\n",i);

	  j=0;
	  do{
		  lens->source_x[0]=xrange[0] + ran2(&seed)*(xrange[1]-xrange[0]);
		  lens->source_x[1]=yrange[0] + ran2(&seed)*(yrange[1]-yrange[0]);
		  //printf("source = %e %e\n",source[0],source[1]);

		  printf("finding images\n");

		  // find source centers
		  lens->source_r=1.0e-6;
		  find_images(lens->source_x,lens->source_r,s_tree,i_tree,&Nimages
				  ,imageinfo,NimageMax,&NImagePoints,0,False,2,verbose,True);
		  ++j;
	  }while(Nimages < 4);
	  printf("%li source position discarded before more than 4 images found\n",j);
	  printf("\n %li  %e %e\n",i,lens->source_x[0],lens->source_x[1]);
	  fprintf(file,"%li  %e %e\n",i,lens->source_x[0],lens->source_x[1]);

	  centers=NewPointArray(Nimages,True);

	  lens->source_r=1.0e-9;
	  find_images(lens->source_x,lens->source_r,s_tree,i_tree,&Nimages
			  ,imageinfo,NimageMax,&NImagePoints,0,False,2,verbose,True);

	  // set star regions
	  for(k=0,lens->star_Nregions=0;k<Nimages;++k){
		  if(imageinfo[k].area/(pi*pow(lens->source_r,2)) > 1.0e-3){
			  for(j=0,rmin=1.0e99;j<imageinfo[k].Npoints;++j){
				  r=sqrt( pow(imageinfo[k].points[j].image->x[0] - lens->source_x[0],2)
					  + pow(imageinfo[k].points[j].image->x[1] - lens->source_x[1],2) );
				  if(rmin>r){
					  rmin=r;
					  PointCopyData(&centers[k],&(imageinfo[k].points[j]));
				  }
			  }
			  ++(lens->star_Nregions);
		  }
	  }

	  //exit(0);
	  //write out extra parameters here and in printAnaLens, fraction in stars etc
	  // overlapping regions

	  for(m=0;m<3;++m){

		  if(m==0){
			  // re-randomize star positions so that host can be kept
			  //  the same while the stars are changed
			  //time(&seed);

			  // populate lens with stars
			  rayshooterInternal(lens->star_Nregions,centers,i_tree,False);
			  implant_stars(lens,centers,lens->star_Nregions,&seed);
			  //PrintAnaLens(lens,False,False);
			  //printf("star 1000 x = %e %e\n",lens->stars_xp[1000][0],lens->stars_xp[1000][1]);
			  //exit(0);
		  }
		  if(m==1){
			  lens->stars_N=0;
		  }
		  if(m==2) lens->sub_N=0;

			  //////////////////////////////
			  // redo grid with stars in it
			  // free old tree to speed up image finding
		  emptyTree(i_tree);
		  emptyTree(s_tree);

			  // build new initial grid
		  i_points=NewPointArray(Ngrid*Ngrid,True);
		  xygridpoints(i_points,range,center,Ngrid,0);
		  s_points=LinkToSourcePoints(i_points,Ngrid*Ngrid);
		  rayshooterInternal(Ngrid*Ngrid,i_points,i_tree,True);
			  // fill trees
		  FillTree(i_tree,i_points,Ngrid*Ngrid);
		  FillTree(s_tree,s_points,Ngrid*Ngrid);


		  fprintf(file,"%i %i\n",m,8*15 + 1);
		  // calculate the magnifications starting with
		  for(r_source_phys = 1.0e-2;r_source_phys >= 1.0e-7*0.99999e-3
		        ;r_source_phys /= pow(10,1.0/15.0) ){
 
 			  lens->source_r = r_source_phys*angDist(0,lens->zlens)/angDist(0,lens->zsource);

			  if(lens->source_r > 1.0e-6)
				  find_images(lens->source_x,lens->source_r,s_tree,i_tree,&Nimages
					  ,imageinfo,NimageMax,&NImagePoints,0,True,2,verbose,True);
			  else
				  find_images(lens->source_x,lens->source_r,s_tree,i_tree,&Nimages
					  ,imageinfo,NimageMax,&NImagePoints,0,False,2,verbose,True);

			  combineCloseImages(min_image_seporation*angDist(0,lens->zlens)*pi/60/60/180
					  ,imageinfo,&Nimages,&NewNimages);

			  printf("%i  %i  %e\n",Nimages,NewNimages,r_source_phys);
			  fprintf(file,"%i  %i  %e\n",Nimages,NewNimages,r_source_phys);
			  for(k=0;k<NewNimages;++k){

				  // find the magnification of the point closest to the source
				  for(j=0,rmin=1.0e99;j<imageinfo[k].Npoints;++j){
					  r=sqrt( pow(imageinfo[k].points[j].image->x[0] - lens->source_x[0],2)
						  + pow(imageinfo[k].points[j].image->x[1] - lens->source_x[1],2) );
					  if(rmin>r){
						  rmin=r;
						  point=&(imageinfo[k].points[j]);
					  }
				  }

				  // find the magnification at the image centroid
				  tmp_point->x[0]=imageinfo[k].centroid[0];
				  tmp_point->x[1]=imageinfo[k].centroid[1];

				  rayshooterInternal(1,tmp_point,i_tree,False);
				  rayshooterInternal(1,point,i_tree,False);

				  printf("%e %e     %e %e %e\n",imageinfo[k].centroid[0],imageinfo[k].centroid[1]
				      ,imageinfo[k].area/(pi*pow(lens->source_r,2)),1.0/tmp_point->invmag
				      ,1.0/point->invmag);
				  fprintf(file,"%e %e     %e %e %e\n",imageinfo[k].centroid[0],imageinfo[k].centroid[1]
				      ,imageinfo[k].area/(pi*pow(lens->source_r,2)),1.0/tmp_point->invmag
				      ,1.0/point->invmag);
				  /*printf("    %li   %e %e   %e %e %e\n",imageinfo[k].Npoints,imageinfo[k].area
				 		  ,imageinfo[k].area_error,imageinfo[k].gridrange[0],imageinfo[k].gridrange[1]
						                                          ,imageinfo[k].gridrange[2]);
				  fprintf(file,"    %li   %e %e   %e %e %e\n",imageinfo[k].Npoints,imageinfo[k].area
						  ,imageinfo[k].area_error,imageinfo[k].gridrange[0],imageinfo[k].gridrange[1]
						                                          ,imageinfo[k].gridrange[2]);
						                                          */
			  }
			  //exit(0);
		  }
		  //PrintList(i_tree->pointlist);
		  //exit(0);

		  /* print all points in images */
		  if(!verbose && !just_mags){
			  for(k=0;k<Nimages;++k){
				  printf(" %li\n",imageinfo[k].Npoints);
				  for(j=0;j<imageinfo[k].Npoints;++j){
					  printf(" %e %e %e\n",imageinfo[k].points[j].x[0]
					                      ,imageinfo[k].points[j].x[1]
					                      ,imageinfo[k].points[j].gridsize);
				  }

				  //PrintList(imageinfo[k].outerborder);
			  }
		  }
		  time(&t4);
	  }
  }

  //free(critical->points);
  //freeImageInfo(critical,Ncrit);
  //FreePointArray(caustic_points);
  fprintf(file,"\n %li  %e %e source\n",i,lens->source_x[0],lens->source_x[1]);
  fprintf(file,"\nlens->host_sigma lens->host_axis_ratio lens->zlens lens->zsource lens->sub_N lens->source_r");
  fprintf(file,"\n%.2f  %.3f  %.3f  %.3f   %i  %e\n",lens->host_sigma,lens->host_axis_ratio,lens->zlens,lens->zsource
 		  ,lens->sub_N,lens->source_r);


  printf("\n\n");
  PrintAnaLens(lens,False,False);

  printf("\nNumber of lenses =  %li average number of sources per lens =  %0.2f\n",m,Ntotal_lenses*1.0/m);
  fprintf(file,"\nNumber of lenses =  %li average number of sources per lens =  %0.2f\n",m,Ntotal_lenses*1.0/m);

  time(&now);
  printf("all images found in %f min\n     %f sec per source    %f sec per lens average\n"
			  ,difftime(now,to)/60.
			  ,difftime(now,to)/Ntotal_lenses
			  ,difftime(now,to)*m/Ntotal_lenses
  );
  fprintf(file,"all images found in %f min\n     %f sec per source    %f sec per lens average\n"
			  ,difftime(now,to)/60.
			  ,difftime(now,to)/Ntotal_lenses
			  ,difftime(now,to)*m/Ntotal_lenses
  );

  printf("\nNumber of rays shot for last source:  %li\n",i_tree->pointlist->Npoints);
  //printf("param file : %s\n",paramfile);

  fclose(file);

  free(imageinfo->points);
  freeImageInfo(imageinfo,NimageMax);
  //FreePointArray(caustic_points);
  FreePointArray(tmp_point->image);
  FreePointArray(tmp_point);
  freeTree(i_tree);
  freeTree(s_tree);

  free(times);
  //fclose(file);

  return 1;
}
