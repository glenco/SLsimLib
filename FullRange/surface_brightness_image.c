
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <omp.h>
#include <nr.h>
#include "../Library/Recipes/nrutil.h"
#include "../Library/Recipes/nrutil.c"
#include "../Library/Recipes/ran2.c"
#include <nrD.h>
#include "../Library/RecipesD/locateD.c"
#include "../Library/RecipesD/powellD.c"
#include "../Library/RecipesD/odeintD.c"
#include "../Library/RecipesD/bsstepD.c"
#include "../Library/RecipesD/mmidD.c"
#include "../Library/RecipesD/pzextrD.c"
#include "../Library/RecipesD/polintD.c"
#include "../Library/Recipes/poidev.c"
#include "../Library/Recipes/gammln.c"
#include "../Library/RecipesD/dfridrD.c"
#include "../Library/Recipes/gasdev.c"

#include "../Library/cosmo.h"
#include "../Library/powerCDM.c"
#include "../Library/cosmo.c"

#include "../SLsim/AnalyticNSIE/analytic_lens.h"
#include "../SLsim/TreeCode_link/Tree.h"
//#include "../TreeCode/TreeNB.h"
#include "../SLsim/TreeCode_link/List/List.h"
#include "../SLsim/Kist/Kist.h"
#include "../SLsim/TreeCode_link/map_images.h"
#include "../SLsim/ImageProcessing/image_processing.h"

COSMOLOGY cosmo;
AnaLens *lens=0;

int main(int arg,char **argv){
  TreeHndl i_tree,s_tree;
  Point *i_points,*s_points,*caustic_points,*tmp_point,*centers;
  double range,r_source_phys,center[2]
                          ,xrange[2],yrange[2],rmin,r,*times;
  unsigned long i,j,m,Ngrid,NImagePoints,Ncritpoints,Npixels=512*2;
  int Nimages;
  long seed=282923;
  double *map,map_range,map_center[2],tmp,area=0;
  //long seed=28374;
  //time(&seed);
  long k;

  printf("seed = %li\n",seed);

  double min_image_seporation;
  ImageInfo *imageinfo,*critical;
  FILE *file;
  int Ncrit,Nsources,Ntotal_lenses;
  time_t to,now,t4;
  Boolean verbose,nocrit,just_mags,success;
  char *paramfile,*string=(char *)malloc(12*sizeof(char)),*filename=(char *)malloc(50*sizeof(char));

  const int Nimagemax=1000;

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

  imageinfo=NewImageInfo(Nimagemax);

  assert(imageinfo->imagekist->Nunits == 0);
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

  assert(imageinfo->imagekist->Nunits == 0);

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
				  ,imageinfo,&NImagePoints,0,False,2,verbose,True);
		  ++j;
	  }while(Nimages < 4);
	  printf("%li source position discarded before more than 4 images found\n",j);
	  printf("\n %li  %e %e\n",i,lens->source_x[0],lens->source_x[1]);

	  centers=NewPointArray(Nimages,True);

	  lens->source_r=1.0e-9;
	  find_images(lens->source_x,lens->source_r,s_tree,i_tree,&Nimages
			  ,imageinfo,&NImagePoints,0,False,2,verbose,True);

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

	  //write out extra parameters here and in printAnaLens, fraction in stars etc
	  // overlapping regions

	  m = 0;
	  switch(m)
	  {
	  case 0:
		  // populate lens with stars
		  //*** why do we need this
		  rayshooterInternal(lens->star_Nregions,centers,i_tree,False);
		  implant_stars(lens,centers,lens->star_Nregions,&seed);
		  break;
	  case 1:
		  lens->stars_N = 0;
		  break;
	  case 2:
		  lens->sub_N = 0;
		  break;
	  }


	  map = (double *)malloc(Npixels*Npixels*sizeof(double));
	  map_range = 3.5e-7;

	  r_source_phys = 1.0e-7;
	  lens->source_r = r_source_phys*angDist(0,lens->zlens)/angDist(0,lens->zsource);
	  lens->source_r2 = pow(0.3*lens->source_r,2);

	  lens->source_r=0.423e-6;  // r_max

	  double r_in = 2.3884456e-9,r_out = 4.2992021e-7; // Mpc
	  int NmacroImages = Nimages,n = 0;
	  double **macro_images = 0;

	  macro_images = dmatrix(0,NmacroImages-1,0,1);


	  for(n=0;n<NmacroImages;++n){
		  macro_images[n][0] = imageinfo[n].centroid[0];
		  macro_images[n][1] = imageinfo[n].centroid[1];
	  }
	  map_center[0] = imageinfo[2].centroid[0];
	  map_center[1] = imageinfo[2].centroid[1];
	  // find image that is closest to map center
	  rmin = min_image_seporation*angDist(0,lens->zlens)*pi/60/60/180;

	  printf("\n\n ***** start *****\n\n");
/*
	  ListHndl neighborlist = NewList();

	  time(&to);
	  for(i=0;i<imageinfo[3].Npoints;++i) PointsWithin(s_tree,imageinfo[3].points[i].image->x,lens->source_r,neighborlist,0);
	  time(&now);
	  printf("%e\n",difftime(now,to));
	  printf("N = %li\n",neighborlist->Npoints);

	  time(&to);
	  for(i=0;i<imageinfo[3].Npoints;++i) PointsWithin_iter(s_tree,imageinfo[3].points[i].image->x,lens->source_r,neighborlist,0);
	  time(&now);
	  printf("%e\n",difftime(now,to));
	  printf("N = %li\n",neighborlist->Npoints);
	  //PrintList(neighborlist);

	  exit(0);
*/

	  //for(lens->source_tau = 100 ; lens->source_tau > 0.009 ; lens->source_tau /= pow(100/0.02,1.0/20) ){
	  //for(k = 202; k > -1 ; --k ){
	  for(k = 64; k < 203 ; ++k ){

		  lens->source_tau = 0.1 + k*15./200.;

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

		  // refine inner ring
		  lens->source_r = sqrt(2*r_out*lens->source_tau* 8.39428142e-10 -
				  pow(lens->source_tau* 8.39428142e-10,2));  // speed of light in Mpc/day

		  map_images(lens,s_tree,i_tree,&Nimages,imageinfo,Nimagemax,range/Ngrid,True,FillHoles,True);
		  pixelize(map,Npixels,map_range,map_center,imageinfo->imagekist,False,True);

		  // open file
		  strcpy(filename,lens->outputfile);
		  sprintf(string,"%ld",k);
		  strcat(filename,string);
		  file=fopen(filename,"w");

		  // print magnification of macro images
		  printf("%i\n",NmacroImages);
		  for(n=0;n<NmacroImages;++n){
			  for(j=0,area=0.0;j<Nimages;++j){
				  tmp = pow(imageinfo[j].centroid[0]-macro_images[n][0],2)
					  + pow(imageinfo[j].centroid[1]-macro_images[n][1],2);
				  if(tmp < rmin*rmin){ area += imageinfo[j].area; }
			  }
			  fprintf(file,"%e  ",area);
			  printf("%e  ",area);
		  }
		  fprintf(file,"\n");
		  printf("\n");

		  // print all points in images
		  //printf("%i %e\n",Nimages,r_source_phys);
		  fprintf(file,"1 %e %e\n",lens->source_r,lens->source_tau);
		  fprintf(file,"%e %e  %e  %e %li\n",map_center[0],map_center[1],area,map_range,Npixels);
		  for(j=0 ; j < Npixels*Npixels ; ++j ) fprintf(file,"%e\n",map[j]);
		  fclose(file);
	  }

	  time(&t4);
  }


  //free(critical->points);
  //freeImageInfo(critical,Ncrit);
  //FreePointArray(caustic_points);

  printf("\n\n");
  PrintAnaLens(lens,False,False);

  printf("\nNumber of lenses =  %li average number of sources per lens =  %0.2f\n",m,Ntotal_lenses*1.0/m);

  time(&now);
  printf("all images found in %f min\n     %f sec per source    %f sec per lens average\n"
			  ,difftime(now,to)/60.
			  ,difftime(now,to)/Ntotal_lenses
			  ,difftime(now,to)*m/Ntotal_lenses
  );

  printf("\nNumber of rays shot for last source:  %li\n",i_tree->pointlist->Npoints);
  //printf("param file : %s\n",paramfile);

  free(imageinfo->points);
  freeImageInfo(imageinfo,Nimagemax);
  //FreePointArray(caustic_points);
  FreePointArray(tmp_point->image);
  FreePointArray(tmp_point);
  freeTree(i_tree);
  freeTree(s_tree);

  free(times);
  //fclose(file);

  return 1;
}
