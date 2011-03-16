/* fixedimages.c
 * driver program that is designed to calculate the variations in
 * lensing properties with fixed image positions
 */

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <string.h>
#include <nr.h>
#include "../Library/Recipes/nrutil.h"
#include <nrD.h>
/*#include "../Library/Recipes/nrutil.c"
#include "../Library/Recipes/ran2.c"
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
#include "../Library/Recipes/gasdev.c"*/

#include "../Library/cosmo.h"
/*#include "../Library/powerCDM.c"
#include "../Library/cosmo.c"*/

#include "../SLsim/AnalyticNSIE/analytic_lens.h"
#include "../SLsim/TreeCode_link/Tree.h"
#include "../SLsim/TreeCode/TreeNB.h"
#include "../SLsim/FitLens/fitlens.h"

COSMOLOGY cosmo;
AnaLens *lens=0;

/*
 *   fullrange.c calculates the magnifications of circular images as a function of
 *   source size including substructures and microlensing by stars
 */

int main(int arg,char **argv){
  TreeHndl i_tree,s_tree;
  Point *i_points,*s_points,*caustic_points,*tmp_point,*point=0,*centers;
  double range,r_source_phys,center[2],kappa=0
                          ,xrange[2],yrange[2],rmin,r,*times;
  unsigned long i,j,k,m,Ngrid,NImagePoints,Ncritpoints;
  int Nimages,NewNimages;
  long seed=282923;
  //long seed=28374;
  //time(&seed);
  const int NimageMax = 750;


  printf("seed = %li\n",seed);

  double min_image_seporation,**subdeflect,alpha[2];
  ImageInfo *imageinfo,*critical;
  //FILE *file;
  int Ncrit,Nsources,Ntotal_lenses,refresh;
  time_t to,now,t4;
  Boolean verbose,nocrit,just_mags,success;
  char *paramfile;
  FILE *file;

  Boolean FiniteSource = False;

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
  Nsources = 1000;
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

   // read parameter file in
  if(paramfile==NULL){
	  paramfile=(char *)malloc(60*sizeof(char));
	  sprintf(paramfile,"FullRange/paramfile");
  }
  lens=(AnaLens *)calloc(1,sizeof(AnaLens));
  readparams_ana(paramfile,&cosmo,lens);

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
   //lens->NSubstruct=10;
   //PrintAnaLens(lens,False,False);


  //strcpy(filename,lens->outputfile);
  //sprintf(string,"%ld",k);
  //strcat(filename,string);

  file=fopen(lens->outputfile,"w");

  j=0;
  do{
	  lens->source_x[0] = xrange[0] + ran2(&seed)*(xrange[1]-xrange[0]);
	  lens->source_x[1] = yrange[0] + ran2(&seed)*(yrange[1]-yrange[0]);
	  //printf("source = %e %e\n",source[0],source[1]);

	  printf("finding images\n");

	  // find source centers
	  lens->source_r = 1.0e-6;
	  find_images(lens->source_x,lens->source_r,s_tree,i_tree,&Nimages
			  ,imageinfo,NimageMax,&NImagePoints,0,False,2,verbose,True);

	  ++j;
  }while(Nimages != 4);

  subdeflect = dmatrix(0,Nimages-1,0,1);
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
				  centers[k].x[0] = imageinfo[k].points[j].x[0];
				  centers[k].x[1] = imageinfo[k].points[j].x[1];

				  //PointCopyData(&centers[k],&(imageinfo[k].points[j]));
			  }
		  }
		  ++(lens->star_Nregions);
	  }
  }

   // image positions for 1422
    centers[0].x[0] = (0.385 - 0.742)/lens->MpcToAsec;
    centers[0].x[1] = (0.317 + 0.656)/lens->MpcToAsec;

    centers[1].x[0] = (0 - 0.742)/lens->MpcToAsec;
    centers[1].x[1] = (0 + 0.656)/lens->MpcToAsec;

    centers[2].x[0] = (-0.336 - 0.742)/lens->MpcToAsec;
    centers[2].x[1] = (-0.750 + 0.656)/lens->MpcToAsec;

    centers[3].x[0] = ( 0.948 - 0.742)/lens->MpcToAsec;
    centers[3].x[1] = (-0.802 + 0.656)/lens->MpcToAsec;

  time(&to);

  kappa = lens->sub_Ndensity;
  for(lens->sub_Mmin = 1.0e5;lens->sub_Mmin <= 1.001e9;lens->sub_Mmin *= pow(10,4.0/10.0)){

	  // keep substructure surface mass density constant
	  lens->sub_Ndensity = kappa * lens->Sigma_crit /averageSubMass(lens);

  printf("\n %i\n",Nsources);
  for(i=0; i < Nsources ;++i){

	  // randomize substructure
	  if(i>0){

		  RandomizeSubstructure3(lens,2,&seed);

		  // find deflections caused by substructure at images
		  for(k=0;k < 4;++k){
			  subdeflect[k][0] = 0;
			  subdeflect[k][1] = 0;
			  for(j=0;j<lens->sub_N;++j){
				  lens->sub_alpha_func(alpha,centers[k].x,lens->sub_Rcut[j]
				                 ,lens->sub_mass[j],lens->sub_beta
			                     ,lens->sub_x[j],lens->Sigma_crit);
				  subdeflect[k][0] -= alpha[0];  // sign checked
				  subdeflect[k][1] -= alpha[1];
			  }
			  //printf("subdeflect[%i] = %e %e\n",k,subdeflect[k][0],subdeflect[k][1]);
		  }

		  FindLensSimple(lens,4,centers,lens->source_x,subdeflect);
	  }

	  //printf("%li source position discarded before more than 4 images found\n",j);
	  printf("\n %li  %e %e\n",i,lens->source_x[0],lens->source_x[1]);
	  fprintf(file,"%li  %e %e  ",i,lens->source_x[0],lens->source_x[1]);

	  //exit(0);
	  //write out extra parameters here and in printAnaLens, fraction in stars etc
	  // overlapping regions

	  if(FiniteSource){
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
	  }
		  //fprintf(file,"%i %i\n",m,8*40 + 1);
		  // calculate the magnifications starting with
		  // for(r_source_phys = 1.0e-2, refresh = 0;r_source_phys >= 1.0e-7*0.99999e-3
		 //       ;r_source_phys /= pow(10,1.0/40.0), ++refresh ){
	  for(r_source_phys = 1.0e-7, refresh = 0;r_source_phys >= 0.99999e-7
			        ;r_source_phys /= pow(10,1.0/40.0), ++refresh ){

			  if(FiniteSource){
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
			  }

			  printf("%e  %e  %e  %e  %e\n",lens->sub_alpha,lens->sub_Ndensity
					  ,lens->sub_Ndensity*averageSubMass(lens)/lens->Sigma_crit
					  ,lens->sub_Mmax,lens->sub_Mmin);
			  fprintf(file,"%e  %e  %e  %e  ",lens->sub_alpha,lens->sub_Ndensity
					  *averageSubMass(lens)/lens->Sigma_crit
					  ,lens->sub_Mmax,lens->sub_Mmin);

			  if(!FiniteSource) NewNimages = Nimages;

			  for(k=0;k<NewNimages;++k){

				  if(FiniteSource){
					  // find closest target image and print displacement
					  for(j=0,r=1.0e100;j<4;++j){
						  r = MIN(r,sqrt(pow(imageinfo[k].centroid[0] - centers[j].x[0],2)
						       	  + pow(imageinfo[k].centroid[1] - centers[j].x[1],2)
							  ) );

					  }
					  printf("%e  ",r);

					  // find the magnification of the point closest to the source
					  for(j=0,rmin=1.0e99;j<imageinfo[k].Npoints;++j){
						  r=sqrt( pow(imageinfo[k].points[j].image->x[0] - lens->source_x[0],2)
						  + pow(imageinfo[k].points[j].image->x[1] - lens->source_x[1],2) );
						  if(rmin>r){
							  rmin=r;
							  point=&(imageinfo[k].points[j]);
						  }
					  }
				  }


				  // find the magnification at the image centers
				  tmp_point->x[0] = centers[k].x[0];
				  tmp_point->x[1] = centers[k].x[1];

				  rayshooterInternal(1,tmp_point,i_tree,False);

				  if(FiniteSource){
					  rayshooterInternal(1,point,i_tree,False);

					  printf("%e %e     %e %e %f\n",imageinfo[k].centroid[0],imageinfo[k].centroid[1]
                                  ,imageinfo[k].area/(pi*pow(lens->source_r,2)),1.0/tmp_point->invmag
                                  ,1.0/point->invmag);
					  fprintf(file,"%e %e     %e %e %e\n",imageinfo[k].centroid[0],imageinfo[k].centroid[1]
					               ,imageinfo[k].area/(pi*pow(lens->source_r,2)),1.0/tmp_point->invmag
				                   ,1.0/point->invmag);
				  }else{
					  printf("%e %e     %e  %f\n",centers[k].x[0]*lens->MpcToAsec,centers[k].x[1]*lens->MpcToAsec
							  ,1.0/tmp_point->invmag,tmp_point->dt);
					  fprintf(file,"%e %e  %e %e    ",centers[k].x[0]*lens->MpcToAsec
							  ,centers[k].x[1]*lens->MpcToAsec,1.0/tmp_point->invmag,tmp_point->dt);
				  }
			  }
			  //exit(0);
		  }
		  fprintf(file,"\n");
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
