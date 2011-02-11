
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

AnaLens *lens = 0;

int main(int arg,char **argv){
  TreeHndl i_tree,s_tree;
  Point *i_points,*s_points,*caustic_points,*tmp_point,*point;
  double range,r_source_phys,center[2],tmp,tmp2
                          ,xrange[2],yrange[2],rmin,r,*times,xtmp[2];
  unsigned long i,j,k,m,Ngrid,NImagePoints,Ncritpoints;
  int Nimages,NewNimages;
  long seed=282923,i_total=0;
  //long seed=282925,i_total=0;
  //long seed=28374;

  printf("seed = %li\n",seed);

  double min_image_seporation,invmag_bg;
  ImageInfo *imageinfo,*critical;
  //FILE *file;
  int Ncrit,tange_caustic=0,Nsources,Ntotal_lenses;
  time_t to,now,t3,t4;
  Boolean verbose,nocrit,just_mags,success;
  const int NimageMax = 100;
  char *paramfile;

  verbose=False;
  //verbose=True;
  nocrit=False;

  just_mags=True;

  if(arg==2) paramfile=argv[1];
  else paramfile=NULL;

  time(&to);
  range=0.5;
  Ngrid=64;
  Ntotal_lenses=100000;
  Nsources=2000;
//  Ntotal_lenses=1000;

  //Ntotal_lenses=5;
  //Nsources=50;

  // select source
  r_source_phys=1.0e-5;
  //r_source_phys=1.0e-6;
  min_image_seporation=0.1;  // in arcsec

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
 
  // initialize tree and lens model
  i_points=NewPointArray(Ngrid*Ngrid,True);
  xygridpoints(i_points,range,center,Ngrid,0);
  s_points=LinkToSourcePoints(i_points,Ngrid*Ngrid);
  i_tree=NULL;
  rayshooterInternal(Ngrid*Ngrid,i_points,i_tree,False);

  // build trees
  i_tree=BuildTree(i_points,Ngrid*Ngrid);
  s_tree=BuildTree(s_points,Ngrid*Ngrid);

  //file=fopen(outputfile,"w");

   // loop through lenses
  i_total=-1;
  m=0;
  while(i_total < Ntotal_lenses){
	  ++m;
	  //for(m=0,i_total=0;m<Ntotal_lenses/Nsources;++m){

	  do{
		  // randomize lens
		  RandomizeHost(lens,r_source_phys,&seed,True);
		  RandomizeSubstructure2(lens,2,&seed);

		  // free old tree to speed up image finding
		  emptyTree(i_tree);
		  emptyTree(s_tree);
		  // build new inital grid
		  i_points=NewPointArray(Ngrid*Ngrid,True);
		  xygridpoints(i_points,range,center,Ngrid,0);
		  s_points=LinkToSourcePoints(i_points,Ngrid*Ngrid);
		  rayshooterInternal(Ngrid*Ngrid,i_points,i_tree,False);
		  // fill trees
		  FillTree(i_tree,i_points,Ngrid*Ngrid);
		  FillTree(s_tree,s_points,Ngrid*Ngrid);

		  // find critical curves

		  critical=find_crit(s_tree,i_tree,&Ncrit,5.0e-5,&success,False,verbose);

		  if(Ncrit == 0){
			  //printf("error: Ncrit=0 \n");
			  //PrintAnaLens(lens,False,False);
			  //printf("\n %li\n",m);
		  }
		  printf("\n\n");

		  /* select critical curve with largest area
		   *     this will generally be the tangential critical curve if
		   */

		  xtmp[0]=0.0;
		  xtmp[1]=0.0;
		  //printf(" %i\n",Ncrit);
		  for(i=0,tmp=0,tange_caustic=0;i<Ncrit;++i){
			  if( tmp < critical[i].area &&
				  0 < abs(windings(xtmp,critical[i].points,critical[i].Npoints,&tmp2,0))
			  ){
				  tmp = critical[i].area;
				  tange_caustic=i;
			  }
		  }

		  if(tmp <= 0){
			  /*
			  printf("error: no tangential caustic found\n");
			  for(i=0;i<Ncrit;++i){
				  printf("%li area=%e windings=%i\n",i,critical[i].area,
					  windings(xtmp,critical[i].points,critical[i].Npoints,&tmp2,0) );
			  }
			  for(i=0;i<Ncrit;++i){
				  printf("%li\n",critical[i].Npoints);
				  for(j=0;j<critical[i].Npoints;++j) printf("%e %e\n",critical[i].points[j].x[0]
				                                                  ,critical[i].points[j].x[1] );
			  }

			  if(Ncrit != 0){
				  PrintAnaLens(lens,False,False);
				  printf("\n %li\n",m);
			  }
			  */
			  //exit(0);
		  }

	  }while(Ncrit == 0 || tmp <= 0.0);

	  if(i_total==-1) PrintAnaLens(lens,False,False);

   // re-assign caustic curve
	  Ncritpoints=critical[tange_caustic].Npoints;
	  caustic_points=NewPointArray(Ncritpoints,True);
	  for(j=0;j<Ncritpoints;++j){
		  caustic_points[j].x[0] = critical[tange_caustic].points[j].image->x[0];
		  caustic_points[j].x[1] = critical[tange_caustic].points[j].image->x[1];
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


  // calculate the area within the caustic
	  windings(caustic_points->x,caustic_points,Ncritpoints,&tmp,0);
	  //fprintf(file,"%e\n",tmp);
	  // area is made negative when find_caustic has failed to order critcurve properly

	  invmag_bg=fabs(pow(1-lens->perturb_modes[0],2)-pow(lens->perturb_modes[1],2)-pow(lens->perturb_modes[2],2));
	  invmag_bg=1.0;

	  if(success) printf("%e %e\n",tmp,pi*pow(lens->host_ro,2)/invmag_bg);
	  else printf("%e %e\n",-tmp,pi*pow(lens->host_ro,2)/invmag_bg);

	  printf("%.2f  %.3f  %.3f  %.3f   %i  %e\n",lens->host_sigma,lens->host_axis_ratio,lens->zlens,lens->zsource
			  ,lens->sub_N,lens->source_r);

	  //fprintf(file,"\n%i\n",Nsources);

	  printf("\n %i\n",(int)(Nsources*tmp/pi/pow(lens->host_ro,2)/invmag_bg + 0.5));
	  for(i=0; (i < (int)(Nsources*tmp/pi/pow(lens->host_ro,2)/invmag_bg + 0.5)) ;++i){
		  //printf("source %i\n",i);

		  do{
			  lens->source_x[0]=xrange[0] + ran2(&seed)*(xrange[1]-xrange[0]);
			  lens->source_x[1]=yrange[0] + ran2(&seed)*(yrange[1]-yrange[0]);
			  //printf("source = %e %e\n",source[0],source[1]);
		  }while(abs(windings(lens->source_x,caustic_points,Ncritpoints,&tmp,0)) < 1);

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

		  time(&t3);

		  ++i_total;
		  //fprintf(file,"\n%i  %e %e  %e\n",i,lens->source_x[0],lens->source_x[1],r_source_phys);
		  printf("\n %li %li  %e %e  %e\n",i_total,i,lens->source_x[0],lens->source_x[1],r_source_phys);

		  // find images
		  if(verbose) printf("finding images\n");

		  //find_images(lens->source_x,lens->source_r,s_tree,i_tree,&Nimages
		 //		  ,imageinfo,&NImagePoints,0,True,2,verbose,True);
		  find_images(lens->source_x,lens->source_r,s_tree,i_tree,&Nimages
			  ,imageinfo,NimageMax,&NImagePoints,0,False,2,verbose,True);

		  // find the point source magnification for each image
		  if(verbose) printf(" %i images\n",Nimages);
		  //else //fprintf(file,"%i\n",Nimages);
			  //printf(" %li\n",Nimages);


		  /*
		  for(k=0;k<Nimages;++k){
			  //printf("%i\n",imageinfo[k].Npoints);


			  tmp_point->x[0]=imageinfo[k].centroid[0];
			  tmp_point->x[1]=imageinfo[k].centroid[1];

			  if(tmp==0){
				  PrintPoint(tmp_point);
				  //fprintf(file,"tmp=%e  Npoints in image = %i\n",tmp,imageinfo[k].Npoints);
				  printf("tmp=%e  Npoints in image =  %li\n",tmp,imageinfo[k].Npoints);
			  }else{
			  rayshooterInternal(1,tmp_point,i_tree,False);

				  //fprintf(file,"%e %e  %e   %e %e     %e\n",x[0]/tmp,x[1]/tmp
				  //           ,imageinfo[k].area/(pi*pow(lens->r_source))
				  //,1.0/point->invmag,1.0/tmp_point->invmag
				  //,1-imageinfo[k].area*fabs(point->invmag)/(pi*r_source*r_source));


			  printf("%e %e  %e   %e %e     %e\n"
					  ,imageinfo[k].centroid[0],imageinfo[k].centroid[1]
						  ,imageinfo[k].area/(pi*pow(lens->r_source,2))
						  ,1.0/point->invmag,1.0/tmp_point->invmag
						  ,1-imageinfo[k].area*fabs(tmp_point->invmag)/(pi*pow(lens->r_source,2)));


			  }
		  }
		  */

			  // remove images that are too close together
		  combineCloseImages(min_image_seporation*angDist(0,lens->zlens)*pi/60/60/180
		          ,imageinfo,&Nimages,&NewNimages);
		  //NewNimages=Nimages;
		  printf("%i  %i\n",Nimages,NewNimages);
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
			  //printf("%e %e   %e %e\n",point->x[0],point->x[1],point->image->x[0]
			  //   ,point->image->x[1]);

			  // find the magnification at the image centroid
			  tmp_point->x[0]=imageinfo[k].centroid[0];
			  tmp_point->x[1]=imageinfo[k].centroid[1];
			  /*if(tmp==0){
				  PrintPoint(tmp_point);
				  //fprintf(file,"tmp=%e  Npoints in image = %i\n",tmp,imageinfo[k].Npoints);
				  printf("tmp=%e  Npoints in image =  %li\n",tmp,imageinfo[k].Npoints);
			  }else{*/
			  rayshooterInternal(1,tmp_point,i_tree,False);
			  rayshooterInternal(1,point,i_tree,False);

			  printf("%e %e     %e %e %e\n",imageinfo[k].centroid[0],imageinfo[k].centroid[1]
						  ,imageinfo[k].area/(pi*pow(lens->source_r,2)),1.0/tmp_point->invmag,1.0/point->invmag);
			  //,(imageinfo[k].area/(pi*pow(lens->r_source,2))-1.0/fabs(point->invmag))*fabs(point->invmag));
		  }

		  // test lines
		  /*if(Nimages > 4){
			  //PrintImages(imageinfo,Nimages);
			  exit(0);
		  }*/

		  //PrintList(i_tree->pointlist);
		  //exit(0);

		  /* print all points in images */
		  if(!verbose && !just_mags){
			  for(k=0;k<Nimages;++k){
				  printf(" %li\n",imageinfo[k].Npoints);
				  for(j=0;j<imageinfo[k].Npoints;++j){
					  printf(" %e %e %e\n",imageinfo[k].points[j].x[0]
					                      ,imageinfo[k].points[j].x[1],imageinfo[k].points[j].gridsize);
				  }

				  //PrintList(imageinfo[k].outerborder);
			  }
		  }
		  time(&t4);
		  times[i]=difftime(t4,t3);
	  }
	  free(critical->points);
	  freeImageInfo(critical,Ncrit);
	  FreePointArray(caustic_points);
  }
  //printf(file,"\n%i\n",NImagePoints);

//  for(i=0;i<Nimages;++i){
//	printf("%i\n",imageinfo[i].Npoints);
//    for(j=0;j<imageinfo[i].Npoints;++j){
//      printf("%e  %e  %e\n",imageinfo[i].points[j].x[0],imageinfo[i].points[j].x[1],imageinfo[i].points[j].invmag);
//    }
//  }

  printf("\nNumber of lenses =  %li average number of sources per lens =  %0.2f\n",m,Ntotal_lenses*1.0/m);

  time(&now);
  printf("all images found in %f min\n     %f sec per source    %f sec per lens average\n"
	 ,difftime(now,to)/60.
	 ,difftime(now,to)/Ntotal_lenses
	 ,difftime(now,to)*m/Ntotal_lenses
  );
  /*
  printf("times for each source\n");
  for(i=1;i<Nsources;++i){
	  //printf("    %f sec\n",times[i]);
	  times[0]+=times[i]/Nsources;
  }
  printf("average time for one source %f sec",times[0]);
	*/

  printf("\nNumber of rays shot for last source:  %li\n",i_tree->pointlist->Npoints);
  //printf("param file : %s\n",paramfile);

  free(imageinfo->points);
  freeImageInfo(imageinfo,200);
  //FreePointArray(caustic_points);
  FreePointArray(tmp_point->image);
  FreePointArray(tmp_point);
  freeTree(i_tree);
  freeTree(s_tree);
  free(times);
  //fclose(file);

  return 1;
}
