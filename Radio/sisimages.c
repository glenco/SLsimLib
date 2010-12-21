/* sismages.x datafile outfile */
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "../../Library/Include/nrutil.h"
#include "../../Library/Recipes/nrutil.c"
/*#include "../Include/nr.h"*/
#include "../../Library/locateD.c"
#include "../../Library/odeintD.c"
#include "../../Library/bsstepD.c"
#include "../../Library/polintD.c"
#include "../../Library/pzextrD.c"
#include "../../Library/mmidD.c"
#include "../../Library/powellD.c"

#include "../../Library/cosmo.h"
#include "../../Library/cosmo.c"

#include "../TreeCode_link/utill.c"
#include "../TreeCode_link/Tree.h"
#include "../TreeCode_link/double_sort.c"
#include "../TreeCode_link/Tree.c"
#include "../TreeCode_link/TreeDriver.c"
#include "../TreeCode_link/image_finder.c"

#include "../../NSIEmodel/nsie.c"/**/

struct cosmology cosmo;

static double sigma_lens,axis_ratio,orientation,core_size,zlens; /* lens parameters */
static double zsource,angle,amax,amin;  /* source parameters */

int main(int arg,char **argv){
  TreeHndl i_tree,s_tree;
  ListHndl spointlist,ipointlist;
  Point *i_points,*s_points,*point;
  double range,map_range,flux,source[2],r_source,center[2],linkinglength,area,Dlens,imax,smax,resolution,zold=-1.0;
  unsigned long imagemarkers[500],i,j,k,m,Ngrid,Nimages,NImagePoints,Nmap;
  ImageInfo *imageinfo;
  short moved;
  FILE *file;
  long Nframes=0,ix,iy;
  int idnumber,refine,refinement_level;
  float *map;
  double  i_151, i_610, i_1400, i_4860, i_18000;

  range=4.0;  /* in arcmin */
  Ngrid=50;
  center[0]=0.0; center[1]=0.0;
  refinement_level=3;
  map_range=3.5;   /* in arcmin */

  /* lens */
  zlens=0.2;
  sigma_lens=1.5e3;
  axis_ratio=0.8;
  orientation=pi/4;  /* in radians */
  core_size=0.01;

  /* cosmology */
  cosmo.h=0.7;
  cosmo.Omo=0.3;
  cosmo.Oml=1.0-cosmo.Omo;
  cosmo.w=-1.0;
  cosmo.w1=0.0;

  printf("resolution %f arcsec\n",60*range/(Ngrid-1)/pow(3,refinement_level));

  Dlens=angDist(0,zlens);
  range*=Dlens*pi/60/180;
  map_range*=Dlens*pi/60/180;

  file=fopen(argv[1],"r");

  resolution=range/(Ngrid-1)/pow(3.0,refinement_level+1);
  Nmap=(Ngrid-1)*pow(3,refinement_level+1);
  map=(float *) malloc(Nmap*Nmap*sizeof(float));
  for(i=0;i<Nmap*Nmap;++i) map[i]=0.0;
  printf("%i\n",Nmap);
 
  /**********************************/
  /*** make initial grid ***/
  imageinfo=NewImageInfo(200);
  linkinglength=1.1*sqrt(2.)*range/(Ngrid-1);

  i_points=NewPointArray(Ngrid*Ngrid,True);
  xygridpoints(i_points,range,center,Ngrid,0);
  s_points=LinkToSourcePoints(i_points,Ngrid*Ngrid);

  /* for(i=0;i<Ngrid*Ngrid;++i) printf("%e %e\n",s_points[i].x[0],s_points[i].x[1]);*/
  /*  printTree(i_tree); */

  for(m=0;m<833;++m){   /*loop through sources */
    /*for(m=0;m<100;++m){   /* loop through sources*/

	  printf("Lensing source %i\n",m);
    /*read in source information */
      /*fscanf(file,"%i %le %le %le %le %le %le",&idnumber,&source[0],&source[1],&zsource,&angle,&amax,&amin);*/
	  fscanf(file,"%i %le %le %le %le %le %le %le %le %le %le %le",&idnumber,&source[0],&source[1],&zsource
			  ,&angle,&amax,&amin,&i_151,&i_610,&i_1400,&i_4860,&i_18000);

	  flux=1.0;
	  flux=pow(10,i_151);

    /*printf("%i %f %f %f %f %f %f\n",idnumber,source[0],source[1],zsource,angle,amax,amin);/**/

 
    /* convert to distance on the lens plane */
	  amax*=Dlens*pi/60/60/180/(1+zlens);
	  amin*=Dlens*pi/60/60/180/(1+zlens);
	  source[0]*=Dlens*pi/180/(1+zlens);
	  source[1]*=Dlens*pi/180/(1+zlens);
    /*printf("%i %f %f %f %f %f %f\n",idnumber,source[0],source[1],zsource,angle,amax,amin);/**/

	  if( fabs(source[0]) > range/2 || fabs(source[1]) > range/2 ){
		  printf("ERROR: source is out of grid range in main source = %e %e index = %i range/2=%e\n",source[0],source[1],idnumber,range/2);
		  exit(0);
	  }


	  if(zsource != zold){
		  /*** free refined points so they don't need to be lensed ***/
		  if(m>0){
			  spointlist=freeTree(s_tree);
			  ipointlist=freeTree(i_tree);

			  /* Note: an intentional memory leak is produced here because the subgrids for */
			  /*       previous sources are simply ignored */

			  /* free points not on coarse grid */
			  /* 	MoveToTopList(spointlist); */
			  /* 	MoveToTopList(ipointlist); */
			  /* 	Nframes=spointlist->Npoints; */

			  /* 	for(i=0;i<Nframes;++i){ */
			  /* 	  if(spointlist->current->id > (Ngrid*Ngrid-1) ){ */
			  /* 	    point=TakeOutCurrent(spointlist); */

			  /* /\* 	    if( ( point->id - Ngrid*Ngrid )%9 == 0 ){  *\/ */
			  /* /\* 	      if(point->id == 2500){ *\/ */
			  /* /\* 		printf(" point=%i \n",point); *\/ */
			  /* /	\* 		PrintPoint(point); *\/ */
			  /* /\* 		exit(0); *\/ */
			  /* /\* 	      } *\/ */
			  /* /\* 	      free(point); } *\/ */
			  /* 	  } */
			  /* 	  if(ipointlist->current->id > (Ngrid*Ngrid-1) ){ */
			  /* 	    point=TakeOutCurrent(ipointlist); */
			  /* 	    /\*if( ( point->id - Ngrid*Ngrid )%9 == 0 ) free(point);*\/ */
			  /* 	  } */
			  /* 	  MoveDownList(spointlist); */
			  /* 	  MoveDownList(ipointlist); */
			  /* 	} */

			  free(spointlist);
			  free(ipointlist);

			  /*** reset gridsize to original ***/
			  for(i=0;i<Ngrid*Ngrid;++i){
				  i_points[i].gridsize=range/(Ngrid-1);
				  s_points[i].gridsize=range/(Ngrid-1);
			  }

		  }

		  /*     for(i=0,imax=1.0e99,smax=1.0e99;i<Ngrid*Ngrid;++i){ */
		  /*       if(i_points[i].gridsize < imax) imax=i_points[i].gridsize; */
		  /*       if(s_points[i].gridsize < smax) smax=s_points[i].gridsize; */
      /*     } */
      /*     printf("imax=%f  smax=%f  %f\n",imax,smax,range/(Ngrid-1)); */


		  /* re-lens original grid and rebuild source tree */
		  rayshooterInternal(Ngrid*Ngrid,i_points,i_tree);

		  /** build trees image and source planes **/
		  /*printf("%i\n",s_points[Ngrid*Ngrid-1].id);*/
		  i_tree=BuildTree(i_points,Ngrid*Ngrid);
		  s_tree=BuildTree(s_points,Ngrid*Ngrid);

		  zold=zsource;
	  } else zold=zsource;

	  modify to work with fixed resolution

	  find_images(source,amax,s_tree,i_tree,&Nimages
			  ,imageinfo,&NImagePoints,0,False,verbose);

	  do{

		  /* find points on grid */
		  if( amin > 0.0){
			  moved=image_finder(source,amax,linkinglength,s_tree,i_tree
					  ,&Nimages,imageinfo,&NImagePoints,-1,0);
		  }else{
			  moved=image_finder(source,0.0,linkinglength,s_tree,i_tree
					  ,&Nimages,imageinfo,&NImagePoints,1,0);
		  }

		  /* refine grid */

	  }while( refine_grid(i_tree,s_tree,imageinfo,Nimages,range/(Ngrid-1)/pow(3.01,refinement_level))
			  || moved );

	  /* record data in the image structure */
 
	  /*   if( amin > 0.0){
		  printf("\n%i image points\n%i images\n",NImagePoints,Nimages); /**/

	  for(i=0;i<Nimages;++i){

		  printf("  Nimages=%i %e\n",Nimages,imageinfo[i].Npoints*pow(imageinfo[i].points[0].gridsize,2));
		  for(j=0;j<imageinfo[i].Npoints;++j){
			  ++Nframes;
			  /*printf("%e  %e  %e\n",imageinfo[i].points[j].x[0]*180*60/Dlens/pi
				  ,imageinfo[i].points[j].x[1]*180*60/Dlens/pi,imageinfo[i].points[j].invmag);*/

			  ix=(long)( (imageinfo[i].points[j].x[0]/resolution- (center[0] - map_range/2 )/resolution) + 0.5 );
			  iy=(long)( (imageinfo[i].points[j].x[1]/resolution- (center[1] - map_range/2 )/resolution) + 0.5 );

			  /* 	ix=(long)( (imageinfo[i].points[j].x[0]/imageinfo[0].points[0].gridsize */
			  /* 		    -(center[0] - map_range/2 )/imageinfo[0].points[0].gridsize) + 0.5 ); */
			  /* 	iy=(long)( (imageinfo[i].points[j].x[1]/imageinfo[0].points[0].gridsize */
			  /* 		    -(center[1] - map_range/2 )/imageinfo[0].points[0].gridsize) + 0.5 ); */

			  /* 	  ix=(long)( (imageinfo[i].points[j].x[0] - center[0] + map_range/2 )/imageinfo[i].points[0].gridsize + 0.5 ); */
			  /* 	  iy=(long)( (imageinfo[i].points[j].x[1] - center[1] + map_range/2 )/imageinfo[i].points[0].gridsize + 0.5 ); */
	
			  if( (ix>-1)*(ix<Nmap) && (iy>-1)*(iy<Nmap) ){
	  
				  if(flux==1.0){
					  map[ix+Nmap*iy] += 1.0;
					  }else{
						  if(amin>0.0 && imageinfo[i].Npoints > 1 ) map[ix+Nmap*iy] += flux/(pi*amax*amin) * pow(imageinfo[i].points[j].gridsize,2);
						  else map[ix+Nmap*iy] += flux/fabs(imageinfo[i].points[j].invmag);
					  }

			  }

		  }
	  }
  }
  printf(" total pixel values %i\n",Nframes);
  fclose(file);

  file=fopen(argv[2],"w");
  resolution*=180*60*60/Dlens/pi; /* convert to arcseconds */
  fwrite(&resolution,sizeof(double),1,file);
  fwrite(&Nmap,sizeof(unsigned long),1,file);
  fwrite(&Nmap,sizeof(unsigned long),1,file);
  fwrite(map,sizeof(float),Nmap*Nmap,file);
  fclose(file);

}
