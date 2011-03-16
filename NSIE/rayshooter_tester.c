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
#include "../TreeCode/TreeNB.h"

char *paramfile;

int main(int arg,char **argv){
  TreeHndl i_tree;
  Point *i_points,*s_points,*ave_point,*sig_point;
  double r,center[2],rmin,rmax,rad[2],cos2theta,sin2theta,old_r;
  unsigned long i,j,Ngrid;
  time_t to,t1;

  if(arg==2) paramfile=argv[1];
  else paramfile=NULL;

  Ngrid=70;
  //   make grid
  center[0]=0.0; center[1]=0.0;

  rmax=3.0e-2;
  rmin=1.0e-5;

  i_points=NewPointArray(Ngrid*Ngrid,True);
  log_polar_grid(i_points,rmax,rmin,center,Ngrid);
  s_points=LinkToSourcePoints(i_points,Ngrid*Ngrid);
  time(&to);
  rayshooterInternal(Ngrid*Ngrid,i_points,i_tree,False);
  time(&t1);
  printf("time for calculating %i points %f sec\n       %f sec per point\n",Ngrid*Ngrid
		  ,difftime(t1,to),difftime(t1,to)/Ngrid/Ngrid);

  ave_point=NewPointArray(1,True);
  sig_point=NewPointArray(1,True);

  // calculate variance and average of radial rings
  for(j=0;j<Ngrid;++j){
	  old_r=sqrt( pow(i_points[j*Ngrid].x[0]-center[0],2)
			  + pow(i_points[j*Ngrid].x[1]-center[1],2) );

	  ave_point->x[0]=0;
	  ave_point->x[1]=0;
	  ave_point->kappa=0;
	  ave_point->invmag=0;
	  ave_point->gamma[0]=0;
	  ave_point->gamma[1]=0;

	  sig_point->x[0]=0;
	  sig_point->x[1]=0;
	  sig_point->kappa=0;
	  sig_point->invmag=0;
	  sig_point->gamma[0]=0;
	  sig_point->gamma[1]=0;

	  for(i=0;i<Ngrid;++i){
		  r=sqrt( pow(i_points[j*Ngrid+i].x[0]-center[0],2)
		  			  + pow(i_points[j*Ngrid+i].x[1]-center[1],2) );
		  //if(r != old_r) printf("r = %e  old_r=%e\n",r,old_r);
		  rad[0]=(i_points[Ngrid*j+i].x[0]-center[0])/r;
		  rad[1]=(i_points[Ngrid*j+i].x[1]-center[1])/r;

		  cos2theta=rad[0]*rad[0]-rad[1]*rad[1];
		  sin2theta=2*rad[0]*rad[1];

		  // radial deflection
		  ave_point->x[0]+=( (i_points[Ngrid*j+i].x[0]-i_points[Ngrid*j+i].image->x[0])*rad[0]
					+ (i_points[Ngrid*j+i].x[1]-i_points[Ngrid*j+i].image->x[1])*rad[1] )/Ngrid;
		  // tangential deflection
		  ave_point->x[1]+=( (i_points[Ngrid*j+i].x[0]-i_points[Ngrid*j+i].image->x[0])*rad[1]
					- (i_points[Ngrid*j+i].x[1]-i_points[Ngrid*j+i].image->x[1])*rad[0] )/Ngrid;
		  ave_point->kappa+=i_points[Ngrid*j+i].kappa/Ngrid;
		  ave_point->invmag+=i_points[Ngrid*j+i].invmag/Ngrid;
		  // tangential shear
		  ave_point->gamma[0]+=(i_points[Ngrid*j+i].gamma[0]*cos2theta
				  + i_points[Ngrid*j+i].gamma[1]*sin2theta)/Ngrid;
		  ave_point->gamma[1]+= (-i_points[Ngrid*j+i].gamma[0]*sin2theta
				  + i_points[Ngrid*j+i].gamma[1]*cos2theta)/Ngrid;

		  // radial deflection
		  sig_point->x[0]+=pow( (i_points[Ngrid*j+i].x[0]-s_points[Ngrid*j+i].x[0])*rad[0]
					+ (i_points[Ngrid*j+i].x[1]-s_points[Ngrid*j+i].x[1])*rad[1] ,2)/(Ngrid-1);
		  // tangential deflection
		  sig_point->x[1]+=pow( (i_points[Ngrid*j+i].x[0]-s_points[Ngrid*j+i].x[0])*rad[1]
					- (i_points[Ngrid*j+i].x[1]-s_points[Ngrid*j+i].x[1])*rad[0] ,2)/(Ngrid-1);
		  sig_point->kappa+=pow(i_points[Ngrid*j+i].kappa,2)/(Ngrid-1);
		  sig_point->invmag+=pow(i_points[Ngrid*j+i].invmag,2)/(Ngrid-1);
		  // tangential shear
		  sig_point->gamma[0]+=pow(i_points[Ngrid*j+i].gamma[0]*cos2theta
				  + i_points[Ngrid*j+i].gamma[1]*sin2theta,2)/(Ngrid-1);
		  sig_point->gamma[1]+= pow(-i_points[Ngrid*j+i].gamma[0]*sin2theta
				  + i_points[Ngrid*j+i].gamma[1]*cos2theta,2)/(Ngrid-1);
	  }


	  sig_point->x[0]= sqrt(sig_point->x[0] - Ngrid*pow(ave_point->x[0],2)/(Ngrid-1));
	  sig_point->x[1]= sqrt(sig_point->x[1] - Ngrid*pow(ave_point->x[1],2)/(Ngrid-1));
	  sig_point->kappa= sqrt(sig_point->kappa - Ngrid*pow(ave_point->kappa,2)/(Ngrid-1));
	  sig_point->invmag= sqrt(sig_point->invmag - Ngrid*pow(ave_point->invmag,2)/(Ngrid-1));
	  sig_point->gamma[0]= sqrt(sig_point->gamma[0] - Ngrid*pow(ave_point->gamma[0],2)/(Ngrid-1));
	  sig_point->gamma[1]= sqrt(sig_point->gamma[1] - Ngrid*pow(ave_point->gamma[1],2)/(Ngrid-1));

	  printf("%e   %e %e   %e %e   %e %e   %e %e   %e %e   %e %e\n",r
			  ,ave_point->x[0],sig_point->x[0]
			  ,ave_point->x[1],sig_point->x[1]
			  ,ave_point->kappa,sig_point->kappa
			  ,ave_point->invmag,sig_point->invmag
			  ,ave_point->gamma[0],sig_point->gamma[0]
			  ,ave_point->gamma[1],sig_point->gamma[1]
             );
  }
  printf("radius (Mpc)         radial alpha               tangential alpha                 kapppa                    invmag                      tangential gamma                  cross gamma\n");
}
