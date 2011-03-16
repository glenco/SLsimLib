#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "TreeNB.h"
#include "../../Library/cosmo.h"
#include "../../Library/Recipes/nrutil.h"

SimLens *readparams(char *filename,CosmoHndl cosmo){
  FILE *file;
  char label[20];
  SimLens lens, *lenses;
  int i,j;

  file=fopen(filename,"r");
  fscanf(file,"%s %s",label,lens.simfilename);
  printf("%s %s\n",label,lens.simfilename);
  fscanf(file,"%s %s",label,lens.treefilenames);
  printf("%s %s\n",label,lens.treefilenames);

  fscanf(file,"%s %i",label,&(lens.Nsph));
  printf("%s %i\n",label,lens.Nsph);
  fscanf(file,"%s %i",label,&(lens.Nspecies));
  printf("%s %i\n",label,lens.Nspecies);

  fscanf(file,"%s %le",label,&(lens.theta_force));
  printf("%s %f\n",label,lens.theta_force);
  fscanf(file,"%s %le",label,&(lens.interpolation_scale));
  printf("%s %f\n",label,lens.interpolation_scale);
  fscanf(file,"%s %le",label,&(lens.theta));
  printf("%s %f\n",label,lens.theta);
  fscanf(file,"%s %le",label,&(lens.phi));
  printf("%s %f\n",label,lens.phi);

  fscanf(file,"%s %le",label,&(lens.mass_units));
  printf("%s %.3e\n",label,lens.mass_units);

  fscanf(file,"%s %le",label,&(lens.zlens));
  printf("%s %f\n",label,lens.zlens);
  fscanf(file,"%s %le",label,&(lens.zsource));
  printf("%s %f\n",label,lens.zsource);

  SetConcordenceCosmology(cosmo);
  cosmo->physical=0;

  fscanf(file,"%s %le",label,&(cosmo->Omo));
  printf("%s %f\n",label,cosmo->Omo);
  fscanf(file,"%s %le",label,&(cosmo->Oml));
  printf("%s %f\n",label,cosmo->Oml);
  fscanf(file,"%s %le",label,&(cosmo->h));
  printf("%s %f\n",label,cosmo->h);
  fclose(file);



  // make copies for other particle species

   /********************************/
  /* set some internal parameters */
  /********************************/

  /* matrix for rotation of the coordinates */
  /* 1st around x-axis by theta and then around y-axis by phi */
  lens.coord=dmatrix(0,2,0,2);

  lens.coord[0][0]=cos(lens.phi);
  lens.coord[0][1]=sin(lens.phi)*sin(lens.theta);
  lens.coord[0][2]=sin(lens.phi)*cos(lens.theta);

  lens.coord[1][0]=0.0;
  lens.coord[1][1]=cos(lens.theta);
  lens.coord[1][2]=-sin(lens.theta);

  lens.coord[2][0]=-sin(lens.phi);
  lens.coord[2][1]=cos(lens.phi)*sin(lens.theta);
  lens.coord[2][2]=cos(lens.phi)*cos(lens.theta);

  /* find critical density */
  lens.Sigma_crit=angDist(0,lens.zsource,cosmo)/angDist(lens.zlens,lens.zsource,cosmo)
		  /angDist(0,lens.zlens,cosmo)/4/pi/Grav;
  printf("critical density is %e Msun/Mpc^2\n",lens.Sigma_crit);

  lenses=(SimLens *)malloc(lens.Nspecies*sizeof(SimLens));

   for(i=0;i<lens.Nspecies;++i){
 	  lenses[i].zlens=lens.zlens;
 	  lenses[i].zsource=lens.zsource;
 	  lenses[i].Nspecies=lens.Nspecies;
 	  lenses[i].phi=lens.phi;
 	  lenses[i].theta=lens.theta;
 	  lenses[i].mass_units=lens.mass_units;
 	  lenses[i].theta_force=lens.theta_force;
 	  lenses[i].Nsph=lens.Nsph;
 	  lenses[i].interpolation_scale=lens.interpolation_scale;
 	  for(j=0;j<50;++j) lenses[i].treefilenames[j]=lens.treefilenames[j];
 	  for(j=0;j<50;++j) lenses[i].simfilename[j]=lens.simfilename[j];

	  lenses[i].Sigma_crit=lens.Sigma_crit;

	  lenses[i].coord=dmatrix(0,2,0,2);

	  lenses[i].coord[0][0]= lens.coord[0][0];
	  lenses[i].coord[0][1]= lens.coord[0][1];
	  lenses[i].coord[0][2]= lens.coord[0][2];

	  lenses[i].coord[1][0]= lens.coord[1][0];
	  lenses[i].coord[1][1]= lens.coord[1][1];
	  lenses[i].coord[1][2]= lens.coord[1][2];

	  lenses[i].coord[2][0]= lens.coord[2][0];
	  lenses[i].coord[2][1]= lens.coord[2][1];
	  lenses[i].coord[2][2]= lens.coord[2][2];
   }
   return lenses;
}
void PrintSimLens(SimLens *lens){

  printf("param file %s\n",lens->simfilename);
  printf("tree file %s\n",lens->treefilenames);

  printf("Nsph: %i\n",lens->Nsph);
  printf("Nspecies: %i\n",lens->Nspecies);

  printf("theta_force %f\n",lens->theta_force);
  printf("interpolations_scale: %f\n",lens->interpolation_scale);
  printf("theta: %f\n",lens->theta);
  printf("phi: %f\n",lens->phi);

  printf("mass_units %.3e\n",lens->mass_units);

  printf("zlens %f\n",lens->zlens);
  printf("zsource %f\n",lens->zsource);

  printf("critical density is %e Msun/Mpc^2\n",lens->Sigma_crit);

}

void readpositions(SimLens *lens){
  FILE *file;
  IndexType i;

  file=fopen(lens->simfilename,"r");
  fread(&(lens->Nparticles),sizeof(IndexType),1,file);
  printf("    number of particles to be read:  %li million\n",lens->Nparticles/1000000);
  lens->particles=(IndexType *)malloc(lens->Nparticles*sizeof(IndexType));
  lens->xp=PosTypeMatrix(0,lens->Nparticles-1,0,2);
  for(i=0;i<lens->Nparticles;++i){
	  fread(lens->xp[i],sizeof(PosType),3,file);
	  lens->particles[i]=i;
  }

  fclose(file);
  lens->rsph=(float *)malloc(lens->Nparticles*sizeof(float));
}

#define NR_END 1

PosType **PosTypeMatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a PosType matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	PosType **m;

	/* allocate pointers to rows */
	m=(PosType **) malloc((size_t)((nrow+NR_END)*sizeof(PosType*)));
	if (!m) {ERROR_MESSAGE(); printf("ERROR: PosTypeMatrix\n  allocation failure 1\n\n"); exit(0);}
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(PosType *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(PosType)));
	if (!m[nrl]) {ERROR_MESSAGE(); printf("ERROR: PosTypeMatrix\n  allocation failure 2\n\n"); exit(0);}
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}
void free_PosTypeMatrix(PosType **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((char*) (m[nrl]+ncl-NR_END));
	free((char*) (m+nrl-NR_END));
}

#undef NR_END
