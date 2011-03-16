/* makeSIS.c makes a singular isothermal sphere made out */
/* of particles */
/* makeSIS.x Nparticles size fy fz filename */
#define pi  3.1415926
#define Grav 4.7788e-20

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <nr.h>
#include "../../Library/Recipes/nrutil.h"
#include "../../Library/Recipes/nrutil.c"
#include <nrD.h>
#include "../../Library/RecipesD/ran2D.c"

int main(int arg,char **argv){
  FILE *file;
  double **xp,r,costheta,phi;
  unsigned long i,Nparticles;
  long seed;
  float size;

  Nparticles=(unsigned long)(atol(argv[1]));
  printf("Nparticles=%i\n",Nparticles);
  size=atof(argv[2]);

  xp=dmatrix(0,Nparticles-1,0,2);

  file=fopen(argv[5],"w");
  fwrite(&Nparticles,sizeof(unsigned long),1,file);

  for(i=0;i<Nparticles;++i){

    r=ran2D(&seed)*size;
    phi=ran2D(&seed)*2*pi;
    costheta=2*ran2D(&seed)-1;

    xp[i][0]=r*sqrt(1-costheta*costheta)*cos(phi)*costheta/fabs(costheta);
    xp[i][1]=r*sqrt(1-costheta*costheta)*sin(phi)*costheta/fabs(costheta);
    xp[i][2]=r*costheta;

    xp[i][1]*=atof(argv[3]);
    xp[i][2]*=atof(argv[4]);

    /*printf("%e  %e  %e\n",xp[i][0],xp[i][1],xp[i][2]);*/
    fwrite(xp[i],sizeof(double),3,file);
  }
  fclose(file);

  free_dmatrix(xp,0,Nparticles-1,0,2);

  return 1;
}
