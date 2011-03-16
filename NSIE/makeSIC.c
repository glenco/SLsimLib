/* makeSIS.c makes a singular isothermal sphere made out */
/* of particles */
/* makeSIS.x Nparticles size fy filename */
#define pi  3.1415926
#define Grav 4.7788e-20

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <nr.h>
#include "../../Library/Recipes/nrutil.h"
//#include "../../Library/Recipes/nrutil.c"
#include <nrD.h>
//#include "../../Library/RecipesD/ran2D.c"

int main(int arg,char **argv){
  FILE *file;
  double **xp,r,theta,phi;
  unsigned long i,Nparticles;
  long seed;
  float size;

  Nparticles=(unsigned long)(atol(argv[1]));
  printf("Nparticles=%i\n",Nparticles);
  size=atof(argv[2]);

  xp=dmatrix(0,Nparticles-1,0,2);

  file=fopen(argv[4],"w");
  fwrite(&Nparticles,sizeof(unsigned long),1,file);

  for(i=0;i<Nparticles;++i){

    r=ran2D(&seed)*size;
    phi=ran2D(&seed)*2*pi;

    xp[i][0]=r*cos(phi);
    xp[i][1]=r*sin(phi)*atof(argv[3]);
    xp[i][2]=ran2D(&seed);

    /*printf("%e  %e  %e\n",xp[i][0],xp[i][1],xp[i][2]);*/
    fwrite(xp[i],sizeof(double),3,file);
  }
  fclose(file);
  free_dmatrix(xp,0,Nparticles-1,0,2);

  return 1;
}
