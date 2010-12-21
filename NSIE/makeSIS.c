/* makeSIS.c makes a singular isothermal sphere made out */
/* of particles */
/* makeSIS.x Nparticles size filename */
#define pi  3.1415926
#define Grav 4.7788e-20

#include <math.h>
#include <stdio.h>
#include <time.h>
#include "../../Library/Recipes/nr.h"
#include "../../Library/Recipes/nrutil.h"
#include "../../Library/Recipes/nrutil.c"
#include "../../Library/RecipesD/nrD.h"
#include "../../Library/RecipesD/ran2D.c"
#include "../../Library/RecipesD/bsstepD.c"
#include "../../Library/RecipesD/dfridrD.c"
#include "../../Library/RecipesD/locateD.c"
#include "../../Library/RecipesD/polintD.c"
#include "../../Library/RecipesD/odeintD.c"
#include "../../Library/RecipesD/pzextrD.c"
#include "../../Library/RecipesD/mmidD.c"
#include "../../Library/cosmo.h"
#include "../../Library/cosmo.c"
#include "../../Library/powerCDM.c"

int main(int arg,char **argv){
  FILE *file;
  double *xp,r,theta,phi,costheta;
  unsigned long i,Nparticles;
  long seed;
  float size;

  Nparticles=(unsigned long)(atol(argv[1]));
  printf("Nparticles=%i\n",Nparticles);
  size=atof(argv[2]);

  xp=(double *)malloc(3*sizeof(double));

  printf("writing to file %s",argv[3]);
  file=fopen(argv[3],"w");
  fwrite(&Nparticles,sizeof(unsigned long),1,file);

  for(i=0;i<Nparticles;++i){

    r=ran2D(&seed)*size;
    phi=ran2D(&seed)*2*pi;
    costheta=2*ran2D(&seed)-1;

    xp[0]=r*sqrt(1-costheta*costheta)*cos(phi)*costheta/fabs(costheta);
    xp[1]=r*sqrt(1-costheta*costheta)*sin(phi)*costheta/fabs(costheta);
    xp[2]=r*costheta;

    /*printf("%e  %e  %e\n",xp[0],xp[1],xp[2]);*/
    fwrite(xp,sizeof(double),3,file);
  }
  fclose(file);

  return 1;
}
