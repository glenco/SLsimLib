/* makeStarField.c makes random field of stars */
/* makeStarField.x Nparticles size filename */
#define pi  3.1415926
#define Grav 4.7788e-20

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <nr.h>
#include "../../Library/Recipes/nrutil.h"
#include "../../Library/Recipes/nrutil.c"
#include <nrD.h>
#include "../../Library/RecipesD/locateD.c"
#include "../../Library/Recipes/ran2.c"
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

int main(int arg,char **argv){
  FILE *file;
  double **xp;
  unsigned long i,Nparticles;
  long seed;
  float size;

  Nparticles=(unsigned long)(atol(argv[1]));
  printf("Nparticles=%i\n",Nparticles);
  size=atof(argv[2]);

  xp=dmatrix(0,Nparticles-1,0,2);

  file=fopen(argv[3],"w");
  fwrite(&Nparticles,sizeof(unsigned long),1,file);

  for(i=0;i<Nparticles;++i){

    xp[i][0]=size*(ran2(&seed)-0.5);
    xp[i][1]=size*(ran2(&seed)-0.5);
    xp[i][2]=size*(ran2(&seed)-0.5);

    /*printf("%e  %e  %e\n",xp[i][0],xp[i][1],xp[i][2]);*/
    fwrite(xp[i],sizeof(double),3,file);
  }
  fclose(file);

  free_dmatrix(xp,0,Nparticles-1,0,2);

  return 1;
}
