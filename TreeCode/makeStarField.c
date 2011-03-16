/* makeStarField.c makes random field of stars */
/* makeStarField.x Nparticles size filename */
#define pi  3.1415926
#define Grav 4.7788e-20

#include <math.h>
#include <stdio.h>
#include <time.h>
/*#include "../Library/Include/nr.h"*/
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

  xp=PosTypeMatrix(0,Nparticles-1,0,2);

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
