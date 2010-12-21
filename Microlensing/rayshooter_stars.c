#define Nstars 240
#define re 0.1
#define kappa_background 0
#define gamma_background 0

void rayshooter(double *ray,double *alpha,double *gamma,double *kappa,double *invmag){
  static long i,seed=19293;
  static short count=0;
  static double **xstars;
  double radius2,length=4.0;

  if(count == 0){
    xstars=dmatrix(0,Nstars-1,0,1);
    for(i=0;i<Nstars;++i){
      xstars[i][0]=length*(ran2(&seed)-0.5);
      xstars[i][1]=length*(ran2(&seed)-0.5);
    }
    ++count;
  }

  alpha[0] = ray[0]*(kappa_background+gamma_background);
  alpha[1] = ray[1]*(kappa_background-gamma_background);

  *kappa = kappa_background;
  gamma[0] = gamma_background;
  gamma[1] = 0;

  for(i=0;i<Nstars;++i){
    radius2 = pow(ray[0]-xstars[i][0],2) + pow(ray[1]-xstars[i][1],2);
    alpha[0] += re*re*(ray[0]-xstars[i][0])/radius2;
    alpha[1] += re*re*(ray[1]-xstars[i][1])/radius2;

    gamma[0] += -re*re*( pow(ray[0]-xstars[i][0],2)-pow(ray[1]-xstars[i][1],2) )/radius2/radius2;
    gamma[1] += -2*re*re*(ray[0]-xstars[i][0])*(ray[1]-xstars[i][1])/radius2/radius2;
  }

  *invmag = pow(1-*kappa,2) - gamma[0]*gamma[0] - gamma[1]*gamma[1];
}

#undef Nstars
#undef re
#undef kappa_background
#undef gamma_background
