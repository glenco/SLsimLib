/*******************************************************
   ellipticalens calculates the image possitions and magnification 
     ratios for a given elliptical power-law lens model
********************************************************/
#define pi  3.141593
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "../../Library/Recipes/nr.h"
#include "../../Library/Recipes/nrutil.h"


struct lensmodel{
  double b;        /* scale size of lens */
  double gamma;    /* magnitude of external shear */
  double thetag;   /* angle of external shear */
  double thetae;   /* angle of ellipticity */
  double e;        /* ellipticity SQUARED, e=1 is circular */
  double n;        /* radial potential power law n=0.5 is isothermal */
  double xl[2];    /* position of lens center */
  double y[2];     /* position of source */
  double a[2];     /* scale lengths of disk, a[0]>a[1] */
  double kappa_d;     /* central kappa of disk */
  double theta_disk; /* orientation of disk, long axis is horizontal for theta_disk=0 */
  double xd[2];    /* center of disk */
} modelg;

static double ys[2],xg[2],size_s,rgl,psf,sigt2,thetat,limit[2];

static int sis=1,disk=0;  /*** if sis=1 then a SIS model is used ***/
                          /*** if disk=0 then a disk is not included ***/

void ellip_exp_lens(struct lensmodel modelIn,double xtol,double **x,double **xo
		   ,double *q,double **Amag,int *check){
  double mago,magnification(double *,double *);
  int i,j,k,kmin;
  void positions(double **,double,int *);
  double s,smin,test[2];
  void dpotdx(double *,double *);

  for(i=0;i<4;++i){ xo[i][0]=x[i][0]; xo[i][1]=x[i][1]; }
  modelg=modelIn;
  if(sis==1) modelg.n=0.5;
  positions(x,1.0e-9,check); /* check=2 if it can't find a root */

  /** check to see if solutions are the same **/
  for(i=0;i<3;++i){
    for(j=i+1;j<4;++j){
      if( (pow(x[i][0]-x[j][0],2)+pow(x[i][1]-x[j][1],2)) < xtol*xtol*1.0e-2){
	/*printf("ERROR: two images are the same  ");/**/
	*check=3;
	/*return;*/
	k=i;
	kmin=j;
      }
    }
  }
  /*  dpotdx(x[3]-1,test-1);
  printf("dy= %e %e\n",test[0],test[0]);
  exit(0);*/

  if(*check==3){
    *check=0;
    i=k;
    j=kmin;
  /** find which image the double is closest to **/
    smin=pow(x[i][0]-xo[0][0],2)+pow(x[i][1]-xo[0][1],2);
    kmin=0;
    for(k=1;k<4;++k){  
      s=pow(x[i][0]-xo[k][0],2)+pow(x[i][1]-xo[k][1],2);
      if(s<smin){
	smin=s;
	kmin=k;
      }
    }
    /*printf("kmin=%i\n",kmin);
    printf("j=%i\n",j);
    printf("i=%i\n",i);*/
    if(kmin!=j) i=j;
    x[i][0]=xo[i][0]+0.02*(xo[i][0]-x[kmin][0]);
    x[i][1]=xo[i][1]+0.02*(xo[i][1]-x[kmin][1]);
    positions(x,1.0e-6,check);

    for(i=0;i<3;++i){
      for(j=i+1;j<4;++j){
	if( sqrt(pow(x[i][0]-x[j][0],2)+pow(x[i][1]-x[j][1],2)) < 1.0e-9){
	  /*printf("ERROR: two images are the same  ");/**/
	  *check=3;
	  return;
	}
      }
    }
  }

  mago=magnification(x[3],Amag[3]);
  for(i=0;i<4;++i){
    q[i]=magnification(x[i],Amag[i])/mago;/**/
  }
}

void positions(double **x,double xtol,int *check){
  int i,iter;
  double fret,potential(double *),test[2];
  void dpotdx2(int,double *,double *),dpotdx(double *,double *),magmatrix(double *,double **),
    rootfinder(double x[],int,int *,double,double,
		void (*vecfunc)(int,double [],double []),  /* the vector function */
		void (*jacob)(double [],double **));

  /** loop through images **/

  for(i=0;i<4;++i){ 
    rootfinder2D(x[i]-1,2,check,1.0e-6,1.0e-9,dpotdx2,magmatrix); /**/
    dpotdx(x[i]-1,test-1);
    if(*check==1 || *check==2){
      ERROR_MESSAGE();
      printf("bad model: ");/**/
      printf("\nERROR: rootfinder failed: x[%i]=%e %e\nf=%e %e\n",i,x[i][0],x[i][1]
	,test[0],test[1]);/**/
      /*exit(0);*/
      break;
    }else if( sqrt(test[0]*test[0]+test[1]*test[1]) > 1.0e-6 ){
      /*printf("\nERROR: rootfinder failed: x[%i]=%e %e\nf=%e %e check=%i\n",i,x[i][0],x[i][1]
	     ,test[0],test[1],*check);/**/
      *check=4;
      break;
      /*exit(0);*/
    }
  }
}

/*** the magnification at possition x ***/
double magnification(double *x,double *Amag){
  double **AA,mu;
  void magmatrix(double *,double **);

  AA=dmatrix(1,2,1,2);
  magmatrix(x-1,AA);

  mu=1.0/(AA[1][1]*AA[2][2]-AA[1][2]*AA[1][2]);

  Amag[0]=AA[1][1];
  Amag[1]=AA[1][2];
  Amag[2]=AA[2][2];

  /*printf("*AA[1][1]=%e AA[2][2]=%e AA[1][2]=%e AA[2][1]=%e\n",AA[1][1]
    ,AA[2][2],AA[1][2],AA[2][1]);/**/
  /*printf("mu=%e\n",mu);/**/
  free_dmatrix(AA,1,2,1,2);

  return mu;
}

void dpotdx2(int n,double *x,double *dpot){ /** to fit rootfinder format **/
  void dpotdx(double *,double *);
  double temp;

  temp=*dpot;
  dpotdx(x,dpot);
}

/*** derivative of potential for use in positions, returns dpot 
    - the source position ***/
void dpotdx(double *x,double *dpot){
  double r,xtemp[2],temp;
  double re2,qx,qy;
  double disk_pot(double *);
  int sigx,sigy;

  xtemp[0]=x[1]-modelg.xl[0];
  xtemp[1]=x[2]-modelg.xl[1];

  /** isothermal sphere **
  r=sqrt(xtemp[0]*xtemp[0]+xtemp[1]*xtemp[1]);
  dpot[1]=x[1]-modelg.y[0]-xtemp[0]*modelg.b/r;
  dpot[2]=x[2]-modelg.y[1]-xtemp[1]*modelg.b/r;
  */

  qx=pow(cos(modelg.thetae),2)+modelg.e*pow(sin(modelg.thetae),2);
  qy=pow(sin(modelg.thetae),2)+modelg.e*pow(cos(modelg.thetae),2);

  re2=qx*xtemp[0]*xtemp[0]+qy*xtemp[1]*xtemp[1]
      +(1-modelg.e)*sin(2*modelg.thetae)*xtemp[0]*xtemp[1];

  /** elliptical power law model **/
  dpot[1]=x[1]-modelg.y[0]-modelg.n*modelg.b*pow(re2,modelg.n-1)*(2*qx*xtemp[0]+(1-modelg.e)*sin(2*modelg.thetae)*xtemp[1]);
  dpot[2]=x[2]-modelg.y[1]-modelg.n*modelg.b*pow(re2,modelg.n-1)*(2*qy*xtemp[1]+(1-modelg.e)*sin(2*modelg.thetae)*xtemp[0]);

  /*** add external shear ***/
  dpot[1]+=-modelg.gamma*(xtemp[0]*cos(2*modelg.thetag)+xtemp[1]*sin(2*modelg.thetag));
  dpot[2]+=-modelg.gamma*(xtemp[0]*sin(2*modelg.thetag)-xtemp[1]*cos(2*modelg.thetag));

  /*** add disk ***/
  if(disk != 0){
    xtemp[0]=x[1]-modelg.xd[0];
    xtemp[1]=x[2]-modelg.xd[1];
    r=disk_pot(xtemp); /* potential */

    sigx=1;
    sigy=1;
    if(xtemp[0]<0) sigx=-1;
    if(xtemp[1]<0) sigy=-1;

    dpot[1]+=( -sigx*cos(modelg.theta_disk)/modelg.a[0]
	       +sigy*sin(modelg.theta_disk)/modelg.a[1])*r;

    dpot[2]+=-( sigx*sin(modelg.theta_disk)/modelg.a[0]
		+sigy*cos(modelg.theta_disk)/modelg.a[1])*r;
  }
  /*************************************************/
}

/** the second derivatives of the potential for use in positions **/
void magmatrix(double *x,double **AA){
  double r,xtemp[2],temp;
  double re2,qx,qy;
  double disk_pot(double *);
  void PrintModel(struct lensmodel,int);
  int sigx,sigy;

  xtemp[0]=x[1]-modelg.xl[0];
  xtemp[1]=x[2]-modelg.xl[1];

  /** isothermal sphere **
  r=sqrt(xtemp[0]*xtemp[0]+xtemp[1]*xtemp[1]);
  AA[1][1]=1.0-modelg.b*(1.0-pow(xtemp[0]/r,2))/r;
  AA[2][2]=1.0-modelg.b*(1.0-pow(xtemp[1]/r,2))/r;
  AA[1][2]=modelg.b*xtemp[0]*xtemp[1]/pow(r,3);
  AA[2][1]=AA[1][2]; */
  /** elliptical power law model **/

  qx=pow(cos(modelg.thetae),2)+modelg.e*pow(sin(modelg.thetae),2);
  qy=pow(sin(modelg.thetae),2)+modelg.e*pow(cos(modelg.thetae),2);

  re2=qx*xtemp[0]*xtemp[0]+qy*xtemp[1]*xtemp[1]
      +(1-modelg.e)*sin(2*modelg.thetae)*xtemp[0]*xtemp[1];
  /*printf("re2=%e qx=%e qy=%e modelg.e=%e\n",re2,qx,qy,modelg.e);/**/
  AA[1][1]=1.0-modelg.n*modelg.b*pow(re2,modelg.n-2)*
    ((modelg.n-1)*pow(2*qx*xtemp[0]+(1-modelg.e)*sin(2*modelg.thetae)*xtemp[1],2)
     +re2*2*qx);

  AA[2][2]=1.0-modelg.n*modelg.b*pow(re2,modelg.n-2)*
    ((modelg.n-1)*pow(2*qy*xtemp[1]+(1-modelg.e)*sin(2*modelg.thetae)*xtemp[0],2)
     +re2*2*qy);

  AA[1][2]=-modelg.n*modelg.b*pow(re2,modelg.n-2)*
    ((modelg.n-1)*(2*qy*xtemp[1]+(1-modelg.e)*sin(2*modelg.thetae)*xtemp[0])*
     (2*qx*xtemp[0]+(1-modelg.e)*sin(2*modelg.thetae)*xtemp[1])
     +re2*(1-modelg.e)*sin(2*modelg.thetae));

  /*** add external shear ***/
  AA[1][1]+=-modelg.gamma*cos(2*modelg.thetag);
  AA[2][2]+= modelg.gamma*cos(2*modelg.thetag);
  AA[1][2]+=-modelg.gamma*sin(2*modelg.thetag);

  /*** add disk ***/
  if(disk != 0){
    xtemp[0]=x[1]-modelg.xd[0];
    xtemp[1]=x[2]-modelg.xd[1];
    r=disk_pot(xtemp); /* potential */
    sigx=1;
    sigy=1;
    if(xtemp[0]<0) sigx=-1;
    if(xtemp[1]<0) sigy=-1;

    AA[1][1]+=(pow(cos(modelg.theta_disk)/modelg.a[0],2)
	       +pow(sin(modelg.theta_disk)/modelg.a[1],2)
	       -sigx*sigy*0.5*sin(2*modelg.theta_disk)/modelg.a[0]/ modelg.a[1] )*r;
    AA[2][2]+=(pow(sin(modelg.theta_disk)/modelg.a[0],2)
	       +pow(cos(modelg.theta_disk)/modelg.a[1],2)
	       +sigx*sigy*0.5*sin(2*modelg.theta_disk)/modelg.a[0]/ modelg.a[1] )*r;

    AA[1][2]+=(0.5*sin(2*modelg.theta_disk)*(1.0/modelg.a[0]/modelg.a[0]
					     +1.0/modelg.a[1]/modelg.a[1])
	       +sigx*sigy*cos(2*modelg.theta_disk)/modelg.a[0]/ modelg.a[1] )*r;
  }
  /***************************************/

  AA[2][1]=AA[1][2];

  /*printf("AA[1][1]=%e AA[2][2]=%e AA[1][2]=%e AA[2][1]=%e\n",AA[1][1],AA[2][2],AA[1][2],AA[2][1]);/**/
  /*PrintModel(modelg,1);*/
}

void findy(struct lensmodel modelIn,double *x,double *y){
  modelg=modelIn;
  if(sis==1) modelg.n=0.5;
  modelg.y[0]=0.0;
  modelg.y[1]=0.0;
  dpotdx(x-1,y-1);
}

// disk potential
double disk_pot(double *x){

  return 2*modelg.kappa_d*pow(modelg.a[0]*modelg.a[1],2)
    *exp(-fabs(x[0]*cos(modelg.theta_disk)+x[1]*sin(modelg.theta_disk))/modelg.a[0]
	 -fabs(x[1]*cos(modelg.theta_disk)-x[0]*sin(modelg.theta_disk))/modelg.a[1]
	 )/(modelg.a[0]*modelg.a[0]+modelg.a[1]*modelg.a[1]);
}


double finite_mag(struct lensmodel modelIn,double ap_rad,double source_rad,double *xo,double *xs,double *muo){
  double A[3];
  double dfinite_dr(double),dfinite_dtheta(double);
  double dfinite_dr2(double),dfinite_dtheta2(double);
  double nintegrateD(double (*func)(double),double,double,double);
  void dpotdx(double *x,double *dpot);
  
  double theta,z[2],y[2];
  /*printf("xo==%e %e\n",xo[0],xo[1]);*/

  modelg=modelIn;
  if(sis==1) modelg.n=0.5;
  psf=ap_rad;
  size_s=source_rad;
  xg[0]=xo[0]; xg[1]=xo[1];
  dpotdx(xs-1,ys-1); /* ys is measured relative to modelg.y */

  *muo=magnification(xs,A);
  /*printf("mu_o= %e ap_rad= %e\n",fabs(magnification(xo,A)),ap_rad);/**/

  sigt2=1.0/(1.0/ap_rad/ap_rad + 1.0/source_rad/source_rad);
  /*
for(theta=0;theta<=2*pi;theta+=pi/100.){
z[0]=xg[0]+ap_rad*cos(theta);
z[1]=xg[1]+ap_rad*sin(theta);

dpotdx(z-1,y-1);
printf("%e %e\n",y[0],y[1]);
}
printf("\n");
return ;
  */
  /*return nintegrateD(dfinite_dtheta2,0.0,6.28319,1.0e-8);*/
  /*return nintegrateD(dfinite_dr,sigt2,0,1.0e-7);/**/
  return 0.5*nintegrateD(dfinite_dr,0.0,ap_rad*ap_rad,1.0e-7);
  /*(2*pi*ap_rad*ap_rad);*/
}

double dfinite_dr(double r2){
  double dfinite_dtheta(double);
  double nintegrate2D(double (*func)(double),double,double,double);
  int i;

  /*rgl=sqrt(-2*sigt2*log(x/sigt2));*/
    rgl=sqrt(r2);
  /*printf("rgl=%e \n",rgl);/**/
  /*return pow(rgl,4);/**/

  return nintegrate2D(dfinite_dtheta,0.0,6.28319,1.0e-7);
  /*return exp(-0.5*rgl*rgl/psf/psf + 0.5*rgl*rgl/sigt2)*nintegrate2D(dfinite_dtheta,0.0,6.28319,1.0e-9);*/
}

double dfinite_dtheta(double theta){
  double z[2],y[2];
  void dpotdx(double *z,double *dpot);
  double source_profile(double *y),new;

  z[0]=xg[0]+rgl*cos(theta);
  z[1]=xg[1]+rgl*sin(theta);

  dpotdx(z-1,y-1);
  y[0]-=ys[0];  y[1]-=ys[1];

  /*if(rgl<size_s/10.) printf("source_profile(y)=%e\n",source_profile(y));/**/
  /*printf("y= %e %e dz= %e %e\n",y[0],y[1],z[0]-xg[0],z[1]-xg[1]);/**/
  /*printf("y= %e %e z= %e %e rgl=%e theta=%f\n",y[0],y[1],z[0],z[1],rgl,theta*180/pi);/**/
  /*printf("y= %e %e\n",y[0],y[1]);/**/
  return source_profile(y);
}
double dfinite_dtheta2(double theta){
  double dfinite_dr2(double);
  double nintegrate2D(double (*func)(double),double,double,double);

  thetat=theta;
  return nintegrate2D(dfinite_dr2,psf*psf*exp(-0.5*pow(3,2)),psf*psf,1.0e-7);
  return nintegrate2D(dfinite_dr2,0,4*psf*psf,1.0e-7);
}

double dfinite_dr2(double x){
  double z[2],y[2];
  void dpotdx(double *z,double *dpot);
  double source_profile(double *y);

  rgl=psf*sqrt(-2*log(x/psf/psf));/**/
  /*rgl=sqrt(r2); /**/
  z[0]=xg[0]+rgl*cos(thetat);
  z[1]=xg[1]+rgl*sin(thetat);

  dpotdx(z-1,y-1);
  y[0]-=ys[0];  y[1]-=ys[1];

  return source_profile(y);
  return exp(-0.5*rgl*rgl/psf/psf)*source_profile(y);
}

double source_profile(double *y){
  double rad,g;

  rad=sqrt(y[0]*y[0]+y[1]*y[1]);

  /*return pow(1+rad/size_s,-2);/**

  if(rad>size_s) return 0.0;
  else return 1.0/(pi*size_s*size_s);/**/

  return exp(-0.5*rad*rad/size_s/size_s)/(2.0*pi*size_s*size_s);/**/
  return exp(-0.5*rad/size_s)/(8*pi*size_s*size_s);/**/
  g=-3.0;  /* g<-2 */
  return (g+1)*(g+2)*pow(1+rad/size_s,g)/(2.0*pi*size_s*size_s);
}

/******************************************************************************
Model printing routines
/*****************************************************************************/
void PrintModel(struct lensmodel modelx,int label){
  if(label==1) printf("\n");
  printf("%.10e %.10e %.10e %.10e %.10e %.10e\n",modelx.b,modelx.gamma,modelx.thetag
	 ,modelx.thetae,modelx.e,modelx.n);
  if(label==1) printf("     b               gammma             thetag             thetae              e                 n\n");

  printf("%.10e %.10e %.10e %.10e\n",modelx.xl[0],modelx.xl[1]
	 ,modelx.y[0],modelx.y[1]);
  if(label==1) printf("    xl[0]             xl[1]             y[0]             y[1]\n");

  printf("%.10e %.10e %.10e %.10e %.10e %.10e\n",modelx.a[0],modelx.a[1],modelx.kappa_d
	 ,modelx.theta_disk,modelx.xd[0],modelx.xd[1]);
  if(label==1) printf("     a[0]             a[1]           kappa_d         theta_disk           xd[0]               xd[1]\n\n");
}


