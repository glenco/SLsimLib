#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <float.h>
#include <assert.h>

#include <nr.h>
#include <nrutil.h>
#include <nrD.h>

//#include "../../Library/RecipesD/powellD.c"
#include "../../Library/RecipesD/rfD.c"
#include "../../Library/RecipesD/rdD.c"
#include "../../Library/RecipesD/svdcmpD.c"
#include "../../Library/RecipesD/svbksbD.c"
#include "../../Library/RecipesD/zriddrD.c"
#include "../../Library/RecipesD/pythagD.c"

#include "../AnalyticNSIE/analytic_lens.h"
#include "fitlens.h"

static double betaT,*modT,**xobT,**dx_subT,sigGT,*modTT,*modoT,**vT,x_centerT[2],**xgT,**dx_subTt;
static int NmodT,NsourcesT,NimagesT,*pairingT,degenT,Nmin;
static double oldsm;//,tang[2],length,yot[2],radsourceT;

void FindLensSimple(AnaLens *lens,int Nimages,Point *image_positions,double *y,double **dx_sub){
	/* rapper that allows simple lens to be found with a single
	 * lens and translates result into data structures used in the other code
	 *
	 * the lens is centered on [0,0]
	 * source position in lens is updated along with all the modes
	 */

	if(lens->perturb_Nmodes <= 0 || Nimages <= 0){
		ERROR_MESSAGE();
		printf("must set perturb_Nmodes lens->perturb_Nmodes = %li Nimages = %i \n"
				,lens->perturb_Nmodes,Nimages);
		exit(1);
	}

	int i;

	//PrintAnaLens(lens,False,False);
	//for(i=0;i<Nimages;++i) printf("  x = %e %e\n",image_positions[i].x[0],image_positions[i].x[1]);

	if(Nimages == 1){
		for(i=1;i<lens->perturb_Nmodes;++i) lens->perturb_modes[i] = 0.0;
		y[0] = image_positions[0].x[0];
		y[1] = image_positions[1].x[1];

		return ;
	}

	int pairing[Nimages],Nsources = 1;
	double **xob,**xg,q[6],*mods;
	double re2 = 0,x_center[2],scale;

	xob = dmatrix(0,Nimages-1,0,1);
	xg = dmatrix(0,1,0,1);
	mods=dvector(0,lens->perturb_Nmodes + 2*Nsources );


	xg[0][0] = xg[0][1] = 0.0;
	x_center[0] = x_center[1] = 0.0;

	// calculate scale to re-normalize.  Otherwise the linear algebra routines will fail.
	for(i=0,scale=0;i<Nimages;++i)
		scale = DMAX(scale,sqrt( pow(image_positions[0].x[0] - image_positions[i].x[0],2)
		    	 + pow(image_positions[0].x[1] - image_positions[i].x[1],2) ) );

	for(i=0;i<Nimages;++i){

		pairing[i] = 1;
		xob[i][0] = image_positions[i].x[0]/scale;
		xob[i][1] = image_positions[i].x[1]/scale;

		dx_sub[i][0] /= scale;
		dx_sub[i][1] /= scale;
		//printf("%e %e\n",xob[i][0],xob[i][1]);
	}

	ElliptisizeLens(Nimages,Nsources,1,pairing,xob,x_center,xg,0,lens->perturb_beta
			,lens->perturb_Nmodes-1,mods,dx_sub,&re2,q);

	for(i=1;i<lens->perturb_Nmodes;++i) lens->perturb_modes[i] = mods[i];

	// source position
	y[0] = mods[i]*scale;
	y[1] = mods[i+1]*scale;

	lens->perturb_modes[0] = 0.0;
	lens->perturb_modes[1] *= -1;  // checked
	lens->perturb_modes[2] *= -1;  // checked
	for(i=3;i<lens->perturb_Nmodes;++i) lens->perturb_modes[i] *= scale; // checked

	for(i=0;i<Nimages;++i){
		dx_sub[i][0] *= scale;
		dx_sub[i][1] *= scale;
	}

	lens->host_ro = 0.0;
	lens->host_sigma = 0.0;

	free_dmatrix(xob,0,Nimages-1,0,1);
	free_dmatrix(xg,0,1,0,1);
	free_dvector(mods,0,lens->perturb_Nmodes + 2*Nsources );

	return ;
}

double ElliptisizeLens(int Nimages,int Nsources,int Nlenses,int *pairing,double **xob
		       ,double *x_center,double **xg,double sigG,double beta,int Nmod
		       ,double *mod,double **dx_sub,double *re2,double *q){


	/*************************************
	****** find most elliptical lens *****
	******
	***  Input:
	***
	****** Nimages  - number of images
	****** Nsources - number of sources
	****** Nlenses  - number of lens centers
	****** pairing  - [0...Nimages-1] which source each image belongs to
	******             , indexed 1 through Nsources
	****** xob[][]      - [0...Nimages-1][0...1]  observed image positions
	****** x_center[][] - ??? expected center of lens
	****** xg       - ??? [0...Nlenses-1][0...1] centers of additional lenses
	****** sigG     - amount by which the center of the lens is allowed to vary
	****** beta     - slope of density profile  kappa propto r^-beta
	****** Nmod     - number of modes used in lens fitting,  this must be greater then ???
	****** dx_sub[][]    - [0,Nimages-1][0,1] deflection offset for each image caused by any
	******             masses not included in the host model (ex substructure)
	******
	***  Output:
	***
	****** mod[]     - [1...Nmod+2*Nsources]  modes of final model
	****** re2       - Einstein radius of additional lens
	****** q[6]
	******      q[1]     - ellipticity of nearest elliptical lens
	******  	q[2]     - position angle of nearest elliptical lens
	****** 	  if sigG == 0
	****** 	    q[3]     - re2 of additional lens if Nlenses = 2, copied to re2
	******    if sigG > 0
	******      q[3,4]   - center of lens, copied to x_center
	******      q[5]     - re2 of additional lens if Nlenses = 2, copied to re2
	******
	*************************************
	**************************************/

	int iter,i;
	double **xi,sm;
	static int count=0;
	double s;

	assert(Nimages > 1);
	assert(Nsources > 0);
	assert(Nlenses > 0);
	for(i=0;i<Nimages;++i) assert(pairing[i] > 0 && pairing[i] <= Nimages);

	//printf("x_center = %e %e\n",x_center[0],x_center[1]);
	//printf("sigG = %e\n",sigG);

	++count;

	// allocate memory
	pairingT=ivector(0,Nimages-1);
	xobT=dmatrix(0,Nimages-1,0,1);
	dx_subT=dmatrix(0,Nimages-1,0,1);
	modT=dvector(1,Nmod+2*Nsources);
	modTT=dvector(1,Nmod+2*Nsources);
	modoT=dvector(1,Nmod+2*Nsources);
	vT=dmatrix(1,Nmod+2*Nsources,1,Nmod+2*Nsources);
	xgT=dmatrix(0,Nlenses-1,0,1);
	if(Nlenses>1) dx_subTt=dmatrix(0,Nimages-1,0,1);


	// copy input parameters into global variables
	NimagesT=Nimages;
	NsourcesT=Nsources;
	betaT=beta;
	NmodT=Nmod;
	sigGT=sigG;

	for(i=1;i<=Nmod+2*Nsources;++i) mod[i]=0.0;

	for(i=0;i<Nimages;++i){
		xobT[i][0]=xob[i][0];
		xobT[i][1]=xob[i][1];
		pairingT[i]=pairing[i];
		dx_subT[i][0]=dx_sub[i][0];
		dx_subT[i][1]=dx_sub[i][1];
		if(Nlenses>1 && count>1 ){
			s=atan2(xg[1][1]-xob[i][1],xg[1][0]-xob[i][0]);
			dx_subT[i][0]=dx_sub[i][0] + *re2*cos(s);
			dx_subT[i][1]=dx_sub[i][1] + *re2*sin(s);
		}
	}
	for(i=0;i<Nlenses;++i){
		xgT[i][0]=xg[i][0];
		xgT[i][1]=xg[i][1];
	}

	// find lens with degeneracy information
	find_lens(NimagesT,NsourcesT,pairingT,xobT,x_center,betaT,NmodT,&degenT,modT,vT,dx_subT);

	//printf("found model\n");
	for(i=1;i<=Nmod + 2*Nsources;++i) mod[i] = modT[i];

	/*
	 * find the most elliptical model amongst the degenerate models that
	 * fit the image positions
	 */

	//if(count == 1){
		q[1] = 0.9;
		q[2] = 0.0; /*find_axis(modT,NmodT);*/
		//}

	x_centerT[0] = x_center[0];
	x_centerT[1] = x_center[1];

	xi=dmatrix(1,5,1,5);

	xi[1][1]=0.01;  xi[1][2]=0.00;    xi[1][3]=0.0;       xi[1][4]=0.0;        xi[1][5]=0.0;
	xi[2][1]=0.00;  xi[2][2]=1.0e-3;  xi[2][3]=0.0;       xi[2][4]=0.0;        xi[2][5]=0.0;
	xi[3][1]=0.00;  xi[3][2]=0.0;     xi[3][3]=1.0e-3;    xi[3][4]=0.0;        xi[3][5]=0.0;
	xi[4][1]=0.00;  xi[4][2]=0.0;     xi[4][3]=0.0;       xi[4][4]=1.0e-3;     xi[4][5]=0.0;
	xi[5][1]=0.00;  xi[5][2]=0.0;     xi[5][3]=0.0;       xi[5][4]=0.0;        xi[5][5]=1.0e-3;

	oldsm=-1;
	if(sigG == 0.0){   /// keep center of lens fixed
		if(Nlenses==1 || count>1){
			Nmin=2;
			powellD(q,xi,Nmin,1.0e-18,&iter,&sm,minEllip);
		}else{
			Nmin=3;
			if(count == 1) q[3]=1.0e-5;
			powellD(q,xi,Nmin,1.0e-18,&iter,&sm,minEllip);
			*re2=q[3];
		}
	}else{
		q[3] = x_center[0];
		q[4] = x_center[1];
		if(Nlenses==1){
			Nmin=4;
			powellD(q,xi,Nmin,1.0e-18,&iter,&sm,minEllip);
			*re2=0.0;
		}else{
			Nmin=5;
			if(count == 1) q[5]=0.001;
			powellD(q,xi,Nmin,1.0e-18,&iter,&sm,minEllip);
			*re2=q[5];
		}
		x_center[0] = q[3];
		x_center[1] = q[4];
	}
	for(i=1;i<=Nmod+2*Nsources;++i) mod[i]=modTT[i];

	//printf("iter = %i\n",iter);
	free_ivector(pairingT,0,Nimages-1);
	free_dmatrix(xobT,0,Nimages-1,0,1);
	free_dmatrix(dx_subT,0,Nimages-1,0,1);
	free_dvector(modT,1,Nmod+2*Nsources);
	free_dvector(modTT,1,Nmod+2*Nsources);
	free_dvector(modoT,1,Nmod+2*Nsources);
	free_dmatrix(vT,1,Nmod+2*Nsources,1,Nmod+2*Nsources);
	free_dmatrix(xgT,0,Nlenses-1,0,1);
	if(Nlenses>1) free_dmatrix(dx_subTt,0,Nimages-1,0,1);
	free_dmatrix(xi,1,5,1,5);

	return sm;
}

double minEllip(double *par){
	// minimized to find model closest to elliptical
	double K,E,q,theta,sm,r,s;
	int i;
	static double x_center[2];

	q=par[1];   // ellipticity
	theta=par[2];  // position angle

	// set modo to elliptical model
	for(i=1;i<=NmodT;++i){
		modoT[i]=0;
	}

  // elliptical integrals
	K = rfD(0,1./q/q,1);
	E = K - (1-1./q/q)*rdD(0,1./q/q,1)/3;

	// fill in modes with their values for an elliptical lens
	modoT[3]=4*K/pi;
	if(q != 1.0){
		if(NmodT>3) modoT[4] =4*( (1+q*q)*K-2*q*q*E )/(1-q*q)/pi/(1-4);
		if(NmodT>7) modoT[8] =4*( (3*q*q+1)*(q*q+3)*K-8*q*q*(1+q*q)*E )
    		  /( 3*pi*pow(1-q*q,2) )/(1-16);
		if(NmodT>11) modoT[12] =4*( (1+q*q)*(15+98*q*q+15*q*q*q*q)*K-2*q*q*(23+82*q*q+23*q*q*q*q)*E )
    		  /( 15*pi*pow(1-q*q,3) )/(1-36);
		if(NmodT>15) modoT[16]=4*( -32*q*q*(1+q*q)*(11+74*q*q+11*q*q*q*q)*E
				+(105+1436*q*q+3062*q*q*q*q+1436*pow(q,6)+105*pow(q,8))*K )
				/(105*pi*pow(1-q*q,4))/(1-64);
	}

	// rotate model
	RotateModel(theta,modoT,NmodT,NsourcesT);   // xobT is used as a dumby variable

	x_center[0]=x_centerT[0];
	x_center[1]=x_centerT[1];

	if(Nmin == 3){
		for(i=1;i<=NimagesT;++i){  // add in second lens
			s=atan2(xgT[1][1]-xobT[i-1][1],xgT[1][0]-xobT[i-1][0]);
			dx_subTt[i-1][0]=dx_subT[i-1][0] + par[3]*cos(s);
			dx_subTt[i-1][1]=dx_subT[i-1][1] + par[3]*sin(s);
		}
		find_lens(NimagesT,NsourcesT,pairingT,xobT,x_centerT,betaT,NmodT,&degenT,modT,vT,dx_subTt);
	}
	if(Nmin == 4){
		if(par[3] != x_center[0] || par[4] != x_center[1]){
			x_center[0]=par[3];
			x_center[1]=par[4];
			find_lens(NimagesT,NsourcesT,pairingT,xobT,x_center,betaT,NmodT,&degenT,modT,vT,dx_subT);
		}
	}
	if(Nmin == 5){  /** add in second lens **/
		if(par[3] != x_center[0] || par[4] != x_center[1]){
			x_center[0]=par[3];
			x_center[1]=par[4];
		}
		for(i=1;i<=NimagesT;++i){  /** add in second lens **/
			s=atan2(xgT[1][1]-xobT[i-1][1],xgT[1][0]-xobT[i-1][0]);
			dx_subTt[i-1][0]=dx_subT[i-1][0] + par[5]*cos(s);
			dx_subTt[i-1][1]=dx_subT[i-1][1] + par[5]*sin(s);
		}
		find_lens(NimagesT,NsourcesT,pairingT,xobT,x_center,betaT,NmodT,&degenT,modT,vT,dx_subTt);
	}

	// find most elliptical model
	sm=regularize(NmodT,3,NmodT,NsourcesT,degenT,modT,vT,modoT);

	if( sm < oldsm || oldsm < 0.0){
		oldsm=sm;
		for(i=1;i<=NmodT+2*NsourcesT;++i) modTT[i]=modT[i];
	}
	//printf("q=%e   theta=%e gamma=%e %e\n",q,theta*180/pi,modT[1],modT[2]);
	//printf("%e %e %e %e %e %e\n",modoT[3],modoT[4],modoT[5],modoT[6],modoT[7],modoT[8]);
	//printf("%e %e %e %e %e %e\n\n",modT[3],modT[4],modT[5],modT[6],modT[7],modT[8]);

	// keep orientation within range
	if(theta > pi/2) sm *= exp( (theta - pi/2)/1.0e-5 );
	if(theta < -pi/2 ) sm *= exp( -(theta - pi/2)/1.0e-5 );
	// keep center of lens within sigG
	if(Nmin == 4 || Nmin == 5){
		r=pow(x_center[0]-x_centerT[0],2) + pow(x_center[1]-x_centerT[1],2);
		if(sqrt(r) > sigGT ) sm *= exp( 10*(r-sigGT*sigGT)/sigGT/sigGT );
	}
	if(fabs(log10(q)) > 2) sm *=exp( pow( (fabs(log10(q))- 2)/0.0001   ,2) );

	//sm *= exp( (modT[1]*modT[1] + modT[2]*modT[2])/pow(0.01,2) );
	return sm;
}

void find_lens(int Nimages,int Nsources,int *pairing,double **xob,double *x_center,double beta
	       ,int Nmodes,int *degen,double *mod,double **v,double **dx_sub){
	/*******************************************************
	calculate lens el that fits the image positions
	  dx_sub[] is a perturbation to the deflection angle calculated elsewhere

	 [1]=gamma1
	 [2]=gamma2
	 [3]=ao

	*******************************************************/

	double **c,*b,*w,r,theta,wmax,**a,*y,*temp,**x;
	int i,k,j;

	if(Nmodes+2*Nsources < 2*Nimages){ printf("ERROR: too few parameters in el\n"); exit(0);}

	y=dvector(0,1);
	temp=dvector(0,1);
    
	x=dmatrix(0,Nimages-1,0,1);
	c=dmatrix(1,2*Nimages,1,Nmodes+2*Nsources);
	a=dmatrix(1,2*Nimages,1,Nmodes+2*Nsources);
	b=dvector(1,2*Nimages);
	w=dvector(1,Nmodes+2*Nsources);

	betaT=beta;
	NmodT=Nmodes;

	// recenter on lens center
	for(i=0;i<Nimages;++i){
		x[i][0]=xob[i][0]-x_center[0]; x[i][1]=xob[i][1]-x_center[1];
		 //printf("%e  %e\n",x[i][0],x[i][1]);
	}

	// fill in data matrix
	for(i=1;i<=Nimages;++i){
		r = sqrt( x[i-1][0]*x[i-1][0] + x[i-1][1]*x[i-1][1] );
		theta = atan2(x[i-1][1],x[i-1][0]);

		c[2*i-1][1]= -x[i-1][0]; c[2*i-1][2]= -x[i-1][1];
		c[2*i][1]  =  x[i-1][1]; c[2*i][2]  = -x[i-1][0];
    
		c[2*i-1][3]= 0.5*beta*pow(r,beta-1)*cos(theta);
		c[2*i][3]  = 0.5*beta*pow(r,beta-1)*sin(theta);

		for(j=4;j<=Nmodes-1;j+=2){
			k=j/2;
			c[2*i-1][j]  = (beta*cos(theta)*cos(k*theta)+k*sin(theta)*sin(k*theta))*pow(r,beta-1);
			c[2*i-1][j+1]= (beta*cos(theta)*sin(k*theta)-k*sin(theta)*cos(k*theta))*pow(r,beta-1);

			c[2*i][j]  = (beta*sin(theta)*cos(k*theta)-k*cos(theta)*sin(k*theta))*pow(r,beta-1);
			c[2*i][j+1]= (beta*sin(theta)*sin(k*theta)+k*cos(theta)*cos(k*theta))*pow(r,beta-1);
		}
    
		// assign images to sources
		for(j=0;j<Nsources;++j){
			if(pairing[i-1]==j+1){
				c[2*i-1][2*j+Nmodes+1] = 1; c[2*i-1][2*j+Nmodes+2] = 0;
				c[2*i][2*j+Nmodes+1]   = 0; c[2*i][2*j+Nmodes+2]   = 1;
			}else{
				c[2*i-1][2*j+Nmodes+1] = 0; c[2*i-1][2*j+Nmodes+2] = 0;
				c[2*i][2*j+Nmodes+1]   = 0; c[2*i][2*j+Nmodes+2]   = 0;
			}
		}
		b[2*i-1] = x[i-1][0] + dx_sub[i-1][0];
		b[2*i]   = x[i-1][1] + dx_sub[i-1][1];
		/*printf("dx_sub = %e %e\n",dx_sub[i-1][0],dx_sub[i-1][1]);*/
	}

	// preserve image data matrix
	for(i=1;i<=2*Nimages;++i){
		for(j=1;j<=Nmodes+2*Nsources;++j){
			a[i][j] = c[i][j];
		}
	}

	//**** fit largest modes to first four images  ***
	// SVD decomposition
	svdcmpD(c,2*Nimages,Nmodes+2*Nsources,w,v);
	// weed out small w's
	wmax=0.0;
	*degen=0;
	for(i=1;i<=Nmodes+2*Nsources;++i) if (w[i] > wmax) wmax=w[i];
	wmax=1.0e-6*wmax;
	for(i=1;i<=Nmodes+2*Nsources;++i) if (w[i] < wmax){ w[i]=0.0; ++*degen;}
	svbksbD(c,w,v,2*Nimages,Nmodes+2*Nsources,b,mod);

	//for(i=1;i<=Nmodes + 2*Nsources;++i) printf("mod[%i]=%e\n",i,mod[i]);
	//for(i=1;i<=Nmodes + 2*Nsources;++i) printf("w[%i]=%e\n",i,w[i]);

	/* find degeneracy vectors
	* and make them the first *degen columns of v
	* mod[i] + a_1 v[i][1] + ... + a_degen v[i][degen] is a solution
	**/

	for(i=1,j=0;i<=Nmodes+2*Nsources;++i){
		if (w[i] == 0.0){
			++j;
			for(k=1;k<=Nmodes+2*Nsources;++k){
				v[k][j]=v[k][i];
			}
		}
	}

	//********************************************************************
	//for(i=1;i<=N+2*Nsources;++i) printf("%e\n",modo[i]);

	free_dvector(y,0,1);
	free_dvector(temp,0,1);
	free_dmatrix(x,0,Nimages-1,0,1);
	free_dmatrix(c,1,2*Nimages,1,Nmodes+2*Nsources);
	free_dmatrix(a,1,2*Nimages,1,Nmodes+2*Nsources);
	free_dvector(b,1,2*Nimages);
	free_dvector(w,1,Nmodes+2*Nsources);

	return ;
}


double deflect_translated(double beta,double *mod,double *x,double *y,double *mag,int Nmodes
		  ,int Nlenses,double Re2,double *x2){
	double kappa,gamma[2],dt;

  /***************************************************************
    calculate the sources position, surface density and magnification at x
          given lens model mod
          Re2 - Einstein radius of second lens
          x2 - position of second lens relative to center of lens
  ***************************************************************/

  if(mod[0] != 0.0){ERROR_MESSAGE(); printf("background kappa should be zero\n"); exit(0);}
  assert(Nlenses < 3);

  // use deflection calculator to reduce code duplication
  kappa = lens_expand(beta,mod,Nmodes,x,y,gamma,&dt);

  // translate result to convention used here

  /// changed from alpha to y,  also convention on shear is opposite
  y[0] = x[0] + 2*(x[0]*mod[1] + x[1]*mod[2]) - y[0];
  y[1] = x[1] - 2*(x[1]*mod[1] + x[0]*mod[2]) - y[1];

  mag[0] = 1 - kappa - gamma[0];
  mag[1] = 1 - kappa + gamma[0];
  mag[3] = -gamma[1];

  /*
  theta=atan2(x[1],x[0]);
  r=sqrt(x[0]*x[0] + x[1]*x[1]);
  cosx=x[0]/r;
  sinx=x[1]/r;
   */

  /*
  F=0.5*mod[3];
  F1=0;
  F2=0;
  for(i=4;i<=Nmodes-1;i+=2){
    k=i/2;
    F += mod[i]*cos(k*theta)     + mod[i+1]*sin(k*theta);
    F1+=-mod[i]*k*sin(k*theta)   + mod[i+1]*k*cos(k*theta);
    F2+=-mod[i]*k*k*cos(k*theta) - mod[i+1]*k*k*sin(k*theta);
  }


  y[0] = -pow(r,beta-1)*(beta*cosx*F - sinx*F1);
  y[1] = -pow(r,beta-1)*(beta*sinx*F + cosx*F1);

  // add shear
  y[0] += x[0] + x[0]*mod[1] + x[1]*mod[2];
  y[1] += x[1] - x[1]*mod[1] + x[0]*mod[2];
  */


  /*
  // magnification matrix in polar coordinates
  dxdr= (1+mod[1])*cosx + mod[2]*sinx 
    - (beta-1)*pow(r,beta-2)*(beta*cosx*F-sinx*F1);
  dxda=-(1+mod[1])*r*sinx + mod[2]*r*cosx
    + pow(r,beta-1)*(beta*sinx*F + (1-beta)*cosx*F1 + sinx*F2);
  dydr= (1-mod[1])*sinx + mod[2]*cosx
    - (beta-1)*pow(r,beta-2)*(beta*sinx*F+cosx*F1);
  dyda= (1-mod[1])*r*cosx - mod[2]*r*sinx 
    + pow(r,beta-1)*(-beta*cosx*F + (1-beta)*sinx*F1 - cosx*F2);
*/


/*   if(Nlenses>1){ */
/*     dxdr= cosx -Re2*( sinx2*sinx2*cosx - sinx2*cosx2*sinx )/r2; */
/*     dydr= sinx -Re2*( cosx2*cosx2*sinx - sinx2*cosx2*cosx )/r2; */
/*     dxda= -r*sinx + r*Re2*( sinx2*sinx2*sinx + sinx2*cosx2*cosx )/r2; */
/*     dyda=  r*cosx - r*Re2*( cosx2*cosx2*cosx + sinx2*cosx2*sinx )/r2; */
/*   } */

  /** actually the inverse magnification **/
  /* *mag=(dxdr*dyda - dxda*dydr)/r; */

  // add additional lens
  if(Nlenses > 1){
	  double dxdr,dxda,dydr,dyda;
	  double theta2,r2,cosx2,sinx2,r,cosx,sinx;

	  //theta=atan2(x[1],x[0]);
	  r=sqrt(x[0]*x[0] + x[1]*x[1]);
	  cosx=x[0]/r;
	  sinx=x[1]/r;

	  theta2=atan2(x2[1]-x[1],x2[0]-x[0]);
	  r2=sqrt( pow(x[0]-x2[0],2) + pow(x[1]-x2[1],2));
	  cosx2=cos(theta2);
	  sinx2=sin(theta2);

	  y[0] += Re2*cos(theta2);
	  y[1] += Re2*sin(theta2);

	  dxdr =  -Re2*( sinx2*sinx2*cosx - sinx2*cosx2*sinx )/r2;
	  dydr =  -Re2*( cosx2*cosx2*sinx - sinx2*cosx2*cosx )/r2;
	  dxda = r*Re2*( sinx2*sinx2*sinx + sinx2*cosx2*cosx )/r2;
	  dyda =-r*Re2*( cosx2*cosx2*cosx + sinx2*cosx2*sinx )/r2;

	  mag[0] += ( -x[1]*dxda + r*x[0]*dxdr )/r/r;  /** xx **/
	  mag[1] += (  x[0]*dyda + r*x[1]*dydr )/r/r;  /** yy **/
	  mag[2] += (  x[0]*dxda + r*x[1]*dxdr )/r/r;  /** xy **/

	  return kappa + 0.5*Re2/r2;
  }

  return kappa;
}

double regularize(int Nmax,int Nmin,int N,int Nsources,int degen
		  ,double *mod,double **v,double *modo){
  double Dsum,sum=0,sumold,aa,*weights;

  /*****************************************************************/
  /** find degenerate model most like modo modulo a normalization **/
  /*****************************************************************/

  int i,j;

  /*
  for(i=Nmin,sum=0.0;i<=Nmax;i+=1){
    printf("%i %e %e\n",i,mod[i],modo[i]);
  }
   */

  weights=dvector(1,degen);

/*   aa=mod[3]/modo[3]; */
/*   for(i=3;i<=Nmax;++i) modo[i]*=aa; */

  for(i=Nmin,Dsum=0;i<=Nmax;i+=1){
    if(i<=3) Dsum += pow( mod[i]-modo[i],2);
    else Dsum += pow(1-pow(i/2,2),2)*pow( mod[i]-modo[i],2);
  }
  sumold=Dsum;

  while(Dsum > 1.0e-6*sumold){

    if(modo[3] != 0.0){
    	/** renormalize model **/
    	for(i=3,sum=0.0;i<=Nmax;++i){
    		if(i<=3) sum += mod[i]*modo[i];
    		else  sum += pow(1-pow(i/2,2),2)*mod[i]*modo[i];
    	}
    	aa=sum;
    	for(i=3,sum=0.0;i<=Nmax;++i){
    		if(i<=3) sum += modo[i]*modo[i];
    		else  sum += pow(1-pow(i/2,2),2)*modo[i]*modo[i];
    	}
    	aa/=sum;

    	if(aa > 0.0) for(i=3;i<=Nmax;++i) modo[i]*=aa;
    }
    /** move in degenerate space to find best model **/
	for(j=1;j<=degen;++j){
      for(i=Nmin,sum=0.0;i<=Nmax;i+=1){
    	  if(i<=3) sum += ( mod[i] - modo[i]  )*v[i][j];
    	  else sum += pow(1-pow(i/2,2),2)*( mod[i] - modo[i]  )*v[i][j];
      }
      weights[j]=-sum;
      for(i=Nmin,sum=0.0;i<=Nmax;i+=1){
    	  if(i<=3) sum += v[i][j]*v[i][j];
    	  else sum += pow(1-pow(i/2,2),2)*v[i][j]*v[i][j];
      }
      weights[j] /= sum;

      for(i=1;i<=N+2*Nsources;++i){
    	  //printf("i=%i j=%i w=%e v=%e\n",i,j,weights[j],v[i][j]);
    	  if( !isnan(weights[j]) ) mod[i] += weights[j]*v[i][j];
      }
      if(sum < 0) printf("max found\n");
      /*printf("weights[%i]=%e\n",j,weights[j]);*/
	}

    for(i=Nmin,sum=0.0;i<=Nmax;i+=1){
      if(i<=3) sum += pow(mod[i]-modo[i],2);
      else sum += pow(1-pow(i/2,2),2)*pow(mod[i]-modo[i],2);
    }
    Dsum=fabs(sumold-sum);
    /*printf("%e\n",sumold-sum);*/
    sumold=sum;
    /*printf("Dsum=%e %e\n",Dsum,sumold);*/
  }

  /* test result */
  free_dvector(weights,1,degen);

  return sum;
}

/*
double critmin(double r){
   double kappa,mag[3],x[2],y[2],AA[4],angle[2],ang_lens[2];

  ang_lens[0]=0;
  ang_lens[1]=0;
  x[0]=r*cos(thetaT)+xt[0];
  x[1]=r*sin(thetaT)+xt[1];

  kappa=deflect_translated(1.0,modT,x,y,mag,NmodT,NlensesT,Re2T,x2T);
  if(withsub==1){
    deflection_total(x,angle,1,1,ang_lens,1,AA);
    return (mag[0]+AA[0]-1)*(mag[1]+AA[1]-1)-(mag[2]+AA[2])*(mag[2]+AA[3]);
  }else{ return mag[0]*mag[1]-mag[2]*mag[2];}
}
*/

/**************************************************
************** rotate lens model *****************
**************************************************/
void RotateModel(double thetaX,double *mod,int N,int Nsources){
  double temp[2];
  int i,k;
  /*  printf("rotate: theta = %e\n",thetaX);*/

  temp[0]=mod[1];  temp[1]=mod[2];
  mod[1]= temp[0]*cos(2*thetaX)+temp[1]*sin(2*thetaX);
  mod[2]=-temp[0]*sin(2*thetaX)+temp[1]*cos(2*thetaX);
  for(i=4;i<=N-1;i+=2){
    k=i/2;
    temp[0]=mod[i]; temp[1]=mod[i+1];
    mod[i]  = temp[0]*cos(k*thetaX)+temp[1]*sin(k*thetaX);
    mod[i+1]=-temp[0]*sin(k*thetaX)+temp[1]*cos(k*thetaX);
  }

  /** source positions **/
  for(i=0;i<Nsources;++i){
    temp[0]=mod[N+1+2*i];  temp[1]=mod[N+2+2*i];
    mod[N+1+2*i]= temp[0]*cos(thetaX)+temp[1]*sin(thetaX);
    mod[N+2+2*i]=-temp[0]*sin(thetaX)+temp[1]*cos(thetaX);
  }
}


double find_axis(double *mod,int Nmod){
	/********************************************
	************* find rotation axis ***********
	********************************************/
  double theta,ans;
  int i,k;

  NmodT=Nmod;
  for(i=1;i<=NmodT;++i) modT[i] = mod[i];
  theta = zriddrD(minaxis,0,pi/4.,1.0e-9);

  /* calc 2nd deriv */
  for(i=4,ans=0.0;i<=5;i+=2){
    k=i/2;
    /*printf("ans=%e   %e %e\n",ans,modT[i],modT[i+1]);*/
    ans += -4*k*k*modT[i]*modT[i+1]*sin(2*k*theta)
      - 2*k*k*(modT[i]*modT[i] - modT[i+1]*modT[i+1])*cos(2*k*theta);
  }

  if(ans < 0) theta += pi/4;  /* make answer a minimum */
  return theta;
}

double minaxis(double thetaX){
  double ans;
  int i,k;

  for(i=4,ans=0.0;i<=5;i+=2){
    k=i/2;
    /*printf("ans=%e   %e %e\n",ans,modT[i],modT[i+1]);*/
    ans += 2*k*modT[i]*modT[i+1]*cos(2*k*thetaX)
      - k*(modT[i]*modT[i] - modT[i+1]*modT[i+1])*sin(2*k*thetaX);
  }
  return ans;
}


int check_model(int Nimages,int Nsources,int Nlenses,int *pairing,double **xob,double *x_center
		,int Nmod,double *mod,double **xg,double Re2,double **dx_sub,double **Amag,double ytol){
  int i;
  double dy,dymax,tmp[2],kappa,mag[3],x2[2],par;


  // check parity of images
  for(i=0,par=0;i<Nimages;++i) par+=(Amag[i][0]*Amag[i][1]-Amag[i][2]*Amag[i][2])
				  /fabs(Amag[i][0]*Amag[i][1]-Amag[i][2]*Amag[i][2]);

  //printf("par=%e\n",par);
  if(fabs(par) > 1.0e-6 ) return 1;

  if(Nlenses>1){
    x2[0]=xg[1][0]-x_center[0]; x2[1]=xg[1][1]-x_center[1];
  }else{ Re2=0.0;}

  // check that there is only one source
  double **y = dmatrix(0,Nimages-1,0,1);
 for(i=0;i<Nimages;++i){
    tmp[0]=xob[i][0]-x_center[0];
    tmp[1]=xob[i][1]-x_center[1];
    kappa=deflect_translated(1.0,mod,tmp,y[i],mag,Nmod,Nlenses,Re2,x2);
    y[i][0]+=dx_sub[i][0];
    y[i][1]+=dx_sub[i][1];
    //printf("y= %e %e    %e %e  %i\n",y[i][0],y[i][1],kappa,1.0/(mag[0]*mag[1]-mag[2]*mag[2]),pairing[i-1]);
  }

  for(i=1,dymax=0;i<Nimages;++i){
    dy=sqrt( pow(y[0][0]-y[i][0],2)+pow(y[0][1]-y[i][1],2) );
    if(dy > dymax) dymax=dy;
  }

  free_dmatrix(y,0,Nimages-1,0,1);

  if(dymax > ytol) return 1;

  // check number of images
  tmp[0]=xob[1][0]-x_center[0];
  tmp[1]=xob[1][1]-x_center[1];
  //printf("number of images: %i\n",find_image_number(y[0],tmp,mod,Nmod,Nlenses,Re2,x2) );
  //printf("number of images: %i\n",magandcrit(y[0],tmp,mod,Nmod,Nlenses,Re2,x2) );
  //printf("finite mag: %e\n",finiteMag(1.0e-3,xob[0],mod,Nmod,Nlenses,Re2,x2) );
  //printf("finite mag: %e\n",finiteMag(1.0e-3,xob[1],mod,Nmod,Nlenses,Re2,x2) );
  //printf("finite mag: %e\n",finiteMag(1.0e-3,xob[2],mod,Nmod,Nlenses,Re2,x2) );
  //printf("finite mag: %e\n",finiteMag(1.0e-3,xob[3],mod,Nmod,Nlenses,Re2,x2) );

   return 0;
}
