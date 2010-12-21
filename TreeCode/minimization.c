#define NRANSI
#define TINY 1.0e-25
#define ITMAX 5000
#define ITMAX1 500
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY1 1.0e-20
#define TOL 2.0e-6
#include "TreeNB.h"

int ncom;
double *pcom,*xicom,(*nrfunc)(SimLens *,double []);

void min2d(SimLens *lens,double p[], double **xi, int n, double ftol, int *iter, double *fret,
	   double (*func)SimLens *lens,double []))
{
  void linmind(SimLens *,double p[], double xi[], int n, double *fret,
	       double (*func)(SimLens *,double []));
  int i,ibig,j;
  double del,fp,fptt,t,*pt,*ptt,*xit;

  pt=dvector(1,n);
  ptt=dvector(1,n);
  xit=dvector(1,n);
  *fret=(*func)(lens,p);
  for (j=1;j<=n;j++) pt[j]=p[j];
  for (*iter=1;;++(*iter)) {
    fp=(*fret);
    ibig=0;
    del=0.0;
    for (i=1;i<=n;i++) {
      for (j=1;j<=n;j++) xit[j]=xi[j][i];
      fptt=(*fret);
      linmind(lens,p,xit,n,fret,func);
      if (fptt-(*fret) > del) {
	del=fptt-(*fret);
	ibig=i;
      }
    }
    if (2.0*(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))+TINY) {
      free_dvector(xit,1,n);
      free_dvector(ptt,1,n);
      free_dvector(pt,1,n);
      return;
    }
    if (*iter == ITMAX) nrerror("powell exceeding maximum iterations.");
    for (j=1;j<=n;j++) {
      ptt[j]=2.0*p[j]-pt[j];
      xit[j]=p[j]-pt[j];
      pt[j]=p[j];
    }
    fptt=(*func)(lens,ptt);
    if (fptt < fp) {
      t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
      if (t < 0.0) {
	linmind(lens,p,xit,n,fret,func);
	for (j=1;j<=n;j++) {
	  xi[j][ibig]=xi[j][n];
	  xi[j][n]=xit[j];
	}
      }
    }
  }
}

/* note #undef's at end of file */

void linmind(SimLens *lens,double p[], double xi[], int n, double *fret, double (*func)(SimLens *,double []))
{
  double brentd(SimLens *,double ax, double bx, double cx,
		double (*f)(SimLens *,double), double tol, double *xmin);
  double f1dimd(SimLens *,double x);
  void mnbrakd(SimLens *,double *ax, double *bx, double *cx, double *fa, double *fb,
	       double *fc, double (*func)(SimLens *,double));
  int j;
  double xx,xmin,fx,fb,fa,bx,ax;

  ncom=n;
  pcom=dvector(1,n);
  xicom=dvector(1,n);
  nrfunc=func;
  for (j=1;j<=n;j++) {
    pcom[j]=p[j];
    xicom[j]=xi[j];
  }
  ax=0.0;
  xx=1.0;
  mnbrakd(lens,&ax,&xx,&bx,&fa,&fx,&fb,f1dimd);
  *fret=brentd(lens,ax,xx,bx,f1dimd,TOL,&xmin);
  for (j=1;j<=n;j++) {
    xi[j] *= xmin;
    p[j] += xi[j];
  }
  free_dvector(xicom,1,n);
  free_dvector(pcom,1,n);
}

void mnbrakd(SimLens *lens,double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
	     double (*func)(SimLens *,double))
{
  double ulim,u,r,q,fu,dum;
  
  *fa=(*func)(lens,*ax);
  *fb=(*func)(lens,*bx);
  if (*fb > *fa) {
    SHFT(dum,*ax,*bx,dum)
      SHFT(dum,*fb,*fa,dum)
      }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc=(*func)(lens,*cx);
  while (*fb > *fc) {
    r=(*bx-*ax)*(*fb-*fc);
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
      (2.0*SIGN(DMAX(fabs(q-r),TINY1),q-r));
    ulim=(*bx)+GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0) {
      fu=(*func)(lens,u);
      if (fu < *fc) {
	*ax=(*bx);
	*bx=u;
	*fa=(*fb);
	*fb=fu;
	return;
      } else if (fu > *fb) {
	*cx=u;
	*fc=fu;
	return;
      }
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(lens,u);
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu=(*func)(lens,u);
      if (fu < *fc) {
	SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
	  SHFT(*fb,*fc,fu,(*func)(lens,u))
	  }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u=ulim;
      fu=(*func)(lens,u);
    } else {
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(lens,u);
    }
    SHFT(*ax,*bx,*cx,u)
      SHFT(*fa,*fb,*fc,fu)
      }
}

double brentd(Lens *lens,double ax, double bx, double cx, double (*f)(Lens *,double), double tol,
	double *xmin)
{
	int iter;
	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(lens,x);
	for (iter=1;iter<=ITMAX1;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*f)(lens,u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	nrerror("Too many iterations in brentd");
	*xmin=x;
	return fx;
}
/* note #undef's at end of file */

double f1dimd(SimLens *lens,double x)
{
	int j;
	double f,*xt;

	xt=dvector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(lens,xt);
	free_dvector(xt,1,ncom);
	return f;
}

#undef ITMAX
#undef ITMAX1
#undef NRANSI
#undef TOL
#undef GOLD
#undef GLIMIT
#undef TINY
#undef TINY1
#undef SHFT
#undef CGOLD
#undef ZEPS

