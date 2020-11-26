#include <math.h>
#include <nrutil.h>
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

//shift in pointer that goes into functor ?

template <typename T,typename F>
T f1dim_tp(T x,std::vector<T> &pcom,std::vector<T> &xicom,F &functor)
{
  int j;
  T f;
  int ncom = pcom.size()-1;
  std::vector<T> xt(ncom+1);
  for (j=1;j<=ncom;j++) xt[j] = pcom[j] + x * xicom[j];
  f = functor(xt.data() + 1);
  
  return f;
}

template <typename T,typename F>
T brent_tp(T ax, T bx, T cx, T tol,T *xmin,std::vector<T> &pcom,std::vector<T> &xicom,F &functor)
{
  int iter;
  T a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  T e=0.0;

  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=f1dim_tp(x,pcom,xicom,functor);
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
    fu=f1dim_tp(u,pcom,xicom,functor);
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
  throw std::runtime_error("Too many iterations in brent_tpd");
  *xmin=x;
  return fx;
}


template <typename T,typename F>
void mnbrak_tp(T *ax, T *bx, T *cx, T *fa, T *fb, T *fc,std::vector<T> &pcom,std::vector<T> &xicom
                     ,F &functor)
{
  T ulim,u,r,q,fu,dum;

  *fa=f1dim_tp(*ax,pcom,xicom,functor);
  *fb=f1dim_tp(*bx,pcom,xicom,functor);
  if (*fb > *fa) {
    SHFT(dum,*ax,*bx,dum)
    SHFT(dum,*fb,*fa,dum)
  }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc=f1dim_tp(*cx,pcom,xicom,functor);
  while (*fb > *fc) {
    r=(*bx-*ax)*(*fb-*fc);
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
      (2.0*SIGN(DMAX(fabs(q-r),TINY1),q-r));
    ulim=(*bx)+GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0) {
      fu=f1dim_tp(u,pcom,xicom,functor);
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
      fu=f1dim_tp(u,pcom,xicom,functor);
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu=f1dim_tp(u,pcom,xicom,functor);
      if (fu < *fc) {
        SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
        SHFT(*fb,*fc,fu,f1dim_tp(u,pcom,xicom,functor))
      }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u=ulim;
      fu=f1dim_tp(u,pcom,xicom,functor);
    } else {
      u=(*cx)+GOLD*(*cx-*bx);
      fu=f1dim_tp(u,pcom,xicom,functor);
    }
    SHFT(*ax,*bx,*cx,u)
    SHFT(*fa,*fb,*fc,fu)
  }
}

template <typename T,typename F>
void linmind_tp(T p[], T xi[], int n, T *fret, F &functor)
{

  int j;
  T xx,xmin,fx,fb,fa,bx,ax;

  std::vector<T> pcom(n+1);
  std::vector<T> xicom(n+1);
  //nrfunc=func;
  for (j=1;j<=n;j++) {
    pcom[j]=p[j];
    xicom[j]=xi[j];
  }
  ax=0.0;
  xx=1.0;
  mnbrak_tp(&ax,&xx,&bx,&fa,&fx,&fb,pcom,xicom,functor);
  *fret=brent_tp(ax,xx,bx,TOL,&xmin,pcom,xicom,functor);
  for (j=1;j<=n;j++) {
    xi[j] *= xmin;
    p[j] += xi[j];
  }
}

/** \brief Minimization routine that does not require derivatives.
 
 This is evolved from the NR routine.  It differs in that it is C++ complient, thread safe and uses 0 indexing, ie pp and
 the input to functor starts at [0] not [1].
 
 */

template <typename T,typename F>
void powell_tp(
                     T *pp         /// initially a gues for x and final solution on completion
                     , int n       /// number of dimenstions
                     , T ftol      /// target tollarence in f(x)
                     , int *iter   ///
                     , T *fret     ///
                     ,F &functor   ///  function to be minimized.  Must have a method T F::opertator()(T *)
                     )
{
	int i,ibig,j;
	T del,fp,fptt,t;
  T *p = pp-1;
  
  std::vector<T> simplex((n+1)*(n+1));
  
  T **xi = new T*[n+1];
  xi[0] = simplex.data();
  for (long i = 1; i <= n; ++i) xi[i] = xi[0] + i * (n+1);
  
  for (long i = 1; i <= n; ++i){
    xi[i][i] = 1;
    for (long j = i+1; j <= n; ++j) {
      xi[i][j] = xi[j][i] = 0;
    }
  }
  
	std::vector<T> pt(n+1);
	std::vector<T> ptt(n+1);
	std::vector<T> xit(n+1);
	*fret = functor(p+1);
	for (j=1;j<=n;j++) pt[j]=p[j];
	for (*iter=1;;++(*iter)) {
		fp=(*fret);
		ibig=0;
		del=0.0;
		for (i=1;i<=n;i++) {
			for (j=1;j<=n;j++) xit[j]=xi[j][i];
			fptt=(*fret);
      linmind_tp(p,xit.data(),n,fret,functor);
			if (fptt-(*fret) > del) {
				del=fptt-(*fret);
				ibig=i;
			}
		}
		if (2.0*(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))+TINY) {
			return;
		}
    if (*iter == ITMAX){
      throw std::runtime_error("powell_tp exceeding maximum iterations.");
    };
		for (j=1;j<=n;j++) {
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
    fptt=functor(ptt.data()+1);
		if (fptt < fp) {
			t=2.0*(fp-2.0*(*fret)+fptt)*sqrt(fp-(*fret)-del)-del*sqrt(fp-fptt);
			if (t < 0.0) {
        linmind_tp(p,xit.data(),n,fret,functor);
				for (j=1;j<=n;j++) {
					xi[j][ibig]=xi[j][n];
					xi[j][n]=xit[j];
				}
			}
		}
	}
  delete [] xi;
}

/* note #undef's at end of file */


/* note #undef's at end of file */


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

