//
//  slicer.h
//  GLAMER
//
//  Created by Robert Benton Metcalf on 15/05/2020.
//

#ifndef slicer_h
#define slicer_h

template <typename V,typename P,typename R>
class Slicer{
public:
  Slicer(int Dim,double scale,int kmax = 10):
  w(scale),Kmax(kmax),dim(Dim){
  }

  void run(P &lnprob,std::vector<V> &chain,V xo,R &ran){
    int n = chain.size();
    int k = 0;
    if(lnprob(xo) <= -1.0e6){
      std::cerr << "Slicer: Initial point has 0 probability" << std::endl;
      throw std::runtime_error("bad point");
    }
    
    chain[0] = xo;
    for(int i=1;i<n;++i){
      chain[i] = chain[i-1];
      //step_double(chain[i],lnprob,k,ran);
      step_out(chain[i],lnprob,k,ran);
      k = (k+1)%dim;
    }
  }
  
  V xl,xr,ll,rr,lshrink,rshrink,x_new;
  double w;
  int Kmax;
  int dim;
  
  void step_double(V &x,P &lnprob,int i,R &ran){
    //double y = prob(x) * ran();
    double y = lnprob(x) + log(ran());
  
    // find range
    int k = Kmax;
    
    xl = x; xl[i] = x[i] - w * ran();;
    xr = x; xr[i] = xl[i] + w;
    
    while( k > 0 && ( y < lnprob(xl) || y < lnprob(xr) )){
      double u = ran();
      if(u<0.5){
        xl[i] = 2*xl[i] - xr[i];
      }else{
        xr[i] = 2*xr[i] - xl[i];
      }
      --k;
    }
    
    // shrink the range
    lshrink = xl;
    rshrink = xr;
    x_new = x;
    
    for(;;){
      x_new[i] = lshrink[i] + ran()*( rshrink[i] - lshrink[i] );
      if( y < lnprob(x_new) && accept(y,x,x_new,i,lnprob)) break;
      if(x_new[i] < x[i]){
        lshrink[i] = x_new[i];
      }else{
        rshrink[i] = x_new[i];
      }
    }
    
    swap(x,x_new);
  }

  void step_out(V &x,P &lnprob,int i,R &ran){
     //double y = prob(x) * ran();
    
    double y = lnprob(x) + log(ran());
   
     // find range
    xl = xr = x;
    
    xl[i] = x[i] - w * ran();;
    xr[i] = xl[i] + w;
    int j = floor(Kmax*ran());
    int k = Kmax - 1 - j;
    
    while( j > 0 && y < lnprob(xl)){
      xl[i] = xl[i] - w;
      --j;
    }
    while( k > 0 && y < lnprob(xr)){
      xr[i] = xr[i] + w;
      --k;
    }
     
    // shrink the range
    lshrink = xl;
    rshrink = xr;
    x_new = x;
     
    for(;;){
      x_new[i] = lshrink[i] + ran()*( rshrink[i] - lshrink[i] );
      if( y < lnprob(x_new) ) break;
      if(x_new[i] < x[i]){
        lshrink[i] = x_new[i];
      }else{
        rshrink[i] = x_new[i];
      }
    }
     
    swap(x,x_new);
   }

private:
  
  bool accept(double y,V &xo,V &x1,int i,P &lnprob){
    ll = xl;
    rr = xr;
    bool D = false;
    while(rr[i] - ll[i] > 1.1 * w){
      double m = (ll[i]+rr[i])/2;
      if( (xo[i] < m)*(x1[i] >= m) ||
          (xo[i] >= m)*(x1[i] < m) ) D = true;
      
      if(x_new[i] < m ){
        rr[i] = m;
      }else{
        ll[i] = m;
      }
      if(D && y >= lnprob(ll) && y >= lnprob(rr) ) return false;
    }
    return true;
  }
};


#endif /* slicer_h */
