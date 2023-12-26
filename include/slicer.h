//
//  slicer.h
//  GLAMER
//
//  Created by Robert Benton Metcalf on 15/05/2020.
//

#ifndef slicer_h
#define slicer_h

#include <cmath>
/**
 This is a slice sampler that can be used to do a markov chain without repeated entries in the chain.
 
 types;  V is usuatly a std::vector<float> or  std::vector<double>
      L is the functor type used for the likelihood function, this functure should
                have a opertor()(V p) that will return the log of the likeliwood
          
 */
template <typename V,typename L,typename R>
class Slicer{
public:
  Slicer(int Dim       /// number of parameters
         ,const V &scale  /// initial stepsize in parameter space
         ,int kmax = 10  ///  maximuma number of attempts made in each step
         ):
  w(scale),Kmax(kmax),dim(Dim),n_evals(0){
    assert(scale.size()>=Dim);
  }

 /// run the MC chain
  void run(L &lnprob               /// log likelihood function or funtor
           ,std::vector<V> &chain  /// will contain the MC chain on return. should be initialized to the desired length
           ,V xo                   /// initial point in parameter space
           ,R &ran                 /// rendom number generator
           ,int step_type /// 0 step_out, !=0 step_double
           ,bool verbose=false
           ){
    if(lnprob(xo) <= -1.0e20){
      std::cerr << "Slicer: Initial point has 0 probability" << std::endl;
      throw std::runtime_error("bad point");
    }
    n_evals=0;
    int n = chain.size();
    int k = 0;
    chain[0] = xo;
    if(verbose) std::cout << std::endl;
    for(int i=1;i<n;++i){
      chain[i] = chain[i-1];
      if(step_type) step_double(chain[i],lnprob,k,ran);
      else step_out(chain[i],lnprob,k,ran);
      k = (k+1)%dim;
      if(verbose) std::cout << i << "  " << n_evals << "  " << chain[i] << std::endl;
    }
  }
  
  /// returns the number of evaluations of the posterior during the last run
  size_t number_evaluations(){ return n_evals;}
  
private:
  V xl,xr,ll,rr,lshrink,rshrink,x_new;
  const V w;
  const int Kmax;
  const int dim;
  
  size_t n_evals;
  
  void step_double(V &x,L &lnprob,int i,R &ran){
    //double y = prob(x) * ran();
    double y = lnprob(x) + log(ran()); ++n_evals;
  
    assert(!std::isnan(y));
    // find range
    int k = Kmax;
    
    xl = x; xl[i] = x[i] - w[i] * ran();;
    xr = x; xr[i] = xl[i] + w[i];
    
    double lprob = lnprob(xl); ++n_evals;
    double rprob = lnprob(xr); ++n_evals;
 
    while( k > 0 && ( y < lprob || y < rprob )){
      double u = ran();
      if(u<0.5){
        xl[i] = 2*xl[i] - xr[i];
        lprob = lnprob(xl); ++n_evals;
      }else{
        xr[i] = 2*xr[i] - xl[i];
        rprob = lnprob(xr); ++n_evals;
      }
      --k;
    }
    
    // shrink the range
    lshrink = xl;
    rshrink = xr;
    x_new = x;
    
    for(;;){
      x_new[i] = lshrink[i] + ran()*( rshrink[i] - lshrink[i] );
      ++n_evals;
      if( y < lnprob(x_new) && accept(y,x,x_new,i,lnprob)) break;
      if(x_new[i] < x[i]){
        lshrink[i] = x_new[i];
      }else{
        rshrink[i] = x_new[i];
      }
    }
    
    swap(x,x_new);
  }

  void step_out(V &x,L &lnprob,int i,R &ran){
     //double y = prob(x) * ran();
    
    double y = lnprob(x) + log(ran()); ++n_evals;
   
     // find range
    xl = xr = x;
    
    xl[i] = x[i] - w[i] * ran();;
    xr[i] = xl[i] + w[i];
    int j = floor(Kmax*ran());
    int k = Kmax - 1 - j;
    
    while( j > 0 && y < lnprob(xl)){
      ++n_evals;
      xl[i] = xl[i] - w[i];
      --j;
    }
    while( k > 0 && y < lnprob(xr)){
      ++n_evals;
      xr[i] = xr[i] + w[i];
      --k;
    }
     
    // shrink the range
    lshrink = xl;
    rshrink = xr;
    x_new = x;
     
    for(;;){
      x_new[i] = lshrink[i] + ran()*( rshrink[i] - lshrink[i] );
      ++n_evals;
      if( y < lnprob(x_new) ) break;
      if(x_new[i] < x[i]){
        lshrink[i] = x_new[i];
      }else{
        rshrink[i] = x_new[i];
      }
    }
     
    swap(x,x_new);
   }
  
  bool accept(const double y,V &xo,V &x1,int i,L &lnprob){
    ll = xl;
    rr = xr;
    bool D = false;
    
    double lprob = lnprob(ll); ++n_evals;
    double rprob = lnprob(rr); ++n_evals;

    while( (rr[i] - ll[i]) > 1.1 * w[i]){
      double m = (ll[i]+rr[i])/2;
      if( (xo[i] < m)*(x1[i] >= m) ||
          (xo[i] >= m)*(x1[i] < m) ) D = true;
      
      if(x_new[i] < m ){
        rr[i] = m;
        rprob = lnprob(rr); ++n_evals;
      }else{
        ll[i] = m;
        lprob = lnprob(ll); ++n_evals;
      }
      if(D && y >= lprob && y >= rprob ) return false;
    }
    return true;
  }
};

template <typename V,typename L,typename R>
class SimpleMH_MCMC{
private:
 
  const V w;
  int dim;
  
  size_t n_evals;
  std::vector<int> active;
  
public:
  SimpleMH_MCMC(int Dim       /// number of parameters
         ,const V &scale  /// initial stepsize in parameter space
         ,std::vector<int> active_parameters  = {} /// active variable, empty does first ones
         ):
  w(scale),dim(Dim),n_evals(0),active(active_parameters)
  {
    
    if(active.size() == 0){
      for(int i ; i<dim ; ++i) active.push_back(i);
    }else{
      if(dim > active.size()) dim = active.size();
    }
    
    assert(scale.size() >= dim);
    assert(active.size() >= dim);
  }
  
  /// run the MC chain
   void run(L &lnprob               /// log likelihood function or funtor
            ,std::vector<V> &chain  /// will contain the MC chain on return. should be initialized to the desired length
            ,V xo                   /// initial point in parameter space
            ,R &ran                 /// rendom number generator
            ,bool cyclic = true     /// if true it cycles thruogh the parameters one at a time
            ,bool verbose=false
            ){
     
     double lnp = lnprob(xo);
     n_evals=0;
     int n = chain.size();
     chain[0] = xo;
     if(verbose) std::cout << std::endl;
     
     int k=0;
     for(int i=1;i<n;++i){
       chain[i] = chain[i-1];
     
       if(cyclic){
         int j = active[k];
         chain[i][j] += w[j]*ran.gauss();
         k = (k+1) % dim;
       }else{
         for(int j=0;j<dim;++j){
           k = active[k];
           chain[i][j] += w[j]*ran.gauss();
         }
       }
       double lnp2 = lnprob(chain[i]);
     
       if( lnp2 <= lnp && exp(lnp2 - lnp) < ran() ){
         chain[i] = chain[i-1];
       }else{
         lnp = lnp2;
         ++n_evals;
       }
       if(verbose) std::cout << i << "  " << n_evals << "  " << chain[i] << std::endl;
     }
   }
  
  /// run the MC chain
  /// class L must have method aux_update(R )
   void runAux(L &lnprob               /// log likelihood function or funtor
            ,std::vector<V> &chain  /// will contain the MC chain on return. should be initialized to the desired length
            ,V xo                   /// initial point in parameter space
            ,R &ran                 /// rendom number generator
            ,bool cyclic = true     /// if true it cycles thruogh the parameters one at a time
            ,bool verbose=false
            ){
     
     double lnp = lnprob(xo);
     n_evals=0;
     int n = chain.size();
     chain[0] = xo;
     lnprob.aux_update(ran,lnp);
     
     if(verbose) std::cout << std::endl;
     
     int k=0;
     for(int i=1;i<n;++i){
       chain[i] = chain[i-1];
       
       if(i%dim == 0){
         lnprob.aux_update(ran,lnp);
       }else{
         if(cyclic){
           chain[i][k] += w[k]*ran.gauss();
           k = (k+1) % dim;
         }else{
           for(int j=0;j<dim;++j){
             chain[i][j] += w[j]*ran.gauss();
           }
         }
         double lnp2 = lnprob(chain[i]);
         
         if( lnp2 <= lnp && exp(lnp2 - lnp) < ran() ){
           chain[i] = chain[i-1];
         }else{
           lnp = lnp2;
           ++n_evals;
         }
       }
       if(verbose) std::cout << i << "  " << n_evals << "  " << chain[i] << std::endl;
     }
   }

  /// returns the number of evaluations of the posterior during the last run
  size_t number_evaluations(){ return n_evals;}
  
};
#endif /* slicer_h */
