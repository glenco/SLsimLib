#include <iostream>
#include <math.h>

using namespace std;

const int N_data = 64;
const int N_obj = 100;
const float A = 1.34;
const float B = 0.87;
const int MAX = 10000;

#define UNIFORM (rand() + 0.5)/(RAND_MAX + 1.0) 
#define PLUS(x,y) (x>y ? x+log(1+exp(y-x)) : y+log(1+exp(x-y)))

class my_data{
public:
  double x[N_data], value[N_data], error[N_data];
  my_data();
};

class my_obj{
public:
  double u, v;
  double a, b;
  double lnL, lnWt;
  void setPrior();
  void getlnLhood(double, double, my_data);
  my_obj Explore(double, my_data);
};

void Results(my_obj*, int, double);

my_data::my_data(){
  int i;
  
  for(i = 0; i < N_data; i++){
    x[i] = UNIFORM;
    value[i] = (A * x[i] + B);
    value[i] += (UNIFORM - 0.5) * 0.1 * value[i];
    error[i] = 0.1 * value[i];
  }  
}

void my_obj::getlnLhood(double a, double b, my_data D){
  int i;
  double value, residue, L;
  
  for(i = 0, L = 0.0; i < N_data; i++){
    value = a * D.x[i] + b;
    residue = (value - D.value[i]) / D.error[i];
    L += - 0.5 * residue * residue - log(D.error[i] * sqrt(2.0 * M_PI));
  }

  lnL = L;
}

void my_obj::setPrior(){
  u = UNIFORM;
  v = UNIFORM;
  a = u * (2.0 - 1.0) + 1.0;
  b = v;
}

my_obj my_obj::Explore(double lnLstar, my_data D){
  double  step = 0.1;  
  int     m    = 20;   
  int     accept = 0;  
  int     reject = 0;  
  my_obj  Try, Get;

  for( ; m > 0; m-- ){
	  // Trial object
	  Try.u = u + step * (2.*UNIFORM - 1.);  // |move| < step
	  Try.v = v + step * (2.*UNIFORM - 1.);  // |move| < step
	  Try.u -= floor(Try.u);      // wraparound to stay within (0,1)
      Try.v -= floor(Try.v);      // wraparound to stay within (0,1)                                                    
      Try.a = Try.u * (2.0 - 1.0) + 1.0;  // map to x                                                                           
      Try.b = Try.v;        // map to y                                                                           
      Try.getlnLhood(Try.a, Try.b, D);  // trial likelihood value                                                     
      // Accept if and only if within hard likelihood constraint                                                                
      if(Try.lnL > lnLstar){  
    	  Get = Try;
    	  accept++;
      }
      else
	reject++;
      // Refine step-size to let acceptance ratio converge around 50%                                                           
      if(accept > reject)    step *= exp(1.0 / accept);
      if(accept < reject)    step /= exp(1.0 / reject);
    }

 return Get;
}

void Results(my_obj* sample, int nest, double lnZ){
  double a = 0.0, aa = 0.0;   // 1st and 2nd moments of x                                                               
  double b = 0.0, bb = 0.0;   // 1st and 2nd moments of y                                                               
  double w;                   // Proportional weight                                                                    
  int    i;                   // Sample counter                                                                         

  for(i = 0; i < nest; i++){
    w = exp(sample[i].lnWt - lnZ);
    a  += w * sample[i].a;
    aa += w * sample[i].a * sample[i].a;
    b  += w * sample[i].b;
    bb += w * sample[i].b * sample[i].b;
  }
  
  cout<<"real A = "<< A <<" mean(A) = "<< a <<" stddev(A) = "<< sqrt(aa - a*a) <<"\n";
  cout<<"real B = "<< B <<" mean(B) = "<< b <<" stddev(B) = "<< sqrt(bb - b*b) <<"\n";
}

int mcmc(){
  double lnwidth;      // ln(width in prior mass)                                                                      
  double lnLstar;      // ln(Likelihood constraint)                                                                    
  double H    = 0.0;    // Information, initially 0                                                                     
  double lnZ = 0.0;// ln(Evidence Z, initially 0)                                                                  
  double lnZnew;       // Updated logZ                                                                                 
  int    i;             // Object counter                                                                               
  int    copy;          // Duplicated object                                                                            
  int    worst;         // Worst object                                                                                 
  int    nest;          // Nested sampling iteration count  
  my_obj Obj[N_obj];
  my_obj Sample[MAX];
  my_data D;
                                                                              
  for(i = 0; i < N_obj; i++){
    Obj[i].setPrior();
    Obj[i].getlnLhood(Obj[i].a,Obj[i].b,D);
   } 
                                                                                      
  lnwidth = log(1.0 - exp(-1.0 / N_obj));
                                        
  for(nest = 0; nest < MAX; nest++){

    for(i = 1, worst = 0; i < N_obj; i++)
      if(Obj[i].lnL < Obj[worst].lnL)  
	worst = i;
    
    Obj[worst].lnWt = lnwidth + Obj[worst].lnL;
                                                                                  
    lnZnew = PLUS(lnZ, Obj[worst].lnWt);

    H = exp(Obj[worst].lnWt - lnZnew) * Obj[worst].lnL
      + exp(lnZ - lnZnew) * (H + lnZ) - lnZnew;

    lnZ = lnZnew;
                                                                                         
    Sample[nest] = Obj[worst];
                                                              
    do{
      copy = (int)(N_obj * UNIFORM) % N_obj;                                                         
    }while(copy == worst && N_obj > 1);                                                 
    
    lnLstar = Obj[worst].lnL;                                                    
    Obj[worst] = Obj[copy];                                                       
                                                                                
    Obj[worst] = Obj[worst].Explore(lnLstar,D);
                                                                                                     
    lnwidth -= 1.0 / N_obj;
  }                                            
                                                   
  cout<<"# iterates = "<< nest <<"\n";
  cout<<"Evidence: ln(Z) = "<< lnZ <<" +- "<< sqrt(H/N_obj) <<"\n";
  cout<<"Information: H = "<< H <<" nats = "<< H/log(2.) <<" bits\n";
  Results(&Sample[0], nest, lnZ);

  return 0;
} 

