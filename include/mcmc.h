/*
 * mcmc.h
 *
 *  Created on: Nov 10, 2011
 *      Author: mpetkova
 */

#ifndef MCMC_H_
#define MCMC_H_

typedef struct{
  double x, value, error;
} data;

typedef struct{
  double u, v;
  double a, b;
  double lnL, lnWt;

} obj;

int mcmc(void);
void setPrior(obj*);
double getlnLhood(double, double);
void Explore(obj*, double);
void Results(obj*, int, double);



#endif /* MCMC_H_ */
