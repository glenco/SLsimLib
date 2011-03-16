/*
 * create_elliptical.c
 *
 *  Created on: Feb 28, 2010
 *      Author: R.B. Metcalf
 */
#include <math.h>
#include <nrD.h>
#include "galaxies.h"

void create_sersic(int n,double Ro,double f,double *center,double theta,double **x,long Nsources){
	static long seed=19282;
	double r,a;
	int i;

	for(i=0;i<Nsources;++i){
		if(n==1){ // exponental
			do r=ran2(&seed); while(r==0);
			r=-log(r)*Ro;
		}else{ ERROR_MESSAGE(); printf("ERROR; create_sersic, only does exponential right now\n"); exit(1);}

		a=2*3.141593*ran2(&seed);
		x[i][0]=r*(cos(a)*cos(theta)-f*sin(a)*sin(theta))+center[0];
		x[i][1]=r*(cos(a)*sin(theta)+f*sin(a)*cos(theta))+center[1];
	}
}
