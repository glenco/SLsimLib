/*
 * profiles.cpp
 *
 *  Created on: Mar 6, 2012
 *      Author: bmetcalf
 */
#include <slsimlib.h>

/// tables for the g,f,and g2 NFW functions
//double *ftable, *gtable, *g2table, *mhattable;
const long NTABLE = 1000;
const double maxrm = 100.0;

void QuadTreeNFW::make_tables(){
	int i;
	double x, dx = maxrm/(double)NTABLE;

	xtable = new double[NTABLE];
	ftable = new double[NTABLE];
	gtable = new double[NTABLE];
	g2table = new double[NTABLE];

	for(i = 0 ; i< NTABLE; i++){
	        x = i*dx;
		xtable[i] = x;
		ftable[i] = ffunction(x);
		gtable[i] = gfunction(x);
		g2table[i] = g2function(x);
	}
}


void QuadTreePseudoNFW::make_tables(){
	int i;
	double x, dx = maxrm/(double)NTABLE;

	xtable = new double[NTABLE];
	mhattable = new double[NTABLE];

	for(i = 0 ; i< NTABLE; i++){
		x = i*dx;
		xtable[i] = x;
		mhattable[i] = mhat(x,beta);
	}
}

void ForceTreeNFW::make_tables(){

	int i;
	double x, dx = maxrm/(double)NTABLE;

	xtable = new double[NTABLE];
	ftable = new double[NTABLE];
	gtable = new double[NTABLE];
	g2table = new double[NTABLE];

	for(i = 0 ; i< NTABLE; i++){
		x = i*dx;
		xtable[i] = x;
		ftable[i] = ffunction(x);
		gtable[i] = gfunction(x);
		g2table[i] = g2function(x);
	}
}

void ForceTreePseudoNFW::make_tables(){
	int i;
	double x, dx = maxrm/(double)NTABLE;

	xtable = new double[NTABLE];
	mhattable = new double[NTABLE];

	for(i = 0 ; i< NTABLE; i++){
		x = i*dx;
		xtable[i] = x;
		mhattable[i] = mhat(x,beta);
	}
}

/*void QuadTreeNFW::point_tables(){
	gt = gtable;
	ft = ftable;
	g2t= g2table;
}*/

/*void QuadTreePseudoNFW::point_tables(){
	mhatt = mhattable;
}*/

double rhos(double x){
	double ans;

	ans = 1.0;
	ans /= log(1+x) - x/(1+x);

	return ans;
}
/// Auxiliary function for PseudoNFW profile
double mhat(double y, double beta){
  if(y==0) y=1e-5;
	if(beta == 1.0) return y - log(1+y);
	if(beta == 2.0) return log(1+y) - y/(1+y);
	if(beta>=3.0) return ( (1 - beta)*y + pow(1+y,beta-1) - 1)/(beta-2)/(beta-1)/pow(1+y,beta-1);

	ERROR_MESSAGE();
	std::cout << "Only beta ==1, ==2 and >=3 are valid" << std::endl;
	exit(1);
	return 0.0;
}

double InterpolateFromTable(double *table, double *x, double y){
	int j;
	j=(int)(y/maxrm*NTABLE);
	assert(y>=x[j] && y<=x[j+1]);
	return (table[j+1]-table[j])/(x[j+1]-x[j])*(y-x[j]) + table[j];
}

