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

	ftable = new double[NTABLE];
	gtable = new double[NTABLE];
	g2table = new double[NTABLE];

	for(i = 0 ; i< NTABLE; i++){
		x = (i+1)*dx;
		ftable[i] = ffunction(x);
		gtable[i] = gfunction(x);
		g2table[i] = g2function(x);
	}
}


void QuadTreePseudoNFW::make_tables(){
	int i;
	double x, dx = maxrm/(double)NTABLE;

	mhattable = new double[NTABLE];

	for(i = 0 ; i< NTABLE; i++){
		x = (i+1)*dx;
		mhattable[i] = mhat(x,beta);
	}
}

void ForceTreeNFW::make_tables(){

	int i;
	double x, dx = maxrm/(double)NTABLE;

	ftable = new double[NTABLE];
	gtable = new double[NTABLE];
	g2table = new double[NTABLE];

	for(i = 0 ; i< NTABLE; i++){
		x = (i+1)*dx;
		ftable[i] = ffunction(x);
		gtable[i] = gfunction(x);
		g2table[i] = g2function(x);
	}
}

void ForceTreePseudoNFW::make_tables(){
	int i;
	double x, dx = maxrm/(double)NTABLE;

	mhattable = new double[NTABLE];

	for(i = 0 ; i< NTABLE; i++){
		x = (i+1)*dx;
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
	if(beta == 1.0) return y - log(1+y);
	if(beta == 2.0) return log(1+y) - y/(1+y);
	if(beta>=3.0) return ( (1 - beta)*y + pow(1+y,beta-1) - 1)/(beta-2)/(beta-1)/pow(1+y,beta-1);

	ERROR_MESSAGE();
	std::cout << "Only beta ==1, ==2 and >=3 are valid" << std::endl;
	exit(1);
	return 0.0;
}

double InterpolateFromTable(double *table, double y){
	int j;
	j=(int)(y/maxrm*NTABLE);
	return (table[j+1]-table[j])*(y-j*maxrm/float(NTABLE)) + table[j];
}

