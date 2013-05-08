/*
 * tables.h
 *
 *  Created on: May 8, 2013
 *      Author: mpetkova
 */

#ifndef TABLES_H_
#define TABLES_H_

#include "standard.h"

double ffunction(double x);
double gfunction(double x);
double g2function(double x);
double mhat(double y, double beta);
double InterpolateFromTable(double *table, double *xtable, double y);

/// tables for the g,f,and g2 NFW functions
const long NTABLE = 1000;
const double maxrm = 100.0;

static double *mhattable,*ftable,*gtable,*g2table,*xtable;

void make_nfw_tables();
void make_pnfw_tables(double beta);

#endif /* TABLES_H_ */
