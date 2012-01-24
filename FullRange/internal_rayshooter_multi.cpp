/*
 * internal_rayshooter_multi.cpp
 *
 *  Created on: Jan 19, 2012
 *      Author: mpetkova
 */

#include <slsimlib.h>

void multiLens::rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off){
	int i, j;
	double convert_fac, tmp, Nplanes;
    double alpha[2], gamma[2];

    Nplanes = getNplanes();

    for(i = 0; i< Npoints; i++){
    	i_points[i].image->x[0] = i_points[i].x[0];
    	i_points[i].image->x[1] = i_points[i].x[0];
    	i_points[i].kappa = 0.0;
    	i_points[i].gamma[0] = 0.0;
    	i_points[i].gamma[1] = 0.0;
    }

    for(j = 0; j < Nplanes; j++){

		convert_fac = mass_scale / Sigma_crit[j] * Dl[0] / Dl[j];

		for(i = 0; i< Npoints; i++){

    		halo_tree[j]->force2D(i_points[i].x,&alpha[0],&tmp,&gamma[0],false);

    		i_points[i].image->x[0] += convert_fac*alpha[0];
    		i_points[i].image->x[1] += convert_fac*alpha[1];

    		if(!kappa_off)
    		{
    			i_points[i].kappa += convert_fac*tmp;
    			i_points[i].gamma[0] += convert_fac*gamma[0];
    			i_points[i].gamma[1] += convert_fac*gamma[1];
    		}

    	}
    }

}
