/*
 * internal_rayshooter_multi.cpp
 *
 *  Created on: Jan 19, 2012
 *      Author: mpetkova, bmetcalf
 */

#include <slsimlib.h>

/** \ingroup DeflectionL2
 *
 * \brief This is the function that does the deflection calculation with multiple lens planes.
 *
 * The convergence, shear and rotation will be calculated if kappa_off == false
 *
 * Needs to be calculated before calling:
 *
 *
 * Dl[j = 0...Nplane] - The angular size distance between the observer and the jth plane not counting the observer plane.
 * 	                    Dl[0] is the first plane with mass on it and Dl[Nplane] is the distance to the source plane.
 * dDl[j = 0...Nplane] - The angular size distance between the jth and (j+1)th planes counting the observer plane as j=0.
 *                      dDl[0] = Dl[0], dDl[Nplane] is between the last plane with mass and the source plane.
 *
 * charge = 4*G*mass_scale/c^2 in units of Mpc
 *
 * i_points[].x[] is in angular units.
 *
 * Warning: Is not valid for a non-flat universe.
 */

void MultiLens::rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off){
	unsigned long i;
	int j;
	double kappa,a,b,c;
    double alpha[2], gamma[2];
    double xminus[2],xplus[2];
    double kappa_minus,gamma_minus[3],kappa_plus,gamma_plus[3];

    Nplanes = getNplanes();

    for(i = 0; i< Npoints; i++){
    	i_points[i].image->x[0] = i_points[i].x[0];
    	i_points[i].image->x[1] = i_points[i].x[0];
    	i_points[i].kappa = 0.0;
    	i_points[i].gamma[0] = 0.0;
    	i_points[i].gamma[1] = 0.0;
    }


    for(i = 0; i< Npoints; i++){
     }

	for(i = 0; i< Npoints; i++){

		// find position on first lens plane
		i_points[i].image->x[0] = i_points[i].x[0]*Dl[0];
       	i_points[i].image->x[1] = i_points[i].x[1]*Dl[0];

       	xminus[0] = 0;
		xminus[1] = 0;

		// Set magnification matrix on first plane.  Also the default if kappa_off == false
		kappa_minus = 0;
		gamma_minus[0] = 0;
		gamma_minus[1] = 0;
		gamma_minus[2] = 0;

		i_points[i].kappa = 1;  // This is actually 1-kappa until after the loop through the planes.
		i_points[i].gamma[0] = 0;
		i_points[i].gamma[1] = 0;
		i_points[i].gamma[2] = 0;

		for(j = 1; j < Nplanes+1; j++){

    		halo_tree[j]->force2D(i_points[i].image->x,&alpha[0],&kappa,&gamma[0],false);

       		xplus[0] = i_points[i].image->x[0]*(dDl[j]+dDl[j-1])/dDl[j-1]
        				                 - xminus[0]*dDl[j]/dDl[j-1]
        		                         - charge*dDl[j]*alpha[0];
       		xplus[1] = i_points[i].image->x[1]*(dDl[j]+dDl[j-1])/dDl[j-1]
        				                 - xminus[1]*dDl[j]/dDl[j-1]
        		                         - charge*dDl[j]*alpha[1];

			xminus[0] = i_points[i].image->x[0];
			xminus[1] = i_points[i].image->x[1];

			i_points[i].image->x[0] = xplus[0];
			i_points[i].image->x[1] = xplus[1];

			if(!kappa_off)
    		{

				a = (1+dDl[j]+dDl[j-1])*Dl[j-1]/dDl[j-1]/Dl[j];
				b = dDl[j]*Dl[j-2]/dDl[j-1]/Dl[j];
				c = charge*dDl[j]*Dl[j-1]/Dl[j];

				kappa_plus = a*i_points[i].kappa - b*kappa_minus
						- c*(kappa*i_points[i].kappa + gamma[0]*i_points[i].gamma[0] + gamma[1]*i_points[i].gamma[1]);
				gamma_plus[0] = a*i_points[i].gamma[0] - b*gamma_minus[0]
						- c*(kappa*i_points[i].gamma[0] + gamma[0]*i_points[i].kappa - gamma[1]*i_points[i].gamma[3]);
				gamma_plus[1] = a*i_points[i].gamma[1] - b*gamma_minus[1]
						- c*(kappa*i_points[i].gamma[1] + gamma[1]*i_points[i].kappa - gamma[0]*i_points[i].gamma[3]);
				gamma_plus[2] = a*i_points[i].gamma[2] - b*gamma_minus[2]
						- c*(kappa*i_points[i].gamma[2] + gamma[0]*i_points[i].gamma[1] - gamma[1]*i_points[i].gamma[0]);

				kappa_minus = i_points[i].kappa;
				gamma_minus[0] = i_points[i].gamma[0];
				gamma_minus[1] = i_points[i].gamma[1];
				gamma_minus[2] = i_points[i].gamma[2];

				i_points[i].kappa = kappa_plus;
				i_points[i].gamma[0] = gamma_plus[0];
				i_points[i].gamma[1] = gamma_plus[1];
				i_points[i].gamma[2] = gamma_plus[2];
    		}

		}
		// Convert units back to angles.

		i_points[i].image->x[0] /= Dl[Nplanes];
		i_points[i].image->x[1] /= Dl[Nplanes];

		i_points[i].kappa = 1 - i_points[i].kappa;
    }

    return;
}
