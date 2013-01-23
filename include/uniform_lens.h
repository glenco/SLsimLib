/*
 * uniform_lens.h
 *
 *  Created on: Jan 21, 2013
 *      Author: D. Leier
 */

#ifndef UNIFORM_LENS_H_
#define UNIFORM_LENS_H_

#include <lens.h>
#include <Tree.h>
#include <forceTree.h>
#include <quadTree.h>
#include <source.h>
#include <base_analens.h>

/**
 * \brief This is the uniform lens model
 */
class UniLens: public BaseAnaLens{
public:
	UniLens(InputParams& params);
	~UniLens();

	virtual void assignParams(InputParams& params);
	void PrintLens(bool show_substruct,bool show_stars);
/*	void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off){
		for(unsigned long i=0;i<Npoints;++i){
			i_points[i].image->kappa = i_points[i].kappa = kappa;
			i_points[i].image->gamma[0] = i_points[i].gamma[0] = gamma[0];
			i_points[i].image->gamma[1] = i_points[i].gamma[1] = gamma[1];
			//i_points[i].image->gamma[2] = i_points[i].gamma[2] = 0.0;
			i_points[i].image->x[0] = (1-kappa-gamma[0])*i_points[i].x[0]-gamma[1]*i_points[i].x[1];
			i_points[i].image->x[1] = (1-kappa+gamma[0])*i_points[i].x[1]-gamma[1]*i_points[i].x[0];
		}
	}
	void rayshooterInternal(double *ray, double *alpha, float *gamma, float *kappa, bool kappa_off){
		alpha[0] = alpha[1] = 0.0;
		gamma[0] = gamma[1] = 0.0;
		*kappa = 0.0;
	}*/


private:
   float kappa_uniform;
   float gamma_uniform[2];

};


#endif /* UNIFORM_LENS_H_ */
