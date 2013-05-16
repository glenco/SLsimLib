/*
 * singlelens.h
 *
 *  Created on: May 10, 2013
 *      Author: mpetkova
 */

#ifndef SINGLELENS_H_
#define SINGLELENS_H_

#include "standard.h"
#include "InputParams.h"
#include "lens.h"
#include "analytic_lens.h"
#include "uniform_lens.h"


class SingleLens : public Lens{
public:
	SingleLens(InputParams& params);
	~SingleLens();

	double getZlens();
	void setZlens(CosmoHndl cosmo,double zlens,double zsource = 1000);
	void setInternalParams(CosmoHndl,SourceHndl);
	void setInternalParams(CosmoHndl,double);
	void error_message1(std::string name,std::string filename);
	void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off);
	void rayshooterInternal(double *ray, double *alpha, KappaType *gamma, KappaType *kappa, bool kappa_off);

	IntProfType DM_halo_type;
	IntProfType galaxy_halo_type;
	int Nprof;
	std::vector<LensHalo *> halo;

private:

	void assignParams(InputParams& params);
};


#endif /* SINGLELENS_H_ */
