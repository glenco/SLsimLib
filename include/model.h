/*
 * model.h
 *
 *  Created on: Jan 23, 2012
 *      Author: mpetkova
 */

#include <analytic_lens.h>
#include <multiplane.h>
#include <source.h>

#ifndef MODEL_H_
#define MODEL_H_

class Model{
public:
	LensHndl lens;
	SourceHndl source;
	CosmoHndl cosmo;

	//double Dl, Ds, Dls;

	Model(LensHndl,SourceHndl,CosmoHndl);
	~Model();

	double getZsource(){return source->zsource;}
	double getZlens(){return lens->getZlens();}

    void RandomizeModel(double r_source_physical,long *seed,bool tables, double angle_factor=1.0);
private:
    void setInternal();
    void change_redshifts(TreeHndl i_tree,TreeHndl s_tree,double z_source,double z_lens);
};

typedef Model *ModelHndl;

#endif /* MODEL_H_ */
