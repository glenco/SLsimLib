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

    void setInternal();
    void change_redshifts(TreeHndl i_tree,TreeHndl s_tree,double z_source,double z_lens);
    void RandomizeModel(double r_source_physical,long *seed,bool tables);
};

typedef Model *ModelHndl;

#endif /* MODEL_H_ */
