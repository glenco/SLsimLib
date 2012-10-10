/*
 * model.h
 *
 *  Created on: Jan 23, 2012
 *      Author: mpetkova
 */

#ifndef MODEL_H_
#define MODEL_H_

#include <analytic_lens.h>
#include <multiplane.h>
#include <source.h>
/**
 * \brief A class to group the lens, source and cosmology models into one package and allows for additional
 * initialization to occur that requires a combination of the above structures.
 */
class Model{
public:

	LensHndl lens;
	SourceHndl source;
	CosmoHndl cosmo;

	Model(LensHndl mylens  /// lens
			,SourceHndl mysource /// source
			,CosmoHndl mycosmo);
	~Model();

	double getZsource(){return source->getZ();}
	double getZlens(){return lens->getZlens();}

    void RandomizeModel(double r_source_physical,long *seed,bool tables,bool randomize_host_z=true,bool randomize_source_z=true,bool in_radians=false);

private:

    void setInternal();
    void change_redshifts(TreeHndl i_tree,TreeHndl s_tree,double z_source,double z_lens);
};

typedef Model *ModelHndl;

#endif /* MODEL_H_ */
