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

	MultiLens *multilens;
	AnaLens *analens;

	Model(std::string paramfile, long *seed, bool multi_lens);
	~Model();

	double getZsource(){return source->zsource;}
	double getZlens(){return lens->getZlens();}

    void RandomizeModel(double r_source_physical,long *seed,bool tables, double angle_factor=1.0);

private:

	/// names of clump and sb models
	typedef enum {Uniform,Gaussian,BLR_Disk,BLR_Sph1,BLR_Sph2} SourceSBModel;

    void setInternal();
    void change_redshifts(TreeHndl i_tree,TreeHndl s_tree,double z_source,double z_lens);
    void readParamfile(std::string paramfile);
	SourceSBModel sb_type;
};

typedef Model *ModelHndl;

#endif /* MODEL_H_ */
