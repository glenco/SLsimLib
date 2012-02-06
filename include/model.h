/*
 * model.h
 *
 *  Created on: Jan 23, 2012
 *      Author: mpetkova
 */

#ifndef MODEL_H_
#define MODEL_H_

#include <Tree.h>
#include <analytic_lens.h>
#include <multiplane.h>

class Model{
public:
	LensHndl lens;
	SourceHndl source;
	CosmoHndl cosmo;

	Model(char *filename, char *filename2, int flag);
	~Model();

	// in internal_rayshooter_nfw.c
	double uniform_SB(double *y);
	double gaussian_SB(double *y);
	double BLR_Disk_SB(double *y);
	double BLR_Sph1_SB(double *y);
	double BLR_Sph2_SB(double *y);

    void setInternal();
};

/// pointer to surface brightness function
static double (Model::*source_sb_func)(double *y);

typedef Model *ModelHndl;

void change_redshifts(TreeHndl i_tree,TreeHndl s_tree,ModelHndl model,double z_source
		,double z_lens);

#endif /* MODEL_H_ */
