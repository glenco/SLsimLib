/*
 * model.cpp
 *
 *  Created on: Jan 23, 2012
 *      Author: mpetkova
 */

#include <model.h>
#include <source_models.h>
/// Sets some distances and passes cosmo and source into the initialization of lens.
Model::Model(LensHndl my_lens, SourceHndl my_source, CosmoHndl my_cosmo){
	lens = my_lens;
	source = my_source;
	cosmo = my_cosmo;

	setInternal();
}

Model::~Model(){
	delete lens;
	delete source;
	delete cosmo;
}

void Model::setInternal(){
	lens->setInternalParams(cosmo,source);
	source->setDlDs(cosmo->angDist(0,lens->getZlens()) / cosmo->angDist(0,source->getZ()));
}

