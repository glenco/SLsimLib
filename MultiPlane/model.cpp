/*
 * model.cpp
 *
 *  Created on: Jan 23, 2012
 *      Author: mpetkova
 */

#include <model.h>
#include <source_models.h>

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
	lens->setInternalParams(cosmo,source->zsource);

	double zlens = lens->getZlens();
	source->DlDs = cosmo->angDist(0,zlens) / cosmo->angDist(0,source->zsource);
}

