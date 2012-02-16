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
	/*Dl = cosmo->angDist(0,lens->zlens);
	Ds = cosmo->angDist(0,source->zsource);
	Dls = cosmo->angDist(lens->zlens,source->zsource);*/

	source->DlDs = cosmo->angDist(0,lens->zlens) / cosmo->angDist(0,source->zsource);

	lens->setInternalParams(cosmo,source->zsource);
}

