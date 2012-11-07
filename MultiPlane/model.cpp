/*
 * model.cpp
 *
 *  Created on: Jan 23, 2012
 *      Author: mpetkova
 */

#include "slsimlib.h"

using namespace std;

/*
 * This creates a model and constructs a Lens, a Source, and a Cosmology.
 * The Type of the lens is determined by multi_lens, the type of the source is read in from the parameter file
 * and the the cosmology is the standard one
 */
Model::Model(LensHndl mylens  /// lens
		,SourceHndl mysource /// source
		,CosmoHndl mycosmo	/// cosmology
		){

	lens = mylens;
	source = mysource;
	cosmo = mycosmo;

	/// sets some source and lens paramaters that are interdependent
	/// e.g. source->DlDs depens on source and lens
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
