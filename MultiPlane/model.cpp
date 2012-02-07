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
	Dl = cosmo->angDist(0,lens->zlens);
	Ds = cosmo->angDist(0,source->zsource);
	Dls = cosmo->angDist(lens->zlens,source->zsource);

	source->DlDs = Dl / Ds;

	lens->MpcToAsec = 60*60*180 / pi / Dl;

	  // in Mpc
	lens->host_ro=4*pi*pow(lens->host_sigma/2.99792e5,2)*Dl
			  *Dls/Ds/(1+lens->zlens);

	  // find critical density
	lens->Sigma_crit=Ds/Dls/Dl/4/pi/Grav;

	lens->to = (1+lens->zlens)*Ds/Dls/Dl/8.39428142e-10;
}

