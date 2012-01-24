/*
 * model.h
 *
 *  Created on: Jan 23, 2012
 *      Author: mpetkova
 */

#ifndef MODEL_H_
#define MODEL_H_

#include <lens.h>
#include <cosmo.h>

class Model{
public:
	Lens *lens;
	Source *source;
	CosmoHndl cosmo;

	Model();
	~Model();
};

typedef Model *ModelHndl;

#endif /* MODEL_H_ */
