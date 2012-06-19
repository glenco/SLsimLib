/*
 * MOKAlens.h
 *
 *  Created on: Jun 8, 2012
 *      Author: mpetkova
 */

#include <slsimlib.h>
#include <MOKAfits.h>
#include "profile.h"


#ifndef MOKALENS_H_
#define MOKALENS_H_


/// A class to represents the MOKA lens map
class MOKALens : public Lens{
public:

	MOKALens(std::string);
	~MOKALens();

	bool set;	/// the name of the MOKA input file
	std::string MOKA_input_file;

	void readParamfile(std::string);
	void rayshooterInternal(double *ray, double *alpha, double *gamma, double *kappa, bool kappa_off);
	void setZlens(double zlens);
	double getZlens();
	void setInternalParams(CosmoHndl,SourceHndl);
	void saveImage(GridHndl grid, bool saveprofile=true);
	void saveProfile();

	MOKAmap *map;

};

#endif /* MOKALENS_H_ */
