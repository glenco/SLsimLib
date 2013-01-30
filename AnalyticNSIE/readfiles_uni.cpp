/*
 * readfiles_uni.cpp
 *
 *  Created on: Jan 22, 2013
 *      Author: Leier
 */


#include "slsimlib.h"

using namespace std;

/**
 * \brief Reads in a parameter file and sets up an uniform lens.
 *
 * Sets many parameters within the lens model, source model and
 * force calculation.
 */

UniLens::UniLens(InputParams& params) : BaseAnaLens(params){

  assignParams(params);

  perturb_Nmodes=3;
  perturb_modes = new double[2];

  perturb_modes[0]=kappa_uniform;
  perturb_modes[1]=gamma_uniform[0];
  perturb_modes[2]=gamma_uniform[1];

  // in degrees
  set = true;

  PrintLens(false,false);
}

UniLens::~UniLens(){

}


void UniLens::assignParams(InputParams& params){

	// Distortion of host lens parameters
	if(perturb_Nmodes > 0){
		if(!params.get("kappa_uniform",kappa_uniform)) error_message1("kappa_uniform",params.filename());
		if(!params.get("gamma_uniform_1",gamma_uniform[0])) error_message1("gamma_uniform_0",params.filename());
		if(!params.get("gamma_uniform_2",gamma_uniform[1])) error_message1("gamma_uniform_1",params.filename());
	}

	return;
}

/** \ingroup ImageFinding
 * \brief Prints the parameters of the analytic lens to stdout
 */
void UniLens::PrintLens(bool show_substruct,bool show_stars){
	int i;

	BaseAnaLens::PrintLens(show_substruct,show_stars);

	// uni lens parameters only
	cout << endl << "**Host lens model**" << endl;
	// redshifts
	cout << "zlens " << zlens << endl;
	cout << "kappa " << kappa_uniform << endl;
	cout << "gamma " << gamma_uniform[0] << " " << gamma_uniform[1] << endl;

}

