/*
 * readfiles_uni.cpp
 *
 *  Created on: Jan 22, 2013
 *      Author: Leier
 */


#include "uniform_lens.h"

using namespace std;

/**
 * \brief Reads in a parameter file and sets up an uniform lens.
 *
 * Sets many parameters within the lens model, source model and
 * force calculation.
 */

UniNSIELensHalo::UniNSIELensHalo(InputParams& params) : BaseNSIELensHalo(params){

  assignParams(params);

  perturb_Nmodes=3;
  perturb_modes = new double[3];

  perturb_modes[0]=kappa_uniform;
  perturb_modes[1]=gamma_uniform[0];
  perturb_modes[2]=gamma_uniform[1];

  PrintLens(false,false);
}

UniNSIELensHalo::~UniNSIELensHalo(){

}


void UniNSIELensHalo::assignParams(InputParams& params){

	//if(perturb_Nmodes > 0){
	if(!params.get("kappa_uniform",kappa_uniform)) error_message1("kappa_uniform",params.filename());
	if(!params.get("gamma_uniform_1",gamma_uniform[0])) error_message1("gamma_uniform_1",params.filename());
	if(!params.get("gamma_uniform_2",gamma_uniform[1])) error_message1("gamma_uniform_2",params.filename());

	if(!params.get("stellar_mass_function",imf_type)) error_message1("stellar_mass_function",params.filename());
	if(!params.get("min_mstar",min_mstar)) error_message1("min_mstar",params.filename());
	if(!params.get("max_mstar",max_mstar)) error_message1("max_mstar",params.filename());
	if(!params.get("bending_point",bend_mstar)) error_message1("bending_point",params.filename());
	if(!params.get("slope_1",lo_mass_slope)) error_message1("slope_1",params.filename());
	if(!params.get("slope_2",hi_mass_slope)) error_message1("slope_2",params.filename());

	return;
}

/** \ingroup ImageFinding
 * \brief Prints the parameters of the analytic lens to stdout
 */
void UniNSIELensHalo::PrintLens(bool show_substruct,bool show_stars){
	int i;

	BaseNSIELensHalo::PrintLens(show_substruct,show_stars);

	// uni lens parameters only
	cout << endl << "**Host lens model**" << endl;
	// redshifts
	cout << "zlens " << zlens << endl;
	cout << "kappa " << kappa_uniform << endl;
	cout << "gamma " << gamma_uniform[0] << " " << gamma_uniform[1] << endl;


}


void UniNSIELensHalo::implant_stars(double x, double y, unsigned long Nregions,long *seed, IMFtype type){

	if(Nregions <= 0) return;
	Point *centers;
	gamma_uniform[2]=0.0; // TODO gamma_uniform[2] determines rotation for multiplane lens, how shall it be implemented here?
	centers = NewPointArray(Nregions,true);
	centers[0].x[0]=x;
	centers[0].x[1]=y;
	centers[0].kappa=kappa_uniform;
	centers[0].gamma[0]=gamma_uniform[0];
	centers[0].gamma[1]=gamma_uniform[1];
	centers[0].gamma[2]=gamma_uniform[2];

	BaseNSIELensHalo::implant_stars(centers,Nregions,seed,type);
}


