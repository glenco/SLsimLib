/*
 * readfiles_ana.c
 *
 *  Created on: Dec 8, 2009
 *      Author: R.B. Metcalf
 *
 *      reads parameters for analytic lens model
 */

#include "slsimlib.h"

using namespace std;

/**
 * \brief Reads in a parameter file and sets up an analytic lens.
 *
 * Sets many parameters within the lens model, source model and
 * force calculation.
 */
void AnaNSIELensHalo::assignParams(InputParams& params){

	// Host lens parameters
	if(!params.get("sigma",sigma)) error_message1("sigma",params.filename());
	if(!params.get("core",rcore)) error_message1("core",params.filename());
	if(!params.get("axis_ratio",fratio)) error_message1("axis_ratio",params.filename());
	if(!params.get("pos_angle",pa)) error_message1("pos_angle",params.filename());


	// Distortion of host lens parameters
	if(!params.get("NDistortionModes",perturb_Nmodes)) error_message1("NDistortionModes",params.filename());
	else if(perturb_Nmodes > 0){
		if(!params.get("beta_perturb",perturb_beta)) error_message1("beta_perturb",params.filename());
		else if(perturb_beta <= 0.0) {ERROR_MESSAGE(); cout << "perturb_beta can't be <= 0.0 in file " << params.filename(); }
		if(!params.get("kappa_perturb",perturb_rms[0])) error_message1("kappa_perturb",params.filename());
		if(!params.get("gamma_perturb",perturb_rms[1])) error_message1("gamma_perturb",params.filename());
		if(!params.get("monopole_perturb",perturb_rms[2])) error_message1("monopole_perturb",params.filename());
		if(!params.get("quadrapole_perturb",perturb_rms[3])) error_message1("quadrapole_perturb",params.filename());
		if(!params.get("hexopole_perturb",perturb_rms[4])) error_message1("hexopole_perturb",params.filename());
		if(!params.get("octopole_perturb",perturb_rms[5])) error_message1("octopole_perturb",params.filename());
	}

}

/** \ingroup ImageFinding
 * \brief Prints the parameters of the analytic lens to stdout

cout << endl << "Nstars "<<stars_N << endl << endl;
	if(stars_N>0){
		if(star_Nregions > 0)
			cout << "stars_Nregions "<< star_Nregions << endl;
		cout << "stars_massscale "<<star_massscale << endl;
		cout << "stars_fstars "<<star_fstars << endl;
		cout << "stars_theta_force "<<star_theta_force << endl;
		if(show_stars){
			if(stars_implanted){
			  for(i=0 ; i < stars_N ; ++i) cout << "    x["<<i<<"]="
							    << stars_xp[i][0] << " " << stars_xp[i][1] << endl;
			}else cout << "stars are not implanted yet" << endl;
		}
	}
*/

void AnaNSIELensHalo::PrintLens(bool show_substruct,bool show_stars){
	int i;

	BaseNSIELensHalo::PrintLens(show_substruct,show_stars);

	// parameters of host elliptical
	cout << endl << "**Host lens model**" << endl;
	// redshifts
	cout << "sigma " << sigma << "km/s" << endl;
	cout << "core " << rcore << " Mpc" << endl;
	cout << "axis_ratio " << fratio << endl;
	cout << "position angle " <<pa << endl;

			// parameters of distortion to host elliptical
	cout << endl << "Nmodes " << perturb_Nmodes << endl;
	if(perturb_Nmodes>0){
		cout << "beta = " << perturb_beta << endl;
		cout << "rms" << endl;
		for(i=0;i<6;++i) cout << "  " << perturb_rms[i] << endl;
		cout << "modes" << endl;
		for(i=0;i<perturb_Nmodes;++i) cout << "  " << perturb_modes[i] << endl;
	}


}


AnaNSIELensHalo::AnaNSIELensHalo(InputParams& params) : BaseNSIELensHalo(params){

  assignParams(params);

  if(perturb_Nmodes){
  	perturb_modes = new double[perturb_Nmodes+1];
  	// zero perturbation modes until use BaseNSIELensHalo::RandomlyDistortLens()
  	for(int i=0;i< perturb_Nmodes+1 ;++i) perturb_modes[i] =  0;
  }

  // in degrees
  pa*=pi/180;

  PrintLens(false,false);
}


AnaNSIELensHalo::~AnaNSIELensHalo(){

}
