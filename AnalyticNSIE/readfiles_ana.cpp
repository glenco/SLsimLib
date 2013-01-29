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
void AnaLens::assignParams(InputParams& params){

	// Host lens parameters
	if(!params.get("sigma",host_sigma)) error_message1("sigma",params.filename());
	if(!params.get("core",host_core)) error_message1("core",params.filename());
	if(!params.get("axis_ratio",host_axis_ratio)) error_message1("axis_ratio",params.filename());
	if(!params.get("pos_angle",host_pos_angle)) error_message1("pos_angle",params.filename());


	// Distortion of host lens parameters
	if(!params.get("NDistortionModes",perturb_Nmodes)) error_message1("NDistortionModes",params.filename());
	else if(perturb_Nmodes > 0){
		if(!params.get("beta_perturb",perturb_beta)) error_message1("beta_perturb",params.filename());
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

void AnaLens::PrintLens(bool show_substruct,bool show_stars){
	int i;

	BaseAnaLens::PrintLens(show_substruct,show_stars);

	// parameters of host elliptical
	cout << endl << "**Host lens model**" << endl;
	// redshifts
	cout << "sigma " << host_sigma << "km/s" << endl;
	cout << "core " << host_core << " Mpc" << endl;
	cout << "axis_ratio " << host_axis_ratio << endl;
	cout << "position angle " <<host_pos_angle << endl;

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


AnaLens::AnaLens(InputParams& params) : BaseAnaLens(params){

  assignParams(params);

  if(perturb_Nmodes){
  	perturb_modes = new double[perturb_Nmodes+1];
  	// zero perturbation modes until use BaseAnaLens::RandomlyDistortLens()
  	for(int i=0;i< perturb_Nmodes+1 ;++i) perturb_modes[i] =  0;
  }

  // in degrees
  host_pos_angle*=pi/180;
  set = true;

  PrintLens(false,false);
}


AnaLens::~AnaLens(){

}
