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
void LensHaloAnaNSIE::assignParams(InputParams& params){

	// Host lens parameters
	if(!params.get("main_sigma",sigma)) error_message1("main_sigma",params.filename());
	if(!params.get("main_core",rcore)) error_message1("main_core",params.filename());
	if(!params.get("main_axis_ratio",fratio)) error_message1("main_axis_ratio",params.filename());
  else if(fratio > 1){
    ERROR_MESSAGE();
    std::cout << "parameter main_axis_ratio must be < 1 in file " << params.filename() << ". Use main_pos_angle to rotate the halo." << std::endl;
    exit(1);
  }

	if(!params.get("main_pos_angle",pa)) error_message1("main_pos_angle",params.filename());


	// Distortion of host lens parameters
	if(!params.get("main_NDistortionModes",perturb_Nmodes)) error_message1("main_NDistortionModes",params.filename());
	else if(perturb_Nmodes > 0){
		if(!params.get("main_perturb_beta",perturb_beta)) error_message1("main_perturb_beta",params.filename());
		else if(perturb_beta <= 0.0) {ERROR_MESSAGE(); cout << "perturb_beta can't be <= 0.0 in file " << params.filename(); }
		if(!params.get("main_perturb_kappa",perturb_rms[0])) error_message1("main_perturb_kappa",params.filename());
		if(!params.get("main_perturb_gamma",perturb_rms[1])) error_message1("main_perturb_gamma",params.filename());
		if(!params.get("main_perturb_monopole",perturb_rms[2])) error_message1("main_perturb_monopole",params.filename());
		if(!params.get("main_perturb_quadrapole",perturb_rms[3])) error_message1("main_perturb_quadrapole",params.filename());
		if(!params.get("main_perturb_hexopole",perturb_rms[4])) error_message1("main_perturb_hexopole",params.filename());
		if(!params.get("main_perturb_octopole",perturb_rms[5])) error_message1("main_perturb_octopole",params.filename());
    
	}

}

void LensHaloFit::PrintLens(bool show_substruct,bool show_stars){
	int i;

	LensHaloBaseNSIE::PrintLens(show_substruct,show_stars);

	// parameters of host elliptical
	cout << endl << "**Host lens model**" << endl;

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

void LensHaloAnaNSIE::PrintLens(bool show_substruct,bool show_stars){
  int i;
  
  LensHaloBaseNSIE::PrintLens(show_substruct,show_stars);
  
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


LensHaloAnaNSIE::LensHaloAnaNSIE(InputParams& params, bool verbose) : LensHaloBaseNSIE(params){
  
  assignParams(params);
  
  if(perturb_Nmodes){
    perturb_modes = new PosType[perturb_Nmodes+1];
    // zero perturbation modes until use LensHaloAnaNSIE::RandomlyDistortLens()
    for(int i=0;i< perturb_Nmodes+1 ;++i) perturb_modes[i] =  0;
  }
  
  if(verbose) PrintLens(false,false);
}
LensHaloAnaNSIE::~LensHaloAnaNSIE(){
  if(perturb_Nmodes) delete [] perturb_modes;
}

LensHaloFit::LensHaloFit(const COSMOLOGY& cosmo, int MyNmodes, PosType beta,PosType zlensref, PosType zsourceref) : LensHaloBaseNSIE(){
  
  //assignParams(params);
  
  zlens = zlensref ;
  // zsource = zsourceref ;
  zsource_reference = zsourceref ;
  
  perturb_Nmodes = MyNmodes ;
  
  if(perturb_Nmodes){
    perturb_modes = new PosType[perturb_Nmodes+1];
    // zero perturbation modes until use LensHaloAnaNSIE::RandomlyDistortLens()
    for(int i=0;i< perturb_Nmodes+1 ;++i) perturb_modes[i] =  0;
  }
  
  sigma = 0.0;
  setCosmology(cosmo);

  Sigma_crit = cosmo.SigmaCrit(zlens, zsource_reference);
  
  perturb_beta = beta;
  Dl = cosmo.angDist(zlens);
  Dls = cosmo.angDist(zlens,zsource_reference);
  Ds = cosmo.angDist(zsource_reference);
  //if(verbose) PrintLens(false,false);
}



LensHaloFit::~LensHaloFit(){
}
