/*
 * readfiles_ana.c
 *
 *  Created on: Dec 8, 2009
 *      Author: R.B. Metcalf
 *
 *      reads parameters for analytic lens model
 */

#include <slsimlib.h>
#include <sstream>

using namespace std;

/// added definition of Grav -- needed it setParams(cosmo)

#ifndef Grav
#define Grav 4.7788e-20
#endif

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
	if(!params.get("z_lens",zlens)) error_message1("z_lens",params.filename());


	// Distortion of host lens parameters
	if(!params.get("NDistortionModes",perturb_Nmodes)) error_message1("NDistortionModes",params.filename());
	else if(perturb_Nmodes > 0){
		if(!params.get("beta_perturb",perturb_beta)) error_message1("beta_perturb",params.filename());
		if(!params.get("kappa_peturb",perturb_rms[0])) error_message1("kappa_peturb",params.filename());
		if(!params.get("gamma_peturb",perturb_rms[1])) error_message1("gamma_peturb",params.filename());
		if(!params.get("monopole_peturb",perturb_rms[2])) error_message1("monopole_peturb",params.filename());
		if(!params.get("quadrapole_peturb",perturb_rms[3])) error_message1("quadrapole_peturb",params.filename());
		if(!params.get("hexopole_peturb",perturb_rms[4])) error_message1("hexopole_peturb",params.filename());
		if(!params.get("octopole_peturb",perturb_rms[5])) error_message1("octopole_peturb",params.filename());
	}
    // Substructure parameters
    if(!params.get("NdensitySubstruct",sub_Ndensity)) error_message1("NdensitySubstruct",params.filename());
    else if(sub_Ndensity > 0){
    	if(!params.get("beta_sub",sub_beta)) error_message1("beta_sub",params.filename());
    	if(!params.get("alpha_sub",sub_alpha)) error_message1("alpha_sub",params.filename());
    	if(!params.get("R_submax",sub_Rmax)) error_message1("R_submax",params.filename());
    	if(!params.get("mass_max",sub_Mmax)) error_message1("mass_max",params.filename());
    	if(!params.get("mass_min",sub_Mmin)) error_message1("mass_min",params.filename());
    	if(!params.get("sub_type",sub_type)) error_message1("sub_type",params.filename());
    }
	  // Stars parameters
    if(!params.get("Nstars",stars_N)) error_message1("Nstars",params.filename());
    else if(stars_N){
    	if(!params.get("fstars",star_fstars)) error_message1("fstars",params.filename());
    	if(!params.get("stars_mass",star_massscale)) error_message1("stars_mass",params.filename());
    }

    if(sub_Ndensity > 0){
    	switch(sub_type){
    	case nfw:
		  sub_alpha_func = alphaNFW;
		  sub_kappa_func = kappaNFW;
		  sub_gamma_func = gammaNFW;
		  sub_phi_func = 0;
		  ERROR_MESSAGE();
		  break;
    	case powerlaw:
		  sub_alpha_func = alphaPowLaw;
		  sub_kappa_func = kappaPowLaw;
		  sub_gamma_func = gammaPowLaw;
		  sub_phi_func = phiPowLaw;
		  break;
    	case pointmass:
		  sub_alpha_func = NULL;
		  sub_kappa_func = NULL;
		  sub_gamma_func = NULL;
		  sub_phi_func = NULL;
		  break;
    	default:
		  ERROR_MESSAGE();
		  cout << "ERROR: no submass internal profile chosen" << endl;
		  exit(1);
		  break;
    	}
    }

    // parameters for stars
    stars_implanted = false; // stars are implanted later
    star_theta_force = 0.1;
    sub_theta_force = 0.1;

    sub_sigmaScale = host_sigma;

    if(sub_Ndensity == 0)
	  sub_N = 0;

    Sigma_crit = 0;

    // in degrees
    host_pos_angle*=pi/180;
    if(perturb_Nmodes)
    	perturb_modes = new double[perturb_Nmodes+1];

    PrintAnaLens(false,false);
}

void AnaLens::error_message1(std::string parameter,std::string file){
		  ERROR_MESSAGE();
		  std::cout << "Parameter " << parameter << " is needed to construct a AnaLens.  It needs to be set in parameter file " << file << "!" << endl;
		  exit(0);
}


double AnaLens::getZlens(){
	return zlens;
}

void AnaLens::setZlens(double z){
	zlens = z;
}

/** \ingroup ImageFinding
 * \brief Prints the parameters of the analytic lens to stdout
 */
void AnaLens::PrintAnaLens(bool show_substruct,bool show_stars){
	int i;

	// parameters of host elliptical
	cout << endl << "**Host lens model**" << endl;
	// redshifts
	cout << "zlens " << zlens << endl;

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

	  // parameters of substructures
	cout << endl << "NdensitySubstruct "<< sub_Ndensity << endl;
	if(sub_Ndensity > 0){
		cout << "betaSubstruct "<<sub_beta << endl;
		cout << "alphaSubstruct "<<sub_alpha << endl;
		cout << "RmaxSubstruct "<<sub_Rmax << " Mpc" << endl;
		cout << "MmaxSubstruct "<<sub_Mmax << " Msun" << endl;
		cout << "MminSubstruct "<<sub_Mmin << " Msun\n" << endl;
	}

	if(sub_N > 0){
		cout << endl << "NSubstruct "<< sub_N << endl;
		if(show_substruct){
			if(substruct_implanted || sub_N > 0){
				for(i=0;i<sub_N;++i){
				  cout << "RcutSubstruct "<<i << " " <<sub_Rcut[i] << " Mpc" << endl;
				  cout << "massSubstruct "<<i<<" "<<sub_mass[i] << " Msun" << endl;
				  cout << "xSubstruct "<<i<<" "<<sub_x[i][0]<<" "<<sub_x[i][1] << " Mpc" << endl;
					switch(sub_type){
					case nfw:
						cout << "  NFW clumps" << endl;
						break;
					case powerlaw:
						cout << "  Power Law clumps" << endl;
						break;
					case pointmass:
						cout << "  Point Mass clumps" << endl;
						break;
					default:
						ERROR_MESSAGE();
						cout << "ERROR: no submass internal profile chosen" << endl;
						exit(1);
						break;
					}
				}
			}else cout << "substructures are not implanted yet" << endl;
		}
	}

	cout << endl << "Nstars "<<stars_N << endl << endl;
	if(stars_N>0){
		if(star_Nregions > 0)
			cout << "stars_Nregions "<<star_Nregions << endl;
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

	if(Sigma_crit)
		cout << "critical density is " << Sigma_crit << " Msun/Mpc^2" << endl << endl;
}

void AnaLens::reNormSubstructure(double kappa_sub){
	/* renomalizes substructure so that
	 * the average surface density it kappa_sub
	 */
	  double avem;
	  avem=sub_Mmax*(sub_alpha+1)
	    /(sub_alpha+2)*(1-pow(sub_Mmin/sub_Mmax,sub_alpha+2))/
	    (1-pow(sub_Mmin/sub_Mmax,sub_alpha+1));

	  sub_Ndensity=kappa_sub*Sigma_crit/avem;

	  return ;
}

void AnaLens::setInternalParams(CosmoHndl cosmo, SourceHndl source){
	double Ds, Dls;

	Dl = cosmo->angDist(0,zlens);
	Ds = cosmo->angDist(0,source->getZ());
	Dls = cosmo->angDist(zlens,source->getZ());

	MpcToAsec = 60*60*180 / pi / Dl;
		// in Mpc
	host_ro=4*pi*pow(host_sigma/2.99792e5,2)*Dl
		*Dls/Ds;
	// find critical density
	Sigma_crit=Ds/Dls/Dl/4/pi/Grav;
	to = (1+zlens)*Ds/Dls/Dl/8.39428142e-10;
}

AnaLens::AnaLens(InputParams& params) : Lens(){

  perturb_rms = new double[6];

  assignParams(params);

  set = true;
}


AnaLens::~AnaLens(){
	cout << "deleting lens" << endl;

	delete[] perturb_rms;

	if(perturb_Nmodes > 0){
		cout << "deleting modes" << endl;
		delete[] perturb_modes;
	}
	if(sub_N > 0 && substruct_implanted){
		cout << "deleting subs" << endl;
		free_dmatrix(sub_x,0,sub_N-1,0,1);
		delete[] sub_Rcut;
		delete[] sub_mass;
		delete[] sub_substructures;
	}
	if(stars_N > 0 && stars_implanted){
		cout << "deleting stars" << endl;
		delete[] star_masses;
		delete[] stars;
		free_PosTypeMatrix(stars_xp,0,stars_N-1,0,2);
		delete[] star_region;
		delete[] star_kappa;
		free_dmatrix(star_xdisk,0,star_Nregions-1,0,1);
		delete star_tree;
	}
}
