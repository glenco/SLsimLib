/*
 * base_analens.cpp
 *
 *  Created on: Jan 22, 2013
 *      Author: mpetkova
 */

#include "slsimlib.h"

using namespace std;

/**
 * \brief Reads in a parameter file and sets up an analytic lens.
 *
 * Sets many parameters within the lens model, source model and
 * force calculation.
 */
void BaseAnaLens::assignParams(InputParams& params){

	// Host lens parameters
	if(!params.get("z_lens",zlens)) error_message1("z_lens",params.filename());

    // Substructure parameters
    if(!params.get("NdensitySubstruct",sub_Ndensity)) error_message1("NdensitySubstruct",params.filename());
    else if(sub_Ndensity > 0){
    	if(!params.get("beta_sub",sub_beta)) error_message1("beta_sub",params.filename());
    	if(!params.get("alpha_sub",sub_alpha)) error_message1("alpha_sub",params.filename());
    	if(!params.get("R_submax",sub_Rmax)) error_message1("R_submax",params.filename());
    	if(!params.get("mass_max",sub_Mmax)) error_message1("mass_max",params.filename());
    	if(!params.get("mass_min",sub_Mmin)) error_message1("mass_min",params.filename());
    	if(sub_Mmin < 1.0e3){
    		ERROR_MESSAGE();
    		std::cout << "Are you sure the minimum halo mass should be " << sub_Mmin << " Msun?" << std::endl;
    		exit(1);
    	}
    	if(!params.get("sub_type",sub_type)) error_message1("sub_type",params.filename());
    }
	  // Stars parameters
    if(!params.get("Nstars",stars_N)) error_message1("Nstars",params.filename());
    else if(stars_N){
    	if(!params.get("fstars",star_fstars)) error_message1("fstars",params.filename());
    	if(star_fstars < 0 || star_fstars > 1){
    		ERROR_MESSAGE();
    		cout << "fstars cannot be less than 0 or larger than 1 in file " << params.filename() <<endl;
    		exit(0);
    	}
    	if(!params.get("stars_mass",star_massscale)) error_message1("stars_mass",params.filename());
    }

}

void BaseAnaLens::error_message1(std::string parameter,std::string file){
		  ERROR_MESSAGE();
		  std::cout << "Parameter " << parameter << " is needed to construct a BaseAnaLens.  It needs to be set in parameter file " << file << "!" << endl;
		  exit(0);
}


double BaseAnaLens::getZlens(){
	return zlens;
}
/// resets Zl, Dl, Sigma_crit, MpcToAsec
void BaseAnaLens::setZlens(CosmoHndl cosmo,double zl,double zsource){
	zlens = zl;
	setInternalParams(cosmo, zsource);
}

void BaseAnaLens::reNormSubstructure(double kappa_sub){
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

/// Sets parameters within BaseLens that depend on the source redshift - Dl,Sigma_crit,etc.
void BaseAnaLens::setInternalParams(CosmoHndl cosmo, SourceHndl source){
	setInternalParams(cosmo,source->getZ());
}
void BaseAnaLens::setInternalParams(CosmoHndl cosmo, double zsource){
	double Ds, Dls;

	if(zsource < zlens) zsource = 1000;
	Dl = cosmo->angDist(0,zlens);
	Ds = cosmo->angDist(0,zsource);
	Dls = cosmo->angDist(zlens,zsource);

	MpcToAsec = 60*60*180 / pi / Dl;
		// in Mpc
	host_ro=4*pi*pow(host_sigma/2.99792e5,2)*Dl
		*Dls/Ds;
	// find critical density
	Sigma_crit=Ds/Dls/Dl/4/pi/Grav;
	to = (1+zlens)*Ds/Dls/Dl/8.39428142e-10;
}

BaseAnaLens::BaseAnaLens(InputParams& params) : Lens(){

  perturb_rms = new double[6];

  assignParams(params);

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

  perturb_Nmodes = 0;
  sub_sigmaScale = host_sigma = host_pos_angle = host_ro = host_axis_ratio = host_core = 0.0;

  if(sub_Ndensity == 0)
	  sub_N = 0;

  Sigma_crit = 0;

  stars_implanted = false;
  set = true;
}


void BaseAnaLens::PrintLens(bool show_substruct,bool show_stars){
	int i;
	cout << "zlens " << zlens << endl;

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


BaseAnaLens::~BaseAnaLens(){
	cout << "deleting lens" << endl;

	delete[] perturb_rms;

	if(perturb_Nmodes > 0){
		cout << "deleting modes" << endl;
		delete[] perturb_modes;
	}
	if(sub_N > 0 && substruct_implanted){
		cout << "deleting subs" << endl;
		free_PosTypeMatrix(sub_x,sub_N,2);
		delete[] sub_Rcut;
		delete[] sub_mass;
		delete[] sub_substructures;
	}
	if(stars_N > 0 && stars_implanted){
		cout << "deleting stars" << endl;
		delete[] star_masses;
		delete[] stars;
		free_PosTypeMatrix(stars_xp,stars_N,3);
		delete[] star_region;
		delete[] star_kappa;
		free_PosTypeMatrix(star_xdisk,star_Nregions,2);
		delete star_tree;
	}
}


