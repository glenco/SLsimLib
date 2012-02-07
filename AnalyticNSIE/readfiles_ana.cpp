/*
 * readfiles_ana.c
 *
 *  Created on: Dec 8, 2009
 *      Author: R.B. Metcalf
 *
 *      reads parameters for analytic lens model
 */

#include <slsimlib.h>

using namespace std;

/// added definition of Grav -- needed it setParams(cosmo)

#ifndef Grav
#define Grav 4.7788e-20
#endif

/** \ingroup ImageFinding
 * \brief Reads in a parameter file and sets up an analytic lens.
 *
 * Sets many parameters within the lens model, source model and
 * force calculation.
 */

void AnaLens::ReadParams_AnaLens(char *filename){
  ifstream file_in(filename);
  char label[20];
  int i, type;
  double tmp = 0;
  
  cout << "reading from " << filename << endl;

  if(!file_in){
    cout << "Can't open file " << filename << endl;
    exit(1);
  }

  // output file
  file_in >> label >> outputfile;
  cout << label << " " <<  outputfile << endl << endl;

  // parameters of host elliptical
  cout << "Host lens model" << endl;

  file_in >> label >> host_sigma;
  cout << label << " " <<  host_sigma << " km/s" << endl;

  file_in >> label >> host_core;
  cout << label << " " <<  host_core << " Mpc" << endl;

  file_in >> label >> host_axis_ratio;
  cout << label << " " <<  host_axis_ratio << endl;

  file_in >> label >> host_pos_angle;
  cout << label << " " <<  host_pos_angle << endl;

  // parameters of distortion to host elliptical
  cout << "Distortion mode" << endl;
  
  file_in >> label >> perturb_Nmodes;
  cout << label << " " <<  perturb_Nmodes;

  if(perturb_Nmodes > 0){
	  perturb_modes=(double *)calloc(perturb_Nmodes+1,sizeof(double));
	  perturb_rms=(double *)calloc(6,sizeof(double));

	  file_in >> label >> perturb_beta;
	  cout << label << " " <<  perturb_beta << endl;

	  // kappa
	  file_in >> label >> perturb_rms[0];
	  cout << label << " " <<  perturb_rms[0] << endl;
	  // gamma
	  file_in >> label >> perturb_rms[1];
	  cout << label << " " <<  perturb_rms[1] << endl;
	  // monopole
	  file_in >> label >> perturb_rms[2];
	  cout << label << " " <<  perturb_rms[2] << endl;
	  // quadropole
	  file_in >> label >> perturb_rms[3];
	  cout << label << " " <<  perturb_rms[3] << endl;
	  // hexopole
	  file_in >> label >> perturb_rms[4];
	  cout << label << " " <<  perturb_rms[4] << endl;
	  // octopole
	  file_in >> label >> perturb_rms[5];
	  cout << label << " " <<  perturb_rms[5] << endl;

	 }else{
  for(i=0;i<7;++i) file_in >> label >> tmp;
  }

  // parameters of substructures
  cout << "**Substructures**" << endl;
  substruct_implanted = false;  // substructures are implanted later where mem is allocated

  file_in >> label >> sub_Ndensity;
  cout << label << " " <<  sub_Ndensity << endl;

  file_in >> label >> sub_beta;
  file_in >> label >> sub_alpha;
  file_in >> label >> sub_Rmax;
  file_in >> label >> sub_Mmax;
  file_in >> label >> sub_Mmin;
  file_in >> label >> type;
  sub_type = (ClumpInternal)type;

  if(sub_Ndensity > 0){
    cout << label << " " <<  sub_beta << endl;
    cout << label << " " <<  sub_alpha << endl;
    cout << label << " " <<  sub_Rmax << endl;
    cout << label << " " <<  sub_Mmax << endl;
    cout << label << " " <<  sub_Mmin << endl;
    cout << label << " " <<  sub_type << endl;
  }

  if(sub_Ndensity > 0){
	  switch(sub_type){
	  case NFW:
		  sub_alpha_func = alphaNFW;
		  sub_kappa_func = kappaNFW;
		  sub_gamma_func = gammaNFW;
		  sub_phi_func = 0;
		  ERROR_MESSAGE();
		  cout << "no time delay function defined for NFW" << endl;
		  cout << "NFW clumps" << endl;
		  break;
	  case powerlaw:
		  sub_alpha_func = alphaPowLaw;
		  sub_kappa_func = kappaPowLaw;
		  sub_gamma_func = gammaPowLaw;
		  sub_phi_func = phiPowLaw;
		  cout << "Power Law clumps" << endl;
		  break;
	  case pointmass:
		  sub_alpha_func = 0;
		  sub_kappa_func = 0;
		  sub_gamma_func = 0;
		  sub_phi_func = 0;
		  cout << "Point Mass clumps" << endl;
		  break;
	  default:
		  ERROR_MESSAGE();
		  cout << "ERROR: no submass internal profile chosen" << endl;
		  exit(1);
		  break;
	  }

  }

  // parameters for stars
  cout << "**Stars**" << endl;
  stars_implanted = false; // stars are implanted later

  file_in >> label >> stars_N;
  cout << label << " " <<  stars_N << endl;

  file_in >> label >> star_fstars;
  cout << label << " " <<  star_fstars << endl;

  file_in >> label >> star_massscale;
  cout << label << " " <<  star_massscale << endl;
  // redshifts

  file_in >> label >> zlens;
  cout << label << " " <<  zlens << endl;
  
  file_in.close();

  sub_sigmaScale = host_sigma;

  // in degrees
  host_pos_angle*=pi/180;
}

/** \ingroup ImageFinding
 * \brief Prints the parameters of the analytic lens to stdout
 */
void AnaLens::PrintAnaLens(bool show_substruct,bool show_stars){
	int i;

	cout << "Output file " << outputfile << endl;
	// parameters of host elliptical
	cout << "Host lens model" << endl;
	cout << "sigma " << host_sigma << "km/s" << endl;
	cout << "core " << host_core << " Mpc" << endl;
	cout << "axis_ratio " << host_axis_ratio << endl;
	cout << "position angle " <<host_pos_angle << endl;

			// parameters of distortion to host elliptical
	cout << "Nmodes " << perturb_Nmodes << endl;
	cout << "beta = " << perturb_beta << endl;
	if(perturb_Nmodes>0){
		cout << "rms" << endl;
		for(i=0;i<6;++i) cout << "  " << perturb_rms[i] << endl;
		cout << "modes" << endl;
		for(i=0;i<perturb_Nmodes;++i) cout << "  " << perturb_modes[i] << endl;
	}

	  // parameters of substructures
	cout << "Substructures" << endl;
	cout << "NdensitySubstruct "<<sub_Ndensity << endl;
	cout << "NSubstruct "<<sub_N << endl;

	if(sub_N > 0){
		cout << "betaSubstruct "<<sub_beta << endl;
		cout << "alphaSubstruct "<<sub_alpha << endl;
		cout << "RmaxSubstruct "<<sub_Rmax << " Mpc" << endl;
		cout << "MmaxSubstruct "<<sub_Mmax << " Msun" << endl;
		cout << "MminSubstruct "<<sub_Mmin << " Msun\n" << endl;

		if(show_substruct){
			if(substruct_implanted){
				for(i=0;i<sub_N;++i){
				  cout << "RcutSubstruct "<<i << " " <<sub_Rcut[i] << " Mpc" << endl;
				  cout << "massSubstruct "<<i<<" "<<sub_mass[i] << " Msun" << endl;
				  cout << "xSubstruct "<<i<<" "<<sub_x[i][0]<<" "<<sub_x[i][1] << " Mpc" << endl;
					switch(sub_type){
					case NFW:
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
			}else cout << "  substructures are implanted yet" << endl;
		}
	}
	if(stars_N>0){
		cout << "Nstars="<<stars_N << endl;
		cout << "stars_Nregions "<<star_Nregions << endl;
		cout << "stars_massscale "<<star_massscale << endl;
		cout << "stars_fstars "<<star_fstars << endl;
		cout << "stars_theta_force "<<star_theta_force << endl;
		if(show_stars){
			if(stars_implanted){
			  for(i=0 ; i < stars_N ; ++i) cout << "    x["<<i<<"]="
							    << stars_xp[i][0] << " " << stars_xp[i][1] << endl;
			}else cout << " stars not implanted yet" << endl;
		}
	}

	// redshifts
	cout << "zlens " << zlens << endl;

	cout << "critical density is " << Sigma_crit << " Msun/Mpc^2" << endl;
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


AnaLens::AnaLens(char filename[]) : Lens(){
  ReadParams_AnaLens(filename);

  set = true;
}

AnaLens::~AnaLens(){
	if(perturb_Nmodes > 0){
	free(perturb_modes);
	free(perturb_rms);
	}
	if(sub_N > 0 && substruct_implanted){
	  free_dmatrix(sub_x,0,sub_N-1,0,1);
	  free(sub_Rcut);
	  free(sub_mass);
	  free(sub_substructures);
	  
	}
	if(stars_N > 0 && stars_implanted){
		free(star_masses);
		free(stars);
		free_PosTypeMatrix(stars_xp,0,stars_N-1,0,2);
		free(star_region);
		free(star_kappa);
		free_dmatrix(star_xdisk,0,star_Nregions-1,0,1);
		delete star_tree;
	}
}
