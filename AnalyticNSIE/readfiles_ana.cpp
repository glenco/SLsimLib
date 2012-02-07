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

/** \ingroup ImageFinding
 * \brief Reads in a parameter file and sets up an analytic lens.
 *
 * Sets many parameters within the lens model, source model and
 * force calculation.
 */

void AnaLens::readParamfile(string filename){
  const int MAXPARAM = 50;
  string label[MAXPARAM], rlabel, rvalue;
  void *addr[MAXPARAM];
  int id[MAXPARAM];
  stringstream ss;
  int i, n, type;
  double tmp = 0;
  int myint;
  double mydouble;
  string mystring;
  
  perturb_rms=(double *)calloc(6,sizeof(double));

  n = 0;

  addr[n] = &outputfile;
  id[n] = 2;
  label[n++] = "outputfile";

  addr[n] = &host_sigma;
  id[n] = 0;
  label[n++] = "sigma";

  addr[n] = &host_core;
  id[n] = 0;
  label[n++] = "core";

  addr[n] = &host_axis_ratio;
  id[n] = 0;
  label[n++] = "axis_ratio";

  addr[n] = &host_pos_angle;
  id[n] = 0;
  label[n++] = "pos_angle";

  addr[n] = &perturb_Nmodes;
  id[n] = 1;
  label[n++] = "NDistortionModes";

  addr[n] = &perturb_beta;
  id[n] = 0;
  label[n++] = "beta_perturb";

  addr[n] = perturb_rms++;
  id[n] = 0;
  label[n++] = "kappa_peturb";

  addr[n] = perturb_rms++;
  id[n] = 0;
  label[n++] = "gamma_peturb";

  addr[n] = perturb_rms++;
  id[n] = 0;
  label[n++] = "monopole_peturb";

  addr[n] = perturb_rms++;
  id[n] = 0;
  label[n++] = "quadrapole_peturb";

  addr[n] = perturb_rms++;
  id[n] = 0;
  label[n++] = "hexopole_peturb";

  addr[n] = perturb_rms++;
  id[n] = 0;
  label[n++] = "octopole_peturb";

  addr[n] = &sub_Ndensity;
  id[n] = 0;
  label[n++] = "NdensitySubstruct";

  addr[n] = &sub_beta;
  id[n] = 0;
  label[n++] = "beta_sub";

  addr[n] = &sub_alpha;
  id[n] = 0;
  label[n++] = "alpha_sub";

  addr[n] = &sub_Rmax;
  id[n] = 0;
  label[n++] = "R_submax";

  addr[n] = &sub_Mmax;
  id[n] = 0;
  label[n++] = "mass_max";

  addr[n] = &sub_Mmin;
  id[n] = 0;
  label[n++] = "mass_min";

  addr[n] = &sub_type;
  id[n] = 1;
  label[n++] = "sub_type";

  addr[n] = &stars_N;
  id[n] = 1;
  label[n++] = "Nstars";

  addr[n] = &star_fstars;
  id[n] = 0;
  label[n++] = "fstars";

  addr[n] = &star_massscale;
  id[n] = 0;
  label[n++] = "stars_mass";

  addr[n] = &zlens;
  id[n] = 0;
  label[n++] = "z_lens";


  // opening file
  cout << "analytic lens: reading from " << filename << endl;

  ifstream file_in(filename.c_str());
  if(!file_in){
    cout << "Can't open file " << filename << endl;
    exit(1);
  }

  // output file
  while(!file_in.eof()){
	  file_in >> rlabel >> rvalue;

	  for(i = 0; i < n; i++){
		  if(rlabel == label[i]){

			  ss << rvalue;

			  switch(id[i]){
			  case 0:
				  ss >> mydouble;
				  *((double *)addr[i]) = mydouble;
				  break;
			  case 1:
				  ss >> myint;
				  *((int *)addr[i]) = myint;
				  break;
			  case 2:
				  ss >> mystring;
				  *((string *)addr[i]) = mystring;
				  break;
			  }

			  ss.clear();
			  ss.str(string());

		  }
	  }
  }


  file_in.close();

  if(sub_Ndensity > 0){
	  switch(sub_type){
	  case NFW:
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
		  sub_alpha_func = 0;
		  sub_kappa_func = 0;
		  sub_gamma_func = 0;
		  sub_phi_func = 0;
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
  
  sub_sigmaScale = host_sigma;

  // in degrees
  host_pos_angle*=pi/180;

  PrintAnaLens(true,true);
}

/** \ingroup ImageFinding
 * \brief Prints the parameters of the analytic lens to stdout
 */
void AnaLens::PrintAnaLens(bool show_substruct,bool show_stars){
	int i;

	// parameters of host elliptical
	cout << endl << "**Host lens model**" << endl;

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
			}else cout << "substructures are implanted yet" << endl;
		}
	}
	if(stars_N>0){
		cout << "Nstars "<<stars_N << endl;
		cout << "stars_Nregions "<<star_Nregions << endl;
		cout << "stars_massscale "<<star_massscale << endl;
		cout << "stars_fstars "<<star_fstars << endl;
		cout << "stars_theta_force "<<star_theta_force << endl;
		if(show_stars){
			if(stars_implanted){
			  for(i=0 ; i < stars_N ; ++i) cout << "    x["<<i<<"]="
							    << stars_xp[i][0] << " " << stars_xp[i][1] << endl;
			}else cout << "stars not implanted yet" << endl;
		}
	}

	// redshifts
	cout << "zlens " << zlens << endl;

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


AnaLens::AnaLens(string filename) : Lens(){
  readParamfile(filename);

  set = true;
}

AnaLens::~AnaLens(){
	free(perturb_rms);
	if(perturb_Nmodes > 0){
	free(perturb_modes);
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
