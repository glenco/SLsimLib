/*
 * readfiles_ana.c
 *
 *  Created on: Dec 8, 2009
 *      Author: R.B. Metcalf
 *
 *      reads parameters for analytic lens model
 */

#include <slsimlib.h>

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

void AnaLens::setInternal(CosmoHndl cosmo){
   /********************************/
  /* set some internal parameters */
  /********************************/
  double NSubstructInRe = 0;

  std::cout << "In setInternal(cosmo)\n";

  sub_sigmaScale = host_sigma;
  MpcToAsec = 60*60*180 / pi / cosmo->angDist(0,zlens);
  std::cout << "Arcseconds/Mpc: " << MpcToAsec << std::endl;

  // in degrees
  host_pos_angle*=pi/180;

  // in Mpc
  host_ro=4*pi*pow(host_sigma/2.99792e5,2)*cosmo->angDist(0,zlens)
		  *cosmo->angDist(zlens,zsource)
		  /cosmo->angDist(0,zsource)/(1+zlens);

  // find critical density
  Sigma_crit=cosmo->angDist(0,zsource)/cosmo->angDist(zlens,zsource)
		  /cosmo->angDist(0,zlens)/4/pi/Grav;
  to = (1+zlens)*cosmo->angDist(0,zsource)
		  /cosmo->angDist(zlens,zsource)/cosmo->angDist(0,zlens)
		  /8.39428142e-10;

  std::cout << "critical density is " << Sigma_crit << " Msun/Mpc^2    ro=" << host_ro << " Mpc  D_l = "
       << cosmo->angDist(0,zlens) << " Mpc D_s = " << cosmo->angDist(0,zsource) << " Mpc  to = " << to << " days/Mpc^2" << std::endl;

  NSubstructInRe = pi*pow(2*host_ro,2)*sub_Ndensity;

  if(NSubstructInRe > 0){

	  // include clumps that are beyond 2 Re
	  /*
	  NSubstruct = (int)(NSubstructInRe*FractionWithinRe(lens,2) + 0.5);
	  Rmax = ro*2 + RmaxSubstruct
		          + pow(2*MmaxSubstruct*ro/pi/Sigma_crit/1.0e-3,1./3.);
	  NSubstruct = (int)(NdensitySubstruct*pi*Rmax*Rmax+0.5);
	  */
	  std::cout << "average number of clumps including outside Re: " << pi*pow(host_ro,2)*sub_Ndensity << std::endl;
	  std::cout << "average clumps mass: " << averageSubMass() << " Msun" << std::endl;

	  /*for(i=0;i<NSubstruct;++i){
		  RcutSubstruct[i]=RmaxSubstruct;
		  massSubstruct[i]=MmaxSubstruct;
	  }*/

  }
}

void AnaLens::ReadParams_AnaLens(char *filename){
  ifstream file_in(filename);
  char label[20];
  int i, type;
  double tmp = 0;
  
  std::cout << "reading from " << filename << std::endl;

  if(!file_in){
    std::cout << "Can't open file " << filename << std::endl;
    exit(1);
  }

  // output file
  file_in >> label >> outputfile;
  std::cout << label << outputfile << std::endl << std::endl;

  // parameters of host elliptical
  std::cout << "Host lens model" << std::endl;

  file_in >> label >> host_sigma;
  std::cout << label << host_sigma << " km/s" << std::endl;

  file_in >> label >> host_core;
  std::cout << label << host_core << " Mpc" << std::endl;

  file_in >> label >> host_axis_ratio;
  std::cout << label << host_axis_ratio << std::endl;

  file_in >> label >> host_pos_angle;
  std::cout << label << host_pos_angle << std::endl;

  // parameters of distortion to host elliptical
  std::cout << "Distortion mode" << std::endl;
  
  file_in >> label >> perturb_Nmodes;
  std::cout << label << perturb_Nmodes;

  if(perturb_Nmodes > 0){
	  perturb_modes=(double *)calloc(perturb_Nmodes+1,sizeof(double));
	  perturb_rms=(double *)calloc(6,sizeof(double));

	  file_in >> label >> perturb_beta;
	  std::cout << label << perturb_beta << std::endl;

	  // kappa
	  file_in >> label >> perturb_rms[0];
	  std::cout << label << perturb_rms[0] << std::endl;
	  // gamma
	  file_in >> label >> perturb_rms[1];
	  std::cout << label << perturb_rms[1] << std::endl;
	  // monopole
	  file_in >> label >> perturb_rms[2];
	  std::cout << label << perturb_rms[2] << std::endl;
	  // quadropole
	  file_in >> label >> perturb_rms[3];
	  std::cout << label << perturb_rms[3] << std::endl;
	  // hexopole
	  file_in >> label >> perturb_rms[4];
	  std::cout << label << perturb_rms[4] << std::endl;
	  // octopole
	  file_in >> label >> perturb_rms[5];
	  std::cout << label << perturb_rms[5] << std::endl;

	 }else{
  for(i=0;i<7;++i) file_in >> label >> tmp;
  }

  // parameters of substructures
  std::cout << "**Substructures**" << std::endl;
  substruct_implanted = false;  // substructures are implanted later where mem is allocated

  file_in >> label >> sub_Ndensity;
  std::cout << label << sub_Ndensity << std::endl;

  file_in >> label >> sub_beta;
  file_in >> label >> sub_alpha;
  file_in >> label >> sub_Rmax;
  file_in >> label >> sub_Mmax;
  file_in >> label >> sub_Mmin;
  file_in >> label >> type;
  sub_type = (ClumpInternal)type;

  if(sub_Ndensity > 0){
    std::cout << label << sub_beta << std::endl;
    std::cout << label << sub_alpha << std::endl;
    std::cout << label << sub_Rmax << std::endl;
    std::cout << label << sub_Mmax << std::endl;
    std::cout << label << sub_Mmin << std::endl;
    std::cout << label << sub_type << std::endl;
  }

  if(sub_Ndensity > 0){
	  switch(sub_type){
	  case NFW:
		  sub_alpha_func = alphaNFW;
		  sub_kappa_func = kappaNFW;
		  sub_gamma_func = gammaNFW;
		  sub_phi_func = 0;
		  ERROR_MESSAGE();
		  std::cout << "no time delay function defined for NFW" << std::endl;
		  std::cout << "NFW clumps" << std::endl;
		  break;
	  case powerlaw:
		  sub_alpha_func = alphaPowLaw;
		  sub_kappa_func = kappaPowLaw;
		  sub_gamma_func = gammaPowLaw;
		  sub_phi_func = phiPowLaw;
		  std::cout << "Power Law clumps" << std::endl;
		  break;
	  case pointmass:
		  sub_alpha_func = 0;
		  sub_kappa_func = 0;
		  sub_gamma_func = 0;
		  sub_phi_func = 0;
		  std::cout << "Point Mass clumps" << std::endl;
		  break;
	  default:
		  ERROR_MESSAGE();
		  std::cout << "ERROR: no submass internal profile chosen" << std::endl;
		  exit(1);
		  break;
	  }

  }

  // parameters for stars
  std::cout << "**Stars**" << std::endl;
  stars_implanted = false; // stars are implanted later

  file_in >> label >> stars_N;
  std::cout << label << stars_N << std::endl;

  file_in >> label >> star_fstars;
  std::cout << label << star_fstars << std::endl;

  file_in >> label >> star_massscale;
  std::cout << label << star_massscale << std::endl;

  // source information
  std::cout << "**Source structure**" << std::endl;

  file_in >> label >> type;
  source_sb_type = (SBModel)type;
  std::cout << label << source_sb_type << std::endl;

  if(source_sb_type == Uniform){
		  source_sb_func = uniform_SB;
		  std::cout << "uniform surface brightness source" << std::endl;
  }else if(source_sb_type == Gaussian){
		  source_sb_func = gaussian_SB;
		  std::cout << "Gaussian surface brightness source" << std::endl;
		  file_in >> label >> source_gauss_r2;
		  std::cout << label << source_gauss_r2 << " Mpc" << std::endl;
  }else{

	  std::cout << "BLR surface brightness source" << std::endl;
	  switch(source_sb_type){
	  case BLR_Disk:
		  source_sb_func = BLR_Disk_SB;
		  std::cout << "disk model" << std::endl;
		  break;
	  case BLR_Sph1:
		  source_sb_func = BLR_Sph1_SB;
		  std::cout << "spherical with circular orbits" << std::endl;
		  break;
	  case BLR_Sph2:
		  source_sb_func = BLR_Sph2_SB;
		  std::cout << "spherical with Gaussian velocities" << std::endl;
		  break;
	  default:
		  ERROR_MESSAGE();
		  std::cout << "ERROR: no submass internal profile chosen" << std::endl;
		  exit(1);
		  break;
	  }

	  file_in >> label >> source_BHmass;
	  std::cout << label << source_BHmass << " Msun" << std::endl;

	  file_in >> label >> source_gamma;
	  std::cout << label << source_gamma << std::endl;

	  file_in >> label >> source_inclination;
	  std::cout << label << source_inclination << " deg" << std::endl;

	  file_in >> label >> source_opening_angle;
	  std::cout << label << source_opening_angle << " deg" << std::endl;


	  source_inclination *= pi/180;
	  source_opening_angle *= pi/180;

	  file_in >> label >> source_r_in;
	  std::cout << label << source_r_in << " Mpc" << std::endl;

	  file_in >> label >> source_r_out;
	  std::cout << label << source_r_out << " Mpc" << std::endl;

	  file_in >> label >> source_nuo;
	  std::cout << label << source_nuo << " Hz" << std::endl;

	  file_in >> label >> source_fK;
	  std::cout << label << source_fK << " X V_Kepler" << std::endl;

	  source_monocrome = false;  // default value

  }

  // redshifts

  file_in >> label >> zlens;
  std::cout << label << zlens << std::endl;

  file_in >> label >> zsource;
  std::cout << label << zsource << std::endl;
  
  file_in.close();
}

/** \ingroup ImageFinding
 * \brief Prints the parameters of the analytic lens to stdout
 */
void AnaLens::PrintAnaLens(bool show_substruct,bool show_stars){
	int i;

	std::cout << "Output file" << outputfile << std::endl;
	// parameters of host elliptical
	std::cout << "Host lens model" << std::endl;
	std::cout << "sigma " << host_sigma << "km/s" << std::endl;
	std::cout << "core " << host_core << " Mpc" << std::endl;
	std::cout << "axis_ratio " << host_axis_ratio << std::endl;
	std::cout << "position angle " <<host_pos_angle << std::endl;
	std::cout << "r_source on lens plane " << source_r*1.0e6 << " pc" << std::endl;

			// parameters of distortion to host elliptical
	std::cout << "Nmodes " << perturb_Nmodes << std::endl;
	std::cout << "beta = " << perturb_beta << std::endl;
	if(perturb_Nmodes>0){
		std::cout << "rms" << std::endl;
		for(i=0;i<6;++i) std::cout << "  " << perturb_rms[i] << std::endl;
		std::cout << "modes" << std::endl;
		for(i=0;i<perturb_Nmodes;++i) std::cout << "  " << perturb_modes[i] << std::endl;
	}

	std::cout << "Source" << std::endl;

	if(source_sb_type == Uniform){
		std::cout << "uniform surface brightness source" << std::endl;
	}else if(source_sb_type == Gaussian){
		std::cout << "Gaussian surface brightness source" << std::endl;
		std::cout << "sigma^2 = " << source_gauss_r2 << " Mpc^2" << std::endl;

	  }else{

	  std::cout << "BLR surface brightness source" << std::endl;
		  switch(source_sb_type){
		  std::cout << "      BH mass "<< source_BHmass << " Msun " << std::endl;
		  std::cout << "      gamma "<< source_gamma << std::endl;
		  std::cout << "      inner radius "<< source_r_in*1.0e6 << " pc" << std::endl;
		  std::cout << "      outer radius "<< source_r_out*1.0e6 << " pc" << std::endl;
		  case BLR_Disk:
		    std::cout << "    disk model" << std::endl;
			  std::cout << "      inclination "<<source_inclination << " rads" << std::endl;
			  std::cout << "      disk opening angle "<<source_opening_angle << " rads" << std::endl;
			  std::cout << "      turbulent/thermal dispersion "<<source_fK << " X V_Kepler" << std::endl;
			  break;
		  case BLR_Sph1:
			  std::cout << "    spherical with circular orbits" << std::endl;
			  break;
		  case BLR_Sph2:
			  std::cout << "    spherical with Gaussian velocities" << std::endl;
			  std::cout << "      turbulent/thermal dispersion "<<source_fK << " x V_kepler" << std::endl;
			  break;
		  default:
			  ERROR_MESSAGE();
			  std::cout << "ERROR: no submass internal profile chosen" << std::endl;
			  break;
		  }


		  if(source_monocrome) std::cout << "      monocromatic" << std::endl;
		  else std::cout << "      center of line "<<source_nuo << " Hz" << std::endl;
	  }


	  // parameters of substructures
	std::cout << "Substructures" << std::endl;
	std::cout << "NdensitySubstruct "<<sub_Ndensity << std::endl;
	std::cout << "NSubstruct "<<sub_N << std::endl;

	if(sub_N > 0){
		std::cout << "betaSubstruct "<<sub_beta << std::endl;
		std::cout << "alphaSubstruct "<<sub_alpha << std::endl;
		std::cout << "RmaxSubstruct "<<sub_Rmax << " Mpc" << std::endl;
		std::cout << "MmaxSubstruct "<<sub_Mmax << " Msun" << std::endl;
		std::cout << "MminSubstruct "<<sub_Mmin << " Msun\n" << std::endl;

		if(show_substruct){
			if(substruct_implanted){
				for(i=0;i<sub_N;++i){
				  std::cout << "RcutSubstruct "<<i << " " <<sub_Rcut[i] << " Mpc" << std::endl;
				  std::cout << "massSubstruct "<<i<<" "<<sub_mass[i] << " Msun" << std::endl;
				  std::cout << "xSubstruct "<<i<<" "<<sub_x[i][0]<<" "<<sub_x[i][1] << " Mpc" << std::endl;
					switch(sub_type){
					case NFW:
						std::cout << "  NFW clumps" << std::endl;
						break;
					case powerlaw:
						std::cout << "  Power Law clumps" << std::endl;
						break;
					case pointmass:
						std::cout << "  Point Mass clumps" << std::endl;
						break;
					default:
						ERROR_MESSAGE();
						std::cout << "ERROR: no submass internal profile chosen" << std::endl;
						exit(1);
						break;
					}
				}
			}else std::cout << "  substructures are implanted yet" << std::endl;
		}
	}
	if(stars_N>0){
		std::cout << "Nstars="<<stars_N << std::endl;
		std::cout << "stars_Nregions "<<star_Nregions << std::endl;
		std::cout << "stars_massscale "<<star_massscale << std::endl;
		std::cout << "stars_fstars "<<star_fstars << std::endl;
		std::cout << "stars_theta_force "<<star_theta_force << std::endl;
		if(show_stars){
			if(stars_implanted){
			  for(i=0 ; i < stars_N ; ++i) std::cout << "    x["<<i<<"]="
							    << stars_xp[i][0] << " " << stars_xp[i][1] << std::endl;
			}else std::cout << " stars not implanted yet" << std::endl;
		}
	}

	// redshifts
	std::cout << "zlens " << zlens << std::endl;
	std::cout << "zsource " << zsource << std::endl;

	std::cout << "critical density is " << Sigma_crit << " Msun/Mpc^2" << std::endl;
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


AnaLens::AnaLens(char filename[],CosmoHndl cosmo){
  ReadParams_AnaLens(filename);
  setInternal(cosmo);

  set = true;
  }

/** \ingroup Constructor
 *
 */
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
