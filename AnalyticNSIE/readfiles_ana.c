/*
 * readfiles_ana.c
 *
 *  Created on: Dec 8, 2009
 *      Author: R.B. Metcalf
 *
 *      reads parameters for analytic lens model
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <nr.h>
#include <nrutil.h>
#include <nrD.h>
#include <cosmo.h>
#include "analytic_lens.h"

/** \ingroup ImageFinding
 * \brief Reads in a parameter file and sets up an analytic lens.
 *
 * Sets many parameters within the lens model, source model and
 * force calculation.
 */

void ReadParams_AnaLens(char *filename,CosmoHndl cosmo,AnaLens *lens){
  FILE *file;
  char label[20];
  int i;
  double tmp=0,NSubstructInRe=0;

  printf(">reading from %s\n>",filename);

  file=fopen(filename,"r");

  // output file

  //lens->outputfile = (char *)malloc(50*sizeof(char));
  fscanf(file,"%s %s",label,&(lens->outputfile));
  printf(">%s %s\n\n",label,lens->outputfile);
  // parameters of host elliptical
  printf(">  Host lens model\n");
  fscanf(file,"%s %le",label,&(lens->host_sigma));
  printf(">%s %f km/s\n",label,lens->host_sigma);
  fscanf(file,"%s %le",label,&(lens->host_core));
  printf(">%s %f Mpc\n",label,lens->host_core);
  fscanf(file,"%s %le",label,&(lens->host_axis_ratio));
  printf(">%s %f\n",label,lens->host_axis_ratio);
  fscanf(file,"%s %le",label,&(lens->host_pos_angle));
  printf(">%s %f\n",label,lens->host_pos_angle);

  // parameters of distortion to host elliptical
  printf("\n>  Distortion mode\n");
  fscanf(file,"%s %li",label,&(lens->perturb_Nmodes));
  printf(">%s %li\n",label,lens->perturb_Nmodes);
  if(lens->perturb_Nmodes > 0){
	  lens->perturb_modes=(double *)calloc(lens->perturb_Nmodes+1,sizeof(double));
	  lens->perturb_rms=(double *)calloc(6,sizeof(double));

	  fscanf(file,"%s %le",label,&(lens->perturb_beta));
	  printf(">%s %.3f\n",label,lens->perturb_beta);

	  // kappa
	  fscanf(file,"%s %le",label,&(lens->perturb_rms[0]));
	  printf(">%s %.3e\n",label,lens->perturb_rms[0]);
	  // gamma
	  fscanf(file,"%s %le",label,&(lens->perturb_rms[1]));
	  printf(">%s %.3e\n",label,lens->perturb_rms[1]);
	  // monopole
	  fscanf(file,"%s %le",label,&(lens->perturb_rms[2]));
	  printf(">%s %.3e\n",label,lens->perturb_rms[2]);
	  // quadropole
	  fscanf(file,"%s %le",label,&(lens->perturb_rms[3]));
	  printf(">%s %.3e\n",label,lens->perturb_rms[3]);
	  // hexopole
	  fscanf(file,"%s %le",label,&(lens->perturb_rms[4]));
	  printf(">%s %.3e\n",label,lens->perturb_rms[4]);
	  // octopole
	  fscanf(file,"%s %le",label,&(lens->perturb_rms[5]));
	  printf(">%s %.3e\n",label,lens->perturb_rms[5]);

  }else{
	  for(i=0;i<7;++i) fscanf(file,"%s %le",label,&tmp);
	  //printf(">%s %.3e\n",label,tmp);
  }

  // parameters of substructures
  printf("\n>  **Substructures**\n");
  lens->substruct_implanted=false;  // substructures are implanted later where mem is allocated
  fscanf(file,"%s %le",label,&(lens->sub_Ndensity));
  printf(">%s %e\n",label,lens->sub_Ndensity);

  fscanf(file,"%s %le",label,&(lens->sub_beta));
  if(lens->sub_Ndensity > 0) printf(">%s %f\n",label,lens->sub_beta);
  fscanf(file,"%s %le",label,&(lens->sub_alpha));
  if(lens->sub_Ndensity > 0) printf(">%s %f\n",label,lens->sub_alpha);

  fscanf(file,"%s %le",label,&lens->sub_Rmax);
  if(lens->sub_Ndensity > 0) printf(">%s %.3e Mpc\n",label,lens->sub_Rmax);
  fscanf(file,"%s %le",label,&lens->sub_Mmax);
  if(lens->sub_Ndensity > 0) printf(">%s %.3e Msun\n",label,lens->sub_Mmax);
  fscanf(file,"%s %le",label,&lens->sub_Mmin);
  if(lens->sub_Ndensity > 0) printf(">%s %.3e Msun\n",label,lens->sub_Mmin);
  fscanf(file,"%s %i",label,&lens->sub_type);
  if(lens->sub_Ndensity > 0) printf(">%s %i ",label,lens->sub_type);

  if(lens->sub_Ndensity > 0){
	  switch(lens->sub_type){
	  case NFW:
		  lens->sub_alpha_func = alphaNFW;
		  lens->sub_kappa_func = kappaNFW;
		  lens->sub_gamma_func = gammaNFW;
		  lens->sub_phi_func = 0;
		  ERROR_MESSAGE();
		  printf(">no time delay function defined for NFW\n");
		  printf(">  NFW clumps\n");
		  break;
	  case powerlaw:
		  lens->sub_alpha_func = alphaPowLaw;
		  lens->sub_kappa_func = kappaPowLaw;
		  lens->sub_gamma_func = gammaPowLaw;
		  lens->sub_phi_func = phiPowLaw;
		  printf(">  Power Law clumps\n");
		  break;
	  case pointmass:
		  lens->sub_alpha_func = 0;
		  lens->sub_kappa_func = 0;
		  lens->sub_gamma_func = 0;
		  lens->sub_phi_func = 0;
		  printf(">  Point Mass clumps\n");
		  break;
	  default:
		  ERROR_MESSAGE();
		  printf(">ERROR: no submass internal profile chosen\n");
		  exit(1);
	  }

  }

  // parameters for stars
  printf("\n>  **Stars**\n");
  lens->stars_implanted=false; // stars are implanted later
  fscanf(file,"%s %li",label,&(lens->stars_N));
  printf(">%s %li\n",label,lens->stars_N);
  fscanf(file,"%s %le",label,&(lens->star_fstars));
  printf(">%s %.3f\n",label,lens->star_fstars);
  fscanf(file,"%s %le",label,&(lens->star_massscale));
  printf(">%s %e\n",label,lens->star_massscale);

  // source information
  printf("\n>  **Source structure**\n");
  fscanf(file,"%s %i",label,&(lens->source_sb_type));
  printf(">%s %i \n",label,lens->source_sb_type);
  if(lens->source_sb_type == Uniform){
		  lens->source_sb_func = uniform_SB;
		  printf(">  uniform surface brightness source\n");
  }else if(lens->source_sb_type == Gaussian){
		  lens->source_sb_func = gaussian_SB;
		  printf(">  Gaussian surface brightness source\n");
		  fscanf(file,"%s %le",label,&(lens->source_gauss_r2));
		  printf(">%s %.3e Mpc\n",label,lens->source_gauss_r2);
  }else{

	  printf(">  BLR surface brightness source\n");
	  switch(lens->source_sb_type){
	  case BLR_Disk:
		  lens->source_sb_func = BLR_Disk_SB;
		  printf(">    disk model\n");
		  break;
	  case BLR_Sph1:
		  lens->source_sb_func = BLR_Sph1_SB;
		  printf(">    spherical with circular orbits\n");
		  break;
	  case BLR_Sph2:
		  lens->source_sb_func = BLR_Sph2_SB;
		  printf(">    spherical with Gaussian velocities\n");
		  break;
	  default:
		  ERROR_MESSAGE();
		  printf(">ERROR: no submass internal profile chosen\n");
		  exit(1);
	  }

	  fscanf(file,"%s %e",label,&(lens->source_BHmass));
	  printf(">    %s %.3e Msun\n",label,lens->source_BHmass);
	  fscanf(file,"%s %e",label,&(lens->source_gamma));
	  printf(">    %s %.3f\n",label,lens->source_gamma);
	  fscanf(file,"%s %e",label,&(lens->source_inclination));
	  printf(">    %s %.3f deg\n",label,lens->source_inclination);
	  lens->source_inclination *= pi/180;
	  fscanf(file,"%s %e",label,&(lens->source_opening_angle));
	  printf(">    %s %.3f deg\n",label,lens->source_opening_angle);
	  lens->source_opening_angle *= pi/180;

	  fscanf(file,"%s %e",label,&(lens->source_r_in));
	  printf(">    %s %.3e Mpc\n",label,lens->source_r_in);
	  fscanf(file,"%s %e",label,&(lens->source_r_out));
	  printf(">    %s %.3e Mpc\n",label,lens->source_r_out);
	  fscanf(file,"%s %e",label,&(lens->source_nuo));
	  printf(">    %s %.5e Hz\n",label,lens->source_nuo);
	  fscanf(file,"%s %e",label,&(lens->source_fK));
	  printf(">    %s %.4f X V_Kepler\n",label,lens->source_fK);
	  lens->source_monocrome = false;  // default value

  }

  // redshifts
  fscanf(file,"%s %le",label,&(lens->zlens));
  printf("\n>%s %f\n",label,lens->zlens);
  fscanf(file,"%s %le",label,&(lens->zsource));
  printf(">%s %f\n",label,lens->zsource);

  // cosmology
  SetConcordenceCosmology(cosmo);
  cosmo->physical=0;

  fscanf(file,"%s %le",label,&(cosmo->Omo));
  printf(">%s %f\n",label,cosmo->Omo);
  fscanf(file,"%s %le",label,&(cosmo->Oml));
  printf(">%s %f\n",label,cosmo->Oml);
  fscanf(file,"%s %le",label,&(cosmo->h));
  printf(">%s %f\n",label,cosmo->h);
  fclose(file);

   /********************************/
  /* set some internal parameters */
  /********************************/

  lens->sub_sigmaScale=lens->host_sigma;
  lens->MpcToAsec=60*60*180/pi/angDist(0,lens->zlens,cosmo);
  printf(">Arcseconds/Mpc: %e\n",lens->MpcToAsec);
  // in degrees
  lens->host_pos_angle*=pi/180;
  // in Mpc
  lens->host_ro=4*pi*pow(lens->host_sigma/2.99792e5,2)*angDist(0,lens->zlens,cosmo)
		  *angDist(lens->zlens,lens->zsource,cosmo)
		  /angDist(0,lens->zsource,cosmo)/(1+lens->zlens);

  // find critical density
  lens->Sigma_crit=angDist(0,lens->zsource,cosmo)/angDist(lens->zlens,lens->zsource,cosmo)
		  /angDist(0,lens->zlens,cosmo)/4/pi/Grav;
  lens->to = (1+lens->zlens)*angDist(0,lens->zsource,cosmo)
		  /angDist(lens->zlens,lens->zsource,cosmo)/angDist(0,lens->zlens,cosmo)
		  /8.39428142e-10;

  printf(">critical density is %e Msun/Mpc^2    ro=%e Mpc  D_l = %e Mpc D_s = %e Mpc  to = %e days/Mpc^2\n"
		  ,lens->Sigma_crit,lens->host_ro
		  ,angDist(0,lens->zlens,cosmo),angDist(0,lens->zsource,cosmo),lens->to);

  NSubstructInRe = pi*pow(2*lens->host_ro,2)*lens->sub_Ndensity;

  if(NSubstructInRe > 0){

	  // include clumps that are beyond 2 Re
	  /*
	  lens->NSubstruct = (int)(NSubstructInRe*FractionWithinRe(lens,2) + 0.5);
	  Rmax = lens->ro*2 + lens->RmaxSubstruct
		          + pow(2*lens->MmaxSubstruct*lens->ro/pi/lens->Sigma_crit/1.0e-3,1./3.);
	  lens->NSubstruct = (int)(lens->NdensitySubstruct*pi*Rmax*Rmax+0.5);
	  */
	  printf(">average number of clumps including outside Re: %f\n",pi*pow(lens->host_ro,2)*lens->sub_Ndensity);
	  printf(">average clumps mass: %e Msun\n",averageSubMass(lens));

	  /*for(i=0;i<lens->NSubstruct;++i){
		  lens->RcutSubstruct[i]=lens->RmaxSubstruct;
		  lens->massSubstruct[i]=lens->MmaxSubstruct;
	  }*/

  }

  lens->set=true;
  printf(">\n");
}

/** \ingroup Constructor
 *
 */
void free_AnaLens(AnaLens *lens){
	free(lens->perturb_modes);
	if(lens->sub_N > 0 && lens->substruct_implanted){
	  free_dmatrix(lens->sub_x,0,lens->sub_N-1,0,1);
	  free(lens->sub_Rcut);
	  free(lens->sub_mass);
	  free(lens->perturb_rms);
	  free(lens->perturb_modes);
	}
	if(lens->stars_N > 0 && lens->stars_implanted){
		free(lens->star_masses);
		free(lens->stars);
		free_PosTypeMatrix(lens->stars_xp,0,lens->stars_N-1,0,2);
		free(lens->star_region);
		free(lens->star_kappa);
		free_dmatrix(lens->star_xdisk,0,lens->star_Nregions-1,0,1);
	}
}

/** \ingroup ImageFinding
 * \brief Prints the parameters of the analytic lens to stdout
 */
void PrintAnaLens(AnaLens *lens,bool show_substruct,bool show_stars){
	int i;

	printf(">Output file %s\n\n",lens->outputfile);
	// parameters of host elliptical
	printf(">  Host lens model\n");
	printf(">sigma %f km/s\n",lens->host_sigma);
	printf(">core %f Mpc\n",lens->host_core);
	printf(">axis_ratio %f\n",lens->host_axis_ratio);
	printf(">position angle %f\n",lens->host_pos_angle);
	printf(">r_source on lens plane %e pc\n",lens->source_r*1.0e6);

			// parameters of distortion to host elliptical
	printf("\n>Nmodes %li\n",lens->perturb_Nmodes);
	printf(">   beta=%e\n",lens->perturb_beta);
	if(lens->perturb_Nmodes>0){
		printf(">rms\n");
		for(i=0;i<6;++i) printf(">  %e\n",lens->perturb_rms[i]);
		printf(">modes\n");
		for(i=0;i<lens->perturb_Nmodes;++i) printf(">  %e\n",lens->perturb_modes[i]);
	}

	printf(">  Source\n");

	if(lens->source_sb_type == Uniform){
		printf(">  uniform surface brightness source\n");
	}else if(lens->source_sb_type == Gaussian){
		printf(">  Gaussian surface brightness source\n");
		printf(">     sigma^2 = %e Mpc^2",lens->source_gauss_r2);

	  }else{

		  printf(">  BLR surface brightness source\n");
		  switch(lens->source_sb_type){
		  printf(">      BH mass %e Msun\n",lens->source_BHmass);
		  printf(">      gamma %.3f \n",lens->source_gamma);
		  printf(">      inner radius %e pc\n",lens->source_r_in*1.0e6);
		  printf(">      outer radius %e pc\n",lens->source_r_out*1.0e6);
		  case BLR_Disk:
			  printf(">    disk model\n");
			  printf(">      inclination %.3f rads\n",lens->source_inclination);
			  printf(">      disk opening angle %.3f rads\n",lens->source_opening_angle);
			  printf(">      turbulent/thermal dispersion %e x V_kepler\n",lens->source_fK);
			  break;
		  case BLR_Sph1:
			  printf(">    spherical with circular orbits\n");
			  break;
		  case BLR_Sph2:
			  printf(">    spherical with Gaussian velocities\n");
			  printf(">      turbulent/thermal dispersion %e x V_kepler\n",lens->source_fK);
			  break;
		  default:
			  ERROR_MESSAGE();
			  printf(">ERROR: no submass internal profile chosen\n");
		  }


		  if(lens->source_monocrome) printf(">      monocromatic\n");
		  else printf(">      center of line %e Hz\n",lens->source_nuo);
	  }


	  // parameters of substructures
	printf(">  Substructures\n");
	printf(">NdensitySubstruct %e\n",lens->sub_Ndensity);
	printf(">NSubstruct %i\n",lens->sub_N);

	if(lens->sub_N > 0){
		printf(">betaSubstruct %.4f\n",lens->sub_beta);
		printf(">alphaSubstruct %.4f\n",lens->sub_alpha);
		printf(">RmaxSubstruct %.3e Mpc\n",lens->sub_Rmax);
		printf(">MmaxSubstruct %.3e Msun\n",lens->sub_Mmax);
		printf(">MminSubstruct %.3e Msun\n\n",lens->sub_Mmin);

		if(show_substruct){
			if(lens->substruct_implanted){
				for(i=0;i<lens->sub_N;++i){
					printf(">RcutSubstruct[%i] %.3e Mpc\n",i,lens->sub_Rcut[i]);
					printf(">massSubstruct[%i] %.3e Msun\n",i,lens->sub_mass[i]);
					printf(">xSubstruct[%i] %e %e Mpc\n\n",i,lens->sub_x[i][0],lens->sub_x[i][1]);
					switch(lens->sub_type){
					case NFW:
						printf(">  NFW clumps\n");
						break;
					case powerlaw:
						printf(">  Power Law clumps\n");
						break;
					case pointmass:
						printf(">  Point Mass clumps\n");
						break;
					default:
						ERROR_MESSAGE();
						printf(">ERROR: no submass internal profile chosen\n");
						exit(1);
					}
				}
			}else printf(">  substructures are implanted yet\n");
		}
	}
	if(lens->stars_N>0){
		printf("\nNstars=%li\n",lens->stars_N);
		printf(">stars_Nregions %i\n",lens->star_Nregions);
		printf(">stars_massscale %e\n",lens->star_massscale);
		printf(">stars_fstars %e\n",lens->star_fstars);
		printf(">stars_theta_force %e\n",lens->star_theta_force);
		if(show_stars){
			if(lens->stars_implanted){
				for(i=0 ; i < lens->stars_N ; ++i) printf(">    x[%li] = %e %e\n",i,lens->stars_xp[i][0],lens->stars_xp[i][1]);
			}else printf("> stars not implanted yet\n");
		}
	}

	// redshifts
	printf("\n>zlens %f\n",lens->zlens);
	printf(">zsource %f\n",lens->zsource);

	printf(">critical density is %e Msun/Mpc^2\n",lens->Sigma_crit);
}

void reNormSubstructure(AnaLens *lens,double kappa_sub){
	/* renomalizes substructure so that
	 * the average surface density it kappa_sub
	 */
	  double avem;
	  avem=lens->sub_Mmax*(lens->sub_alpha+1)
	  			  /(lens->sub_alpha+2)*(1-pow(lens->sub_Mmin/lens->sub_Mmax,lens->sub_alpha+2))/
	  			  (1-pow(lens->sub_Mmin/lens->sub_Mmax,lens->sub_alpha+1));

	  lens->sub_Ndensity=kappa_sub*lens->Sigma_crit/avem;

	  return ;
}
