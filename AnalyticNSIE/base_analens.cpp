/*
 * base_analens.cpp
 *
 *  Created on: Jan 22, 2013
 *      Author: mpetkova
 */

#include "slsimlib.h"

using namespace std;

void BaseNSIELensHalo::force_halo(
		double *alpha     /// mass/Mpc
		,KappaType *kappa
		,KappaType *gamma
		,double *xcm
		,bool no_kappa
		){
     double x_rescale[2];
     long j;
     double alpha_tmp[2];
     KappaType kappa_tmp = 0.0, gamma_tmp[3], dt = 0,tmp = 0;

     gamma_tmp[0] = gamma_tmp[1] = gamma_tmp[2] = 0.0;
     alpha_tmp[0] = alpha_tmp[1] = 0.0;

     alpha[0] = alpha[1] = 0.0;
     gamma[0] = gamma[1] = gamma[2] = 0.0;
     *kappa = 0.0;

     double convert_factor = star_massscale/Sigma_crit;

     if(Einstein_ro > 0){
    	 x_rescale[0] = xcm[0]/Einstein_ro;
    	 x_rescale[1] = xcm[1]/Einstein_ro;

    	 alphaNSIE(alpha,x_rescale,fratio,rcore/Einstein_ro,pa);

    	 if(!no_kappa){
    		 gammaNSIE(gamma,x_rescale,fratio,rcore/Einstein_ro,pa);
    		 *kappa=kappaNSIE(x_rescale,fratio,rcore/Einstein_ro,pa);
    	 }

    	 alpha[0] *= Einstein_ro;
    	 alpha[1] *= Einstein_ro;
     }

  // perturbations of host lens
     if(perturb_Nmodes > 0){
    	 *kappa += lens_expand(perturb_beta,perturb_modes
    			 ,perturb_Nmodes,xcm,alpha_tmp,gamma_tmp,&dt);

    	 alpha[0] += alpha_tmp[0];
    	 alpha[1] += alpha_tmp[1];

   	      if(!no_kappa){
   	    	  gamma[0] += gamma_tmp[0];
   	    	  gamma[1] += gamma_tmp[1];
   	      }

   	     gamma_tmp[0] = gamma_tmp[1] = gamma_tmp[2] = 0.0;
   	     alpha_tmp[0] = alpha_tmp[1] = 0.0;
     }

     // add substructure
     if(substruct_implanted){
    	 for(j=0;j<sub_N;++j){

    		 subs[j].force_halo(alpha_tmp,&kappa_tmp,gamma_tmp,xcm,no_kappa);

    		 alpha[0] += alpha_tmp[0];
    		 alpha[1] += alpha_tmp[1];

    		 if(!no_kappa){
    			 *kappa += kappa_tmp;
    			 gamma[0] += gamma_tmp[0];
    			 gamma[1] += gamma_tmp[1];
    		 }
    	 }

         gamma_tmp[0] = gamma_tmp[1] = gamma_tmp[2] = 0.0;
         alpha_tmp[0] = alpha_tmp[1] = 0.0;
     }

     // add stars for microlensing
     if(stars_N > 0 && stars_implanted){

    	 substract_stars_disks(xcm,alpha,kappa,gamma);

    	 // do stars with tree code
    	 star_tree->force2D_recur(xcm,alpha_tmp,&tmp,gamma_tmp,no_kappa);

    	 alpha[0] -= convert_factor*alpha_tmp[0];
    	 alpha[1] -= convert_factor*alpha_tmp[1];

    	 if(!no_kappa){
    		 *kappa += convert_factor*tmp;
    		 gamma[0] += convert_factor*gamma_tmp[0];
    		 gamma[1] += convert_factor*gamma_tmp[1];
    	 }
     }

     // final operations on results
     convert_factor = 4*pi*Grav*Sigma_crit;

     // convert from physical distance on the lens plane to an angle
	 alpha[0] *= convert_factor;
	 alpha[1] *= convert_factor;

	 // in the multi-plane formalism G^i=partial deflection_angle^i / partial x^i
	 // therefore the quantities need to be in units (1/physical_distance)
	 // --> convert from unitless quantity to (1/physical_distance)
	 *kappa *= convert_factor;
	 gamma[0] *= convert_factor;
	 gamma[1] *= convert_factor;
	 gamma[2] *= convert_factor;

     return ;
}
/**
 * \brief Reads in a parameter file and sets up an analytic lens.
 *
 * Sets many parameters within the lens model, source model and
 * force calculation.
 */
void BaseNSIELensHalo::assignParams(InputParams& params){

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

void BaseNSIELensHalo::error_message1(std::string parameter,std::string file){
		  ERROR_MESSAGE();
		  std::cout << "Parameter " << parameter << " is needed to construct a BaseAnaLens.  It needs to be set in parameter file " << file << "!" << endl;
		  exit(0);
}

/// resets Zl, Dl, Sigma_crit, MpcToAsec
void BaseNSIELensHalo::setZlens(CosmoHndl cosmo,double zl,double zsource){
	zlens = zl;
	setInternalParams(cosmo, zsource);
}

void BaseNSIELensHalo::reNormSubstructure(double kappa_sub){
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
void BaseNSIELensHalo::setInternalParams(CosmoHndl cosmo, SourceHndl source){
	setInternalParams(cosmo,source->getZ());
}

void BaseNSIELensHalo::setInternalParams(CosmoHndl cosmo, double zsource){
	double Ds, Dls;

	if(zsource < zlens) zsource = 1000;
	Dl = cosmo->angDist(0,zlens);
	Ds = cosmo->angDist(0,zsource);
	Dls = cosmo->angDist(zlens,zsource);

	MpcToAsec = 60*60*180 / pi / Dl;
		// in Mpc
	Einstein_ro=4*pi*pow(sigma/2.99792e5,2)*Dl
		*Dls/Ds;
	// find critical density
	Sigma_crit=Ds/Dls/Dl/4/pi/Grav;
	to = (1+zlens)*Ds/Dls/Dl/8.39428142e-10;
}

BaseNSIELensHalo::BaseNSIELensHalo(InputParams& params) : SimpleNSIELensHalo(){

  perturb_rms = new double[6];

  assignParams(params);

  // parameters for stars
  stars_implanted = false; // stars are implanted later
  star_theta_force = 0.1;
  sub_theta_force = 0.1;

  perturb_Nmodes = 0;
  sub_sigmaScale = sigma = pa = Einstein_ro = fratio = rcore = 0.0;

  if(sub_Ndensity == 0)
	  sub_N = 0;

  Sigma_crit = 0;

  stars_implanted = false;
}


void BaseNSIELensHalo::PrintLens(bool show_substruct,bool show_stars){
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
				  cout << "RcutSubstruct "<<i << " " <<subs[i].get_Rmax() << " Mpc" << endl;
				  cout << "massSubstruct "<<i<<" "<<subs[i].get_mass() << " Msun" << endl;
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


BaseNSIELensHalo::~BaseNSIELensHalo(){
	cout << "deleting lens" << endl;

	delete[] perturb_rms;

	if(perturb_Nmodes > 0){
		cout << "deleting modes" << endl;
		delete[] perturb_modes;
	}
	if(sub_N > 0 && substruct_implanted){
		cout << "deleting subs" << endl;
		Utilities::free_PosTypeMatrix(sub_x,sub_N,2);
		delete[] subs;
		delete[] sub_substructures;
	}
	if(stars_N > 0 && stars_implanted){
		cout << "deleting stars" << endl;
		delete[] star_masses;
		delete[] stars;
		Utilities::free_PosTypeMatrix(stars_xp,stars_N,3);
		delete[] star_region;
		delete[] star_kappa;
		Utilities::free_PosTypeMatrix(star_xdisk,star_Nregions,2);
		delete star_tree;
	}
}

/******************************************************
 Below are routines for calculating the deflection etc. 
 for asymetric halos
 ******************************************************/


/// TODO This needs to be thought about more. sets axial modes to reproduce a near elliptically shaped surface density
void LensHalo::setModesToEllip(double q,double theta){
  // elliptical integrals
	double K = rfD(0,1./q/q,1);
	double E = K - (1-1./q/q)*rdD(0,1./q/q,1)/3;
  assert(Nmod == 18);
  
  // set modo to elliptical model
	for(int i=1;i<=Nmod;++i){
		mod[i]=0;
	}
	// fill in modes with their values for an elliptical lens
	mod[3]=4*K/pi;
	if(q != 1.0){
		mod[4] = 4*( (1+q*q)*K-2*q*q*E )/(1-q*q)/pi/(1-4);
		mod[8] = 4*( (3*q*q+1)*(q*q+3)*K-8*q*q*(1+q*q)*E )
      /( 3*pi*pow(1-q*q,2) )/(1-16);
		mod[12] = 4*( (1+q*q)*(15+98*q*q+15*q*q*q*q)*K-2*q*q*(23+82*q*q+23*q*q*q*q)*E )
      /( 15*pi*pow(1-q*q,3) )/(1-36);
		mod[16]= 4*( -32*q*q*(1+q*q)*(11+74*q*q+11*q*q*q*q)*E
                           +(105+1436*q*q+3062*q*q*q*q+1436*pow(q,6)+105*pow(q,8))*K )
      /(105*pi*pow(1-q*q,4))/(1-64);
	}
  
	// rotate model
	RotateModel(theta,mod,Nmod,0);
  
  return;
}

/// Derivatives of the axial potential factor with respect to theta
void LensHalo::faxial(double theta,double f[]){
  int i,k;
  
  //f[0] = 0.5*mod[3];
  f[0] = 1;
  f[1] = f[2] = 0;
  for(i=4;i<Nmod;i+=2){
    k=i/2;
    f[0] +=  mod[i]*cos(k*theta)     + mod[i+1]*sin(k*theta);
    f[1] += -mod[i]*k*sin(k*theta)   + mod[i+1]*k*cos(k*theta);
    f[2] += -mod[i]*k*k*cos(k*theta) - mod[i+1]*k*k*sin(k*theta);
  }

}

