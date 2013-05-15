/*
 * lens_halos.cpp
 *
 *  Created on: 06.05.2013
 *      Author: mpetkova
 */

#include "slsimlib.h"

LensHalo::LensHalo(){

}

LensHalo::LensHalo(InputParams& params){
	assignParams(params);
}

void LensHalo::error_message1(std::string parameter,std::string file){
		  ERROR_MESSAGE();
		  std::cout << "Parameter " << parameter << " is needed to construct a LensHalo.  It needs to be set in parameter file " << file << "!" << std::endl;
		  exit(0);
}

void LensHalo::assignParams(InputParams& params){
	if(!params.get("mass",mass)) error_message1("mass",params.filename());
	if(!params.get("Rmax",Rmax)) error_message1("Rmax",params.filename());
	if(!params.get("z_lens",zlens)) error_message1("z_lens",params.filename());
}

LensHalo::~LensHalo(){

}

NFWLensHalo::NFWLensHalo() : LensHalo(){
}

NFWLensHalo::NFWLensHalo(InputParams& params){
	assignParams(params);
}

void NFWLensHalo::assignParams(InputParams& params){
	if(!params.get("mass",mass)) error_message1("mass",params.filename());
	if(!params.get("Rmax",Rmax)) error_message1("Rmax",params.filename());
	if(!params.get("z_lens",zlens)) error_message1("z_lens",params.filename());
	if(!params.get("concentration",rscale)) error_message1("concentration",params.filename());
	rscale = rscale*Rmax; // was the concentration
}

NFWLensHalo::~NFWLensHalo(){

}

void NFWLensHalo::set_internal(long *seed, float vmax, float r_halfmass){

	NFW_Utility nfw_util;

	// Find the NFW profile with the same mass, Vmax and R_halfmass
	nfw_util.match_nfw(vmax,r_halfmass,mass,&rscale,&Rmax);

	rscale = Rmax/rscale; // Was the concentration
}

PseudoNFWLensHalo::PseudoNFWLensHalo() : LensHalo(){

}

PseudoNFWLensHalo::PseudoNFWLensHalo(InputParams& params){
	assignParams(params);
}

void PseudoNFWLensHalo::assignParams(InputParams& params){
	if(!params.get("mass",mass)) error_message1("mass",params.filename());
	if(!params.get("Rmax",Rmax)) error_message1("Rmax",params.filename());
	if(!params.get("z_lens",zlens)) error_message1("z_lens",params.filename());
	if(!params.get("concentration",rscale)) error_message1("concentration",params.filename());
	if(!params.get("slope_pnfw",beta)) error_message1("slope_pnfw",params.filename());
	rscale = rscale*Rmax; // was the concentration
}

PseudoNFWLensHalo::~PseudoNFWLensHalo(){

}



PowerLawLensHalo::PowerLawLensHalo() : LensHalo(){
	rscale = 1.0;
}

PowerLawLensHalo::PowerLawLensHalo(InputParams& params){
	assignParams(params);
}

void PowerLawLensHalo::assignParams(InputParams& params){
	if(!params.get("mass",mass)) error_message1("mass",params.filename());
	if(!params.get("Rmax",Rmax)) error_message1("Rmax",params.filename());
	if(!params.get("z_lens",zlens)) error_message1("z_lens",params.filename());
	if(!params.get("slope_pl",beta)) error_message1("slope_pl",params.filename());
}

PowerLawLensHalo::~PowerLawLensHalo(){

}


BaseNSIELensHalo::BaseNSIELensHalo() : LensHalo(){
	rscale = 1.0;
}

BaseNSIELensHalo::BaseNSIELensHalo(InputParams& params){
	assignParams(params);
}

void BaseNSIELensHalo::assignParams(InputParams& params){
	if(!params.get("mass",mass)) error_message1("mass",params.filename());
	if(!params.get("Rmax",Rmax)) error_message1("Rmax",params.filename());
	if(!params.get("z_lens",zlens)) error_message1("z_lens",params.filename());

	if(!params.get("sigma",sigma)) error_message1("sigma",params.filename());
	if(!params.get("core",rcore)) error_message1("core",params.filename());
	if(!params.get("axis_ratio",fratio)) error_message1("axis_ratio",params.filename());
	if(!params.get("pos_angle",pa)) error_message1("pos_angle",params.filename());

	Rsize = rmaxNSIE(sigma,mass,fratio,rcore);
	Rmax = MAX(1.0,1.0/fratio)*Rsize;  // redefine

	assert(Rmax >= Rsize);
}

BaseNSIELensHalo::~BaseNSIELensHalo(){

}
void BaseNSIELensHalo::set_internal(long *seed, float vmax, float r_halfmass){
	rcore = 0.0;
	sigma = 126*pow(mass/1.0e10,0.25); // From Tully-Fisher and Bell & de Jong 2001
	fratio = (ran2(seed)+1)*0.5;  //TODO This is a kluge.
	pa = 2*pi*ran2(seed);  //TODO This is a kluge.
	Rsize = rmaxNSIE(sigma,mass,fratio,rcore);
	Rmax = MAX(1.0,1.0/fratio)*Rsize;  // redefine

	assert(Rmax >= Rsize);
}

NSIELensHalo::NSIELensHalo(InputParams& params) : BaseNSIELensHalo(){
	perturb_rms = new double[6];
	substruct_implanted = false;
	stars_implanted = false;
	sub_Ndensity = 0.0;
	perturb_Nmodes = 0;

	assignParams(params);

	if(perturb_Nmodes){
		perturb_modes = new double[perturb_Nmodes+1];
		for(int i=0;i< perturb_Nmodes+1 ;++i) perturb_modes[i] =  0;
	}
}

NSIELensHalo::~NSIELensHalo(){
	std::cout << "deleting NSIELensHalo" << std::endl;

	delete[] perturb_rms;

	if(perturb_Nmodes > 0){
		std::cout << "deleting modes" << std::endl;
		delete[] perturb_modes;
	}
	if(sub_N > 0 && substruct_implanted){
		std::cout << "deleting subs" << std::endl;
		Utilities::free_PosTypeMatrix(sub_x,sub_N,2);
		delete[] subs;
		delete[] sub_substructures;
	}
	if(stars_N > 0 && stars_implanted){
		std::cout << "deleting stars" << std::endl;
		delete[] star_masses;
		delete[] stars;
		Utilities::free_PosTypeMatrix(stars_xp,stars_N,3);
		delete[] star_region;
		delete[] star_kappa;
		Utilities::free_PosTypeMatrix(star_xdisk,star_Nregions,2);
		delete star_tree;
	}
}

void NSIELensHalo::setInternalParams(CosmoHndl cosmo, double zsource){
	double Ds, Dls, Dl;

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

void NSIELensHalo::error_message1(std::string parameter,std::string file){
		  ERROR_MESSAGE();
		  std::cout << "Parameter " << parameter << " is needed to construct a NSIELensHalo.  It needs to be set in parameter file " << file << "!" << std::endl;
		  exit(0);
}

void NSIELensHalo::force_halo(
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


void NSIELensHalo::assignParams(InputParams& params){
	//  lens parameters
	if(!params.get("sigma",sigma)) error_message1("sigma",params.filename());
	if(!params.get("core",rcore)) error_message1("core",params.filename());
	if(!params.get("axis_ratio",fratio)) error_message1("axis_ratio",params.filename());
	if(!params.get("pos_angle",pa)) error_message1("pos_angle",params.filename());
	if(!params.get("z_lens",zlens)) error_message1("z_lens",params.filename());

	// Distortion of host lens parameters
	if(!params.get("NDistortionModes",perturb_Nmodes)) error_message1("NDistortionModes",params.filename());
	else if(perturb_Nmodes > 0){
		if(!params.get("beta_perturb",perturb_beta)) error_message1("beta_perturb",params.filename());
		else if(perturb_beta <= 0.0) {ERROR_MESSAGE(); std::cout << "perturb_beta can't be <= 0.0 in file " << params.filename() << std::endl;}
		if(!params.get("kappa_perturb",perturb_rms[0])) error_message1("kappa_perturb",params.filename());
		if(!params.get("gamma_perturb",perturb_rms[1])) error_message1("gamma_perturb",params.filename());
		if(!params.get("monopole_perturb",perturb_rms[2])) error_message1("monopole_perturb",params.filename());
		if(!params.get("quadrapole_perturb",perturb_rms[3])) error_message1("quadrapole_perturb",params.filename());
		if(!params.get("hexopole_perturb",perturb_rms[4])) error_message1("hexopole_perturb",params.filename());
		if(!params.get("octopole_perturb",perturb_rms[5])) error_message1("octopole_perturb",params.filename());
	}

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
    		std::cout << "fstars cannot be less than 0 or larger than 1 in file " << params.filename() << std::endl;
    		exit(0);
    	}
    	if(!params.get("stars_mass",star_massscale)) error_message1("stars_mass",params.filename());
    }
}

AnaNSIELensHalo::AnaNSIELensHalo(InputParams& params) : NSIELensHalo(params){

  assignParams(params);

  if(perturb_Nmodes){
  	perturb_modes = new double[perturb_Nmodes+1];
  	// zero perturbation modes until use BaseAnaLens::RandomlyDistortLens()
  	for(int i=0;i< perturb_Nmodes+1 ;++i) perturb_modes[i] =  0;
  }

  // in degrees
  pa*=pi/180;
}


AnaNSIELensHalo::~AnaNSIELensHalo(){

}


void AnaNSIELensHalo::assignParams(InputParams& params){

	// Host lens parameters
	if(!params.get("sigma",sigma)) error_message1("sigma",params.filename());
	if(!params.get("core",rcore)) error_message1("core",params.filename());
	if(!params.get("axis_ratio",fratio)) error_message1("axis_ratio",params.filename());
	if(!params.get("pos_angle",pa)) error_message1("pos_angle",params.filename());
	if(!params.get("z_lens",zlens)) error_message1("z_lens",params.filename());


	// Distortion of host lens parameters
	if(!params.get("NDistortionModes",perturb_Nmodes)) error_message1("NDistortionModes",params.filename());
	else if(perturb_Nmodes > 0){
		if(!params.get("beta_perturb",perturb_beta)) error_message1("beta_perturb",params.filename());
		else if(perturb_beta <= 0.0) {ERROR_MESSAGE(); std::cout << "perturb_beta can't be <= 0.0 in file " << params.filename() << std::endl; }
		if(!params.get("kappa_perturb",perturb_rms[0])) error_message1("kappa_perturb",params.filename());
		if(!params.get("gamma_perturb",perturb_rms[1])) error_message1("gamma_perturb",params.filename());
		if(!params.get("monopole_perturb",perturb_rms[2])) error_message1("monopole_perturb",params.filename());
		if(!params.get("quadrapole_perturb",perturb_rms[3])) error_message1("quadrapole_perturb",params.filename());
		if(!params.get("hexopole_perturb",perturb_rms[4])) error_message1("hexopole_perturb",params.filename());
		if(!params.get("octopole_perturb",perturb_rms[5])) error_message1("octopole_perturb",params.filename());
	}

}


UniNSIELensHalo::UniNSIELensHalo(InputParams& params) : NSIELensHalo(params){

  assignParams(params);

  perturb_Nmodes=3;
  perturb_modes = new double[3];

  perturb_modes[0]=kappa_uniform;
  perturb_modes[1]=gamma_uniform[0];
  perturb_modes[2]=gamma_uniform[1];
}

UniNSIELensHalo::~UniNSIELensHalo(){

}

void UniNSIELensHalo::assignParams(InputParams& params){

	if(!params.get("kappa_uniform",kappa_uniform)) error_message1("kappa_uniform",params.filename());
	if(!params.get("gamma_uniform_1",gamma_uniform[0])) error_message1("gamma_uniform_1",params.filename());
	if(!params.get("gamma_uniform_2",gamma_uniform[1])) error_message1("gamma_uniform_2",params.filename());
	if(!params.get("z_lens",zlens)) error_message1("z_lens",params.filename());

	return;
}


void LensHalo::force_halo(
		double *alpha     /// mass/Mpc
		,KappaType *kappa
		,KappaType *gamma
		,double *xcm
		,bool no_kappa
		){

	double rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
	if(rcm2 < 1e-20) rcm2 = 1e-20;

	/// intersecting, subtract the point particle
	if(rcm2 < Rmax*Rmax){
		double prefac = mass/rcm2/pi;
		double arg1 = sqrt(rcm2)/rscale;
		double arg2 = Rmax/rscale;

		double tmp = (alpha_h(arg1,arg2) + 1.0)*prefac;
		alpha[0] += tmp*xcm[0];
		alpha[1] += tmp*xcm[1];

		// can turn off kappa and gamma calculations to save times
		if(!no_kappa){
			*kappa += kappa_h(arg1,arg2)*prefac;

			tmp = (gamma_h(arg1,arg2) + 2.0)*prefac/rcm2;

			gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*tmp;
			gamma[1] += xcm[0]*xcm[1]*tmp;
		}
	}

	return;
}

void BaseNSIELensHalo::force_halo(
		double *alpha
		,KappaType *kappa
		,KappaType *gamma
		,double *xcm
		,bool no_kappa
){

	double rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
	if(rcm2 < 1e-20) rcm2 = 1e-20;

	double ellipR = ellipticRadiusNSIE(xcm,fratio,pa);
	if(ellipR > Rsize){
		double rout = Rsize*MAX(1.0,1.0/fratio);
		// This is the case when the ray is within the NSIE's circular region of influence but outside its elliptical truncation
		// if the ray misses the halo treat it as a point mass
		double prefac = -1.0*mass/rcm2/pi;

		if(rcm2 > rout*rout){
			alpha[0] += prefac*xcm[0];
			alpha[1] += prefac*xcm[1];
		}else{
			double alpha_out[2],alpha_in[2],rin,x_in[2];
			double prefac = -1.0*mass/rout/pi;
			double r = sqrt(rcm2);

			alpha_out[0] = prefac*xcm[0]/r;
			alpha_out[1] = prefac*xcm[1]/r;

			Utilities::rotation(x_in,xcm,pa);
			rin = r*Rsize
					/sqrt( x_in[0]*x_in[0] + pow(fratio*x_in[1],2) );
			//rin = Rsize;

			x_in[0] = rin*xcm[0]/r;
			x_in[1] = rin*xcm[1]/r;

			alpha_in[0] = alpha_in[1] = 0;
			float units = pow(sigma/lightspeed,2)/Grav/sqrt(fratio); // mass/distance(physical)
			alphaNSIE(alpha_in,x_in,fratio,rcore,pa);
			alpha_in[0] *= -units;
			alpha_in[1] *= -units;

			alpha[0] += (r - rin)*(alpha_out[0] - alpha_in[0])/(rout - rin) + alpha_in[0];
			alpha[1] += (r - rin)*(alpha_out[1] - alpha_in[1])/(rout - rin) + alpha_in[1];
		}

		// can turn off kappa and gamma calculations to save times
		if(!no_kappa){
			prefac *= 2.0/rcm2;

			gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*prefac;
			gamma[1] += xcm[0]*xcm[1]*prefac;
		}

	}else{
		double xt[2]={0,0},tmp[2]={0,0};
		float units = pow(sigma/lightspeed,2)/Grav/sqrt(fratio); // mass/distance(physical)
		xt[0]=xcm[0];
		xt[1]=xcm[1];
		alphaNSIE(tmp,xt,fratio,rcore,pa);
		alpha[0] -= units*tmp[0];
		alpha[1] -= units*tmp[1];
		if(!no_kappa){
			KappaType tmp[2]={0,0};
			*kappa += units*kappaNSIE(xt,fratio,rcore,pa);
			gammaNSIE(tmp,xt,fratio,rcore,pa);
			gamma[0] += units*tmp[0];
			gamma[1] += units*tmp[1];
		}
	}
	return;
}
