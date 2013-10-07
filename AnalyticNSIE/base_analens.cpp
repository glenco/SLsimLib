/*
 * base_analens.cpp
 *
 *  Created on: Jan 22, 2013
 *      Author: mpetkova
 */

#include "slsimlib.h"

using namespace std;

void LensHaloBaseNSIE::force_halo(
		double *alpha       /// mass/Mpc
		,KappaType *kappa   /// surface mass density
		,KappaType *gamma
		,double *xcm
		,bool no_kappa
		,bool subtract_point /// if true contribution from a point mass is subtracted
		){
     long j;
     double alpha_tmp[2];
     KappaType kappa_tmp = 0.0, gamma_tmp[3], dt = 0;

     gamma_tmp[0] = gamma_tmp[1] = gamma_tmp[2] = 0.0;
     alpha_tmp[0] = alpha_tmp[1] = 0.0;

     alpha[0] = alpha[1] = 0.0;
     gamma[0] = gamma[1] = gamma[2] = 0.0;
     *kappa = 0.0;

	 double xt[2]={0,0};
	 float units = pow(sigma/lightspeed,2)/Grav;///sqrt(fratio); // mass/distance(physical)
	 xt[0]=xcm[0];
	 xt[1]=xcm[1];
     alphaNSIE(alpha,xt,fratio,rcore,pa);
	 alpha[0] *= units;
	 alpha[1] *= units;

	 if(!no_kappa){
    	gammaNSIE(gamma,xcm,fratio,rcore,pa);
    	*kappa=kappaNSIE(xcm,fratio,rcore,pa);
    	*kappa *= units;
    	gamma[0] *= units;
    	gamma[1] *= units;
		gamma[2] *= units;
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
    	 force_stars(alpha,kappa,gamma,xcm,no_kappa);
     }
     return ;
}
/**
 * \brief Reads in a parameter file and sets up an analytic lens.
 *
 * Sets many parameters within the lens model, source model and
 * force calculation.
 */
void LensHaloBaseNSIE::assignParams(InputParams& params){
	if(!params.get("main_mass",mass)) error_message1("main_mass",params.filename());
	if(!params.get("main_zlens",zlens)) error_message1("main_zlens",params.filename());

	if(!params.get("main_sigma",sigma)) error_message1("main_sigma",params.filename());
	if(!params.get("main_core",rcore)) error_message1("main_core",params.filename());
	if(!params.get("main_axis_ratio",fratio)) error_message1("main_axis_ratio",params.filename());
	if(!params.get("main_pos_angle",pa)) error_message1("main_pos_angle",params.filename());

	Rsize = rmaxNSIE(sigma,mass,fratio,rcore);
	Rmax = MAX(1.0,1.0/fratio)*Rsize;  // redefine

	assert(Rmax >= Rsize);

	if(!params.get("zsource_reference",zsource_reference)) error_message1("zsource_reference",params.filename());

    // Substructure parameters
    if(!params.get("main_sub_Ndensity",sub_Ndensity)) error_message1("main_sub_Ndensity",params.filename());

    else if(sub_Ndensity > 0){
    	if(!params.get("main_sub_beta",sub_beta)) error_message1("main_sub_beta",params.filename());
    	if(!params.get("main_sub_alpha",sub_alpha)) error_message1("main_sub_alpha",params.filename());
    	if(!params.get("main_sub_Rmax",sub_Rmax)) error_message1("main_sub_Rmax",params.filename());
    	if(!params.get("main_sub_mass_max",sub_Mmax)) error_message1("main_sub_mass_max",params.filename());
    	if(!params.get("main_sub_mass_min",sub_Mmin)) error_message1("main_sub_mass_min",params.filename());
    	if(sub_Mmin < 1.0e3){
    		ERROR_MESSAGE();
    		std::cout << "Are you sure the minimum halo mass should be " << sub_Mmin << " Msun?" << std::endl;
    		exit(1);
    	}
    	if(!params.get("main_sub_type",main_sub_type)) error_message1("main_sub_type",params.filename());

    }
	  // Stars parameters
    if(!params.get("main_stars_N",stars_N)) error_message1("main_stars_N",params.filename());

    else if(stars_N){
    	assignParams_stars(params);
    }

}

void LensHaloBaseNSIE::error_message1(std::string parameter,std::string file){
		  ERROR_MESSAGE();
		  std::cout << "Parameter " << parameter << " is needed to construct a BaseAnaLens.  It needs to be set in parameter file " << file << "!" << endl;
		  exit(0);
}

/// resets Zl, Dl, Sigma_crit, MpcToAsec
void LensHaloBaseNSIE::setZlens(double zl){
	zlens = zl;
//	setInternalParams(cosmo, zsource);
}

void LensHaloBaseNSIE::reNormSubstructure(double kappa_sub){
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
void LensHaloAnaNSIE::setCosmology(const COSMOLOGY& cosmo)
{
	Dl = cosmo.angDist(0,zlens);
	Ds = cosmo.angDist(0,zsource_reference);
	Dls = cosmo.angDist(zlens,zsource_reference);
	MpcToAsec = 60*60*180 / pi / Dl;
		// in Mpc
	Einstein_ro=4*pi*pow(sigma/lightspeed,2)*Dl
		*Dls/Ds;
	// find critical density
	Sigma_crit=Ds/Dls/Dl/4/pi/Grav;
	to = (1+zlens)*Ds/Dls/Dl/8.39428142e-10;
}


LensHaloBaseNSIE::LensHaloBaseNSIE(InputParams& params) : LensHalo(){

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

  substruct_implanted = false;

}


void LensHaloBaseNSIE::PrintLens(bool show_substruct,bool show_stars){
	int i;
	cout << "zlens " << zlens << endl;

	 // parameters of substructures
	cout << endl << "main_sub_Ndensity "<< sub_Ndensity << endl;
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
					switch(main_sub_type){
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

	if(Sigma_crit)
		cout << "critical density is " << Sigma_crit << " Msun/Mpc^2" << endl << endl;

	if (stars_implanted) PrintStars(show_stars);
}

std::size_t LensHaloBaseNSIE::Nrandomize() const
{
	return LensHalo::Nrandomize() + 3;
}

Utilities::Any LensHaloBaseNSIE::randomize(std::size_t i, double step, long* seed)
{
	Utilities::Any old;
	
	if(i < LensHalo::Nrandomize())
	{
		old = LensHalo::randomize(i, step, seed);
	}
	else
	{
		switch(i - LensHalo::Nrandomize())
		{
			case 0:
				old = sigma;
				sigma += step*sigma*gasdev(seed);
				break;
			case 1:
				old = fratio;
				fratio += step*gasdev(seed);
				break;
			case 2:
				old = pa;
				pa += step*pi*gasdev(seed);
				break;
			default:
				throw std::invalid_argument("bad parameter index for randomize()");
		}
	}
	
	return old;
}

void LensHaloBaseNSIE::unrandomize(std::size_t i, const Utilities::Any& old)
{
	if(i < LensHalo::Nrandomize())
	{
		LensHalo::unrandomize(i, old);
	}
	else
	{
		switch(i - LensHalo::Nrandomize())
		{
			case 0:
				sigma = Utilities::AnyCast<double>(old);
				break;
			case 1:
				fratio = Utilities::AnyCast<double>(old);
				break;
			case 2:
				pa = Utilities::AnyCast<double>(old);
				break;
			default:
				throw std::invalid_argument("bad parameter index for unrandomize()");
		}
	}
}

void LensHaloBaseNSIE::printCSV(std::ostream& out, bool header) const
{
	if(header)
	{
		out
		<< "sigma" << ","
		<< "fratio" << ","
		<< "PA" << std::endl;
	}
	else
	{
		out
		<< sigma << ","
		<< fratio << ","
		<< pa << std::endl;
	}
}

LensHaloBaseNSIE::~LensHaloBaseNSIE(){
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
		delete[] star_Sigma;
		Utilities::free_PosTypeMatrix(star_xdisk,star_Nregions,2);
		delete star_tree;
	}
}

/******************************************************
 Below are routines for calculating the deflection etc. 
 for asymetric halos
 ******************************************************/


/// TODO: This needs to be thought about more. sets axial modes to reproduce a near elliptically shaped surface density
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
	if(q != 1.0){
    mod[3]=4*K/pi;
		mod[4] = 4*( (1+q*q)*K-2*q*q*E )/(1-q*q)/pi/mod[3];
		mod[8] = 4*( (3*q*q+1)*(q*q+3)*K-8*q*q*(1+q*q)*E )
      /( 3*pi*pow(1-q*q,2) )/mod[3];
		mod[12] = 4*( (1+q*q)*(15+98*q*q+15*q*q*q*q)*K-2*q*q*(23+82*q*q+23*q*q*q*q)*E )
      /( 15*pi*pow(1-q*q,3) )/mod[3];
		mod[16]= 4*( -32*q*q*(1+q*q)*(11+74*q*q+11*q*q*q*q)*E
                           +(105+1436*q*q+3062*q*q*q*q+1436*pow(q,6)+105*pow(q,8))*K )
      /(105*pi*pow(1-q*q,4))/mod[3];
	}
  mod[3]=1.0;
  
	// rotate model
	RotateModel(theta,mod,Nmod,0);
  
  return;
}

/// Derivatives of the axial potential factor with respect to theta
void LensHalo::faxial(double theta,double f[]){
  int i,k;
  
  //f[0] = 0.5*mod[3];
  f[0] = f[1] = f[2] = 0;
  for(i=4;i<Nmod;i+=2){
    k=i/2;
    f[0] +=  mod[i]*cos(k*theta)     + mod[i+1]*sin(k*theta);
    f[1] += -mod[i]*k*sin(k*theta)   + mod[i+1]*k*cos(k*theta);
    f[2] += -mod[i]*k*k*cos(k*theta) - mod[i+1]*k*k*sin(k*theta);
  }

}

/// Derivatives of the axial potential factor with respect to theta
void LensHalo::gradial(double r,double g[]){
  double x = (1+r/r_eps);
  
  g[0] = 1.0/x/x;
  g[1] = -2.0*g[0]/x/r_eps;
  g[2] = -3.0*g[1]/x/r_eps;
}

/** \brief This function returns the lensing quantities for an asymmetric version of the symmetric baseclass halo.
 *  
 *  This function should only be used by the second generation of classes derived from LensHalo.
 *
 *  The math needs to be double checked and the sign convention checked.  
 *  The method used to make the lenses asymmetric is laid out in http://metcalf1.bo.astro.it/~bmetcalf/ExtraNotes/notes_elliptical.pdf
 */
void LensHalo::desymmeterize(double r,double theta,double *alpha,double *kappa,double *gamma){
  double f[3],g[3];
  
  double alpha_iso = alpha_h(r/rscale),phi_iso = phi_h(r/rscale)
  ,kappa_iso = kappa_h(r/rscale),gamma_iso = gamma_h(r/rscale);
  
  double alpha_r,alpha_theta,F;

  faxial(theta,f);
  gradial(r,g);
  
  F = (1+g[0]*f[0]);
  
  alpha_r = (F + g[1]*f[0])*alpha_iso;
  alpha_theta = g[0]*f[1]*phi_iso/r;
  
  alpha[0] = alpha_r*cos(theta) - alpha_theta*sin(theta);
  alpha[1] = alpha_r*sin(theta) + alpha_theta*cos(theta);
  
  *kappa = F*kappa_iso + ( (g[2] + g[1]/r)*f[0] + g[0]*f[2]/r/r )*phi_iso;
  
  double gt = F*gamma_iso + g[1]*f[0]*alpha_iso + 0.5*( g[2]*f[0] - g[0]*f[2]/r/r)*phi_iso;
  double g45 = f[1]*(alpha_iso*g[0]/r + (g[1]-g[0]/r/r)*phi_iso);
  
  gamma[0] = cos(2*theta)*gt + sin(2*theta)*g45;
  gamma[1] = -sin(2*theta)*gt + cos(2*theta)*g45;
  
}

