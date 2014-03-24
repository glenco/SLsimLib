/*
 * base_analens.cpp
 *
 *  Created on: Jan 22, 2013
 *      Author: mpetkova
 */

#include "slsimlib.h"

using namespace std;

void LensHaloBaseNSIE::force_halo(
		PosType *alpha       /// mass/Mpc
		,KappaType *kappa   /// surface mass density
		,KappaType *gamma
        ,KappaType *phi     // PHI BY Fabien
		,PosType const *xcm
		,bool no_kappa
		,bool subtract_point /// if true contribution from a point mass is subtracted
		)
{
     long j;
     PosType alpha_tmp[2];
     KappaType kappa_tmp = 0.0, gamma_tmp[3], dt = 0;
     KappaType phi_tmp = 0.0 ; // PHI BY Fabien
            
            
     gamma_tmp[0] = gamma_tmp[1] = gamma_tmp[2] = 0.0;
     alpha_tmp[0] = alpha_tmp[1] = 0.0;

     alpha[0] = alpha[1] = 0.0;
     gamma[0] = gamma[1] = gamma[2] = 0.0;
     *kappa = 0.0;
     // PHI BY Fabien
     *phi = 0.0 ;

            
	 PosType xt[2]={0,0};
	 float units = pow(sigma/lightspeed,2)/Grav; ///sqrt(fratio); // mass/distance(physical)
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
         
        // PHI BY Fabien
        // *phi = phiNSIE(xcm,fratio,rcore,pa);
        // *phi *= units ; // Fabien : is this necessary for the potential ?
	 }

     // perturbations of host lens
     if(perturb_Nmodes > 0)
     {
    	 *kappa += lens_expand(perturb_beta,perturb_modes
    			 ,perturb_Nmodes,xcm,alpha_tmp,gamma_tmp,&dt);

        // PHI BY Fabien : should I put the computation of the potential somewhere here ?
         
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
     if(substruct_implanted)
     {
    	 for(j=0;j<sub_N;++j)
         {
             
             // PHI BY Fabien
    		 subs[j].force_halo(alpha_tmp,&kappa_tmp,gamma_tmp,&phi_tmp,xcm,no_kappa);
             // subs[j].force_halo(alpha_tmp,&kappa_tmp,gamma_tmp,xcm,no_kappa);

    		 alpha[0] += alpha_tmp[0];
    		 alpha[1] += alpha_tmp[1];

    		 if(!no_kappa)
             {
    			 *kappa += kappa_tmp;
    			 gamma[0] += gamma_tmp[0];
    			 gamma[1] += gamma_tmp[1];
                 
                 // PHY BY Fabien : add something for potential here ?
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
  else if(fratio > 1){
    ERROR_MESSAGE();
    std::cout << "parameter main_axis_ratio must be < 1 in file " << params.filename() << ". Use main_pos_angle to rotate the halo." << std::endl;
    exit(1);
  }
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

/*/ resets Zl, Dl, Sigma_crit, MpcToAsec
void LensHaloBaseNSIE::setZlens(PosType zl){
  
}*/

void LensHaloBaseNSIE::reNormSubstructure(PosType kappa_sub){
	/* renomalizes substructure so that
	 * the average surface density it kappa_sub
	 */
	  PosType avem;
	  avem=sub_Mmax*(sub_alpha+1)
	    /(sub_alpha+2)*(1-pow(sub_Mmin/sub_Mmax,sub_alpha+2))/
	    (1-pow(sub_Mmin/sub_Mmax,sub_alpha+1));

	  sub_Ndensity=kappa_sub*Sigma_crit/avem;

	  return ;
}

/// Sets parameters within BaseLens that depend on the source redshift - Dl,Sigma_crit,etc.
void LensHaloAnaNSIE::setCosmology(const COSMOLOGY& cosmo)
{
  PosType zlens = getZlens();
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

  perturb_rms = new PosType[6];

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
	cout << "zlens " << getZlens() << endl;

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

std::size_t LensHaloBaseNSIE::Nparams() const
{
	return LensHalo::Nparams() + 3;
}

PosType LensHaloBaseNSIE::getParam(std::size_t p) const
{
	if(p < LensHalo::Nparams())
		return LensHalo::getParam(p);
	
	switch(p - LensHalo::Nparams())
	{
		case 0:
//			return std::log10(sigma);
            return sigma*0.02;
		case 1:
			return fratio;
		case 2:
			return pa/pi;
		default:
			throw std::invalid_argument("bad parameter index for getParam()");
	}
}

PosType LensHaloBaseNSIE::setParam(std::size_t p, PosType val)
{
	using Utilities::between;
	
	if(p < LensHalo::Nparams())
		return LensHalo::setParam(p, val);
	
	switch(p - LensHalo::Nparams())
	{
		case 0:
//          return (sigma = std::pow(10., val));
            return (sigma = (float)between<PosType>(val*50., 100., 1000.));
		case 1:
			return (fratio = (float)between<PosType>(val, 1e-10, 1.));
		case 2:
			return (pa = (float)between<PosType>(pi*val, -pi/2, pi/2));
		default:
			throw std::invalid_argument("bad parameter index for setParam()");
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
	// std::cout << "deleting lens" << endl;

	delete[] perturb_rms;

	if(perturb_Nmodes > 0){
		// std::cout << "deleting modes" << endl;
		delete[] perturb_modes;
	}
	if(sub_N > 0 && substruct_implanted){
		// std::cout << "deleting subs" << endl;
		Utilities::free_PosTypeMatrix(sub_x,sub_N,2);
		delete[] subs;
		delete[] sub_substructures;
	}
	if(stars_N > 0 && stars_implanted){
		// std::cout << "deleting stars" << endl;
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


/** \brief This function returns the lensing quantities for an asymmetric version of the symmetric baseclass halo.
 *
 *  This function should only be used by the second generation of classes derived from LensHalo.
 *
 *  The math needs to be PosType checked and the sign convention checked.
 *  The method used to make the lenses asymmetric is laid out in http://metcalf1.bo.astro.it/~bmetcalf/ExtraNotes/notes_elliptical.pdf
 */

/// Derivatives of the potential factor with respect to theta
void LensHalo::faxial0(PosType theta,PosType f[]){
    int i,k;
    //std::cout<< mod[4] << std::endl;
    f[0] = mod[0]; // why is it commented out?
    f[1] = f[2] = 0;
    for(i=4;i<Nmod;i+=2){
        k=i/2;
        f[0] +=  mod[i]*cos(k*theta)   + mod[i+1]*sin(k*theta);
        f[1] += -mod[i]*k*sin(k*theta) + mod[i+1]*k*cos(k*theta);
        f[2] += -mod[i]*k*k*cos(k*theta) - mod[i+1]*k*k*sin(k*theta);
    }
}

void LensHalo::faxial1(PosType theta,PosType f[]){
    int i,k;
    //std::cout<< mod[4] << std::endl;
    f[0] = mod1[0]; // why is it commented out?
    f[1] = f[2] = 0;
    for(i=4;i<Nmod;i+=2){
        k=i/2;
        f[0] +=  mod1[i]*cos(k*theta)   + mod1[i+1]*sin(k*theta);
        f[1] += -mod1[i]*k*sin(k*theta) + mod1[i+1]*k*cos(k*theta);
        f[2] += -mod1[i]*k*k*cos(k*theta) - mod1[i+1]*k*k*sin(k*theta);
    }
}

void LensHalo::faxial2(PosType theta,PosType f[]){
    int i,k;
    //std::cout<< mod[4] << std::endl;
    f[0] = mod2[0]; // why is it commented out?
    f[1] = f[2] = 0;
    for(i=4;i<Nmod;i+=2){
        k=i/2;
        f[0] +=  mod2[i]*cos(k*theta)   + mod2[i+1]*sin(k*theta);
        f[1] += -mod2[i]*k*sin(k*theta) + mod2[i+1]*k*cos(k*theta);
        f[2] += -mod2[i]*k*k*cos(k*theta) - mod2[i+1]*k*k*sin(k*theta);
    }
}

/// Derivatives of the potential damping factor with respect to r ... TODO: come up with a better damping faction
void LensHalo::gradial(PosType r,PosType g[]){
  double r_eps=0.5*Rmax;

  PosType x = (1+r/r_eps);
  g[0] = 1.0/x/x;
  g[1] = -2.0*g[0]/x/r_eps;
  g[2] = -3.0*g[1]/x/r_eps;
}

void LensHalo::gradial2(PosType r,PosType mu, PosType sigma, PosType g[]){
    // gaussian
    //PosType x = r;
    //sigma=sigma*rscale;
    //g[0] = (exp(-1.0*(x-mu)*(x-mu)/2/sigma/sigma));
    //g[1] = g[0]*(-1.0*(x-mu)/sigma/sigma);
    //g[2] = g[0]*((x-mu)*(x-mu)-sigma*sigma)/pow(sigma,4);

    double r_eps=Rmax;
    // former tweaked
    PosType x = (1+exp(-sigma*(r-mu)/r_eps));
    g[0] = 1.0-1.0/x;
    g[1] = -1.0*sigma/x/x *(x-1);
    g[2] = -1.0*sigma*sigma*(x-1)/x/x*(2*(x-1)-1);
    //if(r<mu){
    //    std::cout << r << " " << mu << " " << g[0] << std::endl;
    //}
}

/// Calculates fourier-coefficients for power law halo
double LensHalo::fourier_coeff(double n, double q, double beta){
    struct fourier_func f(n,q,beta);
    return Utilities::nintegrate<fourier_func>(f,0.0,2*pi,1.0e-6);
}

/// Calculates potential (phi_int) from alpha_h. If flag is_alphah_a_table is True it takes and integrates directly the gfunction instead of alpha_h. The gfunction is used for the InterpolationTable used in alpha_h. Setting the flag to False speeds up the calculation of phi_h.
PosType LensHalo::alpha_int(PosType x){
    struct Ig_func g(*this);
    return Utilities::nintegrate<Ig_func>(g,1E-8,x,1.0e-6);
}

/// Calculates the modes for fourier expansion of power law halo. All the modes are relative to the zero mode to conserve mass throughout the calculation of kappa etc.
void LensHalo::calcModes(double q, double beta, double rottheta, PosType my_mod[]){
    int i,k;
    for(int i=1;i<Nmod;++i){
		my_mod[i]=0;
	}
	// fill in modes with their values for an elliptical lens
	if(q != 1.0){
        my_mod[0] = fourier_coeff(0, q, beta)/pi/2.;
        for(i=4;i<Nmod;i+=2){
            k=i/2;
            assert(i<=Nmod);
            my_mod[i] = beta*beta*fourier_coeff(k, q, beta)/pi/(beta*beta-k*k)/my_mod[0];
        }
        my_mod[0]=1.0;
    }
	else{
		cout << "circular case" << endl;
		my_mod[0]=1.0;
	}
  
    std::cout << "calcModes for beta=" << beta << " " << my_mod[0] << " " << my_mod[4] << " " << my_mod[8] << " " << std::endl;
    //for(int i=1;i<Nmod;++i){
    //    std::cout << i << " " << my_mod[i] << std::endl;
    //}
    // rotate model
    RotateModel(rottheta,my_mod,Nmod,0);
}

/// In alpha_asym phi_int(x) is used, which calculates the potential unlike phi_h from alpha_h.

void LensHalo::alpha_asym(PosType x,PosType theta, PosType alpha[]){
	PosType F,f[3],g[3],alpha_r,alpha_theta;
    PosType phi=phi_int(x);

    faxial1(theta,f);
    F=f[0]-1;
    gradial(x,g);
    beta=get_slope();
    double fac=1.0/(beta*beta/(2-beta)/(2-beta));
    
    //alpha_r=alpha_h(x)*f[0]; // w/o damping
    //alpha_theta=f[1]*phi/x; //  w/o damping
    
    alpha_r=alpha_h(x)*(1+F*g[0])+phi*fac*F*g[1]; // with damping
    alpha_theta=f[1]*g[0]*phi*fac/x; //  with damping
    
    //std::cout << "in alpha_asym: " << beta << " " << alpha_theta << std::endl;
    
	alpha[0] = (alpha_r*cos(theta) - alpha_theta*sin(theta))/cos(theta);
	alpha[1] = (alpha_r*sin(theta) + alpha_theta*cos(theta))/sin(theta);
	return;
}

/// In kappa phi_int(x) is used, which calculates the potential unlike phi_h from alpha_h.

PosType LensHalo::kappa_asym(PosType x,PosType theta){
	PosType F, f[3],g[3], kappa;
    PosType phi=phi_int(x);
    PosType f0[3], f1[3], f2[3];
    
    faxial0(theta,f0);
    faxial1(theta,f1);
    faxial2(theta,f2);
    //gradial2(x,1,2,g);
    gradial(x,g);
    g[0]=1;
    g[1]=0;
    g[2]=0;
    /*if(x<rscale){
        std::cout << x << " " << rscale << std::endl;
    
     f[0]=(f1[0]-f2[0])*x+f2[0];
     f[1]=(f1[1]-f2[1])*x+f2[1];
     f[2]=(f1[2]-f2[2])*x+f2[2];

    };*/
    
    
    double x0=0.2, x1=1.0, x2=2.5, x3=5;
//double y0=f0[0], y1=f1[0], y2=f0[0], y3=0.0;
    double y0=f0[0], y1=f1[0], y2=f2[0], y3=0.0;
    
    
    f[0]=(x-x1)*(x-x2)*(x-x3)/((x0-x1)*(x0-x2)*(x0-x3))*y0+(x-x0)*(x-x2)*(x-x3)/((x1-x0)*(x1-x2)*(x1-x3))*y1+(x-x0)*(x-x1)*(x-x3)/((x2-x0)*(x2-x1)*(x2-x3))*y2+(x-x0)*(x-x1)*(x-x2)/((x3-x0)*(x3-x1)*(x3-x2))*y3;
    
    y0=f0[1]; y1=f1[1]; y2=f2[1]; y3=0.0;
    
    f[1]=(x-x1)*(x-x2)*(x-x3)/((x0-x1)*(x0-x2)*(x0-x3))*y0+(x-x0)*(x-x2)*(x-x3)/((x1-x0)*(x1-x2)*(x1-x3))*y1+(x-x0)*(x-x1)*(x-x3)/((x2-x0)*(x2-x1)*(x2-x3))*y2+(x-x0)*(x-x1)*(x-x2)/((x3-x0)*(x3-x1)*(x3-x2))*y3;
    
    y0=f0[2]; y1=f1[2]; y2=f2[2]; y3=0.0;
    
    f[2]=(x-x1)*(x-x2)*(x-x3)/((x0-x1)*(x0-x2)*(x0-x3))*y0*(1.0/(0.5*0.5/(2-0.5)/(2-0.5)))+(x-x0)*(x-x2)*(x-x3)/((x1-x0)*(x1-x2)*(x1-x3))*y1*(1.0/(1/(2-1)/(2-1)))+(x-x0)*(x-x1)*(x-x3)/((x2-x0)*(x2-x1)*(x2-x3))*y2*(1.0/(1.5*1.5/(2-1.5)/(2-1.5)))+(x-x0)*(x-x1)*(x-x2)/((x3-x0)*(x3-x1)*(x3-x2))*y3*(1.0/(1.5*1.5/(2-1.5)/(2-1.5)));
    
    
    
    //faxial(theta,f);
    
    F=f[0]-1;

   // PosType F2, f2[3],g2[3];
   // faxial2(theta,f2);
   // F2=f2[0]-1;
   // gradial2(x,1.0,3,g);
    
    
    beta=get_slope();
    double fac=1.0/(beta*beta/(2-beta)/(2-beta));
    
    
    kappa=f[0]*kappa_h(x)-0.5*f[2]*phi;//  w/o damping
    
    // kappa=(1+(F*g[0]+F2*g2[0])/(g[0]+g2[0]))*kappa_h(x)-0.5*phi*fac*((F*g[1]+F2*g2[1])/(g[1]+g2[1])/x+(F*g[2]+F2*g2[2])/(g[2]+g2[2])+(f[2]*g[0]+f2[2]*g2[0])/(g[0]+g2[0])/x/x)*x*x-(F*g[1]+F2*g2[1])/(g[1]+g2[1])*alpha_h(x)*x*x; /// with bi-gaussian damping
	
    //kappa=(1+F*g[0])*kappa_h(x)-0.5*phi*fac*(F*g[1]/x+F*g[2]+f[2]*g[0]/x/x)*x*x-F*g[1]*alpha_h(x)*x*x; /// with damping
	return kappa;
}

/// In gamma_asym phi_int(x) is used, which calculates the potential unlike phi_h from alpha_h.

void LensHalo::gamma_asym(PosType x,PosType theta, PosType gamma[]){
	PosType F, f[3],g[3];
    PosType phi=phi_int(x);
    
    
    faxial1(theta,f);
    F=f[0]-1;
    gradial(x,g);
    beta=get_slope();
    double fac=1.0/(beta*beta/(2-beta)/(2-beta));
    //double gt = f[0]*gamma_h(x)+0.5*phi*fac*f[2];// w/o damping
    //double g45 = (-alpha_h(x)*f[1]*g[0])*x+(phi*fac*f[1]);// w/o damping
    
    PosType gt = (1+F*g[0])*gamma_h(x)+0.5*phi*fac*(-F*g[1]/x+F*g[2]-f[2]*g[0]/x/x)*x*x-F*g[1]*alpha_h(x)*x*x;// with damping
    PosType g45 = (-alpha_h(x)*f[1]*g[0]-phi*fac*f[1]*g[1])*x+(phi*fac*f[1]*g[0]);// with damping
    
    gt *= 0.5*pow(x/xmax,2);
    g45 *= 0.5*pow(x/xmax,2);
    
	gamma[0] = cos(2*theta)*gt-sin(2*theta)*g45;
    gamma[1] = sin(2*theta)*gt+cos(2*theta)*g45;
    
	return;
}


/* makes beta=-2 Power Law Elliptical

 
 // for Ansatz I
 void LensHalo::felliptical(double x, double q, double theta, double f[], double g[]){ // q not a function of r
 double A;
 //reps=rmax
 //q=r/reps+q0*(1-r/reps) // q=q0 for small radii, q=1 for large radii
 A=1/q/q;
 f[0]=pow((cos(theta)*cos(theta)+A*sin(theta)*sin(theta)),-0.5);
 g[0]=1; //f[0]/x;
 f[1]=-((A-1.)*cos(theta)*sin(theta))*pow(f[0],3);
 f[2]=(A-1.)*(4*(A+1.)*cos(2.*theta)+(A-1)*(cos(4*theta)-5.))/(pow(2,0.5)*pow(1.+A-(A-1.)*cos(2*theta),2.5));
 g[1]=0.; //f[0]/x;
 g[2]=0.;
 }
 
 void LensHalo::felliptical(double x, double q, double theta, double f[], double g[]){ // q not a function of r
 double A;
 //reps=rmax
 //q=r/reps+q0*(1-r/reps) // q=q0 for small radii, q=1 for large radii
 A=1/q/q;
 f[0]=pow((cos(theta)*cos(theta)+A*sin(theta)*sin(theta)),-0.5);
 g[0]=1; //f[0]/x;
 f[1]=-((A-1.)*cos(theta)*sin(theta))*pow(f[0],3);
 f[2]=(A-1.)*(4*(A+1.)*cos(2.*theta)+(A-1)*(cos(4*theta)-5.))/(pow(2,0.5)*pow(1.+A-(A-1.)*cos(2*theta),2.5));
 g[1]=0.; //f[0]/x;
 g[2]=0.;
 }
 
 double LensHalo::kappa_asym(double x,double theta){
 double F, f[3],g[3], kappa;
 
 felliptical(x,0.4,theta,f,g);
 g[0]=1;
 g[1]=g[2]=0;
 
 F=f[0];
 kappa=F*kappa_h(x)-0.5*f[2]/x/x*phi_h(x);
 
 return kappa/x;
 }
 */

/* ANSATZ II
 // eq 34
 void LensHalo::felliptical(double x, double q, double theta, double f[], double g[]){
 double A[3];
 //reps=rmax
 //q=r/reps+q0*(1-r/reps) // q=q0 for small radii, q=1 for large radii
 A[0]=1/q/q;
 A[1]=0.;
 A[2]=0.;
 f[0]=x*pow((cos(theta)*cos(theta)+A[0]*sin(theta)*sin(theta)),0.5);
 g[0]=f[0]/x;
 f[1]=x*(A[0]-1)*(cos(theta)*sin(theta))/(g[0]);
 f[2]=(x*(A[0]-1)*(pow(cos(theta),4)-A[0]*pow(sin(theta),4)))/(g[0]*g[0]*g[0]);
 g[1]=g[0]; // (g[0]*g[0]+0.5*x*A[1]*sin(theta)*sin(theta))/g[0];
 g[2]=0; //sin(theta)*sin(theta)*(A[1]*(4*g[0]*g[0]-x*A[1]*sin(theta)*sin(theta))+2*f[0]*g[0]*A[2] )/(4*g[0]*g[0]*g[0]);
 //g[3]=cos(theta)*sin(theta)*(2*(A[0]*A[0]-1-(A[0]-1)*(A[0]-1)*cos(2*theta))+x*A[1]*(3+cos(2*theta)+2*A[0]*sin(theta)*sin(theta)))/(4*g[0]*g[0]*g[0]);
 }
 
 double LensHalo::kappa_asym(double x,double theta){
 double f[3],g[3];
 double G,Gr,Grr, Gt,Gtt, kappa;
 felliptical(x,0.4,theta,f,g);
 G = f[0];
 Gt = f[1];
 Gtt= f[2];
 Gr=g[1];
 Grr=g[2];
 double phitwo=(2*kappa_h(G)/G/G - alpha_h(G)/G/G);
 kappa=-0.5*alpha_h(G)/G*(Gr/x+Grr+Gtt/x/x)+0.5*phitwo*(Gr*Gr+Gt*Gt/x/x);
 kappa/=G*G;
 return kappa;
 }
 */

