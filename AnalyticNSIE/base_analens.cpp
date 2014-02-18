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
		,PosType *xcm
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
			return std::log10(sigma);
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
			return (sigma = std::pow(10., val));
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


/// TODO: This needs to be thought about more. sets axial modes to reproduce a near elliptically shaped surface density
void LensHalo::setModesToEllip(PosType q,PosType theta){
  // elliptical integrals

	PosType K = rfD(0,1./q/q,1);
	PosType E = K - (1-1./q/q)*rdD(0,1./q/q,1)/3;
  assert(Nmod == 18);
  
    
    
  // set modo to elliptical model
	for(int i=1;i<=Nmod;++i){
		mod[i]=0.0;
	}
	// fill in modes with their values for an elliptical lens
	if(q != 1.0){
		mod[3] = 4*K/pi; // /2?
		mod[4] = 4*( (1+q*q)*K-2*q*q*E )/(1-q*q)/pi/mod[3];
		mod[8] = 4*( (3*q*q+1)*(q*q+3)*K-8*q*q*(1+q*q)*E )/( 3*pi*pow(1-q*q,2) )/mod[3];
		mod[12]= 4*( (1+q*q)*(15+98*q*q+15*q*q*q*q)*K-2*q*q*(23+82*q*q+23*q*q*q*q)*E )/( 15*pi*pow(1-q*q,3) )/mod[3];
		mod[16]= 4*( -32*q*q*(1+q*q)*(11+74*q*q+11*q*q*q*q)*E+(105+1436*q*q+3062*q*q*q*q+1436*pow(q,6)+105*pow(q,8))*K )/(105*pi*pow(1-q*q,4))/mod[3];
	}
    mod[3]=1.0;
  
	// rotate model
	RotateModel(theta,mod,Nmod,0);
  
  return;
}


/// Derivatives of the axial potential factor with respect to theta

void LensHalo::fangular(PosType theta,PosType f[]){
  int i,k;

  for(i=0;i<=Nmod/2;++i){
	k=2*i;
    f[0] +=  mod[k]*cos(i*theta)     + mod[k+1]*sin(i*theta);
    f[1] += -mod[k]*i*sin(i*theta)   + mod[k+1]*i*cos(i*theta);
    f[2] += -mod[k]*i*i*cos(i*theta) - mod[k+1]*i*i*sin(i*theta);
    //std::cout << f[0] << " " << k << " " << mod[k] << " " << mod[k+1] << endl;
  }
	//cout << "fangular=" << "\n" << f[0] << "\n" << f[1] << "\n" << f[2] << "\n" << endl;

//	cout << mod[0] << " " << mod[1] << " " << mod[2] << " " << mod[3] << " " << mod[4] << " " << mod[5] << " " << mod[6] << " " << mod[7] << " " << mod[8] << " " << mod[9] << endl;

//	cout << mod[10] << " " << mod[11] << " " << mod[12] << " " << mod[13] << " " << mod[14] << " " << mod[15] << " " << mod[16] << " " << mod[17] << " " << mod[18] << endl;
}



/// Derivatives of the axial potential factor with respect to theta
void LensHalo::gradial(PosType r,PosType g[]){
  double r_eps=0.1*Rmax; // TODO: r_eps = Rmax for now, but must be thought about later
  PosType x = (1+r/r_eps);
  g[0] = 1.0/x/x;
  g[1] = -2.0*g[0]/x/r_eps;
  g[2] = -3.0*g[1]/x/r_eps;
  //cout << "ginside: rmax " << Rmax  << " "<< g[0] << " " << g[1] << " " << g[2] << endl;
}

/** \brief This function returns the lensing quantities for an asymmetric version of the symmetric baseclass halo.
 *  
 *  This function should only be used by the second generation of classes derived from LensHalo.
 *
 *  The math needs to be PosType checked and the sign convention checked.  
 *  The method used to make the lenses asymmetric is laid out in http://metcalf1.bo.astro.it/~bmetcalf/ExtraNotes/notes_elliptical.pdf
 */
void LensHalo::desymmeterize(PosType r,PosType theta,PosType *alpha,PosType *kappa,PosType *gamma){
  PosType f[3],g[3];
  
  PosType alpha_iso = alpha_h(r/rscale),phi_iso = phi_h(r/rscale)
  ,kappa_iso = kappa_h(r/rscale),gamma_iso = gamma_h(r/rscale);
  
  PosType alpha_r,alpha_theta,F;

  faxial(theta,f);
  gradial(r,g);
  
  F = (1+g[0]*f[0]);
  
  alpha_r = (F + g[1]*f[0])*alpha_iso;
  alpha_theta = g[0]*f[1]*phi_iso/r;
  
  alpha[0] = alpha_r*cos(theta) - alpha_theta*sin(theta);
  alpha[1] = alpha_r*sin(theta) + alpha_theta*cos(theta);
  
  *kappa = F*kappa_iso + ( (g[2] + g[1]/r)*f[0] + g[0]*f[2]/r/r )*phi_iso;
  
  PosType gt = F*gamma_iso + g[1]*f[0]*alpha_iso + 0.5*( g[2]*f[0] - g[0]*f[2]/r/r)*phi_iso;
  PosType g45 = f[1]*(alpha_iso*g[0]/r + (g[1]-g[0]/r/r)*phi_iso);
  
  gamma[0] = cos(2*theta)*gt + sin(2*theta)*g45;
  gamma[1] = -sin(2*theta)*gt + cos(2*theta)*g45;
  
}

double LensHalo::fourier_coeff(double n, double q, double beta){
    struct fourier_func f(n,q,beta);
    //f.n = n;
    //f.q = q;
    //f.beta = beta;
    return Utilities::nintegrate<fourier_func>(f,0.0,2*pi,1.0e-6);
}


void LensHalo::setEllipModes(double q,double theta){
  // elliptical integrals
    int i,k;
	PosType K = rfD(0,1./q/q,1);
	PosType E = K - (1-1./q/q)*rdD(0,1./q/q,1)/3;
    assert(Nmod == 32);
  //assert(q == 0.3);
    

    
  // set modo to elliptical model
	for(int i=1;i<=Nmod;++i){
		mod[i]=0;
	}
    PosType beta=1;
    
	// fill in modes with their values for an elliptical lens
	if(q != 1.0){
        mod[0] = fourier_coeff(0, q, beta)/pi;
        for(i=4;i<Nmod;i+=2){
            k=i/2;
            mod[i] = fourier_coeff(k, q, beta)/pi/(beta*beta-k*k);
        }
        //mod[0]=1;
        /*
		mod[0] = 4*K/pi; //2?;
		mod[4] = 4*( (1+q*q)*K-2*q*q*E )/(1-q*q)/pi/(beta*beta-2*2); ///mod[0];fourier_coeff(2)/pi/(beta*beta-2*2); //fourier_coeff(2, q, beta)/(beta*beta-2*2);//
        mod[4] = fourier_coeff(2, q, beta)/pi/(-3);
		mod[8] = 4*( (3*q*q+1)*(q*q+3)*K-8*q*q*(1+q*q)*E )/( 3*pi*pow(1-q*q,2) )/(beta*beta-4*4); // fourier_coeff(4, q, beta)/(beta*beta-4*4); // ///mod[0];
		mod[12]=4*( (1+q*q)*(15+98*q*q+15*q*q*q*q)*K-2*q*q*(23+82*q*q+23*q*q*q*q)*E )/( 15*pi*pow(1-q*q,3) )/(beta*beta-6*6); //fourier_coeff(6, q, beta)/(beta*beta-6*6); // ///mod[0];
		mod[16]=4*( -32*q*q*(1+q*q)*(11+74*q*q+11*q*q*q*q)*E+(105+1436*q*q+3062*q*q*q*q+1436*pow(q,6)+105*pow(q,8))*K )/(105*pi*pow(1-q*q,4))/(beta*beta-8*8); //fourier_coeff(8, q, beta)/(beta*beta-8*8); //  ///mod[0]; */
		//mod[18]=0.;
	}
	else{
		cout << "here in setEllipModes" << endl;

		mod[0]=1.0;

	}
    //mod[0]=1.0;
	// rotate model
	//RotateModel(theta,mod,Nmod,0);

  return;
}

void LensHalo::calcModes(double q, double beta, double rottheta, PosType newmod[]){
    int i,k;
	//assert(Nmod == 32);
    for(int i=1;i<=Nmod;++i){
		mod[i]=0;
	}
	// fill in modes with their values for an elliptical lens
	if(q != 1.0){
        mod[0] = fourier_coeff(0, q, beta)/pi/2.;
        for(i=4;i<Nmod;i+=2){
            k=i/2;
            mod[i] = fourier_coeff(k, q, beta)/pi/(beta*beta-k*k);
        }
    }
	else{
		cout << "here in setEllipModes" << endl;
		mod[0]=1.0;
	}
    // rotate model
    std::cout << "calcModes for beta=" << beta << " " << mod[0] << " " << mod[4] << " " << mod[8] << " " << std::endl;
	RotateModel(rottheta,mod,Nmod,0);
}

/// Derivatives of the axial potential factor with respect to theta
void LensHalo::faxial(PosType theta,PosType f[]){
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


void LensHalo::alpha_asym(PosType x,PosType theta, PosType alpha[]){
	double f[3],g[3];
	double alpha_r,alpha_theta,F;

    faxial(theta,f);
    F=f[0]-1;
    gradial(x,g);

    PosType beta=get_slope();
	
    //alpha_r=-beta*(1+F*g[0])*pow(x,-beta-1)-F*g[1]*pow(x,-beta); // for power law halo with damping
    alpha_r=alpha_h(x)*pow(x,-2*beta-3)*beta*(1+F*g[0])-phi_h(x)*(2-beta)*pow(x,-2)*(F*g[1]); // for arbitrary halos with damping
    
    //alpha_theta=-f[1]*g[0]*pow(x,-beta-1); // for power law halo with damping
    alpha_theta=-f[1]*g[0]/x*phi_h(x)*(2-beta)*pow(x,-2); // for arbitrary halos with damping

    alpha_r*= pow(x,3);
    alpha_theta*= pow(x,3);

	alpha[0] = (alpha_r*cos(theta) - alpha_theta*sin(theta))/cos(theta);
	alpha[1] = (alpha_r*sin(theta) + alpha_theta*cos(theta))/sin(theta);
	return;
}


double LensHalo::kappa_asym(PosType x,PosType theta){
	PosType F, f[3],g[3], kappa;
    PosType beta=get_slope();

    faxial(theta,f);
    gradial(x,g);
    F=f[0]-1;
    g[0]=1;
    g[1]=0;
    g[2]=0;
    //kappa=f[0]*kappa_h(x)*beta*beta/(2+beta)/pow(x,4)-0.5*f[2]/x/x*phi_h(x)*(2-beta)*pow(x,2*beta-2); // for arbitrary halos w/o damping
    
    //kappa=0.5*(1+F*g[0])*beta*beta*pow(x,beta-2)+0.5*(F*g[1]/x+F*g[2]+f[2]*g[0]/x/x)*pow(x,beta)                      +2*F*g[1]*beta*pow(x,beta-1); // for powerlaw halo with damping
    
    //kappa=0.5*((beta*beta*(1+F*g[0])+g[0]*f[2])*pow(x,beta-2)+((2*beta+1)*F*g[1])*pow(x,beta-1)+(F*g[2])*pow(x,beta)); // equivalent to latter expression: for powerlaw halo with damping
    
    //kappa=0.5*(1+F*g[0])*beta*beta*2*kappa_h(x)/(2+beta)/pow(x,4)-0.5*(F*g[1]/x+F*g[2]+f[2]*g[0]/x/x)*phi_h(x)*(2-beta)*pow(x,2*beta-2)-2*F*g[1]*beta*alpha_h(x)/pow(x,3); // for arbitrary halos with damping
    //std::cout <<  beta << " " << kappa_h(x) << std::endl;
   
    //kappa=0.5*(1+F*g[0])*beta*beta*2*kappa_h(x)-0.5*(f[2]*g[0]/x/x)*phi_h(x)*4*pi; // for arbitrary halos with damping
    
    //if(x<4.185 && x>4.18){
      //  std::cout << "in kappa_asym: " <<  x << " " << kappa_h(x) << std::endl;
    //}
    
    double phi_iso=-x*x*phi_h(x);
    double alpha_iso=-1*alpha_h(x);
    kappa=kappa_h(x)*f[0]+0.5*phi_iso*f[2]/x/x; // w/o damping function NFW
    //kappa=kappa_h(x)*(1+F*g[0])+0.5*F*phi_iso*x*x*(g[1]/x+g[2]+f[2]/x/x);
    
	return kappa; //*(1/(beta*beta))*pow(x,-2*beta+4);
}



void LensHalo::gamma_asym(PosType x,PosType theta, PosType gamma[]){
	double f[3],g[3];
	double F;
	PosType beta=get_slope();

    faxial(theta,f);
    F=f[0]-1;
    gradial(x,g);

    
    //double gt = 0.5*(beta*(beta-2)*f[0]-f[2])*pow(x,beta-2); // w/o damping
    //double g45 = (beta-1)*f[1]*pow(x,beta-2) ; // w/o damping
    
    //double gt = 0.5*( (beta*(beta-2)*(1+F*g[0]) -g[0]*f[2])*pow(x,beta-2)+(2*beta-1)*F*g[1]*pow(x,beta-1)+F*g[2]*pow(x,beta)); // direct calculation for powerlaw with damping
    double gt = 0.5*((1+F*g[0])*gamma_h(x)*pow(x,-4)*(beta-2)-0.5*(-F*g[1]/x+F*g[2]-f[2]*g[0]/x/x)*phi_h(x)*pow(x,2*beta-2)*(2-beta)-2*F*g[1]*beta*alpha_h(x)/pow(x,3)); // for arbitrary halos
    
    //double g45 = f[1]/x*((beta-1)*g[0]*pow(x,beta-1)+g[1]*pow(x,beta)) ; // direct calculation for powerlaw with damping
    double g45 = 0.5*f[1]/x*( -alpha_h(x)/pow(x,3)*g[0]-phi_h(x)*pow(x,2*beta-2)*(2-beta)*(g[1]-g[0]/x) ) ; //for arbitrary halos
    
    gt *= pow(x,2*beta+2);
    g45 *= pow(x,2*beta+2);
    
	gamma[0] = cos(2*theta)*gt-sin(2*theta)*g45;
    gamma[1] = sin(2*theta)*gt+cos(2*theta)*g45;
    
	return;
}









/* makes beta=-2 Power Law Elliptical
 
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

