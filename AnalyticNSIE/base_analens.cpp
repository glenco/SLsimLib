/*
 * base_analens.cpp
 *
 *  Created on: Jan 22, 2013
 *      Author: mpetkova
 */

//#include "slsimlib.h"
#include "lens_halos.h"
#include "base_analens.h"
#include "analytic_lens.h"
#include "fitlens.h"



using namespace std;

void LensHaloBaseNSIE::force_halo(
                                  PosType *alpha       /// mass/PhysMpc
                                  ,KappaType *kappa    /// surface mass density , mass / /PhysMpc/PhysMpc
                                  ,KappaType *gamma    /// mass / /PhysMpc/PhysMpc
                                  ,KappaType *phi      /// mass
                                  ,PosType const *xcm  /// Position in PhysMpc
                                  ,bool subtract_point /// if true contribution from a point mass is subtracted
                                  ,PosType screening   /// the factor by which to scale the mass for screening of the point mass subtraction
                                  )
{

  PosType rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
  if(force_point(alpha,kappa,gamma,phi,xcm,rcm2
                 ,subtract_point,screening)) return;


  PosType alpha_tmp[2];
  KappaType gamma_tmp[2];
  KappaType phi_tmp,kappa_tmp ;
  
  gamma_tmp[0] = gamma_tmp[1] = 0.0;
  alpha_tmp[0] = alpha_tmp[1] = 0.0;
  phi_tmp = kappa_tmp = 0.0 ;
  
  
  
  if(sigma > 0.0){
    PosType xt[2]={0,0};
    float units = pow(sigma/lightspeed,2)/Grav ; /// sqrt(fratio); // in mass / PhysMpc
    units *= 2. * Rsize / PI ; // units now in mass /// Multiplying by 2*Rsize/PI to match with Power Law
      
    xt[0]=xcm[0]; // in PhysMpc
    xt[1]=xcm[1];
    
    alphaNSIE(alpha_tmp,xt,fratio,rcore,pa);
    alpha_tmp[0] *= units;
    alpha_tmp[1] *= units;
    
    {
      gammaNSIE(gamma_tmp,xcm,fratio,rcore,pa);
      kappa_tmp=kappaNSIE(xcm,fratio,rcore,pa);
      kappa_tmp *= units ;
      gamma_tmp[0] *= units ;
      gamma_tmp[1] *= units ;
      
      phi_tmp = -1.0 * phiNSIE(xcm,fratio,rcore,pa) ;  // phi in Mpc (physical)
      phi_tmp *= units ;                               // phi now in Msun * Mpc
    }
  }
  
  // As before :
  alpha[0] += alpha_tmp[0];
  alpha[1] += alpha_tmp[1];
  gamma[0] += gamma_tmp[0];
  gamma[1] += gamma_tmp[1];
  *phi += phi_tmp ;
  *kappa += kappa_tmp;
    
  // perturbations of host lens
  if(perturb_Nmodes > 0)
  {
    PosType xt[2]={0,0};
    xt[0]=xcm[0] ; // in PhysMpc
    xt[1]=xcm[1] ;
    
    // Calling lens_expand : xt in PhysMpc, mod[0,1,2] in mass / PhysMpc^2, mod[3,4,5,...] in mass / PhysMpc :
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    *kappa += lens_expand(perturb_beta,perturb_modes,perturb_Nmodes-1,xt,alpha_tmp,gamma_tmp,&phi_tmp);

    // For consistency with FindLensSimple and other calls of lens_expand we had removed the -1 after perturb_Nmodes.
      
    // We put it back because otherwise the modes mod[i+1] in lens_expand would be used up to mod[imax+1] with :
    // - if Nmodes (in lens_expand) is odd  : imax = Nmodes - 1 ;
    // - if Nmodes (in lens_expand) is even : imax = Nmodes - 2 ;
    // Hence, with Nmodes = perturb_Nmodes - 1 we have for example :
    // perturb_Nmodes = 7  =>  Nmodes = perturb_Nmodes - 1 = 6  =>  imax = 4  => max mode used is mod[5].
    // Without the -1 we would have :
    // perturb_Nmodes = 7  =>  Nmodes = perturb_Nmodes = 7  =>  imax = 6  => max mode used is mod[7].
    // This is not possible because 7 modes mean we use mod[0], mod[1], ... , mod[6] !
      
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    // std::cout << "alpha_tmp just after lens_expand : " << alpha_tmp[0] << " " << alpha_tmp[1] << std::endl;
    // alpha_tmp is in mass / PhyMpc here !
    
    // Printing quantities for control :
    // std::cout << "         xt in force_halo : @@@ " << xt[0]/Dls << " " << xt[1]/Dls << " @@@" << std::endl ;
    // std::cout << "alpha in force_halo : " << alpha_tmp[0] * (4*PI*Grav) << " " << alpha_tmp[1]* (4*PI*Grav) << std::endl ;
    // std::cout << "xt - alpha in force_halo : !!! " << xt[0] - alpha_tmp[0] << " " << xt[1] - alpha_tmp[1] << " !!!" << std::endl ;

    // Changing the sign because there is a change of sign in LensPlaneSingular::force :
    alpha_tmp[0] *= -1. ;
    alpha_tmp[1] *= -1. ;
    
    // As before :
    alpha[0] += alpha_tmp[0];
    alpha[1] += alpha_tmp[1];
    
    {
      gamma[0] += gamma_tmp[0];
      gamma[1] += gamma_tmp[1];

      *phi += phi_tmp ;
    }
    gamma_tmp[0] = gamma_tmp[1] = 0.0;
    alpha_tmp[0] = alpha_tmp[1] = 0.0;
    phi_tmp = 0.0;
  }
  
   // std::cout << "alpha at the end of force_halo : " << alpha[0] << " " << alpha[1] << std::endl;
    
    
  // add substructure
//  if(substruct_implanted)
//  {
//    for(j=0;j<sub_N;++j)
//    {
//
//      subs[j].force_halo(alpha_tmp,&kappa_tmp,gamma_tmp,&phi_tmp,xcm);
//
//      alpha[0] += alpha_tmp[0];
//      alpha[1] += alpha_tmp[1];
//
//      {
//        *kappa += kappa_tmp;
//        gamma[0] += gamma_tmp[0];
//        gamma[1] += gamma_tmp[1];
//        *phi += phi_tmp;
//      }
//    }
//
//    gamma_tmp[0] = gamma_tmp[1] = 0.0;
//    alpha_tmp[0] = alpha_tmp[1] = 0.0;
//    phi_tmp = 0.0;
//  }
  
  // add stars for microlensing
//  if(stars_N > 0 && stars_implanted){
//    force_stars(alpha,kappa,gamma,xcm);
//  }
  
  //assert(alpha[0] == alpha[0] && alpha[1] == alpha[1]);

  return ;
}


/**
 * \brief Reads in a parameter file and sets up an analytic lens.
 *
 * Sets many parameters within the lens model, source model and
 * force calculation.
 */
void LensHaloBaseNSIE::assignParams(InputParams& params,const COSMOLOGY &cosmo){
  double tmp;
	if(!params.get("main_mass",tmp)) error_message1("main_mass",params.filename());
  LensHalo::setMass(tmp);
	if(!params.get("main_zlens",tmp)) error_message1("main_zlens",params.filename());
  LensHalo::setZlens(tmp,cosmo);
	if(!params.get("main_sigma",sigma)) error_message1("main_sigma",params.filename());
	if(!params.get("main_core",rcore)) error_message1("main_core",params.filename());
	if(!params.get("main_axis_ratio",fratio)) error_message1("main_axis_ratio",params.filename());
  else if(fratio > 1){
    ERROR_MESSAGE();
    std::cerr << "parameter main_axis_ratio must be < 1 in file " << params.filename() << ". Use main_pos_angle to rotate the halo." << std::endl;
    exit(1);
  }
	if(!params.get("main_pos_angle",pa)) error_message1("main_pos_angle",params.filename());
  
	Rsize = rmaxNSIE(sigma,get_mass(),fratio,rcore);
  
	Rmax = MAX(1.0,1.0/fratio)*Rsize;  // redefine
  
	assert(Rmax >= Rsize);
  
	if(!params.get("zsource_reference",zsource_reference)) error_message1("zsource_reference",params.filename());
  
  // Substructure parameters
//  if(!params.get("main_sub_Ndensity",sub_Ndensity)) error_message1("main_sub_Ndensity",params.filename());
//
//  else if(sub_Ndensity > 0){
//    if(!params.get("main_sub_beta",sub_beta)) error_message1("main_sub_beta",params.filename());
//    if(!params.get("main_sub_alpha",sub_alpha)) error_message1("main_sub_alpha",params.filename());
//    if(!params.get("main_sub_Rsize",sub_Rsize)) error_message1("main_sub_Rsize",params.filename());
//    if(!params.get("main_sub_mass_max",sub_Mmax)) error_message1("main_sub_mass_max",params.filename());
//    if(!params.get("main_sub_mass_min",sub_Mmin)) error_message1("main_sub_mass_min",params.filename());
//    if(sub_Mmin < 1.0e3){
//      ERROR_MESSAGE();
//      std::cerr << "Are you sure the minimum halo mass should be " << sub_Mmin << " Msun?" << std::endl;
//      exit(1);
//    }
//    if(!params.get("main_sub_type",main_sub_type)) error_message1("main_sub_type",params.filename());
//    
//  }
  // Stars parameters
  //if(!params.get("main_stars_N",stars_N)) error_message1("main_stars_N",params.filename());
  
//  else if(stars_N){
 //   assignParams_stars(params);
 // }
  
}

void LensHaloBaseNSIE::error_message1(std::string parameter,std::string file){
		  ERROR_MESSAGE();
		  std::cout << "Parameter " << parameter << " is needed to construct a BaseAnaLens.  It needs to be set in parameter file " << file << "!" << endl;
		  exit(0);
}

/*/ resets Zl, Dl, Sigma_crit, MpcToAsec
void LensHaloBaseNSIE::setZlens(PosType zl){
  
}*/

//void LensHaloBaseNSIE::reNormSubstructure(PosType kappa_sub){
//	// renomalizes substructure so that
//  //the average surface density it kappa_sub
//	  PosType avem;
//	  avem=sub_Mmax*(sub_alpha+1)
//	    /(sub_alpha+2)*(1-pow(sub_Mmin/sub_Mmax,sub_alpha+2))/
//	    (1-pow(sub_Mmin/sub_Mmax,sub_alpha+1));
//
//	  sub_Ndensity=kappa_sub*Sigma_crit/avem;
//
//	  return ;
//}

/// Sets parameters within BaseLens that depend on the source redshift - Dl,Sigma_crit,etc.
void LensHaloAnaNSIE::setCosmology(const COSMOLOGY& cosmo)
{
  PosType zlens = getZlens();
	Dl = cosmo.angDist(0,zlens);
	Ds = cosmo.angDist(0,zsource_reference);
	Dls = cosmo.angDist(zlens,zsource_reference);
	MpcToAsec = 60*60*180 / PI / Dl;
		// in Mpc
	Einstein_ro=4*PI*pow(sigma/lightspeed,2)*Dl
		*Dls/Ds;
	// find critical density
	Sigma_crit=Ds/Dls/Dl/4/PI/Grav;
	to = (1+zlens)*Ds/Dls/Dl/8.39428142e-10;
}
void LensHaloFit::setCosmology(const COSMOLOGY& cosmo)
{
  PosType zlens = getZlens();
  Dl = cosmo.angDist(0,zlens);
  Ds = cosmo.angDist(0,zsource_reference);
  Dls = cosmo.angDist(zlens,zsource_reference);
  MpcToAsec = 60*60*180 / PI / Dl;

  // find critical density
  Sigma_crit=Ds/Dls/Dl/4/PI/Grav;
  to = (1+zlens)*Ds/Dls/Dl/8.39428142e-10;
}


//LensHaloBaseNSIE::LensHaloBaseNSIE(InputParams& params,const COSMOLOGY &cosmo) : LensHalo(){
//
//  perturb_rms = new PosType[6];
//
//  assignParams(params,cosmo);
//
//  // parameters for stars
//  stars_implanted = false; // stars are implanted later
//  star_theta_force = 0.1;
//  sub_theta_force = 0.1;
//
//  perturb_Nmodes = 0;
//  //sub_sigmaScale = sigma = pa = Einstein_ro = fratio = rcore = 0.0;
//
//  if(sub_Ndensity == 0)
//	  sub_N = 0;
//
//  Sigma_crit = 0;
//
//  substruct_implanted = false;
//
//}
LensHaloBaseNSIE::LensHaloBaseNSIE() : LensHalo(){
  
  perturb_rms = new PosType[6];
  
  // parameters for stars
  //stars_implanted = false; // stars are implanted later
  //star_theta_force = 0.1;
  //sub_theta_force = 0.1;
  
  perturb_Nmodes = 0;
  //sub_sigmaScale = sigma = pa = Einstein_ro = fratio = rcore = 0.0;
  
  //if(sub_Ndensity == 0)
  //  sub_N = 0;
  
  Sigma_crit = 0;
  
  substruct_implanted = false;
}

void LensHaloBaseNSIE::PrintLens(bool show_substruct,bool show_stars){
	cout << "zlens " << getZlens() << endl;

	 // parameters of substructures
	//cout << endl << "main_sub_Ndensity "<< sub_Ndensity << endl;
	//if(sub_Ndensity > 0){
		//cout << "betaSubstruct "<<sub_beta << endl;
		//cout << "alphaSubstruct "<<sub_alpha << endl;
		//cout << "RsizeSubstruct "<<sub_Rsize << " Mpc" << endl;
		//cout << "MmaxSubstruct "<<sub_Mmax << " Msun" << endl;
		//cout << "MminSubstruct "<<sub_Mmin << " Msun\n" << endl;
	//}

//	if(sub_N > 0){
//		//cout << endl << "NSubstruct "<< sub_N << endl;
//		if(show_substruct){
//			if(substruct_implanted || sub_N > 0){
//				for(i=0;i<sub_N;++i){
//				  cout << "RcutSubstruct "<<i << " " <<subs[i].getRsize() << " Mpc" << endl;
//				  cout << "massSubstruct "<<i<<" "<<subs[i].get_mass() << " Msun" << endl;
//				  cout << "xSubstruct "<<i<<" "<<sub_x[i][0]<<" "<<sub_x[i][1] << " Mpc" << endl;
//					switch(main_sub_type){
//					case nfw:
//						cout << "  NFW clumps" << endl;
//						break;
//					case powerlaw:
//						cout << "  Power Law clumps" << endl;
//						break;
//					case pointmass:
//						cout << "  Point Mass clumps" << endl;
//						break;
//					default:
//						ERROR_MESSAGE();
//						cerr << "ERROR: no submass internal profile chosen" << endl;
//						exit(1);
//						break;
//					}
//				}
//			}else cout << "substructures are not implanted yet" << endl;
//		}
//	}

	if(Sigma_crit)
		cout << "critical density is " << Sigma_crit << " Msun/Mpc^2" << endl << endl;

	//if (stars_implanted) PrintStars(show_stars);
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
			return pa/PI;
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
			return (pa = (float)between<PosType>(PI*val, -PI/2, PI/2));
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
	//if(sub_N > 0 && substruct_implanted){
		// std::cout << "deleting subs" << endl;
		//Utilities::free_PosTypeMatrix(sub_x,sub_N,2);
		//delete[] subs;
		//delete[] sub_substructures;
	//}
//	if(stars_N > 0 && stars_implanted){
//		// std::cout << "deleting stars" << endl;
//		//delete[] star_masses;
//		star_xdisk.clear();
//		delete star_tree;
//	}
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


void LensHalo::faxial(PosType x, PosType theta,PosType f[]){
    int i,k;
    //for(int i=0;i<Nmod;i++){
        //mod1[i]=modfunc(i, bfunction(x), fratio);
        //mod1[i]=hi_order_modfunc(x, i, bfunction(x), fratio);
        //cout << "comp: " << mod1[i] << " " << hi_order_modfunc(x, i, bfunction(x), fratio) << endl;
    //}
    
    //calcModesC(dhfunction(x), fratio, pa, mod1); // used for Ansatz IV
    //calcModesB(x, fratio, pa, mod1); // used for Ansatz IIIb plus derivatives of Fourier modes
    //calcModes(fratio, beta, pa, mod1);
    //std::cout << "fratio:: " << fratio << " " << beta <<" " << pa <<  std::endl;
    //calcModes(fratio, bfunction(x), pa, mod1); // used for Ansatz IIIa w/o derivatives of Fourier modes
    
    
    //if((x<1 && x>0.99)||(bfunction(x*0.61)<1)){
    
    /*if((1.0<bfunction(x) && bfunction(x)<1.02)){
        //std::cout<< x << " " << bfunction(x) << std::endl;
        
        mod1[4]=0;
        mod1[8]=0;
        mod1[12]=0;
    }
    if((1.5<=bfunction(x) && bfunction(x)<1.51)){
        //std::cout<< x << " " << bfunction(x) << std::endl;
        
        mod1[4]=0;
        mod1[8]=0;
        mod1[12]=0;
    }*/
    
    
    f[0] = mod1[0];
    f[1] = f[2] = 0;
    for(i=4;i<Nmod;i+=2){
        k=i/2;
        f[0] +=  mod1[i]*cos(k*theta)   + mod1[i+1]*sin(k*theta);
        f[1] += -mod1[i]*k*sin(k*theta) + mod1[i+1]*k*cos(k*theta);
        f[2] += -mod1[i]*k*k*cos(k*theta) - mod1[i+1]*k*k*sin(k*theta);
    }
}


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
    f[0] = mod1[0]; // why is it commented out?
    f[1] = f[2] = 0;
    for(i=4;i<Nmod;i+=2){
        k=i/2;
        //std::cout<< mod1[i] << std::endl;
        f[0] +=  mod1[i]*cos(k*theta)   + mod1[i+1]*sin(k*theta);
        f[1] += -mod1[i]*k*sin(k*theta) + mod1[i+1]*k*cos(k*theta);
        f[2] += -mod1[i]*k*k*cos(k*theta) - mod1[i+1]*k*k*sin(k*theta);
    }
}

void LensHalo::faxial2(PosType theta,PosType f[]){
    int i,k;
    //std::cout<< mod[4] << std::endl;
    f[0] = mod2[0];
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
  double r_eps=0.3*Rsize;

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

    double r_eps=Rsize;
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
    return Utilities::nintegrate<fourier_func>(f,0.0,2*PI,1.0e-6);
}

/// Calculates potential (phi_int) from alpha_h. If flag is_alphah_a_table is True it takes and integrates directly the gfunction instead of alpha_h. The gfunction is used for the InterpolationTable used in alpha_h. Setting the flag to False speeds up the calculation of phi_h.
PosType LensHalo::alpha_int(PosType x) const{
    struct Ig_func g(*this);
    return Utilities::nintegrate<Ig_func>(g,1E-8,x,1.0e-6);
}

PosType LensHalo::norm_int(PosType r_max){
    struct norm_func g(*this,r_max);
    return Utilities::nintegrate<norm_func>(g,0.0,2.*PI,1.0e-8);
}

//PosType LensHalo::norm_intt(PosType theta){
//    struct norm_funct g(*this,theta);
//    return Utilities::nintegrate<norm_funct>(g,0.0,1.0,1.0e-8);
//}

/// Calculates the modes for fourier expansion of power law halo. All the modes are relative to the zero mode to conserve mass throughout the calculation of kappa etc.
void LensHalo::calcModes(double q, double beta, double rottheta, PosType my_mod[]){
    int i,k;
    for(int i=1;i<Nmod;++i){
		my_mod[i]=0;
	}
	// fill in modes with their values for an elliptical lens
	if(q != 1.0){
        my_mod[0] = fourier_coeff(0, q, beta)/PI/2.;
        for(i=4;i<Nmod;i+=2){
            k=i/2;
            assert(i<=Nmod);
            my_mod[i] = beta*beta*fourier_coeff(k, q, beta)/PI/(beta*beta-k*k)/my_mod[0];
        }
        my_mod[0]=1.0;
    }
	else{
		cout << "circular case" << endl;
		my_mod[0]=1.0;
	}
  
    //std::cout << "calcModes for beta=" << beta << " " << my_mod[0] << " " << my_mod[4] << " " << my_mod[8] << " " << std::endl;
    //for(int i=1;i<Nmod;++i){
    //    std::cout << i << " " << my_mod[i] << std::endl;
    //}
    // rotate model
    RotateModel(rottheta,my_mod,Nmod,0);
}

void LensHalo::calcModesC(PosType beta_r,double q, double rottheta, PosType my_mod[]){
    int i,k;
    //double beta_r=dhfunction(x); // TODO check if the dhfunction (slope of the potential) for the correct halo type is taken
    for(int i=1;i<Nmod;++i){
		my_mod[i]=0;
	}
	// fill in modes with their values for an elliptical lens
	if(q != 1.0){
        my_mod[0] = fourier_coeff(0, q, beta_r)/PI/2.;
        for(i=4;i<Nmod;i+=2){
            k=i/2;
            assert(i<=Nmod);
            my_mod[i] = fourier_coeff(k, q, beta_r)/PI/my_mod[0];
        }
        my_mod[0]=1.0;
    }
	else{
		cout << "circular case" << endl;
		my_mod[0]=1.0;
	}
    //std::cout << "calcModesC for beta=" << beta_r << " " << my_mod[0] << " " << my_mod[4] << " " << my_mod[8] << " " << std::endl;
    RotateModel(rottheta,my_mod,Nmod,0);
}



/*void LensHalo::calcModesB(PosType x,double q, double beta, double rottheta, PosType my_mod[]){ // was used for Ansatz III w derivatives of the Fourier modes
    int i,k;
    PosType dla, ddla;
    for(int i=1;i<Nmod;++i){
		my_mod[i]=0;
	}
	// fill in modes with their values for an elliptical lens
	if(q != 1.0){
        my_mod[0] = fourier_coeff(0, q, beta)/PI/2.;
        for(i=4;i<Nmod;i+=2){
            dla=dlnmod_dr(x,i, beta, q);
            ddla=ddlnmod_dr(x,i, beta, q);
            //std::cout << dla << " dla ddla " << ddla << std::endl;
            k=i/2;
            assert(i<=Nmod);
            //if(x<1){
            my_mod[i] = beta*beta*fourier_coeff(k, q, beta)/PI/(beta*beta-k*k+0*0.666*(2.*beta+1.0)*x*dla-0*log(10.0)*x*x*(ddla+dla*dla))/my_mod[0];
        
            
    //if(i==4 && x<2.0 && x>0.1 ){
      //     cout << x << " " << my_mod[i] << " "<< beta*beta-k*k << " " << dbfunction(x) << " " << ddla+dla*dla << endl;
      //}
            
        }
        my_mod[0]=1.0;
    }
	else{
		cout << "circular case" << endl;
		my_mod[0]=1.0;
	}
    
    RotateModel(rottheta,my_mod,Nmod,0);
}*/





PosType LensHalo::alpha_ell(PosType x,PosType theta){ // used only for calculation in the mass norm integral for E(A)
    PosType G, Gr,Gt, f[3],g[3];
    felliptical(x,fratio,theta,f,g);
    G = f[0];
    Gr= g[1];
    Gt = f[1];
   // std::cout << alpha_h(G) << " " << G << std::endl;
   return -alpha_h(G)/G*Gr;
   //return -alpha_h(G)/G*(pow(G/x,0.3*pow(get_slope(),4./3.)))*Gr;
}



/*void LensHalo::alphakappagamma_asym(
      PosType r         /// Radius in Mpc (not scale lengths)
      ,PosType theta    /// angle of ray
      ,PosType alpha[]  /// output deflection
      ,PosType *kappa   /// output kappa
      ,PosType gamma[]  /// outpot gamma
      ,PosType *phi
                          ){ // According to Ansatz II
  PosType f[3],g[4],alpha_r,alpha_theta;
    
  felliptical(r,fratio,theta,f,g);
  
  PosType alpha_isoG = mass*alpha_h(f[0]/rscale)/f[0]/PI;
    
  alpha_r = alpha_isoG * g[1]; // with damping
  alpha_theta =  alpha_isoG * f[1] / r; //  with damping
    
	//alpha[0] = (alpha_r*cos(theta) - alpha_theta*sin(theta))/cos(theta);
	//alpha[1] = (alpha_r*sin(theta) + alpha_theta*cos(theta))/sin(theta);

	alpha[0] = (alpha_r*cos(theta) - alpha_theta*sin(theta));
	alpha[1] = alpha_r*sin(theta) + alpha_theta*cos(theta);
    
  double kappa_isoG = mass*kappa_h(f[0]/rscale)/f[0]/f[0]/PI;
    
  double phitwo = ( 2*kappa_isoG + alpha_isoG/f[0] );
  
  *kappa = -0.5*alpha_isoG*(g[1]/r+g[2]+f[2]/r/r) + 0.5*phitwo*(g[1]*g[1]+f[1]*f[1]/r/r);
  
  double gt = -0.5*alpha_isoG*(-g[1]/r+g[2]-f[2]/r/r) + 0.5*phitwo*(g[1]*g[1]-f[1]*f[1]/r/r);
  double g45 = -alpha_isoG*(g[3]/r-f[1]/r/r) + phitwo*g[1]*f[1]/r;
  
  gamma[0] = cos(2*theta)*gt + sin(2*theta)*g45;
  gamma[1] = -sin(2*theta)*gt + cos(2*theta)*g45;

  *phi= 0; //mass*phi_h(f[0]/rscale); // <- THIS IS JUST A PLACEHOLDER
  
	return;
}*/


/** \brief Pseudo-elliptical profiles by Phi(G)-Ansatz */
void LensHalo::alphakappagamma_asym(
 PosType r         /// Radius in Mpc (not scale lengths)
 ,PosType theta    /// angle of ray
 ,PosType alpha[]  /// output deflection
 ,PosType *kappa   /// output kappa
 ,PosType gamma[]  /// outpot gamma
 ,PosType *phi
 ){ // According to Ansatz II
 PosType f[3],g[4],alpha_r,alpha_theta;
 PosType x=r/rscale;
 
 //Utilities::rotation(x,xt,theta);
  
 felliptical(x,fratio,theta+pa,f,g);
 PosType alpha_isoG = mass*alpha_h(f[0])/f[0]/PI;
 PosType kappa_isoG = mass*kappa_h(f[0])/f[0]/f[0]/PI;
 PosType phi_isoG = mass*phi_int(f[0])/PI;
 
 alpha_r = alpha_isoG * g[1]; // with damping
 alpha_theta = alpha_isoG * f[1]/x; //  with damping
 
 alpha[0] = (alpha_r*cos(theta) - alpha_theta*sin(theta));
 alpha[1] = (alpha_r*sin(theta) + alpha_theta*cos(theta));

 PosType G=f[0];
 
 /* To check for the correctness of alpha_isoG being the first derivatice of phi
  PosType h=1e-4;
  double phione=mass/PI*(phi_int(G+h)-phi_int(G-h))/2/h;
  std::cout << "I: " << alpha_isoG << " "  <<  phione << " " <<std::endl;
  */
  
  double phitwo=(2.0*kappa_isoG + alpha_isoG/G);
  *kappa = -0.5*alpha_isoG*(g[2]+g[1]/x+f[2]/x/x)+0.5*phitwo*(g[1]*g[1]+f[1]*f[1]/x/x);
 
 double gt =2.0*( -0.5*alpha_isoG*(g[2]-g[1]/x-f[2]/x/x) + 0.5*phitwo*(g[1]*g[1]-f[1]*f[1]/x/x));
 double g45 = 2.0*(-alpha_isoG*(g[3]/x-f[1]/x/x) + phitwo*g[1]*f[1]*g[1]);
 
 gamma[0] = cos(2*theta)*gt - sin(2*theta)*g45;
 gamma[1] = sin(2*theta)*gt + cos(2*theta)*g45;
 
  *phi= -1.0*phi_isoG;
 return;
 }

/** \brief Elliptical profiles by Fourier-Ansatz */
void LensHalo::alphakappagamma1asym(
       PosType r        /// Radius in Mpc (not scale lensgths)
      ,PosType theta    /// angle of ray
      ,PosType alpha[]  /// output deflection
      ,PosType *kappa   /// output kappa
      ,PosType gamma[]
      ,PosType *phi
                                    ){
  PosType f[3],g[3],alpha_r,alpha_theta;
  PosType F;
  PosType x=r/rscale;
  PosType phi_iso=mass*(phi_int(x))/PI;
  //PosType phi_iso=-1.0*mass*(phi_h(r/rscale)-log(Rsize)+1./(2-beta))/PI;    //  ( pow(x/xmax,2-beta) - 1 )/(2-beta) + log(Rsize)
  PosType alpha_iso=mass*alpha_h(x)/PI/x;  //-pow(x/xmax,-beta+1)
  PosType kappa_iso=mass*kappa_h(x)/PI/x/x;
  PosType gamma_iso=mass*gamma_h(x)/PI/x/x; // -beta*pow(x/xmax,-beta)
  
  theta=theta-PI/2.;
  faxial0(theta,f);

  //gradial(x,g);
  g[0]=1;
  g[1]=0;
  g[2]=0;
  
  
  F=f[0]-1.;
  beta=get_slope(); // only for fixed beta, i.e. PowerLaw
  //beta=bfunction(x); // only for NFW
  PosType fac=1.0/(beta*beta/(2.-beta)/(2.-beta));
  
  //alpha_r = (alpha_iso*f[0]); // w/o damping
  //alpha_theta = (phi_iso/x*f[1]); //  w/o damping
  
  alpha_r = (alpha_iso*(1+F*g[0])+phi_iso*F*g[1]); // with damping
  alpha_theta = (phi_iso*g[0]*f[1]/x); //  with damping
  
  alpha[0] = (alpha_r*cos(theta+PI/2.) - alpha_theta*sin(theta+PI/2.));
	alpha[1] = (alpha_r*sin(theta+PI/2.) + alpha_theta*cos(theta+PI/2.));
  
  //*kappa = (f[0]*kappa_iso-0.5*f[2]*phi_iso/x/x);//  w/o damping
  *kappa=(1+F*g[0])*kappa_iso-0.5*phi_iso*fac*(F*g[1]/x+F*g[2]+f[2]*g[0]/x/x)-F*g[1]*alpha_iso*x*x; /// with damping

  //PosType gt = (f[0]*gamma_iso+0.5*f[2]*phi_iso/x/x);// w/o damPIng
  //PosType g45 = ((alpha_iso*f[1])/x+(phi_iso/x/x*f[1]));// w/o damping

  
  PosType gt = (1+F*g[0])*gamma_iso+0.5*phi_iso*fac*(F*g[1]/x-F*g[2]+f[2]*g[0]/x/x)-F*g[1]*alpha_iso*x*x;// with damping
  PosType g45 = alpha_iso*f[1]*g[0]/x-phi_iso*fac*f[1]*g[1]/x+phi_iso/x/x*fac*f[1]*g[0];// with damping
  PosType tp = PI;
  gamma[0] = (cos(2.*theta-tp)*gt + sin(2.*theta-tp)*g45);
  gamma[1] = (sin(2.*theta-tp)*gt + cos(2.*theta-tp)*g45);
  
  *phi= f[0]*phi_iso; // w/o damping
  
  
  return;
}


void LensHalo::alphakappagamma2asym( // Schramm 1990
                                    PosType r        /// Radius in Mpc (not scale lensgths)
                                    ,PosType theta    /// angle of ray
                                    ,PosType alpha[]  /// output deflection
                                    ,PosType *kappa   /// output kappa
                                    ,PosType gamma[]
                                    ,PosType *phi
                                    ){
  double a = Rsize*sqrt(fratio);
  double b = Rsize/sqrt(fratio);
  double a2=a*a,b2 = b*b;
  double xm[2];
  xm[0]=r*cos(theta);
  xm[1]=r*sin(theta);
  
  //std::cout << "inside akg2: " << pa << " " << xm[0] << " " << xm[1] <<  std::endl;
  double tmp = (a2 + b2 - xm[0]*xm[0] - xm[1]*xm[1]);
  double c = cos(pa),s = sin(pa);
  double xtmp[2] = {xm[0]*c - xm[1]*s,xm[0]*s + xm[1]*c};
  //double xtmp[2] = {xm[0],xm[1]};
  
  
  double lambda = (tmp + sqrt(tmp*tmp + 4*(xm[0]*xm[0]*b2 + xm[1]*xm[1]*a2 - a2*b2 )) )/2;
  assert(lambda == lambda);
  
  PosType mo = sqrt(xtmp[0]*xtmp[0]/a2 + xtmp[1]*xtmp[1]/b2);
  //std::cout<< mo << " " << lambda << " " << xtmp[0] << " " << xtmp[1] << " " << r << " " << theta << std::endl;
  //DALPHAXDM funcX(lambda,a2,b2,xtmp,this);
  alpha[0] = -8*a*b*xtmp[0]*IDAXDM(lambda,a2,b2,xtmp,Rsize,mo)/PI;
  //alpha[0] = -8*a*b*xtmp[0]*Utilities::nintegrate<DALPHAXDM,PosType>(funcX,0,MIN(mo,1.0),1.0e-1)/PI;
  //DALPHAYDM funcY(lambda,a2,b2,xtmp,this);
  alpha[1] = -8*a*b*xtmp[1]*IDAYDM(lambda,a2,b2,xtmp,Rsize,mo)/PI;
  //alpha[1] = -8*a*b*xtmp[1]*Utilities::nintegrate<DALPHAYDM,PosType>(funcY,0,MIN(mo,1.0),1.0e-1)/PI;
  double temp = alpha[0];
  //std::cout << c << " " << s <<  std::endl;
  alpha[0] = temp*c + alpha[1]*s;
  alpha[1] = -temp*s + alpha[1]*c;
  
  double xi=sqrt(xm[0]*xm[0]+xm[1]*xm[1]/fratio/fratio);
  *kappa = mass*kappa_h(xi)/PI/xi/xi;
  
  gamma[0] = 0.0;
  gamma[1] = 0.0;
  
  *phi = 0;
  
};


void LensHalo::alphakappagamma3asym( // Keeton's 2001 adaption of Schramm 1990
                                    PosType r        /// Radius in Mpc (not scale lensgths)
                                    ,PosType theta    /// angle of ray
                                    ,PosType alpha[]  /// output deflection
                                    ,PosType *kappa   /// output kappa
                                    ,PosType gamma[]
                                    ,PosType *phi
                                    ){
  //double a = Rsize*sqrt(fratio);
  //double b = Rsize/sqrt(fratio);
  //double a2=a*a,b2 = b*b;
  double xm[2];
  xm[0]=r*cos(theta)/Rsize;
  xm[1]=r*sin(theta)/Rsize;
  
  //std::cout << "inside akg4: " << pa << " " << xm[0] << " " << Rsize <<  std::endl;
  double c = cos(pa),s = sin(pa);
  double xtmp[2] = {xm[0]*c - xm[1]*s,xm[0]*s + xm[1]*c};
  
  double j0=SCHRAMMJN(0,xtmp,Rsize);
  double j1=SCHRAMMJN(1,xtmp,Rsize);
  
  
  alpha[0] = -fratio*xtmp[0]*j0*2.75;
	alpha[1] = -fratio*xtmp[1]*j1*2.25;
  
  double temp = alpha[0];
  //std::cout << c << " " << s <<  std::endl;
  alpha[0] = temp*c + alpha[1]*s;
  alpha[1] = -temp*s + alpha[1]*c;
  
  double xi=sqrt(xtmp[0]*xtmp[0]+xtmp[1]*xtmp[1]/fratio/fratio);
  *kappa = mass*kappa_h(xi)/PI/xi/xi;
  double gt = 0.5*(2.0*fratio*xtmp[0]*xtmp[0]*SCHRAMMKN(0,xtmp,Rsize)+fratio*j0*2.75 - 2.0*fratio*xtmp[1]*xtmp[1]*SCHRAMMKN(2,xtmp,Rsize)-fratio*j1*2.25);
  double g45 = -2.0*fratio*xtmp[0]*xtmp[1]*SCHRAMMKN(1,xtmp,Rsize);
  gamma[0] = (cos(2*theta)*gt - sin(2*theta)*g45);
  gamma[1] = (sin(2*theta)*gt + cos(2*theta)*g45);
  
  *phi = -0.5*fratio*SCHRAMMI(xtmp,Rsize);
  
};

double LensHalo::SCHRAMMKN(double n, double x[], double rmax){
  struct SCHRAMMK_func f(n,x,rmax,this);
  return Utilities::nintegrate<SCHRAMMK_func>(f,0.0,1.0,1.0e-3);
}

double LensHalo::SCHRAMMJN(double n, double x[], double rmax){
  struct SCHRAMMJ_func f(n,x,rmax,this);
  return Utilities::nintegrate<SCHRAMMJ_func>(f,0.0,1.0,1.0e-3);
}

double LensHalo::SCHRAMMI(double x[], double rmax){
  struct SCHRAMMI_func f(x,rmax,this);
  return Utilities::nintegrate<SCHRAMMI_func>(f,0.0,1.0,1.0e-3);
}

double LensHalo::IDAXDM(double lambda, double a2, double b2, double x[], double rmax, double mo){
  struct IDAXDM_func f(lambda,a2,b2,x,rmax,this);
  return Utilities::nintegrate<IDAXDM_func>(f,0.0,MIN(mo,1.0),1.0e-4);
}

double LensHalo::IDAYDM(double lambda, double a2, double b2, double x[], double rmax, double mo){
  struct IDAYDM_func f(lambda,a2,b2,x,rmax,this);
  return Utilities::nintegrate<IDAYDM_func>(f,0.0,MIN(mo,1.0),1.0e-4);
}



/**
 \brief Calculate the derivatives of the G function = r*sqrt(cos(theta)^2 + q(r)^2 sin(theta))
 *
 * output: f[0] = G, f[1] = dG/dtheta, f[2] = ddg/dtheta^2, g[0] = R/r, g[1] = dG/dr
 *   , g[2] = ddG/dr^2, g[3] = ddG/dr/dtheta
 */
void LensHalo::felliptical(double r, double q, double theta, double f[], double g[]){ // According to Ansatz II
    double A[3];
    // dampening:
    //reps=rmax
    //q=r/reps+q0*(1-r/reps) // q=q0 for small radii, q=1 for large radii
    //double eps=1.-1./q/q;
    //double reps=0.5*Rsize;
    //double dampslope=2;
  
    // no damping
    A[0]=q*q; 
    A[1]=0;
    A[2]=0;
     

/*    double sig2 = r_eps*r_eps;
    double tmp = exp(-0.5*(r-Rsize)*(r-Rsize)/sig2);
    A[0] = q*q + (1-q*q)*tmp; // Compare the two cases: 1/q/q or q*q
    A[1] = -(1-q*q)*tmp*(r-Rsize)/sig2; // no r dependence thus set to 0
    A[2] = -(1-q*q)*tmp*(1 - (r-Rsize)*(r-Rsize)/sig2)/sig2;;
*/ 
    //A[0]=1-eps*(1-pow(x/reps,dampslope));
    //A[1]=eps*dampslope/reps*pow((x/reps),dampslope-1.);
    //A[2]=eps*dampslope*(dampslope-1.)/reps/reps*pow((x/reps),dampslope-2.);
  
/*    // Equations (41) - (46)
    f[0]=x*sqrt( cos(theta)*cos(theta)+A[0]*sin(theta)*sin(theta) );
    g[0]=f[0]/x;
    f[1]=x*(A[0]-1)*(cos(theta)*sin(theta))/(g[0]); // Gt
    f[2]=(x*(A[0]-1)*(pow(cos(theta),4)-A[0]*pow(sin(theta),4)))/(g[0]*g[0]*g[0]); // Gtt
    g[1]=(g[0]*g[0]+0.5*x*A[1]*sin(theta)*sin(theta))/g[0]; // Gr
    g[2]=sin(theta)*sin(theta)*(A[1]*(4*g[0]*g[0]-x*A[1]*sin(theta)*sin(theta))+2*f[0]*g[0]*A[2] )/(4*g[0]*g[0]*g[0]); // Grr
    //g[3]=cos(theta)*sin(theta)*(2*(A[0]*A[0]-1-(A[0]-1)*(A[0]-1)*cos(2*theta))+x*A[1]*(3+cos(2*theta)+2*A[0]*sin(theta)*sin(theta)))/(4*g[0]*g[0]*g[0]); // Gtr not needed here
*/
  
  double s2 = sin(theta)*sin(theta),c2 = cos(theta)*cos(theta),sc=sin(theta)*cos(theta);
  
  f[0] = r*sqrt( c2 + A[0]*s2 );
  g[0] = f[0]/r;
  f[1] = r*(A[0]-1.)*sc/g[0]; // Gt
  f[2] = r*(A[0]-1.)*( c2 - s2 - (A[0]-1)*s2*c2/g[0]/g[0] )/g[0]; // Gtt
  g[1] = g[0] + 0.5*r*A[1]*s2/g[0]; // Gr
  g[2] = 0.5*s2*( A[1]*(2 - 0.5*r*s2*A[1]/g[0]/g[0]) + r*A[2] )/g[0];  // Grr
  g[3] = sc*( A[0]-1 + r*A[1] - 0.5*r*(A[0]-1)*A[1]*s2/g[0]/g[0] )/g[0]; //Grt
}

PosType LensHalo::renormalization(PosType r_max){ // Calculates renormalization factor in the constructor of PowerLaw and NFW only (for now)
  double fac=1;
  return norm_int(fac*r_max)*fac*r_max/2.0/PI/(-1.0*alpha_h(fac*r_max)/r_max/fac);
}

PosType LensHalo::kappa_asym(PosType r, PosType theta){ /// radius in Mpc (Not x = r/rscale)
  PosType F, f[3],g[3], kappa, x;
  x=r;
  PosType phi=phi_int(x);
  faxial(x,theta,f); // for this to work calculate modes in faxial with calcModesB !!
  gradial(x,g);
  F=f[0]-1;
  beta=get_slope(); // only for fixed beta, i.e. PowerLaw
  //beta=bfunction(x); // only for NFW
  double fac=1.0/(beta*beta/(2.-beta)/(2.-beta));
  kappa=f[0]*kappa_h(x)-0.5*f[2]*fac*phi;//  w/o damping
  //kappa=(1+F*g[0])*kappa_h(x)-0.5*phi*fac*(F*g[1]/x+F*g[2]+f[2]*g[0]/x/x)*x*x-F*g[1]*alpha_h(x)*x*x; /// with damping
  return kappa;
}




/* Phi(G(r)) ansatz (IV in notes)  with r-dependent modes a_n(r) [requires renormalization?]
PosType LensHalo::renormalization(PosType r_max){ // Calculates renormalization factor in the constructor of PowerLaw and NFW only (for now)
  double fac=1;
  return norm_int(fac*r_max)*fac*r_max/2.0/PI/(-1.0*alpha_h(fac*r_max)/r_max/fac);
}

PosType LensHalo::kappa_asym(PosType r, PosType theta){ /// radius in Mpc (Not x = r/rscale)
  // According to Ansatz II
  double f[3],g[3];
  //double kappa;
  felliptical(r,fratio,theta,f,g);
  double alpha_isoG = mass*alpha_h(f[0]/rscale)/f[0]/PI;
  double kappa_isoG = mass*kappa_h(f[0]/rscale)/f[0]/f[0]/PI;
  
  // double b=dbfunction(G);
  // std::cout<< b << std::endl;
  // double fac=1.0/(b*b/(2.-b)/(2.-b));
  
  // we agreed on the following to be wrong - this was an attempt to express phitwo in terms of iso kappa and alpha
  //double phitwo=(2.0*kappa_h(f[0]/rscale) - alpha_h(f[0]/rscale)/f[0]/f[0]*(g[1]/r+g[2]))/g[1]/g[1];
  
  double phitwo = ( 2*kappa_isoG + alpha_isoG/f[0] );
  
  
  //double phitwo=(2.0*kappa_h(G) - alpha_h(G)/G*(Gr/x+Grr))/Gr/Gr;
  // double phitwo = (2.0*kappa_h(x)+alpha_h(x)/x/x+alpha_h(G)/G*Grr)/Gr/Gr;
  
  //double phitwo=-1.0*ddhfunction(G,true); // correct phi'' for NFW, the flag distinguishes numerical calculation (true) from an analytic one (false), caution the false-case might contain an error.
  
  
  //double beta=get_slope(); // needed for PowerLawHalos
  //double phitwo = mass*(1-beta)*pow(f[0]/Rsize,-beta)/PI/Rsize/Rsize; // works for PowerLawHalos
  
  //phitwo=dgfunctiondx(G); // used for a check of values, leave it here for now.
  
  return -0.5*alpha_isoG*(g[1]/r+g[2]+f[2]/r/r) + 0.5*phitwo*(g[1]*g[1]+f[1]*f[1]/r/r); //#print1
  //return 0.5*kappa*x*x*(2.-beta) ; ///pow(mnorm,4./3.); // if A(r)=1/q/q  E(A)^(4/3) gives sort-of the correct normalization, i.e mnorm^(4/3.)  for now I tool the normalization out again.
}

*/


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



// Ansatz IV, issues a curl in kappa at larger radii.
/*
PosType LensHalo::kappa_asym(PosType x,PosType theta){
	PosType F, f[4],g[3], kappa;
    PosType amod[3],amodp1[3];
    
    PosType imod, imodp1, dmod_db, ddmod_db, dmod_dq, ddmod_dq, dbdr, ddbdr, dqdr, ddqdr;
    PosType dmod_dbp1, ddmod_dbp1, dmod_dqp1, ddmod_dqp1;
    
    
    PosType phi=phi_int(x);
    
 
    // the following for loop is used to check if InterpolateModes yields the correct Fourier modes
     //for(int l=1; l<99; l++ ){
        //PosType q=0.01*float(l);
        //for(int i = 1; i < 199 ; i++){
            //calcModesC(i*0.01, q, pa, mod1);
            //std::cout << i*0.01 << " " <<  mod1[4] << " " << mod1[8] << " " << mod1[12] << " " << mod1[16] << " " << mod1[20] << " " << mod1[24] << " " << mod1[28] << " " << mod1[32] << " " << q << std::endl;
            //std::cout << " " <<  InterpolateModes(4, q+0.005, 1.005) << " " <<  std::endl;
            //std::cout << i*0.01 << " " <<  InterpolateModes(4, q, i*0.01) << " " << InterpolateModes(8, q, i*0.01) << " " << InterpolateModes(12, q, i*0.01) << " " <<InterpolateModes(16, q, i*0.01) << " " << InterpolateModes(20, q, i*0.01) << " " << InterpolateModes(24, q, i*0.01) << " " << InterpolateModes(28, q, i*0.01) << " " << InterpolateModes(32, q, i*0.01) << " " << q << std::endl;
     // }
     //}



    PosType b=dhfunction(x);
    for(int i=4;i<Nmod;i+=2){
        int k=i/2;
        imod=InterpolateModes(i, fratio, b);
        imodp1=0.0;//InterpolateModes(i+1, fratio, b);
        dmod_db=dmoddb(i, fratio, b);
        dmod_dbp1=0.0;//dmoddb(i+1, fratio, b);
        ddmod_db=ddmoddb(i, fratio, b);
        ddmod_dbp1=0.0;//ddmoddb(i+1, fratio, b);
        dmod_dq=dmoddq(i, fratio, b);
        dmod_dqp1=0.0;//dmoddq(i+1, fratio, b);
        ddmod_dq=ddmoddq(i, fratio, b);
        ddmod_dqp1=0.0; //ddmoddq(i+1, fratio, b);
        dbdr=ddhfunction(x,true);
        ddbdr=dddhfunction(x,true);
        dqdr=0;
        ddqdr=0;
        amod[0]=imod;
        amodp1[0]=0.0; //imodp1;
        amod[1]=dmod_db*dbdr+dmod_dq*dqdr;
        amodp1[1]=0.0; //dmod_dbp1*dbdr+dmod_dqp1*dqdr;
        amod[2]=ddmod_db*dbdr*dbdr+dmod_db*ddbdr+ddmod_dq*dqdr*dqdr+dmod_dq*ddqdr;
        amodp1[2]=0.0;// ddmod_dbp1*dbdr*dbdr+dmod_dbp1*ddbdr+ddmod_dqp1*dqdr*dqdr+dmod_dqp1*ddqdr;
        PosType coskt=cos(k*theta);
        PosType sinkt=sin(k*theta);
        PosType aux=amod[0]*coskt+amodp1[0]*sinkt;
        f[0] +=  aux;
        f[1] +=  amod[1]*coskt + amodp1[1]*sinkt;
        f[2] +=  amod[2]*coskt + amodp1[2]*sinkt;
        f[3] += -k*k*aux;
    }
    
    double bk=dhfunction(x);
    double fac=0.5/(bk*bk/(2.-bk)/(2.-bk));
    
    kappa=kappa_h(x);
    for(int i=4;i<Nmod;i+=2){
        int k=i/2;
        amod[0]=InterpolateModes(i, fratio, b);
        dmod_db=dmoddb(i, fratio, b);
        dbdr=ddhfunction(x,true);
        amod[1]=dmod_db*dbdr;
        ddmod_db=ddmoddb(i, fratio, b);
        ddbdr=dddhfunction(x,true);
        amod[2]=ddmod_db*dbdr*dbdr+dmod_db*ddbdr;
        kappa+=(amod[0]*(kappa_h(x)+k*k/2*phi*fac)+amod[1]*(-1.0*alpha_h(x)*x*x-phi*fac*x/2)-amod[2]*phi*fac*x*x/2)*cos(k*theta);
    }
 
    return kappa;
}
*/



// NOTE THAT GAMMA_ASYM does not use Ansatz II yet!
// In gamma_asym phi_int(x) is used, which calculates the potential unlike phi_h from alpha_h.

void LensHalo::gamma_asym(PosType x,PosType theta, PosType gamma[]){
	PosType F, f[3],g[3];
  PosType phi=phi_int(x);
  
  
  faxial0(theta,f);
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



