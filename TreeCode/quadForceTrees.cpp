/*
 * forceTree.cpp
 *
 *  Created on: Oct 18, 2011
 *      Author: bmetcalf
 */

#include "slsimlib.h"

QuadTreePowerLaw::QuadTreePowerLaw(
		float beta                 /// slop of mass profile \f$ \Sigma \propto r^\beta \f$
		,PosType **xp              /// positions of the halos xp[0..Npoints-1][0..1 or 2]
		,IndexType Npoints         /// number of halos
		,HaloStructure *h_params   /// array with internal properties of halos
		,double my_sigma_bk        /// Background convergence to be subtracted
		,int bucket                /// maximum number of halos in a leaf of the tree
		,PosType theta             /// Opening angle used in tree force calculation, default 0.1
		) :
		QuadTree(xp,h_params,Npoints,my_sigma_bk,bucket,theta), beta(beta)
{

}

QuadTreePowerLaw::~QuadTreePowerLaw(){
}

long QuadTreeNFW::ob_count = 0;
double *QuadTreeNFW::ftable = NULL,*QuadTreeNFW::gtable = NULL,*QuadTreeNFW::g2table = NULL,*QuadTreeNFW::xtable = NULL;

QuadTreeNFW::QuadTreeNFW(
		PosType **xp               /// positions of the halos xp[0..Npoints-1][0..1 or 2]
		,IndexType Npoints         /// number of halos
		,HaloStructure *h_params   /// array with internal properties of halos
		,double my_sigma_bk        /// Background convergence to be subtracted
		,int bucket                /// maximum number of halos in a leaf of the tree
		,PosType theta             /// Opening angle used in tree force calculation, default 0.1
		) :
		QuadTree(xp,h_params,Npoints,my_sigma_bk,bucket,theta)
{
	for(unsigned long i=0;i<Npoints;++i){
		if(h_params[i].Rmax <= 0.0 || h_params[i].rscale <= 0.0){
			ERROR_MESSAGE();
			std::cout << "Illegal values for halo internal valuables." << std::endl;
			exit(1);
		}
	}

	// make halo profile lookup tables if this is the first instance of a QuadTreeNFW
	if(ob_count == 0) make_tables();
	++ob_count;
}

QuadTreeNFW::~QuadTreeNFW(){
	--ob_count;
	if(ob_count == 0){
		delete[] xtable;
		delete[] gtable;
		delete[] ftable;
		delete[] g2table;
	}
}

long QuadTreePseudoNFW::ob_count = 0;
double * QuadTreePseudoNFW::mhattable = NULL,*QuadTreePseudoNFW::xtable = NULL;

QuadTreePseudoNFW::QuadTreePseudoNFW(
		double my_beta                 /// outer slope of profile is \f$ \Sigma \propto r^{-\beta} \f$
		,PosType **xp              /// positions of the halos xp[0..Npoints-1][0..1 or 2]
		,IndexType Npoints         /// number of halos
		,HaloStructure *h_params   /// array with internal properties of halos
		,double my_sigma_bk       /// Background convergence to be subtracted
		,int bucket                /// maximum number of halos in a leaf of the tree
		,PosType theta             /// Opening angle used in tree force calculation, default 0.1
		) :
		QuadTree(xp,h_params,Npoints,my_sigma_bk,bucket,theta), beta(my_beta)
{

	if(beta == 0.0){
		std::cout << "QuadTreePseudoNFW: The slope can not be zero!" << std::endl;
		exit(1);
	}

	// Check for values that would make the rayshooter return nan.
	for(unsigned long i=0;i<Npoints;++i){
		if(h_params[i].Rmax <= 0.0 || h_params[i].rscale <= 0.0){
			ERROR_MESSAGE();
			std::cout << "Illegal values for halo internal valuables." << std::endl;
			exit(1);
		}
	}

	// make halo profile lookup tables if this is the first instance of a QuadTreePseudoNFW
	if(ob_count == 0) make_tables();
	++ob_count;
}

QuadTreePseudoNFW::~QuadTreePseudoNFW(){
	--ob_count;
	if(ob_count == 0){
		delete[] xtable;
		delete[] mhattable;
	}
}

QuadTreeNSIE::QuadTreeNSIE(
		PosType **xp              /// positions of the halos xp[0..Npoints-1][0..1 or 2]
		,IndexType Npoints         /// number of halos
		,HaloStructure *h_params   /// array with internal properties of halos
		,double my_sigma_bk       /// Background convergence to be subtracted
		,int bucket                /// maximum number of halos in a leaf of the tree
		,PosType theta             /// Opening angle used in tree force calculation, default 0.1
		) :
		QuadTree(xp,h_params,Npoints,my_sigma_bk,bucket,theta,true)
{

	//for(unsigned long i=0;i<Npoints;++i){
	//	h_params[i].Rsize_nsie = rmaxNSIE(h_params[i].sigma_nsie,h_params[i].mass,h_params[i].fratio_nsie,h_params[i].rscale);
	//}
}
QuadTreeNSIE::~QuadTreeNSIE(){
}


/// Testing function
void QuadTreeNSIE::test_force_halo(
	  	HaloStructure &halo_params
		){
	double alpha[2];
	float kappa, gamma[2];
	double xcm[2];
	int N = 200;

	std::cout << N << std::endl;

	for(int i = 0; i<N ; ++i){
		xcm[0] = 2.5*halo_params.Rsize_nsie*MAX(1.0,1.0/halo_params.fratio_nsie)*(i*1.0/(N-1) - 0.5);
		for(int j = 0; j<N ; ++j){
			xcm[1] = 2.5*halo_params.Rsize_nsie*MAX(1.0,1.0/halo_params.fratio_nsie)
			*(j*1.0/(N-1) - 0.5);

			alpha[0]=alpha[1] = 0.0;
			gamma[0]=gamma[1] = 0.0;
			force_halo(alpha,&kappa,gamma,xcm,halo_params,true);

			std::cout << xcm[0] << " " << xcm[1] << "     " << alpha[0] << "  " << alpha[1] << std::endl;
		}
	}

}

/// This is the function that override QuadTree::force_halo to make it 2D.
void QuadTreeNSIE::force_halo(
		double *alpha
		,float *kappa
		,float *gamma
		,double *xcm
	  	,HaloStructure &halo_params
	  	,bool no_kappa
	  	){

	double rcm2 = xcm[0]*xcm[0] + xcm[1]*xcm[1];
	if(rcm2 < 1e-20) rcm2 = 1e-20;

	double ellipR = ellipticRadiusNSIE(xcm,halo_params.fratio_nsie,halo_params.pa_nsie);
	if(ellipR > halo_params.Rsize_nsie){
		double rout = halo_params.Rsize_nsie*MAX(1.0,1.0/halo_params.fratio_nsie);
		// This is the case when the ray is within the NSIE's circular region of influence but outside its elliptical truncation
		  // if the ray misses the halo treat it as a point mass
		double prefac = -1.0*halo_params.mass_nsie/rcm2/pi;

		if(rcm2 > rout*rout){
			alpha[0] += prefac*xcm[0];
			alpha[1] += prefac*xcm[1];
		}else{
			double alpha_out[2],alpha_in[2],rin,x_in[2];
			double prefac = -1.0*halo_params.mass_nsie/rout/pi;
			double r = sqrt(rcm2);

			alpha_out[0] = prefac*xcm[0]/r;
			alpha_out[1] = prefac*xcm[1]/r;

			rotation(x_in,xcm,halo_params.pa_nsie);
			rin = r*halo_params.Rsize_nsie
					  /sqrt( x_in[0]*x_in[0] + pow(halo_params.fratio_nsie*x_in[1],2) );
			//rin = halo_params.Rsize_nsie;

			x_in[0] = rin*xcm[0]/r;
			x_in[1] = rin*xcm[1]/r;

			alpha_in[0] = alpha_in[1] = 0;
			float units = pow(halo_params.sigma_nsie/lightspeed,2)/Grav/sqrt(halo_params.fratio_nsie); // mass/distance(physical)
			alphaNSIE(alpha_in,x_in,halo_params.fratio_nsie,halo_params.rcore_nsie,halo_params.pa_nsie);
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
		float units = pow(halo_params.sigma_nsie/lightspeed,2)/Grav/sqrt(halo_params.fratio_nsie); // mass/distance(physical)
		xt[0]=xcm[0];
		xt[1]=xcm[1];
		alphaNSIE(tmp,xt,halo_params.fratio_nsie,halo_params.rcore_nsie,halo_params.pa_nsie);
		alpha[0] -= units*tmp[0];
		alpha[1] -= units*tmp[1];
		if(!no_kappa){
			float tmp[2]={0,0};
			*kappa += units*kappaNSIE(xt,halo_params.fratio_nsie,halo_params.rcore_nsie,halo_params.pa_nsie);
			gammaNSIE(tmp,xt,halo_params.fratio_nsie,halo_params.rcore_nsie,halo_params.pa_nsie);
			gamma[0] += units*tmp[0];
			gamma[1] += units*tmp[1];
		}
	}
	return;
}

QuadTreeNFW_NSIE::QuadTreeNFW_NSIE(
		PosType **xp              /// positions of the halos xp[0..Npoints-1][0..1 or 2]
		,IndexType Npoints         /// number of halos
		,HaloStructure *h_params   /// array with internal properties of halos
		,double my_sigma_bk       /// Background convergence to be subtracted
		,int bucket                /// maximum number of halos in a leaf of the tree
		,PosType theta             /// Opening angle used in tree force calculation, default 0.1
		) :
		QuadTreeNFW(xp,Npoints,h_params,my_sigma_bk,bucket,theta)
{
	// The background should be taken care of in the QuadTreeNFW part
	qtreensie = new QuadTreeNSIE(xp,Npoints,h_params,0.0,bucket,theta);
}
QuadTreeNFW_NSIE::~QuadTreeNFW_NSIE(){
	delete qtreensie;
}

//virtual void force2D(double *ray,double *alpha,float *kappa,float *gamma,bool no_kappa);
//virtual void force2D_recur();

void QuadTreeNFW_NSIE::force2D_recur(
		double *ray
		,double *alpha
		,float *kappa
		,float *gamma
		,bool no_kappa
		){

	float tmp_kappa,tmp_gamma[2];
	double tmp_alpha[2];

	qtreensie->force2D_recur(ray,tmp_alpha,&tmp_kappa,tmp_gamma,no_kappa);
	QuadTree::force2D_recur(ray,alpha,kappa,gamma,no_kappa);

	alpha[0] += tmp_alpha[0];
	alpha[1] += tmp_alpha[1];

	gamma[0] += tmp_gamma[0];
	gamma[1] += tmp_gamma[1];

	*kappa += tmp_kappa;
}
void QuadTreeNFW_NSIE::force2D(
		double *ray
		,double *alpha
		,float *kappa
		,float *gamma
		,bool no_kappa
		){

	float tmp_kappa,tmp_gamma[2];
	double tmp_alpha[2];

	qtreensie->force2D(ray,tmp_alpha,&tmp_kappa,tmp_gamma,no_kappa);
	QuadTree::force2D(ray,alpha,kappa,gamma,no_kappa);

	alpha[0] += tmp_alpha[0];
	alpha[1] += tmp_alpha[1];

	gamma[0] += tmp_gamma[0];
	gamma[1] += tmp_gamma[1];

	*kappa += tmp_kappa;
}
