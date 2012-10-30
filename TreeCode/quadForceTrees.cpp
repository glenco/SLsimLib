/*
 * forceTree.cpp
 *
 *  Created on: Oct 18, 2011
 *      Author: bmetcalf
 */
#include <assert.h>
#include <math.h>
#include <slsimlib.h>
//#include <quadTree.h>
#include <iostream>


QuadTreePowerLaw::QuadTreePowerLaw(
		float beta                 /// slop of mass profile \f$ \Sigma \propto r^\beta \f$
		,PosType **xp              /// positions of the halos xp[0..Npoints-1][0..1 or 2]
		,IndexType Npoints         /// number of halos
		,HaloStructure *h_params   /// array with internal properties of halos
		,double my_kappa_bk        /// Background convergence to be subtracted
		,int bucket                /// maximum number of halos in a leaf of the tree
		,PosType theta             /// Opening angle used in tree force calculation, default 0.1
		) :
		QuadTree(xp,h_params,Npoints,my_kappa_bk,bucket,theta), beta(beta)
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
		,double my_kappa_bk        /// Background convergence to be subtracted
		,int bucket                /// maximum number of halos in a leaf of the tree
		,PosType theta             /// Opening angle used in tree force calculation, default 0.1
		) :
		QuadTree(xp,h_params,Npoints,my_kappa_bk,bucket,theta)
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
		,double my_kappa_bk       /// Background convergence to be subtracted
		,int bucket                /// maximum number of halos in a leaf of the tree
		,PosType theta             /// Opening angle used in tree force calculation, default 0.1
		) :
		QuadTree(xp,h_params,Npoints,my_kappa_bk,bucket,theta), beta(my_beta)
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
		,double my_kappa_bk       /// Background convergence to be subtracted
		,int bucket                /// maximum number of halos in a leaf of the tree
		,PosType theta             /// Opening angle used in tree force calculation, default 0.1
		) :
		QuadTree(xp,h_params,Npoints,my_kappa_bk,bucket,theta)
{

	for(long i=0;i<Npoints;++i){
		h_params[i].Rsize = rmaxNSIE(h_params[i].sigma,h_params[i].mass,h_params[i].fratio,h_params[i].rscale);
	}
}
QuadTreeNSIE::~QuadTreeNSIE(){
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

	double ellipR = ellipticRadiusNSIE(xcm,halo_params.fratio,halo_params.pa);
	if(ellipR > halo_params.Rsize){
		// if the ray misses the halo treat it as a point mass
		  double prefac = -1.0*halo_params.mass/rcm2/pi;

		  alpha[0] += prefac*xcm[0];
		  alpha[1] += prefac*xcm[1];

		  // can turn off kappa and gamma calculations to save times
		  if(!no_kappa){
			  prefac *= 2.0/rcm2;

			  gamma[0] += 0.5*(xcm[0]*xcm[0]-xcm[1]*xcm[1])*prefac;
			  gamma[1] += xcm[0]*xcm[1]*prefac;
		  }
	}else{
		double xt[2]={0,0},tmp[2]={0,0};
		float units = 0.5*pow(halo_params.sigma/lightspeed,2)/Grav; // mass/distance(physical)
		xt[0]=xcm[0];
		xt[1]=xcm[1];
		alphaNSIE(tmp,xt,halo_params.fratio,halo_params.rscale,halo_params.pa);
		alpha[0] -= units*tmp[0];
		alpha[1] -= units*tmp[1];
		if(!no_kappa){
			float tmp[2]={0,0};
			*kappa += units*kappaNSIE(xt,halo_params.fratio,halo_params.rscale,halo_params.pa);
			gammaNSIE(tmp,xt,halo_params.fratio,halo_params.rscale,halo_params.pa);
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
		,double my_kappa_bk       /// Background convergence to be subtracted
		,int bucket                /// maximum number of halos in a leaf of the tree
		,PosType theta             /// Opening angle used in tree force calculation, default 0.1
		) :
		QuadTreeNFW(xp,Npoints,h_params,my_kappa_bk,bucket,theta),
		QuadTreeNSIE(xp,Npoints,h_params,my_kappa_bk,bucket,theta)
{

}
QuadTreeNFW_NSIE::~QuadTreeNFW_NSIE(){
}
void QuadTreeNFW_NSIE::force_halo(
		double *alpha
		,float *kappa
		,float *gamma
		,double *xcm
	  	,HaloStructure &halo_params
	  	,bool no_kappa
	  	){
	float tmp_kappa,tmp_gamma[2];
	double tmp_alpha[2];

	QuadTreeNSIE::force_halo(tmp_alpha,&tmp_kappa,tmp_gamma,xcm,halo_parms,no_kappa);
	QuadTreeNFW::force_halo(alpha,kappa,gamma,xcm,halo_parms,no_kappa);

	alpha[0] += tmp_alpha[0];
	alpha[1] += tmp_alpha[1];

	gamma[0] += tmp_gamma[0];
	gamma[1] += tmp_gamma[1];

	*kappa += tmp_kappa;
}
