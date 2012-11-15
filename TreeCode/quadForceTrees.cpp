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

	// TODO BEN Worry about background subtraction so that it isn't redundant.

	qtreensie->force2D(ray,tmp_alpha,&tmp_kappa,tmp_gamma,no_kappa);
	QuadTree::force2D(ray,alpha,kappa,gamma,no_kappa);

	alpha[0] += tmp_alpha[0];
	alpha[1] += tmp_alpha[1];

	gamma[0] += tmp_gamma[0];
	gamma[1] += tmp_gamma[1];

	*kappa += tmp_kappa;
}
