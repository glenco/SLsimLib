/*
 * forceTree.cpp
 *
 *  Created on: Oct 18, 2011
 *      Author: bmetcalf
 */
#include <assert.h>
#include <math.h>
#include <forceTree.h>
#include <quadTree.h>
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
}

QuadTreeNFW::~QuadTreeNFW(){
}

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
}

QuadTreePseudoNFW::~QuadTreePseudoNFW(){
}
