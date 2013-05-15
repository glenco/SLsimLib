/*
 * singlelens.cpp
 *
 *  Created on: May 10, 2013
 *      Author: mpetkova
 */

#include "singlelens.h"

SingleLens::SingleLens(InputParams& params) : Lens(){
	assignParams(params);

	switch(DM_halo_type){
	case PowerLaw:
		halo.push_back(new PowerLawLensHalo(params));
		break;
	case NFW:
		halo.push_back(new NFWLensHalo(params));
		break;
	case PseudoNFW:
		halo.push_back(new PseudoNFWLensHalo(params));
		break;
	case NSIE:
		halo.push_back(new BaseNSIELensHalo(params));
		break;
	case AnaNSIE:
		halo.push_back(new AnaNSIELensHalo(params));
		break;
	case UniNSIE:
		halo.push_back(new UniNSIELensHalo(params));
		break;
	case PointMass:
		halo.push_back(new LensHalo(params));
		break;
	default:
		ERROR_MESSAGE();
		std::cout << "Incorrect halo_type selected! Please choose from:" << std::endl;
		std::cout << "0: PowerLaw, 1: NFW, 2: PseudoNFW, 3: NSIE, 4: AnaNSIE, 5: UniNSIE, 6: PointMass" << std::endl;
		exit(1);
		break;
	}

	if(Nprof == 2){
		switch(galaxy_halo_type){
		case PowerLaw:
			halo.push_back(new PowerLawLensHalo(params));
			break;
		case NFW:
			halo.push_back(new NFWLensHalo(params));
			break;
		case PseudoNFW:
			halo.push_back(new PseudoNFWLensHalo(params));
			break;
		case NSIE:
			halo.push_back(new BaseNSIELensHalo(params));
			break;
		case AnaNSIE:
			halo.push_back(new AnaNSIELensHalo(params));
			break;
		case UniNSIE:
			halo.push_back(new UniNSIELensHalo(params));
			break;
		case PointMass:
			halo.push_back(new LensHalo(params));
			break;
		default:
			ERROR_MESSAGE();
			std::cout << "Incorrect halo_type selected! Please choose from:" << std::endl;
			std::cout << "0: PowerLaw, 1: NFW, 2: PseudoNFW, 3: NSIE, 4: AnaNSIE, 5: UniNSIE, 6: PointMass" << std::endl;
			exit(1);
			break;
		}
	}
}


SingleLens::~SingleLens(){
	std::cout << "deleting SingleLens" << std::endl;
}

void SingleLens::assignParams(InputParams& params){
	if(!params.get("DM_halo_type",DM_halo_type)) error_message1("DM_halo_type",params.filename());
	if(!params.get("DM_halo_type",galaxy_halo_type)){
		Nprof = 1;
	}
	else{
		Nprof = 2;
	}
}

void SingleLens::error_message1(std::string parameter,std::string file){
	ERROR_MESSAGE();
	std::cout << "Parameter " << parameter << " is needed to construct a SingleLens.  It needs to be set in parameter file " <<
			file << "!" << std::endl;
	exit(0);
}

double SingleLens::getZlens(){
	return halo[0]->get_zlens();
}

/// resets Zl, Dl, Sigma_crit, MpcToAsec
void SingleLens::setZlens(CosmoHndl cosmo,double zl,double zsource){
	halo->set_zlens(zl);
	halo->setInternalParams(cosmo,zsource);
}

/// Sets parameters within BaseLens that depend on the source redshift - Dl,Sigma_crit,etc.
void SingleLens::setInternalParams(CosmoHndl cosmo, SourceHndl source){
	setInternalParams(cosmo,source->getZ());
}

void SingleLens::setInternalParams(CosmoHndl cosmo, double zsource){
	halo->setInternalParams(cosmo,zsource);
}
