/*
 * MOKAlens.h
 *
 *  Created on: Jun 8, 2012
 *      Author: mpetkova
 */


#ifndef MOKALENS_H_
#define MOKALENS_H_

#include "standard.h"
#include "profile.h"
#include "InputParams.h"
#include "lens_halos.h"
#include "model.h"

/**
 * \brief The MOKA map structure, containing all quantities that define it
 *
 * The MOKA map, that is read in from the fits file. Its components include the
 * lensing properties, as well as the cosmology, the size of the field of view,
 * the redshifts of the lens and source, the properties of the cluster, etc.
 *
 * Note: To use this class requires setting the ENABLE_FITS compiler flag and linking
 * the cfits library.
 */
struct MOKAmap{
	/// values for the map
	std::valarray<float> convergence;
	std::valarray<float> alpha1;
	std::valarray<float> alpha2;
	std::valarray<float> gamma1;
	std::valarray<float> gamma2;
	std::valarray<float> gamma3;
	std::valarray<float> Signlambdar;
	std::valarray<float> Signlambdat;
	std:: vector<double> x;	 
    int nx,ny;
    // boxlMpc is Mpc/h for MOKA
	/// lens and source properties
    double zlens,m,zsource,DL,DLS,DS,c,cS,fsub,mstar,minsubmass;
    double boxlarcsec,boxlMpc,boxlrad;
    /// cosmology
    double omegam,omegal,h,wq;
	double inarcsec;
	double center[2];
};

/**
 *  \brief A class that includes the MOKA lens map
 *
 * A class, where the lens is represented by a MOKA cluster in the form of a
 * MOKAmap object (see the description of MOKAmap). It can either be used with the
 * Lens model or by itself.
 *
 * Note: To use this class requires setting the ENABLE_FITS compiler flag and linking
 * the cfits library.
 */


class LensHaloMOKA : public LensHalo{
public:

	LensHaloMOKA(InputParams& params);

	~LensHaloMOKA();

	bool set;	/// the name of the MOKA input file
	std::string MOKA_input_file;
	/// if >=1 (true), do analyzis only; if = 0 (false) change units to internal GLAMER units and prepare for ray-shooting
	int flag_MOKA_analyze;
	int flag_background_field;

	void assignParams(InputParams& params);
	void setup(CosmoHndl,SourceHndl);
	void saveImage(bool saveprofile=true);
	void saveKappaProfile();
	void saveGammaProfile();
	void saveProfiles(double &RE3, double &xxc, double &yyc);
	void initMap();
	void force_halo(double *alpha,KappaType *kappa,KappaType *gamma,double *xcm,bool no_kappa,bool subtract_point=false);

	MOKAmap *map;

	void estSignLambdas();
	void EinsteinRadii(double &RE1, double &RE2, double &xxc, double &yyc);

	void getDims();
	void readImage();
	void writeImage(std::string fn);

};

/**
 * \brief The MOKA model that contains the lens, source and cosmology.
 * 
 * This is a Model subclass which reads in a MOKA file and updates the
 * cosmology and lens according to the values read from the FITS file.
 */
template<typename SourceT = SourceUniform>
class ModelMOKA : public Model<SourceT>
{
public:
	ModelMOKA(std::string paramfile, long* seed) : Model<SourceT>(paramfile, seed)
	{
		// make sure the lens contains a LensHaloMOKA
		if(this->lens->main_halos.template size<LensHaloMOKA>() != 1)
			throw std::runtime_error("ModelMOKA needs to contain exactly one main halo of type LensHaloMOKA!");
		
		// update cosmology and source from halo
		this->lens->main_halos.template get<LensHaloMOKA>(0)->setup(this->cosmo, this->source);
	}
	
private:
	LensHaloMOKA* halo;
};

void make_friendship(int ii,int ji,int np,std:: vector<int> &friends, std:: vector<double> &pointdist);

int fof(double l,std:: vector<double> xci, std:: vector<double> yci, std:: vector<int> &groupid);
#endif /* MOKALENS_H_ */

