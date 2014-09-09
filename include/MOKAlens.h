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
#include "grid_maintenance.h"

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
	std::valarray<double> convergence;
	std::valarray<double> alpha1;
	std::valarray<double> alpha2;
	std::valarray<double> gamma1;
	std::valarray<double> gamma2;
	std::valarray<double> gamma3;
  std::valarray<double> phi;
	std::valarray<double> Signlambdar;
	std::valarray<double> Signlambdat;
	std:: vector<double> x;
  int nx,ny;
  // boxlMpc is Mpc/h for MOKA
	/// lens and source properties
  double zlens,m,zsource,Dlens,DLS,DS,c,cS,fsub,mstar,minsubmass;
  /// range in x direction, pixels are square
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
class LensHaloMOKA : public LensHalo
{
public:
	LensHaloMOKA(const std::string& filename, LensHaloType maptype, const COSMOLOGY& lenscosmo);
	LensHaloMOKA(InputParams& params, const COSMOLOGY& lenscosmo);
	
	~LensHaloMOKA();
	
	std::string MOKA_input_file;
	/// if >=1 (true), do analyzis only; if = 0 (false) change units to internal GLAMER units and prepare for ray-shooting
	int flag_MOKA_analyze;
	int flag_background_field;
	
	void assignParams(InputParams& params);

	void checkCosmology();
	
	void saveImage(bool saveprofile=true);
	void saveKappaProfile();
	void saveGammaProfile();
	void saveProfiles(double &RE3, double &xxc, double &yyc);
    
	void force_halo(double *alpha,KappaType *kappa,KappaType *gamma,KappaType *phi,double const *xcm,bool subtract_point=false);
    
	void saveImage(GridHndl grid,bool saveprofiles);
	
	void estSignLambdas();
	void EinsteinRadii(double &RE1, double &RE2, double &xxc, double &yyc);
	
	void getDims();
	void readImage();
	void writeImage(std::string fn);
  
	/// return center in physical Mpc
	const double* getCenter() const { return map->center; }
	/// return range of input map in rad
	double getRangeRad() const { return map->boxlrad; }
	/// return range of input map in physical Mpc
	double getRangeMpc() const { return map->boxlMpc; }
	/// return number of pixels on a side in original map
	size_t getN() const { return map->nx; }
	 
private:
	LensHaloType maptype;
	void initMap();
	void convertmap(MOKAmap *map,LensHaloType maptype);
	MOKAmap* map;
	const COSMOLOGY& cosmo;
	void PreProcessFFTWMap();
};

void make_friendship(int ii,int ji,int np,std:: vector<int> &friends, std:: vector<double> &pointdist);

int fof(double l,std:: vector<double> xci, std:: vector<double> yci, std:: vector<int> &groupid);
#endif /* MOKALENS_H_ */

