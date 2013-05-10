/*
 * multiplane.h
 *
 *  Created on: Jan 12, 2012
 *      Author: mpetkova
 */

#ifndef MULTIPLANE_H_
#define MULTIPLANE_H_

#include <analytic_lens.h>
#include <MOKAlens.h>
#include <quadTree.h>
#include <sourceAnaGalaxy.h>
#include "tables.h"

/**
 * \brief A class to represents a lens with multiple planes.
 *
 *<pre>
 *	The rays are traced through multiple deflections.  On each plane there is a deflection
 *	solver.  An AnaLens or MOKALens can be put on one of the planes.  The other planes can be
 *	populated with random halos drawn from a mass function or they can be retrieved from an
 *	external catalog.
 *
 *
 *	Lens plane indexing scheme
 *
 *              --------------------------------  i = Nplanes-1 = source plane, No mass
 *
 *              --------------------------------  i = Nplanes-2 last plane with mass on it
 *
 *
 *              --------------------------------  i = j == (flag_input_lens % Nplanes)
 *
 *
 *              --------------------------------  i = 0 first plane with mass on it at finite distance from observer
 *
 *   Input Parameters:
 *
 *	Nplanes                     Number of lens planes
 *	flag_input_lens             Implanted lens - 0: no lens, 1: AnaLens, 2: MOKALens, The redshifts and internal parameters are the same as for constructing these lenses separately, see classes for each type of lens
 *	fov                         Field of view for the simulation region (not nessisarily the grided region)
 *	internal_profile            The internal profile type for the halos, 0 or PowerLaw,1 or NFW,2 or PseudoNFW, 3 or NSIE, 4 or NFW_NSIE .
 *  halo_to_galaxy_ratio        If NFW_NSIE is chosen this must be set to the ratio of the mass put in each component.
 *	z_source                    The "source" redshift, but actually the redshift of the last plane in the simulation.  The source can be at higher or lower redshift (see ResetSourcePlane)
 *	input_simulation_file       File that contains the input catalog of halos.  If it is missing a random set of halos will be generated.
 *	mass_func_type              The mass function used to generate random halos 0 through 2 or PS (Press & Schechter), ST (Sheth & Torman) or PowLaw (Power-law).  Not needed if input_simulation_file is provided.
 *	min_mass                    Minimum mass of halo when mass function is used (solar masses).  Not used when catalog is used.
 *	mass_scale                  The conversion between the mass units used and solar masses.  Usually 1.
 *	field_buffer                Field of view buffer in physical, rest frame Mpc.  Default is 0. Set to provide a buffer to the field of view so that halos that are centered outside the conical field of view but overlap it will be included.
 *	deflection_off              If true turns deflection off for testing purposes, default if false.
 *  background_off              If true turns deflection caused by background surface density off for testing purposes, default if false
 *
 * </pre>
 */

class MultiLens : public Lens{
public:

  void unusedHalos();

	MultiLens(InputParams& params,long *seed);
	~MultiLens();

	std::string outputfile;

	void resetNplanes(CosmoHndl cosmo, int Np);
	void resetHalos(CosmoHndl cosmo);
	
	void calc_error_test(unsigned long Npoints,Point *point,bool kappa_off);
	void calc_error_test_multi(unsigned long Npoints,Point *i_points,bool kappa_off,CosmoHndl cosmo);

	void buildHaloTrees(CosmoHndl cosmo);
	void buildHaloTrees_test(CosmoHndl cosmo);
	template <class LH> void createHaloData(CosmoHndl cosmo,long *seed);
	template <class LH> void createHaloData_test(CosmoHndl cosmo,long *seed);
	void RandomizeHost(long *seed,bool tables);
	void RandomizeSigma(long *seed,bool tables);
	double getZlens();
	void setZlens(CosmoHndl cosmo,double zlens,double zsource);
	void printMultiLens();

	void setInternalParams(CosmoHndl,SourceHndl);
	void rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off);
	void rayshooterInternal(double *ray, double *alpha, KappaType *gamma, KappaType *kappa, bool kappa_off);
	void rayshooterInternal_halos(unsigned long Npoints, Point *i_points, bool kappa_off, double *Dl_halos, double *dDl_halos);
	
	LensHndl input_lens;
	AnaLens *analens;
	MOKALens *mokalens;
	/// field of view in square degrees
	double fieldofview;

	// methods used for use with implanted sources

	short ResetSourcePlane(CosmoHndl cosmo,double z,bool nearest, unsigned long GalID=0, double *xx=NULL);

	/// Revert the source redshift to the value it was when the MultiLens was created.
	void RevertSourcePlane(){ toggle_source_plane = false;}
	//void ImplantSource(unsigned long index,CosmoHndl cosmo);
	double getSourceZ(){
		if(toggle_source_plane){
			return zs_implant;
		}else{
			return plane_redshifts[Nplanes-1];
		}
	}
	double getZmax(){return plane_redshifts[Nplanes-1];}

	/// if = 0 there is no input lens, if = 1 there is an analytic lens, if = 2 there is a MOKA lens
	InputLens flag_input_lens;
	/// Dl[j = 0...]angular diameter distances
	double *Dl;
	/// Redshifts of lens planes, 0...Nplanes.  Last one is the source redshift.
	double *plane_redshifts;
	/// dDl[j] is the distance between plane j-1 and j plane
	double *dDl;
	/// charge for the tree force solver (4*pi*G*mass_scale)
	double charge;
	/// an array of smart pointers to halo trees on each plane, uses the haloModel in the construction
	std::auto_ptr<QuadTree> *halo_tree;
	/// if >= 1, deflection in the rayshooting is switched if
	bool flag_switch_deflection_off;
	/// if >= 1, the background is switched of and only the main lens is present
	bool flag_switch_background_off;

private:

	void setCoorDist(CosmoHndl cosmo);
	
	double *coorDist_table;
	double *redshift_table;
	unsigned long NTABLE;
	bool table_set;
	void make_table(CosmoHndl cosmo);

	long *seed;

	void assignParams(InputParams& params);

	double r_print_halos;

	/* the following parameters are read in from the parameter file */

	/// type of mass function PS (0), ST (1), and power law (2) default is ST
	MassFuncType mass_func_type;
	/// slope of the mass function is mass_func_type == 2
	double pw_alpha;
	/// mass scale
	double mass_scale;
	/// min mass for the halo model
	double min_mass;
	/// internal profile type, 0=powerlaw,1=nfw,2=pseudoNfw, 3=NSIE
	IntProfType internal_profile;
	/// power law or pseudo NFW internal profile slope
	double beta;
	double galaxy_mass_fraction;


	/// read particle/halo data in from a file
	template <class DM_halo, class galaxy_halo> void readInputSimFile(CosmoHndl cosmo);

	std::string input_sim_file;
	bool sim_input_flag;
	//std::string input_gal_file;
	//bool gal_input_flag;
	bool read_sim_file;
    // TODO Margarita:  This parameter is undocumented!!!
	bool partial_cone;

	/// enables to two plane test
	bool flag_run_twop_test;
	/// enables the multi planes halos test
	bool flag_run_multip_test;

	/// pointer to first of all the halo internal structures
	std::vector<LensHalo *> halos;
	/// number of halos on all the planes
	IndexType Nhalos;
	double *halo_zs;
	double **halo_pos_Mpc;
	double **halo_pos;
	unsigned long *halo_id;

	// Variables for implanted source
	//std::auto_ptr<MultiSourceAnaGalaxy> anasource;
	/// turns source plane on and off
	bool toggle_source_plane;
	/// the distance from the source to the next plane
	double dDs_implant;
	double zs_implant,Ds_implant;
	/// This is the index of the plane at one larger distance than the new source distance
	int index_of_new_sourceplane;

	/// This is the source redshift that is read in from the parameter file and becomes the maximum redshift
	double zsource;

	double field_buffer;

	bool nfw_table_set;
	bool pnfw_table_set;

	bool second_halo;

	void quicksort(std::vector<LensHalo *> *halos,double **brr,double *arr,unsigned long *id,unsigned long N);
};

typedef  MultiLens* MultiLensHndl;


template <class LH> void MultiLens::createHaloData(
		CosmoHndl cosmo     /// cosmology
		,long *seed
	){

	if(internal_profile == NSIE || internal_profile == NFW_NSIE){
		std::cout << "ERROR: MultiLens Cann't make NSIE halos without input file yet" << std::endl;
		ERROR_MESSAGE();
		exit(1);
	}

	const int Nzbins=64;
	const int Nmassbin=64;
	int NZSamples = 50;
	std::vector<double> zbins,Nhalosbin(Nzbins);
	unsigned long i,k,j_max,k1,k2;
	std::vector<double> Logm;
	double pos_max[2], z_max;
	const double MaxLogm=16.;
	double z1, z2, mass_max,mass_tot,Nhaloestot;
	int np;
	double rr,theta,maxr;
	HALO *halo_calc = new HALO(cosmo,min_mass*mass_scale,0.0);

	double aveNhalos = cosmo->haloNumberInBufferedCone(min_mass*mass_scale,0,zsource,fieldofview*pow(pi/180,2),field_buffer,mass_func_type,pw_alpha);

	fill_linear(zbins,Nzbins,0.0,zsource);
	// construct redshift distribution table
	Nhalosbin[0] = 1;
	zbins[0] = 0;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(k)
#endif
	for(k=1;k<Nzbins-1;++k){
		Nhalosbin[k] = cosmo->haloNumberInBufferedCone(min_mass*mass_scale,zbins[k],zsource,fieldofview*pow(pi/180,2),field_buffer,mass_func_type,pw_alpha)/aveNhalos;
	}
	zbins[Nzbins-1] = zsource;
	Nhalosbin[Nzbins-1] = 0.0;

	Nhalos = (long)(poidev(float(aveNhalos), seed) );

	// allocate memory for halos
	halo_zs = new double[Nhalos];
	halo_id = new unsigned long[Nhalos];
	halo_pos = Utilities::PosTypeMatrix(Nhalos,3);

	// assign redsshifts to halos and sort them

	for(i=0;i < Nhalos;++i){
		halo_zs[i] = InterpolateYvec(Nhalosbin,zbins,ran2(seed));
	}
	std::sort(halo_zs,halo_zs + Nhalos);

	assert(halo_zs[0] < halo_zs[1]);
	assert(halo_zs[0] < halo_zs[Nhalos-1]);

	std::cout << halo_zs[Nhalos-1] << std::endl;;

	// fill the log(mass) vector
	Logm.resize(Nmassbin);
	Nhalosbin.resize(Nmassbin);
	fill_linear(Logm,Nmassbin,log10(min_mass*mass_scale),MaxLogm);

	k2 = 0;
	for(np=0,mass_max=0;np<NZSamples;np++){

		z1 = np*zsource/(NZSamples);
		z2 = (np+1)*zsource/(NZSamples);

		locateD(halo_zs-1,Nhalos,z1,&k1);
		if(k1 > Nhalos) k1 = Nhalos;
		assert(k1 == k2);
		locateD(halo_zs-1,Nhalos,z2,&k2);
		if(k2 > Nhalos) k2 = Nhalos;

		Nhaloestot = cosmo->haloNumberInBufferedCone(pow(10,Logm[0]),z1,z2,fieldofview*pow(pi/180,2),field_buffer,mass_func_type,pw_alpha);

		Nhalosbin[0] = 1;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(k)
#endif
		for(k=1;k<Nmassbin-1;k++){
			// cumulative number density in one square degree
			Nhalosbin[k] = cosmo->haloNumberInBufferedCone(pow(10,Logm[k]),z1,z2,fieldofview*pow(pi/180,2),field_buffer,mass_func_type,pw_alpha)
					/Nhaloestot;
		}
		Nhalosbin[Nmassbin-1] = 0;

		for(i = k1; i < k2; i++){
// TODO Ben - this will never work for NSIE or NFW+NSIE models fix it
			/// positions need to be in radians initially

			double Ds = cosmo->angDist(0,halo_zs[i]);

			maxr = pi*sqrt(fieldofview/pi)/180. + field_buffer/Ds; // fov is a circle
			rr = maxr*sqrt(ran2(seed));

			assert(rr == rr);

			theta = 2*pi*ran2(seed);


			halo_pos[i][0] = rr*cos(theta)*Ds;
			halo_pos[i][1] = rr*sin(theta)*Ds;
			halo_pos[i][2] = 0;

			halos.push_back(new LH);

			halos[i]->set_slope(beta);

			float mass = pow(10,InterpolateYvec(Nhalosbin,Logm,ran2 (seed)));

			halos[i]->set_mass(mass/mass_scale);

			halo_calc->reset(mass,halo_zs[i]);

			halos[i]->set_Rmax(halo_calc->getRvir());
			halos[i]->set_rscale(halos[i]->get_Rmax()/halo_calc->getConcentration(0));

			if(halos[i]->get_mass() > mass_max) {
				mass_max = halos[i]->get_mass();
				j_max = i;
				pos_max[0] = halo_pos[i][0];
				pos_max[1] = halo_pos[i][1];
				z_max = halo_zs[i];
			}

			halo_id[i] = i;

			mass_tot += halos[i]->get_mass();
		}

		Nhalosbin.empty();
	}

	assert(k2 == Nhalos);
	delete halo_calc;

	std::cout << Nhalos << " halos created." << std::endl
	    << "Max input mass = " << mass_max << "  R max = " << halos[j_max]->get_Rmax()
	    << " at z = " << z_max << std::endl;

	std::cout << "leaving MultiLens::createHaloData_buffered()" << std::endl;
}

/*
 * A test function that creates two halos, at random positions on the sky, but in
 * such a way that they will be onto two separate lensing planes.
 * This is used to test the convergence of the ray shooter.
 *
 */
template <class LH> void MultiLens::createHaloData_test(
		CosmoHndl cosmo     /// cosmology
		,long *seed
	){

	HALO *halo_calc = new HALO(cosmo,min_mass*mass_scale,0.0);

	Nhalos = 2;

	// allocate memory for halos
	halo_zs = new double[Nhalos];
	halo_id = new unsigned long[Nhalos];
	halo_pos = Utilities::PosTypeMatrix(Nhalos,3);

	for(int i=0;i<Nhalos;i++){
		double maxr = pi*sqrt(fieldofview/pi)/180.; // fov is a circle
		double rr = 0.5*maxr;

		assert(rr == rr);

		double theta = 2*pi*ran2(seed);

		halo_pos[i][0] = rr*cos(theta);
		halo_pos[i][1] = rr*sin(theta);
		halo_pos[i][2] = 0;

		halo_zs[i] = plane_redshifts[i]+0.1;

		halos.push_back(new LH);

		halos[i]->set_slope(beta);
		halos[i]->set_mass(ran2(seed)*1e12/mass_scale);
		halo_calc->reset(halos[i]->get_mass()*mass_scale,halo_zs[i]);

		halos[i]->set_Rmax(halo_calc->getRvir());
		halos[i]->set_rscale(halos[i]->get_Rmax()/halo_calc->getConcentration(0));

		std::cout<< "halo z " << halo_zs[i] << std::endl;
		std::cout<< "halo pos " << halo_pos[i][0] << " " << halo_pos[i][1] << std::endl;
		std::cout<< "halo mass " << halos[i]->get_mass()*mass_scale << std::endl;
		std::cout<< "halo Rmax and rscale " << halos[i]->get_Rmax() << " " << halos[i]->get_rscale() << std::endl;
	}

	delete halo_calc;

	std::cout << "leaving MultiLens::createHaloData_test()" << std::endl;
}

/**
 * \brief Read in information from a Virgo Millennium Data Base http://gavo.mpa-garching.mpg.de/MyMillennium/
 *
 * query select * from MockLensing.dbo.m21_20_39_021_bc03_Ben_halos
 *
 * This is information on the dark matter halos only.  There are 13 entries in each line separated by commas.
 * The comments must be removed from the beginning of the data file and the total number of halos must be added
 * as the first line.
 */
template <class DM_halo, class galaxy_halo> void MultiLens::readInputSimFile(CosmoHndl cosmo){

	double ra,dec,z,vmax,vdisp,r_halfmass;
	unsigned long i,j;
	unsigned long haloid,idd,np;
	double mo=7.3113e10,M1=2.8575e10,gam1=7.17,gam2=0.201,be=0.557;

	double rmax=0,rtmp=0;

	//int index;

	if(internal_profile == PseudoNFW){
		ERROR_MESSAGE();
		std::cout << "Input to MultiLens from a simulation is not yet enabled for PseudoNFW profiles"
				<< std::endl << "Change this is parameter file" << std::endl;
		exit(1);
	}
	if(internal_profile == PowerLaw){
		ERROR_MESSAGE();
		std::cout << "Input to MultiLens from a simulation is not yet enabled for PowerLaw profiles"
				<< std::endl << "Change this is parameter file" << std::endl;
		exit(1);
	}

	std::ifstream file_in(input_sim_file.c_str());
	if(!file_in){
		std::cout << "Can't open file " << input_sim_file << std::endl;
		exit(1);
	}

	// skip through header information in data file
	i=0;
	while(file_in.peek() == '#'){
		file_in.ignore(10000,'\n');
		++i;
	}
	std::cout << "skipped "<< i << " comment lines in file " << input_sim_file << std::endl;

	std::vector<double> halo_zs_vec;
	std::vector<double *> halo_pos_vec;
	std::vector<unsigned long> halo_id_vec;

	// read in data
	int j_max;
	double mass_max=0,R_max=0,V_max=0,minmass=1e30;
	double *theta;
	int ncolumns = 9;

	void *addr[ncolumns];
	addr[0] = &haloid;
	addr[1] = &idd;
	addr[2] = &ra;
	addr[3] = &dec;
	addr[4] = &z;
	addr[5] = &np;
	addr[6] = &vdisp;
	addr[7] = &vmax;
	addr[8] = &r_halfmass;

	unsigned long myint;
	double mydouble;
	std::string myline;
	std::string strg;
	std::string f=",";
	std::stringstream buffer;

	for(i=0,j=0 ; ; ++i){
		// read a line of data
		myline.clear();
		std::getline(file_in,myline);

		if(myline[0] == '#')
			break;
		for(int l=0;l<ncolumns; l++){
			int pos = myline.find(f);
			strg.assign(myline,0,pos);
			buffer << strg;
			if(l == 0 || l == 1 || l == 5){
				buffer >> myint;
				*((unsigned long *)addr[l]) = myint;
			}
			else{
				buffer >> mydouble;
				*((double *)addr[l]) = mydouble;
			}
			myline.erase(0,pos+1);
			strg.clear();
			buffer.clear();
			buffer.str(std::string());
		}

		/// pos in radians
		theta = new double[2];
		/// pos in physical radians
		double Ds = cosmo->angDist(0,z);
		theta[0] = ra*pi/180.*Ds;
		theta[1] = dec*pi/180.*Ds;

        if(rmax < (rtmp = theta[0]*theta[0]+theta[1]*theta[1])) rmax = rtmp;

		if(partial_cone){
			double r = sqrt(theta[0]*theta[0]+theta[1]*theta[1]);
			if(r > 1.5*mokalens->map->boxlMpc)
				continue;
		}

		if(np > 0.0 && vdisp > 0.0 && z <= zsource){

			halo_pos_vec.push_back(theta);

			halos.push_back(new DM_halo);

			float mass = np*8.6e8/cosmo->gethubble();

			halos[j]->set_mass(mass);
			halos[j]->set_Rmax(mass*Grav/2/pow(vmax/lightspeed,2));  // SIS value

			if(halos[j]->get_mass() > mass_max) {
				mass_max = halos[j]->get_mass();
				j_max = j;
			}
			if(halos[j]->get_mass() < minmass) {
				minmass = halos[j]->get_mass();
			}

			/*
			 * set the other properties of the LensHalo, such as rscale or the NSIE properties
			 */

			halos[j]->set_slope(beta);
			halos[j]->set_internal(seed,vmax,r_halfmass*cosmo->gethubble());

			if(halos[j]->get_Rmax() > R_max) R_max = halos[j]->get_Rmax();
			if(vdisp > V_max) V_max = vdisp;


			halo_zs_vec.push_back(z);
			halo_id_vec.push_back(haloid);

			halos[j]->set_mass(halos[j]->get_mass()/mass_scale);

			++j;

			if(second_halo){
				galaxy_mass_fraction = 2*mo*pow(mass/M1,gam1)
				  /pow(1+pow(mass/M1,be),(gam1-gam2)/be)/mass;
				if(galaxy_mass_fraction > 1.0) galaxy_mass_fraction = 1;

				halo_pos_vec.push_back(theta);

				halos.push_back(new galaxy_halo);

				halos[j]->set_mass(mass*galaxy_mass_fraction);
				halos[j-1]->set_mass(mass*(1-galaxy_mass_fraction)/mass_scale);

				halos[j]->set_Rmax(mass*galaxy_mass_fraction*Grav/2/pow(vmax/lightspeed,2));  // SIS value


				if(halos[j]->get_mass() > mass_max) {
					mass_max = halos[j]->get_mass();
					j_max = j;
				}
				if(halos[j]->get_mass() < minmass) {
					minmass = halos[j]->get_mass();
				}

				/*
				 * set the other properties of the LensHalo, such as rscale or the NSIE properties
				 */

				halos[j]->set_slope(beta);
				halos[j]->set_internal(seed,vmax,r_halfmass*cosmo->gethubble());

				if(halos[j]->get_Rmax() > R_max) R_max = halos[j]->get_Rmax();
				if(vdisp > V_max) V_max = vdisp;


				halo_zs_vec.push_back(z);
				halo_id_vec.push_back(haloid);

				halos[j]->set_mass(halos[j]->get_mass()/mass_scale);

				++j;
			}

		}
	}
	file_in.close();
	std::cout << halos.size() << " halos read in."<< std::endl
			<< "Max input mass = " << mass_max << "  R max = " << R_max << "  V max = " << V_max << std::endl;

	Nhalos = halos.size();

	/// setting the minimum halo mass in the simulation
	min_mass = minmass;
	if(field_buffer > 0.0){
		std::cout << "Overiding field_buffer to make it 0 because halos are read in." << std::endl;
		field_buffer = 0.0;
	}

	if(partial_cone == false)
		fieldofview = pi*rmax*pow(180/pi,2);  // Resets field of view to estimate of inputed one

	halo_zs = new double[Nhalos];
	halo_id = new unsigned long[Nhalos];
	halo_pos = Utilities::PosTypeMatrix(Nhalos,3);

	for(i=0;i<Nhalos;++i){
		halo_id[i] = halo_id_vec[i];
		halo_zs[i] = halo_zs_vec[i];
		halo_pos[i] = halo_pos_vec[i];
	}

	std::cout << "sorting in MultiLens::readInputSimFile()" << std::endl;
	// sort the halos by readshift
	MultiLens::quicksort(&halos,halo_pos,halo_zs,halo_id,Nhalos);

	std::cout << "leaving MultiLens::readInputSimFile()" << std::endl;

	read_sim_file = true;
}


#endif /* MULTIPLANE_H_ */
