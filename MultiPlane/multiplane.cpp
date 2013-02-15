/*
 * halo_model.cpp
 *
 *  Created on: Jan 14, 2012
 *      Author: mpetkova, bmetcalf
 */

#include "slsimlib.h"

using namespace std;

HaloData::HaloData(HaloStructure *halostrucs,double sigma_back,PosType **positions, double *zz, unsigned long *id,unsigned long Nhaloss,double Dl):
	pos(positions), halos(halostrucs), Nhalos(Nhaloss),z(zz),haloID(id),sigma_background(sigma_back)
{
	//convert to physical Mpc on the plane
	int i;
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
	for(i=0; i<Nhalos; i++){
		pos[i][0]*=Dl;
		pos[i][1]*=Dl;
	}
}

HaloData::~HaloData(){
}

void MultiLens::make_table(CosmoHndl cosmo){
	int i;
	double x, dx = (zsource+1.0)/(double)NTABLE;

	coorDist_table = new double[NTABLE];
	redshift_table = new double[NTABLE];
	
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,x)
#endif
	for(i = 0 ; i< NTABLE; i++){
		x = i*dx;
		redshift_table[i] = x;
		coorDist_table[i] = cosmo->coorDist(0,x);
	}
	table_set=true;
}

void MultiLens::resetNplanes(CosmoHndl cosmo, int Np){
	int i,j;
	///revert back to radians
	for(j=0; j<Nplanes-1; j++){
		if(flag_input_lens && j == (flag_input_lens % Nplanes))
			continue;
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
		for(i=0; i<halo_data[j]->Nhalos; i++){
			halo_data[j]->pos[i][0]/=Dl[j]/(1+plane_redshifts[j]);
			halo_data[j]->pos[i][1]/=Dl[j]/(1+plane_redshifts[j]);
		}
#ifdef _OPENMP
#pragma omp barrier
#endif
	}

  delete[] halo_tree;
  delete[] halo_data;
  
  delete[] Dl;
  delete[] plane_redshifts;
  delete[] dDl;
  
  Nplanes = Np;

  plane_redshifts = new double[Nplanes];
  Dl = new double[Nplanes];
  dDl = new double[Nplanes];
  
  halo_tree = new auto_ptr<QuadTree>[Nplanes-1];
  halo_data = new auto_ptr<HaloData>[Nplanes-1];
   
  setCoorDist(cosmo);
  buildHaloTrees(cosmo);
}

void MultiLens::resetHalos(CosmoHndl cosmo){
  delete[] halo_tree;
  delete[] halo_data;
  delete[] halos;

  delete[] halo_zs;
  delete[] halo_id;
  free_PosTypeMatrix(halo_pos,Nhalos,3);
  if(flag_run_multip_test)
	  free_PosTypeMatrix(halo_pos_Mpc,Nhalos,3);

  halo_tree = new auto_ptr<QuadTree>[Nplanes-1];
  halo_data = new auto_ptr<HaloData>[Nplanes-1];

  if(field_buffer > 0.0) createHaloData_buffered(cosmo,seed);
  else{
	  if(flag_run_twop_test) createHaloData_test(cosmo,seed);
	  else createHaloData(cosmo,seed);
  }

  if(flag_run_twop_test) buildHaloTrees_test(cosmo);
  else buildHaloTrees(cosmo);
}

/*
 * \ingroup Constructor
 * allocates space for the halo trees and the inout lens, if there is any
 */
MultiLens::MultiLens(InputParams& params,long *my_seed) : Lens(){
	assignParams(params);

	NTABLE=1000;
	table_set = false;

	// flag to determine if halos are created randomly or read in from a external simulation.
	if(input_sim_file.size() < 1) sim_input_flag = false;
	else sim_input_flag = true;

	std::cout << input_sim_file.c_str() << std::endl;

	read_sim_file = false;

	plane_redshifts = new double[Nplanes];
	Dl = new double[Nplanes];
	dDl = new double[Nplanes];

	charge = 4*pi*Grav*mass_scale;
	std::cout << "charge: " << charge << std::endl;

	halo_tree = new auto_ptr<QuadTree>[Nplanes-1];
	halo_data = new auto_ptr<HaloData>[Nplanes-1];

	switch(flag_input_lens){
	case null:
		input_lens = NULL;
		break;
	case ana_lens:
		input_lens = new AnaLens(params);
		analens = static_cast<AnaLens*>(input_lens);
		break;
	case moka_lens:
		input_lens = new MOKALens(params);
		mokalens = static_cast<MOKALens*>(input_lens);
		fieldofview = pow(1.5*mokalens->map->boxlrad*180/pi,2.0);
		break;
	default:
		ERROR_MESSAGE();
		cout << "Incorrect flag_input_lens selected! Please choose from:" << endl;
		cout << "0: no lens, 1: AnaLens, 2: MOKALens" << endl;
		exit(1);
		break;
	}

	// initially let source be the one inputed from parameter file
	index_of_new_sourceplane = -1;
	toggle_source_plane = false;

	seed = my_seed;
}

MultiLens::~MultiLens(){
	delete[] halo_tree;
	delete[] halo_data;

	delete[] Dl;
	delete[] plane_redshifts;
	delete[] dDl;

	delete[] halos;
	delete[] halo_zs;
	delete[] halo_id;
	free_PosTypeMatrix(halo_pos,Nhalos,3);
	if(flag_run_multip_test)
		free_PosTypeMatrix(halo_pos_Mpc,Nhalos,3);
	
	if(flag_input_lens)
		delete input_lens;

	delete[] coorDist_table;
	delete[] redshift_table;
}

/// Retrieve input parameters for construction
void MultiLens::assignParams(InputParams& params){

  if(!params.get("outputfile",outputfile)){
		  ERROR_MESSAGE();
		  cout << "parameter outputfile needs to be set in the parameter file " << params.filename() << endl;
		  exit(0);
	}	
	if(!params.get("Nplanes",Nplanes)){
		  ERROR_MESSAGE();
		  cout << "parameter Nplanes needs to be set in the parameter file " << params.filename() << endl;
		  exit(0);
	}
	if(!params.get("flag_input_lens",flag_input_lens)){
		  ERROR_MESSAGE();
		  cout << "parameter flag_input_lens needs to be set in the parameter file " << params.filename() << endl;
		  exit(0);
	}
	if(!params.get("internal_profile",internal_profile)){
		  ERROR_MESSAGE();
		  cout << "parameter internal_profile needs to be set in the parameter file " << params.filename() << endl;
		  exit(0);
	}
	if(internal_profile == NFW_NSIE){
		if(!params.get("galaxy_mass_fraction",galaxy_mass_fraction)){
			  ERROR_MESSAGE();
			  cout << "to construct a NFW + NSIE model the parameter halo_to_galaxy_ratio needs to be set in the parameter file " << params.filename() << endl;
			  exit(0);
		}
	}
	if(!params.get("z_source",zsource)){
		  ERROR_MESSAGE();
		  cout << "parameter z_source needs to be set in the parameter file " << params.filename() << endl;
		  exit(0);
	}

	if(!params.get("input_simulation_file",input_sim_file)){

		// No simulation input file provided
		if(!params.get("mass_func_type",mass_func_type)){
			  ERROR_MESSAGE();
			  cout << "parameter mass_func_type needs to be set in the parameter file " << params.filename() << endl;
			  exit(0);
		}
		if(!params.get("min_mass",min_mass)){
			  ERROR_MESSAGE();
			  cout << "parameter min_mass needs to be set in the parameter file " << params.filename() << endl;
			  exit(0);
		}
		if(!params.get("fov",fieldofview)){
			  ERROR_MESSAGE();
			  cout << "parameter fov needs to be set in the parameter file " << params.filename() << endl;
			  exit(0);
		}
		if(!params.get("field_buffer",field_buffer)){
			field_buffer = 0.0;
			cout << "default field buffer of 0 Mpc is being used." << endl;
		}
		partial_cone=false;
	}else{
		if(!params.get("partial_cone",partial_cone))	partial_cone=false;
		min_mass = 0.0;
	}

	if(!params.get("mass_scale",mass_scale)){
		  ERROR_MESSAGE();
		  cout << "parameter mass_scale needs to be set in the parameter file " << params.filename() << endl;
		  exit(0);
	}
	// parameters with default values
	if(!params.get("alpha",pw_alpha))                pw_alpha = 1./6.;
	if(!params.get("internal_slope_pw",pw_beta))     pw_beta = -1.0;
	if(!params.get("internal_slope_pnfw",pnfw_beta)) pnfw_beta = 2.0;
	if(!params.get("deflection_off",flag_switch_deflection_off)) flag_switch_deflection_off = false;
	if(!params.get("background_off",flag_switch_background_off)) flag_switch_background_off = false;
	if(!params.get("twop_test",flag_run_twop_test)) flag_run_twop_test = false;
	if(!params.get("multip_test",flag_run_multip_test)) flag_run_multip_test = false;
	if(!params.get("print_halos",r_print_halos)) r_print_halos = 0.0;

	// Some checks for valid parameters
	  if(pw_beta >= 0){
		  ERROR_MESSAGE();
		  cout << "Internal slope >=0 not possible." << endl;
		  exit(1);
	  }

	  if(pnfw_beta <= 0){
		  ERROR_MESSAGE();
		  cout << "Internal slope <=0 not possible." << endl;
		  exit(1);
	  }

	  if(pnfw_beta / floor(pnfw_beta) > 1.0){
		  ERROR_MESSAGE();
		  cout << "Internal slope needs to be a whole number." << endl;
		  exit(1);
	  }

	  if(input_sim_file.size() < 1 && internal_profile == NSIE){
		  ERROR_MESSAGE();
		  cout << "The NSIE internal profile works only for Millenium DM simulations for now." << endl;
		  cout << "Set input_simulation_file in sample_paramfile." << endl;
		  exit(1);
	  }

	  // to compensate for the last plane, which is the source plane
	  Nplanes++;

	  // to compensate for additional lens planes
	  if(flag_input_lens)
		  Nplanes++;

	  // convert to square degrees
	  fieldofview /= 3600. * 3600.;

	  if(partial_cone == true && flag_input_lens != 2){
		  ERROR_MESSAGE();
		  cout << "Partial cone selection is possible only with a MOKA lens!" << endl;
		  exit(1);
	  }

	  printMultiLens();
}

void MultiLens::printMultiLens(){

	cout << endl << "**multi lens model**" << endl;

	cout << "Nplanes " << Nplanes << endl;

	cout << "mass scale " << mass_scale << endl;

	cout << "min mass " << min_mass << endl;

	cout << "flag input lens " << flag_input_lens << endl;
	switch(flag_input_lens){
	case null:
		cout << "  No input lens specified " << endl;
		break;
	case ana_lens:
		cout << "  AnaLens " << endl;
		break;
	case moka_lens:
		cout << "  MOKALens " << endl;
		break;
	}

	cout << "field of view " << fieldofview << endl;

	cout << "internal profile type " << internal_profile << endl;
	switch(internal_profile){
	case PowerLaw:
		cout << "  Power law internal profile " << endl;
		cout << "  slope: " << pw_beta << endl;
		break;
	case NFW:
		cout << "  NFW internal profile " << endl;
		break;
	case PseudoNFW:
		cout << "  Pseudo NFW internal profile " << endl;
		cout << "  slope: " << pnfw_beta << endl;
		break;
	case NSIE:
		cout << "  NonSingular Isothermal Ellipsoid internal profile " << endl;
		break;
	case NFW_NSIE:
		cout << "  NFW plus NonSingular Isothermal Ellipsoid internal profile " << endl;
		break;
	}

	cout << "mass function type " << mass_func_type << endl;
	switch(mass_func_type){
		case PS:
			cout << "  Press-Schechter mass function " << endl;
			break;
		case ST:
			cout << "  Sheth-Tormen mass function " << endl;
			break;
		case PL:
			cout << "  Power law mass function " << endl;
			cout << "  slope: " << pw_alpha << endl;
			break;
		}

	cout << endl;
}

/**
 *  \brief In this constructor the halos are created from a mass function
 *
 *  The lightcone is populated by dividing it into small redshifts bins and creating halos,
 *  which then are sorted according to redshift
 */
void MultiLens::createHaloData(
		CosmoHndl cosmo     /// cosmology
		,long *seed
	){

	if(internal_profile == NSIE || internal_profile == NFW_NSIE){
		std::cout << "ERROR: MultiLens Cann't make NSIE halos without input file yet" << std::endl;
		ERROR_MESSAGE();
		exit(1);
	}
  HALO *halo_calc = new HALO(cosmo,min_mass*mass_scale,0.0);
  
  std::vector<double> Logm,Nhalosbin;
  std::vector<HaloStructure> halo_vec;
  std::vector<double> halo_zs_vec;
  std::vector<double *> halo_pos_vec;
  std::vector<unsigned long> halo_id_vec;

  double *pos, pos_max[2], z_max;
  const int Nmassbin=32;
  const double MaxLogm=17.;

  
  Logm.resize(Nmassbin);
  Nhalosbin.resize(Nmassbin);
  
  /* fill the log(mass) vector */
  
  fill_linear(Logm,Nmassbin,log10(min_mass*mass_scale),MaxLogm);

  int Nsample = 50;
  
  double tmp_dDl, Dl1, Dl2, z1, z2, mass_max,mass_tot;
  tmp_dDl = cosmo->coorDist(0,zsource)/(Nsample);
  int np;
  unsigned long h_index=0,j_max,k;
  for(np=0,mass_max=0;np<Nsample;np++){

	mass_tot = 0.0;
    double Nhaloestot;
    if(np == 0){
      z1 = 0.0;
      Dl1 = 0.0;
    }
    else{
      Dl1 = np*tmp_dDl;
      locateD(coorDist_table-1,NTABLE,Dl1,&k);
      z1 = redshift_table[k];
    }

    Dl2 = Dl1+tmp_dDl;
    
    locateD(coorDist_table-1,NTABLE,Dl2,&k);
    z2 = redshift_table[k];

    double tailarea = cosmo->haloNumberDensityOnSky(pow(10,MaxLogm),z1,z2,mass_func_type,pw_alpha)*fieldofview;
    Nhalosbin[0] = cosmo->haloNumberDensityOnSky(pow(10,Logm[0]),z1,z2,mass_func_type,pw_alpha)*fieldofview;
    
    /*
    std::cout << "tail = " << tailarea << "  " << tailarea/Nhalosbin[0] << " % "<< std::endl;
    std::cout << "number tail above 1.0e16 = " << 100*cosmo->haloNumberDensityOnSky(1.0e16,z1,z2,mass_func_type,pw_alpha)*fieldofview/Nhalosbin[0]
     		<< " % "<< std::endl;
    std::cout << "number tail above 1.0e15 = " << 100*cosmo->haloNumberDensityOnSky(1.0e15,z1,z2,mass_func_type,pw_alpha)*fieldofview/Nhalosbin[0]
     		<< " % "<< std::endl;
    std::cout << "number tail above 1.0e14 = " << 100*cosmo->haloNumberDensityOnSky(1.0e14,z1,z2,mass_func_type,pw_alpha)*fieldofview/Nhalosbin[0]
     		<< " % "<< std::endl;
     */

    Nhaloestot = Nhalosbin[0]-tailarea;
    Nhalosbin[0] = 1;
    int k;
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(k)
#endif
    for(k=1;k<Nmassbin-1;k++){
      // cumulative number density in one square degree
      Nhalosbin[k] = cosmo->haloNumberDensityOnSky(pow(10,Logm[k]),z1,z2,mass_func_type,pw_alpha)*fieldofview;
      // normalize the cumulative distribution to one
      Nhalosbin[k] = (Nhalosbin[k]-tailarea)/Nhaloestot;
    }
    Nhalosbin[Nmassbin-1] = 0;
    
    long Nh = (long)(poidev(float(Nhaloestot), seed) + 0.5);
    
    double rr,theta,maxr,zi;
    unsigned long i;
    for(i = 0,mass_max=0; i < Nh; i++){
      HaloStructure halo;
      
      zi = z1+(z2-z1)*ran2 (seed); 
       /// positions need to be in radians initially
      maxr = pi*sqrt(fieldofview/pi)/180.; // fov is a circle
      rr = maxr*sqrt(ran2(seed));
      
      theta = 2*pi*ran2(seed);
      
      pos = new double[3];
      pos[0] = rr*cos(theta);
      pos[1] = rr*sin(theta);
      pos[2] = 0.0;
      
      // TODO Ben - this will never work for NSIE or NFW+NSIE models fix it
     halo.mass = pow(10,InterpolateYvec(Nhalosbin,Logm,ran2 (seed)));
      halo_calc->reset(halo.mass,zi);
      halo.mass /= mass_scale;
      halo.Rmax = halo_calc->getRvir();
      halo.rscale = halo.Rmax/halo_calc->getConcentration(0);

      if(halo.mass > mass_max) {
    	  mass_max = halo.mass;
    	  j_max = h_index;
    	  pos_max[0] = pos[0];
    	  pos_max[1] = pos[1];
    	  z_max = zi;
      }
      
      halo_vec.push_back(halo);
      halo_zs_vec.push_back(zi);
      halo_pos_vec.push_back(pos);
      halo_id_vec.push_back(h_index);
      h_index++;
      
      mass_tot += halo.mass;
    }
    
    Nhalosbin.empty();


    double mftot = cosmo->totalMassDensityinHalos(mass_func_type,pw_alpha,pow(10,Logm[0]),(z1+z2)/2,z1,z2)
       	       		*pow(cosmo->angDist(0,(z1+z2)/2),2)*fieldofview*pow(pi/180,2);
    // Test lines
    /*
     std::cout <<
       mass_tot << "  "<< (mftot - mass_tot) << std::endl
        << " in tail above 1.0e16 = " << cosmo->totalMassDensityinHalos(mass_func_type,pw_alpha,1.0e16,(z1+z2)/2,z1,z2)
  		*pow(cosmo->angDist(0,(z1+z2)/2),2)*fieldofview*pow(pi/180,2)/mftot
        << " in tail above 1.0e15 = " << cosmo->totalMassDensityinHalos(mass_func_type,pw_alpha,1.0e15,(z1+z2)/2,z1,z2)
   		*pow(cosmo->angDist(0,(z1+z2)/2),2)*fieldofview*pow(pi/180,2)/mftot
   		<< " in tail above 1.0e14 = " << cosmo->totalMassDensityinHalos(mass_func_type,pw_alpha,1.0e14,(z1+z2)/2,z1,z2)
		*pow(cosmo->angDist(0,(z1+z2)/2),2)*fieldofview*pow(pi/180,2)/mftot
		<< " Msun  number of halo = " << Nh << std::endl;
     */
  }
  
  delete halo_calc;
  
  Nhalos = halo_vec.size();
  
  std::cout << Nhalos << " halos created."<< std::endl
	    << "Max input mass = " << mass_max*mass_scale << "  R max = " << halo_vec[j_max].Rmax
	    << " at z = " << z_max << std::endl;
  
  halos = new HaloStructure[Nhalos];
  halo_zs = new double[Nhalos];
  halo_id = new unsigned long[Nhalos];
  halo_pos = PosTypeMatrix(Nhalos,3);
  if(flag_run_multip_test)
	  halo_pos_Mpc = PosTypeMatrix(Nhalos,3);
  
  for(int i=0;i<Nhalos;++i){
    halo_id[i] = halo_id_vec[i];
    halo_zs[i] = halo_zs_vec[i];
    halo_pos[i] = halo_pos_vec[i];
    halos[i] = halo_vec[i];
  }

  std::cout << "sorting in MultiLens::createHaloData()" << std::endl;
  // sort the halos by readshift
  MultiLens::quicksort(halos,halo_pos,halo_zs,halo_id,Nhalos);
  
  if(flag_run_multip_test){
	  for(int i=0;i<Nhalos;++i){
		  double Dh = cosmo->angDist(0.0,halo_zs[i]);
		  halo_pos_Mpc[i][0] = halo_pos[i][0]*Dh;
		  halo_pos_Mpc[i][1] = halo_pos[i][1]*Dh;
	  }
  }

  std::cout << "leaving MultiLens::createHaloData()" << std::endl;
}

void MultiLens::createHaloData_buffered(
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
	halos = new HaloStructure[Nhalos];
	halo_id = new unsigned long[Nhalos];
	halo_pos = PosTypeMatrix(Nhalos,3);

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
			maxr = pi*sqrt(fieldofview/pi)/180. + field_buffer/cosmo->angDist(0,halo_zs[i]); // fov is a circle
			rr = maxr*sqrt(ran2(seed));

			assert(rr == rr);

			theta = 2*pi*ran2(seed);

			halo_pos[i][0] = rr*cos(theta);
			halo_pos[i][1] = rr*sin(theta);
			halo_pos[i][2] = 0;

			halos[i].mass = pow(10,InterpolateYvec(Nhalosbin,Logm,ran2 (seed)));
			halo_calc->reset(halos[i].mass,halo_zs[i]);
			halos[i].mass /= mass_scale;
			halos[i].Rmax = halo_calc->getRvir();
			halos[i].rscale = halos[i].Rmax/halo_calc->getConcentration(0);

			if(halos[i].mass > mass_max) {
				mass_max = halos[i].mass;
				j_max = i;
				pos_max[0] = halo_pos[i][0];
				pos_max[1] = halo_pos[i][1];
				z_max = halo_zs[i];
			}

			halo_id[i] = i;

			mass_tot += halos[i].mass;
		}

		Nhalosbin.empty();

/*
		double mftot = cosmo->totalMassDensityinHalos(mass_func_type,pw_alpha,pow(10,Logm[0]),(z1+z2)/2,z1,z2)
       	       		*pow(cosmo->angDist(0,(z1+z2)/2),2)*fieldofview*pow(pi/180,2);
		// Test lines
		std::cout <<
				mass_tot << "  "<< (mftot - mass_tot) << std::endl
				<< " in tail above 1.0e16 = " << cosmo->totalMassDensityinHalos(mass_func_type,pw_alpha,1.0e16,(z1+z2)/2,z1,z2)
				*pow(cosmo->angDist(0,(z1+z2)/2),2)*fieldofview*pow(pi/180,2)/mftot
				<< " in tail above 1.0e15 = " << cosmo->totalMassDensityinHalos(mass_func_type,pw_alpha,1.0e15,(z1+z2)/2,z1,z2)
				*pow(cosmo->angDist(0,(z1+z2)/2),2)*fieldofview*pow(pi/180,2)/mftot
				<< " in tail above 1.0e14 = " << cosmo->totalMassDensityinHalos(mass_func_type,pw_alpha,1.0e14,(z1+z2)/2,z1,z2)
				*pow(cosmo->angDist(0,(z1+z2)/2),2)*fieldofview*pow(pi/180,2)/mftot
				<< " Msun  number of halo in bin = " << k2-k1 << std::endl;
				*/
	}

	assert(k2 == Nhalos);
	delete halo_calc;

	std::cout << Nhalos << " halos created." << std::endl
	    << "Max input mass = " << mass_max << "  R max = " << halos[j_max].Rmax
	    << " at z = " << z_max << std::endl;

	std::cout << "leaving MultiLens::createHaloData_buffered()" << std::endl;
}

/*
 * A test function that creates two halos, at random positions on the sky, but in
 * such a way that they will be onto two separate lensing planes.
 * This is used to test the convergence of the ray shooter.
 *
 */
void MultiLens::createHaloData_test(
		CosmoHndl cosmo     /// cosmology
		,long *seed
	){

	HALO *halo_calc = new HALO(cosmo,min_mass*mass_scale,0.0);

	Nhalos = 2;

	// allocate memory for halos
	halo_zs = new double[Nhalos];
	halos = new HaloStructure[Nhalos];
	halo_id = new unsigned long[Nhalos];
	halo_pos = PosTypeMatrix(Nhalos,3);

	for(int i=0;i<Nhalos;i++){
		double maxr = pi*sqrt(fieldofview/pi)/180.; // fov is a circle
		double rr = 0.5*maxr;

		assert(rr == rr);

		double theta = 2*pi*ran2(seed);

		halo_pos[i][0] = rr*cos(theta);
		halo_pos[i][1] = rr*sin(theta);
		halo_pos[i][2] = 0;

		halo_zs[i] = plane_redshifts[i]+0.1;
		halos[i].mass = ran2(seed)*1e12;
		halo_calc->reset(halos[i].mass,halo_zs[i]);
		halos[i].mass /= mass_scale;
		halos[i].Rmax = halo_calc->getRvir();
		halos[i].rscale = halos[i].Rmax/halo_calc->getConcentration(0);

		std::cout<< "halo z " << halo_zs[i] << std::endl;
		std::cout<< "halo pos " << halo_pos[i][0] << " " << halo_pos[i][1] << std::endl;
		std::cout<< "halo mass " << halos[i].mass*mass_scale << std::endl;
		std::cout<< "halo Rmax and rscale " << halos[i].Rmax << " " << halos[i].rscale << std::endl;
	}

	delete halo_calc;

	std::cout << "leaving MultiLens::createHaloData_test()" << std::endl;
}

/**
 * Populates the planes with halos by dividing the space around the planes into
 * equal redshift distances, where the plane with the input lens is excluded
 * since it will not contain any halos
 *
 * Then the halo trees are built, depending on the internal profile model that
 * has been chosen in the parameter file
 */
void MultiLens::buildHaloTrees(
		CosmoHndl cosmo /// the cosmology
		){
	int j, Ntot;
	unsigned long ind;
	double z1, z2;
	unsigned long j1,j2;

	std::cout << "MultiLens::buildHaloTrees zsource = " << zsource << std::endl;


	ofstream file_area("halos.dat");
	if(r_print_halos){
		if(!file_area){
			cout << "unable to create file " << "halos.dat" << endl;
			exit(1);
		}
		file_area << "redshift      mass        pos[0]       pos[1]" <<endl;
	}

	assert(plane_redshifts[Nplanes-1] == zsource);

	for(j=0,Ntot=0;j<Nplanes-1;j++){
		if(flag_input_lens && j == (flag_input_lens % Nplanes))
			continue;

		/*
		 * Setting the redshift range
		 * If there is a plane with an input lens on it, it is skipped over
		 * since it will not contain any halos
		 */
		if(j == 0) z1 = 0.0;
		else{
			locateD(coorDist_table-1,NTABLE,(Dl[j]-0.5*dDl[j]),&ind);
			z1 = redshift_table[ind];
		}

		if(flag_input_lens && j-1 == (flag_input_lens % Nplanes)){
			locateD(coorDist_table-1,NTABLE,(Dl[j] - 0.5*(Dl[j] - Dl[j-2])),&ind);
			z1 = redshift_table[ind];
		}

		if(j == Nplanes-2) z2 = zsource;
		else{
			locateD(coorDist_table-1,NTABLE,(Dl[j] + 0.5*dDl[j+1]),&ind);
			z2 = redshift_table[ind];
		}

		if(flag_input_lens && j+1 == (flag_input_lens % Nplanes)){
			locateD(coorDist_table-1,NTABLE,(Dl[j] + 0.5*(Dl[j+2] - Dl[j])),&ind);
			z2 = redshift_table[ind];
		}

		/// Find which halos are in redshift range
		locateD(halo_zs-1,Nhalos,z1,&j1);
		if(j1 == Nhalos) j1 = Nhalos-1;
		locateD(halo_zs-1,Nhalos,z2,&j2);
		if(j2 == Nhalos) j2 = Nhalos-1;

		/*
		 * finding the average mass surface density in halos
		 */

		// TODO Ben: test this
		double sigma_back;
		if(field_buffer > 0.0)
			sigma_back = cosmo->haloMassInBufferedCone(min_mass*mass_scale,z1,z2,fieldofview*pow(pi/180,2),field_buffer,mass_func_type,pw_alpha)
			             /(pi*pow(sqrt(fieldofview/pi)*pi*Dl[j]/180/(1+plane_redshifts[j]) + field_buffer,2))/mass_scale;
		else sigma_back = cosmo->totalMassDensityinHalos(mass_func_type,pw_alpha,min_mass*mass_scale,plane_redshifts[j],z1,z2)/mass_scale;

		double sb=0.0;
		for(int m=0;m<j2-j1;m++){
		  sb+=halos[j1+m].mass;
		  if(r_print_halos){
			  // halo_pos is still in radians
			  double r=sqrt(halo_pos[j1+m][0]*halo_pos[j1+m][0]+halo_pos[j1+m][1]*halo_pos[j1+m][1]);
			  if(r < r_print_halos)
				  file_area << halo_zs[j1+m] << " " << halos[j1+m].mass << " " << halo_pos[j1+m][0] << " " << halo_pos[j1+m][0] <<endl;
		  }
		}
		if(field_buffer > 0.0) sb /= (pi*pow(sqrt(fieldofview/pi)*pi*Dl[j]/180/(1+plane_redshifts[j]) + field_buffer,2));
		else sb /= fieldofview*pow(pi/180.*Dl[j]/(1+plane_redshifts[j]),2);
		//else sb /= fieldofview*pow(pi/180.*Dl[j],2);

		std::cout << sigma_back << " " << sb << " " << sb/sigma_back - 1 << std::endl;
		if(sim_input_flag) sigma_back = sb;
		if(flag_run_multip_test) sigma_back = 0.0;

		halo_data[j].reset(new HaloData(&halos[j1],sigma_back,&halo_pos[j1],&halo_zs[j1],&halo_id[j1],j2-j1,Dl[j]/(1+plane_redshifts[j])));
		//halo_data[j].reset(new HaloData(&halos[j1],sigma_back,&halo_pos_Mpc[j1],&halo_zs[j1],&halo_id[j1],j2-j1,1.0));

		/// Use other constructor to create halo data
		std::cout << "  Building tree on plane " << j << " number of halos: " << halo_data[j]->Nhalos << std::endl;

		switch(internal_profile){
		case PowerLaw:
			halo_tree[j].reset(new QuadTreePowerLaw(pw_beta,&halo_data[j]->pos[0],halo_data[j]->Nhalos
							,halo_data[j]->halos,halo_data[j]->sigma_background));
			break;
		case NFW:
			halo_tree[j].reset(new QuadTreeNFW(&halo_data[j]->pos[0],halo_data[j]->Nhalos
							,halo_data[j]->halos,halo_data[j]->sigma_background));
			break;
		case PseudoNFW:
			halo_tree[j].reset(new QuadTreePseudoNFW(pnfw_beta,&halo_data[j]->pos[0],halo_data[j]->Nhalos
							,halo_data[j]->halos,halo_data[j]->sigma_background));
			break;
		case NSIE:
			halo_tree[j].reset(new QuadTreeNSIE(&halo_data[j]->pos[0],halo_data[j]->Nhalos
							,halo_data[j]->halos,halo_data[j]->sigma_background));
			break;
		case NFW_NSIE:
			halo_tree[j].reset(new QuadTreeNFW_NSIE(&halo_data[j]->pos[0],halo_data[j]->Nhalos
							,halo_data[j]->halos,halo_data[j]->sigma_background));
			break;
		default:
			ERROR_MESSAGE();
			cout << "There is no such case for the halo trees. Please choose from:" << endl;
			cout << "0: PowerLaw, 1: NFW, 2: PseudoNFW, 3: NSIE, 4: NFW_NSIE" << endl;
			exit(1);
			break;
		}

	}

	if(r_print_halos){
		 file_area.close();
	}

	cout << "constructed " << Nhalos << " halos" << endl;
}

/*
 * Builds the halo trees for a test without an AnaLens and
 * with just two lensing planes, used to test the convergence
 * of the ray-shooter.
 */
void MultiLens::buildHaloTrees_test(
		CosmoHndl cosmo /// the cosmology
		){
	int j, Ntot;
	unsigned long ind;
	double z1, z2;
	unsigned long j1,j2;

	std::cout << "MultiLens::buildHaloTrees zsource = " << zsource << std::endl;

	assert(plane_redshifts[Nplanes-1] == zsource);

	for(j=0,Ntot=0;j<Nplanes-1;j++){
		/*
		 * Setting the redshift range
		 * If there is a plane with an input lens on it, it is skipped over
		 * since it will not contain any halos
		 */
		if(j == 0) z1 = 0.0;
		else{
			locateD(coorDist_table-1,NTABLE,(Dl[j]-0.5*dDl[j]),&ind);
			z1 = redshift_table[ind];
		}

		if(j == Nplanes-2) z2 = zsource;
		else{
			locateD(coorDist_table-1,NTABLE,(Dl[j] + 0.5*dDl[j+1]),&ind);
			z2 = redshift_table[ind];
		}

		/// Find which halos are in redshift range
		locateD(halo_zs-1,Nhalos,z1,&j1);
		locateD(halo_zs-1,Nhalos,z2,&j2);

		double sigma_back = 0.0;

		halo_data[j].reset(new HaloData(&halos[j1],sigma_back,&halo_pos[j1],&halo_zs[j1],&halo_id[j1],j2-j1,Dl[j]/(1+plane_redshifts[j])));

		/// Use other constructor to create halo data
		std::cout << "  Building tree on plane " << j << " number of halos: " << halo_data[j]->Nhalos << std::endl;

		switch(internal_profile){
		case PowerLaw:
			halo_tree[j].reset(new QuadTreePowerLaw(pw_beta,&halo_data[j]->pos[0],halo_data[j]->Nhalos
							,halo_data[j]->halos,halo_data[j]->sigma_background));
			break;
		case NFW:
			halo_tree[j].reset(new QuadTreeNFW(&halo_data[j]->pos[0],halo_data[j]->Nhalos
							,halo_data[j]->halos,halo_data[j]->sigma_background));
			break;
		case PseudoNFW:
			halo_tree[j].reset(new QuadTreePseudoNFW(pnfw_beta,&halo_data[j]->pos[0],halo_data[j]->Nhalos
							,halo_data[j]->halos,halo_data[j]->sigma_background));
			break;
		case NSIE:
			halo_tree[j].reset(new QuadTreeNSIE(&halo_data[j]->pos[0],halo_data[j]->Nhalos
							,halo_data[j]->halos,halo_data[j]->sigma_background));
			break;
		case NFW_NSIE:
			halo_tree[j].reset(new QuadTreeNFW_NSIE(&halo_data[j]->pos[0],halo_data[j]->Nhalos
							,halo_data[j]->halos,halo_data[j]->sigma_background));
			break;
		default:
			ERROR_MESSAGE();
			cout << "There is no such case for the halo trees. Please choose from:" << endl;
			cout << "0: PowerLaw, 1: NFW, 2: PseudoNFW, 3: NSIE, 4: NFW_NSIE" << endl;
			exit(1);
			break;
		}

	}

	cout << "constructed " << Nhalos << " halos" << endl;
}

/*
 * Calculates the relative error of the rayshooter to the
 * analytical solution for a system with two lensing planes.
 * The error values are saved instead of the real physical values.
 */
void MultiLens::calc_error_test(
		unsigned long Npoints
		,Point *point
		,bool kappa_off
		){
	double alpha[2];
	KappaType kappa,gamma[3];
	double alpha0[2];
	KappaType kappa0,gamma0[3];
	double alpha1[2];
	KappaType kappa1,gamma1[3];

	double aa = charge*dDl[2]*Dl[1]/Dl[2];
	double bb = charge*(dDl[2]+dDl[1])*Dl[0]/Dl[2];
	double cc = charge*charge*dDl[2]*dDl[1]*Dl[0]/Dl[2];

	double xx[2], xx0[2], xx1[2], xx2[2];

	std::cout << aa << " " << bb << " " << cc << std::endl;

	for(unsigned long i=0; i<Npoints; i++){
		xx[0] = point[i].x[0];
		xx[1] = point[i].x[1];

		rayshooterInternal(xx,alpha,gamma,&kappa,kappa_off);

		xx0[0] = xx[0]*Dl[0]/(1+plane_redshifts[0]);
		xx0[1] = xx[1]*Dl[0]/(1+plane_redshifts[0]);

		halo_tree[0]->force2D_recur(xx0,alpha0,&kappa0,gamma0,kappa_off);
		double fac = 1/(1+plane_redshifts[0]);
		kappa0*=fac;
		gamma0[0]*=fac;
		gamma0[1]*=fac;
		gamma0[2]*=fac;

		if(flag_switch_deflection_off)
			alpha0[0] = alpha0[1] = 0.0;

		xx1[0] = (dDl[1]+dDl[0])/dDl[0]*xx0[0]*(1+plane_redshifts[0]) - charge*dDl[1]*alpha0[0];
		xx1[1] = (dDl[1]+dDl[0])/dDl[0]*xx0[1]*(1+plane_redshifts[0]) - charge*dDl[1]*alpha0[1];

		xx1[0] /= (1+plane_redshifts[1]);
		xx1[1] /= (1+plane_redshifts[1]);

		halo_tree[1]->force2D_recur(xx1,alpha1,&kappa1,gamma1,kappa_off);
		fac = 1/(1+plane_redshifts[1]);
		kappa1*=fac;
		gamma1[0]*=fac;
		gamma1[1]*=fac;
		gamma1[2]*=fac;

		if(flag_switch_deflection_off)
			alpha1[0] = alpha1[1] = 0.0;

		xx2[0] = (dDl[2]+dDl[1])/dDl[1]*xx1[0]*(1+plane_redshifts[1]) - dDl[2]/dDl[1]*xx0[0]*(1+plane_redshifts[0]) - charge*dDl[2]*alpha1[0];
		xx2[1] = (dDl[2]+dDl[1])/dDl[1]*xx1[1]*(1+plane_redshifts[1]) - dDl[2]/dDl[1]*xx0[1]*(1+plane_redshifts[0]) - charge*dDl[2]*alpha1[1];

		xx2[0] /= Dl[2];
		xx2[1] /= Dl[2];

		point[i].image->x[0] = (xx2[0]/alpha[0] - 1.0);

		point[i].image->x[1] = (xx2[1]/alpha[1] - 1.0);

		point[i].kappa = ((aa*kappa1+bb*kappa0-cc*(kappa0*kappa1+gamma0[0]*gamma1[0]+gamma0[1]*gamma1[1]))/kappa - 1.0);

		point[i].gamma[0] = ((-aa*gamma1[0]-bb*gamma0[0]+cc*(kappa0*gamma1[0]+gamma0[0]*kappa1))/gamma[0] - 1.0);

		point[i].gamma[1] = ((-aa*gamma1[1]-bb*gamma0[1]+cc*(kappa0*gamma1[1]+gamma0[1]*kappa1))/gamma[1] - 1.0);

		point[i].gamma[2] = (cc*(-gamma0[0]*gamma1[1]+gamma0[1]*gamma1[0])/gamma[2] - 1.0);

	}
}


/*
 * Calculates the lensing properties of two sets of images planes, once
 * using the standard multiplane rayshooter, and once using the halo by halo rayshooter.
 * This way one can see the relative error between the plane approximation and the
 * "analytical" solution of the lens equation.
 */
void MultiLens::calc_error_test_multi(
		unsigned long Npoints
		,Point *i_points
		,bool kappa_off
		,CosmoHndl cosmo
		){

	unsigned long i,j;

	double *Dl_halos = new double[Nhalos+1];
	double *dDl_halos = new double[Nhalos+1];

	for(j=0; j<Nhalos; j++)
		Dl_halos[j] = cosmo->coorDist(0.0,halo_zs[j]);
	Dl_halos[Nhalos] = cosmo->coorDist(0.0,zsource);
	dDl_halos[0] = Dl_halos[0];
	for(j=1; j<Nhalos+1; j++)
		dDl_halos[j] = Dl_halos[j] - Dl_halos[j-1];

	rayshooterInternal(Npoints,i_points,kappa_off);
	rayshooterInternal_halos_diff(Npoints,i_points,kappa_off,Dl_halos,dDl_halos);

	delete[] Dl_halos;
	delete[] dDl_halos;
}


/**
 * Set the coordinate distances of the planes by dividing the coordinate distance space into equal intervals
 * and then plugging the analytic input plane in between.
 *
 * After this flag_input_lens will hold the index of the plane it is on
 * In case it is on the first plane, it will hold the index Nplanes, to make
 * sure that it is not zero (i.e. not set)
 */
void MultiLens::setCoorDist(CosmoHndl cosmo){
	std:: vector<double> lD;
	int i, Np;

	if(flag_input_lens)
		Np = Nplanes-1;
	else
		Np = Nplanes;

	double Ds = cosmo->coorDist(0,zsource);

	double Dlens;
	if(flag_input_lens) Dlens = cosmo->coorDist(0,input_lens->getZlens());
	else Dlens = Ds;

	/// spaces lD equally up to the source, including 0 and Ds
	/// therefore we need Nplanes+1 values
	/// however, if there is an input plane, we will need Nplanes values, since the input plane will take up a value itself
	fill_linear(lD,Np,0.,Ds);

	/// ensures that the first plane and the last before the source plane have the same volume
	/// as all ther planes
	double dlD = lD[1]-lD[0];
	for(i=0; i<Np; i++){
	  lD[i] -= 0.5*dlD;
	}

	/// puts the input plane first if the case
	int j=0, flag=0;
	if(flag_input_lens && Dlens < lD[1]){
		Dl[j] = Dlens;
		flag_input_lens = (InputLens)Nplanes;
		flag = 1;
		j++;
	}

	/// assigns the redshifts and plugs in the input plane
	for(i=1; i<Np; i++){
		Dl[j] = lD[i];

		if(flag_input_lens && flag == 0)
			if(Dlens > lD[i] && Dlens <= lD[i+1]){
				Dl[j] = lD[i];
				Dl[++j] = Dlens;
				flag_input_lens = (InputLens)j;
				flag = 1;
			}
		j++;
	}

	Dl[Nplanes-1] = Ds;

	dDl[0] = Dl[0];  // distance between jth plane and the previous plane
	for(j = 1; j < Nplanes; j++){
		dDl[j] = Dl[j] - Dl[j-1]; // distance between jth plane and the previous plane
	}

	if(flag_input_lens)
		cout << "zlens " << input_lens->getZlens() << " on plane number " << (flag_input_lens % Nplanes) << endl;

	cout << "Dl: ";
	for(j = 0; j < Nplanes; j++)
		cout << Dl[j] << " ";
	cout << endl;

	cout << "dDl: ";
	for(j = 0; j < Nplanes; j++)
		cout << dDl[j] << " ";
	cout << endl;

	unsigned long k;
	// assigns the redshifts and plugs in the input plane
	cout << "z: ";
	for(i=0; i<Nplanes-1; i++){
		locateD(coorDist_table-1,NTABLE,Dl[i],&k);
		plane_redshifts[i] = redshift_table[k];
		cout << plane_redshifts[i] << " ";
	}
	plane_redshifts[Nplanes-1] = zsource;
	cout << plane_redshifts[i] << " " << std::endl;
}

double MultiLens::getZlens(){
	if(flag_input_lens)
		return input_lens->getZlens();
	else
		return plane_redshifts[0];
}

void MultiLens::setZlens(CosmoHndl cosmo,double z,double zsource){
	if(flag_input_lens)
		input_lens->setZlens(cosmo,z,zsource);
	else{
		cout << "It is not possible to reset the redshift of the input AnaLens plane a MultiLens" << endl;
		ERROR_MESSAGE();
		exit(1);
	}
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
void MultiLens::readInputSimFile(CosmoHndl cosmo){

	double ra,dec,z,vmax,vdisp,r_halfmass;
	unsigned long i,j;
	unsigned long haloid,idd,np;
	double mo=7.3113e10,M1=2.8575e10,gam1=7.17,gam2=0.201,beta=0.557;

	double rmax=0,rtmp=0;

	//int index;

	if(internal_profile == PseudoNFW){
		ERROR_MESSAGE();
		std::cout << "Input to MultiLens is not yet enabled for PseudoNFW profiles"
				<< std::endl << "Change this is parameter file" << std::endl;
		exit(1);
	}

	ifstream file_in(input_sim_file.c_str());
	if(!file_in){
		cout << "Can't open file " << input_sim_file << endl;
		exit(1);
	}

	// skip through header information in data file
	i=0;
	while(file_in.peek() == '#'){
		file_in.ignore(10000,'\n');
		++i;
	}
	std::cout << "skipped "<< i << " comment lines in file " << input_sim_file << std::endl;

	std::vector<HaloStructure> halo_vec;
	std::vector<double> halo_zs_vec;
	std::vector<double *> halo_pos_vec;
	std::vector<unsigned long> halo_id_vec;

	// read in data
	int j_max;
	double mass_max=0,R_max=0,V_max=0,minmass=1e30;
	double *theta;
	HaloStructure halo;
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
		getline(file_in,myline);

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
		theta[0] = ra*pi/180.;
		theta[1] = dec*pi/180.;

		if(partial_cone){
			double r = sqrt(theta[0]*theta[0]+theta[1]*theta[1]);
			if(r > 1.5*mokalens->map->boxlrad)
				continue;
		}

		//file_in >> haloid >>  idd >>  ra >>  dec >>  z
		//		 >>  np >>  vdisp >>  vmax >>  r_halfmass;
		//file_in >> c >> haloid >> c >> idd >> c >> ra >> c >> dec >> c >> z
		//		 >> c >> np >> c >> vdisp >> c >> vmax >> c >> r_halfmass >> c;
		//cout << haloid << c << idd << c << ra << c << dec << c << z
		//				 << c << np << c << vdisp << c << vmax << c << r_halfmass << std::endl;
		//std::cout << i << "  z: " << z << " np: " << np << " vmax :" << vmax << " vdisp: " << vdisp << "  " << file_in.peek() << std::endl;

		if(np > 0.0 && vdisp > 0.0 && z <= zsource){

			halo_pos_vec.push_back(theta);
			
			halo_vec.push_back(halo);

			halo_vec[j].mass = np*8.6e8/cosmo->gethubble();
			halo_vec[j].Rmax = halo_vec[j].mass*Grav/2/pow(vmax/lightspeed,2);  // SIS value

			if(halo_vec[j].mass > mass_max) {
				mass_max = halo_vec[j].mass;
				j_max = j;
			}
			if(halo_vec[j].mass < minmass) {
				minmass = halo_vec[j].mass;
			}

			if(internal_profile == NFW || internal_profile == NFW_NSIE){
				NFW_Utility nfw_util;

				// Find the NFW profile with the same mass, Vmax and R_halfmass
				nfw_util.match_nfw(vmax,r_halfmass*cosmo->gethubble(),halo_vec[j].mass
						   ,&(halo_vec[j].rscale),&(halo_vec[j].Rmax));
				halo_vec[j].rscale = halo_vec[j].Rmax/halo_vec[j].rscale; // Was the concentration
				//std::cout << "z: " << z << " vmax:" << vmax << " Rmax: " << halo_vec[j].Rmax << " rscale: " << halo_vec[j].rscale << std::endl;
			}
			
			if(internal_profile == NSIE || internal_profile == NFW_NSIE){
				NFW_Utility nfw_util;

				if(internal_profile == NFW_NSIE){
					// Stellar mass fraction in from Moster et al. 2010
					galaxy_mass_fraction = 2*mo*pow(halo_vec[j].mass/M1,gam1)
					  /pow(1+pow(halo_vec[j].mass/M1,beta),(gam1-gam2)/beta)/halo_vec[j].mass;
					if(galaxy_mass_fraction > 1.0) galaxy_mass_fraction = 1;

					halo_vec[j].mass_nsie = halo_vec[j].mass*galaxy_mass_fraction;   //TODO This is a kluge. A mass dependent ratio would be better
					halo_vec[j].mass *= (1-galaxy_mass_fraction);
					halo_vec[j].rcore_nsie = 0.0;
				}else{
					halo_vec[j].mass_nsie = halo_vec[j].mass;
					halo_vec[j].rscale = 0.0;
					halo_vec[j].mass = 0.0;
					halo_vec[j].rcore_nsie = 0.0;
				}

//				halo_vec[j].sigma_nsie = vmax/sqrt(2.0);   //TODO This is a kluge.
				halo_vec[j].sigma_nsie = 126*pow(halo_vec[j].mass_nsie/1.0e10,0.25); // From Tully-Fisher and Bell & de Jong 2001
				halo_vec[j].fratio_nsie = (ran2(seed)+1)*0.5;  //TODO This is a kluge.
				halo_vec[j].pa_nsie = 2*pi*ran2(seed);  //TODO This is a kluge.
				halo_vec[j].Rsize_nsie = rmaxNSIE(halo_vec[j].sigma_nsie,halo_vec[j].mass_nsie
						,halo_vec[j].fratio_nsie,halo_vec[j].rcore_nsie);
				cout << halo_vec[j].Rmax << " " << halo_vec[j].Rsize_nsie << " " << vmax << " " << halo_vec[j].mass_nsie << " " << halo_vec[j].fratio_nsie << " " << halo_vec[j].rcore_nsie << endl;
				//std::cout << "sigma_nsie: " << halo_vec[j].sigma_nsie << " fratio_nsie: " << halo_vec[j].fratio_nsie
				//		<< " pa_nsie: " << halo_vec[j].pa_nsie << " Rsize_nsie: " << halo_vec[j].Rsize_nsie << std::endl;

				if(internal_profile == NSIE)
					halo_vec[j].Rmax = MAX(1.0,1.0/halo_vec[j].fratio_nsie)*halo_vec[j].Rsize_nsie;  // redefine
				assert(halo_vec[j].Rmax >= halo_vec[j].Rsize_nsie); // If this should allowable then the tree Rmax and the NFW profile Rmax need to be separated.

			}else{
				// initialize unused variables to harmless values
				halo_vec[j].Rsize_nsie = halo_vec[j].Rmax;
				halo_vec[j].rcore_nsie = halo_vec[j].Rmax;
				halo_vec[j].fratio_nsie = 0;
				halo_vec[j].pa_nsie = 0;
				halo_vec[j].sigma_nsie = 0;
			}

			// initialize unused variables to harmless values in PowerLaw case
			if(internal_profile == PowerLaw) halo_vec[j].rscale = 0.0;

			if(halo_vec[j].Rmax > R_max) R_max = halo_vec[j].Rmax;
			if(vdisp > V_max) V_max = vdisp;
			/*
			halo_vec[j].mass = mass*1.0e10*cosmo->gethubble();
			halo_vec[j].Rmax = cosmo->R200(z,mass*1.0e10*cosmo->gethubble());
			assert(halo_vec[j].Rmax > 0.0);
			cout << "Rmax:" << halo_vec[j].Rmax << endl;
			halo_vec[j].rscale = halo_vec[j].Rmax/cosmo->NFW_Concentration(vmax,halo_vec[j].mass,halo_vec[j].Rmax);
			 */

			halo_zs_vec.push_back(z);
			halo_id_vec.push_back(haloid);

			halo_vec[j].mass /= mass_scale;

			if(rmax < (rtmp = theta[0]*theta[0]+theta[1]*theta[1])) rmax = rtmp;
			++j;
		
		}
	}
	file_in.close();
	std::cout << halo_vec.size() << " halos read in."<< std::endl
			<< "Max input mass = " << mass_max << "  R max = " << R_max << "  V max = " << V_max << std::endl;

	Nhalos = halo_vec.size();

	/// setting the minimum halo mass in the simulation
	min_mass = minmass;
	if(field_buffer > 0.0){
		std::cout << "Overiding field_buffer to make it 0 because halos are read in." << endl;
		field_buffer = 0.0;
	}

	if(partial_cone == false)
		fieldofview = pi*rmax*pow(180/pi,2);  // Resets field of view to estimate of inputed one

	halos = new HaloStructure[Nhalos];
	halo_zs = new double[Nhalos];
	halo_id = new unsigned long[Nhalos];
	halo_pos = PosTypeMatrix(Nhalos,3);

	for(i=0;i<Nhalos;++i){
		halo_id[i] = halo_id_vec[i];
		halo_zs[i] = halo_zs_vec[i];
		halo_pos[i] = halo_pos_vec[i];
		halos[i] = halo_vec[i];
	}

	std::cout << "sorting in MultiLens::readInputSimFile()" << std::endl;
	// sort the halos by readshift
	MultiLens::quicksort(halos,halo_pos,halo_zs,halo_id,Nhalos);

	std::cout << "leaving MultiLens::readInputSimFile()" << std::endl;

	read_sim_file = true;
}


/**
 * Sets the internal parameters of the multiple lens model
 * first the redshifts of the planes are calculated
 * then the coordinate distances to the different planes
 * the planes are populated by halos and the halo trees are built
 * and finally the internal parameters of the input plane are set, in case there is one
 */
void MultiLens::setInternalParams(CosmoHndl cosmo, SourceHndl source){

	if( (cosmo->getOmega_matter() + cosmo->getOmega_lambda()) != 1.0 ){
		printf("ERROR: MultiLens can only handle flat universes at present.  Must change cosmology.\n");
		exit(1);
	}

	/// makes the oordinate distance table for the calculation of the redshifts of the different planes
	if(table_set == false) {std::cout << "making tables" << std::endl; make_table(cosmo);}
	
	if(flag_input_lens)
	  input_lens->setInternalParams(cosmo,source);
	
	setCoorDist(cosmo);

	if(sim_input_flag){
		if(read_sim_file == false) readInputSimFile(cosmo);
	}
	else{
		// TODO Ben swap function here or provide toggle
		if(field_buffer > 0.0) createHaloData_buffered(cosmo,seed);
		else{
			if(flag_run_twop_test) createHaloData_test(cosmo,seed);
			else createHaloData(cosmo,seed);
		}
	}

	if(flag_run_twop_test) buildHaloTrees_test(cosmo);
	else buildHaloTrees(cosmo);
	std:: cout << " done " << std:: endl;
}

/// Sort halos[], brr[][], and id[] by content off arr[]
void MultiLens::quicksort(HaloStructure *halos,double **brr,double *arr,unsigned long  *id,unsigned long N){
	double pivotvalue;
	unsigned long pivotindex,newpivotindex,i;

	if(N <= 1) return ;

	// pick pivot as the median of the first, last and middle values
	if ((arr[0] >= arr[N/2] && arr[0] <= arr[N-1])
			|| (arr[0] >= arr[N-1] && arr[0] <= arr[N/2])) pivotindex = 0;
	else if ((arr[N/2] >= arr[0] && arr[N/2] <= arr[N-1])
			|| (arr[N/2] >= arr[N-1] && arr[N/2] <= arr[0])) pivotindex = N/2;
	else pivotindex = N-1;
	pivotvalue=arr[pivotindex];

	// move pivet to end of array
	std::swap(arr[pivotindex],arr[N-1]);
	std::swap(halos[pivotindex],halos[N-1]);
	std::swap(brr[pivotindex],brr[N-1]);
	std::swap(id[pivotindex],id[N-1]);
	newpivotindex=0;

	// partition list and array
	for(i=0;i<N;++i){
		if(arr[i] <= pivotvalue){
			std::swap(arr[newpivotindex],arr[i]);
			std::swap(halos[newpivotindex],halos[i]);
			std::swap(brr[newpivotindex],brr[i]);
			std::swap(id[newpivotindex],id[i]);
			++newpivotindex;
		}
	}
	--newpivotindex;

	quicksort(halos,brr,arr,id,newpivotindex);
	quicksort(&halos[newpivotindex+1],&brr[newpivotindex+1],&arr[newpivotindex+1],&id[newpivotindex+1],N-newpivotindex-1);

	return ;
}

/**
 * \brief Changes the maximum redshift that the rays are shot to.
 *
 * The multilens must have been initially constructed with a source redshift that is higher
 * than this redshift.  This is used to rayshoot to a source whose line of sight passes through the
 * simulation volume.  The source can be at higher redshift than the simulation volume.
 *
 * To revert the source redshift to its original value use MultiLens::RevertSourcePlane().
 *
 */
short MultiLens::ResetSourcePlane(
		CosmoHndl cosmo           /// cosmology
		,double z                 /// redshift of implanted source
		,bool nearest             /** If true, set the source plane to the nearest (in coordinate distance)
			                      * lensing plane that was created already.  This can be used to avoid self-lensing
			                      * by the halo of the source.  If the source is at higher redshift than the simulation
			                      * volume the source will be at its real redshift.
			                      */
		,unsigned long GalID
		,double *xx
		){
	unsigned long j;
	short out=0;

	toggle_source_plane = true;

	if(z<=0.0){
		cout << "Warning: Source redshift can't be set to " << z << " in MultiLens::ResetSourcePlane." << endl;
		return out;
	}

	// j is the index of the next plane at higher redshift, This plane will be temporarily replaced and used as a source plane

	double Ds = cosmo->coorDist(0,z);
	locateD(Dl-1,Nplanes,Ds,&j);
	assert(j <= Nplanes && j >=0);

	if(j >= Nplanes-1){
	  j--;
	}
	else if(j > 0){
		unsigned long ind;

		locateD(coorDist_table-1,NTABLE,(Dl[j]-0.5*dDl[j]),&ind);
		double z1 = redshift_table[ind];
		if(nearest) j = (z>=z1) ? j : j-1;  // reset j to the nearest plane
	}

	if(nearest && (j < Nplanes-1) ){
		zs_implant = plane_redshifts[j];
		Ds_implant = Dl[j];
		if(j > 0) dDs_implant = dDl[j];
		else  dDs_implant = Ds_implant;
	}else{
		// if nearest==false or the source is at higher redshift than the last plane use the real redshift
		Ds_implant = Ds;
		zs_implant = z;
		if(j > 0) dDs_implant = cosmo->coorDist(plane_redshifts[j-1],z);
		else  dDs_implant = Ds;
	}

	std::cout << "Source on plane " << j << " zs " << zs_implant << " Ds " << Ds << " dDs " << dDs_implant << std::endl;
	index_of_new_sourceplane = j;

	out = j;
	return out;
}

void MultiLens::unusedHalos(){
  ofstream file_area("halos_with.data.txt");
  if(!file_area){
    cout << "unable to create file " << endl;
    exit(1);
  }
  
  for(int l=0; l < Nplanes-1; l++){
    for(int i=0; i<halo_data[l]->Nhalos; i++){
      if(halo_data[l]->haloID[i] == 0){
	halo_data[l]->haloID[i] = 0;
	file_area << halo_data[l]->halos[i].mass << " " << halo_data[l]->z[i] 
		  << " " << halo_data[l]->pos[i][0] << " " << halo_data[l]->pos[i][1] << endl;
      }
    }
  }

  file_area.close();

}


