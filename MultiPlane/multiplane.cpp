/*
 * halo_model.cpp
 *
 *  Created on: Jan 14, 2012
 *      Author: mpetkova, bmetcalf
 */

#include "slsimlib.h"

using namespace std;


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
  delete[] halo_tree;

  delete[] Dl;
  delete[] plane_redshifts;
  delete[] dDl;

  Nplanes = Np;

  plane_redshifts = new double[Nplanes];
  Dl = new double[Nplanes];
  dDl = new double[Nplanes];

  halo_tree = new auto_ptr<QuadTree>[Nplanes-1];

  setCoorDist(cosmo);
  buildHaloTrees(cosmo);
}

void MultiLens::resetHalos(CosmoHndl cosmo){
  delete[] halo_tree;
  halos.resize(0);

  delete[] halo_zs;
  delete[] halo_id;
  Utilities::free_PosTypeMatrix(halo_pos,Nhalos,3);
  if(flag_run_multip_test)
	  Utilities::free_PosTypeMatrix(halo_pos_Mpc,Nhalos,3);

  halo_tree = new auto_ptr<QuadTree>[Nplanes-1];

  switch(internal_profile){
  	case PowerLaw:
  		if(flag_run_twop_test) createHaloData_test<PowerLawLensHalo>(cosmo,seed);
  		else createHaloData<PowerLawLensHalo>(cosmo,seed);
  		break;
  	case NFW:
  		if(nfw_table_set == false){
  			make_nfw_tables();
  			nfw_table_set = true;
  		}
  		if(flag_run_twop_test) createHaloData_test<NFWLensHalo>(cosmo,seed);
  		else createHaloData<NFWLensHalo>(cosmo,seed);
  		break;
  	case PseudoNFW:
  		if(nfw_table_set == false){
  			make_pnfw_tables(beta);
  			pnfw_table_set = true;
  		}
  		if(flag_run_twop_test) createHaloData_test<PseudoNFWLensHalo>(cosmo,seed);
  		else createHaloData<PseudoNFWLensHalo>(cosmo,seed);
  		break;
  	case NSIE:
  		if(flag_run_twop_test) createHaloData_test<NSIELensHalo>(cosmo,seed);
  		else createHaloData<NSIELensHalo>(cosmo,seed);
  		break;
  	case NFW_NSIE:
  		if(nfw_table_set == false){
  			make_nfw_tables();
  			nfw_table_set = true;
  		}
  		if(flag_run_twop_test) createHaloData_test<NFW_NSIELensHalo>(cosmo,seed);
  		else createHaloData<NFW_NSIELensHalo>(cosmo,seed);
  		break;
  }

  if(flag_run_twop_test) buildHaloTrees_test(cosmo);
  else buildHaloTrees(cosmo);
}

/*
 * \ingroup Constructor
 * allocates space for the halo trees and the inout lens, if there is any
 */
MultiLens::MultiLens(InputParams& params,long *my_seed) : Lens(){

	second_halo = false;

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

	nfw_table_set = false;
	pnfw_table_set = false;

	seed = my_seed;
}

MultiLens::~MultiLens(){
	delete[] halo_tree;

	delete[] Dl;
	delete[] plane_redshifts;
	delete[] dDl;

	halos.resize(0);
	delete[] halo_zs;
	delete[] halo_id;
	Utilities::free_PosTypeMatrix(halo_pos,Nhalos,3);
	if(flag_run_multip_test)
		Utilities::free_PosTypeMatrix(halo_pos_Mpc,Nhalos,3);

	if(flag_input_lens)
		delete input_lens;

	delete[] coorDist_table;
	delete[] redshift_table;

	if(nfw_table_set){
		delete[] xtable;
		delete[] gtable;
		delete[] ftable;
		delete[] g2table;
	}

	if(pnfw_table_set){
		delete[] xtable;
		delete[] mhattable;
	}
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
		second_halo = true;
		if(!params.get("galaxy_mass_fraction",galaxy_mass_fraction)){
			  ERROR_MESSAGE();
			  cout << "to construct a NFW + NSIE model the parameter galaxy_mass_fraction needs to be set in the parameter file " << params.filename() << endl;
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
	if(!params.get("internal_slope",beta) && internal_profile == PowerLaw)     beta = -1.0;
	if(!params.get("internal_slope",beta) && internal_profile == PseudoNFW) beta = 2.0;
	if(!params.get("deflection_off",flag_switch_deflection_off)) flag_switch_deflection_off = false;
	if(!params.get("background_off",flag_switch_background_off)) flag_switch_background_off = false;
	if(!params.get("twop_test",flag_run_twop_test)) flag_run_twop_test = false;
	if(!params.get("multip_test",flag_run_multip_test)) flag_run_multip_test = false;
	if(!params.get("print_halos",r_print_halos)) r_print_halos = 0.0;

	// Some checks for valid parameters
	  if(internal_profile == PowerLaw && beta >= 0){
		  ERROR_MESSAGE();
		  cout << "Power Law internal slope >=0 not possible." << endl;
		  exit(1);
	  }

	  if(internal_profile == PseudoNFW && beta <= 0){
		  ERROR_MESSAGE();
		  cout << "Pseudo NFW internal slope <=0 not possible." << endl;
		  exit(1);
	  }

	  if(internal_profile == PseudoNFW && (beta / floor(beta) > 1.0)){
		  ERROR_MESSAGE();
		  cout << "Pseudo NFW internal slope needs to be a whole number." << endl;
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
		cout << "  slope: " << beta << endl;
		break;
	case NFW:
		cout << "  NFW internal profile " << endl;
		break;
	case PseudoNFW:
		cout << "  Pseudo NFW internal profile " << endl;
		cout << "  slope: " << beta << endl;
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
		double sigma_back = cosmo->haloMassInBufferedCone(min_mass*mass_scale,z1,z2,fieldofview*pow(pi/180,2),field_buffer,mass_func_type,pw_alpha)
			             /(pi*pow(sqrt(fieldofview/pi)*pi*Dl[j]/180/(1+plane_redshifts[j]) + field_buffer,2))/mass_scale;

		double sb=0.0;
		for(int m=0;m<j2-j1;m++){
		  sb+=halos[j1+m]->get_mass();
		  if(r_print_halos){
			  // halo_pos is still in radians
			  double r=sqrt(halo_pos[j1+m][0]*halo_pos[j1+m][0]+halo_pos[j1+m][1]*halo_pos[j1+m][1]);
			  if(r < r_print_halos)
				  file_area << halo_zs[j1+m] << " " << halos[j1+m]->get_mass() << " " << halo_pos[j1+m][0] << " " << halo_pos[j1+m][0] <<endl;
		  }
		}
		sb /= (pi*pow(sqrt(fieldofview/pi)*pi*Dl[j]/180/(1+plane_redshifts[j]) + field_buffer,2));

		std::cout << sigma_back << " " << sb << " " << sb/sigma_back - 1 << std::endl;
		if(sim_input_flag) sigma_back = sb;
		if(flag_run_multip_test) sigma_back = 0.0;

		/// Use other constructor to create halo data
		std::cout << "  Building tree on plane " << j << " number of halos: " << j2-j1 << std::endl;

		halo_tree[j].reset(new QuadTree(&halo_pos[j1],halos[j1],j2-j1,sigma_back));
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

		/// Use other constructor to create halo data
		std::cout << "  Building tree on plane " << j << " number of halos: " << j2-j1 << std::endl;

		halo_tree[j].reset(new QuadTree(&halo_pos[j1],halos[j1],j2-j1,sigma_back));
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

	rayshooterInternal_halos(Npoints,i_points,kappa_off,Dl_halos,dDl_halos);

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

	switch(internal_profile){
	case PowerLaw:
		if(sim_input_flag){
			if(read_sim_file == false) readInputSimFile<PowerLawLensHalo,NSIELensHalo>(cosmo);
		}
		else{
			if(flag_run_twop_test) createHaloData_test<PowerLawLensHalo>(cosmo,seed);
			else createHaloData<PowerLawLensHalo>(cosmo,seed);
		}
		break;
	case NFW:
  		if(nfw_table_set == false){
  			make_nfw_tables();
  			nfw_table_set = true;
  		}
		if(sim_input_flag){
			if(read_sim_file == false) readInputSimFile<NFWLensHalo,NSIELensHalo>(cosmo);
		}
		else{
			if(flag_run_twop_test) createHaloData_test<NFWLensHalo>(cosmo,seed);
			else createHaloData<NFWLensHalo>(cosmo,seed);
		}
		break;
	case PseudoNFW:
  		if(nfw_table_set == false){
  			make_pnfw_tables(beta);
  			pnfw_table_set = true;
  		}
		if(sim_input_flag){
			if(read_sim_file == false) readInputSimFile<PseudoNFWLensHalo,NSIELensHalo>(cosmo);
		}
		else{
			if(flag_run_twop_test) createHaloData_test<PseudoNFWLensHalo>(cosmo,seed);
			else createHaloData<PseudoNFWLensHalo>(cosmo,seed);
		}
		break;
	case NSIE:
		if(sim_input_flag){
			if(read_sim_file == false) readInputSimFile<NSIELensHalo,NSIELensHalo>(cosmo);
		}
		else{
			if(flag_run_twop_test) createHaloData_test<NSIELensHalo>(cosmo,seed);
			else createHaloData<NSIELensHalo>(cosmo,seed);
		}
		break;
	case NFW_NSIE:
  		if(nfw_table_set == false){
  			make_nfw_tables();
  			nfw_table_set = true;
  		}
		if(sim_input_flag){
			if(read_sim_file == false) readInputSimFile<NFWLensHalo,NSIELensHalo>(cosmo);
		}
		else{
			if(flag_run_twop_test) createHaloData_test<NFW_NSIELensHalo>(cosmo,seed);
			else createHaloData<NFW_NSIELensHalo>(cosmo,seed);
		}
		break;
	}

	if(flag_run_twop_test) buildHaloTrees_test(cosmo);
	else buildHaloTrees(cosmo);
	std:: cout << " done " << std:: endl;
}

/// Sort halos[], brr[][], and id[] by content off arr[]
void MultiLens::quicksort(LensHaloHndl *halos,double **brr,double *arr,unsigned long  *id,unsigned long N){
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

	quicksort(&halos[0],brr,arr,id,newpivotindex);
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



