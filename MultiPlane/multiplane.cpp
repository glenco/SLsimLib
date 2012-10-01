/*
 * halo_model.cpp
 *
 *  Created on: Jan 14, 2012
 *      Author: mpetkova
 */

#include <slsimlib.h>
#include <sstream>
#include <string>
#include <utilities.h>
#include <source.h>
#include <sourceAnaGalaxy.h>

using namespace std;

HaloData::HaloData(CosmoHndl cosmo
		,double mass
		,double zlens){
	allocation_flag = true;
	kappa_background = 0.0;
	Nhalos = 1;
	halos = new HaloStructure[Nhalos];
	pos = PosTypeMatrix(0,Nhalos-1,0,2);

	HALO *ha = new HALO(cosmo,mass,zlens);

	int i;
	for(i = 0; i < Nhalos; i++){
		pos[i][0] = 0.0;
		pos[i][1] = 0.0;
		pos[i][2] = 0.0;

		halos[i].mass = mass;
		halos[i].Rmax = ha->getRvir()*cosmo->gethubble();
		halos[i].rscale = halos[i].Rmax/ha->getConcentration(0);

		cout << endl << "mass: " << halos[i].mass  << ", Rmax: " << halos[i].Rmax << ", rscale: " << halos[i].rscale << endl;
	}

	delete ha;
}

HaloData::HaloData(CosmoHndl cosmo,HaloStructure *halostrucs,PosType **positions
		   , double *zz, unsigned long *id, unsigned long Nhaloss,double z1,double z2,double mass_scale):
	pos(positions), halos(halostrucs), Nhalos(Nhaloss),z(zz)
{
	allocation_flag = false;

	// Set background convergence
	double min_mass = halos[0].mass;
	for(long i=0;i<Nhalos;++i)
		if(halos[i].mass < min_mass) min_mass = halos[i].mass;
	double zave =(z1+z2)/2;
	HALO *ha = new HALO(cosmo,min_mass,zave);

	// calculate the mass density on the plane
	// and convert it to physical 1/physical_distance^2
	kappa_background = ha->totalMassDensityinHalos(0,1)*pow(cosmo->gethubble(),2)*pow(1+zave,3);
	kappa_background *= cosmo->getOmega_matter()*cosmo->angDist(z1,z2)/mass_scale;
}

/*HaloData::HaloData(NSIEstructure *halostrucs,PosType **positions,unsigned long Nhaloss):
	pos(positions), halos(NULL), nsiehalos(halostrucs), Nhalos(Nhaloss)
{
	allocation_flag = false;
	kappa_background = 0.0;  //TODO This should be set properly at some point.
}*/

/**
 *  \brief In this constructor the halos are created from a mass function and given random positions on the lens plane
 */

HaloData::HaloData(
		double fov            /// field of view in square degrees
		,double min_mass      /// Minimum mass of a halo
		,double mass_scale    /// mass scale
		,double z1            /// lowest redshift
		,double z2            /// highest redshift
		,int mass_func_type   /// mass function type: 0 Press-Schechter, 1 Sheth-Tormen, 2 Power Law
		,double alpha		/// slope of the Power Law
		,CosmoHndl cosmo     /// cosmology
		,long *seed
	)
	{

	allocation_flag = true;

    HALO *ha = new HALO(cosmo,min_mass,(z1+z2)/2);

	// calculate the mass density on the plane
    // and convert it to physical 1/physical_distance^2
	kappa_background = ha->totalMassDensityinHalos(0,1,alpha)*pow(cosmo->gethubble(),2)*pow(1+(z1+z2)/2,3);
	kappa_background *= cosmo->getOmega_matter()*cosmo->angDist(z1,z2)/mass_scale;

    //int iterator;
    //std::vector<double> vmasses,vsizes,vscale,vz;
    std::vector<int> vindex;
	std::vector<double> Logm,Nhalosbin;
	std::vector<double> Dli;

    Logm.resize(Nmassbin);
    Nhalosbin.resize(Nmassbin);

	/* fill the log(mass) vector */

	fill_linear(Logm,Nmassbin,min_mass,MaxLogm);

	double Nhaloestot;

	Nhalosbin[0] = cosmo->haloNumberDensityOnSky(pow(10,Logm[0])/cosmo->gethubble(),z1,z2,mass_func_type,alpha)*fov;

	Nhaloestot = Nhalosbin[0];
	Nhalosbin[0] = 1;
	int k;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(k)
#endif
	for(k=1;k<Nmassbin;k++){
		// cumulative number density in one square degree
		Nhalosbin[k] = cosmo->haloNumberDensityOnSky(pow(10,Logm[k])/cosmo->gethubble(),z1,z2,mass_func_type,alpha)*fov;
		// normalize the cumulative distribution to one
		Nhalosbin[k] = Nhalosbin[k]/Nhaloestot;
	}


	Nhalos = (long)(poidev(float(Nhaloestot), seed) + 0.5);

	halos = new HaloStructure[Nhalos];
	pos = PosTypeMatrix(0,Nhalos-1,0,2);
	double rr,theta,maxr,zi;
	unsigned long i;
	for(i = 0; i < Nhalos; i++){
		zi = z1+(z2-z1)*ran2 (seed);
		Dli.push_back(cosmo->angDist(0,zi));
		maxr = pi*sqrt(fov/pi)/180*Dli[i]; // fov is a circle
		rr = maxr*sqrt(ran2(seed));
		theta = 2*pi*ran2(seed);
		pos[i][0] = rr*cos(theta);
		pos[i][1] = rr*sin(theta);

		halos[i].mass = pow(10,InterpolateYvec(Nhalosbin,Logm,ran2 (seed)));
		ha->reset(halos[i].mass,zi);
		halos[i].mass /= mass_scale;
		//halos[i].mass = vmasses[i];
		//halos[i].Rmax = vsizes[i];
		halos[i].Rmax = ha->getRvir()*cosmo->gethubble();
		//halos[i].rscale = vsizes[i]/vscale[i]; // get the Rscale=Rmax/c
		halos[i].rscale = halos[i].Rmax/ha->getConcentration(0); // get the Rscale=Rmax/c
		pos[i][2] = 0.0;//halos[i].Rmax;
		/* -> if you want to print the masses and redshift
		std:: cout << halos[i].mass << "   " << zi << std:: endl; */		
	}

	delete ha;
}

HaloData::~HaloData(){

	if(allocation_flag){
		free_PosTypeMatrix(pos,0,Nhalos-1,0,2);
		delete[] halos;
	}
}

void MultiLens::make_table(CosmoHndl cosmo){
	int i;
	double x, dx = 10.0/(double)NTABLE;

	coorDist_table = new double[NTABLE];
	redshift_table = new double[NTABLE];
	
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,x)
#endif
	for(i = 0 ; i< NTABLE; i++){
		x = (i+1)*dx;
		redshift_table[i] = x;
		coorDist_table[i] = cosmo->coorDist(0,x);
	}
	table_set=true;
}

/*
 * \ingroup Constructor
 * allocates space for the halo trees and the inout lens, if there is any
 */
MultiLens::MultiLens(string filename,long *my_seed) : Lens(){
	readParamfile(filename);

	NTABLE=1000;
	table_set = false;

	// flag to determine if halos are created randomly or read in from a external simulation.
	if(input_sim_file.size() < 1) sim_input_flag = false;
	else sim_input_flag = true;

	std::cout << input_sim_file.c_str() << std::endl;

	if(input_gal_file.size() < 1) gal_input_flag = false;
	else gal_input_flag = true;

	std::cout << input_gal_file.c_str() << std::endl;

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
		input_lens = new AnaLens(filename);
		analens = static_cast<AnaLens*>(input_lens);
		break;
	case moka_lens:
		input_lens = new MOKALens(filename);
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

	if(sim_input_flag){  // Otherwise these deallocations are done in the destructor of halo_data's
		delete[] halos;
		delete[] halo_zs;
		delete[] halo_id;
		free_PosTypeMatrix(halo_pos,0,Nhalos-1,0,2);
	}

	if(flag_input_lens)
		delete input_lens;

	delete[] coorDist_table;
	delete[] redshift_table;
}

void MultiLens::readParamfile(string filename){
      const int MAXPARAM = 50;
	  string label[MAXPARAM], rlabel, rvalue;
	  void *addr[MAXPARAM];
	  int id[MAXPARAM];
	  stringstream ss;
	  int i ,n;
	  int myint;
	  double mydouble;
	  string mystring;
	  char dummy[100];
	  string escape = "#";
	  int flag;

	  n = 0;

	  /// id[] = 2 = string, 1 = int, 0 = double
	  addr[n] = &outputfile;
	  id[n] = 2;
	  label[n++] = "outputfile";

	  addr[n] = &Nplanes;
	  id[n] = 1;
	  label[n++] = "Nplanes";

	  addr[n] = &min_mass;
	  id[n] = 0;
	  label[n++] = "min_mass";

	  addr[n] = &mass_scale;
	  id[n] = 0;
	  label[n++] = "mass_scale";

	  addr[n] = &flag_input_lens;
	  id[n] = 1;
	  label[n++] = "flag_input_lens";

	  addr[n] = &fieldofview;
	  id[n] = 0;
	  label[n++] = "fov";

	  addr[n] = &mass_func_type;
	  id[n] = 1;
	  label[n++] = "mass_func_type";

	  addr[n] = &internal_profile;
	  id[n] = 1;
	  label[n++] = "internal_profile";

	  addr[n] = &input_sim_file;
	  id[n] = 2;
	  label[n++] = "input_simulation_file";

	  addr[n] = &input_gal_file;
	  id[n] = 2;
	  label[n++] = "input_galaxy_file";

	  addr[n] = &pw_alpha;
	  id[n] = 0;
	  label[n++] = "alpha";

	  addr[n] = &pw_beta;
	  id[n] = 0;
	  label[n++] = "internal_slope_pw";

	  addr[n] = &pnfw_beta;
	  id[n] = 0;
	  label[n++] = "internal_slope_pnfw";

	  addr[n] = &flag_switch_deflection_off;
	  id[n] = 1;
	  label[n++] = "deflection_off";

	  addr[n] = &zsource;
	  id[n] = 0;
	  label[n++] = "z_source";

	  assert(n < MAXPARAM);

	  cout << "Multi lens: reading from " << filename << endl;

	  ifstream file_in(filename.c_str());
	  if(!file_in){
	    cout << "Can't open file " << filename << endl;
	    exit(1);
	  }

	  // output file
	  while(!file_in.eof()){
		  file_in >> rlabel >> rvalue;
		  file_in.getline(dummy,100);

		  if(rlabel[0] == escape[0])
			  continue;

		  flag = 0;

		  for(i = 0; i < n; i++){
			  if(rlabel == label[i]){

				  flag = 1;
				  ss << rvalue;

				  switch(id[i]){
				  case 0:
					  ss >> mydouble;
					  *((double *)addr[i]) = mydouble;
					  break;
				  case 1:
					  ss >> myint;
					  *((int *)addr[i]) = myint;
					  break;
				  case 2:
					  ss >> mystring;
					  *((string *)addr[i]) = mystring;
					  break;
				  }

				  ss.clear();
				  ss.str(string());

				  id[i] = -1;
			  }
		  }
	  }

	  for(i = 0; i < n; i++){
		  if(id[i] >= 0 && addr[i] != &input_sim_file && addr[i] != &input_gal_file &&
				  addr[i] != &pw_alpha && addr[i] != &pw_beta && addr[i] != &pnfw_beta &&
				  addr[i] != &flag_switch_deflection_off){
			  ERROR_MESSAGE();
			  cout << "parameter " << label[i] << " needs to be set in the parameter file " << filename << endl;
			  exit(0);
		  }

		  /// DEFAULT VALUES
		  /// in case they are not set in the parameter file
		  if(id[i] >= 0 && addr[i] == &pw_alpha){
			  pw_alpha = 1./6.;
		  }
		  if(id[i] >= 0 && addr[i] == &pw_beta){
			  pw_beta = -1.0;
		  }
		  if(id[i] >= 0 && addr[i] == &pnfw_beta){
			  pnfw_beta = 2.0;
		  }
		  if(id[i] >= 0 && addr[i] == &flag_switch_deflection_off){
			  flag_switch_deflection_off = 0; //false, deflection is on
		  }
	  }

	  file_in.close();

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

	  // to compenstate for additional lens planes
	  if(flag_input_lens)
		  Nplanes++;

	  // convert to square degrees
	  fieldofview /= 3600. * 3600.;

	  printMultiLens();
}



void MultiLens::printMultiLens(){
	cout << endl << "outputfile " << outputfile << endl;

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
	int i, j, Ntot;
	unsigned long ind;
	double z1, z2;

	std::cout << "MultiLens::buildHaloTrees zsource = " << zsource << std::endl;

	if(!sim_input_flag){   /// If no input file is provided synthesize halos

		for(j=0,Ntot=0;j<Nplanes-1;j++){
			if(flag_input_lens && j == (flag_input_lens % Nplanes))
				continue;

			/*
			 * Setting the redshift range
			 * If there is a plane with an input lens on it, it is skipped over
			 * since it will not contain any halos
			 */
			if(j == 0) z1 = 0.0;
			//else z1 = QuickFindFromTable(Dl[j] - 0.5*dDl[j]);
			//else z1 = plane_redshifts[j] - 0.5*(plane_redshifts[j] - plane_redshifts[j-1]);
			else{
				locateD(coorDist_table-1,NTABLE,(Dl[j]-0.5*dDl[j]),&ind);
				z1 = redshift_table[ind];
			}

			if(flag_input_lens && j-1 == (flag_input_lens % Nplanes)){
				//z1 = plane_redshifts[j] - 0.5*(plane_redshifts[j] - plane_redshifts[j-2]);
				//z1 = QuickFindFromTable(Dl[j] - 0.5*(Dl[j] - Dl[j-2]));
				locateD(coorDist_table-1,NTABLE,(Dl[j] - 0.5*(Dl[j] - Dl[j-2])),&ind);
				z1 = redshift_table[ind];
			}

			if(j == Nplanes-2) z2 = zsource;
			//else z2 = QuickFindFromTable(Dl[j] + 0.5*dDl[j+1]);
			//else z2 = plane_redshifts[j] + 0.5*(plane_redshifts[j+1] - plane_redshifts[j]);
			else{
				locateD(coorDist_table-1,NTABLE,(Dl[j] + 0.5*dDl[j+1]),&ind);
				z2 = redshift_table[ind];
			}

			if(flag_input_lens && j+1 == (flag_input_lens % Nplanes)){
				//z2 = plane_redshifts[j] + 0.5*(plane_redshifts[j+2] - plane_redshifts[j]);
				//z2 = QuickFindFromTable(Dl[j] + 0.5*(Dl[j+2] - Dl[j]));
				locateD(coorDist_table-1,NTABLE,(Dl[j] + 0.5*(Dl[j+2] - Dl[j])),&ind);
				z2 = redshift_table[ind];
			}

			halo_data[j].reset(new HaloData(fieldofview,min_mass,mass_scale,z1,z2,mass_func_type,pw_alpha,cosmo,seed));

			Ntot+=halo_data[j]->Nhalos;
		}

	}else{
		// If input file is provided, read it in and put the halos onto lens planes.

		readInputSimFile(cosmo);
		unsigned long j1=0,j2=0;

		for(j=0,Ntot=0;j<Nplanes-1;j++){
			if(flag_input_lens && j == (flag_input_lens % Nplanes))
				continue;

			/*
			 * Setting the redshift range
			 * If there is a plane with an input lens on it, it is skipped over
			 * since it will not contain any halos
			 */
			if(j == 0) z1 = 0.0;
			//else z1 = QuickFindFromTable(Dl[j] - 0.5*dDl[j]);
			//else z1 = plane_redshifts[j] - 0.5*(plane_redshifts[j] - plane_redshifts[j-1]);
			else{
				locateD(coorDist_table-1,NTABLE,(Dl[j]-0.5*dDl[j]),&ind);
				z1 = redshift_table[ind];
			}

			if(flag_input_lens && j-1 == (flag_input_lens % Nplanes)){
				//z1 = plane_redshifts[j] - 0.5*(plane_redshifts[j] - plane_redshifts[j-2]);
				//z1 = QuickFindFromTable(Dl[j] - 0.5*(Dl[j] - Dl[j-2]));
				locateD(coorDist_table-1,NTABLE,(Dl[j] - 0.5*(Dl[j] - Dl[j-2])),&ind);
				z1 = redshift_table[ind];
			}

			if(j == Nplanes-2) z2 = zsource;
			//else z2 = QuickFindFromTable(Dl[j] + 0.5*dDl[j+1]);
			//else z2 = plane_redshifts[j] + 0.5*(plane_redshifts[j+1] - plane_redshifts[j]);
			else{
				locateD(coorDist_table-1,NTABLE,(Dl[j] + 0.5*dDl[j+1]),&ind);
				z2 = redshift_table[ind];
			}

			if(flag_input_lens && j+1 == (flag_input_lens % Nplanes)){
				//z2 = plane_redshifts[j] + 0.5*(plane_redshifts[j+2] - plane_redshifts[j]);
				//z2 = QuickFindFromTable(Dl[j] + 0.5*(Dl[j+2] - Dl[j]));
				locateD(coorDist_table-1,NTABLE,(Dl[j] + 0.5*(Dl[j+2] - Dl[j])),&ind);
				z2 = redshift_table[ind];
			}

			/// Find which halos are in redshift range
			locateD(halo_zs-1,Nhalos,z1,&j1);
			if(j1 == Nhalos) j1 = Nhalos-1;
			locateD(halo_zs-1,Nhalos,z2,&j2);
			if(j2 == Nhalos) j2 = Nhalos-1;

			/// Use other constructor to create halo data

			halo_data[j].reset(new HaloData(cosmo,&halos[j1],&halo_pos[j1],&halo_zs[j1],&halo_id[j1],j2-j1,z1,z2,mass_scale));

			//for(int i = 0; i<10 ;++i) cout << "Rmax:" << halos[j1+i].Rmax << "mass:" << halos[j1+i].mass << "rscale:" << halos[j1+i].rscale << "x = " << halo_pos[j1+i][0] << " " << halo_pos[j1+i][1] << endl;

			Ntot += halo_data[j]->Nhalos;
			///TODO: MARGARITA/BEN can be removed when the self-lensing problem is fixed 100%
			cout << "j: " << j << " z1: " << z1 << " j1: " << j1 << " z2: " << z2 << " j2: " << j2 << endl;
			//cout << "j: " << j << " z1: " << Dl[j] - 0.5*dDl[j] << " j1: " << j1 << " z2: " << Dl[j] + 0.5*dDl[j+1] << " j2: " << j2 << endl;

		}

	}

	for(j=0;j<Nplanes-1;j++){
		if(flag_input_lens && j == (flag_input_lens % Nplanes))
			continue;

		std::cout << "  Building tree on plane " << j << " number of halos: " << halo_data[j]->Nhalos << std::endl;

		switch(internal_profile){
		case PowerLaw:
			halo_tree[j].reset(new QuadTreePowerLaw(pw_beta,&halo_data[j]->pos[0],halo_data[j]->Nhalos
							,halo_data[j]->halos,halo_data[j]->kappa_background));
			break;
		case NFW:
			halo_tree[j].reset(new QuadTreeNFW(&halo_data[j]->pos[0],halo_data[j]->Nhalos
							,halo_data[j]->halos,halo_data[j]->kappa_background));
			break;
		case PseudoNFW:
			halo_tree[j].reset(new QuadTreePseudoNFW(pnfw_beta,&halo_data[j]->pos[0],halo_data[j]->Nhalos
							,halo_data[j]->halos,halo_data[j]->kappa_background));
			break;
		case NSIE:
			halo_tree[j].reset(new QuadTreeNSIE(&halo_data[j]->pos[0],halo_data[j]->Nhalos
							,halo_data[j]->halos,halo_data[j]->kappa_background));
			break;
		default:
			ERROR_MESSAGE();
			cout << "There is no such case for the halo trees. Please choose from:" << endl;
			cout << "0: PowerLaw, 1: NFW, 2: PseudoNFW, 3: NSIE" << endl;
			exit(1);
			break;
		}

	}

	cout << "constructed " << Ntot << " halos" << endl;

	/*
	stringstream f;
	f << "halos_" << zsource << ".data";
	string filename = f.str();
	ofstream file_area(filename.c_str());
	if(!file_area){
		cout << "unable to create file " << filename << endl;
		exit(1);
	}

	for(j = 0; j < Nplanes-1; j++){

		if(flag_analens && j == (flag_analens % Nplanes))
			continue;

		for(i = 0; i < halo_data[j]->Nhalos; i++){

			double fac = 180/pi/3600./Dl[j]*(1+plane_redshifts[j]);

			file_area << plane_redshifts[j] << " ";
			file_area << i << " " << halo_data[j]->halos[i].mass << " " << fac*halo_data[j]->halos[i].Rmax << " " << fac*halo_data[j]->halos[i].rscale << " ";
			file_area << fac*halo_data[j]->pos[i][0] << " " << fac*halo_data[j]->pos[i][1] << endl;

		}
	}


	file_area.close();
	*/
	/*
	for(int l=0; l < 100; l++){
	double ray[2];

	double fac = l/100.;

	ray[0] = fac*sqrt(fieldofview);
	ray[1] = ray[0];

	ray[0] = ray[0]*pi/180.;
	ray[1] = ray[1]*pi/180.;

	double xx[2]={0.0,0.0};

	int halos;

	for(j = 0, halos = 0; j < Nplanes-1; j++){

		if(flag_analens && j == (flag_analens % Nplanes))
			continue;

		xx[0] = ray[0]*Dl[j]/(1+plane_redshifts[j]);
		xx[1] = ray[1]*Dl[j]/(1+plane_redshifts[j]);

		for(i = 0; i < halo_data[j]->Nhalos; i++){
			double r2 = (halo_data[j]->pos[i][0] - xx[0])*(halo_data[j]->pos[i][0] - xx[0])
						+ (halo_data[j]->pos[i][1] - xx[1])*(halo_data[j]->pos[i][1] - xx[1]);

			if(r2 <= halo_data[j]->halos[i].Rmax*halo_data[j]->halos[i].Rmax)
				halos++;
		}

	}


	cout << "ray: x " << ray[0]*180/pi << " y " << ray[1]*180/pi << ", halos: " << halos/float(Ntot) << endl;
	}

	 */

	/*
	double xp,x,xo;
	double yp,y,yo;
	double alpha_x_one, alpha_x_two;
	double alpha_y_one, alpha_y_two;

	double fac_one = (Dl[Nplanes-1]-Dl[0])/Dl[Nplanes-1];
	double fac_two = (Dl[Nplanes-1]-Dl[1])/Dl[Nplanes-1];

	x=xo=0; y=yo=0;

	xp = x*Dl[0]/(1+redshift[0]);
	yp = y*Dl[0]/(1+redshift[0]);

	xp = haloModel[0]->pos[0][0]-xp;
	yp = haloModel[0]->pos[0][1]-yp;

	double r2 = xp*xp+yp*yp;

	double ratio = 0.5*r2/haloModel[0]->halos[0].Rmax/haloModel[0]->halos[0].Rmax;

	alpha_x_one = 4*Grav*haloModel[0]->halos[0].mass*xp*(exp(-ratio)-1)/r2;
	alpha_y_one = 4*Grav*haloModel[0]->halos[0].mass*yp*(exp(-ratio)-1)/r2;

	double fac = (Dl[1]-Dl[0])/Dl[1];

	// angle on the second plane
	x = xo - fac*alpha_x_one;
	y = yo - fac*alpha_y_one;

	xp = x*Dl[1]/(1+redshift[1]);
	yp = y*Dl[1]/(1+redshift[1]);

	xp = haloModel[1]->pos[0][0]-xp;
	yp = haloModel[1]->pos[0][1]-yp;

	r2 = xp*xp+yp*yp;

	ratio = 0.5*r2/haloModel[1]->halos[0].Rmax/haloModel[1]->halos[0].Rmax;

	alpha_x_two = 4*Grav*haloModel[1]->halos[0].mass*xp*(exp(-ratio)-1)/r2;
	alpha_y_two = 4*Grav*haloModel[1]->halos[0].mass*yp*(exp(-ratio)-1)/r2;

	x = xo - fac_one*alpha_x_one - fac_two*alpha_x_two;
	y = yo - fac_one*alpha_y_one - fac_two*alpha_y_two;

	cout << "0 " << alpha_x_one << " " << alpha_y_one << endl;
	cout << "1 " << alpha_x_two << " " << alpha_y_two << endl;
	cout << xo << " " << yo << " " << x << " " << y << endl << endl;
	*/

}

/**
 * Set the redshifts of the planes by mapping the correct
 * redshift by using the coordinate distance table
 */
void MultiLens::setRedshifts(){
	int i;
	unsigned long j;
	// assigns the redshifts and plugs in the input plane
	cout << "z: ";
	for(i=0; i<Nplanes; i++){
		locateD(coorDist_table-1,NTABLE,Dl[i],&j);
		plane_redshifts[i] = redshift_table[j];
		cout << plane_redshifts[i] << " ";
	}
	cout << endl;
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
		Np = Nplanes;
	else
		Np = Nplanes+1;

	double Ds = cosmo->coorDist(0,zsource);

	double Dlens;
	if(flag_input_lens) Dlens = cosmo->coorDist(0,input_lens->getZlens());
	else Dlens = Ds;

	/// spaces lD equally up to the source, including 0 and Ds
	/// therefore we need Nplanes+1 values
	/// however, if there is an input plane, we will need Nplanes values, since the input plane will take up a value itself
	fill_linear(lD,Np,0.,Ds);

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
	cout << endl << endl;

}

double MultiLens::getZlens(){
	if(flag_input_lens)
		return input_lens->getZlens();
	else
		return plane_redshifts[0];
}

void MultiLens::setZlens(double z){
	if(flag_input_lens)
		input_lens->setZlens(z);
	else{
		cout << "It is not possible to reset the redshift of the input AnaLens plane a MultiLens" << endl;
		ERROR_MESSAGE();
		exit(1);
	}
}

/// read in halos from a simulation file
void MultiLens::readInputSimFile(CosmoHndl cosmo){

	char c;
	double ra,dec,z,vmax,vdisp,r_halfmass;
	string strg;
	unsigned long i,j;
	unsigned long haloid,id,np;

	//int index;

	if(internal_profile != NSIE && internal_profile != PowerLaw){
		std::cout << "The internal profile of the halos while using simulation input files must be a Power Law."
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

	//file_in >> Nhalos;
	//Nhalos = 10;
	//cout << Nhalos << endl;

	//std::vector<HaloStructure> halo_vec;
	std::vector<HaloStructure> halo_vec;
	std::vector<double> halo_zs_vec;
	std::vector<double *> halo_pos_vec;
	std::vector<unsigned long> halo_id_vec;

	// read in data
	int j_max;
	double mass_max=0,R_max=0,V_max=0;
	double *theta;
	HaloStructure halo;

	//for(i=0,j=0 ; i<21152; ++i){
	for(i=0,j=0 ; c!='#'; ++i){
		// read a line of data
		//file_in >> haloid >>  id >>  ra >>  dec >>  z
			//	 >>  np >>  vdisp >>  vmax >>  r_halfmass;
		file_in >> c >> haloid >> c >> id >> c >> ra >> c >> dec >> c >> z
				 >> c >> np >> c >> vdisp >> c >> vmax >> c >> r_halfmass >> c;  //TODO the GalID will miss the first digit using this method.  No other method stops at the end of file.
		//std::cout << id << c << id << c << ra << c << dec << c << z
		//				 << c << np << c << vdisp << c << vmax << c << r_halfmass << std::endl;
		//cout << i << "  z: " << z << " np: " << np << " vmax:" << vmax << "  " << file_in.peek() << endl;


		if(np > 0.0 && vdisp > 0.0){

			halo_vec.push_back(halo);

			halo_vec[j].mass = np*8.6e8/cosmo->gethubble();
			halo_vec[j].Rmax = halo_vec[j].mass*Grav/2/pow(vdisp/lightspeed,2);  // SIS value
			if(internal_profile == NSIE){
				halo_vec[j].sigma = vdisp;
				halo_vec[j].fratio = (ran2(seed)+1)*0.5;  //TODO This is a kluge.
				halo_vec[j].pa = 2*pi*ran2(seed);  //TODO This is a kluge.
				  // TODO Needs to be changed.
				//if(halo_vec[j].mass > 1.0e14) halo_vec[j].rscale = halo_vec[j].Rmax*0.1*pow(halo_vec[j].mass/1.0e14,0.25);
				//else
				halo_vec[j].rscale = 0.0;
				halo_vec[j].Rsize = rmaxNSIE(halo_vec[j].sigma,halo_vec[j].mass,halo_vec[j].fratio,halo_vec[j].rscale);
				halo_vec[j].Rmax = MAX(1.0,1.0/halo_vec[j].fratio)*halo_vec[j].Rsize;  // redefine
			}else{
				// initialize unused variables to harmless values in PowerLaw case
				halo_vec[j].Rsize = halo_vec[j].Rmax;
				halo_vec[j].rscale = halo_vec[j].Rmax;
				halo_vec[j].fratio = 0;
				halo_vec[j].pa = 0;
				halo_vec[j].sigma = 0;
			}
			if(halo_vec[j].mass > mass_max) {
				mass_max = halo_vec[j].mass;
				j_max = j;
			}
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

			theta = new double[2];
			theta[0] = ra*pi/180.;
			theta[1] = dec*pi/180.;
			halo_pos_vec.push_back(theta);

			++j;
		}
	}
	file_in.close();
	std::cout << halo_vec.size() << " halos read in."<< std::endl
			<< "Max input mass = " << mass_max << " zmax " << j_max << "  R max = " << R_max << "  V max = " << V_max << std::endl;

	Nhalos = halo_vec.size();

	halos = new HaloStructure[Nhalos];
	halo_zs = new double[Nhalos];
	halo_id = new unsigned long[Nhalos];
	halo_pos = PosTypeMatrix(0,Nhalos-1,0,2);

	for(i=0;i<Nhalos;++i){
		halo_id[i] = halo_id_vec[i];
		halo_zs[i] = halo_zs_vec[i];
		halo_pos[i] = halo_pos_vec[i];
		halos[i] = halo_vec[i];
	}

	std::cout << "sourting in MultiLens::readInputSimFile()" << std::endl;
	// sort the halos by readshift
	MultiLens::quicksort(halos,halo_pos,halo_zs,halo_id,Nhalos);

	std::cout << "leaving MultiLens::readInputSimFile()" << std::endl;
}


/** Sets the internal parameters of the multiple lens model
 * first the redshifts of the planes are calculated
 * then the coordinate distances to the different planes
 * the planes are populated by halos and the halo trees are built
 * and finally the internal parameters of the input plane are set, in case there is one
 */
void MultiLens::setInternalParams(CosmoHndl cosmo, SourceHndl source){

	if(flag_input_lens)
		input_lens->setInternalParams(cosmo,source);

	if( (cosmo->getOmega_matter() + cosmo->getOmega_lambda()) != 1.0 ){
		printf("ERROR: MultiLens can only handle flat universes at present.  Must change cosmology.\n");
		exit(1);
	}

	/// makes the oordinate distance table for the calculation of the redshifts of the different planes
	if(table_set == false) {std::cout << "making tables" << std::endl; make_table(cosmo);}
	setCoorDist(cosmo);
	setRedshifts();

	buildHaloTrees(cosmo);
	std:: cout << " done " << std:: endl;
}

/** \brief Read in information from a Virgo Millennium Data Base http://gavo.mpa-garching.mpg.de/MyMillennium/
 *
 * query select * from MockLensing.dbo.m21_20_39_021_bc03_Ben_halos
 *
 * This is information on the dark matter halos only.  There are 13 entries in each line separated by commas.
 * The comments must be removed from the beginning of the data file and the total number of halos must be added
 * as the first line.
 */
/// Sort halos[] and brr[][] by content off arr[]
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
	swap(&arr[pivotindex],&arr[N-1]);
	//SwapPointsInArray(&pointarray[pivotindex],&pointarray[N-1]);
	swap(&halos[pivotindex],&halos[N-1]);
	swap(&brr[pivotindex],&brr[N-1]);
	swap(&id[pivotindex],&id[N-1]);
	newpivotindex=0;

	// partition list and array
	for(i=0;i<N;++i){
		if(arr[i] <= pivotvalue){
			swap(&arr[newpivotindex],&arr[i]);
			//SwapPointsInArray(&pointarray[newpivotindex],&pointarray[i]);
			swap(&halos[newpivotindex],&halos[i]);
			swap(&brr[newpivotindex],&brr[i]);
			swap(&id[newpivotindex],&id[i]);
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

		){
	unsigned long j;
	short out = 1;

	toggle_source_plane = true;

	if(z<=0.0){
		cout << "Warning: Source redshift cann't be set to " << z << " in MultiLens::ResetSourcePlane." << endl;
		return 0;
	}
/*
	if(z > plane_redshifts[Nplanes-1]){
		cout << "Warning: Implanted source is at higher redshift than simulation was constructed for." << endl
		<< "It is not being added. " << plane_redshifts[Nplanes-1] << " < " << z << endl;

		return;
	}
*/
	// j is the index of the next plane at higher redshift, This plane will be temporarily replaced and used as a source plane

	double Ds = cosmo->coorDist(0,z);

	//locateD(plane_redshifts-1,Nplanes,z,&j);
	locateD(Dl-1,Nplanes,Ds,&j);
	assert(j <= Nplanes && j >=0);
	if(j >= Nplanes-1){
		j--;
	}
	else if(j > 0){
		unsigned long ind;

		locateD(coorDist_table-1,NTABLE,(Dl[j]-0.5*dDl[j]),&ind);
		double z1 = redshift_table[ind];

		//if(nearest) j = ((z-plane_redshifts[j-1]) > (plane_redshifts[j]-z)) ? j : j-1;
		//if(nearest) j = ((Ds-Dl[j-1]) > (Dl[j]-Ds)) ? j : j-1;

		//if(nearest)
		j = (z>=z1) ? j : j-1;
	}

	if(nearest && (j < Nplanes-1) ){
		zs_implant = plane_redshifts[j];
		Ds_implant = Dl[j];
		if(j > 0) dDs_implant = dDl[j];
		else  dDs_implant = Ds_implant;
	}else{
		Ds_implant = Ds;
		zs_implant = z;
		if(j > 0) dDs_implant = cosmo->coorDist(plane_redshifts[j-1],z);
		else  dDs_implant = Ds;
	}

	//std::cout << "Source on plane " << j << " zs " << zs_implant << " Ds " << Ds << " dDs " << dDs_implant << std::endl;
	index_of_new_sourceplane = j;

	return out;
}


/**
 * \brief Change the implanted source in the Multilens.
 *
 * When rays are subsequently shot through the simulation the surface brightness of the source will
 * be added to the point->surface_brightness.  Sources can be implanted without altering the existing
 * lens or rays.  The rays need to be re-shot after the index is re-shot.
 *
 * This is meant for use when the internal AnaSource has already been initialized with multiple sources.
 *
void MultiLens::ImplantSource(
		unsigned long index        /// the index of the galaxy to be made the current galaxy
		,CosmoHndl cosmo           /// cosmology
		){

	if(!gal_input_flag){
		ERROR_MESSAGE();
		std::cout << "The AnaSource has not been constructed within MultiLens" << std::endl;
		exit(1);
	}
	unsigned long j;
	double z;

	anasource->setIndex(index);
	z = anasource->getZ();

	if(anasource->getZ() > plane_redshifts[Nplanes-1]){
		cout << "Warning: Implanted source is at higher redshift than simulation was constructed for." << endl
		<< "It is not being added." << endl;
		return;
	}

	Ds_implant = cosmo->angDist(0,z);
	zs_implant = z;


	//ys_implant[0] = Ds_implant*anasource->get_theta()[0];
	//ys_implant[1] = Ds_implant*anasource->get_theta()[1];

	ys_implant[0] = ys_implant[1] = 0.0;


	locateD(plane_redshifts-1,Nplanes,zs_implant,&j);

	dDs_implant = cosmo->coorDist(z,plane_redshifts[j]);

	index_of_new_sourceplane = j;
}
*/

void swap(double **a,double **b){
	double *tmp;
	tmp=*a;
	*a=*b;
	*b=tmp;
}
void swap(unsigned long *a,unsigned long *b){
	unsigned long tmp;
	tmp=*a;
	*a=*b;
	*b=tmp;
}


