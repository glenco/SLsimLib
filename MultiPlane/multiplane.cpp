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

using namespace std;

HaloData::HaloData(HaloStructure *halostrucs,PosType **positions,unsigned long Nhaloss):
	pos(positions), halos(halostrucs), Nhalos(Nhaloss)
{
	allocation_flag = false;
}

/**
 *  \brief In this constructor the halos are created from a mass function and given random positions on the lens plane
 */

HaloData::HaloData(
		double fov            /// field of view in square degrees
		,double min_mass      /// Minimum mass of a halo
		,double mass_scale      /// mass scale
		,double z1            /// lowest redshift
		,double z2            /// highest redshift
		,int mass_func_type   /// mass function type: 0 Press-Schechter, 1 Sheth-Tormen
		,CosmoHndl cosmo     /// cosmology
		,long *seed
	)
	{

	allocation_flag = true;

    HALO ha(cosmo,min_mass,0.);

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

	Nhalosbin[0] = cosmo->haloNumberDensityOnSky(pow(10,Logm[0])/cosmo->gethubble(),z1,z2,mass_func_type)*fov;

	Nhaloestot = Nhalosbin[0];
	Nhalosbin[0] = 1;
	int k;

	for(k=1;k<Nmassbin;k++){
		// cumulative number density in one square degree
		Nhalosbin[k] = cosmo->haloNumberDensityOnSky(pow(10,Logm[k])/cosmo->gethubble(),z1,z2,mass_func_type)*fov;
		// normalize the cumulative distribution to one
		Nhalosbin[k] = Nhalosbin[k]/Nhaloestot;
	}

	Nhalos = (long)(poidev(float(Nhaloestot), seed) + 0.5);

	halos = new HaloStructure[Nhalos];
	pos = PosTypeMatrix(0,Nhalos-1,0,2);
	double rr,theta,maxr,zi;
	unsigned long i;
	for(i = 0, kappa = 0.0; i < Nhalos; i++){
		//maxr = pi*sqrt(fov*0.5/pi)/180*Dli[i]; // fov is a circle
		zi = z1+(z2-z1)*ran2 (seed);
		Dli.push_back(cosmo->angDist(0,zi));
		maxr = pi*sqrt(fov/pi)/180*Dli[i]; // fov is a circle
		rr = maxr*sqrt(ran2(seed));
		theta = 2*pi*ran2(seed);
		/*do{
			x = 2*maxr*ran2(seed) - maxr;
			y = 2*maxr*ran2(seed) - maxr;
		}while((x*x+y*y) > (maxr*maxr));*/
		pos[i][0] = rr*cos(theta);
		pos[i][1] = rr*sin(theta);

		halos[i].mass = pow(10,InterpolateYvec(Nhalosbin,Logm,ran2 (seed)));
		ha.reset(halos[i].mass,zi);
		halos[i].mass /= mass_scale;
		//halos[i].mass = vmasses[i];
		//halos[i].Rmax = vsizes[i];
		halos[i].Rmax = ha.getRvir()*cosmo->gethubble();
		//halos[i].rscale = vsizes[i]/vscale[i]; // get the Rscale=Rmax/c
		halos[i].rscale = halos[i].Rmax/ha.getConcentration(0); // get the Rscale=Rmax/c
		pos[i][2] = 0.0;//halos[i].Rmax;


		/* calculate the mass density on the plane */
		kappa += halos[i].mass;
	}

}

HaloData::~HaloData(){

	if(allocation_flag){
		free_PosTypeMatrix(pos,0,Nhalos-1,0,2);
		delete[] halos;
	}
}


MultiLens::MultiLens(string filename,long *my_seed) : Lens(){
	readParamfile(filename);

	// flag to determine if halos are created randomly or read in from a external simulation.
	if(input_sim_file.size() < 1) sim_input_flag = false;
	else sim_input_flag = true;

	plane_redshifts = new double[Nplanes];
	Dl = new double[Nplanes];
	dDl = new double[Nplanes];

	charge = 4*pi*Grav*mass_scale;

	//halo_tree = new ForceTreeHndl[Nplanes-1];
	halo_tree = new auto_ptr<ForceTree>[Nplanes-1];

	//halodata = new HaloDataHndl[Nplanes-1];
	halodata = new auto_ptr<HaloData>[Nplanes-1];

	if(flag_analens){
		analens = new AnaLens(filename);
	}

	seed = my_seed;
}

MultiLens::~MultiLens(){

	delete[] halo_tree;
	for(int j=0;j<Nplanes-1;j++) delete &(halo_tree[j]);
	delete[] Dl;
	delete[] plane_redshifts;
	delete[] dDl;
	for(int j=0;j<Nplanes-1;j++) delete &(halodata[j]);
	delete[] halodata;

	if(sim_input_flag){  // Otherwise these deallocations are done in the destructor of halodata's
		delete[] halos;
		delete[] halo_zs;
		free_PosTypeMatrix(halo_pos,0,Nhalos-1,0,2);
	}


	if(flag_analens)
		delete analens;
}

void MultiLens::readParamfile(string filename){
      const int MAXPARAM = 11;
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

	  addr[n] = &flag_analens;
	  id[n] = 1;
	  label[n++] = "flag_analens";

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

	  addr[n] = &MOKA_input_file;
	  id[n] = 2;
	  label[n++] = "MOKA_input_file";

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
		  if(id[i] >= 0 && addr[i] != &input_sim_file && addr[i] != &MOKA_input_file){
			  ERROR_MESSAGE();
			  cout << "parameter " << label[i] << " needs to be set!" << endl;
			  exit(0);
		  }
	  }

	  file_in.close();

	  // to compensate for the lst plane, which is the source plane
	  Nplanes++;

	  // to compenstate for analytic lens plane
	  if(flag_analens)
		  Nplanes++;

	  printMultiLens();

	  // convert to square degrees
	  fieldofview /= 3600. * 3600.;
}



void MultiLens::printMultiLens(){
	cout << endl << "outputfile " << outputfile << endl;

	cout << endl << "**multi lens model**" << endl;

	cout << "Nplanes " << Nplanes << endl;

	cout << "mass scale " << mass_scale << endl;

	cout << "min mass " << min_mass << endl;

	cout << "flag analens " << flag_analens << endl;

	cout << "field of view " << fieldofview << endl;

	cout << "internal profile type " << internal_profile << endl;
	switch(internal_profile){
	case 0:
		cout << "  Gaussian internal profile " << endl;
		break;
	case 1:
		cout << "  Power law internal profile " << endl;
		break;
	case 2:
		cout << "  NFW internal profile " << endl;
		break;
	case 3:
		cout << "  Pseudo NFW internal profile " << endl;
		break;
	}

	cout << "mass function type " << mass_func_type << endl << endl;
}

/// Populates the lensing planes with halos and builds force trees
void MultiLens::buildHaloTrees(
		CosmoHndl cosmo /// the cosmology
		,double zsource /// the source resdhift
		){
	int i, j, Ntot;
	double z1, z2;

	if(!sim_input_flag){   /// If no input file is provided synthesize halos

		for(j=0,Ntot=0;j<Nplanes-1;j++){
			if(flag_analens && j == (flag_analens % Nplanes))
				continue;

			// redshift range
			if(j == 0) z1 = 0.0;
			else z1 = plane_redshifts[j] - 0.5*(plane_redshifts[j] - plane_redshifts[j-1]);

			if(j-1 == (flag_analens % Nplanes)) z1 = plane_redshifts[j] - 0.5*(plane_redshifts[j] - plane_redshifts[j-2]);

			if(j == Nplanes-2) z2 = zsource;
			else z2 = plane_redshifts[j] + 0.5*(plane_redshifts[j+1] - plane_redshifts[j]);

			if(j+1 == (flag_analens % Nplanes)) z2 = plane_redshifts[j] + 0.5*(plane_redshifts[j+2] - plane_redshifts[j]);

			//halodata[j] = new HaloData(fieldofview,min_mass,z1,z2,mass_func_type,cosmo,seed);
			halodata[j] = auto_ptr<HaloData>(new HaloData(fieldofview,min_mass,mass_scale,z1,z2,mass_func_type,cosmo,seed));

			Ntot+=halodata[j]->Nhalos;
		}

	}else{
		// If input file is provided, read it in and put the halos onto lens planes.

		readInputSimFile(cosmo);
		unsigned long j1,j2;

		for(j=0,Ntot=0;j<Nplanes-1;j++){
			if(flag_analens && j == (flag_analens % Nplanes))
				continue;

			// redshift range
			if(j == 0) z1 = 0.0;
			else z1 = plane_redshifts[j] - 0.5*(plane_redshifts[j] - plane_redshifts[j-1]);

			if(j-1 == (flag_analens % Nplanes)) z1 = plane_redshifts[j] - 0.5*(plane_redshifts[j] - plane_redshifts[j-2]);

			if(j == Nplanes-2) z2 = zsource;
			else z2 = plane_redshifts[j] + 0.5*(plane_redshifts[j+1] - plane_redshifts[j]);

			if(j+1 == (flag_analens % Nplanes)) z2 = plane_redshifts[j] + 0.5*(plane_redshifts[j+2] - plane_redshifts[j]);

			/// Find which halos are in redshift range

			locateD(halo_zs-1,Nhalos,z1,&j1);
			if(j1 == Nhalos) j1 = Nhalos-1;
			locateD(halo_zs-1,Nhalos,z2,&j2);
			if(j2 == Nhalos) j2 = Nhalos-1;

			/// Use other constructor to create halo data

			//halodata[j] = new HaloData(&halos[j1],&halo_pos[j1],j2-j1);
			halodata[j] = auto_ptr<HaloData>(new HaloData(&halos[j1],&halo_pos[j1],j2-j1));

			//for(int i = 0; i<10 ;++i) cout << "Rmax:" << halos[j1+i].Rmax << "mass:" << halos[j1+i].mass << "rscale:" << halos[j1+i].rscale << "x = " << halo_pos[j1+i][0] << " " << halo_pos[j1+i][1] << endl;

			Ntot += halodata[j]->Nhalos;
		}

	}

	// *** these are not freed if this routine in used repeatedly
	for(j=0;j<Nplanes-1;j++){
		if(flag_analens && j == (flag_analens % Nplanes))
			continue;

		switch(internal_profile){
		case 1:
			//halo_tree[j] = new ForceTreePowerLaw(1.9,&halodata[j]->pos[0],halodata[j]->Nhalos,halodata[j]->halos);
			halo_tree[j] = auto_ptr<ForceTree>(new ForceTreePowerLaw(1.9,&halodata[j]->pos[0],halodata[j]->Nhalos,halodata[j]->halos));
			break;
		case 2:
			//halo_tree[j] = new ForceTreeNFW(&halodata[j]->pos[0],halodata[j]->Nhalos,halodata[j]->halos);
			halo_tree[j] = auto_ptr<ForceTree>(new ForceTreeNFW(&halodata[j]->pos[0],halodata[j]->Nhalos,halodata[j]->halos));
			break;
		case 3:
			//halo_tree[j] = new ForceTreePseudoNFW(1.9,&halodata[j]->pos[0],halodata[j]->Nhalos,halodata[j]->halos);
			//halo_tree[j] = auto_ptr<ForceTree>(new ForceTreePseudoNFW(2.0,&halodata[j]->pos[0],halodata[j]->Nhalos,halodata[j]->halos,true,5,3,false,0.1));
			halo_tree[j] = auto_ptr<ForceTree>(new ForceTreePseudoNFW(2.0,&halodata[j]->pos[0],halodata[j]->Nhalos,halodata[j]->halos));
			break;
		case 0:
			//halo_tree[j] = new ForceTreeGauss(&halodata[j]->pos[0],halodata[j]->Nhalos,halodata[j]->halos);
			halo_tree[j] = auto_ptr<ForceTree>(new ForceTreeGauss(&halodata[j]->pos[0],halodata[j]->Nhalos,halodata[j]->halos));
			break;
		}

		/* to obtain 1/physical_distance^2
		 * to be compatible with the rayshooter*/
		double area = fieldofview * pow(pi/180*Dl[j]/(1+plane_redshifts[j]),2.0);
		halodata[j]->kappa /= area;

		cout << halodata[j]->kappa << endl;

		halo_tree[j]->kappa = halodata[j]->kappa;
	}


	cout << "constructed " << Ntot << " halos" << endl;

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

		for(i = 0; i < halodata[j]->Nhalos; i++){

			double fac = 180/pi/3600./Dl[j]*(1+plane_redshifts[j]);

			file_area << plane_redshifts[j] << " ";
			file_area << i << " " << halodata[j]->halos[i].mass << " " << fac*halodata[j]->halos[i].Rmax << " " << fac*halodata[j]->halos[i].rscale << " ";
			file_area << fac*halodata[j]->pos[i][0] << " " << fac*halodata[j]->pos[i][1] << endl;

		}
	}


	file_area.close();

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

		for(i = 0; i < halodata[j]->Nhalos; i++){
			double r2 = (halodata[j]->pos[i][0] - xx[0])*(halodata[j]->pos[i][0] - xx[0])
						+ (halodata[j]->pos[i][1] - xx[1])*(halodata[j]->pos[i][1] - xx[1]);

			if(r2 <= halodata[j]->halos[i].Rmax*halodata[j]->halos[i].Rmax)
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

void MultiLens::setRedshift(double zsource){
	std:: vector<double> lz;
	int Np;
	if(flag_analens)
		Np = Nplanes;
	else
		Np = Nplanes+1;

	fill_linear(lz,Np,0.,zsource);

	int j=0, flag=0;
	if(flag_analens && analens->zlens < lz[1]){
		plane_redshifts[j] = analens->zlens;
		flag_analens = Nplanes;
		flag = 1;
		j++;
	}
	for(int i=1; i<Np; i++){
		plane_redshifts[j] = lz[i];

		if(flag_analens && flag == 0)
			if(analens->zlens > lz[i] && analens->zlens <= lz[i+1]){
				plane_redshifts[j] = lz[i];
				plane_redshifts[++j] = analens->zlens;
				flag_analens = j;
				flag = 1;
			}
		j++;
	}

	if(flag_analens)
		cout << "zlens " << analens->zlens << " on plane number " << (flag_analens % Nplanes) << endl;

	cout << "z: ";
	for(int i = 0; i < Nplanes; i++)
		cout << plane_redshifts[i] << " ";
	cout << endl;

}

double MultiLens::getZlens(){
	if(flag_analens)
		return analens->zlens;
	else
		return plane_redshifts[0];
}

void MultiLens::setZlens(double z){
	if(flag_analens)
		analens->zlens = z;
}

// sets the redshifts and distances for the lens planes
void MultiLens::setInternalParams(CosmoHndl cosmo, double zsource){
	int j;

	if( (cosmo->getOmega_matter() + cosmo->getOmega_lambda()) != 1.0 ){
		printf("ERROR: MultiLens can only handle flat universes at present.  Must change cosmology.\n");
		exit(1);
	}

	setRedshift(zsource);

	Dl[0] = cosmo->coorDist(0,plane_redshifts[0]);
	dDl[0] = Dl[0];  // distance between jth plane and the previous plane
	for(j = 1; j < Nplanes; j++){

		Dl[j] = cosmo->coorDist(0,plane_redshifts[j]);
		dDl[j] = Dl[j] - Dl[j-1]; // distance between jth plane and the previous plane
	}

	cout << "Dl: ";
	for(j = 0; j < Nplanes; j++)
		cout << Dl[j] << " ";
	cout << endl;

	cout << "dDl: ";
	for(j = 0; j < Nplanes; j++)
		cout << dDl[j] << " ";
	cout << endl << endl;

	buildHaloTrees(cosmo,zsource);

	if(flag_analens)
		analens->setInternalParams(cosmo,zsource);
}

/** \brief Read in information from a Virgo Millennium Data Base http://gavo.mpa-garching.mpg.de/MyMillennium/
 *
 * query select * from MockLensing.dbo.m21_20_39_021_bc03_Ben_halos
 *
 * This is information on the dark matter halos only.  There are 13 entries in each line separated by commas.
 * The comments must be removed from the beginning of the data file and the total number of halos must be added
 * as the first line.
 */
void MultiLens::readInputSimFile(CosmoHndl cosmo){

	char c;
	int type;
	double ra,dec,z,zob,mass,massct,vmax,vdisp,r_halfmass;
	string strg;
	unsigned long i,j;
	unsigned long id,np;

	//int index;

	ifstream file_in(input_sim_file.c_str());
	if(!file_in){
		cout << "Can't open file " << input_sim_file << endl;
		exit(1);
	}


	file_in >> Nhalos;

	cout << Nhalos << endl;

	halos = new HaloStructure[Nhalos];
	halo_zs = new double[Nhalos];
	halo_pos = PosTypeMatrix(0,Nhalos-1,0,2);

	// read in data
	for(i=0,j=0 ; i < Nhalos && !file_in.eof() ; ++i){

		// rerad a line of data
		file_in >> id >> c >> id >> c >> type >> c >> ra >> c >> dec >> c >> z >> c >> zob
				 >> c >> np >> c >> mass >> c >> massct >> c >> vmax >> c >> vdisp >> c >> r_halfmass;
		//cout << id << c << id << c << type << c << ra << c << dec << c << z << c << zob
		  //		 << c << np << c << r200 << c << mass << c << vmax << c << vdisp << c << r_halfmass << endl;
		cout << "z:" << z << " np:" << mass*1.0e10/np << " mass:" << mass*1.0e10 << " vmax:" << vmax << endl;

		if(mass > 0.0){
			halos[j].mass = mass*1.0e10*cosmo->gethubble();
			halos[j].Rmax = cosmo->R200(z,mass*1.0e10*cosmo->gethubble());
			assert(halos[j].Rmax > 0.0);
			cout << "Rmax:" << halos[j].Rmax << endl;
			halos[j].rscale = halos[j].Rmax/cosmo->NFW_Concentration(vmax,halos[j].mass,halos[j].Rmax);
			halo_zs[j] = z;

			halos[j].mass /= mass_scale;

			halo_pos[j][0] = ra;
			halo_pos[j][1] = dec;
			++j;
		}
	}

	Nhalos = j;  // There is some waisted memory here which would have contained the halos with zero mass.

	// sort the halos by readshift
	MultiLens::quicksort(halos,halo_pos,halo_zs,Nhalos);

}

/// Sort halos[] and brr[][] by content off arr[]
void MultiLens::quicksort(HaloStructure *halos,double **brr,double *arr,unsigned long N){
	double pivotvalue;
	unsigned long pivotindex,newpivotindex,i;
	void swap(double *a,double *b);
	void swap(double **a,double **b);
	//void swapLong(unsigned long *a,unsigned long *b);

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
	newpivotindex=0;

	// partition list and array
	for(i=0;i<N;++i){
		if(arr[i] <= pivotvalue){
			swap(&arr[newpivotindex],&arr[i]);
			//SwapPointsInArray(&pointarray[newpivotindex],&pointarray[i]);
			swap(&halos[newpivotindex],&halos[i]);
			swap(&brr[newpivotindex],&brr[i]);
			++newpivotindex;
		}
	}
	--newpivotindex;

	quicksort(halos,brr,arr,newpivotindex);
	quicksort(&halos[newpivotindex+1],&brr[newpivotindex+1],&arr[newpivotindex+1],N-newpivotindex-1);

	return ;
}

void swap(double **a,double **b){
	double *tmp;
	tmp=*a;
	*a=*b;
	*b=tmp;
}
