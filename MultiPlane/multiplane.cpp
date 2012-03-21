/*
 * halo_model.cpp
 *
 *  Created on: Jan 14, 2012
 *      Author: mpetkova
 */

#include <slsimlib.h>
#include <sstream>
#include <utilities.h>

haloM::haloM(int jplane /// the index of the plane that holds the halos
		,double zsource  /// source redshift
		,CosmoHndl cosmo     /// cosmology
		,MultiLens *lens     /// lens
		){

	double fov = lens->fieldofview;
	int mfty = lens->mass_func_type;

    long seed = 2203;

    HALO ha(cosmo,lens->min_mass,0.);

    int iterator;
    std::vector<double> vmasses,vsizes,vscale;
    std::vector<int> vindex;
    Logm.resize(Nmassbin);
    Nhalosbin.resize(Nmassbin);

	std:: vector<double> Dli;
	/* fill the log(mass) vector */
	fill_linear(Logm,Nmassbin,lens->min_mass,MaxLogm);

	int Nplanes = lens->getNplanes();

	double Nhaloestot;
	float Nhaloestotf;

	double z1, z2;

	if(jplane == 0) z1 = 0.0;
	else z1 = lens->redshift[jplane] - 0.5*(lens->redshift[jplane] - lens->redshift[jplane-1]);

	if(jplane == Nplanes-2) z2 = zsource;
	else z2 = lens->redshift[jplane] + 0.5*(lens->redshift[jplane+1] - lens->redshift[jplane]);

	Nhalosbin[0]=cosmo->haloNumberDensityOnSky(pow(10,Logm[0]),z1,z2,mfty)*fov;
	Nhaloestot = Nhalosbin[0];
	int k;

	for(k=1;k<Nmassbin;k++){
		// cumulative number density in one square degree
		Nhalosbin[k]=cosmo->haloNumberDensityOnSky(pow(10,Logm[k]),z1,z2,mfty)*fov;
		// normalize the cumulative distribution to one
		Nhalosbin[k] = Nhalosbin[k]/Nhaloestot;
	}

	Nhaloestotf=poidev(float(Nhaloestot), &seed);
	Nhalos = Nhaloestotf;

	for(int k=0;k<Nhaloestotf;k++){
		iterator++;
		double ni = ran2 (&seed);
		// compute the mass inverting the cumulative distribution
		double logmi = getY(Nhalosbin,Logm,ni);
		double mi = pow(10.,logmi);
		vmasses.push_back(mi);
		vindex.push_back(iterator);
		double zi = z1+(z2-z1)*ran2 (&seed);
		ha.reset(mi,zi);
		double Rvir = ha.getRvir();
		vsizes.push_back(Rvir);
		double scale = ha.getConcentration(0);
		vscale.push_back(scale);
		Dli.push_back(cosmo->angDist(0,zi));
	}

	halos = new HaloStructure[Nhalos];
	pos = PosTypeMatrix(0,Nhalos-1,0,2);
	double x,y;
	for(int i = 0; i < Nhalos; i++){
		double maxr = pi*sqrt(fov*0.5/pi)/180*Dli[i]; // fov is a circle
		do{
			x = 2*maxr*ran2(&seed) - maxr;
			y = 2*maxr*ran2(&seed) - maxr;
		}while((x*x+y*y) > (maxr*maxr));
		pos[i][0] = x;
		pos[i][1] = y;
		pos[i][2] = 0.0;

		halos[i].mass = vmasses[i];
		halos[i].Rmax = vsizes[i];
		halos[i].rscale = vsizes[i]/vscale[i]; // get the Rscale=Rmax/c
	}

	/*
	stringstream f;
	f << "halos_" << jplane << ".data";
	string filename = f.str();
	ofstream file_area(filename.c_str());
	if(!file_area){
		cout << "unable to create file " << filename << endl;
		exit(1);
	}

	for(int i = 0; i < Nhalos; i++){
		file_area << i << " " << halos[i].mass << " " << halos[i].Rmax << " ";
		file_area << pos[i][0] << " " << pos[i][1] << endl;
	}
	file_area.close();
	*/

}

haloM::~haloM(){
	free_PosTypeMatrix(pos,0,Nhalos-1,0,2);
	delete[] halos;
}


MultiLens::MultiLens(string filename) : Lens(){
	readParamfile(filename);

	redshift = new double[Nplanes];
	Dl = new double[Nplanes];
	dDl = new double[Nplanes];

	charge = 4*pi*Grav*mass_scale;

	halo_tree = new ForceTreeHndl[Nplanes-1];
	haloModel = new haloMHndl[Nplanes-1];

	if(flag_analens){
		analens = new AnaLens(filename);
	}

}

MultiLens::~MultiLens(){
	delete[] halo_tree;
	delete[] Dl;
	delete[] redshift;
	delete[] dDl;
	delete[] haloModel;

	if(flag_analens)
		delete analens;
}

void MultiLens::readParamfile(string filename){
      const int MAXPARAM = 10;
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
		  if(id[i] >= 0){
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

	cout << "mass function type " << mass_func_type << endl << endl;
}

void MultiLens::buildHaloTree(CosmoHndl cosmo /// the cosmology
		,double zsource /// the source resdhift
		){
	int j, Ntot;

	for(j=0,Ntot=0;j<Nplanes-1;j++){
		if(flag_analens && j == (flag_analens % Nplanes))
			continue;
		haloModel[j] = new haloM(j,zsource,cosmo,this);

		switch(internal_profile){
		case 1:
			halo_tree[j] = new ForceTreePowerLaw(1.9,&haloModel[j]->pos[0],haloModel[j]->Nhalos,haloModel[j]->halos);
			break;
		case 2:
			halo_tree[j] = new ForceTreeNFW(&haloModel[j]->pos[0],haloModel[j]->Nhalos,haloModel[j]->halos);
			break;
		case 3:
			halo_tree[j] = new ForceTreePseudoNFW(1.9,&haloModel[j]->pos[0],haloModel[j]->Nhalos,haloModel[j]->halos);
			break;
		case 0:
			halo_tree[j] = new ForceTreeGauss(&haloModel[j]->pos[0],haloModel[j]->Nhalos,haloModel[j]->halos);
			break;
		}

		Ntot+=haloModel[j]->Nhalos;
	}

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

	cout << "constructed " << Ntot << " halos" << endl;
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
		redshift[j] = analens->zlens;
		flag_analens = Nplanes;
		flag = 1;
		j++;
	}
	for(int i=1; i<Np; i++){
		redshift[j] = lz[i];

		if(flag_analens && flag == 0)
			if(analens->zlens > lz[i] && analens->zlens <= lz[i+1]){
				redshift[j] = lz[i];
				redshift[++j] = analens->zlens;
				flag_analens = j;
				flag = 1;
			}
		j++;
	}

	/*
	cout << "z: ";
	for(int i = 0; i < Nplanes; i++)
		cout << redshift[i] << " ";
	cout << endl;
	*/
}

double MultiLens::getZlens(){
	if(flag_analens)
		return analens->zlens;
	else
		return redshift[0];
}

void MultiLens::setZlens(double z){
	if(flag_analens)
		analens->zlens = z;
	else{
		cout << "ERROR in setZlens() -- there is no analytic lens!" << endl;
		exit(1);
	}
}

void MultiLens::setInternalParams(CosmoHndl cosmo, double zsource){
	int j;

	if( (cosmo->getOmega_matter() + cosmo->getOmega_lambda()) != 1.0 ){
		printf("ERROR: MultiLens can only handle flat universes at present.  Must change cosmology.\n");
		exit(1);
	}

	setRedshift(zsource);

	Dl[0] = cosmo->coorDist(0,redshift[0]);
	dDl[0] = Dl[0];  // distance between jth plane and the previous plane
	for(j = 1; j < Nplanes; j++){

		Dl[j] = cosmo->coorDist(0,redshift[j]);
		dDl[j] = Dl[j] - Dl[j-1]; // distance between jth plane and the previous plane
	}

	/*
	cout << "Dl: ";
	for(j = 0; j < Nplanes; j++)
		cout << Dl[j] << " ";
	cout << endl;

	cout << "dDl: ";
	for(j = 0; j < Nplanes; j++)
		cout << dDl[j] << " ";
	cout << endl << endl;
	*/

	buildHaloTree(cosmo,zsource);

	if(flag_analens)
		analens->setInternalParams(cosmo,zsource);
}
