/*
 * halo_model.cpp
 *
 *  Created on: Jan 14, 2012
 *      Author: mpetkova
 */

#include <slsimlib.h>
#include <sstream>
#include <utilities.h>

haloM::haloM(double zsource  /// source redshift
		,CosmoHndl cosmo     /// cosmology
		,MultiLens *lens     /// lens
		,double fieldofview  /// field of view in square degree
		,int mfty            /// type of mass function PS (0) and ST (1) default is ST
		){

    long seed = 2203;

    HALO ha(cosmo,lens->min_mass,0.);

    int iterator;
    std::vector<double> vmasses,vsizes,vredshifts;
    std::vector<int> vindex;
    Logm.resize(Nmassbin);
    Nhaloes.resize(Nmassbin);

	std:: vector<double> Dli;
	/* fill the log(mass) vector */
	fill_linear(Logm,Nmassbin,lens->min_mass,MaxLogm);

	int Nplanes = lens->getNplanes();

	for(int i=0;i<Nplanes-1;i++){

		if(lens->flag_analens && i==lens->flag_analens)
			continue;

		double Nhaloestot;
		float Nhaloestotf;

		double z1, z2;

		if(i == 0) z1 = 0.0;
		else z1 = lens->redshift[i] - 0.5*(lens->redshift[i] - lens->redshift[i-1]);

		if(i == Nplanes-2) z2 = zsource;
		else z2 = lens->redshift[i] + 0.5*(lens->redshift[i+1] - lens->redshift[i]);

		Nhaloes[0]=cosmo->haloNumberDensityOnSky(pow(10,Logm[0]),z1,z2,mfty)*fieldofview;
		Nhaloestot = Nhaloes[0];
		int k;

		for(k=1;k<Nmassbin;k++){
			//if(omp_get_thread_num() == 0) cout << "k " << k << endl;
			// cumulative number density in one square degree
			Nhaloes[k]=cosmo->haloNumberDensityOnSky(pow(10,Logm[k]),z1,z2,mfty)*fieldofview;
			// normalize the cumulative distribution to one
			Nhaloes[k] = Nhaloes[k]/Nhaloestot;
		}

		Nhaloestotf=poidev(float(Nhaloestot), &seed);
		lens->NhalosinPlane[i] = Nhaloestotf;
		//cout << "Plane " << i << " z1 = " << z1 << " z2 = " << z2 << " Nhalos in plane = " << Nhaloestotf << endl;

		/*stringstream filename;
		filename << "output_" << i << ".dat";
		string f = filename.str();
		ofstream file_test(f.c_str());*/

		for(int k=0;k<Nhaloestotf;k++){
			iterator++;
			double ni = ran2 (&seed);
			// compute the mass inverting the cumulative distribution
		    double logmi = getY(Nhaloes,Logm,ni);
		    double mi = pow(10.,logmi);
		    vmasses.push_back(mi);
		    vindex.push_back(iterator);
		    double zi = z1+(z2-z1)*ran2 (&seed);
		    vredshifts.push_back(zi);
		    ha.reset(mi,zi);
		    double Rvir = ha.getRvir();
		    vsizes.push_back(Rvir);
		    Dli.push_back(cosmo->angDist(0,zi));

			//file_test << mi << endl;
		}

		//file_test.close();
	}
	N = vmasses.size();

	masses = new float[N];
	sizes = new float[N];
	redshifts = new float[N];

	pos = PosTypeMatrix(0,N-1,0,2);
	for(int i = 0; i < N; i++){
		double maxr = sqrt(fieldofview/M_PI)*Dli[i]*M_PI/180;
		double r = maxr*ran2(&seed);
		double theta=2*pi*ran2(&seed);

		pos[i][0] = r*cos(theta);
		pos[i][1] = r*sin(theta);
		pos[i][2] = 0.0;

		masses[i] = vmasses[i];
		sizes[i] = vsizes[i];
		redshifts[i] = vredshifts[i];
	}

	vmasses.clear();
	vsizes.clear();
	vredshifts.clear();
}

haloM::~haloM(){
	free_PosTypeMatrix(pos,0,N-1,0,2);
	delete[] masses;
	delete[] sizes;
	delete[] redshifts;
}


MultiLens::MultiLens(string filename) : Lens(){
	readParamfile(filename);

	redshift = new double[Nplanes];
	Dl = new double[Nplanes];
	dDl = new double[Nplanes];

	charge = 4*Grav*mass_scale;

	halo_tree = new ForceTreeHndl[Nplanes-1];
	NhalosinPlane = new IndexType[Nplanes-1];

	if(flag_analens){
		analens = new AnaLens(filename);
	}

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

	cout << "flag analens " << flag_analens << endl << endl;
}

MultiLens::~MultiLens(){
	delete[] halo_tree;
	delete[] NhalosinPlane;
	delete[] Dl;
	delete[] redshift;
	delete[] dDl;
	delete haloModel;

	if(flag_analens)
		delete analens;
}

void MultiLens::buildHaloTree(CosmoHndl cosmo /// the cosmology
		,double zsource /// the source resdhift
		,double fieldofview // the field of view in square degrees
		){
	IndexType N, N_last;

    haloModel = new haloM(zsource,cosmo,this,fieldofview,1);

	for(int j = 0, N_last = 0; j < Nplanes-1; j++){
		if(flag_analens && j==flag_analens)
			continue;

		N = NhalosinPlane[j];
		halo_tree[j] = new ForceTree(&haloModel->pos[N_last + N],N,&haloModel->masses[N_last + N],&haloModel->sizes[N_last + N],true,true,5,2,true,0.1);

		N_last = N;
	}

}

void MultiLens::setRedshift(double zsource){
	std:: vector<double> lz;
	/* fill the redshift vector logarithmically */
	int Np;
	if(flag_analens)
		Np = Nplanes;
	else
		Np = Nplanes+1;

	fill_linear(lz,Np,0.,zsource);

	int j=0, flag=0;
	for(int i = 1; i < Np; i++){
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

	cout << "z: ";
	for(int i = 0; i < Nplanes; i++)
		cout << redshift[i] << " ";
	cout << endl;
}

double MultiLens::getZlens(){
	return redshift[0];
}

void MultiLens::setZlens(double z){

}

void MultiLens::setInternalParams(CosmoHndl cosmo, double zsource){
	int j;

	if( (cosmo->getOmega_matter() + cosmo->getOmega_lambda()) != 1.0 ){
		printf("ERROR: MultiLens can only handle flat universes at present.  Must change cosmology.\n");
		exit(1);
	}

	setRedshift(zsource);

	Dl[0] = cosmo->coorDist(0,redshift[0]);
	dDl[0] = Dl[0];  // distance between jth plane and the next plane
	for(j = 1; j < Nplanes; j++){

		Dl[j] = cosmo->coorDist(0,redshift[j]);
		dDl[j] = Dl[j] - Dl[j-1]; // distance between jth plane and the next plane
	}

	cout << "Dl: ";
	for(j = 0; j < Nplanes; j++)
		cout << Dl[j] << " ";
	cout << endl;

	cout << "dDl: ";
	for(j = 0; j < Nplanes; j++)
		cout << dDl[j] << " ";
	cout << endl << endl;


	buildHaloTree(cosmo,zsource,0.01);

	if(flag_analens)
		analens->setInternalParams(cosmo,zsource);
}
