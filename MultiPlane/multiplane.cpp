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

    HALO ha(cosmo,lens->mass_resolution,0.);

    int iterator;
    std::vector<double> vmasses,vsizes,vredshifts;
    std::vector<int> vindex;
    Logm.resize(Nmassbin);
    Nhaloes.resize(Nmassbin);

	std:: vector<double> Dli;
	/* fill the log(mass) vector */
	fill_linear(Logm,Nmassbin,9.0,MaxLogm);

	int Nplanes = lens->getNplanes();

	for(int i=0;i<Nplanes;i++){

		double Nhaloestot;
		float Nhaloestotf;

		double z1, z2;

		if(i == 0) z1 = 0.0;
		else z1 = 0.5*(lens->redshift[i] - lens->redshift[i-1]);

		if(i == Nplanes) z2 = zsource;
		else z2 = 0.5*(lens->redshift[i+1] - lens->redshift[i]);

		for(int k=0;k<Nmassbin;k++){
			// cumulative number density in one square degree
			Nhaloes[k]=cosmo->number(pow(10,Logm[k]),z1,z2,mfty)*fieldofview;
			if(k == 0) Nhaloestot = Nhaloes[k];
			// normalize the cumulative distribution to one
			Nhaloes[k] = Nhaloes[k]/Nhaloestot;
		}
		Nhaloestotf=poidev(float(Nhaloestot), &seed);
		lens->NhalosinPlane[i] = Nhaloestotf;
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
		}
	}
	N = vmasses.size();

	masses = new float[N];
	sizes = new float[N];
	redshifts = new float[N];

	double r,theta,maxr;
	pos = PosTypeMatrix(0,N-1,0,2);
	for(int i = 0; i < N; i++){
		maxr = sqrt(fieldofview/M_PI)*Dli[i]*M_PI/180;
		r = maxr*ran2(&seed);
		theta=2*pi*ran2(&seed);
		pos[i][0] = r*cos(theta);
		pos[i][1] = r*sin(theta);
		pos[i][2] = 0.0;

		masses[i] = vmasses[i];
		sizes[i] = vsizes[i];
		redshifts[i] = vredshifts[i];

		cout << redshifts[i] << " " << masses[i] << endl;
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

	redshift = new double[Nplanes+1];

	Dl = new double[Nplanes+1];

	NhalosinPlane = new IndexType[Nplanes];

	dDl = new double[Nplanes+1];

	charge = mass_scale/4/pi/Grav;

	halo_tree = new ForceTreeHndl[Nplanes];
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

	  addr[n] = &zlens;
	  id[n] = 0;
	  label[n++] = "z_lens";

	  addr[n] = &mass_resolution;
	  id[n] = 0;
	  label[n++] = "mass_resolution";

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

	  printMultiLens();
}



void MultiLens::printMultiLens(){
	cout << endl << "outputfile " << outputfile << endl;

	cout << endl << "**multi lens model**" << endl;

	cout << "Nplanes " << Nplanes << endl;

	cout << "z_lens " << zlens << endl << endl;

	cout << "mass_resolution " << mass_resolution << endl << endl;
}

MultiLens::~MultiLens(){
	delete[] halo_tree;
	delete[] NhalosinPlane;
	delete[] Dl;
	delete[] redshift;
	delete[] dDl;
}

void buildHaloTree(MultiLens *lens,CosmoHndl cosmo, double zsource,double fieldofview){
	IndexType N, N_last;
	haloHndl halo;

    halo = new haloM(zsource,cosmo,lens,fieldofview,1);

	for(int j = 0, N_last = 0; j < lens->getNplanes(); j++){
		N = lens->NhalosinPlane[j];
		lens->halo_tree[j] = new ForceTree(&halo->pos[N_last + N],N,&halo->masses[N_last + N],&halo->sizes[N_last + N],true,true,5,2,true,0.1);

		N_last = N;
	}

}

void MultiLens::setRedshift(double zsource){
	std:: vector<double> lz;
	/* fill the redshift vector logarithmically */
	fill_linear(lz,Nplanes+2,0.,zsource);

	for(int i = 1; i < Nplanes+2; i++)
		redshift[i-1] = -1 + pow(10.,lz[i]);
}

void MultiLens::setInternalParams(CosmoHndl cosmo, double zsource){
	int j;

	mass_scale = 1.0;

	setRedshift(zsource);

	zlens = redshift[0];

	for(j = 0; j < Nplanes; j++){
		Dl[j] = cosmo->angDist(0,redshift[j]);
		dDl[j] = cosmo->angDist(redshift[j],redshift[j+1]);  // distance between jth plane and the next plane
	}
}

/*
void swap(float *a,float *b){
	float tmp;
	tmp=*a;
	*a=*b;
	*b=tmp;
}
/*
/*void swap(PosType *a,PosType *b){
	PosType tmp;
	tmp=*a;
	*a=*b;
	*b=tmp;
}*/
/*
void swap(IndexType a,IndexType b){
	IndexType tmp;
	tmp=a;
	a=b;
	b=tmp;
}
void swap(IndexType *a,IndexType *b){
	IndexType tmp;
	tmp=*a;
	*a=*b;
	*b=tmp;
}

void quicksort(IndexType *particles,float *redshifts,PosType **pos,float *sizes,float *masses,IndexType N){
	double pivotvalue;
	unsigned long pivotindex,newpivotindex,i;

	if(N <= 1) return ;

	// pick pivot as the median of the first, last and middle values
	if ((redshifts[0] >= redshifts[N/2] && redshifts[0] <= redshifts[N-1])
			|| (redshifts[0] >= redshifts[N-1] && redshifts[0] <= redshifts[N/2])) pivotindex = 0;
	else if ((redshifts[N/2] >= redshifts[0] && redshifts[N/2] <= redshifts[N-1])
			|| (redshifts[N/2] >= redshifts[N-1] && redshifts[N/2] <= redshifts[0])) pivotindex = N/2;
	else pivotindex = N-1;
	pivotvalue=redshifts[pivotindex];

	// move pivet to end of array
	swap(&redshifts[pivotindex],&redshifts[N-1]);
	swap(&particles[pivotindex],&particles[N-1]);
	swap(pos[pivotindex],pos[N-1]);
	swap(&sizes[pivotindex],&sizes[N-1]);
	swap(&masses[pivotindex],&masses[N-1]);

	newpivotindex=0;

	// partition list and array
	for(i=0;i<N;++i){
		if(redshifts[i] <= pivotvalue){
			swap(&redshifts[newpivotindex],&redshifts[i]);
			swap(&particles[newpivotindex],&particles[i]);
			swap(pos[pivotindex],pos[i]);
			swap(&sizes[pivotindex],&sizes[i]);
			swap(&masses[pivotindex],&masses[i]);
			++newpivotindex;
		}
	}
	--newpivotindex;

	quicksort(particles,redshifts,pos,sizes,masses,newpivotindex);
	quicksort(&particles[newpivotindex+1],&redshifts[newpivotindex+1],&pos[newpivotindex+1],&sizes[newpivotindex+1],&masses[newpivotindex+1],N-newpivotindex-1);

	return ;
}
*/
