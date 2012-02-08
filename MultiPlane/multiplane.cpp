/*
 * halo_model.cpp
 *
 *  Created on: Jan 14, 2012
 *      Author: mpetkova
 */

#include <slsimlib.h>
#include <sstream>

haloM::haloM(double maxr, double zmax){
	int i, n;
	long seed;
	double r, theta;

	seed = 2203;

	N = 100;

	masses = new float[N];
	sizes = new float[N];
	redshifts = new float[N];
	index = new IndexType[N];

	pos = PosTypeMatrix(0,N-1,0,2);

	for(i = 0; i < N; i++){
		r = maxr*ran2(&seed);
		theta=2*pi*ran2(&seed);
		pos[i][0] = r*cos(theta);
		pos[i][1] = r*sin(theta);
		pos[i][2] = 0.0;

		redshifts[i] = ran2(&seed) * zmax;

		sizes[i] = ran2(&seed) * 0.3;

		masses[i] = ran2(&seed);

		index[i] = i;
	}

	quicksort(index,redshifts,pos,sizes,masses,N);
}

haloM::~haloM(){
	free_PosTypeMatrix(pos,0,N-1,0,2);
	delete[] redshifts;
	delete[] sizes;
	delete[] masses;
}

multiLens::multiLens(string filename) : Lens(){
	readParamfile(filename);

	redshift = new double[Nplanes];
	Ds = new double[Nplanes];
	Dl = new double[Nplanes];
	Dls = new double[Nplanes];
	Sigma_crit = new double[Nplanes];
	halo_tree = new ForceTreeHndl[Nplanes];
}

void multiLens::readParamfile(string filename){
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

	  cout << "multi lens: reading from " << filename << endl;

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
		  if(id[i] > 0){
			  ERROR_MESSAGE();
			  cout << "parameter " << label[i] << " needs to be set!" << endl;
			  exit(0);
		  }
	  }

	  file_in.close();

	  printMultiLens();
}



void multiLens::printMultiLens(){
	cout << endl << "outputfile "<< outputfile << endl;

	cout << endl << "**multi lens model**" << endl;

	cout << "Nplanes " << Nplanes << endl;
	cout << "z_lens " << zlens << endl << endl;
}

multiLens::~multiLens(){
	delete halo;
	delete[] halo_tree;
	delete[] Sigma_crit;
	delete[] Dls;
	delete[] Dl;
	delete[] Ds;
	delete[] redshift;
}

void buildHaloTree(multiLens *lens){
	IndexType N, N_last, Nplanes;
	double dz;
	int i, j, n;

    lens->halo = new haloM(0.4, 3.64);

    dz = lens->redshift[1] - lens->redshift[0];

    Nplanes = lens->getNplanes();

	for(j = 0, N_last = 0; j < Nplanes; j++){

		for(i = 0, N = 0; i < lens->halo->N; i++)
			if(lens->halo->redshifts[i] >= lens->redshift[j] && lens->halo->redshifts[i] < (lens->redshift[j]+dz))
				N++;

		lens->halo_tree[j] = new ForceTree(&lens->halo->pos[N_last + N],N,&lens->halo->masses[N_last + N],&lens->halo->sizes[N_last + N],true,true,5,2,true,0.1);

		N_last = N;
	}
}

void multiLens::setRedshift(double zsource){
	int i;
	double dz;

	dz = zsource / (Nplanes+1);

	for(i = 0; i < Nplanes; i++)
		redshift[i] = (i + 1)*dz;
}

void multiLens::setInternalParams(CosmoHndl cosmo, double zsource){
	int j;

	setRedshift(zsource);

	mass_scale = 1.0;

	for(j = 0; j < Nplanes; j++){
		Ds[j] = cosmo->angDist(0,zsource);
		Dl[j] = cosmo->angDist(0,redshift[j]);
		Dls[j] = cosmo->angDist(redshift[j],zsource);

		Sigma_crit[j]=cosmo->angDist(0,zsource)/cosmo->angDist(redshift[j],zsource)
					/cosmo->angDist(0,redshift[j])/4/pi/Grav;
	}
}

void swap(float *a,float *b){
	float tmp;
	tmp=*a;
	*a=*b;
	*b=tmp;
}
/*void swap(PosType *a,PosType *b){
	PosType tmp;
	tmp=*a;
	*a=*b;
	*b=tmp;
}*/
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
