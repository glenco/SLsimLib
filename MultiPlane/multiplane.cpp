/*
 * halo_model.cpp
 *
 *  Created on: Jan 14, 2012
 *      Author: mpetkova
 */

#include <slsimlib.h>

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

multiLens::multiLens(char filename[]) : Lens(filename){
	redshift = new double[Nplanes];
	Ds = new double[Nplanes];
	Dl = new double[Nplanes];
	Dls = new double[Nplanes];
	Sigma_crit = new double[Nplanes];
	halo_tree = new ForceTreeHndl[Nplanes];
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

void multiLens::buildHaloTree(){
	IndexType N, N_last;
	double dz;
	int i, j, n;

    halo = new haloM(0.4, 3.64);

    dz = redshift[1] - redshift[0];

	for(j = 0, N_last = 0; j < Nplanes; j++){

		for(i = 0, N = 0; i < halo->N; i++)
			if(halo->redshifts[i] >= redshift[j] && halo->redshifts[i] < (redshift[j]+dz))
				N++;

		halo_tree[j] = new ForceTree(&halo->pos[N_last + N],N,&halo->masses[N_last + N],&halo->sizes[N_last + N],true,true,5,2,true,0.1);

		N_last = N;
	}
}

void multiLens::setInternalParams(CosmoHndl cosmo, double zsource){
	int j;

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
