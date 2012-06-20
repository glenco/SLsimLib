/*
 * MOKAlens.cpp
 *
 *  Created on: Jun 8, 2012
 *      Author: mpetkova
 */

#include <MOKAlens.h>
#include <string>
#include <sstream>
#include <fstream>

using namespace std;

/*
 * \ingroup Constructor
 * \brief allocates and reads the MOKA map in
 */
MOKALens::MOKALens(std::string paramfile) : Lens(){
	map = new MOKAmap;

	readParamfile(paramfile);

	getDims(MOKA_input_file,&(map->nx),&(map->ny));

	map->convergence.resize(map->nx*map->ny);
	map->alpha1.resize(map->nx*map->ny);
	map->alpha2.resize(map->nx*map->ny);
	map->gamma1.resize(map->nx*map->ny);
	map->gamma2.resize(map->nx*map->ny);
	map->gamma3.resize(map->nx*map->ny);

	readImage(MOKA_input_file
			,&map->convergence
			,&map->alpha1
			,&map->alpha2
			,&map->gamma1
			,&map->gamma2
			,&(map->boxl)
			,&(map->boxlMpc)
			,&(map->zlens)
			,&(map->zsource)
			,&(map->omegam)
			,&(map->omegal)
			,&(map->h));

	map->center[0] = map->center[1] = 0.0;
	map->boxlMpc /= map->h;

	/// to radians
	map->boxl *= pi/180/3600.;
}

MOKALens::~MOKALens(){
	map->convergence.resize(0);
	map->alpha1.resize(0);
	map->alpha2.resize(0);
	map->gamma1.resize(0);
	map->gamma2.resize(0);
	map->gamma3.resize(0);
	delete map;
}

/*
 * sets the cosmology and the lens and the source according to the MOKA map parameters
 */
void MOKALens::setInternalParams(CosmoHndl cosmo, SourceHndl source){
	cosmo->setOmega_matter(map->omegam,true);
	cosmo->sethubble(map->h);
	setZlens(map->zlens);
	source->zsource = map->zsource;

	double Ds = cosmo->angDist(0,map->zsource);
	double Dl = cosmo->angDist(0,map->zlens);
	double Dls = cosmo->angDist(map->zlens,map->zsource);
	double fac = Ds/Dls/Dl;

	/// converts to the code units
	int i, j;
	for(i=0;i<map->nx;i++)
		for(j=0;j<map->ny;j++){
			int index = i+map->ny*j;
			map->convergence[index] *= fac;
			map->gamma1[index] *= fac;
			map->gamma2[index] *= fac;
		}
}


/** \ingroup ImageFinding
 * \brief Reads in a parameter file and sets up a MOKA lens map.
 *
 * Sets many parameters within the MOKA lens model
 */

void MOKALens::readParamfile(std::string filename){
  const int MAXPARAM = 50;
  string label[MAXPARAM], rlabel, rvalue;
  void *addr[MAXPARAM];
  int id[MAXPARAM];
  std::stringstream ss;
  int i, n;
  int myint;
  double mydouble;
  string mystring;
  string escape = "#";
  char dummy[300];
  int flag;

  n = 0;

  // id[] = 0 double, 1 int, 2 string

  addr[n] = &zlens;
  id[n] = 0;
  label[n++] = "z_lens";

  addr[n] = &MOKA_input_file;
  id[n] = 2;
  label[n++] = "MOKA_input_file";

  std::ifstream file_in(filename.c_str());
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
			  ss.str(std::string());

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

  set = true;

}


double MOKALens::getZlens(){
	return zlens;
}

void MOKALens::setZlens(double z){
	zlens = z;
}

/*
 * saves the image, by rading off the calues from the image tree
 * and then saving to a fits file and computing the radial profile
 * of the convergence
 */
void MOKALens::saveImage(GridHndl grid,bool saveprofiles){
	std::stringstream f;
	std::string filename;

	f << "noisy_" << MOKA_input_file;
	filename = f.str();

	MoveToTopList(grid->i_tree->pointlist);
	do{
		long index = IndexFromPosition(grid->i_tree->pointlist->current->x,map->nx,map->boxl,map->center);
		if(index > -1){

			map->convergence[index] = grid->i_tree->pointlist->current->kappa;
			map->gamma1[index] = grid->i_tree->pointlist->current->gamma[0];
			map->gamma2[index] = grid->i_tree->pointlist->current->gamma[1];
			map->gamma3[index] = grid->i_tree->pointlist->current->gamma[2];
		}
	}while(MoveDownList(grid->i_tree->pointlist)==true);

	map->boxl *= 180/pi*3600;

	writeImage(filename
			,map->convergence
			,map->gamma1
			,map->gamma2
			,map->gamma3
			,map->nx
			,map->ny
			,map->boxl
			,map->boxlMpc
			,map->zlens
			,map->zsource
			,map->omegam
			,map->omegal
			,map->h);

	if(saveprofiles == true){
		saveKappaProfile();
		saveGammaProfile();
	}
}

/*
 * computing and saving the radial profile of the convergence
 */
void MOKALens::saveKappaProfile(){

	/* measuring the differential and cumulative profile*/
	double xmin = -map->boxlMpc*0.5;
	double xmax =  map->boxlMpc*0.5;
	double drpix = map->boxlMpc/map->nx;
	std::valarray<float> pxdist(map->nx*map->ny);
	for(int i=0; i<map->nx; i++ ) for(int j=0; j<map->ny; j++ ){
		pxdist[i+map->ny*j]= sqrt(pow((xmin+(drpix*0.5)+i*drpix),2) +
				pow((xmin+(drpix*0.5)+j*drpix),2));
	}
	double dr0 = 8.*(0.5*map->boxlMpc)/(map->nx/2.);
	int nbin = int(0.5*map->boxlMpc/dr0);
	//
	std:: cout << "   " << std:: endl;
	std:: cout << " nbins = " << nbin << "  dr0 = " << dr0 << std:: endl;
	std:: cout << " ______________________________________________________ " << std:: endl;
	std:: cout << " computing profiles assuming spherical symmetry";
	// - - - - - - - - - - - - - - - - -
	double *kprofr = estprof(map->convergence,map->nx,map->ny,pxdist,dr0,xmax);
	double *sigmakprof = estsigmaprof(map->convergence,map->nx,map->ny,pxdist,dr0,xmax,kprofr);
	double *ckprofr = estcprof(map->convergence,map->nx,map->ny,pxdist,dr0,xmax);
	double *sigmackprof = estsigmacprof(map->convergence,map->nx,map->ny,pxdist,dr0,xmax,kprofr);
	std::ostringstream fprof;
	fprof << "MAP_radial_prof_" << MOKA_input_file << ".dat";
	std:: ofstream filoutprof;
	std:: string filenameprof = fprof.str();
	filoutprof.open(filenameprof.c_str());
	filoutprof <<"# r      kappa     sig_k     ckappa     sig_ck" << std:: endl;
	int l;
	for(l=0;l<nbin;l++){
		filoutprof << dr0*l + dr0/2. << "  "
				<< kprofr[l] << "  " << sigmakprof[l] << "  "
				<< ckprofr[l] << "  " << sigmackprof[l] <<
				std:: endl;
	}
	filoutprof.close();
}

void MOKALens::saveGammaProfile(){

	/* measuring the differential and cumulative profile*/
	double xmin = -map->boxlMpc*0.5;
	double xmax =  map->boxlMpc*0.5;
	double drpix = map->boxlMpc/map->nx;
	std::valarray<float> pxdist(map->nx*map->ny);
	for(int i=0; i<map->nx; i++ ) for(int j=0; j<map->ny; j++ ){
		pxdist[i+map->ny*j]= sqrt(pow((xmin+(drpix*0.5)+i*drpix),2) +
				pow((xmin+(drpix*0.5)+j*drpix),2));
	}
	double dr0 = 8.*(0.5*map->boxlMpc)/(map->nx/2.);
	int nbin = int(0.5*map->boxlMpc/dr0);
	//
	std:: cout << "   " << std:: endl;
	std:: cout << " nbins = " << nbin << "  dr0 = " << dr0 << std:: endl;
	std:: cout << " ______________________________________________________ " << std:: endl;
	std:: cout << " computing profiles assuming spherical symmetry";
	// - - - - - - - - - - - - - - - - -
	double *kprofr = estprof(map->gamma2,map->nx,map->ny,pxdist,dr0,xmax);
	double *sigmakprof = estsigmaprof(map->gamma2,map->nx,map->ny,pxdist,dr0,xmax,kprofr);
	double *ckprofr = estcprof(map->gamma2,map->nx,map->ny,pxdist,dr0,xmax);
	double *sigmackprof = estsigmacprof(map->gamma2,map->nx,map->ny,pxdist,dr0,xmax,kprofr);
	std::ostringstream fprof;
	fprof << "MAP_radial_prof_" << MOKA_input_file << ".dat";
	std:: ofstream filoutprof;
	std:: string filenameprof = fprof.str();
	filoutprof.open(filenameprof.c_str());
	filoutprof <<"# r      kappa     sig_k     ckappa     sig_ck" << std:: endl;
	int l;
	for(l=0;l<nbin;l++){
		filoutprof << dr0*l + dr0/2. << "  "
				<< kprofr[l] << "  " << sigmakprof[l] << "  "
				<< ckprofr[l] << "  " << sigmackprof[l] <<
				std:: endl;
	}
	filoutprof.close();
}

/** \ingroup DeflectionL2
   *
   * \brief Routine for obtaining the deflection and other lensing quantities for
   * a MOKA map (MOKALens), for just one ray!!
   *
*/
void MOKALens::rayshooterInternal(double *xx, double *alpha, double *gamma, double *kappa, bool kappa_off){

	long index = IndexFromPosition(xx,map->nx,map->boxlMpc,map->center);

	if(index > -1){
		alpha[0] = map->alpha1[index];
		alpha[1] = map->alpha2[index];
		gamma[0] = map->gamma1[index];
		gamma[1] = map->gamma2[index];
		gamma[2] = 0.0;
		*kappa = map->convergence[index];
	}
	else{
		alpha[0] = alpha[1] = 0.0;
		gamma[0] = gamma[1] = gamma[2] = 0.0;
		*kappa = 0.0;
	}


	return;
}
